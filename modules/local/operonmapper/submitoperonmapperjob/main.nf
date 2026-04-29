process OPERONMAPPER_SUBMITOPERONMAPPERJOB {
    errorStrategy "ignore"
    tag "${meta.id}"
    label 'process_medium'

    container "oras://community.wave.seqera.io/library/curl_requests:e1dc553ba087bf18"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path('*')
    val email

    output:
    tuple val(meta), path("${prefix}_operons.gff"), emit: gff, optional: true
    tuple val(meta), path("${prefix}_operons.tsv"), emit: tsv, optional: true

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    # Extract CDS from GFF for submission
    awk -F'\\t' '\$3 == "CDS"{ OFS="\\t"; print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9 }' ${prefix}.gff* > tmpcds.gff

    # Remove > and < symbols from columns 4 and 5
    awk -F'\\t' '
    BEGIN { OFS="\\t" }
    {
        gsub(/^[><]/, "", \$4)
        gsub(/^[><]/, "", \$5)
        print
    }
    ' tmpcds.gff > ${prefix}_clean_CDS.gff

python <<'PYCODE'
import requests, re

strain = "${prefix}"
email = "${email}"

session = requests.Session()

headers = {
    "User-Agent": "Mozilla/5.0 (X11; Linux x86_64; rv:133.0) Gecko/20100101 Firefox/133.0",
    "Referer": "https://biocomputo.ibt.unam.mx/operon_mapper/",
    "Origin": "https://biocomputo.ibt.unam.mx",
}

files = {
    'fastafile': open(f"{strain}.fna", 'rb'),
    'gfffile': open(f"{strain}_clean_CDS.gff", 'rb')
}

userdata = {
    "email1": email,
    "descri": strain + " operons",
    "genepairs": "si",
    "operons": "si",
    "cogs": "si",
    "orfsdescri": "si"
}

url = "https://biocomputo.ibt.unam.mx/operon_mapper/capta_forma_01.pl"
resp = session.post(url, headers=headers, data=userdata, files=files, verify=False)
resp.raise_for_status()

# Extract the results URL (out_XXXXX.html)
m = re.search(r'(https?://biocomputo\\.ibt\\.unam\\.mx/operon_mapper/out/out_\\d+\\.html)', resp.text)
if not m:
    raise RuntimeError("Could not locate results URL from submission.")

results_url = m.group(1)

# Write the results URL and status to file for the bash wait-loop
with open("operon-mapper_results_url", "w") as f:
    f.write(str(resp.status_code) + "\\n")
    f.write(results_url + "\\n")

print("Operon Mapper job submitted successfully. Waiting for results in bash loop.")
PYCODE

    # --- Bash wait loop ---
    max_wait_time=120  # minutes
    sleep_interval=5   # minutes

    URL=`tail -n1 operon-mapper_results_url`
    submitstatus=`head -n1 operon-mapper_results_url`
    wait_count=0

    sleep 10s

    while true; do
        response=\$(curl -s -L \$URL)
        if echo "\$response" | grep -qF 'Download'; then
            echo "Operon Mapper job finished. Downloading results..."
            genepairsurl=\$(echo "\$response" | grep 'operonic_gene_pairs_' | sed 's/^.*href="//g;s/".*//g;s#^\\.\\.#https://biocomputo.ibt.unam.mx/operon_mapper#g')
            operonsurl=\$(echo "\$response" | grep 'list_of_operons_' | sed 's/^.*href="//g;s/".*//g;s#^\\.\\.#https://biocomputo.ibt.unam.mx/operon_mapper#g')
            curl -s -L -o ${prefix}_operons.gff "\$operonsurl"
            curl -s -L -o ${prefix}_operons.tsv "\$operonsurl"
            echo "Operon Mapper files downloaded."
            break
        elif echo "\$response" | grep -iqF 'errors'; then
            echo "Error in Operon Mapper job. Exiting."
            touch ${prefix}_operons.tsv
            break
        else
            wait_count=\$((wait_count+sleep_interval))
            if [[ \$wait_count -ge \$max_wait_time ]]; then
                echo "Exceeded maximum wait time (\$max_wait_time minutes). Skipping Operon Mapper."
                echo "Check your email(${email}) for results."
                touch ${prefix}_operons.tsv
                break
            else
                echo "Job still running. Waiting \${sleep_interval} minutes..."
                sleep \${sleep_interval}m
                continue
            fi
        fi
    done
    """
}