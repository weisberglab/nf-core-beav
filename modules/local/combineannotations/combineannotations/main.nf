process COMBINE_ANNOTATIONS{
    tag "${meta.id}"
    label 'process_medium'
    
    publishDir "${params.outdir}/${meta.id}_results/final_annotation_files", mode: 'copy'

    container "oras://community.wave.seqera.io/library/bcbio-gff_pip_biopython:1c8350bea1249708"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path('*')

    output:
    tuple val(meta), path("${prefix}_final.gbk"), emit: finalgbk
    tuple val(meta), path("${prefix}_final.faa"), emit: finalfaa
    tuple val(meta), path("${prefix}_final.ffn"), emit: finalffn
    tuple val(meta), path("${prefix}_annotation_summary.html"), emit: summary
    val(true), emit: done

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mv ${prefix}.faa ${prefix}_final.faa
    mv ${prefix}.ffn ${prefix}_final.ffn
    if [[ -f "${prefix}.gbff" ]]; then
	    mv ${prefix}.gbff ${prefix}.gbk
    fi

    echo "Combining annotations for sample ${prefix}"

    # Cut macsyfinder table
    if [[ -f "${prefix}_best_solution.tsv" ]]; then
        cut -f 2,3 ${prefix}_best_solution.tsv \\
        | grep -v '^#' \\
        | grep -v '^\$' \\
        > ${prefix}_macsyfinder.tsv
    fi

    # Reorder nhmmer columns, replace + with +1 and - with -1 in last column(5). Then, go through each line and swap start and end position numbers if end > start
    if [[ -f "${prefix}.nhmmer.all.tsv" ]]; then
        awk -F'\\t' -v OFS='\\t' '
        {
            contig = \$1
            start = \$2
            end = \$3
            type = \$7
            strand = \$4
            
            # Replace + with +1 and - with -1
            strand = (strand == "+") ? "+1" : (strand == "-") ? "-1" : strand
            
            # Make start < end
            if (start > end) {
                tmp = start
                start = end
                end = tmp
            }

            print contig, start, end, type, strand
        }' ${prefix}.nhmmer.all.tsv > ${prefix}_nhmmer_reduced.tsv 
    fi

    # Reorder phispy columns
    if [[ -f "${prefix}_phispy.tsv" ]]; then
        awk -F'\\t' -v OFS='\\t' '{ print \$2, \$3, \$4, \$1, \$11 }' ${prefix}_phispy.tsv > ${prefix}_phispy_reduced.tsv
    fi

    # Cut defensefinder file
    if [[ -f "${prefix}_defense_finder_genes.tsv" ]]; then
        cut -f 2,3 ${prefix}_defense_finder_genes.tsv > ${prefix}_defense_reduced.tsv
    fi

    # Mobsuite file name change for consistency
    if [[ -f "contig_report.txt" ]]; then
        cp contig_report.txt ${prefix}_contig_report.txt
    fi
    
python - <<'PYCODE'

import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
#from BCBio import GFF
from html import escape
import os
import os.path

#Check if input exists
bakta_path = f"${prefix}.gbk"
macsyfinder_path = f"${prefix}_macsyfinder.tsv"
defensefinder_path = f"${prefix}_defense_reduced.tsv"
hmmdb_path = f"${prefix}_nhmmer_reduced.tsv"
phispy_path = f"${prefix}_phispy_reduced.tsv"
antismash_path = f"${prefix}_antismash_final.tsv.beav.subset"
tiger_path = f"${prefix}_TIGER2_final.table.out"
integron_path = f"${prefix}.integrons"
integron_gene_path = f"${prefix}.summary"
gapmind_path = f"${prefix}_combined_GapMind_results.tab"
operon_path = f"${prefix}_operons.tsv"
oriT_path = f"${prefix}_oriT.table"
dif_path = f"${prefix}_dif.table"
mobsuite_path = f"${prefix}_contig_report.txt"
ani_path = f"${prefix}.ani.txt"
outfile_path = f"${prefix}_final.gbk"
faa_path = f"${prefix}_final.faa"
ffn_path = f"${prefix}_final.ffn"

#read in annotation tables as dictionary
locus_dict = {}
if os.path.isfile(macsyfinder_path) == True:
    with open (f"{macsyfinder_path}", "r") as f:
        for line in f:
            (key,values) = line.split()
            locus_dict[key] = values 

plasmid_contigs = set()
if os.path.isfile(mobsuite_path):
    with open(mobsuite_path) as mob_table:
        header = mob_table.readline().rstrip().split('\\t')
        idx_contig = header.index("contig_id")
        idx_type = header.index("molecule_type")

        for line in mob_table:
            cols = line.rstrip().split('\\t')
            contig = cols[idx_contig].split()[0]
            mol_type = cols[idx_type]
            if mol_type == "plasmid":
                plasmid_contigs.add(contig)

oriT = {}
if os.path.isfile(oriT_path) == True:
    with open (f"{oriT_path}", "r") as ori_table:
         for line in ori_table:
             contig,annot,start,end,strand,reference = line.strip().split('\\t')
             if contig in oriT:
                 oriT[contig].append((contig,annot,start,end,strand,reference))
             else:
                 oriT[contig] = [(contig,annot,start,end,strand,reference)]

defense_dict = {}
if os.path.isfile(defensefinder_path) == True:
    with open (f"{defensefinder_path}", 'r') as k:
        for line in k:
            (key,values) = line.split()
            defense_dict[key] = values

gapmind_dict = {}
if os.path.isfile(gapmind_path) == True:
    with open (gapmind_path, 'r') as gp:
        for line in gp:
            (key,values) = line.rstrip().split("\\t")
            gapmind_dict[key] = values

operon_dict = {}
if os.path.isfile(operon_path) == True:
    with open (operon_path, 'r') as opp:
        for line in opp:
            (key,values) = line.rstrip().split("\\t")
            operon_dict[key] = values

hmmdb = {}
if os.path.isfile(hmmdb_path) == True:
    with open(f"{hmmdb_path}", 'r') as hmmtable_file:
        for l in hmmtable_file:
            replicon,start,end,annot,strand = l.strip('\\n').split('\\t')
            if replicon in hmmdb:
                hmmdb[replicon].append((replicon,start,end,annot,strand))
            else:
                hmmdb[replicon]=[(replicon,start,end,annot,strand)]

phispy = {}
if os.path.isfile(phispy_path) == True:     
    with open (f"{phispy_path}", 'r') as prophage_file:
        for l in prophage_file:
            replicon,start,end,annot,category = l.strip().split('\\t')
            if replicon in phispy:
                phispy[replicon].append((start,end,annot,category))
            else:
                phispy[replicon] = [(start,end,annot,category)] 

antismash_dict = {}
if os.path.isfile(antismash_path) == True:
    with open (f"{antismash_path}", 'r') as antismash_file:
        next(antismash_file) # Skip header
        for line in antismash_file:
            locus,cluster,function,nrps,domain = line.rstrip('\\n').split('\\t')
            if locus in antismash_dict:
                antismash_dict[locus].append((cluster,function,nrps,domain))
            else:
                antismash_dict[locus] = [(cluster,function,nrps,domain)]

tiger_dict = {}
if os.path.isfile(tiger_path) == True:          
    with open (f"{tiger_path}", 'r') as tiger_file:
        for line in tiger_file:
            replicon,start,end,annot = line.strip().split(None, 3)
            if replicon in tiger_dict:
                tiger_dict[replicon].append((start,end,annot))
            else:
                tiger_dict[replicon] = [(start,end,annot)]

integron = {}
if os.path.isfile(integron_path) == True:                
    with open(integron_path, 'r') as file:
        for l in file:
            if l.startswith("#"):
                continue
            cols = l.strip().split("\\t")
            if len(cols) < 6:
                continue  # skip summary or "No integron found" lines

            replicon, element, start, end, strand, complete = cols

            if replicon not in integron:
                integron[replicon] = []
            integron[replicon].append((element, start, end, strand, complete))

integron_gene = {}
if os.path.isfile(integron_gene_path) == True:
    protein = "IntI"
    with open (f"{integron_gene_path}", 'r') as integron_gene_file:
        for locus in integron_gene_file:
            key = locus.strip()
            integron_gene[key] = protein    
               
#Parse genbank and match locus tags to annotation                        
new_records = []
for record in SeqIO.parse(f"{bakta_path}","gb"):
    accession = record.annotations.get("accessions")
    for feature in record.features:
        locus_tags = feature.qualifiers.get("locus_tag")
        if locus_tags is not None and feature.type == "CDS":
            if locus_tags[0] in locus_dict:
                if "note" in feature.qualifiers:
                    feature.qualifiers["note"].append("MacSyFinder: " + locus_dict[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = ["MacSyFinder: " + locus_dict[locus_tags[0]]]                

            if locus_tags[0] in defense_dict:
                if "note" in feature.qualifiers:
                    feature.qualifiers["note"].append("DefenseFinder: " + defense_dict[locus_tags[0]])
                else: 
                    feature.qualifiers["note"] = ["DefenseFinder: " + defense_dict[locus_tags[0]]]
                    
            if locus_tags[0] in antismash_dict:
                cluster = antismash_dict[locus_tags[0]][0][0]
                gene_product = antismash_dict[locus_tags[0]][0][1]
                NRPS_PKS = antismash_dict[locus_tags[0]][0][2]
                domains = antismash_dict[locus_tags[0]][0][3]
                if cluster != "":
                    if  "note" in feature.qualifiers:
                        feature.qualifiers["note"].append("antiSMASH cluster: " +  cluster)
                    else:
                        feature.qualifiers["note"] = ["antiSMASH cluster: " + cluster]
                if gene_product != "":
                    if "note" in feature.qualifiers:
                        feature.qualifiers["note"].append("antiSMASH gene product: " +  gene_product)
                    else:
                        feature.qualifiers["note"] = ["antiSMASH gene product: " + gene_product]

                if NRPS_PKS != "":
                    if "note" in feature.qualifiers:
                        feature.qualifiers["note"].append("antiSMASH NRPS/PKS: " +  NRPS_PKS)
                    else:
                        feature.qualifiers["note"] = ["antiSMASH NRPS/PKS: " + NRPS_PKS]

                if domains != "":
                    if "note" in feature.qualifiers:
                        feature.qualifiers["note"].append("antiSMASH domain: " +  domains)
                    else:
                        feature.qualifiers["note"] = ["antiSMASH domain: " + domains]

            if locus_tags[0] in integron_gene:
                if "note" in feature.qualifiers:
                    feature.qualifiers["note"].append("IntegronFinder: " + integron_gene[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = ["IntegronFinder: " + integron_gene[locus_tags[0]]]

            if locus_tags[0] in gapmind_dict:
                if "note" in feature.qualifiers:
                    feature.qualifiers["note"].append(gapmind_dict[locus_tags[0]])
                else:
                    feature.qualifiers["note"] = [gapmind_dict[locus_tags[0]]]
            
            if locus_tags[0] in operon_dict:
                feature.qualifiers["operon"] = operon_dict[locus_tags[0]]

#adding new features
    if record.id in plasmid_contigs:
        plasmid_feat = SeqFeature(
            FeatureLocation(0, len(record.seq)),
            type="plasmid",
            qualifiers={
                "note": ["MOB-suite predicted plasmid"],
                "inference": ["MOB-suite contig_report"]
            }
        )
        record.features.append(plasmid_feat)

    if record.id in oriT:
        for oriT_annot in oriT[record.id]:
            oriT_new_feat = SeqFeature((FeatureLocation(int(oriT_annot[2]), int(oriT_annot[3]), strand = int(oriT_annot[4]))), type="oriT", qualifiers = {"reference": [oriT_annot[5]], "inference" : "blastn"})
            record.features.append(oriT_new_feat)
    
    if record.id in hmmdb:
        for current_border in hmmdb[record.id]:
            new_feat = SeqFeature((FeatureLocation(int(current_border[1]), int(current_border[2]), strand = int(current_border[4]))),type="misc_feature", qualifiers= {"note": [current_border[3]], "inference" : "BEAV"})
            record.features.append(new_feat)
        
    if record.id in integron:
        for integrons in integron[record.id]:
            integron_newfeat = SeqFeature(FeatureLocation(int(integrons[1]), int(integrons[2])), type="misc_feature", qualifiers= {"note": [integrons[0]], "inference" : "MacSyFinder TXSS models"}, strand = int(integrons[3]))
            integron_newfeat.qualifiers["note"].append("Integron Status: " +integrons[4])
            record.features.append(integron_newfeat)

    if record.id in tiger_dict:
        for ice in tiger_dict[record.id]:
            if "Phage" in ice[2]:
                tiger_newfeat = SeqFeature(FeatureLocation(int(ice[0]),int(ice[1])), type="mobile_element", qualifiers={"mobile_element_type": "phage", "note": ice[2] , "inference" : "TIGER2"})
            else:
                tiger_newfeat = SeqFeature(FeatureLocation(int(ice[0]),int(ice[1])), type="mobile_element", qualifiers={"mobile_element_type": "integrative element", "note": ice[2] , "inference" : "TIGER2"})
            record.features.append(tiger_newfeat)
    for number in accession:
        if number in phispy:
            for phage in phispy[number]:
                newfeat = SeqFeature(FeatureLocation(int(phage[0]),int(phage[1])), type="mobile_element", qualifiers={"mobile_element_type": "phage", "note": phage[2] + " " + phage[3] , "inference" : "PHISPY"})
                record.features.append(newfeat)
    record.features.sort(key=lambda x: x.location.start,reverse=False)  
    new_records.append(record)
output_gbk_handle = open(outfile_path, 'w')
SeqIO.write(new_records,output_gbk_handle, "genbank")
output_gbk_handle.close()

#writing output into different formats
output_faa_handle = open(faa_path, 'w')
for record in new_records:
        for feature in record.features:
                if feature.type=="CDS" and 'pseudogene' not in feature.qualifiers:
                        product = feature.qualifiers.get("product")[0]
                        gene = feature.qualifiers.get("gene")
                        if gene is not None:
                                gene = feature.qualifiers.get("gene")[0]
                        else:
                                gene = ""

                        translation = feature.qualifiers.get("translation")
                        if translation is not None:
                                translation = feature.qualifiers.get("translation")[0]
                        else:
                                translation = ""
                        output_faa_handle.write(">%s %s %s\\n%s\\n" % (
                                feature.qualifiers['locus_tag'][0],
                                product,
                                gene,
                                translation))
output_faa_handle.close()


output_ffn_handle = open(ffn_path, 'w')

for record in new_records:
        for feature in record.features:
                if feature.type=="CDS":
                        product = feature.qualifiers.get("product")[0]
                        gene = feature.qualifiers.get("gene")
                        if gene is not None:
                                gene = feature.qualifiers.get("gene")[0]
                        else:
                                gene = ""

                        locus = feature.qualifiers.get("locus_tag")
                        if locus is not None:
                                locus = feature.qualifiers.get("locus_tag")[0]
                        else:
                                locus = ""
                        nucleotide = feature.extract(record).seq
                        locus = feature.qualifiers.get("locus_tag")[0]
                        output_ffn_handle.write(">%s %s %s\\n%s\\n" % (
                                locus,
                                product,
                                gene,
                                nucleotide))
output_ffn_handle.close()

def get_feature_label(feature):
    q = feature.qualifiers

    # Mobile elements
    if feature.type == "mobile_element":
        met = q.get("mobile_element_type", [])
        if met:
            if "phage" in met[0].lower():
                return "Prophage"
            if "integrative" in met[0].lower():
                return "Integrative element"
        return "Mobile element"

    # Plasmid
    if feature.type == "plasmid":
        return "Plasmid"

    # antiSMASH
    for n in q.get("note", []):
        if "antismash cluster" in n.lower():
            return "Secondary metabolite cluster"

    # Integrons
    for n in q.get("note", []):
        if "integron" in n.lower():
            return "Integron"

    # Defense systems
    for n in q.get("note", []):
        if "defensefinder" in n.lower():
            return "Defense system"

    # Transposition
    for n in q.get("product", []) + q.get("note", []):
        v = n.lower()
        if "transposase" in v:
            return "Transposase"
        if "integrase" in v:
            return "Integrase"
        if "recombinase" in v:
            return "Recombinase"

    # CDS naming
    if "product" in q:
        return q["product"][0]

    if "gene" in q:
        return q["gene"][0]

    if "locus_tag" in q:
        return q["locus_tag"][0]

    return feature.type

# create an output HTML with a summary of the annotations
html_path = f"${prefix}_annotation_summary.html"

highlight_features = [
    "mobile_element", "CDS", "rep_origin", "oriT", "misc_feature"
]

highlight_keywords = [
    "phage", "integrative element", "T3SS", "T4SS", "T6SS",
    "vir", "tra", "trb", "repA", "repB", "repC",
    "acs", "ags", "chsA", "cus", "mas1", "mas2", "nos", "ocs", "sus", "vis", "unk",
    "accA", "accB", "accC", "accD", "accE", "accF", "accF_2", "accG", "accG_2", "accR",
    "integrase", "transposase", "recombinase", "SMASH", "secondary metabolite", "biosynthetic",
    "CRISPR", "restriction-modification", "RM system", "anti-phage", "repeat_region",
    "vgrG"
]

label = get_feature_label(feature)

html_lines = [
    "<html><head><style>",
    "table {border-collapse: collapse; width: 100%;}",
    "th, td {border: 1px solid black; padding: 4px; text-align: left; font-family: Arial, sans-serif; font-size: 12px;}",
    "th {background-color: #f2f2f2;}",
    "</style></head><body>",
    "<h2>Highlighted Genomic Features</h2>",
    "<table>",
    "<tr><th>Contig</th><th>Start</th><th>End</th><th>Strand</th><th>Feature Type</th><th>Annotation</th></tr>"
]

for record in SeqIO.parse(f"{outfile_path}", "genbank"):
    for feature in record.features:

        if feature.type not in highlight_features:
            continue

        # Gather searchable text
        notes = " ".join(feature.qualifiers.get("note", []))
        genes = " ".join(feature.qualifiers.get("gene", []))
        products = " ".join(feature.qualifiers.get("product", []))
        searchable = " ".join([notes, genes, products]).lower()

        if not any(k.lower() in searchable for k in highlight_keywords):
            continue

        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand

        label = get_feature_label(feature)

        # Full text for hover
        hover = " | ".join(
            feature.qualifiers.get("note", []) +
            feature.qualifiers.get("product", []) +
            feature.qualifiers.get("gene", [])
        )

        html_lines.append(
            f"<tr>"
            f"<td>{escape(record.id)}</td>"
            f"<td>{start}</td>"
            f"<td>{end}</td>"
            f"<td>{strand}</td>"
            f"<td>{escape(feature.type)}</td>"
            f"<td title='{escape(hover)}'>{escape(label)}</td>"
            f"</tr>"
        )

html_lines.append("</table></body></html>")

with open(f"{html_path}", "w") as f:
    f.write("\\n".join(html_lines))


PYCODE
    """
}
