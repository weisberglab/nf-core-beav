process ANTISMASH_GENBANKS_TO_TABLE{
    tag "${meta.id}"
    label 'process_medium'

    container "oras://community.wave.seqera.io/library/bcbio-gff_pandas_pip_biopython_ordereddict:f4b3ca75388bef64"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(gbk_results_dir)

    output:
    tuple val(meta), path("${meta.id}_antismash_final.tsv"), emit: asfinal
    tuple val(meta), path("${meta.id}_antismash_final.tsv.beav.subset"), emit: subset
    tuple val(meta), path("${meta.id}_synopsis.tsv"), emit: synopsis
    val(true)                                               , emit: done

    script:
    """
python - <<'PYCODE'
#downloaded from https://raw.githubusercontent.com/jolespin/veba/main/src/scripts/antismash_genbanks_to_table.py
import sys, os, argparse, gzip, glob
import pandas as pd
from Bio import SeqIO

sys.argv = ["script.py",
    "-i", "${gbk_results_dir}",
    "-o", "${meta.id}_antismash_final.tsv",
    "-s", "${meta.id}_synopsis.tsv"
]

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2023.01.08"

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = "Running: {} v{} via Python v{} | {}".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -i ${gbk_results_dir} -o ${meta.id}_antismash_final.tsv".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--antismash_directory", type=str, help = "path/to/antismash_directory")
    parser.add_argument("-o","--output", type=str, default="stdout", help = "path/to/output.tsv [Default: stdout]")
    parser.add_argument("-e","--exclude_contig_edges", action="store_true", help = "Exclude clusters that are on contig edges")
    parser.add_argument("-s","--synopsis", type=str, help = "Summary file output")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output
    if opts.output == "stdout":
        opts.output = sys.stdout 

    # Wildcard genbank files
    output = list()
    for fp in glob.glob(os.path.join(opts.antismash_directory, "*region*.gbk")):
        # Get genome and region identifiers
        id_genome = fp.split("/")[-2]
        id_region = fp.split(".")[-2]
        # Iterate through sequence records
        for seq_record in SeqIO.parse(fp, "genbank"):
            # Iterate through sequence features
            id_contig = seq_record.id
            cds_features = list()
            product = None
            contig_edge = None
            for feature in seq_record.features:
                if feature.type == "CDS":
                    data = {"genome_id":id_genome, "region_id":id_region, "start":int(feature.location.start), "end":int(feature.location.end), "strand":feature.location.strand}
                    for k,v in feature.qualifiers.items():
                        if isinstance(v,list):
                            if len(v) == 1:
                                v = v[0]
                        data[k] = v
                    cds_features.append(pd.Series(data))
                elif feature.type == "region":                    
                    product = feature.qualifiers["product"][0]
                    contig_edge = feature.qualifiers["contig_edge"][0]
            df = pd.DataFrame(cds_features)
            
            # NCBI Genbank
            # Add gene id preferentially if one doesn't exist: gene_id, Name, ID, and locus_tag in that order.
            if "gene_id" not in df.columns:
                genes = None
                if "Name" in df.columns:
                    genes = df["Name"]
                elif "ID" in df.columns:
                    genes = df["ID"]
                elif "locus_tag" in df.columns:
                    genes = df["locus_tag"]
                df["gene_id"] = genes
            # If no contig identifiers, add it.  This prioritizes gff fields.
            if "contig_id" not in df.columns:
                df["contig_id"] = id_contig

            
            df.insert(0, "bgc_type", product)
            df.insert(1, "cluster_on_contig_edge", contig_edge)
            output.append(df)
    df_bgcs = pd.concat(output, axis=0).sort_values(["genome_id", "contig_id", "start", "end"])
    # for field in ["start","end"]:
    #     df_bgcs[field] = df_bgcs[field].astype(int)
    df_bgcs = df_bgcs.set_index(["genome_id", "contig_id", "region_id", "gene_id"]).sort_index()
    if "translation" in df_bgcs.columns:
        df_bgcs["translation"] = df_bgcs.pop("translation")

    if opts.synopsis:
        df_tmp = df_bgcs.reset_index()

        # Ensure required columns exist
        for col in ["genome_id", "bgc_type", "cluster_on_contig_edge"]:
            if col not in df_tmp.columns:
                df_tmp[col] = None

        # Count all BGCs
        df_all = (
            df_tmp.groupby(["genome_id", "bgc_type"])
            .size()
            .reset_index(name="number_of_bgcs")
            .sort_values("number_of_bgcs", ascending=False)
        )

        # Count excluding edges
        mask = df_tmp["cluster_on_contig_edge"].map(lambda x: str(x).lower() == "false")
        df_noedges = (
            df_tmp[mask]
            .groupby(["genome_id", "bgc_type"])
            .size()
            .reset_index(name="number_of_bgcs(not_on_edge)")
            .sort_values("number_of_bgcs(not_on_edge)", ascending=False)
        )

        df_synopsis = pd.merge(
            df_all, df_noedges,
            on=["genome_id", "bgc_type"],
            how="outer"
        ).fillna(0)

        df_synopsis.to_csv(opts.synopsis, sep="\t", index=False)

    # Exclude contig edges
    if opts.exclude_contig_edges:
        df_bgcs = df_bgcs.loc[~df_bgcs["cluster_on_contig_edge"].map(eval).values]

    # Output
    df_bgcs.to_csv(opts.output, sep="\\t")

    tabdf = df_bgcs.reset_index() 
    if set(['NRPS_PKS','sec_met_domain']).issubset(tabdf.columns):
        newdf = tabdf[ ['locus_tag','contig_id','region_id','bgc_type','gene_functions','NRPS_PKS','sec_met_domain'] ].copy()
    elif 'sec_met_domain' in tabdf.columns:
        newdf = tabdf[ ['locus_tag','contig_id','region_id','bgc_type','gene_functions','sec_met_domain'] ].copy()
        newdf['NRPS_PKS'] = ''
    else:
        newdf = tabdf[ ['locus_tag','contig_id','region_id','bgc_type','gene_functions'] ].copy()
        newdf['NRPS_PKS'] = ''
        newdf['sec_met_domain'] = ''
    newdf['cluster'] = newdf['contig_id'] + "_" + newdf['region_id'] + ":" + newdf['bgc_type']
    newdf = newdf[ ['locus_tag','cluster','gene_functions','NRPS_PKS','sec_met_domain' ] ]
    newdf.to_csv(opts.output + ".beav.subset", sep="\\t", index=False)

if __name__ == "__main__":
    main()
    
PYCODE
    """
}