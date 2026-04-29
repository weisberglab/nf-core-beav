process GENERATE_GENBANK_FILES{
    tag "${meta.id}"
    label 'process_medium'

    container "oras://community.wave.seqera.io/library/bcbio-gff_biopython:43dd3826ecc32a41"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input), val(type)

    output:
    tuple val(meta), path("${meta.id}.gbk")   , emit: gbk, optional: true
    tuple val(meta), path("${meta.id}.gbff")  , emit: gbff, optional: true
    tuple val(meta), path("${meta.id}.fna")   , emit: fna, optional: true
    tuple val(meta), path("${meta.id}.faa")   , emit: faa, optional: true
    tuple val(meta), path("${meta.id}.gff3")  , emit: gff, optional: true
    tuple val(meta), path("${meta.id}.ffn")   , emit: ffn, optional: true
    tuple val(meta), path("${meta.id}_B.fna") , emit: fnaB, optional: true

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
python - <<'PYCODE'
import sys
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from BCBio import GFF
import shutil

infile = "${input}"
ext = os.path.splitext(infile)[1].lower() # Should only be ".gbff", ".gbk", ".fna", ".fa", or ".fasta"
print(ext)
bakta_state = "${type}"

if bakta_state == "NB":
    # Needs Bakta annotation

    if ext == ".gbff" or ext == ".gbk":
        # Convert to fna for Bakta
        SeqIO.convert(infile, "genbank", "${prefix}_B.fna", "fasta")
        os.remove(infile)
    else:
        # Ready for Bakta, just match output style for parsing
        shutil.copy(infile, f"${prefix}_B{ext}")
        os.remove(infile)

else:
    # Does not need Bakta annotation

    if ext == ".gbff" or ext == ".gbk":
        # Generate files from genbank
        gbk_filename = "${input}"

        fna_filename = "${prefix}.fna"
        faa_filename = "${prefix}.faa"
        gff3_filename = "${prefix}.gff3"
        ffn_filename = "${prefix}.ffn"

        gbk_input_handle = open(gbk_filename, "r")

        #gbk to fna
        fna_output_handle = open(fna_filename, "w")
        for seq_record in SeqIO.parse(gbk_input_handle, "genbank"):
        #       print ("Dealing with GenBank record %s" % seq_record.id)
                fna_output_handle.write(">%s %s\\n%s\\n" % (
                        seq_record.id,
                        seq_record.description,
                        seq_record.seq))
        fna_output_handle.close()
        gbk_input_handle.close()

        #gbk to faa
        gbk_input_handle = open(gbk_filename, "r")
        faa_output_handle = open(faa_filename, "w")
        for seq_record in SeqIO.parse(gbk_input_handle, "genbank") :
        #       print ("Dealing with GenBank record %s" % seq_record.id)
                for seq_feature in seq_record.features:
                        if seq_feature.type=="CDS" and 'translation' in seq_feature.qualifiers and 'pseudogene' not in seq_feature.qualifiers:
                                assert len(seq_feature.qualifiers['translation'])==1
                                faa_output_handle.write(">%s %s\\n%s\\n" % (
                                        seq_feature.qualifiers['locus_tag'][0],
                                        seq_feature.qualifiers['product'][0],
                                        seq_feature.qualifiers['translation'][0]))

        faa_output_handle.close()
        gbk_input_handle.close()

        #gbk to gff3
        gbk_input_handle = open(gbk_filename, "r")
        gff3_output_handle = open(gff3_filename, "w")
        GFF.write(SeqIO.parse(gbk_input_handle, "genbank"), gff3_output_handle)
        gff3_output_handle.close()
        gbk_input_handle.close()

        #gbk to ffn
        gbk_input_handle = open(gbk_filename, "r")
        ffn_input_handle = open(ffn_filename, 'w')
        for seq_record in SeqIO.parse(gbk_input_handle, "genbank") :
        #        print ("Dealing with GenBank record %s" % seq_record.id)
            for seq_feature in seq_record.features:
                if seq_feature.type=="CDS":
                    product = seq_feature.qualifiers.get("product")[0]
                    gene = seq_feature.qualifiers.get("gene")
                    if gene is not None:
                        gene = seq_feature.qualifiers.get("gene")[0]
                    else:
                        gene = ""
                    locus = seq_feature.qualifiers.get("locus_tag")
                    if locus is not None:
                        locus = seq_feature.qualifiers.get("locus_tag")[0]
                    else:
                        locus = ""
                    nucleotide = seq_feature.extract(seq_record).seq
                    locus = seq_feature.qualifiers.get("locus_tag")[0]
                    ffn_input_handle.write(">%s %s %s\\n%s\\n" % (
                                        locus,
                                        product,
                                        gene,
                                        nucleotide))
        ffn_input_handle.close()
        gbk_input_handle.close()
    else:
        print("Already annotated non genbank files cannot be used since they are missing required information, send to Bakta.")
        shutil.copy(infile, f"${prefix}_B.{ext}")
        os.remove(infile)
PYCODE
    """
}
