process IDENTIFY_TI_RI_PLASMIDS{
    tag "${meta.id}"
    label 'process_medium'

    container "oras://community.wave.seqera.io/library/bbmap_git:2fc3d4caccbdba08"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(fna)
    path(db)

    output:
    tuple val(meta), path("${prefix}.oncogenic_plasmid_final.out")          , emit: plasmidfinal
    tuple val(meta), path("${prefix}.oncogenic_plasmid_type.sketch.out")    , emit: sketch
    tuple val(meta), path("${prefix}.oncogenic_plasmid_contigs.sketch.out") , emit: contigs
    tuple val(meta), path("${prefix}.oncogenic_plasmid_final.out.contiglist"), emit: contiglist
    val(true), emit: done

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """

    echo "Identifying oncogenic (Ti/Ri) plasmid contigs"

    INREF=`ls -1 ${db}/Weisberg2022PhilTransB/*.sketch ${db}/burr/*.sketch | tr '\\n' ','`

    comparesketch.sh in=${fna} ref=\$INREF minwkid=0.5 out=${prefix}.oncogenic_plasmid_type.sketch.out 1>&2
    comparesketch.sh in=${fna} ref=\$INREF mode=sequence minwkid=0.5 out=${prefix}.oncogenic_plasmid_contigs.sketch.out 1>&2

    echo "oncogenic plasmid type: (best hit reference prefix)" | tee ${prefix}.oncogenic_plasmid_final.out
    grep 'WKID' -A 1 ${prefix}.oncogenic_plasmid_type.sketch.out | tail -n 1 | sed 's/^.*\\s//g' | sed 's/__/\\t/g' | tee -a ${prefix}.oncogenic_plasmid_final.out || true
    echo "" | tee -a ${prefix}.oncogenic_plasmid_final.out
    echo "oncogenic plasmid-like contigs: (best hit reference prefix)" | tee -a ${prefix}.oncogenic_plasmid_final.out
    echo "contig	length	weighted_k-mer_id	k-mer_id	plasmid_type	best_hit_ref_prefix" | tee -a ${prefix}.oncogenic_plasmid_final.out

    grep 'WKID' ${prefix}.oncogenic_plasmid_contigs.sketch.out -A1 -B1 | grep -v '^--\$' | grep -v 'WKID' | sed 's/^Query: //g' | sed 's/\\s.*Bases: /___BASES___/g' | sed 's/\\s.*File:.*\$//g' | paste -d " " - - | sed 's/\\s\\+/\\t/g' | cut -f 1,2,3,14 | sed 's/___BASES___/\\t/g' | sed 's/__\\(\\S\\+\\)\$/\\t\\1/g' | tee -a ${prefix}.oncogenic_plasmid_final.out | tee ${prefix}.oncogenic_plasmid_final.out.contiglist || true

    """

}