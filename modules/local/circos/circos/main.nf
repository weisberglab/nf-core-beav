process CIRCOS_CIRCOS {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/${meta.id}_results/circos_plots", mode: 'copy'

    container "oras://community.wave.seqera.io/library/biopython_matplotlib_numpy_pycirclize_python:e2d6ec959830e694"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path('*')

    output:
    path("${prefix}.circos.png"), emit: png
    path("${prefix}.circos.pdf"), emit: pdf
    path("${prefix}.oncogenes.pdf"), emit: oncpdf, optional: true
    path("${prefix}.oncogenes.png"), emit: oncpng, optional: true

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
python - <<'PYCODE'

# Dictionaries to define gene of interest in Ti Plasmid

# agrocinopine transport/catabolism
agrocinopine_dict = {
    'accA': 'agrocinopine ABC transporter, substrate binding protein accA',
    'accB': 'agrocinopine ABC transporter, nucleotide binding/ATPase accB',
    'accC': 'agrocinopine ABC transporter, nucleotide binding/ATPase accC',
    'accD': 'agrocinopine ABC transporter, membrane spanning protein accD',
    'accE': 'agrocinopine ABC transporter, membrane spanning protein accE',
    'accF_2': 'agrocinopine phosphodiesterase accF_2',
    'accF': 'agrocinopine phosphodiesterase accF',
    'accG_2': 'arabinose phosphate phosphatase accG_2',
    'accG': 'arabinose phosphate phosphatase accG',
    'accR': 'Repressor of the acc and arc operons accR'
}

# opine transport/catabolism
opine_tracat_dict = {
    'ocd': 'Ornithine cyclodeaminase',
    'nocQ': 'nopaline ABC transporter permease NocQ',
    'nocT': 'nopaline ABC transporter substrate-binding protein NocT',
    'nocP': 'nopaline ABC transporter ATP-binding protein NocP',
    'nocR': 'Regulatory protein NocR',
    'ooxA': 'Opine oxidase subunit A',
    'ooxB': 'Opine oxidase subunit B',
    'occT': 'Octopine-binding periplasmic protein',
    'occP': 'Octopine permease ATP-binding protein P',
    'occM': 'Octopine transport system permease protein OccM',
    'occQ': 'Octopine transport system permease protein OccQ',
    'occR': 'Octopine catabolism/uptake operon regulatory protein OccR'
}

# opine synthase genes
opine_synth_dict = {
    'acs': 'agrocinopine synthase',
    'ags': 'agropine synthase',
    'chsA': 'chrysopine synthase',
    'cus': 'cucumopine synthase',
    'mas1': 'mannopine synthase 1',
    'mas2': 'mannopine synthase 2',
    'nos': 'nopaline synthase',
    'ocs': 'octopine synthase',
    'sus': 'succinamopine synthase',
    'vis': 'vitopine synthase',
    'unk': 'unknown opine synthase'
}

# virulence genes
vir_dict = {
    'virA': 'vir locus protein VirA',
    'virB10': 'vir locus Type IV secretion system VirB10',
    'virB11': 'vir locus Type IV secretion system VirB11',
    'virB1': 'vir locus Type IV secretion system VirB1',
    'virB2': 'vir locus Type IV secretion system VirB2',
    'virB3': 'vir locus Type IV secretion system VirB3',
    'virB4': 'vir locus Type IV secretion system VirB4',
    'virB5': 'vir locus Type IV secretion system VirB5',
    'virB6': 'vir locus Type IV secretion system VirB6',
    'virB7': 'vir locus Type IV secretion system VirB7',
    'virB8': 'vir locus Type IV secretion system VirB8',
    'virB9': 'vir locus Type IV secretion system VirB9',
    'virC1': 'vir locus protein VirC1',
    'virC2': 'vir locus protein VirC2',
    'virD1': 'vir locus protein VirD1',
    'virD2': 'vir locus protein VirD2',
    'virD3': 'vir locus protein VirD3',
    'virD4': 'vir locus protein VirD4',
    'virD5': 'vir locus protein VirD5',
    'virE1': 'vir locus protein VirE1',
    'virE2': 'vir locus protein VirE2',
    'virE3': 'vir locus protein VirE3',
    'virF': 'exported virulence protein virF',
    'virG': 'vir locus protein VirG',
    'virH1': 'P-450 monooxygenase virH1',
    'virH2': 'P-450 monoxygenase virH2',
    'virK': 'virA/G regulated gene virK',
    'GALLS': 'vir protein GALLS'
}

# plasmid conjugation tra genes
tra_dict = {
    'traA': 'conjugal transfer protein TraA',
    'traB': 'conjugal transfer protein TraB',
    'traC': 'conjugal transfer protein TraC',
    'traD': 'conjugal transfer protein TraD',
    'traF': 'conjugal transfer protein TraF',
    'traG': 'conjugal transfer protein TraG',
    'traH': 'conjugal transfer protein TraH',
    'traI': 'conjugal transfer protein TraI',
    'traM': 'conjugal transfer protein TraM',
    'traR': 'conjugal transfer protein TraR'
}

# plasmid conjugation trb genes
trb_dict = {
    'trbB': 'conjugal transfer protein TrbB',
    'trbC': 'conjugal transfer protein TrbC',
    'trbD': 'conjugal transfer protein TrbD',
    'trbE': 'conjugal transfer protein TrbE',
    'trbF': 'conjugal transfer protein TrbF',
    'trbG': 'conjugal transfer protein TrbG',
    'trbH': 'conjugal transfer protein TrbH',
    'trbI': 'conjugal transfer protein TrbI',
    'trbJ': 'conjugal transfer protein TrbJ',
    'trbK': 'conjugal transfer protein TrbK',
    'trbL': 'conjugal transfer protein TrbL'
}

# replication gene
rep_dict = {
    'repA': 'plasmid partitioning protein RepA',
    'repB': 'plasmid partitioning protein RepB',
    'repC': 'plasmid replication protein RepC'
}

# T-DNA/oncogene genes
oncogene_dict = {
    'iaaH': 'indoleacetamide hydrolase',
    'iaaM': 'Tryptophan 2-monooxygenase',
    'ipt': 'Adenylate dimethylallyltransferase',
    'tzs': 'Adenylate dimethylallyltransferase'
}

oncogene_list = [
    'Tryptophan 2-monooxygenase',
    'Adenylate dimethylallyltransferase',
    '6a protein',
    'Uncharacterized protein 6',
    '6b protein',
    "Gene 4' protein",
    'C protein',
    'D protein',
    'E protein',
    '5 protein',
    'RolB-RolC domain-containing protein',
    'T-DNA oncoprotein'
]

import warnings
import argparse
import sys, os

from pycirclize import Circos
from pycirclize.parser import Genbank

import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore", message=".*Fontconfig.*")

def all_contig_circos(gbk_file, onco_label):
    "INPUT: Genbank file containing beav annotation\\nOUTPUT: Circos plot containing feature distribution in differnt contigs"
    # Load GBK file
    gbk = Genbank(gbk_file)

    # Get contig genome seqid & size, features dict
    seqid2size = gbk.get_seqid2size()
    seqid2features = gbk.get_seqid2features(feature_type=None)


    circos = Circos(seqid2size, space=min(1.0, 150.0/len(seqid2size)), start=15, end=345)

    intPresent = 0 # status will change if integron is present
    plasmidPresent = 0 # status will change if plasmid is present

    # Loop through the contigs
    contig_i = 0
    for sector in circos.sectors:
        #  highlight plasmids
        is_plasmid = False
        for feat in seqid2features[sector.name]:
            if feat.type == "plasmid":
                is_plasmid = True
                plasmidPresent = 1
                break
        # plot sector labels
        kwargs = dict(color="navy")
        sector.text(f"{sector.name}", orientation='vertical', r=110, size=10, **kwargs)
        # Plot outer track
        outer_track = sector.add_track((98, 100))
        if is_plasmid:
            outer_track.axis(fc="gold")
        else:
            outer_track.axis(fc="lightgrey")
        major_interval = 1000000
        minor_interval = int(major_interval / 100)
        if sector.size > minor_interval:
            outer_track.xticks_by_interval(major_interval, label_formatter=lambda v: f"{v / 1000000:.0f} Mb")
            outer_track.xticks_by_interval(minor_interval, tick_length=1, show_label=False)
    
        f_cds_track = sector.add_track((90, 95), r_pad_ratio=0.1)
        r_cds_track = sector.add_track((85, 90), r_pad_ratio=0.1)
        beav_track = sector.add_track((80, 85), r_pad_ratio=0.1)
        extra_feature =  sector.add_track((72, 77))
        gc_content_track = sector.add_track((67, 72))
        gc_skew_track = sector.add_track((62, 67))
        
        # Plot forward/reverse CDS, ICE, origin, beav tracks
        for feature in seqid2features[sector.name]:
            if feature.type == "CDS":
                if feature.location.strand == 1:
                    f_cds_track.genomic_features([feature], fc="tomato")
                elif feature.location.strand == -1:
                    r_cds_track.genomic_features([feature], fc="skyblue")
            elif feature.type == "mobile_element":
                    if 'mobile_element_type' in feature.qualifiers.keys():
                        if feature.qualifiers['mobile_element_type'][0] == 'phage':
                            beav_track.genomic_features([feature], fc='royalblue') # Phage
                        elif feature.qualifiers['mobile_element_type'][0] == 'integrative element':
                            beav_track.genomic_features([feature], fc="turquoise") # Non phage mobile element
            elif feature.type == "rRNA":
                fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                gc_skew_track.xticks([(fx1 + fx2)/2], outer=True, label_size=6, labels=[''], label_orientation="vertical", line_kws={'ec':'yellowgreen'}) #rRNA
            elif feature.type == "tRNA":
                fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                gc_skew_track.xticks([(fx1 + fx2)/2], outer=True, label_size=6, labels=[''], label_orientation="vertical", line_kws={'ec':'orange'}) #tRNA

                
            if feature.type == "CDS":
                
                if "note" in feature.qualifiers.keys(): 
                    if 'MacSy' in feature.qualifiers['note'][0]:
                        sec_sys = feature.qualifiers['note'][0].split(': ')[1]
                        beav_track.genomic_features([feature], fc="darkorange") # Secretion Systems
                        if 'T4SS' in sec_sys:
                            fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                            beav_track.xticks([(fx1 + fx2)/2], outer=False, line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # T4SS
                        elif 'T3SS' in sec_sys:
                            fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                            beav_track.xticks([(fx1 + fx2)/2], outer=False, line_kws={'ec':'plum'}, text_kws={'color':'darkred'}) # T3SS
                        elif 'T6SS' in sec_sys:
                            fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                            beav_track.xticks([(fx1 + fx2)/2], outer=False, line_kws={'ec':'slateblue'}, text_kws={'color':'darkred'}) # T6SS
                    elif 'SMASH' in feature.qualifiers['note'][0]:
                        beav_track.genomic_features([feature], fc="forestgreen") # Secondary metabolite
                    elif 'IntegronFinder' in feature.qualifiers['note'][0]:
                        beav_track.genomic_features([feature], fc="black") # Integron
                        intPresent = 1
                    
                    
                if 'gene' in feature.qualifiers.keys():
                    if feature.qualifiers['gene'][0] == 'vgrG':
                        beav_track.genomic_features([feature], fc="darkorange") # Secretion Systems 
                        fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                        beav_track.xticks([(fx1 + fx2)/2], outer=False, line_kws={'ec':'slateblue'}, text_kws={'color':'darkred'}) # T6SS
            
            if feature.type == "rep_origin":
                fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                gc_skew_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['Origin'], label_orientation="vertical", line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # Origin of replication

            if feature.type == "misc_feature" and feature.qualifiers['note'][0] == 'dif':
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end    ))
                    gc_skew_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['dif'], label_orientation="vertical", line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # Dif site
            
            if feature.type == 'oriT':
                fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                gc_skew_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['oriT'], label_orientation="vertical") # Origin of transfer

            if feature.type == "CDS":
                if 'gene' in feature.qualifiers.keys():
                    if feature.qualifiers['gene'][0] in ['repA', 'repB', 'repC']:
                        extra_feature.genomic_features([feature], fc='black',  plotstyle="arrow")
                                            
        # Plot GC skew
        pos_list, gc_skews = gbk.calc_gc_skew(seq=str(gbk.records[contig_i].seq))
        positive_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
        negative_gc_skews = np.where(gc_skews < 0, gc_skews, 0)
        abs_max_gc_skew = np.max(np.abs(gc_skews))
        vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew
        gc_skew_track.fill_between(
            pos_list, positive_gc_skews, 0, vmin=vmin, vmax=vmax, color="olive"
        )
        gc_skew_track.fill_between(
            pos_list, negative_gc_skews, 0, vmin=vmin, vmax=vmax, color="purple"
        )


        # update contig index
        contig_i += 1


    text_common_kws = dict(ha="center", va="center", size=6)
    circos.text("CDS +strand", r=93, color="dimgrey", **text_common_kws)
    circos.text("CDS -strand", r=87, color="dimgrey", **text_common_kws)
    circos.text("Beav", r=82, color="dimgrey", **text_common_kws)
    circos.text("repABC", r=74, color="dimgrey", **text_common_kws)
    circos.text("RNA", r=68, color="dimgrey", **text_common_kws)
    circos.text("GC skew", r=63, color="dimgrey", **text_common_kws)

    # Prepare legends for circos plot
    legend_list = [
            Patch(color="darkorange", label="Secretion Systems"),
            Line2D([], [], color="plum", label="T3SS", marker='_', ls='None'),
            Line2D([], [], color="darkred", label="T4SS", marker='_', ls='None'),
            Line2D([], [], color="slateblue", label="T6SS", marker='_', ls='None'),
            Patch(color="royalblue", label="Prophages"),
            Patch(color="turquoise", label="ICEs"),
            Patch(color="forestgreen", label="Metabolite Clusters"),
            Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=6, ls="None"),
            Line2D([], [], color="purple", label="Negative GC Skew", marker="v", ms=6, ls="None"),
            Line2D([], [], color="yellowgreen", label="rRNA", marker="_", ms=6, ls="None"),
            Line2D([], [], color="orange", label="tRNA", marker="_", ms=6, ls="None")
        ]

    if intPresent == 1:
            legend_list.insert(6, Patch(color="black", label="Integron"))
    if plasmidPresent == 1:
            legend_list.insert(1, Patch(color="gold", label="Plasmid Region"))

    fig = circos.plotfig()
    _ = circos.ax.legend(
        title = f"{onco_label}\\n",
        handles = legend_list,
        bbox_to_anchor=(0.5, 0.45),
        loc="center",
        fontsize=6
    )

    fig.savefig(f'{onco_label}.circos.png', dpi=300)
    fig.savefig(f'{onco_label}.circos.pdf', dpi=300)
    print(f'Image saved: {onco_label}.circos.png/.pdf')

    pass


def splice_genbank_contig(gbk_file, contig_id, onco_label):
    "INPUT:\\n\\tgbk_file: Genbank file containing beav annotation and has multiple contigs\\n\\tcontig_id: List of contigs\\nOUTPUT: A spliced genbank file that has only the contig of interest"
    with open(gbk_file, 'r') as f:
        gbk_raw = f.read()
    
    gbk_contigs = gbk_raw.strip().split('//\\n')
    
    contig_of_interest = []
    for i in gbk_contigs:
        for j in contig_id:
            if j in i.split('\\n', 1)[0].split():
                contig_of_interest.append(i)

    with open(f'{onco_label}.oncogenic.gbk', 'w') as f:
        for i in range(len(contig_of_interest)):
            f.writelines(contig_of_interest[i])
            if i < len(contig_of_interest):
                f.writelines('//\\n')
            
    

    return f'{onco_label}.oncogenic.gbk'


def oncogenic_circos(gbk_file, onco_label):
    "INPUT:\\ngbk_file: Single contig genbank file containing beav annotation\\nOUTPUT: Circos plot containing oncogenic plasmid feature distribution of that contig"

    #Load GBK file
    gbk = Genbank(gbk_file)

    #Get features
    seqid2size = gbk.get_seqid2size()
    seqid2features = gbk.get_seqid2features(feature_type=None)
    
    #Define sector for single contig
    circos = Circos(seqid2size, space=3)
    #main_kwargs = dict(ha="center", va="top") 
    #circos.text(f"{get_base_file_name(gbk_file)}\\nTi/Ri plasmid", r=5, size=8, **main_kwargs)

    # Extract CDS gene labels
    # List containing CDS product labels if gene names are not present
    product_of_interest = []

    for sector in circos.sectors:
        # plot sector labels
        kwargs = dict(color="navy")
        sector.text(f"{sector.name}", orientation='vertical', r=110, size=10, **kwargs)

        # add tracks
        cds_track = sector.add_track((95, 100))
        scale_track = sector.add_track((95, 95))
        cds_track.axis(fc="#EEEEEE", ec="none") 

        # Get sector specific labels
        ##print(sector, seqid2features[sector.name])

        sector_pos_list, sector_labels = [], []
        for feat in seqid2features[sector.name]:
            start, end = int(str(feat.location.end)), int(str(feat.location.start))
            pos = (start + end) / 2
            label = feat.qualifiers.get("gene", [""])[0]
            product = feat.qualifiers.get("product", [""])[0]
            if label == "":
                if product in product_of_interest:
                    label = product
            if product == "" or product.startswith("hypothetical"):
                continue
            if len(label) > 20:
                label = label[:20] + "..."
            sector_pos_list.append(pos)
            sector_labels.append(label)

        # plot oncogene specific features
        for feature in seqid2features[sector.name]:
            if feature.type == "CDS":
                cds_track.genomic_features([feature], fc="lightgrey", plotstyle="arrow")

            if 'gene' in feature.qualifiers.keys():
                if feature.qualifiers['gene'][0] in vir_dict.keys():
                    cds_track.genomic_features([feature], fc="olive", plotstyle="arrow") #T-DNA transfer
                elif feature.qualifiers['gene'][0] in tra_dict.keys():
                    cds_track.genomic_features([feature], fc="indigo", plotstyle="arrow") #tra genes
                elif feature.qualifiers['gene'][0] in trb_dict.keys():
                    cds_track.genomic_features([feature], fc="purple", plotstyle="arrow") #trb genes
                elif feature.qualifiers['gene'][0] in oncogene_dict.keys():
                    cds_track.genomic_features([feature], fc="orange", plotstyle="arrow") #T-DNA/Oncogene            
                elif feature.qualifiers['gene'][0] in rep_dict.keys():
                    cds_track.genomic_features([feature], fc="salmon", plotstyle="arrow") #repABC
                elif feature.qualifiers['gene'][0] in opine_synth_dict.keys():
                    cds_track.genomic_features([feature], fc="darkgreen", plotstyle="arrow") #Opine synthase genes
                elif feature.qualifiers['gene'][0] in opine_tracat_dict.keys():
                    cds_track.genomic_features([feature], fc="turquoise", plotstyle="arrow") #Opine transport/catabolism genes
                elif feature.qualifiers['gene'][0] in agrocinopine_dict.keys():
                    cds_track.genomic_features([feature], fc="steelblue", plotstyle="arrow") #Agrocinopine transport/catabolism genes

            if feature.type == 'oriT':
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['oriT'], label_orientation="vertical") # Origin of transfer
            
            if feature.type == "rep_origin":
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['Origin'], label_orientation="vertical", line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # Origin of replication
                
            if 'product' in feature.qualifiers.keys():
                if feature.qualifiers['product'][0] in oncogene_list:
                    product_of_interest.append(feature.qualifiers['product'][0]) 
                    cds_track.genomic_features([feature], fc="orange", plotstyle="arrow") #T-DNA/Oncogene
                elif feature.qualifiers['product'][0] in rep_dict.values():
                    product_of_interest.append(feature.qualifiers['product'][0])
                    cds_track.genomic_features([feature], fc="salmon", plotstyle="arrow") #repABC
                elif feature.qualifiers['product'][0] in opine_synth_dict.values():
                    product_of_interest.append(feature.qualifiers['product'][0])
                    cds_track.genomic_features([feature], fc="darkgreen", plotstyle="arrow") #Opine synthase genes
                elif feature.qualifiers['product'][0] in opine_tracat_dict.values():
                    product_of_interest.append(feature.qualifiers['product'][0])
                    cds_track.genomic_features([feature], fc="turquoise", plotstyle="arrow") #Opine transport/catabolism genes
                elif feature.qualifiers['product'][0] in agrocinopine_dict.values():
                    product_of_interest.append(feature.qualifiers['product'][0])
                    cds_track.genomic_features([feature], fc="steelblue", plotstyle="arrow") #Agrocinopine transport/catabolism genes
                elif 'transposase' in feature.qualifiers['product'][0].lower():
                    cds_track.genomic_features([feature], fc="gray", plotstyle="arrow") # Transposase
                elif 'integrase' in feature.qualifiers['product'][0].lower():
                    cds_track.genomic_features([feature], fc="gray", plotstyle="arrow") # Integrase
                elif 'recombinase' in feature.qualifiers['product'][0].lower():
                    cds_track.genomic_features([feature], fc="gray", plotstyle="arrow") # Recombinase
                elif 'ABC transporter' in feature.qualifiers['product'][0]:
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, line_kws={'ec':'darkgreen'}) #ABC Transporter
            
            if 'note' in feature.qualifiers.keys():
                if 'origin_of_replication' in ' '.join(feature.qualifiers['note']):
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['origin'], label_orientation="vertical") # Origin of replication
                elif 'virbox' in ' '.join(feature.qualifiers['note']):
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['virbox'], label_orientation="vertical", line_kws={'ec':'navy'}, text_kws={'color':'navy'}) # virbox
                elif 'trabox' in ' '.join(feature.qualifiers['note']):
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['trabox'], label_orientation="vertical", line_kws={'ec':'navy'}, text_kws={'color':'navy'}) # trabox
                elif 'T-DNA_right_border' in ' '.join(feature.qualifiers['note']):
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['T-DNA Right'], label_orientation="vertical", line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # t-dna right border
                elif 'T-DNA_left_border' in feature.qualifiers['note']:
                    fx1, fx2 = int(str(feature.location.parts[0].start)), int(str(feature.location.parts[-1].end))
                    cds_track.xticks([(fx1 + fx2)/2], outer=False, label_size=6, labels=['T-DNA Left'], label_orientation="vertical", line_kws={'ec':'darkred'}, text_kws={'color':'darkred'}) # t-dna left border
                elif 'integrase' in ' '.join(feature.qualifiers['note']).lower():
                    cds_track.genomic_features([feature], fc="gray", plotstyle="arrow") # Integrase

        # Plot CDS product labels on outer position
        cds_track.xticks(
            sector_pos_list,
            sector_labels,
            label_orientation="vertical",
            show_bottom_line=True,
            outer=True,
            label_size=6,
            line_kws=dict(ec="grey"),
        )
        # Plot xticks & intervals on inner position
        scale_track.xticks_by_interval(
            interval=25000,
            label_size=6,
            outer=False,
            show_bottom_line=True,
            label_formatter=lambda v: f"{v/ 1000:.1f} Kb",
            label_orientation="vertical",
            line_kws=dict(ec="grey"),
        )


    fig = circos.plotfig()
    _ = circos.ax.legend(
        title = f"{onco_label}\\n",
        handles=[
            Line2D([], [], color="lightgrey", label="CDS +strand", marker=">", ms=6, ls="None"),
            Line2D([], [], color="lightgrey", label="CDS -strand", marker="<", ms=6, ls="None"),
            Patch(color="olive", label=r"\$\\it{vir}\$ genes"),
            Patch(color="orange", label="T-DNA/Oncogenes"),
            Line2D([], [], color="darkgreen", label="ABC Transporter", marker="_", ms=6, ls="None"),
            Patch(color="indigo", label=r"\$\\it{tra}\$"),
            Patch(color="purple", label=r"\$\\it{trb}\$"),
            Patch(color="salmon", label=r"\$\\it{repABC}\$"),
            Patch(color="darkgreen", label="Opine synthase "),
            Patch(color="turquoise", label="Opine transport"),
            Patch(color="steelblue", label="Agrocinopine transport"),
            Patch(color="gray", label="Transposase/Recombinase"), 
        ],
        bbox_to_anchor=(0.5, 0.45),
        loc="center",
        fontsize=6
    )
    if '\\n' in onco_label:
        file_name = onco_label.split('\\n')[0]
    else:
        file_name = onco_label
    fig.savefig(f'{file_name}.oncogenes.png', dpi=300)
    fig.savefig(f'{file_name}.oncogenes.pdf', dpi=300)
    print(f'Image saved: {file_name}.oncogenes.png/pdf')

    pass

# Parse arguments from Nextflow

parser = argparse.ArgumentParser()
parser.add_argument('--input', '-i', type=str, required=True, help='Input GenBank file')
parser.add_argument('--contigs', '-c', type=str, nargs='*', required=False, help='Plot selected contig(s)')
parser.add_argument('--pTi', type=str, nargs='*', required=False, help='Plot oncogenic contig(s)')
parser.add_argument('--plasmid', '-p', type=str, required=False, help='Oncogenic plasmid type')
parser.add_argument('--label', '-l', type=str, required=False, help='Label')

contig_path = f"${prefix}.oncogenic_plasmid_final.out.contiglist"

contigs = []

if os.path.exists(contig_path):
    with open(contig_path, 'r') as f:
        contigs = [line.strip().split('\\t')[0] for line in f if line.strip()]

args = parser.parse_args([
    '--input', '${prefix}_final.gbk', *(['--contigs', *contigs] if contigs else [])
])

strain_id = "${prefix}"

# Run genome-wide circos
all_contig_circos(args.input, strain_id)

# If agrobacterium run, make plasmid-specific plots
if contigs:
    for contig_id in contigs:
        plasmid_gbk = splice_genbank_contig(args.input, [contig_id], strain_id)
        oncogenic_circos(plasmid_gbk, f"{strain_id}\\n{contig_id}")

PYCODE
    """
}
