# Import Python Packages
import numpy as np
import pandas as pd
import argparse
import logging

# Get Parameter
parser = argparse.ArgumentParser(description='SV background allele frequence annotation')
parser.add_argument('--spid', type=str, default = "None", metavar = "", help='Sample ID list if sv raw file without columns names. Default is None.')
parser.add_argument('--threshold', type=str, default = 0.8, metavar = "", help='SV overlap threshold to confirm same SV. Default is 0.8.')
parser.add_argument('--rare', type=str, default = "true", metavar = "", help='Deciding whether output rare sv part which background af < 0.01. Default is true.')
parser.add_argument('--background', type=str, required=True, metavar = "", help='Background cohort sv allele frequence set for annotation.')
parser.add_argument('--sv', type=str, required=True, metavar = "", help='SV file which contains CHROM, START, END, SVTYPE at least.')
parser.add_argument('-o', '--output_prefix', type=str, required=True, metavar = "", help='Output file prefix.')
args = parser.parse_args()

threshold = args.threshold
spid_list = args.spid
output_prefix = args.output_prefix

logger = logging.getLogger('SV-ANNO')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# Overlap function
def sv_overlap(sv_source, chrom, start, end, sv_type):
    if type(chrom) == int:
        chrom = str(chrom)
    # same chromosome and sv type
    if "chr" in chrom:
        tmp_bg_sv = sv_source[(sv_source.CHROM == chrom) & (sv_source.SVTYPE == sv_type)].copy()
    else:
        tmp_bg_sv = sv_source[(sv_source.CHROM == "chr" + chrom) & (sv_source.SVTYPE == sv_type)].copy()
    
    # choose background svs overlaped with appointed sv
    tmp_bg_sv["e1e2"] = tmp_bg_sv.END-end
    tmp_bg_sv["e1s2"] = tmp_bg_sv.END-start
    tmp_bg_sv["e2s1"] = end-tmp_bg_sv.START
    tmp_bg_sv["s2s1"] = start-tmp_bg_sv.START
    tmp_bg_sv_overlap = tmp_bg_sv[~(((tmp_bg_sv.e1e2 < 0) & (tmp_bg_sv.e1s2 < 0) & (tmp_bg_sv.s2s1 > 0) & (tmp_bg_sv.e2s1 > 0)) |
                                   ((tmp_bg_sv.e1e2 > 0) & (tmp_bg_sv.e1s2 > 0) & (tmp_bg_sv.s2s1 < 0) & (tmp_bg_sv.e2s1 < 0)))].copy()
    
    if len(tmp_bg_sv_overlap) == 0:
        return 0
    elif len(tmp_bg_sv_overlap) == 1:
        tmp_bg_sv_overlap["SVLEN"] = tmp_bg_sv_overlap.END - tmp_bg_sv_overlap.START
        tmp_bg_sv_overlap_percent = (min(tmp_bg_sv_overlap.END.values[0], end) - max(tmp_bg_sv_overlap.START.values[0], start)) / max(tmp_bg_sv_overlap.SVLEN.values[0], end - start)

        if tmp_bg_sv_overlap_percent > threshold:
            return tmp_bg_sv_overlap.Background_AF.tolist()[0]
        else:
            return 0
    else:
        tmp_bg_sv_overlap["SVLEN"] = tmp_bg_sv_overlap.END - tmp_bg_sv_overlap.START
        # choose the max overlap background sv
        tmp_bg_sv_overlap_percent = []
        
        for i in range(len(tmp_bg_sv_overlap)):
            tmp_bg_sv_overlap_len = (min([tmp_bg_sv_overlap.iloc[i,].END, end]) - max(tmp_bg_sv_overlap.iloc[i,].START, start))
            tmp_bg_sv_overlap_percent = tmp_bg_sv_overlap_percent + [tmp_bg_sv_overlap_len / max(tmp_bg_sv_overlap.iloc[i,].SVLEN, end - start)]
        tmp_bg_sv_overlap["overlap_percent"] = tmp_bg_sv_overlap_percent
        
        if len([x for x in tmp_bg_sv_overlap_percent if x >threshold]) > 0:
            return tmp_bg_sv_overlap.loc[tmp_bg_sv_overlap.overlap_percent > threshold].Background_AF.max()
        else:
            return 0


# Load data
if spid_list == "None":
    sv_file = pd.read_csv(args.sv,sep='\t', low_memory=False)
else:
    sample_list = pd.read_csv(spid_list,sep='\t', names=["spid"]).spid.tolist()
    sv_file = pd.read_csv(args.sv,sep='\t', names=["CHROM","START","END","SVTYPE"] + sample_list, low_memory=False)

background_sv_af = pd.read_csv(args.background,sep='\t')

# Background af annotation
sv_file["Background_AF"] = sv_file.apply(lambda row: sv_overlap(background_sv_af,row["CHROM"],row["START"],row["END"],row["SVTYPE"]), axis=1)
sv_file.to_csv(output_prefix + "_background_af_annotation.txt", sep="\t", index=False, header=True)

# Extract rare sv
if args.rare == "true":
    sv_rare = sv_file.loc[sv_file.Background_AF < 0.01].copy()
    sv_rare.to_csv(output_prefix + "_background_af_annotation_rare.txt", sep="\t", index=False, header=True)
