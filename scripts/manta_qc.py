# working environment is sv-anno

import numpy as np
import pandas as pd
import argparse
#import logging

# Load data
# Get Parameter
parser = argparse.ArgumentParser(description='Manta SV result filter by Lumpy, CNVnator and paragraph result')
parser.add_argument('-m', '--manta', type=str, required=True, metavar = "", help='Manta SV VCF file without vcf head.')
parser.add_argument('-l', '--lumpy', type=str, required=True, metavar = "", help='Lumpy SV bed file for DEL and DUP which contains CHROM, START, END, SVTYPE.')
parser.add_argument('-c', '--cnvnator', type=str, required=True, metavar = "", help='CNVnator SV bed file which contains CHROM, START, END, SVTYPE.')
parser.add_argument('-p', '--paragraph', type=str, required=True, metavar = "", help='Paragraph regenotype file for SVLEN < 200bp in Manta result.')
parser.add_argument('-t', '--threshold', type=str, default = 0.8, metavar = "", help='SV overlap threshold to confirm same SV. Default = 0.8.')
parser.add_argument('-o', '--output', type=str, required=True, metavar = "", help='Output file prefix.')
parser.add_argument('-f', '--filter', type=str, default = "true", metavar = "", help='Output file having removed low quality SV and TRA. Default is true.')
args = parser.parse_args()

# Load Files
threshold = args.threshold
output_prefix = args.output
manta_sv_raw = pd.read_table(args.manta, names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER_RAW","INFO_RAW","FORMAT","GT_RAW"])
lumpy_sv_info = pd.read_table(args.lumpy, names=["CHROM","START","END","SVTYPE"])
cnvnator_sv_info = pd.read_table(args.cnvnator, names=["CHROM","START","END","SVTYPE"])
manta_paragraph_sv = pd.read_table(args.paragraph, names=["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT_RAW"])
manta_paragraph_sv[["GT_NEW", "Others"]] = manta_paragraph_sv["GT_RAW"].str.split(":", n=1, expand =True)
manta_paragraph_sv["SVTYPE"] = manta_paragraph_sv["INFO"].str.split(";",expand = True).iloc[:,3].str.split("=",expand = True).iloc[:,1]
manta_paragraph_sv_gt = manta_paragraph_sv[["ID","GT_NEW"]].copy()

del manta_paragraph_sv

# Overlap function
def sv_overlap(sv_source, chrom, start, end, sv_type):
    # same chromosome and sv type
    tmp_bg_sv = sv_source[(sv_source.CHROM == chrom) & (sv_source.SVTYPE == sv_type)].copy()
    
    # choose background svs overlaped with appointed sv
    tmp_bg_sv["e1e2"] = tmp_bg_sv.END-end
    tmp_bg_sv["e1s2"] = tmp_bg_sv.END-start
    tmp_bg_sv["e2s1"] = end-tmp_bg_sv.START
    tmp_bg_sv["s2s1"] = start-tmp_bg_sv.START
    tmp_bg_sv_overlap = tmp_bg_sv[~(((tmp_bg_sv.e1e2 < 0) & (tmp_bg_sv.e1s2 < 0) & (tmp_bg_sv.s2s1 > 0) & (tmp_bg_sv.e2s1 > 0)) |
                                   ((tmp_bg_sv.e1e2 > 0) & (tmp_bg_sv.e1s2 > 0) & (tmp_bg_sv.s2s1 < 0) & (tmp_bg_sv.e2s1 < 0)))].copy()
    tmp_bg_sv_overlap["SVLEN"] = tmp_bg_sv_overlap.END - tmp_bg_sv_overlap.START
    
    if len(tmp_bg_sv_overlap) == 0:
        return 0
    elif len(tmp_bg_sv_overlap) == 1:
        tmp_bg_sv_overlap_percent = (min(tmp_bg_sv_overlap.END.tolist() + [end]) - max(tmp_bg_sv_overlap.START.tolist() + [start])) / max(tmp_bg_sv_overlap.SVLEN.tolist() + [end - start])
        if tmp_bg_sv_overlap_percent > threshold:
            return 1
        else:
            return 0
    else:
        # the appointed sv was recognized as multiple svs in other svtools
        tmp_bg_sv_overlap_all = 0
        tmp_bg_sv_overlap_percent = []
        
        for i in range(len(tmp_bg_sv_overlap)):
            tmp_bg_sv_overlap_len = (min([tmp_bg_sv_overlap.iloc[i,].END, end]) - max([tmp_bg_sv_overlap.iloc[i,].START, start]))
            tmp_bg_sv_overlap_percent = tmp_bg_sv_overlap_percent + [tmp_bg_sv_overlap_len/max(tmp_bg_sv_overlap.iloc[i,].SVLEN + [end - start])]
            tmp_bg_sv_overlap_all = tmp_bg_sv_overlap_all + tmp_bg_sv_overlap_len
            
        if len([x for x in tmp_bg_sv_overlap_percent if x > threshold]) > 0:
            return 1
        else:
            tmp_bg_sv_overlap_percent_sum = tmp_bg_sv_overlap_all / max(tmp_bg_sv_overlap.SVLEN.tolist() + [end - start])
            if tmp_bg_sv_overlap_percent_sum > threshold:
                return 1
            else:
                return 0




# Splite Manta SV by different classes: SVTYPE and SVLEN
manta_sv_raw["SVTYPE"] = manta_sv_raw["INFO_RAW"].str.split(";",expand = True).iloc[:,1].str.split("=",expand = True).iloc[:,1]

manta_sv_filter_type = manta_sv_raw[(manta_sv_raw.SVTYPE.isin(["DEL","DUP","INS","INV"]))].copy()
manta_sv_filter_type["END"] = manta_sv_filter_type["INFO_RAW"].str.split(";",expand = True).iloc[:,0].str.split("=",expand = True).iloc[:,1].astype(int)
manta_sv_filter_type[["GT_OLD","GT_OTHERS"]] = manta_sv_filter_type["GT_RAW"].str.split(":", n=1,expand = True)
manta_sv_filter_type["SVLEN"] = manta_sv_filter_type.END - manta_sv_filter_type.POS

# Filter part
manta_sv_svlen_min200_max100000_del_dup = manta_sv_filter_type[(manta_sv_filter_type.SVLEN >= 200) & (manta_sv_filter_type.SVLEN < 100000) & (manta_sv_filter_type.SVTYPE.isin(["DEL", "DUP"]))].copy()
manta_sv_svlen_min100000_del_dup = manta_sv_filter_type[(manta_sv_filter_type.SVLEN >= 100000) & (manta_sv_filter_type.SVTYPE.isin(["DEL", "DUP"]))].copy()
manta_sv_svlen_max200 = manta_sv_filter_type[(manta_sv_filter_type.SVLEN < 200) & (manta_sv_filter_type.ALT != "<INS>")].copy()

# Other part
manta_sv_svlen_ins = manta_sv_filter_type[manta_sv_filter_type.ALT == "<INS>"].copy()
manta_sv_svlen_ins.rename(columns={"FILTER_RAW": "FILTER"}, inplace=True)
manta_sv_svlen_ins.rename(columns={"INFO_RAW": "INFO"}, inplace=True)
manta_sv_svlen_ins.rename(columns={"GT_RAW": "GT"}, inplace=True)

manta_sv_svlen_min200_inv = manta_sv_filter_type[(manta_sv_filter_type.SVLEN >= 200) & (manta_sv_filter_type.SVTYPE == "INV")].copy()
manta_sv_svlen_min200_inv.rename(columns={"FILTER_RAW": "FILTER"}, inplace=True)
manta_sv_svlen_min200_inv.rename(columns={"INFO_RAW": "INFO"}, inplace=True)
manta_sv_svlen_min200_inv.rename(columns={"GT_RAW": "GT"}, inplace=True)

manta_sv_other_type = manta_sv_raw[~(manta_sv_raw.SVTYPE.isin(["DEL","DUP","INS","INV"]))].copy()
manta_sv_other_type.rename(columns={"FILTER_RAW": "FILTER"}, inplace=True)
manta_sv_other_type.rename(columns={"INFO_RAW": "INFO"}, inplace=True)
manta_sv_other_type.rename(columns={"GT_RAW": "GT"}, inplace=True)

# SV Filter
# 1. Filter SVLEN <= 200bp by Paragraph
manta_sv_svlen_max200_merge = pd.merge(manta_sv_svlen_max200, manta_paragraph_sv_gt, how = "left", on = "ID")

def record_svtools(x):
    if ((x == "1/1") | (x == "0/1")):
        return "Manta,Paragraph"
    else:
        return "Manta"

manta_sv_svlen_max200_merge["SVTOOL"] = manta_sv_svlen_max200_merge.GT_NEW.apply(record_svtools)
manta_sv_svlen_max200_merge["FILTER"] = manta_sv_svlen_max200_merge.SVTOOL.apply(lambda x: "PASS" if x=="Manta,Paragraph" else "Non-Reproduction")
manta_sv_svlen_max200_merge["INFO"] = manta_sv_svlen_max200_merge.INFO_RAW + ";SVTOOLS=" + manta_sv_svlen_max200_merge.SVTOOL
manta_sv_svlen_max200_merge["GT"] = manta_sv_svlen_max200_merge.GT_NEW + ":" + manta_sv_svlen_max200_merge.GT_OTHERS

# 2. Filter 200bp < SVLEN <= 100000bp & SVTYPE in (DEL, DUP) by overlap with Lumpy or CNVnator
lumpy_overlap =  manta_sv_svlen_min200_max100000_del_dup.apply(lambda row: sv_overlap(lumpy_sv_info,row["CHROM"],row["POS"],row["END"],row["SVTYPE"]), axis=1)
cnvnator_overlap = manta_sv_svlen_min200_max100000_del_dup.apply(lambda row: sv_overlap(cnvnator_sv_info,row["CHROM"],row["POS"],row["END"],row["SVTYPE"]), axis=1)
both_overlap_result = cnvnator_overlap*2 + lumpy_overlap
manta_sv_svlen_min200_max100000_del_dup["OVERLAP_IND"] = both_overlap_result

def record_svtools(x):
    if x == 3:
        return "Manta,CNVnator,Lumpy"
    elif x == 2:
        return "Manta,CNVnator"
    elif x == 1:
        return "Manta,Lumpy"
    else:
        return "Manta"

manta_sv_svlen_min200_max100000_del_dup["SVTOOL"] = manta_sv_svlen_min200_max100000_del_dup.OVERLAP_IND.apply(record_svtools)
manta_sv_svlen_min200_max100000_del_dup["FILTER"] = manta_sv_svlen_min200_max100000_del_dup.OVERLAP_IND.apply(lambda x: "PASS" if x > 0 else "Non-Reproduction")
manta_sv_svlen_min200_max100000_del_dup["INFO"] = manta_sv_svlen_min200_max100000_del_dup.INFO_RAW + ";SVTOOLS=" + manta_sv_svlen_min200_max100000_del_dup.SVTOOL
manta_sv_svlen_min200_max100000_del_dup.rename(columns={"GT_RAW": "GT"}, inplace=True)

# 3. Filter SVLEN > 100000bp & SVTYPE in (DEL, DUP) by overlap with CNVnator
lumpy_overlap =  manta_sv_svlen_min100000_del_dup.apply(lambda row: sv_overlap(lumpy_sv_info,row["CHROM"],row["POS"],row["END"],row["SVTYPE"]), axis=1)
cnvnator_overlap = manta_sv_svlen_min100000_del_dup.apply(lambda row: sv_overlap(cnvnator_sv_info,row["CHROM"],row["POS"],row["END"],row["SVTYPE"]), axis=1)
both_overlap_result = cnvnator_overlap*2 + lumpy_overlap
manta_sv_svlen_min100000_del_dup["OVERLAP_IND"] = both_overlap_result

def record_svtools(x):
    if x == 3:
        return "Manta,CNVnator,Lumpy"
    elif x == 2:
        return "Manta,CNVnator"
    elif x == 1:
        return "Manta,Lumpy"
    else:
        return "Manta"

manta_sv_svlen_min100000_del_dup["SVTOOL"] = manta_sv_svlen_min100000_del_dup.OVERLAP_IND.apply(record_svtools)
manta_sv_svlen_min100000_del_dup["FILTER"] = manta_sv_svlen_min100000_del_dup.OVERLAP_IND.apply(lambda x: "PASS" if x > 1 else "Non-Reproduction")
manta_sv_svlen_min100000_del_dup["INFO"] = manta_sv_svlen_min100000_del_dup.INFO_RAW + ";SVTOOLS=" + manta_sv_svlen_min100000_del_dup.SVTOOL
manta_sv_svlen_min100000_del_dup.rename(columns={"GT_RAW": "GT"}, inplace=True)

# Merge filter result
final_sv_filter_result = pd.concat([manta_sv_svlen_max200_merge[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]],
                                    manta_sv_svlen_min200_max100000_del_dup[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]],
                                    manta_sv_svlen_min100000_del_dup[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]],
                                    manta_sv_svlen_ins[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]],
                                    manta_sv_svlen_min200_inv[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]],
                                    manta_sv_other_type[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]]])
final_sv_filter_result.to_csv(output_prefix + "_qc_by_lumpy_cnvnator_paragraph_mark_no_header.vcf", sep="\t", index=False, header=False)

# Remove low quality sv
final_sv_filter_result_filter = pd.concat([manta_sv_svlen_max200_merge.loc[manta_sv_svlen_max200_merge.FILTER == "PASS"][["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]],
                                           manta_sv_svlen_min200_max100000_del_dup.loc[manta_sv_svlen_min200_max100000_del_dup.FILTER == "PASS"][["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]],
                                           manta_sv_svlen_min100000_del_dup.loc[manta_sv_svlen_min100000_del_dup.FILTER == "PASS"][["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]], 
                                           manta_sv_svlen_min200_inv[["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GT"]]])

final_sv_filter_result_filter.sort_values(by=['CHROM', 'POS'], inplace=True, ascending=True)
final_sv_filter_result_filter.to_csv(output_prefix + "_qc_by_lumpy_cnvnator_paragraph_filtered_no_header.vcf", sep="\t", index=False, header=False)