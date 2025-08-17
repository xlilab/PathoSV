# working environment is sv-anno

# Import Python Packages
import numpy as np
import pandas as pd
import argparse
import logging

# Logging Set
logger = logging.getLogger('SV-ANNO')
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
logger.addHandler(ch)

# Get Parameter
parser = argparse.ArgumentParser(description='SV transcriptome pathogenic annotation')
parser.add_argument('--gene', type=str, default = "None", metavar = "", help='Gene and phenotype annotation file.')
parser.add_argument('--exon', type=str, default = "None", metavar = "", help='Exon annotation file.')
parser.add_argument('--gnocchi', type=str, default = "None", metavar = "", help='Gnocchi score file which is a genomic mutational constraint map calculated from gnomAD dataset.')
parser.add_argument('--clinvar', type=str, default = "None", metavar = "", help='Clinvar SV info file.')
parser.add_argument('--tissue', type=str, default = "None", metavar = "", help='Interesting tissue for annotation. If value is None, output file will not contain transcriptom annotation.')
parser.add_argument('--tpm_trans', type=str, default = "None", metavar = "", help='Transcript TPM matrix file. Default is None.')
parser.add_argument('--mods', type=str, default = "None", metavar = "", help='Path to the Map of Dosage sensitivity (MoDs) file. Default is None.')
parser.add_argument('--sv', type=str, required=True, metavar = "", help='SV file which contains CHROM, START, END, SVTYPE at least.')
parser.add_argument('-o', '--output_prefix', type=str, required=True, metavar = "", help='Output file prefix.')
args = parser.parse_args()

# Load Files
if args.gene != "None":
    hg38_gene = pd.read_table(args.gene)
else:
    hg38_gene = pd.read_table("ref_dir/gene_gencodev26_OMIM_GO_info.txt.gz")
if args.exon != "None":
    hg38_exon = pd.read_table(args.exon)
else:
    hg38_exon = pd.read_table("ref_dir/gencode.v26.annotation_exon_info.txt.gz")
if args.gnocchi != "None":
    gnocchi_info = pd.read_table(args.gnocchi)[["chrom", "start", "end", "z"]]
else:
    gnocchi_info = pd.read_table("ref_dir/constraint_z_genome_1kb.qc.download.txt.gz")[["chrom", "start", "end", "z"]]
if args.clinvar != "None":
    clinvar_info = pd.read_table(args.clinvar, low_memory=False)
else:
    clinvar_info = pd.read_table("ref_dir/clinvar_20241027_sv_info.txt.gz", low_memory=False)

sv_info = pd.read_table(args.sv, low_memory=False)
sv_info["CHROM"] = sv_info["CHROM"].astype("str")
sv_info["IND"] = sv_info["CHROM"].str.cat([sv_info["START"].astype("str"), sv_info["END"].astype("str"), sv_info["SVTYPE"].astype("str")], sep="-")

output_prefix = args.output_prefix
select_tissue = args.tissue

if select_tissue != "None":
    if args.tpm_trans != "None":
        trans_tpm = pd.read_table(args.tpm_trans, low_memory=False)[["ENSG", "ENST", select_tissue]]
    else:
        transcript_tpm = pd.read_table("ref_dir/55tissues_p10_v26_transcript_tpm_mean.txt.gz", low_memory=False)[["ENSG", "ENST", select_tissue]]

# annotation sv to genome
tmp_IND = []
tmp_SYMBOL = []
tmp_GENE_TYPE = []
tmp_Consequence = []
tmp_ENSG = []
tmp_ENST = []
tmp_TRANS_TYPE = []
tmp_Gnocchi = []
tmp_MIM = []
tmp_Phenotypes = []
tmp_GO_MF = []
tmp_GO_BP = []
tmp_GO_CC = []
tmp_ClinVar_CLNSIG = []
tmp_ClinVar_CLNDN = []
tmp_ClinVar_ALLELEID = []

for i in range(len(sv_info)):
    # extract same chrom and svtype annotation
    tmp_sv = sv_info.loc[i, :].copy()
    
    if "chr" in tmp_sv.CHROM:
        tmp_hg38_gene = hg38_gene[hg38_gene.CHROM == tmp_sv.CHROM].copy()
        tmp_gnocchi_info = gnocchi_info[gnocchi_info.chrom == tmp_sv.CHROM].copy()
        tmp_clinvar_info = clinvar_info[clinvar_info.CHROM == tmp_sv.CHROM].copy()
    else:
        tmp_hg38_gene = hg38_gene[hg38_gene.CHROM == "chr" + tmp_sv.CHROM].copy()
        tmp_gnocchi_info = gnocchi_info[gnocchi_info.chrom == "chr" + tmp_sv.CHROM].copy()
        tmp_clinvar_info = clinvar_info[clinvar_info.CHROM == "chr" + tmp_sv.CHROM].copy()
        
    if tmp_sv.SVTYPE in ["INV", "Inversion"]:
        tmp_clinvar_info = tmp_clinvar_info[tmp_clinvar_info.SVTYPE == "Inversion"]
    elif tmp_sv.SVTYPE in ["DEL", "Deletion", "copy_number_loss"]:
        tmp_clinvar_info = tmp_clinvar_info[tmp_clinvar_info.SVTYPE.isin(["Deletion", "copy_number_loss"])]
    elif tmp_sv.SVTYPE in ["DUP", "Duplication", "copy_number_gain"]:
        tmp_clinvar_info = tmp_clinvar_info[tmp_clinvar_info.SVTYPE.isin(["Duplication", "copy_number_gain"])]
    elif tmp_sv.SVTYPE in ["INS", "Insertion"]:
        tmp_clinvar_info = tmp_clinvar_info[tmp_clinvar_info.SVTYPE.isin(["Duplication", "copy_number_gain"])]
    else:
        logger.error("SVTYPE is unrecognizable")
        
    # overlap with gene
    if tmp_sv.SVTYPE == "INS":
        # only annote the start location of insertion variants
        tmp_hg38_gene["e1e2"] = np.sign(tmp_hg38_gene.END-tmp_sv.START)
        tmp_hg38_gene["e1s2"] = np.sign(tmp_hg38_gene.END-tmp_sv.START)
        tmp_hg38_gene["e2s1"] = np.sign(tmp_sv.START-tmp_hg38_gene.START)
        tmp_hg38_gene["s2s1"] = np.sign(tmp_sv.START-tmp_hg38_gene.START)

        tmp_gnocchi_info["e1e2"] = np.sign(tmp_gnocchi_info.end-tmp_sv.END)
        tmp_gnocchi_info["e1s2"] = np.sign(tmp_gnocchi_info.end-tmp_sv.START)
        tmp_gnocchi_info["e2s1"] = np.sign(tmp_sv.END-tmp_gnocchi_info.start)
        tmp_gnocchi_info["s2s1"] = np.sign(tmp_sv.START-tmp_gnocchi_info.start)
    else:
        tmp_hg38_gene["e1e2"] = np.sign(tmp_hg38_gene.END-tmp_sv.END)
        tmp_hg38_gene["e1s2"] = np.sign(tmp_hg38_gene.END-tmp_sv.START)
        tmp_hg38_gene["e2s1"] = np.sign(tmp_sv.END-tmp_hg38_gene.START)
        tmp_hg38_gene["s2s1"] = np.sign(tmp_sv.START-tmp_hg38_gene.START)

        tmp_gnocchi_info["e1e2"] = np.sign(tmp_gnocchi_info.end-tmp_sv.END)
        tmp_gnocchi_info["e1s2"] = np.sign(tmp_gnocchi_info.end-tmp_sv.START)
        tmp_gnocchi_info["e2s1"] = np.sign(tmp_sv.END-tmp_gnocchi_info.start)
        tmp_gnocchi_info["s2s1"] = np.sign(tmp_sv.START-tmp_gnocchi_info.start)

        
    # extract overlap annotation
    tmp_overlap_gene_list = tmp_hg38_gene[~(((tmp_hg38_gene.e1e2 < 0) & (tmp_hg38_gene.e1s2 < 0) &
                                             (tmp_hg38_gene.s2s1 > 0) & (tmp_hg38_gene.e2s1 > 0)) |
                                            ((tmp_hg38_gene.e1e2 > 0) & (tmp_hg38_gene.e1s2 > 0) &
                                             (tmp_hg38_gene.s2s1 < 0) & (tmp_hg38_gene.e2s1 < 0)))].copy()
    
    tmp_gnocchi_info = tmp_gnocchi_info[~(((tmp_gnocchi_info.e1e2 < 0) & (tmp_gnocchi_info.e1s2 < 0) &
                                           (tmp_gnocchi_info.s2s1 > 0) & (tmp_gnocchi_info.e2s1 > 0)) |
                                          ((tmp_gnocchi_info.e1e2 > 0) & (tmp_gnocchi_info.e1s2 > 0) &
                                           (tmp_gnocchi_info.s2s1 < 0) & (tmp_gnocchi_info.e2s1 < 0)))].copy()
    tmp_gnocchi_score = tmp_gnocchi_info.z.max()
    
    if len(tmp_clinvar_info) == 0:
        tmp_clinvar_sig = "-"
        tmp_clinvar_pheno = "-"
        tmp_clinvar_id = "-"
    else:
        tmp_clinvar_info["e1e2"] = np.sign(tmp_clinvar_info.END-tmp_sv.END)
        tmp_clinvar_info["e1s2"] = np.sign(tmp_clinvar_info.END-tmp_sv.START)
        tmp_clinvar_info["e2s1"] = np.sign(tmp_sv.END-tmp_clinvar_info.START)
        tmp_clinvar_info["s2s1"] = np.sign(tmp_sv.START-tmp_clinvar_info.START)
            
        tmp_clinvar_var_list = tmp_clinvar_info[~(((tmp_clinvar_info.e1e2 < 0) & (tmp_clinvar_info.e1s2 < 0) & 
                                                   (tmp_clinvar_info.s2s1 > 0) & (tmp_clinvar_info.e2s1 > 0)) |
                                                  ((tmp_clinvar_info.e1e2 > 0) & (tmp_clinvar_info.e1s2 > 0) & 
                                                   (tmp_clinvar_info.s2s1 < 0) & (tmp_clinvar_info.e2s1 < 0)))].copy()
        if len(tmp_clinvar_var_list) == 0:
            tmp_clinvar_sig = "-"
            tmp_clinvar_pheno = "-"
            tmp_clinvar_id = "-"
        elif len(tmp_clinvar_var_list) > 0:
            def calculate_overlap_percent(tmp_row):
                return (min([tmp_row.END, tmp_sv.END]) - max([tmp_row.START, tmp_sv.START])) / max([(tmp_row.END - tmp_row.START), (tmp_sv.END - tmp_sv.START)])
            tmp_clinvar_var_list["overlap_percent"] = tmp_clinvar_var_list.apply(calculate_overlap_percent, axis=1)
            tmp_clinvar_var_list_filter = tmp_clinvar_var_list.loc[tmp_clinvar_var_list.overlap_percent >0.8].copy()
            if len(tmp_clinvar_var_list_filter) > 0:
                tmp_clinvar_sig = tmp_clinvar_var_list_filter.CLNSIG.tolist()[0]
                tmp_clinvar_pheno = tmp_clinvar_var_list_filter.CLNDN.tolist()[0]
                tmp_clinvar_id = tmp_clinvar_var_list_filter.ALLELEID.tolist()[0]
            else:
                tmp_clinvar_sig = "-"
                tmp_clinvar_pheno = "-"
                tmp_clinvar_id = "-"
    
    # add annotation by transcript
    if len(tmp_overlap_gene_list) > 0:
        for j in tmp_overlap_gene_list.index:
            if ((tmp_overlap_gene_list.loc[j, "e1e2"] < 0) &
                (tmp_overlap_gene_list.loc[j, "e1s2"] > 0) &
                (tmp_overlap_gene_list.loc[j, "s2s1"] < 0) &
                (tmp_overlap_gene_list.loc[j, "e2s1"] > 0)):
                tmp_hg38_exon = hg38_exon[hg38_exon.ENSG == tmp_overlap_gene_list.loc[j, "ENSG"]].copy()
                for k in tmp_hg38_exon.loc[tmp_hg38_exon["ENSG"] == tmp_overlap_gene_list.loc[j, "ENSG"], "ENST"]:
                    tmp_IND.append(tmp_sv.IND)
                    tmp_SYMBOL.append(tmp_overlap_gene_list.loc[j, "SYMBOL"])
                    tmp_GENE_TYPE.append(tmp_overlap_gene_list.loc[j, "GENE_TYPE"])
                    tmp_Consequence.append("Overlap all gene")
                    tmp_ENSG.append(tmp_overlap_gene_list.loc[j, "ENSG"])
                    tmp_ENST.append(k)
                    tmp_TRANS_TYPE.append(tmp_hg38_exon.loc[tmp_hg38_exon.ENST == k, "TRANS_TYPE"].unique().tolist()[0])
                    tmp_Gnocchi.append(tmp_gnocchi_score)
                    tmp_MIM.append(tmp_overlap_gene_list.loc[j, "MIM"])
                    tmp_Phenotypes.append(tmp_overlap_gene_list.loc[j, "Phenotypes"])
                    tmp_GO_MF.append(tmp_overlap_gene_list.loc[j, "GO_MF"])
                    tmp_GO_BP.append(tmp_overlap_gene_list.loc[j, "GO_BP"])
                    tmp_GO_CC.append(tmp_overlap_gene_list.loc[j, "GO_CC"])
                    tmp_ClinVar_CLNSIG.append(tmp_clinvar_sig)
                    tmp_ClinVar_CLNDN.append(tmp_clinvar_pheno)
                    tmp_ClinVar_ALLELEID.append(tmp_clinvar_id)
            else:
                tmp_hg38_exon = hg38_exon[hg38_exon.ENSG == tmp_overlap_gene_list.loc[j, "ENSG"]].copy()
                if tmp_sv.SVTYPE in ["INS","Insertion"]:
                    tmp_hg38_exon["e1s2"] = np.sign(tmp_hg38_exon.END-tmp_sv.START)
                    tmp_hg38_exon["e2s1"] = np.sign(tmp_sv.START-tmp_hg38_exon.START)
                else:
                    tmp_hg38_exon["e1s2"] = np.sign(tmp_hg38_exon.END-tmp_sv.START)
                    tmp_hg38_exon["e2s1"] = np.sign(tmp_sv.END-tmp_hg38_exon.START)
                tmp_hg38_exon["Condition"] = tmp_hg38_exon.e1s2 * tmp_hg38_exon.e2s1
                if len(tmp_hg38_exon.loc[tmp_hg38_exon.Condition > 0,]) > 0:
                    for k in tmp_hg38_exon[tmp_hg38_exon.Condition > 0]["ENST"].unique():
                        tmp_IND.append(tmp_sv.IND)
                        tmp_SYMBOL.append(tmp_overlap_gene_list.loc[j, "SYMBOL"])
                        tmp_GENE_TYPE.append(tmp_overlap_gene_list.loc[j, "GENE_TYPE"])
                        tmp_Consequence.append("Overlap " + str(", ".join(m for m in tmp_hg38_exon.loc[tmp_hg38_exon.ENST == k, "CLASS"].unique())))
                        tmp_ENSG.append(tmp_overlap_gene_list.loc[j, "ENSG"])
                        tmp_ENST.append(k)
                        tmp_TRANS_TYPE.append(tmp_hg38_exon.loc[tmp_hg38_exon.ENST == k, "TRANS_TYPE"].unique().tolist()[0])
                        tmp_Gnocchi.append(tmp_gnocchi_score)
                        tmp_MIM.append(tmp_overlap_gene_list.loc[j, "MIM"])
                        tmp_Phenotypes.append(tmp_overlap_gene_list.loc[j, "Phenotypes"])
                        tmp_GO_MF.append(tmp_overlap_gene_list.loc[j, "GO_MF"])
                        tmp_GO_BP.append(tmp_overlap_gene_list.loc[j, "GO_BP"])
                        tmp_GO_CC.append(tmp_overlap_gene_list.loc[j, "GO_CC"])
                        tmp_ClinVar_CLNSIG.append(tmp_clinvar_sig)
                        tmp_ClinVar_CLNDN.append(tmp_clinvar_pheno)
                        tmp_ClinVar_ALLELEID.append(tmp_clinvar_id)
                else:
                    tmp_IND.append(tmp_sv.IND)
                    tmp_SYMBOL.append(tmp_overlap_gene_list.loc[j, "SYMBOL"])
                    tmp_GENE_TYPE.append(tmp_overlap_gene_list.loc[j, "GENE_TYPE"])
                    tmp_Consequence.append("Overlap intron")
                    tmp_ENSG.append(tmp_overlap_gene_list.loc[j, "ENSG"])
                    tmp_ENST.append("-")
                    tmp_TRANS_TYPE.append("-")
                    tmp_Gnocchi.append(tmp_gnocchi_score)
                    tmp_MIM.append(tmp_overlap_gene_list.loc[j, "MIM"])
                    tmp_Phenotypes.append(tmp_overlap_gene_list.loc[j, "Phenotypes"])
                    tmp_GO_MF.append(tmp_overlap_gene_list.loc[j, "GO_MF"])
                    tmp_GO_BP.append(tmp_overlap_gene_list.loc[j, "GO_BP"])
                    tmp_GO_CC.append(tmp_overlap_gene_list.loc[j, "GO_CC"])
                    tmp_ClinVar_CLNSIG.append(tmp_clinvar_sig)
                    tmp_ClinVar_CLNDN.append(tmp_clinvar_pheno)
                    tmp_ClinVar_ALLELEID.append(tmp_clinvar_id)

    else:
        tmp_IND.append(tmp_sv.IND)
        tmp_SYMBOL.append("-")
        tmp_GENE_TYPE.append("-")
        tmp_Consequence.append("Overlap intergenic")
        tmp_ENSG.append("-")
        tmp_ENST.append("-")
        tmp_TRANS_TYPE.append("-")
        tmp_Gnocchi.append(tmp_gnocchi_score)
        tmp_MIM.append("-")
        tmp_Phenotypes.append("-")
        tmp_GO_MF.append("-")
        tmp_GO_BP.append("-")
        tmp_GO_CC.append("-")
        tmp_ClinVar_CLNSIG.append(tmp_clinvar_sig)
        tmp_ClinVar_CLNDN.append(tmp_clinvar_pheno)
        tmp_ClinVar_ALLELEID.append(tmp_clinvar_id)

result_anno_df = pd.DataFrame({"IND": tmp_IND,"SYMBOL": tmp_SYMBOL, "GENE_TYPE": tmp_GENE_TYPE, "Consequence": tmp_Consequence,
                               "ENSG": tmp_ENSG, "ENST": tmp_ENST, "TRANS_TYPE": tmp_TRANS_TYPE,"Gnocchi": tmp_Gnocchi, 
                               "MIM": tmp_MIM, "Phenotypes": tmp_Phenotypes,"GO_MF": tmp_GO_MF, "GO_BP": tmp_GO_BP, "GO_CC": tmp_GO_CC,
                               "ClinVar_CLNSIG": tmp_ClinVar_CLNSIG, "ClinVar_CLNDN": tmp_ClinVar_CLNDN,
                               "ClinVar_ALLELEID": tmp_ClinVar_ALLELEID})
result_df = pd.merge(left=sv_info, right=result_anno_df, on="IND")
del(result_anno_df)


# add tpm annotation
if select_tissue != "None":
    result_df_tpm = pd.merge(left=result_df, right=transcript_tpm[["ENST", select_tissue]], on="ENST", how="left")

    result_df_tpm.loc[result_df_tpm.Consequence.str.contains("exon"),"Consequence"] = "Overlap exon"
    result_df_tpm_unique = result_df_tpm[sv_info.columns.tolist() + ["SYMBOL", "ENSG", "GENE_TYPE", "Consequence", "Gnocchi","MIM","Phenotypes","GO_MF", "GO_BP", "GO_CC","ClinVar_CLNSIG", "ClinVar_CLNDN", "ClinVar_ALLELEID"]].drop_duplicates()
    result_df_tpm_unique.index = range(len(result_df_tpm_unique))

    tpm_sum_truncated = []
    tpm_sum_all = []
    tpm_percent = []
    tpm_rank = []

    for i in range(len(result_df_tpm_unique)):
        if (result_df_tpm_unique.loc[i,"Consequence"] != "Overlap intron") & (result_df_tpm_unique.loc[i,"Consequence"] != "Overlap intergenic"):
            tmp_trans_tpm_df = result_df_tpm.loc[(result_df_tpm.IND == result_df_tpm_unique.iloc[i,:].IND) & (result_df_tpm.ENSG == result_df_tpm_unique.iloc[i,:].ENSG), ["ENST", select_tissue]].drop_duplicates().copy()
            tmp_gene_tpm_df = transcript_tpm.loc[transcript_tpm.ENSG == result_df_tpm_unique.iloc[i,:].ENSG, ["ENST", select_tissue]]
            tmp_sum_truncated = round(tmp_trans_tpm_df[select_tissue].sum(),3)
            tmp_sum_all = round(tmp_gene_tpm_df[select_tissue].sum(),3)
            tmp_gene_tpm_df["Rank"] = tmp_gene_tpm_df[select_tissue].rank(ascending=False).astype(int)
            if len(tmp_gene_tpm_df.loc[tmp_gene_tpm_df.ENST.isin(tmp_trans_tpm_df.ENST),]) > 0:
                tmp_truncated_top_rank = min(tmp_gene_tpm_df.loc[tmp_gene_tpm_df.ENST.isin(tmp_trans_tpm_df.ENST),"Rank"])
            else:
                tmp_truncated_top_rank = "NAN"
            if ((tmp_sum_all == 0)):
                tpm_sum_truncated.append(tmp_sum_truncated)
                tpm_sum_all.append(tmp_sum_all)
                tpm_percent.append(0)
                tpm_rank.append(tmp_truncated_top_rank)
            else:
                tpm_sum_truncated.append(tmp_sum_truncated)
                tpm_sum_all.append(tmp_sum_all)
                tpm_percent.append(round(tmp_sum_truncated/tmp_sum_all, 2))
                tpm_rank.append(tmp_truncated_top_rank)
        else:
            tpm_sum_truncated.append("-")
            tpm_sum_all.append("-")
            tpm_percent.append("-")
            tpm_rank.append("-")

    tpm_percent_df = pd.DataFrame({"Sum_truncated_trascript_TPM": tpm_sum_truncated,
                                   "Gene_TPM": tpm_sum_all,
                                   "TDR": tpm_percent,
                                   "Top_truncated_trascript_TPM_rank": tpm_rank})
    result_df_tpm_unique.index = range(len(result_df_tpm_unique))
    result_df_tpm_unique_tpm_percent_df = pd.concat([result_df_tpm_unique, tpm_percent_df], axis=1)
    
    if (args.mods != "None"):
        mods_anno = pd.read_table(args.mods)[["ENSG", select_tissue]]
        mods_anno.columns = ["ENSG", "MODS"]
        result_df_tpm_unique_tpm_percent_mods_df = pd.merge(left=result_df_tpm_unique_tpm_percent_df, right=mods_anno, how="left", on="ENSG")
        result_df_tpm_unique_tpm_percent_mods_df.to_csv(output_prefix + "_" +select_tissue + "_transcriptom_annotation.txt", sep="\t", index=False, header=True)
    else:
        result_df_tpm_unique_tpm_percent_df.to_csv(output_prefix + "_" +select_tissue + "_transcriptom_annotation.txt", sep="\t", index=False, header=True)
else:
    result_df.to_csv(output_prefix + "_genome_annotation.txt", sep="\t", index=False, header=True)
