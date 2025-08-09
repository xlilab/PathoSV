## Need tools: bcftools, Paragraph
## Need python package: pandas, numpy

# parameter
sv_raw_dir=$1
manta_vcf_name=$2
cnvnator_vcf_name=$3
lumpy_vcf_name=$4
ref_file=$5
bam_depth_read_info=$6
spid=$7

lumpy_dir=lumpy-sv

# Tools pathway
manta_convertInversion_tool=/home/liuxubing/miniconda3/envs/sv/share/manta-1.6.0-0/libexec/convertInversion.py
samtools_tool=/home/liuxubing/miniconda3/envs/sv/bin/samtools
paragraph_tools=/home/liuxubing/miniconda3/envs/paragraph/bin/multigrmpy.py
merge_tool=/home/liuxubing/work/sv/sv_filter/filter_manta_by_lumpy_cnvnator/manta_qc.py

# Extract Lumpy SV bed
cd $sv_raw_dir"/"$spid"/"$lumpy_dir

if [ ! -f lumpy_del_dup.bed ];then
	bcftools query -f "%CHROM\t%POS\t%END\t%SVTYPE\n" $lumpy_vcf_name | grep -v "BND" | grep -v "INV" | grep -v "INS" >lumpy_del_dup.bed
fi

# Extract cnvnator bed
cd $sv_raw_dir"/"$spid"/cnvnator"

if [ ! -f cnvnator_info.bed ];then
	if [ -f $cnvnator_vcf_name".gz" ];then
		gunzip $cnvnator_vcf_name".gz"
		rm $cnvnator_vcf_name".gz.tbi"
	fi
	bcftools query -f "%CHROM\t%POS\t%END\t%SVTYPE\n" -o cnvnator_info.bed $cnvnator_vcf_name
fi

# Split manta result
cd $sv_raw_dir"/"$spid"/manta"

if [ ! -f manta_pass_autosome_convertInversion.vcf ];then
	grep "#" $manta_vcf_name >manta_pass_autosome.vcf
	for i in {1..22} X Y
	do
	        awk -v chr="chr$i" 'BEGIN{FS="\t";OFS="\t"}{if($1==chr)print $0}' $manta_vcf_name >>manta_pass_autosome.vcf
	done

	$manta_convertInversion_tool $samtools_tool $ref_file manta_pass_autosome.vcf >manta_pass_autosome_convertInversion.vcf
fi

# split sv for paragraph
## need to attention:
## 1. The INS records need SEQ info in INFO column which could get from ALT column(which ALT != <INS>)
## 2. The ALT column of INV records need to overwrite as <INV>

if [ ! -f manta_pass_autosome_max_svlen200_for_paragraph.vcf ];then
	grep "#" manta_pass_autosome_convertInversion.vcf >manta_pass_autosome_header
	grep -v "#" manta_pass_autosome_convertInversion.vcf >manta_pass_autosome_convertInversion_noheader.vcf

# extract SVLEN <200 & SVTYPE != "BND"
	awk 'NR==FNR{ a[$1]; next }FNR in a' <(awk '{if($1 < 200)print NR}' <(bcftools query -f "%END\t%POS\n" manta_pass_autosome_convertInversion.vcf | awk 'BEGIN{FS="\t";OFS="\t"}{print $1-$2}')) <(grep -v "#" manta_pass_autosome_convertInversion.vcf) | grep -v "SVTYPE=BND" | grep -v "<INS>" >manta_pass_autosome_max_svlen200.vcf

# fix INS format
	bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\tEND=%END;SEQ=%ALT;SVTYPE=%SVTYPE;SVLEN=%SVLEN\n" <(cat manta_pass_autosome_header <(grep "SVTYPE=INS" manta_pass_autosome_max_svlen200.vcf)) >tmp_ins.vcf
	cat <(cut -f 1-8 manta_pass_autosome_header) <(grep -v "SVTYPE=INS" manta_pass_autosome_max_svlen200.vcf| cut -f 1-8) tmp_ins.vcf >manta_pass_autosome_max_svlen200_for_paragraph.vcf
	rm tmp_ins.vcf
fi

# run Paragraph
if [ ! -d $sv_raw_dir"/"$spid"/paragraph_max_svlen200" ];then
	mkdir $sv_raw_dir"/"$spid"/paragraph_max_svlen200"
fi

cd $sv_raw_dir"/"$spid"/paragraph_max_svlen200"

if [ ! -f $sv_raw_dir"/"$spid"/paragraph_max_svlen200/manta_pass_autosome_max_svlen200_paragraph_gt.vcf" ];then
	cat <(sed -n '1p' $bam_depth_read_info) <(grep $spid $bam_depth_read_info) >$sv_raw_dir"/"$spid"/paragraph_max_svlen200/bam_depth_length_info.txt"
	python3 $paragraph_tools -i $sv_raw_dir"/"$spid"/manta/manta_pass_autosome_max_svlen200_for_paragraph.vcf" -m $sv_raw_dir"/"$spid"/paragraph_max_svlen200/bam_depth_length_info.txt" -r $ref_file -o $sv_raw_dir"/"$spid"/paragraph_max_svlen200/paragraph_output"
	if [ -f $sv_raw_dir"/"$spid"/paragraph_max_svlen200/paragraph_output/genotypes.vcf.gz" ];then
	    	zcat $sv_raw_dir"/"$spid"/paragraph_max_svlen200/paragraph_output/genotypes.vcf.gz" | grep -v "#" >$sv_raw_dir"/"$spid"/paragraph_max_svlen200/manta_pass_autosome_max_svlen200_paragraph_gt.vcf"
	fi
fi

# Filter Manta SV by Lumpy, CNVnator and Paragraph result
cd $sv_raw_dir"/"$spid"/manta"

if [ -f $sv_raw_dir"/"$spid"/paragraph_max_svlen200/manta_pass_autosome_max_svlen200_paragraph_gt.vcf" ];then
	python3 $merge_tool -m $sv_raw_dir"/"$spid"/manta/manta_pass_autosome_convertInversion_noheader.vcf" -l $sv_raw_dir"/"$spid"/"$lumpy_dir"/lumpy_del_dup.bed" -c $sv_raw_dir"/"$spid"/cnvnator/cnvnator_info.bed" -p $sv_raw_dir"/"$spid"/paragraph_max_svlen200/manta_pass_autosome_max_svlen200_paragraph_gt.vcf" -o manta_pass_autosome

	cat manta_pass_autosome_header manta_pass_autosome_qc_by_lumpy_cnvnator_paragraph_mark_no_header.vcf >manta_pass_autosome_qc_by_lumpy_cnvnator_paragraph_mark.vcf
	cat manta_pass_autosome_header manta_pass_autosome_qc_by_lumpy_cnvnator_paragraph_filtered_no_header.vcf >manta_pass_autosome_qc_by_lumpy_cnvnator_paragraph_filtered.vcf

	rm $sv_raw_dir"/"$spid"/manta/manta_pass_autosome_qc_by_lumpy_cnvnator_paragraph_mark_no_header.vcf"
	rm $sv_raw_dir"/"$spid"/manta/manta_pass_autosome_qc_by_lumpy_cnvnator_paragraph_filtered_no_header.vcf"
fi
