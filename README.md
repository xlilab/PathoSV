# PathoSV: A Transcriptome-Aware Framework for Prioritizing Pathogenic Structural Variants

**PathoSV** is a computational framework to identify and prioritize pathogenic structural variants (SVs) by quantifying their dosage-disrupting impact within a specific tissue context. It is particularly effective for diagnosing rare genetic disorders, such as neurodegenerative conditions.

---

### Live Web Server
For users who prefer a graphical interface or do not wish to install the tool locally, we provide a user-friendly web server:  
**[https://www.biosino.org/pathosv/](https://www.biosino.org/pathosv/)**

The web server includes an advanced feature using a Retrieval-Augmented Generation (RAG) framework with the DeepSeek Large Language Model to synthesize literature evidence and aid in interpreting potential relationships between an SV, the affected gene, and clinical phenotypes.

---

## Table of Contents
- [Introduction](#introduction)  
- [Framework Overview](#framework-overview)  
- [Installation](#installation)  
- [Annotation Data](#annotation-data)  
- [Usage Workflow](#usage-workflow)  
  - [Step 1 (Optional): SV Merging and Pre-processing](#step-1-optional-sv-merging-and-pre-processing)  
  - [Step 2: Background Allele Frequency (AF) Annotation](#step-2-background-allele-frequency-af-annotation)  
  - [Step 3: Pathogenicity Annotation](#step-3-pathogenicity-annotation)  
- [Example Output](#example-output)  
- [Repository Structure](#repository-structure)  
- [Data and Code Availability](#data-and-code-availability)  
- [How to Cite](#how-to-cite)  

---

## Introduction
Structural variants (SVs) are significant contributors to human genetic diversity and disease, but identifying which variants are pathogenic remains a central challenge in genomics. The pathogenicity of an SV is fundamentally tied to its effect on gene dosage within a relevant biological context.

PathoSV addresses this challenge by implementing a "transcriptome-aware" strategy. It moves beyond simple variant calling to a quantitative assessment of dosage disruption. The core innovation is the **Transcript Disruption Ratio (TDR)**, a metric that quantifies the proportion of a gene's functional transcript output (in a specific tissue) that is disrupted by an SV. By integrating tissue-specific expression data, PathoSV can effectively prioritize SVs most likely to be pathogenic.

---

## Framework Overview
The PathoSV workflow is a rigorous, multi-step filtering cascade designed to systematically isolate candidate dosage-disrupting SVs from tens of thousands of raw calls.

1. **High-Sensitivity SV Detection** â The framework begins with enhanced SV detection by integrating evidence from complementary callers (e.g., Manta, CNVnator, Lumpy) that leverage both breakpoint (BP) and read-depth (RD) signals.
2. **Population Filtering** â A stringent filter removes common polymorphisms by leveraging allele frequencies from large-scale, multi-ancestral reference cohorts (1000 Genomes Project East Asian and GTEx European). Rare variants are defined as those with a combined background allele frequency < 0.01.
3. **Transcriptome-Aware Prioritization** â The final and most crucial step uses the tissue-specific **Transcript Disruption Ratio (TDR)** to prioritize variants with the highest predicted functional impact. A TDR > 0.25 is used as the threshold to identify variants with a substantial biological effect.

This cascade dramatically reduces the number of candidates, leaving an average of only 2â5 high-confidence SVs per individual for final clinical interpretation.

---

## Installation
PathoSV is written in Python 3 and requires the `numpy` and `pandas` libraries.

```bash
git clone https://github.com/xlilab/PathoSV.git
cd PathoSV
pip install numpy pandas
```

---

## Annotation Data
The annotation scripts require several reference data files. Place them in a directory named `ref_dir/` within your working directory.

- `1KG_china_gtex_manta_qc_jasmine_merge0.8_max_af.txt` â Background SV AF database from 196 East Asian (1KGP) and 838 European (GTEx) individuals.  
- `gene_gencodev26_OMIM_GO_info.txt` â Gene annotation file combining GENCODE v26, OMIM, and GO.  
- `gencode.v26.annotation_exon_info.txt` â Exon/transcript annotation from GENCODE v26.  
- `constraint_z_genome_1kb.qc.download.txt.gz` â Gnocchi genomic constraint scores.  
- `55tissues_p10_v26_transcript_tpm_mean.txt` â Mean transcript TPM across 55 tissues (GTEx + retina).  
- `clinvar_20241027_sv_info.txt` â ClinVar SV annotations.

You may use custom annotation files, but they must follow the same format.

---

## Usage Workflow

### Step 1 (Optional): SV Merging and Pre-processing
Merge SV callsets from multiple samples/callers with [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR):

```bash
ls *.vcf > sample_files.txt
SURVIVOR merge sample_files.txt 1000 1 1 1 0 50 merged_svs.vcf
```

Filter problematic regions (e.g., HLA, ALT contigs) using `bedtools subtract` with exclusion lists. Scripts `filter_sv.sh` and `filter_sv_add_chr.sh` are provided.

---

### Step 2: Background Allele Frequency (AF) Annotation
```bash
python src/sv_background_af_annotation.py     --backgound ref_dir/1KG_china_gtex_manta_qc_jasmine_merge0.8_max_af.txt     --sv your_input_svs.vcf     --o output_prefix     --threshold 0.01     --rare True
```
Outputs:
- `[output]_background_af_annotation.txt`
- `[output]_background_af_annotation_rare.txt`

---

### Step 3: Pathogenicity Annotation
```bash
python src/sv_pathogenic_annotaion.py     --sv rare_svs.txt     -o final_output     --tissue "Retina"     --gene ref_dir/gene_gencodev26_OMIM_GO_info.txt     --exon ref_dir/gencode.v26.annotation_exon_info.txt     --gnocchi ref_dir/constraint_z_genome_1kb.qc.download.txt.gz     --clinvar ref_dir/clinvar_20241027_sv_info.txt     --tpm_trans ref_dir/55tissues_p10_v26_transcript_tpm_mean.txt
```

---

## Example Output
Columns:
- `SYMBOL` â gene symbol  
- `Consequence` â e.g., exon overlap  
- `Gnocchi` â max constraint score in SV region  
- `Phenotypes` â OMIM terms  
- `GO_BP`, `GO_CC`, `GO_MF` â GO terms  
- `ClinVar_CLNSIG` â ClinVar pathogenicity  
- `Sum_truncated_trascript_tpm` â disrupted transcript TPM sum  
- `Gene_tpm` â total gene TPM in tissue  
- `Truncated_ratio` â TDR (prioritize if > 0.25)  


---

## Data and Code Availability
- **PathoSV Code** â available in this repo.  
- **MoDs Metric** â [https://github.com/xlilab/MoDs](https://github.com/xlilab/MoDs)  
- **Public Data** â GTEx, 1000 Genomes Project  
- **Study Data** â NODE: OEP004860 (IRD), OEZ00021280 (HA)  

---

## License
This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)** â see the [LICENSE](LICENSE) file for details.

---

## How to Cite
Liu, X., Chen, Z., Jiang, Q., *et al.*  
*Transcriptome-aware identification of rare dosage-disrupting structural variants enhances diagnostic yield in neurodegenerative disorders.* (Journal, Year)  

```bibtex
# BibTeX entry will be provided upon publication.
```
