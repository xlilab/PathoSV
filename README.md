# PathoSV: A Transcriptome-Aware Framework for Prioritizing Pathogenic Structural Variants

[![Python Version](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**PathoSV** is a computational framework to identify and prioritize pathogenic structural variants (SVs) by quantifying their dosage-disrupting impact within a specific tissue context. It is particularly effective for diagnosing rare genetic disorders, such as neurodegenerative conditions.

<img width="872" height="497" alt="PathoSV" src="https://github.com/user-attachments/assets/b54acf68-da49-46e4-bee6-67cd44fb2bb0" />

### Live Web Server
For users who prefer a graphical interface or do not wish to install the tool locally, we provide a user-friendly web server:
**[https://www.biosino.org/pathosv/](https://www.biosino.org/pathosv/)**

The web server includes an advanced feature using a Retrieval-Augmented Generation (RAG) framework with the DeepSeek Large Language Model to synthesize literature evidence and aid in interpreting potential relationships between an SV, the affected gene, and clinical phenotypes.

## Table of Contents
- [Introduction](#introduction)
- [Framework Overview](#framework-overview)
- [Installation](#installation)
- [Annotation Data](#annotation-data)
- [Usage Workflow](#usage-workflow)
  - [Step 1: SV Calling and Quality Control (QC)](#step-1-sv-calling-and-quality-control-qc)
  - [Step 2 (Optional): SV Merging and Pre-processing](#step-2-optional-sv-merging-and-pre-processing)
  - [Step 3: Background Allele Frequency (AF) Annotation](#step-3-background-allele-frequency-af-annotation)
  - [Step 4: Pathogenicity Annotation](#step-4-pathogenicity-annotation)
- [Example Output](#example-output)
- [Repository Structure](#repository-structure)
- [Data and Code Availability](#data-and-code-availability)
- [How to Cite](#how-to-cite)

## Introduction

Structural variants (SVs) are significant contributors to human genetic diversity and disease, but identifying which variants are pathogenic remains a critical challenge in genomics. The pathogenicity of an SV is fundamentally tied to its effect on gene dosage within a relevant biological context.

PathoSV addresses this challenge by implementing a "transcriptome-aware" strategy. It moves beyond simple variant calling to a quantitative assessment of dosage disruption. We introduce the **Transcript Disruption Ratio (TDR)**, a metric that quantifies the proportion of a gene's functional transcript output (in a specific tissue) that is disrupted by an SV. By integrating tissue-specific expression data, PathoSV can effectively prioritize SVs most likely to be pathogenic.

## Framework Overview

The PathoSV workflow is a rigorous, multi-step filtering cascade designed to systematically isolate candidate dosage-disrupting SVs from tens of thousands of raw calls.

1.  **High-Sensitivity SV Detection:** The framework begins with enhanced SV detection by integrating evidence from complementary callers (e.g., Manta, CNVnator, Lumpy) that leverage both breakpoint (BP) and read-depth (RD) signals.
2.  **Population Filtering:** A stringent filter removes common polymorphisms by leveraging allele frequencies from large-scale, multi-ancestral reference cohorts (1000 Genomes Project East Asian and GTEx European). Rare variants are defined as those with a combined background allele frequency < 0.01.
3.  **Transcriptome-Aware Prioritization:** The final and most crucial step uses the tissue-specific **Transcript Disruption Ratio (TDR)** to prioritize variants with the highest predicted functional impact. A TDR > 0.25 is used as the threshold to identify variants with a substantial biological effect.

This cascade dramatically reduces the number of candidates, leaving an average of only 2-5 high-confidence SVs per individual for final clinical interpretation.

## Installation

PathoSV is written in Python 3 and requires the `numpy` and `pandas` libraries.

1.  Clone the repository:
    ```bash
    git clone [https://github.com/xlilab/PathoSV.git](https://github.com/xlilab/PathoSV.git)
    cd PathoSV
    ```
2.  Ensure dependencies are installed:
    ```bash
    pip install numpy pandas
    ```

## Annotation Data

The annotation scripts require several reference data files. These files should be placed in a directory named `ref_dir/` within your working directory.

-   **`1KG_china_gtex_manta_qc_jasmine_merge0.8_max_af.txt`**: Background SV allele frequency database from 196 East Asian (1KGP) and 838 European (GTEx) individuals. Contains columns: `CHROM`, `START`, `END`, `SVTYPE`, `SVLEN`, `Background_AF`.
-   **`gene_gencodev26_OMIM_GO_info.txt`**: Gene annotation file combining information from GENCODE v26, OMIM, and Gene Ontology (GO).
-   **`gencode.v26.annotation_exon_info.txt`**: Exon and transcript annotation from GENCODE v26. Contains 9 columns including `CHROM`, `START`, `END`, `ENSG`, `ENST`, `SYMBOL`.
-   **`constraint_z_genome_1kb.qc.download.txt.gz`**: Gnocchi genomic constraint scores (z-scores).
-   **`55tissues_p10_v26_transcript_tpm_mean.txt`**: Mean transcript TPM values across 55 tissues (54 from GTEx + retina).
-   **`clinvar_20241027_sv_info.txt`**: SV information from the ClinVar database, including columns `CHROM`, `START`, `END`, `SVTYPE`, `CLNSIG`, `CLNDN`.

**Note:** You can use your own custom annotation files, but they must adhere to the same format and column structure.

## Usage Workflow

### Step 1: SV Calling and Quality Control (QC)

Accurate downstream analysis begins with high-quality SV calls. While **Manta** is a leading tool known for its speed and precise breakpoint detection, its reliance on paired-end and split-read evidence can lead to **false-positive** deletions and duplications because it doesn't use a read-depth strategy. 

To enhance accuracy, our workflow integrates evidence from complementary tools to filter Manta's output. We use the read-depth based callers **CNVnator** and **Lumpy**, along with the SV re-genotyper **Paragraph**, to remove false-positive variants and produce a high-confidence callset.

<img width="2012" height="1070" alt="FigureS_QC" src="https://github.com/user-attachments/assets/52d01a31-811b-48d1-854b-9b47482e14c0" />


#### QC Workflow Overview

1.  **Initial Calling**: Perform SV calling on each sample using **Manta**, **Lumpy**, and **CNVnator**.
2.  **Information Extraction**: Use `bedtools` to extract key information (chromosome, start, end, SV type, and length) from the Lumpy and CNVnator VCF files.
3.  **Partition Manta Calls**: Split the Manta SV calls into two groups based on length: those less than 200 bp and those 200 bp or greater.
4.  **Re-genotype Short SVs**: Input the short SVs (< 200 bp) into **Paragraph** for re-genotyping to confirm their presence.
5.  **Filter and Annotate**: Process the Manta calls, Lumpy/CNVnator data, and Paragraph output using our custom `manta_qc.py` script. This script applies the filtering logic and produces two headerless VCF files: one with filter annotations (`*_mark_no_header.vcf`) and another with only the filtered, high-confidence calls (`*_filtered_no_header.vcf`).
6.  **Finalize VCFs**: Add the original VCF headers back to the two output files.

The complete bash script for this pipeline can be found in `code/filter_manta_sv_pipeline.sh`, and the core Python logic is in `code/manta_qc.py`.

### Step 2 (Optional): SV Merging and Pre-processing

Before annotation, it is common to merge SV callsets from multiple samples or callers.

**SV Merging with Jasmine**

Jasmine consolidates SVs of the same type that have a significant reciprocal overlap across different samples.

```bash
# Create a file listing all VCF files to be merged
ls *.vcf > sample_files.txt

# Run Jasmine to merge the VCFs
# Key parameters include the file list, output file name, and overlap requirements.
jasmine file_list=sample_files.txt \
        out_file=merged_svs.vcf \
        --spec_overlap=0.8 \
        --threads=8
```



### Filtering Problematic Regions

It is also critical to filter out SVs in known problematic genomic regions (e.g., HLA, ALT contigs). We recommend using the exclusion list from the [SpeedSeq paper](https://github.com/hall-lab/speedseq/blob/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed) and `bedtools subtract`. The scripts `filter_sv.sh` and `filter_sv_add_chr.sh` are provided to assist with this process.

---

### Step 3: Background Allele Frequency (AF) Annotation

This script annotates your input SVs with their allele frequencies from our background population database.

**Usage:**

```bash
python /path/to/sv_background_af_annotation.py \
    --backgound /path/to/ref_dir/1KG_china_gtex_manta_qc_jasmine_merge0.8_max_af.txt \
    --sv your_input_svs.vcf \
    --o your_output_prefix \
    [--threshold 0.01] \
    [--rare True]
```

**Arguments:**

* `--backgound`: (Required) Path to the background AF database file.
* `--sv`: (Required) Input SV file (VCF format). Must have a header.
* `--o`: (Required) Output file prefix.
* `--threshold`: (Optional) AF threshold for defining rare variants. Default is `0.01`.
* `--rare`: (Optional) If `True`, outputs an additional file containing only rare SVs. Default is `True`.

**Output:**

* `[output]_background_af_annotation.txt`: Input SVs annotated with the `Background_AF` column.
* `[output]_background_af_annotation_rare.txt`: A subset of the above file containing only rare SVs.

---
### Step 4: Pathogenicity Annotation
This is the core script that annotates SVs with genomic information, constraint scores, ClinVar data, and the tissue-specific **Transcript Disruption Ratio (TDR)**.

**Usage:**
```bash
python /path/to/sv_pathogenic_annotaion.py \
    --sv your_rare_svs.txt \
    -o your_final_output_prefix \
    --tissue "Retina" \
    --gene /path/to/ref_dir/gene_gencodev26_OMIM_GO_info.txt \
    --exon /path/to/ref_dir/gencode.v26.annotation_exon_info.txt \
    --gnocchi /path/to/ref_dir/constraint_z_genome_1kb.qc.download.txt.gz \
    --clinvar /path/to/ref_dir/clinvar_20241027_sv_info.txt \
    --tpm_trans /path/to/ref_dir/55tissues_p10_v26_transcript_tpm_mean.txt
```

**Arguments:**

* `--sv`: (Required) Input SV file (must have a header with CHR, START, END, SVTYPE).
* `--o`: (Required) Output file prefix.
* `--tissue`: (Optional) Tissue of interest for TDR calculation (e.g., "Retina", "Cerebellum"). If `None`, transcriptome annotation is skipped.
* `--gene`, `--exon`, `--gnocchi`, `--clinvar`, `--tpm_trans`: (Optional) Paths to reference files. Defaults to filenames within a local `ref_dir/`.

---
### Example Output
The primary output is a text file where each row represents an SV-gene annotation. Key columns include:

* **SYMBOL**: Gene symbol.
* **Consequence**: How the SV intersects the gene (e.g., `Overlap exon`, `Overlap all gene`).
* **Gnocchi**: Max genomic constraint score in the SV region.
* **Phenotypes**: OMIM phenotypes and inheritance patterns.
* **GO_BP, GO_CC, GO_MF**: Gene Ontology terms.
* **ClinVar_CLNSIG**: ClinVar pathogenicity.
* **Sum_truncated_trascript_tpm**: Sum of TPMs for transcripts disrupted by the SV in the specified tissue.
* **Gene_tpm**: Total TPM of all transcripts for the gene in that tissue.
* **Truncated_ratio**: The calculated **TDR**. A value > 0.25 is prioritized as potentially pathogenic.
---
### Data and Code Availability
* **PathoSV Code:** The source code for the PathoSV pipeline is available in this repository.
* **MoDs:** The Map of Dosage sensitivity (MoDs) metric and code are available at [https://github.com/xlilab/MoDs](https://github.com/xlilab/MoDs).
* **Public Data:** Data from GTEx and the 1000 Genomes Project are publicly available.
* **Study Data:** Processed data and SV calls from the IRD and HA cohorts have been deposited in the National Omics Data Encyclopedia (NODE) under accession numbers OEP004860 (IRD) and OEZ00021280 (HA). Access is subject to committee approval.

---
### How to Cite
If you use PathoSV in your research, please cite our paper:

```bibtex
# BibTeX entry will be provided upon publication.

