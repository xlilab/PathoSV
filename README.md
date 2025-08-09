
# PathoSV: A Transcriptome-Aware Framework for Prioritizing Pathogenic Structural Variants

[![Python Version](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[cite_start]**PathoSV** is a computational framework to identify and prioritize pathogenic structural variants (SVs) by quantifying their dosage-disrupting impact within a specific tissue context[cite: 18, 28, 46]. [cite_start]It is particularly effective for diagnosing rare genetic disorders, such as neurodegenerative conditions[cite: 1, 29].

### Live Web Server
For users who prefer a graphical interface or do not wish to install the tool locally, we provide a user-friendly web server:
[cite_start]**[https://www.biosino.org/pathosv/](https://www.biosino.org/pathosv/)** [cite: 312]

[cite_start]The web server includes an advanced feature using a Retrieval-Augmented Generation (RAG) framework with the DeepSeek Large Language Model to synthesize literature evidence and aid in interpreting potential relationships between an SV, the affected gene, and clinical phenotypes[cite: 67, 313, 314].

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

## Introduction

[cite_start]Structural variants (SVs) are significant contributors to human genetic diversity and disease, but identifying which variants are pathogenic remains a central challenge in genomics[cite: 31, 32, 35]. [cite_start]The pathogenicity of an SV is fundamentally tied to its effect on gene dosage within a relevant biological context[cite: 36].

[cite_start]PathoSV addresses this challenge by implementing a "transcriptome-aware" strategy[cite: 36, 46]. [cite_start]It moves beyond simple variant calling to a quantitative assessment of dosage disruption[cite: 47]. [cite_start]The core innovation is the **Transcript Disruption Ratio (TDR)**, a metric that quantifies the proportion of a gene's functional transcript output (in a specific tissue) that is disrupted by an SV[cite: 20, 21, 56, 95]. [cite_start]By integrating tissue-specific expression data, PathoSV can effectively prioritize SVs most likely to be pathogenic[cite: 22, 28].

## Framework Overview

[cite_start]The PathoSV workflow is a rigorous, multi-step filtering cascade designed to systematically isolate candidate dosage-disrupting SVs from tens of thousands of raw calls[cite: 52, 121, 122].

1.  [cite_start]**High-Sensitivity SV Detection:** The framework begins with enhanced SV detection by integrating evidence from complementary callers (e.g., Manta, CNVnator, Lumpy) that leverage both breakpoint (BP) and read-depth (RD) signals[cite: 19, 54, 59, 60, 180].
2.  [cite_start]**Population Filtering:** A stringent filter removes common polymorphisms by leveraging allele frequencies from large-scale, multi-ancestral reference cohorts (1000 Genomes Project East Asian and GTEx European)[cite: 19, 117, 182]. [cite_start]Rare variants are defined as those with a combined background allele frequency < 0.01[cite: 119, 123, 298].
3.  [cite_start]**Transcriptome-Aware Prioritization:** The final and most crucial step uses the tissue-specific **Transcript Disruption Ratio (TDR)** to prioritize variants with the highest predicted functional impact[cite: 124, 320]. [cite_start]A TDR > 0.25 is used as the threshold to identify variants with a substantial biological effect[cite: 98, 124, 310].

[cite_start]This cascade dramatically reduces the number of candidates, leaving an average of only 2-5 high-confidence SVs per individual for final clinical interpretation[cite: 125].

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

### Step 1 (Optional): SV Merging and Pre-processing

Before annotation, it is common to merge SV callsets from multiple samples or callers.

**SV Merging with SURVIVOR**

We recommend using [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR) for merging.

```bash
# Create a file listing all VCF files to be merged
ls *.vcf > sample_files.txt

# Run SURVIVOR merge
# Parameters: sample_list, max_distance, min_callers_supporting, type_and_strand_consistency, min_length, output_vcf
SURVIVOR merge sample_files.txt 1000 1 1 1 0 50 merged_svs.vcf
