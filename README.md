# Breast-Cancer-with-Liver-Metastasis
This repository contains the data and codes necessary for analysis of the scRNA-seq data of human breast cancer with liver metastasis presented in the artical titled "•	Single-cell transcriptomics reveals cellular dynamics, immune reprogramming, and metabolic adaptations during breast cancer liver metastasis".

## Data Processing
Raw sequencing data were processed using the Cell Ranger software provided by 10x Genomics to perform sample demultiplexing, barcode processing, and single-cell gene expression quantification(31). The reads were aligned to the human reference genome (GRCh38) using STAR within the Cell Ranger pipeline, and a gene-barcode matrix was generated for each sample. scTCR-seq -seq data were preprocessed by CellRanger for V(D)J sequence assembly and TCR reconstruction using the human reference genome build 38. Only high-confidence and productive TRA/TRB annotations were used for further analysis. Cells with fewer than 200 detected genes or with >20% mitochondrial gene expression were filtered out to exclude low-quality cells and potential doublets.

## Data Downloading
Raw fastq data: GSA-Human Bioproject ID: PRJCA039607
Processed Seurat object preview link: https://data.mendeley.com/preview/s27yn9ykbv?a=70f52439-8947-4f3a-8a2a-cdf9af1d444c

## Data visualization
### Requirements
1. R version 4.3.1
2. R packages
      Seurat
      ggplot2
      data.table
      dplyr
      tidyr
      ArchR
      pheatmap
      ggpubr
