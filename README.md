# Single-cell Transcriptomic Analysis of Differentially Abundant Immune Cell Populations in PTCB and TCB

## Overview

This repository contains the computational workflow developed for the analysis of 
single-cell RNA sequencing (scRNA-seq) data from the study **GSE271413**. 

The project focuses on identifying and characterizing differentially abundant immune 
cell populations between **PTCB** and **TCB** groups using single-cell transcriptomic 
analysis approaches.

The workflow was implemented in **R** using the **Seurat** framework and includes 
quality control, normalization, dimensionality reduction, clustering, cell-type 
annotation, differential abundance analysis, differential gene expression analysis, 
and functional enrichment analysis.

---

## Research Objectives

The main objectives of this project were:

- Develop a reproducible single-cell RNA-seq analysis pipeline.
- Compare immune cell composition between PTCB and TCB groups.
- Identify cell populations with significant abundance differences.
- Characterize cell-type-specific transcriptional changes.
- Investigate biological pathways associated with immune regulation.

---

## Dataset

**Dataset:** GSE271413  
**Technology:** Single-cell RNA sequencing  
**Organism:** Human  
**Analysis focus:** Immune cell populations

The dataset contains multi-patient single-cell transcriptomic profiles that were 
integrated and analyzed to identify differences between experimental groups.

---

## Workflow Overview

### 1. Data Loading and Quality Control

Performed initial preprocessing steps including:

- Importing single-cell expression matrices
- Creating Seurat objects
- Evaluating quality metrics:
  - Number of detected genes per cell
  - Total transcript counts
  - Mitochondrial gene percentage
- Filtering low-quality cells

---

### 2. Normalization and Feature Selection

Applied Seurat-based preprocessing:

- Log-normalization of gene expression data
- Identification of highly variable genes
- Scaling of expression values
- Preparation of datasets for downstream analysis

---

### 3. Multi-patient Dataset Integration

To reduce patient-specific technical variation:

- Applied **Canonical Correlation Analysis (CCA)** integration
- Identified shared biological structures across samples
- Generated an integrated expression dataset for downstream analysis

---

### 4. Dimensionality Reduction and Clustering

Performed:

- Principal Component Analysis (PCA)
- UMAP dimensionality reduction
- Graph-based clustering

Generated low-dimensional representations to visualize cellular heterogeneity.

---

### 5. Cell-Type Annotation

Immune cell populations were identified using:

- Canonical marker gene expression
- Cluster-specific differential expression patterns
- Biological annotation of immune subpopulations

---

### 6. Differential Abundance Analysis

Compared immune cell composition between:

- PTCB group
- TCB group

Identified top differentially abundant immune cell populations and quantified 
changes in cellular distribution.

---

### 7. Differential Gene Expression Analysis

Performed cell-type-specific differential expression analysis to identify:

- Upregulated genes
- Downregulated genes
- Group-specific transcriptional signatures

---

### 8. Functional Enrichment Analysis

Investigated biological significance using:

- Gene Ontology (GO) enrichment analysis
- Immune-related pathway interpretation
- Functional characterization of differentially expressed genes

---

## Repository Structure
