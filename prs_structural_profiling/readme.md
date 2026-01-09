# PRS Summarization & Structural Comparison

This module provides a reproducible pipeline to **summarize and structurally compare multiple polygenic risk scores (PRS)** prior to aggregation.

The goal of this track is to answer a fundamental question:

> **Before aggregating PRS, how similar are they—at the SNP and gene levels?**

---

## Why Summarization Comes First

Multiple PRS often exist for the same phenotype, developed by different studies using different cohorts and methods.  
Although these PRS are often treated as interchangeable, their **underlying variant composition can differ substantially**.

This module focuses on **structural comparison**, not predictive performance, by answering:

- How many SNPs does each PRS contain?
- How much SNP overlap exists across PRS?
- How many SNPs fall within known genes?
- How much gene-level overlap exists when considering regulatory flanks?

These summaries provide **critical context** for downstream PRS selection and aggregation.

---

## What This Module Does

Given two or more harmonized PRS files, this pipeline:

1. Standardizes SNP identifiers using genomic coordinates  
2. Computes SNP- and gene-level summary statistics per PRS  
3. Maps SNPs to genes using gene bodies and optional flanking windows  
4. Visualizes overlap across PRS using UpSet plots  
5. Generates publication-ready tables and figures  

---

## Input Requirements

This module assumes **harmonized PRS inputs**, such as those from the PGS Catalog.

Each input file **must include** the following columns:

- `hm_chr` — chromosome (numeric, no `chr` prefix)
- `hm_pos` — base-pair position (GRCh37 / hg19)
  
Additional columns may be present but are not required for summarization.

> ⚠️ All PRS compared together must use the **same genome build**.

---

## How It Works (High-Level)

The pipeline performs three levels of summarization:

### 1. SNP-Level
- Total number of SNPs per PRS
- SNP overlap across PRS
- Chromosome-level SNP distribution

### 2. Gene-Level
- SNPs overlapping known gene bodies
- SNPs overlapping genes with user-defined flanking windows
- Gene overlap across PRS

### 3. Cross-PRS Comparison
- Presence/absence matrices
- UpSet visualizations for SNP and gene overlap

---

## Running the Pipeline

This module is implemented as a command-line R script.

### Example Command

```bash
Rscript Rscript_PRSAggregator_Summarization.R \
  --files PGS000020_hmPOS_GRCh37.txt.gz,PGS000804_hmPOS_GRCh37.txt.gz,PGS001818_hmPOS_GRCh37.txt.gz \
  --out results \
  --flank 50000
