# PRSAggregator
Visual Exploration and Structural Comparison of Polygenic Risk Scores for Aggregation

## Overview Diagram
<img width="1155" height="702" alt="Screenshot 2026-01-09 at 1 31 51 PM" src="https://github.com/user-attachments/assets/cd6e3d50-b70d-43fe-9956-68d94d4b41de" />




## Background
Polygenic Risk Scores (PRS) are widely used to estimate genetic susceptibility to complex diseases. For many common traits and diseases, multiple PRS have been developed by different studies using diverse cohorts, methodologies, and SNP selection strategies.

Although the PGS Catalog provides harmonized PRS data, researchers still lack practical tools to compare multiple PRS at the structural level and to understand how these scores relate to one another before downstream use.

## Motivation

Aggregating multiple PRS has the potential to improve robustness and generalizability. However, PRS aggregation is challenging because:

- Different PRS often use partially overlapping but non-identical SNP sets  
- Redundancy and complementarity between PRS are unclear  
- PRS selection is often arbitrary and poorly justified  

**Before aggregating PRS, it is essential to understand how they overlap and differ.**

## What This Project Does

This project provides a framework to **summarize, visualize, and explore overlap among multiple PRS** using harmonized data from the PGS Catalog.

Specifically,our project contains three main part:

- **PRS Profiling**: Build a wrapper pipeline to summarize SNP- and gene-level information across multiple PRS. And visualize using Upset Plot. 
- **PRS Locus Viewer**: Enable an interactive tool to explore of SNPs and genes in genomic context
- **PRS Federated Representation**: Establish an representation learning approach for PRS scores across ancestry and trait from different studies.

## 1. PRS Profiling
Detailed information locates under `/PRS_Structural_Profiling/`

### Motivation

This module provides a reproducible pipeline to **summarize and structurally compare multiple polygenic risk scores (PRS)** prior to aggregation.

The goal of this track is to answer a fundamental question:

> **Before aggregating PRS, how similar are they—at the SNP and gene levels?**



### Methods

#### input file format

PRS score files downloaded from the **PGS Catalog** [https://www.pgscatalog.org], with Genome build: **GRCh37**


Each input file **must include** the following columns:

- `hm_chr` — chromosome (numeric, no `chr` prefix)
- `hm_pos` — base-pair position (GRCh37 / hg19)

#### Example Command to run the pipeline

```bash
Rscript Rscript_PRSAggregator_Summarization.R \
  --files PGS000020_hmPOS_GRCh37.txt.gz,PGS000804_hmPOS_GRCh37.txt.gz,PGS001818_hmPOS_GRCh37.txt.gz \
  --out results \
  --flank 50000

```

### Results
We used Type 2 Diabetes as an example and we took 3 PGS scores overall. 

#### Structural summary of PRSs for Type 2 Diabetes (T2D)

| PRS ID | # SNPs | SNPs within genes | # Genes (genic) | SNPs within ±50kb | # Genes (±50kb) |
|------|-------:|------------------:|----------------:|------------------:|----------------:|
| PGS000020 | 7,502 | 3,541 | 2,735 | 5,227 | 7,624 |
| PGS000804 | 578 | 342 | 366 | 496 | 1,183 |
| PGS001818 | 30,745 | 14,137 | 5,084 | 21,379 | 12,805 |
 
#### UpSet plot (SNP overlap)

<img width="1333" height="933" alt="upset_snp" src="https://github.com/user-attachments/assets/f27e6751-4a34-46e0-bad2-523e92a18b45" />

#### UpSet plot (Gene overlap)

<img width="1333" height="933" alt="upset_gene" src="https://github.com/user-attachments/assets/05f53a89-5a9a-4a35-9bd8-c2d131563106" />


## 2. PRS Locus Viewer


### Motivation
### Methods
### Results

<img width="1068" height="561" alt="Screenshot 2026-01-08 at 11 41 16" src="https://github.com/user-attachments/assets/c0c668e8-7580-499f-a8ae-d0eac87d3d3b" />

## 3. Federated Representation Learning for Polygenic Risk Scores
### Motivation

A **representation learning framework for Polygenic Risk Scores (PGS)** that enables systematic analysis of genetic architecture across ancestries and studies, without requiring individual-level genotype data.

The framework:

<img width="1059" height="460" alt="Screenshot 2026-01-09 at 2 03 20 PM" src="https://github.com/user-attachments/assets/e2681e2d-1d03-442b-b996-44549c453bb0" />


- Integrates **heterogeneous PGS scoring files** across ancestries and studies  
- Harmonizes variants at the **locus (genomic position) level**  
- Enriches variants with **biological annotations** (genes, gene regions, mutation types)  

It learns **interpretable embeddings at two levels**:

- **PGS-level embeddings** → cohort / ancestry representations  
- **Variant-level embeddings** → locus representations  

These embeddings are interpreted through the lens of:

- **Ancestry**
- **Disease / trait**
- **Cross-ancestry sharing of variants**

The entire pipeline naturally aligns with a **federated learning perspective**, where each PGS acts as a client and shared vs ancestry-specific signals emerge from the learned representations.

---
### Methods and Results
### 1. Data Harmonization & Feature Construction

<img width="998" height="442" alt="Screenshot 2026-01-09 at 1 42 27 PM" src="https://github.com/user-attachments/assets/91a72a23-3dd0-42a9-9cf2-a494c020bd59" />


- Parsed PGS Catalog scoring files across multiple ancestries
- Unified heterogeneous formats (different genome builds, weight types)
- Defined a canonical **locus identifier** using harmonized coordinates: locus_id = hm_chr : hm_pos
- Normalized effect sizes into a common numeric space (`weight_scaled`)
- Annotated each locus using Ensembl VEP: (genomic information)

### 2. Autoencoder I — PGS-Level Representation Learning

**Objective:**  
Learn compact embeddings that capture how each PGS distributes genetic risk across loci.

**Input:**  
- Sparse matrix: **PGS × loci**
- Values = normalized effect sizes

**Architecture (PGS Autoencoder):**
- Encoder:
- Dense projection (sparse → latent)
- Nonlinear activation
- Decoder:
- Reconstructs original PGS risk profile
- Loss:
- Masked reconstruction loss (emphasizes non-zero effects)
- Evaluation:
- Train/test split
- Cosine similarity in input space

**Output:**  
- One embedding per PGS representing:
- ancestry-specific genetic architecture
- disease-related signal
- similarity to other cohorts

### PGS-Level Embeddings
- **PCA / UMAP plots**  
  *PGS embeddings cluster by ancestry even for the same disease, highlighting population-specific genetic architectures.*

<img width="892" height="691" alt="Screenshot 2026-01-09 at 1 43 22 PM" src="https://github.com/user-attachments/assets/29443c34-77de-4c20-a443-6191d3916ae0" />

- **Cosine distance heatmaps**  
  *PGSs derived from similar ancestries show higher similarity in genetic risk profiles.*

<img width="593" height="496" alt="Screenshot 2026-01-09 at 1 43 49 PM" src="https://github.com/user-attachments/assets/0bc960fd-2034-4148-9906-2cd82d902f10" />


---


### 3. Autoencoder II — Variant-Level Representation Learning

**Key Idea:**  
Transpose the problem and treat each variant as a vector across PGSs.

Variant = [effect in PGS₁, effect in PGS₂, …]

**Objective:**  
Learn embeddings that capture how loci behave across ancestries and studies.

**Architecture (Variant Autoencoder):**
- Input:
  - Variant × PGS effect vectors
- Encoder:
  - Low-dimensional latent space (compact locus representation)
- Decoder:
  - Reconstructs variant effect profile
- Filtering:
  - Optional minimum PGS support (focus on shared variants)

**Output:**  
- One embedding per variant encoding:
  - consistency vs heterogeneity across ancestries
  - shared vs ancestry-specific behavior

### Variant-Level Embeddings
- **Variant PCA / UMAP **  
<img width="947" height="761" alt="Screenshot 2026-01-09 at 1 45 22 PM" src="https://github.com/user-attachments/assets/90c117ca-87b4-4fab-a354-069d6935e98f" />

<img width="947" height="761" alt="Screenshot 2026-01-09 at 1 45 37 PM" src="https://github.com/user-attachments/assets/6265e1b1-1add-49b5-b8d5-d4b094cc20ba" />

<img width="947" height="761" alt="Screenshot 2026-01-09 at 1 45 56 PM" src="https://github.com/user-attachments/assets/d82385a3-7139-4966-b0c0-0cd9ea41a8e2" />

- **Sharedness stratification (shared_2, shared_3, …)**  
  *Widely shared variants show structured dispersion, suggesting consistent but nuanced cross-population effects.*
<img width="947" height="761" alt="Screenshot 2026-01-09 at 1 46 23 PM" src="https://github.com/user-attachments/assets/dc4c0fa2-e6d4-44af-9924-2ebd9cc21039" />

---

## Federated Learning Perspective

- Each **PGS = client**
- Each **trait can have a different model (AE)**
- Each **variant = model parameter**
- **Shared variants** behave like global parameters
- **Ancestry-specific variants** act as client-private signals

This framework enables **privacy-preserving, interpretable analysis of genetic architecture** across populations without pooling individual-level data.

---


## Future directions

## Initial Workflow
<img width="1128" height="633" alt="Screenshot 2026-01-08 at 12 41 09 PM" src="https://github.com/user-attachments/assets/6b59ed8e-de1e-47aa-8bf6-941dc2c11170" />









## Contributers

| Name | Email | ORCID | Institution |
|------|-------|-------|-------------|
| Ashok K. Sharma | ashoks773@gmail.com | https://orcid.org/0000-0002-2264-7628 | Cedars-Sinai Medical Center, LA |
| Dmitriy Ivkov | Divkov@umich.edu | https://orcid.org/0009-0008-4536-3274 | University of Michigan, Michigan |
| Jasmine Baker | jasmine.baker@bcm.edu | https://orcid.org/0000-0001-7545-6086 | Baylor College of Medicine, Houston |
| Mengying Hu | meh251@pitt.edu | https://orcid.org/0000-0003-4827-3051 | University of Pittsburgh, Pittsburgh |
| Qianqian Liang | qil57@pitt.edu | https://orcid.org/0000-0002-1737-5031 | Population Health Sciences, Geisinger, Danville, PA |
| Shivank Sadasivan | ssadasiv@andrew.cmu.edu | https://orcid.org/0009-0004-4699-2129 | Carnegie Mellon University, Pittsburgh |

