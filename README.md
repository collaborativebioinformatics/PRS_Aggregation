# PRSAggregator
Visual Exploration and Structural Comparison of Polygenic Risk Scores for Aggregation

## Overview Diagram

<img width="1271" height="768" alt="Screenshot 2026-01-09 at 1 13 52 PM" src="https://github.com/user-attachments/assets/51fc26c9-5789-4493-8605-cdcbbd0546bb" />




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

## 1 PRS Profiling


### Motivation
### Methods
### Results
<img width="790" height="460" alt="Upset_updated" src="https://github.com/user-attachments/assets/0ac3c48a-72d6-464f-a761-c5f64cd34f2c" />

## 2 PRS Locus Viewer


### Motivation
### Methods
### Results

<img width="1068" height="561" alt="Screenshot 2026-01-08 at 11 41 16" src="https://github.com/user-attachments/assets/c0c668e8-7580-499f-a8ae-d0eac87d3d3b" />

## 3 Federated Representation Learning for Polygenic Risk Scores
### Motivation

A **representation learning framework for Polygenic Risk Scores (PGS)** that enables systematic analysis of genetic architecture across ancestries and studies, without requiring individual-level genotype data.

The framework:

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

- **Cosine distance heatmaps**  
  *PGSs derived from similar ancestries show higher similarity in genetic risk profiles.*

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
- **Variant PCA / UMAP (colored by ancestry sharing)**  
  *Variants shared across multiple ancestries occupy distinct regions, while ancestry-specific variants cluster tightly.*

- **Sharedness stratification (shared_2, shared_3, …)**  
  *Widely shared variants show structured dispersion, suggesting consistent but nuanced cross-population effects.*

---

### Biological Context
- **Coloring by gene region / mutation type**  
  *Embedding structure aligns with biological annotation, supporting interpretability.*

---
## Federated Learning Perspective

- Each **PGS = client**
- Each **trait can have a different model (AE)**
- Each **variant = model parameter**
- **Shared variants** behave like global parameters
- **Ancestry-specific variants** act as client-private signals

This framework enables **privacy-preserving, interpretable analysis of genetic architecture** across populations without pooling individual-level data.

---



<img width="1153" height="329" alt="Screenshot 2026-01-08 at 12 44 08 PM" src="https://github.com/user-attachments/assets/b0ce0f36-14d5-4184-b004-a5097ec1f406" />

<img width="1128" height="633" alt="Screenshot 2026-01-08 at 12 41 09 PM" src="https://github.com/user-attachments/assets/6b59ed8e-de1e-47aa-8bf6-941dc2c11170" />



## Future directions





## Contributers

| Name | Email | ORCID | Institution |
|------|-------|-------|-------------|
| Ashok K. Sharma | ashoks773@gmail.com | https://orcid.org/0000-0002-2264-7628 | Cedars-Sinai Medical Center, LA |
| Dmitriy Ivkov | Divkov@umich.edu | https://orcid.org/0009-0008-4536-3274 | University of Michigan, Michigan |
| Jasmine Baker | jasmine.baker@bcm.edu | https://orcid.org/0000-0001-7545-6086 | Baylor College of Medicine, Houston |
| Mengying Hu | meh251@pitt.edu | https://orcid.org/0000-0003-4827-3051 | University of Pittsburgh, Pittsburgh |
| Qianqian Liang | qil57@pitt.edu | https://orcid.org/0000-0002-1737-5031 | Population Health Sciences, Geisinger, Danville, PA |
| Shivank Sadasivan | ssadasiv@andrew.cmu.edu | https://orcid.org/0009-0004-4699-2129 | Carnegie Mellon University, Pittsburgh |

