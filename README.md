# PRS_Aggregation
Polygenic Risk Score (PRS) Aggregation in common diseases and phenotypes across databases

## Overview Diagram

<img width="2624" height="1600" alt="Gemini_Generated_Image_7ft1ga7ft1ga7ft1" src="https://github.com/user-attachments/assets/4dfa70f3-2da5-45e3-a58e-907b69c4b147" />

## Introduction
Polygenic risk scores (PRS) are a powerful approach for predicting an individual’s likelihood of developing a phenotype or disease and have contributed substantially to our understanding of human health. PRS models are typically developed by multiple research groups, and aggregating PRS across ancestry groups represents an important opportunity to better leverage diverse datasets.
However, PRS data are often highly heterogeneous. For example, different PRS models may be based on distinct sets of SNPs, and SNP identifiers and formats can vary across resources. To address these challenges, we are developing a pipeline to harmonize PRS datasets across studies and generate a centralized PRS that integrates information from multiple sources.
In this work, we utilize the PRS Catalog, which contains PRS models from a wide range of resources. We also generated visualization tools to visualize the overlapping between different scores. We will also zoom into single-nucleotide resolution to visualize the surrounding SNPs. 

## Goals
The project aims to create a practical, reproducible framework for aggregating polygenic risk scores (PRS) across databases, studies, and ancestry groups. 

Specifically, we will:

Aggregate PRS models across resources to produce a centralized, integrated PRS representation for common diseases and phenotypes.

Harmonize heterogeneous PRS inputs (differing SNP sets, identifier conventions, file formats, and metadata schemas) into a consistent, comparable structure.

Enable comparability and quality control, including standardized SNP identifiers, allele alignment checks, coordinate normalization, and consistent weight representation.

Quantify and visualize overlap among PRS models, supporting assessment of shared vs. unique variants and model similarity.

Lay groundwork for scalable aggregation strategies, specifically exploration of federated learning/representation approaches that can help integrate signals across cohorts or sites without compromising sensitive data.

## Description
PRS are widely used to estimate an individual’s genetic predisposition to phenotypes and common diseases. However, PRS models are frequently developed by independent research groups and released through different resources (public and private), leading to significant complications in reusing and aggregating data. Complications typically include: 1) Non-overlapping or partially overlapping SNP sets 2)Inconsistent SNP identifiers (rsIDs vs. chrom:pos vs. internal IDs) 3) Variable formatting (TSV/flat files vs. Hail MatrixTable representations) 4) Divergent allele encoding and genome build conventions 5) Uneven metadata completeness and inconsistent trait definitions. We address these challenges by developing a pipeline that harmonizes PRS datasets across studies and enables the construction of a centralized PRS integrating information from multiple sources. 

The initial implementation focuses on leveraging PGS Catalog (Polygenic Score Catalog) as a primary source of published PRS models and accompanying metadata and 
GWAS Catalog as a complementary resource. In addition to harmonization and aggregation, the repository includes visual analytics to help users understand how different PRS relate at the model level and at the single-variant level. The overarching outcome is a standardized approach to aggregate different PRS across databases in a way that is transparent, extensible, and suitable for downstream analyses across diverse populations.

## Methods
### PRS visualization tool (Mengying/Qianqian)
<img width="1068" height="561" alt="Screenshot 2026-01-08 at 11 41 16" src="https://github.com/user-attachments/assets/c0c668e8-7580-499f-a8ae-d0eac87d3d3b" />

<img width="790" height="460" alt="Upset_updated" src="https://github.com/user-attachments/assets/0ac3c48a-72d6-464f-a761-c5f64cd34f2c" />

## Federated Learning / Representation
<img width="1153" height="329" alt="Screenshot 2026-01-08 at 12 44 08 PM" src="https://github.com/user-attachments/assets/b0ce0f36-14d5-4184-b004-a5097ec1f406" />

<img width="1128" height="633" alt="Screenshot 2026-01-08 at 12 41 09 PM" src="https://github.com/user-attachments/assets/6b59ed8e-de1e-47aa-8bf6-941dc2c11170" />



## Results 









## Contributers

| Name | Email | ORCID | Institution |
|------|-------|-------|-------------|
| Ashok K. Sharma | ashoks773@gmail.com | https://orcid.org/0000-0002-2264-7628 | Cedars-Sinai Medical Center, LA |
| Dmitriy Ivkov | Divkov@umich.edu | https://orcid.org/0009-0008-4536-3274 | University of Michigan, Michigan |
| Jasmine Baker | jasmine.baker@bcm.edu | https://orcid.org/0000-0001-7545-6086 | Baylor College of Medicine, Houston |
| Mengying Hu | meh251@pitt.edu | https://orcid.org/0000-0003-4827-3051 | — |
| Qianqian Liang | qil57@pitt.edu | https://orcid.org/0000-0002-1737-5031 | Population Health Sciences, Danville |
| Shivank Sadasivan | ssadasiv@andrew.cmu.edu | https://orcid.org/0009-0004-4699-2129 | Carnegie Mellon University, Pittsburgh |

