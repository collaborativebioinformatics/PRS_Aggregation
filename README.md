# PRS_Aggregation
Polygenic Risk Aggregation in common diseases and phenotypes

## Contributers

| Name | Email | ORCID | Institution |
|------|-------|-------|-------------|
| Ashok Kumar Sharma | ashoks773@gmail.com | https://orcid.org/0000-0002-2264-7628 | Cedars-Sinai Medical Center |
| Dmitriy Ivkov | Divkov@umich.edu | https://orcid.org/0009-0008-4536-3274 | — |
| Jasmine Baker | jasmine.baker@bcm.edu | https://orcid.org/0000-0001-7545-6086 | Baylor College of Medicine |
| Mengying Hu | meh251@pitt.edu | https://orcid.org/0000-0003-4827-3051 | — |
| Qianqian Liang | qil57@pitt.edu | https://orcid.org/0000-0002-1737-5031 | — |
| Shivank Sadasivan | ssadasiv@andrew.cmu.edu | https://orcid.org/0009-0004-4699-2129 | — |

## Introduction
Polygenic risk scores (PRS) are a powerful approach for predicting an individual’s likelihood of developing a phenotype or disease and have contributed substantially to our understanding of human health. PRS models are typically developed by multiple research groups, and aggregating PRS across ancestry groups represents an important opportunity to better leverage diverse datasets.
However, PRS data are often highly heterogeneous. For example, different PRS models may be based on distinct sets of SNPs, and SNP identifiers and formats can vary across resources. To address these challenges, we are developing a pipeline to harmonize PRS datasets across studies and generate a centralized PRS that integrates information from multiple sources.
In this work, we utilize the PRS Catalog, which contains PRS models from a wide range of resources. We also generated visualization tools to visualize the overlapping between different scores. We will also zoom into single-nucleotide resolution to visualize the surrounding SNPs. 

## Goals

## Description

## Overview Diagram
Data Sets: GWAS Catalog
-comes in TSV and HailMatrixTable option
-ideal to extract SNP weights for PRS 
PGS catalog- polygenic score catalog:https://www.pgscatalog.org/

<img width="617" height="435" alt="Screenshot 2026-01-07 at 12 46 30 PM" src="https://github.com/user-attachments/assets/72ca0ae5-263d-464d-93f4-5c79d3c4ddad" />

<img width="2624" height="1600" alt="Gemini_Generated_Image_7ft1ga7ft1ga7ft1" src="https://github.com/user-attachments/assets/4dfa70f3-2da5-45e3-a58e-907b69c4b147" />

## Methods
### PRS visualization tool (Mengying/Qianqian)
<img width="1068" height="561" alt="Screenshot 2026-01-08 at 11 41 16" src="https://github.com/user-attachments/assets/c0c668e8-7580-499f-a8ae-d0eac87d3d3b" />

<img width="790" height="460" alt="Upset_updated" src="https://github.com/user-attachments/assets/0ac3c48a-72d6-464f-a761-c5f64cd34f2c" />

## Federated Learning / Representation
<img width="1153" height="329" alt="Screenshot 2026-01-08 at 12 44 08 PM" src="https://github.com/user-attachments/assets/b0ce0f36-14d5-4184-b004-a5097ec1f406" />

<img width="1128" height="633" alt="Screenshot 2026-01-08 at 12 41 09 PM" src="https://github.com/user-attachments/assets/6b59ed8e-de1e-47aa-8bf6-941dc2c11170" />



## Results 
