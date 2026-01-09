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
- **PRS Federated Presentations**: Establish an representation tool to visualize PRS scores across ancestry and different locations.

## PRS Profiling


### Motivation
### Methods
### Results
<img width="790" height="460" alt="Upset_updated" src="https://github.com/user-attachments/assets/0ac3c48a-72d6-464f-a761-c5f64cd34f2c" />

## PRS Locus Viewer


### Motivation
### Methods
### Results

<img width="1068" height="561" alt="Screenshot 2026-01-08 at 11 41 16" src="https://github.com/user-attachments/assets/c0c668e8-7580-499f-a8ae-d0eac87d3d3b" />

## Federated Learning/Representation
### Motivation
### Methods
### Results

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

