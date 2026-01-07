# PRS_Aggregation
Polygenic Risk Aggregation in common diseases and phenotypes

Data Sets: GWAS Catalog
-comes in TSV and HailMatrixTable option
-ideal to extract SNP weights for PRS 
Synthetic data from https://biobank.ndph.ox.ac.uk/synthetic_dataset/ 
PGS catalog- polygenic score catalog:https://www.pgscatalog.org/

Goals 

Working Prototype 
Aggregate PRS across biobanks
Tag based on attributes 
Download phenotype TSV 
Variant harmonization & SNP selection
Quality control and dropping ambiguous data
Choose top # of SNPs by pvalue 
create /download synthetic data set
Simulate genotypes 
Simulate phenotype of choice + covariates and account for noise in the data set
Model & Evaluate 
Computation
Train Model: Logistic Regression 
Example: height ~ PRS + age + sex 
Ancestry specific PRS ? (population structure artifact or real thing)
NVFlare
Aggregate model parameters 
Iterate until convergence
Visualize 
Reproducibility
Github repository 


<img width="591" height="363" alt="PRS" src="https://github.com/user-attachments/assets/c4dc46fc-3ee3-40ae-b84b-5a5eacf4655a" />
