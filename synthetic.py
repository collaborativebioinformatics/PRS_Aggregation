"""
Naive synthetic data generation script - doesn't account for LD, marginal distributions

Method:
1) Random individual SNP intiialization
2) rudimentary Gradient Descent toward desired PRS weights

Parameters:
1) Variant (names)
2) Number of individuals
3) H2 Heritability
4) MAF range
"""


import pandas as pd
import numpy as np
import os
from tqdm import tqdm

def load_pgs_variants(pgs_file, max_variants=5e4, sort_variants_when_truncating=True):
    max_variants = int(max_variants)
    pgs = pd.read_csv(pgs_file, sep='\t', comment='#', header=None, low_memory=False)
    
    # this is the format in PGS000729.txt # this is probably a poor way to select for formatting
    if len(pgs.columns) == 6: 
        pgs.columns = ['rsid', 'chr', 'pos', 'effect_allele', 'other_allele', 'effect_weight'] 
    else: 
        pgs.columns = ['rsid', 'chr', 'effect_allele', 'other_allele', 'effect_weight'] + [f'col_{i}' for i in range(5, len(pgs.columns))]
    
    pgs = pgs.dropna(subset=['effect_weight'])
    pgs['effect_weight'] = pd.to_numeric(pgs['effect_weight'], errors='coerce')
    pgs = pgs.dropna(subset=['effect_weight'])
    
    if len(pgs) > max_variants:
        print(f"Subsetting to {max_variants} variants from {len(pgs)} total")
        if sort_variants_when_truncating:
            pgs = pgs.sort_values(by='effect_weight', ascending=False)
        pgs = pgs.iloc[:max_variants]
    
    return pgs[['rsid', 'chr', 'effect_allele', 'effect_weight']]

def generate_realistic_genotypes(variants, n_samples=1e4, h2=0.3, maf_range=(0.05, 0.45)):
    n_variants = len(variants)
    n_samples = int(n_samples)
    np.random.seed(42)
    
    weights = variants['effect_weight'].values
    mafs = np.random.uniform(maf_range[0], maf_range[1], n_variants)
    target_prs = np.random.normal(0, 1, n_samples)
    
    genotypes = np.random.binomial(2, mafs[:, np.newaxis], (n_variants, n_samples)).T.astype(np.float32)
    
    for iteration in tqdm(range(500)):
        current_prs = genotypes @ weights
        prs_corr = np.corrcoef(current_prs, target_prs)[0, 1] if np.std(current_prs) > 0 else 0
        prs_mse = np.mean((current_prs - target_prs)**2)
        
        print(f"Iteration {iteration}: Corr={prs_corr:.4f}, MSE={prs_mse:.4f}")
        
        target_scaled = target_prs * np.std(current_prs) / np.std(target_prs)
        prs_diff = target_scaled - current_prs # direction
        significant_mask = np.abs(weights) > 0.001 # identify only significant weights
        
        for j in range(n_variants):
            if significant_mask[j]:
                adjustment = prs_diff * weights[j] / np.sum(weights**2)
                prob_change = np.clip(adjustment * 0.05, -0.4, 0.4)
                
                random_vals = np.random.random(n_samples)
                
                # increase/decrease mask (where needed)
                increase_mask = (prob_change > 0) & (genotypes[:, j] < 2) & (random_vals < np.abs(prob_change))
                genotypes[increase_mask, j] += 1
                decrease_mask = (prob_change < 0) & (genotypes[:, j] > 0) & (random_vals < np.abs(prob_change))
                genotypes[decrease_mask, j] -= 1
    
    final_prs = genotypes @ weights
    
    final_corr = np.corrcoef(final_prs, target_prs)[0, 1]
    final_mse = np.mean((final_prs - target_prs)**2)
    
    print(f"\nFinal:")
    print(f"PRS Correlation: {final_corr:.4f}")
    print(f"PRS MSE: {final_mse:.4f}")
    
    prs_var = np.var(final_prs)
    noise_var = prs_var * (1 - h2) / h2
    noise = np.random.normal(0, np.sqrt(noise_var), n_samples)
    phenotypes = final_prs + noise

    geno_df = pd.DataFrame(genotypes, 
                          columns=variants['rsid'].values,
                          index=[f"Individual_{i}" for i in range(n_samples)])
    
    return geno_df, mafs, final_prs, phenotypes

def save_genotype_data(geno_df, variants, prs_scores=None, phenotypes=None, output_dir="synthetic_data"):
    os.makedirs(output_dir, exist_ok=True)
    
    geno_file = f"{output_dir}/synthetic_genotypes.tsv"
    geno_df.to_csv(geno_file, sep='\t')
    
    variant_file = f"{output_dir}/variants.tsv"
    variants.to_csv(variant_file, sep='\t', index=False)
    
    if prs_scores is not None and phenotypes is not None:
        pheno_df = pd.DataFrame({
            'Individual_ID': geno_df.index,
            'PRS': prs_scores,
            'Phenotype': phenotypes
        })
        pheno_file = f"{output_dir}/phenotypes.tsv"
        pheno_df.to_csv(pheno_file, sep='\t', index=False)
    
    return geno_file, variant_file

if __name__ == "__main__":
    pgs_file = "data/PGS000729.txt"
    variants = load_pgs_variants(pgs_file)
    
    genotypes, mafs, prs_scores, phenotypes = generate_realistic_genotypes(variants, n_samples=1000)
    
    geno_file, variant_file = save_genotype_data(genotypes, variants, prs_scores, phenotypes)
    
    
    
