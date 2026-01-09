#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import argparse
import pickle
from pathlib import Path
from mpi4py import MPI
import time
import logging


def setup_logging(rank):
    logging.basicConfig(
        level=logging.INFO,
        format=f'[Rank {rank}] %(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(f'synthetic_generation_rank_{rank}.log'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


def load_pgs_variants(pgs_file, max_variants=5e4, sort_variants_when_truncating=True):
    """ Chopped PGS loading - naive and inflexible - matches poorly and only supports PGS000729's formating """
    max_variants = int(max_variants)
    pgs = pd.read_csv(pgs_file, sep='\t', comment='#', header=None, low_memory=False)
    
    if len(pgs.columns) == 6: 
        pgs.columns = ['rsid', 'chr', 'pos', 'effect_allele', 'other_allele', 'effect_weight'] 
    else: 
        pgs.columns = ['rsid', 'chr', 'effect_allele', 'other_allele', 'effect_weight'] + [f'col_{i}' for i in range(5, len(pgs.columns))]
    
    pgs = pgs.dropna(subset=['effect_weight'])
    pgs['effect_weight'] = pd.to_numeric(pgs['effect_weight'], errors='coerce')
    pgs = pgs.dropna(subset=['effect_weight'])
    
    if len(pgs) > max_variants:
        if sort_variants_when_truncating:
            pgs = pgs.sort_values(by='effect_weight', ascending=False)
        pgs = pgs.iloc[:max_variants]
    
    return pgs[['rsid', 'chr', 'effect_allele', 'effect_weight']]


def distribute_samples(n_samples, comm_size, rank):
    """Distribute samples across MPI processes."""
    base_samples = n_samples // comm_size
    extra_samples = n_samples % comm_size
    
    if rank < extra_samples:
        local_samples = base_samples + 1
        start_idx = rank * (base_samples + 1)
    else:
        local_samples = base_samples
        start_idx = extra_samples * (base_samples + 1) + (rank - extra_samples) * base_samples
    
    return local_samples, start_idx


def generate_distributed_genotypes(
    variants, 
    n_samples=1e4, 
    h2=0.3, 
    maf_range=(0.05, 0.45),
    max_iterations=500,
    checkpoint_dir="checkpoints",
    output_dir="synthetic_data_distributed"
):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    logger = setup_logging(rank)
    
    n_samples = int(n_samples)
    n_variants = len(variants)
    weights = variants['effect_weight'].values
    
    local_samples, start_idx = distribute_samples(n_samples, size, rank)
    logger.info(f"Processing {local_samples} samples (indices {start_idx}-{start_idx + local_samples - 1})")
    
    np.random.seed(42 + rank)
    
    if rank == 0:
        mafs = np.random.uniform(maf_range[0], maf_range[1], n_variants)
        target_prs_global = np.random.normal(0, 1, n_samples)
    else:
        mafs = None
        target_prs_global = None
    
    mafs = comm.bcast(mafs, root=0)
    target_prs_global = comm.bcast(target_prs_global, root=0)
    
    target_prs_local = target_prs_global[start_idx:start_idx + local_samples]
    
    local_genotypes = np.random.binomial(
        2, mafs[:, np.newaxis], (n_variants, local_samples)
    ).T.astype(np.float32)
    
    if rank == 0:
        os.makedirs(checkpoint_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
    comm.barrier()
    
    start_time = time.time()
    
    for iteration in range(max_iterations):
        iter_start = time.time()
        
        local_prs = local_genotypes @ weights
        
        if iteration % 10 == 0 or iteration < 5:  # reduce communication overhead
            all_prs = comm.gather(local_prs, root=0)
            
            if rank == 0:
                global_prs = np.concatenate(all_prs)
                prs_corr = np.corrcoef(global_prs, target_prs_global)[0, 1] if np.std(global_prs) > 0 else 0
                prs_mse = np.mean((global_prs - target_prs_global)**2)
                logger.info(f"Iteration {iteration}: Corr={prs_corr:.4f}, MSE={prs_mse:.4f}")
            else:
                prs_corr = None
                prs_mse = None
            
            prs_corr = comm.bcast(prs_corr, root=0)
            prs_mse = comm.bcast(prs_mse, root=0)
        else:
            prs_corr = None
            prs_mse = None
        
        local_std_current = np.std(local_prs)
        local_std_target = np.std(target_prs_local)
        
        if local_std_target > 0 and local_std_current > 0:
            target_scaled = target_prs_local * local_std_current / local_std_target
        else:
            target_scaled = target_prs_local
        
        prs_diff = target_scaled - local_prs
        significant_mask = np.abs(weights) > 0.001
        
        weights_squared_sum = np.sum(weights**2)
        
        significant_indices = np.where(significant_mask)[0]
        
        for j in significant_indices:
            adjustment = prs_diff * weights[j] / weights_squared_sum
            prob_change = np.clip(adjustment * 0.05, -0.4, 0.4)
            
            random_vals = np.random.random(local_samples)
            
            # masking
            increase_mask = (
                (prob_change > 0) & 
                (local_genotypes[:, j] < 2) & 
                (random_vals < np.abs(prob_change))
            )
            decrease_mask = (
                (prob_change < 0) & 
                (local_genotypes[:, j] > 0) & 
                (random_vals < np.abs(prob_change))
            )
            
            # update updates
            local_genotypes[:, j] += increase_mask.astype(np.float32)
            local_genotypes[:, j] -= decrease_mask.astype(np.float32)
        
        # checkpoint
        if iteration % 50 == 0:
            checkpoint_data = {
                'iteration': iteration,
                'local_genotypes': local_genotypes,
                'local_samples': local_samples,
                'start_idx': start_idx,
                'rank': rank
            }
            
            checkpoint_file = Path(checkpoint_dir) / f"checkpoint_rank_{rank}_iter_{iteration}.pkl"
            with open(checkpoint_file, 'wb') as f:
                pickle.dump(checkpoint_data, f)
        
        # check for convergence
        if prs_corr is not None and prs_mse is not None:
            if abs(1.0 - prs_corr) < 1e-6 or prs_mse < 1e-6:
                logger.info(f"Converged at iteration {iteration}")
                break
        
        iter_time = time.time() - iter_start
        if rank == 0 and iteration % 50 == 0:
            logger.info(f"Iteration {iteration} completed in {iter_time:.3f}s")
    
    total_time = time.time() - start_time
    logger.info(f"Optimization completed in {total_time:.2f}s")
    
    final_local_prs = local_genotypes @ weights
    
    prs_var = np.var(final_local_prs)
    noise_var = prs_var * (1 - h2) / h2
    noise = np.random.normal(0, np.sqrt(noise_var), local_samples)
    local_phenotypes = final_local_prs + noise
    
    all_genotypes = comm.gather(local_genotypes, root=0)
    all_prs = comm.gather(final_local_prs, root=0)
    all_phenotypes = comm.gather(local_phenotypes, root=0)
    all_start_indices = comm.gather(start_idx, root=0)
    
    if rank == 0:
        logger.info("comb results from all")
        
        combined_genotypes = np.zeros((n_samples, n_variants), dtype=np.float32)
        combined_prs = np.zeros(n_samples)
        combined_phenotypes = np.zeros(n_samples)
        
        for i, (geno, prs, pheno, start) in enumerate(zip(all_genotypes, all_prs, all_phenotypes, all_start_indices)):
            end = start + len(prs)
            combined_genotypes[start:end] = geno
            combined_prs[start:end] = prs
            combined_phenotypes[start:end] = pheno
        
        geno_df = pd.DataFrame(
            combined_genotypes,
            columns=variants['rsid'].values,
            index=[f"Individual_{i}" for i in range(n_samples)]
        )
        
        # final global statistics
        final_corr = np.corrcoef(combined_prs, target_prs_global)[0, 1]
        final_mse = np.mean((combined_prs - target_prs_global)**2)
        
        logger.info(f"Final global metrics:")
        logger.info(f"PRS Correlation: {final_corr:.4f}")
        logger.info(f"PRS MSE: {final_mse:.4f}")
        
        return geno_df, mafs, combined_prs, combined_phenotypes
    else:
        return None, None, None, None


def save_distributed_results(geno_df, variants, prs_scores, phenotypes, mafs, output_dir):
    if geno_df is None:  # Only rank 0 has the combined results
        return
    
    os.makedirs(output_dir, exist_ok=True)
    
    geno_file = f"{output_dir}/synthetic_genotypes_distributed.tsv"
    geno_df.to_csv(geno_file, sep='\t')
    
    # Phenotypes
    variant_file = f"{output_dir}/variants.tsv"
    variants.to_csv(variant_file, sep='\t', index=False)
    
    pheno_df = pd.DataFrame({
        'Individual_ID': geno_df.index,
        'PRS': prs_scores,
        'Phenotype': phenotypes
    })
    pheno_file = f"{output_dir}/phenotypes_distributed.tsv"
    pheno_df.to_csv(pheno_file, sep='\t', index=False)
    
    # MAFs
    maf_df = pd.DataFrame({
        'rsid': variants['rsid'].values,
        'MAF': mafs
    })
    maf_file = f"{output_dir}/mafs.tsv"
    maf_df.to_csv(maf_file, sep='\t', index=False)
    
    return geno_file, variant_file, pheno_file, maf_file


def main():
    parser = argparse.ArgumentParser(description='SLURM-distributed synthetic genotype generation')
    parser.add_argument('--pgs_file', required=True, help='Path to PGS file')
    parser.add_argument('--n_samples', type=int, default=10000, help='Number of samples to generate')
    parser.add_argument('--max_variants', type=int, default=50000, help='Maximum number of variants')
    parser.add_argument('--h2', type=float, default=0.3, help='Heritability')
    parser.add_argument('--max_iterations', type=int, default=500, help='Maximum optimization iterations')
    parser.add_argument('--output_dir', default='synthetic_data_distributed', help='Output directory')
    parser.add_argument('--checkpoint_dir', default='checkpoints', help='Checkpoint directory')
    
    args = parser.parse_args()
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    if rank == 0:
        print(f"Starting distributed generation with {comm.Get_size()} processes")
        print(f"Loading variants from {args.pgs_file}")
        
        variants = load_pgs_variants(args.pgs_file, args.max_variants)
        print(f"Loaded {len(variants)} variants")
    else:
        variants = None
    
    variants = comm.bcast(variants, root=0)

    geno_df, mafs, prs_scores, phenotypes = generate_distributed_genotypes(
        variants,
        n_samples=args.n_samples,
        h2=args.h2,
        max_iterations=args.max_iterations,
        checkpoint_dir=args.checkpoint_dir,
        output_dir=args.output_dir
    )
    
    if rank == 0:
        save_distributed_results(geno_df, variants, prs_scores, phenotypes, mafs, args.output_dir)
        print(f"Results saved to {args.output_dir}")


if __name__ == "__main__":
    main()