#!/bin/bash
#SBATCH --job-name=synthetic_genotype_large
#SBATCH --output=synthetic_large_%j.out
#SBATCH --error=synthetic_large_%j.err
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=08:00:00
#SBATCH --partition=compute

# For very large-scale generation (1M+ samples)

module load python/3.9
module load openmpi/4.1.0

export OMPI_MCA_btl_vader_single_copy_mechanism=none
export OMP_NUM_THREADS=1

# Enable memory mapping for large datasets
export OMPI_MCA_mpi_leave_pinned=1
export OMPI_MCA_mpi_leave_pinned_pipeline=1

mpirun -np $SLURM_NTASKS python synthetic_slurm.py \
    --pgs_file data/PGS000729.txt \
    --n_samples 1000000 \
    --max_variants 100000 \
    --h2 0.3 \
    --max_iterations 1000 \
    --output_dir synthetic_data_large \
    --checkpoint_dir checkpoints_large

echo "Large-scale job completed at $(date)"