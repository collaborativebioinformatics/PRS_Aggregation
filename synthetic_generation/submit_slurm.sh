#!/bin/bash
#SBATCH --job-name=synthetic_genotype_generation
#SBATCH --output=synthetic_generation_%j.out
#SBATCH --error=synthetic_generation_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --partition=compute

# Load required modules (adjust based on your cluster)
module load python/3.9
module load openmpi/4.1.0

# Activate virtual environment if needed
# source venv/bin/activate

# Set OpenMPI parameters for optimal performance
export OMPI_MCA_btl_vader_single_copy_mechanism=none
export OMP_NUM_THREADS=1

# Run the distributed synthetic data generation
mpirun -np $SLURM_NTASKS python synthetic_slurm.py \
    --pgs_file data/PGS000729.txt \
    --n_samples 100000 \
    --max_variants 50000 \
    --h2 0.3 \
    --max_iterations 500 \
    --output_dir synthetic_data_distributed \
    --checkpoint_dir checkpoints

echo "Job completed at $(date)"