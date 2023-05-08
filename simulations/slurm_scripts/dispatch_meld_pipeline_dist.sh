#!/bin/bash
#SBATCH --job-name="meld"
#SBATCH --output=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/slurm_output/sout_'%J'.txt
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=3G
#SBATCH --chdir=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/src/
#SBATCH --array=0-99
#SBATCH --account=ccmb-condo

###############################################################################
# Need to run dist_meld_jobs.py first to generate the filename2job mapping!!! #
# python dist_meld_jobs.py #directory                                         #
###############################################################################

srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID GxE 10 1
srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID GxE_noCASS 10 1

srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID ld2 10 1


srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID GxAncestry 10 1
srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID GxAncestry_noCASS 10 1

srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID sim1 0.0
srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID sim2 0.0
srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID sim3 0.5
srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID sim4 1.0
srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID sim5 1
srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID sim6 10 0.5
srun time python compute_ldscore_regr_dist_1kg.py $SLURM_ARRAY_TASK_ID sim7 1.0

