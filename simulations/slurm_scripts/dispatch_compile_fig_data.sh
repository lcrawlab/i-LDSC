#!/bin/bash
#SBATCH --job-name="meld"
#SBATCH --output=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/slurm_output/sout_'%J'.txt
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=3G
#SBATCH --chdir=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/src/
#SBATCH --account=ccmb-condo

# python compile.figure.data.py GxE
# python compile.figure.data.py GxE_noCASS
python compile.figure.data.py GxE_sc

# python dist_meld_jobs.py sim1
# python dist_meld_jobs.py sim2
# python dist_meld_jobs.py sim3
# python dist_meld_jobs.py sim4
# python dist_meld_jobs.py sim5
# python dist_meld_jobs.py sim6
# python dist_meld_jobs.py sim7
