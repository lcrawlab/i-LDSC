#!/bin/bash
#SBATCH --job-name="ld2_dg"
#SBATCH --output=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/slurm_output/sout_'%J'.txt
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --chdir=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/src/
#SBATCH --array=73-100
#SBATCH --account=ccmb-condo


# #########
# #  ld2 architecture  #
# #########
#hold other parameters constant and change the magnitude of amplification, generate with genetic pve of 0.6 and 0.3
#generate phenotypes that are determined only bny the SNPs in the highest x percentil of the ld score distribution
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.6 --percentage 1 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.6 --percentage 5 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.6 --percentage 10 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.6 --percentage 25 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.6 --percentage 50 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.6 --percentage 100 --exp-alpha 1

srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.3 --percentage 1 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.3 --percentage 5 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.3 --percentage 10 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.3 --percentage 25 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.3 --percentage 50 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample high --datagen -p 0.3 --percentage 100 --exp-alpha 1

#Now generate data where phenotypes are determined by SNPs in the the lower percentiles of the data
#don't repeat 100 because it is the same 
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.6 --percentage 1 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.6 --percentage 5 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.6 --percentage 10 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.6 --percentage 25 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.6 --percentage 50 --exp-alpha 1

srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.3 --percentage 1 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.3 --percentage 5 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.3 --percentage 10 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.3 --percentage 25 --exp-alpha 1
srun time python TopDriver.py --exppath ld2 --ld2based --l2sample low --datagen -p 0.3 --percentage 50 --exp-alpha 1



