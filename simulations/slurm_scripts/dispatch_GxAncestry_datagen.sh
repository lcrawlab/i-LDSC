#!/bin/bash
#SBATCH --job-name="ancestry_dg"
#SBATCH --output=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/slurm_output/sout_'%J'.txt
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --chdir=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/src/
#SBATCH --array=1-100
#SBATCH --account=ccmb-condo

# #############################
# #  GxAncestry architecture  #
# #############################
# # hold other parameters constant and change the magnitude of amplification, generate with genetic pve of 0.6 and 0.3
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 1 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 2 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 3 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 4 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 5 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 6 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 7 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 8 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 9 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 10 --exp-alpha 1

srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 1 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 2 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 3 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 4 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 5 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 6 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 7 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 8 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 9 --exp-alpha 1
srun time python TopDriver.py --exppath GxAncestry --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.25 -s 5 --ss 5 --pc 10 --exp-alpha 1

# #GxAncestry with no CASS interactions
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 1 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 2 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 3 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 4 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 5 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 6 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 7 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 8 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 9 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.6 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 10 --exp-alpha 1

# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 1 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 2 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 3 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 4 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 5 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 6 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 7 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 8 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 9 --exp-alpha 1
# srun time python TopDriver.py --exppath GxAncestry_noCASS --ancestry --datagen -p 0.3 -r 0.5 --kappa 0.5 -s 5 --ss 5 --pc 10 --exp-alpha 1
