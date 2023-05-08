#!/bin/bash
#SBATCH --job-name="melddg"
#SBATCH --output=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/slurm_output/sout_'%J'.txt
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --chdir=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/src/
#SBATCH --array=1-100


#########
# SIM 1 #
#########
#This first simulation command sets the percent of heritability explained by nonlinear interactions to 0.5, 
#the window size to 10000 kb, the sparsity to 1 percent, and the total pve to 0.6, there is no interaction
#between maf and effect size (exp-alpha)

#Ensuing commands toggle each of these parameters

srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.5 -w 10000 -s 1 -p 0.6 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.5 -w 10000 -s 5 -p 0.6 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.5 -w 10000 -s 10 -p 0.6 --exp-alpha 0.0

srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.8 -w 10000 -s 1 -p 0.6 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.6 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.8 -w 10000 -s 10 -p 0.6 --exp-alpha 0.0

srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.5 -w 10000 -s 1 -p 0.3 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.5 -w 10000 -s 5 -p 0.3 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.5 -w 10000 -s 10 -p 0.3 --exp-alpha 0.0

srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.8 -w 10000 -s 1 -p 0.3 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.3 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 0.8 -w 10000 -s 10 -p 0.3 --exp-alpha 0.0

srun time python TopDriver.py --exppath sim1/ --datagen --no-refdata -r 1.0 -w 10000 -s 1 -p 0.6 --exp-alpha 0.0

#########
# SIM 2 #
#########
#Iterates over same parameter sets as sim1 but with a larger window size of 100kb

srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.5 -w 100000 -s 1 -p 0.6 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.5 -w 100000 -s 5 -p 0.6 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.5 -w 100000 -s 10 -p 0.6 --exp-alpha 0.0

srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.8 -w 100000 -s 1 -p 0.6 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.8 -w 100000 -s 5 -p 0.6 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.8 -w 100000 -s 10 -p 0.6 --exp-alpha 0.0

srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.5 -w 100000 -s 1 -p 0.3 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.5 -w 100000 -s 5 -p 0.3 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.5 -w 100000 -s 10 -p 0.3 --exp-alpha 0.0

srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.8 -w 100000 -s 1 -p 0.3 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.8 -w 100000 -s 5 -p 0.3 --exp-alpha 0.0
srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 0.8 -w 100000 -s 10 -p 0.3 --exp-alpha 0.0

srun time python TopDriver.py --exppath sim2/ --datagen --no-refdata -r 1.0 -w 100000 -s 1 -p 0.6 --exp-alpha 0.0


# srun time python TopDriver.py --exppath GxE_alladd/ --datagen --no-refdata -r 1.0 -w 100000 -s 1 -p 0.3 --exp-alpha 0.0
# srun time python TopDriver.py --exppath GxE_alladd/ --datagen --no-refdata -r 1.0 -w 100000 -s 1 -p 0.6 --exp-alpha 0.0

#########
# SIM 3 #
#########
#Iterates over same parameter sets as sim1 but with a dependency between maf and effect size

srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.5 -w 10000 -s 1 -p 0.6 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.5 -w 10000 -s 5 -p 0.6 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.5 -w 10000 -s 10 -p 0.6 --exp-alpha 0.5

srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.8 -w 10000 -s 1 -p 0.6 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.6 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.8 -w 10000 -s 10 -p 0.6 --exp-alpha 0.5

srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.5 -w 10000 -s 1 -p 0.3 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.5 -w 10000 -s 5 -p 0.3 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.5 -w 10000 -s 10 -p 0.3 --exp-alpha 0.5

srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.8 -w 10000 -s 1 -p 0.3 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.3 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 0.8 -w 10000 -s 10 -p 0.3 --exp-alpha 0.5

srun time python TopDriver.py --exppath sim3/ --datagen --no-refdata -r 1.0 -w 10000 -s 1 -p 0.6 --exp-alpha 0.5

#########
# SIM 4 #
#########
#Iterates over same parameter sets as sim1 but with a total dependency between maf and effect size

srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.5 -w 10000 -s 1 -p 0.6 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.5 -w 10000 -s 5 -p 0.6 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.5 -w 10000 -s 10 -p 0.6 --exp-alpha 1.0

srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.8 -w 10000 -s 1 -p 0.6 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.6 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.8 -w 10000 -s 10 -p 0.6 --exp-alpha 1.0

srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.5 -w 10000 -s 1 -p 0.3 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.5 -w 10000 -s 5 -p 0.3 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.5 -w 10000 -s 10 -p 0.3 --exp-alpha 1.0

srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.8 -w 10000 -s 1 -p 0.3 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.3 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 0.8 -w 10000 -s 10 -p 0.3 --exp-alpha 1.0

srun time python TopDriver.py --exppath sim4/ --datagen --no-refdata -r 1.0 -w 10000 -s 1 -p 0.6 --exp-alpha 1.0


#########
# SIM 5 #
#########
#Testing across a range of rho paramters
# srun time python TopDriver.py --exppath sim5/ --datagen --no-refdata -r 0.2 -w 10000 -s 5 -p 0.6 --exp-alpha 0.0
# srun time python TopDriver.py --exppath sim5/ --datagen --no-refdata -r 0.4 -w 10000 -s 5 -p 0.6 --exp-alpha 0.0
# srun time python TopDriver.py --exppath sim5/ --datagen --no-refdata -r 0.6 -w 10000 -s 5 -p 0.6 --exp-alpha 0.0
# srun time python TopDriver.py --exppath sim5/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.6 --exp-alpha 0.0
# srun time python TopDriver.py --exppath sim5/ --datagen --no-refdata -r 0.2 -w 10000 -s 5 -p 0.3 --exp-alpha 0.0
# srun time python TopDriver.py --exppath sim5/ --datagen --no-refdata -r 0.4 -w 10000 -s 5 -p 0.3 --exp-alpha 0.0
# srun time python TopDriver.py --exppath sim5/ --datagen --no-refdata -r 0.6 -w 10000 -s 5 -p 0.3 --exp-alpha 0.0
# srun time python TopDriver.py --exppath sim5/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.3 --exp-alpha 0.0

# srun time python TopDriver.py --exppath sim5_alpha1/ --datagen --no-refdata -r 0.2 -w 10000 -s 5 -p 0.6 --exp-alpha 1.0
# srun time python TopDriver.py --exppath sim5_alpha1/ --datagen --no-refdata -r 0.4 -w 10000 -s 5 -p 0.6 --exp-alpha 1.0
# srun time python TopDriver.py --exppath sim5_alpha1/ --datagen --no-refdata -r 0.6 -w 10000 -s 5 -p 0.6 --exp-alpha 1.0
# srun time python TopDriver.py --exppath sim5_alpha1/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.6 --exp-alpha 1.0
# srun time python TopDriver.py --exppath sim5_alpha1/ --datagen --no-refdata -r 0.2 -w 10000 -s 5 -p 0.3 --exp-alpha 1.0
# srun time python TopDriver.py --exppath sim5_alpha1/ --datagen --no-refdata -r 0.4 -w 10000 -s 5 -p 0.3 --exp-alpha 1.0
# srun time python TopDriver.py --exppath sim5_alpha1/ --datagen --no-refdata -r 0.6 -w 10000 -s 5 -p 0.3 --exp-alpha 1.0
# srun time python TopDriver.py --exppath sim5_alpha1/ --datagen --no-refdata -r 0.8 -w 10000 -s 5 -p 0.3 --exp-alpha 1.0

# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.1

#########
# SIM 6 #
#########
#100kb window with Schoech's alpha of 0.5
srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.5 -w 100000 -s 1 -p 0.6 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.5 -w 100000 -s 5 -p 0.6 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.5 -w 100000 -s 10 -p 0.6 --exp-alpha 0.5

srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.8 -w 100000 -s 1 -p 0.6 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.8 -w 100000 -s 5 -p 0.6 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.8 -w 100000 -s 10 -p 0.6 --exp-alpha 0.5

srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.5 -w 100000 -s 1 -p 0.3 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.5 -w 100000 -s 5 -p 0.3 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.5 -w 100000 -s 10 -p 0.3 --exp-alpha 0.5

srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.8 -w 100000 -s 1 -p 0.3 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.8 -w 100000 -s 5 -p 0.3 --exp-alpha 0.5
srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 0.8 -w 100000 -s 10 -p 0.3 --exp-alpha 0.5

srun time python TopDriver.py --exppath sim6/ --datagen --no-refdata -r 1.0 -w 100000 -s 1 -p 0.6 --exp-alpha 0.5

#########
# SIM 7 #
#########
#Larger window full alpha dependency
srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.5 -w 100000 -s 1 -p 0.6 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.5 -w 100000 -s 5 -p 0.6 --exp-alpha 1.0 
srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.5 -w 100000 -s 10 -p 0.6 --exp-alpha 1.0

srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.8 -w 100000 -s 1 -p 0.6 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.8 -w 100000 -s 5 -p 0.6 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.8 -w 100000 -s 10 -p 0.6 --exp-alpha 1.0

srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.5 -w 100000 -s 1 -p 0.3 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.5 -w 100000 -s 5 -p 0.3 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.5 -w 100000 -s 10 -p 0.3 --exp-alpha 1.0

srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.8 -w 100000 -s 1 -p 0.3 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.8 -w 100000 -s 5 -p 0.3 --exp-alpha 1.0
srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 0.8 -w 100000 -s 10 -p 0.3 --exp-alpha 1.0

srun time python TopDriver.py --exppath sim7/ --datagen --no-refdata -r 1.0 -w 100000 -s 1 -p 0.6 --exp-alpha 1.0


