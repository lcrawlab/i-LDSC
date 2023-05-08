#!/bin/bash
#SBATCH --job-name="gxe_dg"
#SBATCH --output=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/slurm_output/sout_'%J'.txt
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=6G
#SBATCH --chdir=/gpfs/data/sramacha/ukbiobank_jun17/ssmith/ongoing/meld/meld_sim/src/
#SBATCH --array=1-100
#SBATCH --account=ccmb-condo

# #########
# #  GxE  #
# #########
#hold other parameters constant and change the magnitude of amplification, generate with genetic pve of 0.6 and 0.3
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.1
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.2
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.3
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.4
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.5
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.6
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.7
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.8
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.9
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.6 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 2.0

# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.1
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.2
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.3
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.4
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.5
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.6
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.7
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.8
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.9
# srun time python TopDriver.py --exppath GxE --gxe --datagen -p 0.3 -r 0.5 --phi 0.25 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 2.0

# #hold other parameters constant and change the magnitude of amplification, generate with genetic pve of 0.6 and 0.3
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.1
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.2
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.3
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.4
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.5
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.6
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.7
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.8
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.9
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.6 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 2.0

# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.1
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.2
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.3
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.4
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.5
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.6
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.7
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.8
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 1.9
# srun time python TopDriver.py --exppath GxE_noCASS --gxe --datagen -p 0.3 -r 0.5 --phi 0.5 -w 10000 -s 0 --ss 5 --exp-alpha 1.0 --amplifier 2.0


#sim5test
srun time python TopDriver.py --exppath sim5 --gxe --datagen -p 0.6 -r 0.2 --phi 0.0 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.0
srun time python TopDriver.py --exppath sim5 --gxe --datagen -p 0.6 -r 0.4 --phi 0.0 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.0
srun time python TopDriver.py --exppath sim5 --gxe --datagen -p 0.6 -r 0.6 --phi 0.0 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.0
srun time python TopDriver.py --exppath sim5 --gxe --datagen -p 0.6 -r 0.8 --phi 0.0 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.0

srun time python TopDriver.py --exppath sim5 --gxe --datagen -p 0.3 -r 0.2 --phi 0.0 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.0
srun time python TopDriver.py --exppath sim5 --gxe --datagen -p 0.3 -r 0.4 --phi 0.0 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.0
srun time python TopDriver.py --exppath sim5 --gxe --datagen -p 0.3 -r 0.6 --phi 0.0 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.0
srun time python TopDriver.py --exppath sim5 --gxe --datagen -p 0.3 -r 0.8 --phi 0.0 -w 10000 -s 5 --ss 5 --exp-alpha 1.0 --amplifier 1.0





