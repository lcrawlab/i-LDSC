# Marginal Epistatic LD score regression (MELD)

The inflation of test statistics in genome-wide association (GWA) studies due to confounding such as cryptic relatedness, population stratification, and spurious non-zero genetic effects driven by linkage disequilibrium (LD) has been well characterized in the literature. This repository holds **MELD**: a software package to run **marginal epistatic LD score regression**. MELD is an extended framework which takes in GWA test statistics and accurately partitions true additive genetic variation from confounding non-additive genetic variation, as well as other biases.

## Installation Requirements

We recommend creating an anaconda environment for running MELD, instructions for setting up conda environments can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

## Calculation of MELD Scores (Optional)

Calculation of MELD scores requires a working version of pip and Python 3.9.12 to be installed.

To install the necessary python packages, navigate to the woking directory and submit:

```pip install -r meld.requirements.txt```

If you do not want to use the provided MELD scores for the EUR and EAS populations from the 1000 genomes database, you will need the following files from your own dataset:

* Plink binary format files (bed, bim, fam) partitioned by chromosome for efficient calculation
* A list.txt file with each SNP in the bim files listed on its own line

Once those files are generated, the following command can be run to generate MELD scores:

```python compute_ld_1000G.py --chrom 22 --win 100```

##MELD analysis

MELD requires the same dependencies as the original LD score framework. Instructions for installation and requirements can be found at https://github.com/bulik/ldsc. 

Once the appropriate packages are intalled in the environment the following command will initiate MELD analyses.

```python ldsc.py --h2 $path_to_sumstats$ --ref-ld-chr $path_to_Full_MELD_directory/$ --w-ld-chr $path_to_Full_LD_directory/$ --out $outfile_prefix --print-coefficients```

 ## Tutorials and Examples
 
We provide example code and a toy dataset which illustrate how to use MELD and conduct downstream analyses.
For instance, in order to analyze Triglyceride levels using the Biobank Japan Data, run:

```python ldsc.py --h2 BiobankJapan/sumstats/bbj.Triglyceride.sumstats.gz 
    --ref-ld-chr EAS_1000G/Full_MELD_/R2_0.07/MELD.win_100.@ 
    --w-ld-chr EAS_1000G/Full_LD_/R2_0.07/LD.win_100.@ 
    --out EAS.Triglyceride 
    --print-coefficients```

This will write out an EAS.Triglyceride.log file that contains MELD estimates using an alpha value of 0.07 and a MELD window of 100 SNPs.
 ## RELEVANT CITATIONS

G. Darnell*, S.P. Smith*, D. Udwin, S. Ramachandran, and L. Crawford. Partitioning marginal epistasis distinguishes confounding nonlinear effects from polygenicity in GWA summary statistics. _biorxiv_.

## QUESTIONS AND FEEDBACK
For questions or concerns with the MELD software, please contact [Sam Smith](mailto:samuel_smith1@brown.edu) or [Lorin Crawford](mailto:lcrawford@microsoft.com).

We welcome and appreciate any feedback you may have with our software and/or instructions.

##Directory contents

```BiobankJapan``` contains the sumstats and outputs directories whcih contain the summary statistic files and outputs for the Biobank Japan analyses, respectively. ```UKBiobank``` has the same file system for the UK Biobank analyses. 

```EAS_1000G``` contains the LD score and MELD score files calculated using the EAS superpopulation of the 1000 Genomes. The scores are located in the ```Full_LD_``` and ```Full_MELD_``` subdirectories, respectively. The same system is in ```EUR_1000G``` used for scores calculated in the EUR  superpopulation of the 1000 Genomes. 

```ldscore``` contains the necessary supplemnental scripts to perform MELD analysis and should always be located in the directory where MELD is being initiated. 
