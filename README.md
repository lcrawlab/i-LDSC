# Marginal Epistatic LD score regression (MELD)

The misestimation of test statistics in genome-wide association (GWA) studies due to confounding such as cryptic relatedness, population stratification, and spurious non-zero genetic effects driven by linkage disequilibrium (LD) has been well characterized in the literature. This repository contains code for **MELD**: a software package to run **marginal epistatic LD score regression**. MELD is an extended statistical framework which takes in GWA test statistics and accurately partitions true additive genetic variation from non-additive genetic variation, as well as other biases.

## Directory Contents

`BiobankJapan` contains subdirectories entitlted `sumstats` and `outputs` which, respectively, contain GWA summary statistics and results for the Biobank Japan analyses performed in the paper. 

`UKBiobank` has the same type of files as the `BiobankJapan` directory but for the UK Biobank analyses performed in the paper. 

`EAS_1000G` contains files of both the LD scores and the marginal epistatic LD scores calculated using the EAS superpopulation from the 1000 Genomes Project. The scores are located in the `Full_LD_` and `Full_MELD_` subdirectories, respectively.

`EUR_1000G` has the same type of score files as the `EAS_1000G` directory but computed from the EUR superpopulation of the 1000 Genomes Project.

`ldscore` contains the necessary supplemnental scripts to perform the MELD analysis and should always be located in the directory where MELD software is being initiated.

## Software Installation Requirements

We recommend creating an anaconda environment for running MELD. Detailed instructions for setting up conda environments can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html.

## Calculation of the Marginal Epistatic LD Scores (Optional)

Calculation of marginal epistatic LD scores requires a working version of pip and Python 3.9.12 to be installed.

To install the necessary python packages, navigate to the woking directory and use the command:

```pip install -r meld.requirements.txt```

If you do not want to use the provided MELD scores for the EUR and EAS populations from the 1000 Genomes database, you will need the following files from your own dataset:

* Plink binary format files (bed, bim, fam) partitioned by chromosome for efficient calculation;
* A list.txt file with each SNP in the bim files listed on its own line.

Once those files are generated, the following command can be used to generate MELD scores:

```python compute_ld_1000G.py --chrom 22 --win 100```

## Analyzing GWA Summary Statistics with MELD

MELD requires the same dependencies as the original LD score framework (LDSC) ([(Bulik-Sullivan et al. (2015)](https://www.nature.com/articles/ng.3211)). Instructions and software requirements for LDSC can be found at https://github.com/bulik/ldsc. 

Once the appropriate packages are intalled in the environment, the following command will initiate an analysis using the MELD model:

```python 
ldsc.py --h2 $path_to_sumstats$ 
--ref-ld-chr $path_to_Full_MELD_directory/$ 
--w-ld-chr $path_to_Full_LD_directory/$ 
--out $outfile_prefix 
--print-coefficients
```

 ## Tutorials and Examples
 
Here, we briefly provide example code which illustrates how to use MELD and conduct downstream analyses.
For instance, in order to analyze triglyceride levels using GWA summary statistics from Biobank Japan, one would use the command:

```python 
python ldsc.py --h2 BiobankJapan/sumstats/bbj.Triglyceride.sumstats.gz
--ref-ld-chr EAS_1000G/Full_MELD_/R2_0.07/MELD.win_100.@
--w-ld-chr EAS_1000G/Full_LD_/R2_0.07/LD.win_100.@
--out EAS.Triglyceride
--print-coefficients
 ```

This command will write out an EAS.Triglyceride.log file that contains MELD estimates using an alpha value of 0.07 and a generating window of 100 SNPs for the marginal epistatic LD scores.

 ## RELEVANT CITATIONS

G. Darnell*, S.P. Smith*, D. Udwin, S. Ramachandran, and L. Crawford. Partitioning marginal epistasis distinguishes nonlinear effects from polygenicity in GWA summary statistics. _biorxiv_.

## QUESTIONS AND FEEDBACK
For questions or concerns with the MELD software, please contact [Sam Smith](mailto:samuel_smith1@brown.edu) or [Lorin Crawford](mailto:lcrawford@microsoft.com).

We welcome and appreciate any feedback you may have with our software and/or instructions. 
