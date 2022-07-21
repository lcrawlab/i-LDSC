# Marginal Epistatic LD score regression (MELD)

The misestimation of test statistics in genome-wide association (GWA) studies due to confounding such as cryptic relatedness, population stratification, and spurious non-zero genetic effects driven by linkage disequilibrium (LD) has been well characterized in the literature. This repository contains code for **MELD**: a software package to run **marginal epistatic LD score regression**. MELD is an extended statistical framework which takes in GWA test statistics and accurately partitions true additive genetic variation from non-additive genetic variation, as well as other biases.

## Directory Contents

`BioBankJapan` contains subdirectories entitlted `sumstats` and `outputs` which, respectively, contain GWA summary statistics and results for the BioBank Japan analyses performed in the paper. 

`UKBiobank` has the same type of files as the `BioBankJapan` directory but for the UK Biobank analyses performed in the paper. 

`EAS_1000G` contains files of both the LD scores and the marginal epistatic LD scores calculated using the EAS superpopulation from the 1000 Genomes Project. The scores are located in the `Full_LD_` and `Full_MELD_` subdirectories, respectively.

`EUR_1000G` has the same type of score files as the `EAS_1000G` directory but computed from the EUR superpopulation of the 1000 Genomes Project.

`ldscore` contains the necessary supplemental scripts to perform the MELD analysis and should always be located in the directory where MELD software is being initiated.

`compute_ld_1000G.py` is the python script used to calculate MELD scores from plink formatted genetic data. See description of use below.

`meld.py` is the command line tool for launching MELD analysis. See description of use below.

`meld.score.requirements.txt` contains the python packages necessary for an environment to perform calculation of MELD scores with the `compute_ld_1000G.py` script.

`munge_sumstats.py` is the python script originally developed by [Bulik-Sullivan et al. (2015)](https://www.nature.com/articles/ng.3211) to format GWA summary statistics for analysis. 

## Recommendation: Partially Clone the Repo

Due to the size of the data and summary statistics, we recommend that users consider [partially cloning](https://docs.gitlab.com/ee/topics/git/partial_clone.html) the repo to reduce memory for their working copy and downloading only the subdirectories that are relevant for their purposes. To do this, one would first clone the repo without including any files:

```git clone --filter=blob:none --sparse https://github.com/lcrawlab/MELD.git Test```

Then after changing the directory `cd MELD`, one could copy over all of the python scripts needed to run analyses using the command line:

```git sparse-checkout set --no-cone \*.py \*requirements.txt \*.md```

If one wanted to extract the UK Biobank or BioBank Japan results and corresponding summary statistics, one would use (as an example):

```git sparse-checkout set "UK Biobank"```

If one wanted to extract the MELD and LDSC scores from the EUR and EAS reference panels, one would use (as an example):

```git sparse-checkout set --no-cone \*.ldscore.gz```

For more information on using the partial clone function, please see more details [here](https://docs.gitlab.com/ee/topics/git/partial_clone.html).

## Software Installation Requirements

We recommend creating an anaconda environment for running MELD. Detailed instructions for setting up conda environments can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html.

## Calculation of the Marginal Epistatic LD Scores (Optional)

Calculation of marginal epistatic LD scores requires a working version of pip and Python 3.9.12 to be installed.

To install the necessary python packages, navigate to the woking directory and use the command:

```pip install -r meld.score.requirements.txt```

If you do not want to use the provided MELD scores for the EUR and EAS populations from the 1000 Genomes database, you will need the following files from your own dataset:

* Plink binary format files (bed, bim, fam) partitioned by chromosome for efficient calculation;
* A list.txt file with each SNP in the bim files listed on its own line.

Once those files are generated, commands similar to the following can be used to generate MELD scores:

```python compute_ld_1000G.py --chrom 22 --win 100```

If one would like to use the MELD framework using their own GWA summary statistics, we would highly recommend that they convert their data to the `.sumstats` file format that the LDSC software recognizes (see instructions [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics)).

## Analyzing GWA Summary Statistics with MELD

MELD requires the same dependencies as the original LD score framework (LDSC) ([Bulik-Sullivan et al. 2015](https://www.nature.com/articles/ng.3211)). Instructions and software requirements for LDSC can be found at https://github.com/bulik/ldsc. 

NOTE: MELD analyses are run in python 2 and the computation of MELD scores are run in python 3. You must make a new environment to obtain MELD estimates. 

Once the appropriate packages are intalled in the environment, the following command will initiate an analysis using the MELD model:

```python 
meld.py --h2 $path_to_sumstats$ 
--ref-ld-chr $path_to_Full_MELD_directory/$ 
--w-ld-chr $path_to_Full_LD_directory/$ 
--out $outfile_prefix 
--print-coefficients
```

 ## Tutorials and Examples
 
Here, we briefly provide example code which illustrates how to use MELD and conduct downstream analyses.
For instance, in order to analyze triglyceride levels using GWA summary statistics from BioBank Japan, one would use the command:

```python 
python ldsc.py --h2 BioBankJapan/sumstats/bbj.Triglyceride.sumstats.gz
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
