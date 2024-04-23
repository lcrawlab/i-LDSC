# Interaction-LD Score (i-LDSC) Regression

LD score regression (LDSC) was developed to distinguish true genetic signal from confounding biases such as cryptic relatedness and population stratification. The key concept underlying the LDSC framework is that there is a positive linear relationship between test statistics derived from genome-wide association studies (GWAS) and linkage disequilibrium (LD) when complex traits are generated under the infinitesimal model --- that is, when all genetic variants equally contribute to phenotypic variation and their effects are uniformly distributed along the genome. This repository contains code for **interaction-LD score (i-LDSC) regression**: an extended framework which aims to recover missing heritability from GWAS summary statistics by incorporating an additional score that measures the amount of non-additive genetic variation that is tagged by each variant in the data. After model fitting, these scores produce nonzero regression coefficients when the distribution of genetic effects is non-uniform due to a subset of variants being involved in interactions that also contribute to trait architecture.

## Directory Contents
In addition to the files below, files of LD scores, cis-interaction LD scores, and GWAS summary statistics used for our analysis of the UK Biobank and BioBank Japan can be downloaded from the [Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/W6MA8J&faces-redirect=true). Please be aware that these files are quiet large and will take time to download completely. 

`1000G_EAS` contains cis-interaction LD scores calculated on the East Asian (EAS) superpopulation of the 1000 Genomes Project phase 3. These scores are calculated across a range of alpha values, which account for the trait specific MAF-effect size relationship. 

`EUR_1000G` has the same type of score files as the `EAS_1000G` directory but computed from the European (EUR) superpopulation of the 1000 Genomes Project.

`ldscore` contains the necessary supplemental scripts to perform i-LDSC analysis and should always be located in the directory where i-LDSC software is being initiated.

`compute_ld_cisinteractionscores.py` is the python script used to calculate cis-interaction LD scores from PLINK formatted genetic data. See description of use below.

`ildsc.py` is the command line tool for launching i-LDSC analysis. See description of use below.

`ildsc.score.requirements.txt` contains the python packages necessary for an environment to perform calculation of cis-interaction LD scores with the `compute_ld_1000G.py` script.

`munge_sumstats.py` is the python script originally developed by [Bulik-Sullivan et al. (2015)](https://www.nature.com/articles/ng.3211) to format GWAS summary statistics for analysis. 

`simulations` contains the python scripts necessary to generate the simulated data included in the manuscript, including: various alpha and window size parameters, GxE architecture under an amplification model, and GxAncestry architecture using principal component (PC) loadings.

`figures` contains both the main text and supplemental figures that are presented in the manuscript as well as the scripts and data necessary to reproduce them.

## Recommendation: Partially Clone the Repo

Due to the size of the data and summary statistics, we recommend that users consider [partially cloning](https://docs.gitlab.com/ee/topics/git/partial_clone.html) the repo to reduce memory for their working copy and downloading only the subdirectories that are relevant for their purposes. To do this, one would first clone the repo without including any files:

```git clone --filter=blob:none --sparse https://github.com/lcrawlab/i-LDSC.git Test```

Then after changing the directory `cd i-LDSC`, one could copy over all of the python scripts needed to run analyses using the command line:

```git sparse-checkout set --no-cone \*.py \*requirements.txt \*.md```

If one wanted to extract the cis-interaction score for the EAS superpopulation from the 1000 Genomes phase 3 data, one would use (as an example):

```git sparse-checkout set 1000G_EAS```

If one wanted to extract the cis-interaction and additive LD scores from the EUR and EAS reference panels, one would use (as an example):

```git sparse-checkout set --no-cone \*.ldscore.gz```

For more information on using the partial clone function, please see more details [here](https://docs.gitlab.com/ee/topics/git/partial_clone.html).

## Software Installation Requirements

We recommend creating an anaconda environment for running i-LDSC. Detailed instructions for setting up conda environments can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html.

## Calculation of the cis-interaction LD scores (Optional)

Calculation of cis-interaction LD scores requires a working version of pip and Python 3.9.12 to be installed.

To install the necessary python packages, navigate to the woking directory and use the command:

```pip install -r ildsc.score.requirements.txt```

If you do not want to use the provided cis-interaction LD scores for the EUR and EAS populations from the 1000 Genomes database, you will need the following files from your own dataset:

* Plink binary format files (bed, bim, fam) partitioned by chromosome for efficient calculation;
* A list.txt file with each SNP in the bim files listed on its own line.

Once those files are generated, commands similar to the following can be used to generate cis-interaction LD scores:

```python compute_ld_cisinteractionscores.py --plink_dir $path_to_plink_files/plink_file_prefix$ EAS_1000G/ --chrom 22 --win 100```

They are available upon request. 

If one would like to use the i-LDSC framework using their own GWAS summary statistics, we would highly recommend that they convert their data to the `.sumstats` file format that the LDSC software recognizes (see instructions [here](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics)).

## Analyzing GWAS Summary Statistics with i-LDSC

i-LDSC requires the same dependencies as the original LD score regression framework (LDSC) ([Bulik-Sullivan et al. 2015](https://www.nature.com/articles/ng.3211)). Instructions and software requirements for LDSC can be found at https://github.com/bulik/ldsc. 

_NOTE: i-LDSC analyses are run in python 2, while the computation of the cis-interaction LD scores is done in python 3. Therefore, two separate environments are needed to move between tasks._

Once the appropriate packages are intalled in the environment, the following command will initiate an analysis using the i-LDSC model:

```python 
ildsc.py --h2 $path_to_sumstats$ 
--ref-ld-chr $path_to_cisinteractionscore_directory/$ 
--w-ld-chr $path_to_weight_directory/$ 
--out $outfile_prefix 
--overlap-annot
--print-coefficients
```

 ## Tutorials and Examples
 
Here, we briefly provide example code which illustrates how to use i-LDSC and conduct downstream analyses.
For instance, in order to analyze triglyceride levels using GWAS summary statistics from BioBank Japan, one would use the command:

```python 
python ildsc.py --h2 sumstats/Height_norm.hapmap.1kg.sumstats.gz
--ref-ld-chr cisinteraction_scores/
--w-ld-chr weight_hm3_no_hla/
--out UKB.Triglyceride
--print-coefficients
 ```

This command will write out an UKB.Triglyceride.log file that contains i-LDSC estimates using an alpha value of 0.07 and a generating window of 100 SNPs for the cis-interaction LD scores.

 ## RELEVANT CITATIONS

S.P. Smith*, G. Darnell*, D. Udwin, J. Stamp, A. Harpak, S. Ramachandran, and L. Crawford. Discovering non-additive heritability using additive GWAS summary statistics. _eLife_. In Press.

## QUESTIONS AND FEEDBACK
For questions or concerns with the i-LDSC software, please contact [Sam Smith](mailto:samuel.smith@utexas.edu) or [Lorin Crawford](mailto:lcrawford@microsoft.com).

We welcome and appreciate any feedback you may have with our software and/or instructions. 
