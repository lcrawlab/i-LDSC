# Marginal Epistatic LD score regression (MELD)

The inflation of test statistics in genome-wide association (GWA) studies due to confounding such as cryptic relatedness, population stratification, and spurious non-zero genetic effects driven by linkage disequilibrium (LD) has been well characterized in the literature. This repository holds **MELD**: a software package to run **marginal epistatic LD score regression**. MELD is an extended framework which takes in GWA test statistics and accurately partitions true additive genetic variation from confounding non-additive genetic variation, as well as other biases.

## Installation Requirements

We recommend creating an anaconda environment for running MELD, instructions for setting up conda environments can be found here: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

MELD requires a working version of pip and Python 3.9.12 to be installed.

To install the necessary python packages, navigate to the woking directory and submit:

```pip install -r requirements.txt```

## Calculation of MELD Scores (Optional)

If you do not want to use the provided MELD scores for the EUR and EAS populations from the 1000 genomes database, you will need the following files from your own dataset:

* Plink binary format files (bed, bim, fam) partitioned by chromosome for efficient calculation
* A list.txt file with each SNP in the bim files listed on its own line

Once those files are generated, the following command can be run to generate MELD scores:

```python compute_ld_1000G.py --chrom 22 --win 100```

 ## Tutorials and Examples
 
We provide example code and a toy dataset which illustrate how to use MELD and conduct downstream analyses.
 
 ## RELEVANT CITATIONS

G. Darnell*, S.P. Smith*, D. Udwin, S. Ramachandran, and L. Crawford. Partitioning marginal epistasis distinguishes confounding nonlinear effects from polygenicity in GWA summary statistics. _biorxiv_.

## QUESTIONS AND FEEDBACK
For questions or concerns with the MELD software, please contact [Sam Smith](mailto:samuel_smith1@brown.edu) or [Lorin Crawford](mailto:lcrawford@microsoft.com).

We welcome and appreciate any feedback you may have with our software and/or instructions.
