# Population-Branch-Statistic-PBS
Script for Calculating the Population Branch Statistic



This script calculates the Population Branch Statistic (PBS) using the equation supplied in the supplemental materials in Yi et al., 2010.
This script uses the Hudson's Fst Estimator as recommended by Bhattia et al. 2013. It also outputs a normalization of the branch lengths using Crawford et al. 2017's modification.

The script uses the allele frequency files from plink2, but can be modified to use other inputs. 
It assumes that your files are filtered to keep only biallelic sites, and that all your files have the same REF alleles.

The file can be run using the following commandline: 
python3 PBS.py <input1>.afreq <input2>.afreq <input3>.afreq <output>.txt
