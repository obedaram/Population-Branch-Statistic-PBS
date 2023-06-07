#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import math
import os

# Load the CSV as dataframes, and specify the output files 
csv_file1 = sys.argv[1]
csv_file2 = sys.argv[2]
csv_file3 = sys.argv[3]
output_file = sys.argv[4]
output_info_file = os.path.splitext(output_file)[0] + "_pop_info.txt"

# Read CSV files into dataframes
df1 = pd.read_csv(csv_file1, sep='\t', low_memory=False)
df2 = pd.read_csv(csv_file2, sep='\t', low_memory=False)
df3 = pd.read_csv(csv_file3, sep='\t', low_memory=False)


# Calculations are performed for each row in the dataframe using Hudson's Fst Estimator as given in Bhattia et al. 2013. 
results = []
for index, row in df1.iterrows():
    maf1 = row['ALT_FREQS']
    maf2 = df2.loc[index, 'ALT_FREQS']
    maf3 = df3.loc[index, 'ALT_FREQS']
    nchrobs1 = df1.loc[index, 'OBS_CT']
    nchrobs2 = df2.loc[index, 'OBS_CT']
    nchrobs3 = df3.loc[index, 'OBS_CT']

#Calculates Fst between populations 1 and 2
    hw12 = (((maf1 - maf2) ** 2) - ((maf1 * (1 - maf1)) / (nchrobs1 - 1)) - ((maf2 * (1 - maf2)) / (nchrobs2 - 1)))
    hb12 = (maf1 * (1 - maf2) + maf2 * (1 - maf1))
    fst12 = (hw12 / hb12) if not (math.isnan(hw12) or math.isnan(hb12) or math.isinf(hw12) or math.isinf(hb12) or abs(hb12) < 1e-10) else np.nan

#Calculates Fst between populations 2 and 3
    hw23 = (((maf2 - maf3) ** 2) - ((maf2 * (1 - maf2)) / (nchrobs2 - 1)) - ((maf3 * (1 - maf3)) / (nchrobs3 - 1)))
    hb23 = (maf2 * (1 - maf3) + maf3 * (1 - maf2))
    fst23 = (hw23 / hb23) if not (math.isnan(hw23) or math.isnan(hb23) or math.isinf(hw23) or math.isinf(hb23) or abs(hb23) < 1e-10) else np.nan

 #Calculates Fst between populations 1 and 3 
    hw13 = (((maf1 - maf3) ** 2) - ((maf1 * (1 - maf1)) / (nchrobs1 - 1)) - ((maf3 * (1 - maf3)) / (nchrobs3 - 1)))
    hb13 = (maf1 * (1 - maf3) + maf3 * (1 - maf1))
    fst13 = (hw13 / hb13) if not (math.isnan(hw13) or math.isnan(hb13) or math.isinf(hw13) or math.isinf(hb13) or abs(hb13) < 1e-10) else np.nan

    if math.isnan(fst12) or math.isnan(fst23) or math.isnan(fst13):
        # Skip calculations and set 'nan' for dependent variables
        T12, T23, T13, bl1, bl2, bl3, nbl1, nbl2, nbl3 = [np.nan] * 9
    else:
        if fst12 >= 1 or fst23 >= 1 or fst13 >= 1:
            # Set 'nan' if the value inside the log function is invalid
            T12, T23, T13, bl1, bl2, bl3, nbl1, nbl2, nbl3 = [np.nan] * 9
        else:
            
            # Calculates PBS based on the equation in the supplemental materials in Yi et al., 2010
            
            # Calculate population divergence time T for population pairs
            T12 = -math.log10(1 - fst12)
            T13 = -math.log10(1 - fst13)
            T23 = -math.log10(1 - fst23)

            # Calculate the PBS branch lengths
            bl1 = (T12 + T13 - T23) / 2
            bl2 = (T12 + T23 - T13) / 2
            bl3 = (T23 + T13 - T12) / 2
            
            # The normalization of the branch lengths follows Crawford et al. 2017's method to scale the statistic and avoid artificially high PBS values when differentiation was low or high between all groups. 

            nbl1 = bl1 / (1 + bl1 + bl2 + bl3)
            nbl2 = bl2 / (1 + bl1 + bl2 + bl3)
            nbl3 = bl3 / (1 + bl1 + bl2 + bl3)
    
    # Prints the total NCHROBS for each SNP in each population. You can use this as a check to see how many individuals were genotyped for that marker and in what order they were loaded into the dataframe
    NCHROBStotal = str(nchrobs1) + ',' + str(nchrobs2) + ',' + str(nchrobs3)

    chr_val = df1.loc[index, '#CHROM']
    snp_val = df1.loc[index, 'ID']
    a1_val = df1.loc[index, 'REF']
    a2_val = df1.loc[index, 'ALT']

    results.append((chr_val, snp_val, a1_val, a2_val, fst12, fst23, fst13, bl1, bl2, bl3, nbl1, nbl2, nbl3, NCHROBStotal))

# Column headers for the new dataframe
column_headers = ['CHR', 'ID', 'REF', 'ALT', 'fst12', 'fst23', 'fst13', 'PBS_1', 'PBS_2', 'PBS_3', 'norm_PBS_1', 'norm_PBS_2', 'norm_PBS_3', 'Individuals_Pop']

# Updating the header for the results
output_df = pd.DataFrame(results, columns=column_headers)

# add NaN to all the values where it couldn't calculate the FST
output_df = output_df.fillna("NaN")

# Print the final dataframe to the output tab-delimited text file
output_df.to_csv(output_file, sep='\t', index=False)

# Print the pop info file for the scan
with open(output_info_file, 'w') as info_file:
    info_file.write(f"For {output_file}\n")
    info_file.write(f"pop1 = {csv_file1}\n")
    info_file.write(f"pop2 = {csv_file2}\n")
    info_file.write(f"pop3 = {csv_file3}\n")
