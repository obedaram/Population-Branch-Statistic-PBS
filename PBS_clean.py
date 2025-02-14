#!/usr/bin/env python3

import pandas as pd
import sys
import os

def process_pbs_file(input_file_path, output_file_path):
    df = pd.read_csv(input_file_path, sep='\t')

    # Remove rows with NaN
    df.dropna(inplace=True)

    # Extracting BP from the 'ID' column
    df['BP'] = df['ID'].apply(lambda x: x.split('_')[1])

    # Write new BP column
    cols = df.columns.tolist()
    cols = cols[:1] + ['BP'] + cols[1:-1]
    df = df[cols]

    # Output dataframe
    df.to_csv(output_file_path, sep='\t', index=False)

# Main function to execute the script
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python PBS_clean.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    base_name = os.path.splitext(input_file)[0]
    output_file = f"{base_name}_cleaned.txt"

    process_pbs_file(input_file, output_file)
    print(f"Cleaned results saved as {output_file}")
