#!/usr/bin/env python3

import pandas as pd
import sys
import os

def process_pbs_file(input_file_path, column_title):
    try:
        df = pd.read_csv(input_file_path, sep='\t')

        # Sort requested column
        df.sort_values(by=column_title, ascending=False, inplace=True)

        # Rank values
        rank_column = f'Rank_{column_title}'
        df[rank_column] = df[column_title].rank(ascending=False, method='first')

        # Calculate Empirical P-Value
        pval_column = f'pval_{column_title}'
        df[pval_column] = df[rank_column] / df[rank_column].max()

        # Write file
        output_file_path = f'{os.path.splitext(input_file_path)[0]}_pval_{column_title}.txt'
        df.to_csv(output_file_path, sep='\t', index=False)

        # write new file with p-values
        thresholds = [0.001, 0.01, 0.05]
        threshold_values = {}
        for t in thresholds:
            # Finding the closest value to the p-value threshold
            closest_row = df.iloc[(df[pval_column] - t).abs().argsort()[:1]]
            closest_value = closest_row[column_title].values[0]
            threshold_values[f'p-value_{t}'] = closest_value

        thresholds_output_file_path = f'{os.path.splitext(input_file_path)[0]}_pval_{column_title}_thresholds.txt'
        with open(thresholds_output_file_path, 'w') as file:
            for key, value in threshold_values.items():
                file.write(f'{key}\t{value}\n')

        return output_file_path, thresholds_output_file_path

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

# Main
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python PBS_pval.py <input_file> <column_title>")
        sys.exit(1)

    input_file = sys.argv[1]
    column_title = sys.argv[2]

    output_file, thresholds_file = process_pbs_file(input_file, column_title)
    print(f"Processed file saved as {output_file}")
    print(f"Thresholds file saved as {thresholds_file}")

