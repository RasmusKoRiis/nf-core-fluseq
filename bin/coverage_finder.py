from Bio import SeqIO
import pandas as pd
import os
import sys

fasta = sys.argv[1]
output = sys.argv[2]
name = sys.argv[3]
segment = sys.argv[4]

# Function to calculate coverage
def calculate_coverage(sequence):
    n_count = sequence.count('N')
    return (1 - n_count / len(sequence)) * 100

# Function to extract sample name
def extract_sample_name(filepath):
    filename = os.path.basename(filepath)
    parts = filename.split('_')
    print(parts)
    sample_name = parts[0] + '_N'
    return sample_name

# Read FASTA file
fasta_file = fasta
#sample_name = extract_sample_name(fasta_file)

# Initialize a list to store data
data = []

# Parse each record in the FASTA file
for record in SeqIO.parse(fasta_file, "fasta"):
    coverage = calculate_coverage(str(record.seq))
    data.append({'id': record.id, 'coverage': coverage, 'Sample': name, 'Segment': segment})

# Create DataFrame
df = pd.DataFrame(data)



# Pivot table to have one row per sample and separate columns for each 'id'
df_pivot = df.pivot(index='Sample', columns='Segment', values='coverage')

# Add 'coverage status' column
#df_pivot['coverage status'] = df_pivot.apply(lambda row: 'fail' if any(row < 90) else 'ok', axis=1)

# Select only coverage columns for calculating the average
#coverage_columns = df_pivot.columns[df_pivot.columns != 'coverage status']
#df_pivot['average coverage'] = df_pivot[coverage_columns].mean(axis=1)


# Rename columns that contain 'HA', 'NA', 'PA'
#for col in df_pivot.columns:
  #  if 'HA' in col:
  #      df_pivot.rename(columns={col: 'HA'}, inplace=True)
 #   elif 'NA' in col:
 #       df_pivot.rename(columns={col: 'NA'}, inplace=True)
  #  elif 'PA' in col:
    #    df_pivot.rename(columns={col: 'PA'}, inplace=True)

# Add 'HA', 'NA', 'PA' columns if they are missing
#for segment in ['HA', 'NA', 'PA']:
   # if segment not in df_pivot.columns:
     #   df_pivot[segment] = None

#df_pivot = df_pivot.drop(['HA', 'NA', 'PA'], axis=1)

# Save DataFrame to CSV file
txt_file = df_pivot.loc[name, segment]
filename = name + '_' + segment + '_coverage.txt'

with open(filename, 'w') as file:
    file.write(str(txt_file))

df_pivot.to_csv(output)




