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
    sequence_upper = sequence.upper()  # Convert sequence to uppercase
    n_count = sequence_upper.count('N')  # Count 'N' which now includes both 'N' and 'n'
    return (1 - n_count / len(sequence)) * 100

# Read FASTA file
fasta_file = fasta

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


# Save DataFrame to CSV file
txt_file = df_pivot.loc[name, segment]
filename = name + '_' + segment + '_coverage.txt'

with open(filename, 'w') as file:
    file.write(str(txt_file))

df_pivot.to_csv(output)




