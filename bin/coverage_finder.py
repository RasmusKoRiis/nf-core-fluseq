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
    if 'HA' in sequence:
        length = 1800
    elif 'NA' in sequence:
        length = 1450
    elif 'PB2' in sequence:
        length = 2400
    elif 'PB1' in sequence:
        length = 2400
    elif 'PA' in sequence:
        length = 2300
    elif 'NP' in sequence:
        length = 1600
    elif 'NS' in sequence:
        length = 920
    elif 'M' in sequence:
        length = 1100
    else:
        # Default length if no specific segment is found
        length = len(sequence)
    
    sequence_upper = sequence.upper()  # Convert sequence to uppercase
    n_count = sequence_upper.count('N')  # Count 'N' which now includes both 'N' and 'n'
    return (length - n_count ) / length * 100

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




