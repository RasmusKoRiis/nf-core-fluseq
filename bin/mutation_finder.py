import pandas as pd
import sys
import re
import os.path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner


sequence_file = sys.argv[1]
reference_file = sys.argv[2]
segment = sys.argv[3]
subtype = sys.argv[4]
output_file = sys.argv[5]
type = sys.argv[6]


reference_file = os.path.join(reference_file, type + '/' + subtype + '/' + segment + '.fasta')
print("python reference: {}".format(reference_file))

print("python segment: {}".format(segment))


# Function to align sequences and find differences
def find_differences(reference, seq):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10  # Penalty for opening a gap
    aligner.extend_gap_score = -1  # Penalty for extending a gap
    alignments = aligner.align(reference, seq)
    
    # Take the best alignment
    best_alignment = alignments[0]
    differences = []
    ref_pos = 1  # Biological indexing starts at 1

    for i in range(len(best_alignment[0])):
        ref_char = best_alignment[0][i]
        seq_char = best_alignment[1][i]

        if ref_char != seq_char:
            if ref_char == "-":  # Insertion
                if differences and differences[-1].startswith(f"ins{ref_pos}"):
                    differences[-1] += seq_char
                else:
                    differences.append(f"ins{ref_pos}{seq_char}")
            elif seq_char == "-":  # Deletion
                differences.append(f"del{ref_pos}{ref_char}")
            else:  # Substitution
                differences.append(f"{ref_char}{ref_pos}{seq_char}")

        if ref_char != "-":
            ref_pos += 1

    return ";".join(differences)



def check_frameshift(seq):
    if 'X' in seq:
        return seq.index('X') + 1
    else:
        return 'N.A.'

def process_differences(row):
        return row['Differences']
    

sequences = []

print(sequence_file)
for ref in SeqIO.parse(reference_file, 'fasta'):
    reference = Seq(str(ref.seq))
    for record in SeqIO.parse(sequence_file, 'fasta'):
        sequence = Seq(str(record.seq))
        differences = find_differences(reference, sequence)
        frameshift = check_frameshift(str(sequence))
        sequences.append({'ID': record.id, 'Differences': differences})

# Create a DataFrame from the sequences
df = pd.DataFrame(sequences)

 
df['Differences'] = df.apply(process_differences, axis=1)

# Split "ID" column into "sample" and "Ref_Name" columns
df[['Sample', 'Ref_Name']] = df['ID'].str.split('|', n=1, expand=True)

# Drop the original "ID" column
df.drop(columns=['Ref_Name'], inplace=True)
df.drop(columns=['ID'], inplace=True)


# Reorder the columns
df = df[['Sample', 'Differences']]


# Remove any instance of 'ins...' or 'del...' in the 'Differences' column
df['Differences'] = df['Differences'].str.replace(r'\bins[^\s;]*;?|\bdel[^\s;]*;?', '', regex=True)

# Remove any trailing or leading semicolons that may remain
df['Differences'] = df['Differences'].str.strip(';')

# Replace empty strings with 'No mutations found'
df['Differences'] = df['Differences'].replace('', 'No mutations found')

# Save the final dataframe to a CSV file
df.to_csv(output_file, index=False)

output_file_report = output_file.replace('.csv', '_report.csv')

new_name = segment + ' ' + 'Differences' + ' ' + type
df.rename(columns={'Differences': new_name}, inplace=True)

# Save the final dataframe to a CSV file
df.to_csv(output_file_report, index=False)
