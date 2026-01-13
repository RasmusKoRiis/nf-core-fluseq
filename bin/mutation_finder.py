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
    all_positions = []
    ref_pos = 1  # Biological indexing starts at 1

    for i in range(len(best_alignment[0])):
        ref_char = best_alignment[0][i]
        seq_char = best_alignment[1][i]

        if ref_char == "-":
            if seq_char != "-":
                if differences and differences[-1].startswith(f"ins{ref_pos}"):
                    differences[-1] += seq_char
                else:
                    differences.append(f"ins{ref_pos}{seq_char}")
            continue

        if seq_char == "-":  # Deletion
            differences.append(f"del{ref_pos}{ref_char}")
            all_positions.append(f"{ref_char}{ref_pos}-")
        else:
            token = f"{ref_char}{ref_pos}{seq_char}"
            all_positions.append(token)
            if ref_char != seq_char:
                differences.append(token)

        if ref_char != "-":
            ref_pos += 1

    return ";".join(differences), ";".join(all_positions)



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
        differences, all_positions = find_differences(reference, sequence)
        frameshift = check_frameshift(str(sequence))
        sequences.append({'ID': record.id, 'Differences': differences, 'All_Positions': all_positions})

# Create a DataFrame from the sequences
df = pd.DataFrame(sequences)

# Create a DataFrame from the sequences
df = pd.DataFrame(sequences)

if df.empty:
    print("Warning: No sequences were processed; DataFrame is empty.")
    # Create an empty DataFrame with the expected columns
    empty_df = pd.DataFrame(columns=["Sample", "Differences"])
    # Save the empty CSVs so downstream steps have output files
    empty_df.to_csv(output_file, index=False)
    output_file_report = output_file.replace('.csv', '_report.csv')
    empty_df.to_csv(output_file_report, index=False)
    full_output = output_file.replace('.csv', '_full_mutation_list.csv')
    full_output_report = full_output.replace('.csv', '_report.csv')
    empty_full_df = pd.DataFrame(columns=["Sample", "All_Positions"])
    empty_full_df.to_csv(full_output, index=False)
    empty_full_df.to_csv(full_output_report, index=False)
    sys.exit(0)

df['Differences'] = df.apply(process_differences, axis=1)


 
#df['Differences'] = df.apply(process_differences, axis=1)

# Split "ID" column into "sample" and "Ref_Name" columns
df[['Sample', 'Ref_Name']] = df['ID'].str.split('|', n=1, expand=True)

# Drop the original "ID" column
df.drop(columns=['Ref_Name'], inplace=True)
df.drop(columns=['ID'], inplace=True)


df = df[['Sample', 'Differences', 'All_Positions']]


# Remove any instance of 'ins...' or 'del...' in the 'Differences' column
df['Differences'] = df['Differences'].str.replace(r'\bins[^\s;]*;?|\bdel[^\s;]*;?', '', regex=True)

# Remove any trailing or leading semicolons that may remain
df['Differences'] = df['Differences'].str.strip(';')

# Replace empty strings with 'No mutations found'
df['Differences'] = df['Differences'].replace('', 'No mutations found')

# Save the final dataframe to a CSV file
df_main = df[['Sample', 'Differences']].copy()
df_main.to_csv(output_file, index=False)

output_file_report = output_file.replace('.csv', '_report.csv')
full_output = output_file.replace('.csv', '_full_mutation_list.csv')
full_output_report = full_output.replace('.csv', '_report.csv')

new_name = segment + ' ' + 'Differences' + ' ' + type
df_report = df_main.rename(columns={'Differences': new_name})

# Save the final dataframe to a CSV file
df_report.to_csv(output_file_report, index=False)

full_column_name = f"{segment} {type} full amino acid list"
df_full = df[['Sample', 'All_Positions']].copy()
df_full.rename(columns={'All_Positions': full_column_name}, inplace=True)
df_full.to_csv(full_output, index=False)
df_full.to_csv(full_output_report, index=False)
