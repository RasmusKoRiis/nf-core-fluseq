import pandas as pd
import sys
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


def find_differences(reference, seq):
    # Align the sequences
    from Bio.Align import PairwiseAligner
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    alignments = aligner.align(reference, seq)
    
    best_alignment = max(alignments)
    print(best_alignment)
    
    differences = []
    ref_pos = 1  # Start from 1 to match biological sequence indexing conventions
    
    for i in range(len(best_alignment[0])):
        ref_char = best_alignment[0][i]
        seq_char = best_alignment[1][i]

        if ref_char != seq_char:
            if ref_char == '-':  # Insertion in sequence compared to reference
                if differences and differences[-1].startswith(f"ins{ref_pos}"):
                    # If the last difference was an insertion at the same position, append this char to it
                    differences[-1] += seq_char
                else:
                    # Append the insertion at the current position and adjust ref_pos backwards by one
                    differences.append(f"ins{ref_pos}{seq_char}")
                    ref_pos += 1  # Move back ref_pos as we are inserting, not moving along the reference
            elif seq_char == '-':  # Deletion in sequence compared to reference
                # Note the deletion at the current position without adjustment
                differences.append(f"del{ref_pos}{ref_char}")
                ref_pos -= 1  # Move back ref_pos as we are inserting, not moving along the reference
            else:  # Mismatch
                # Note the mismatch at the current position without adjustment
                differences.append(f"{ref_char}{ref_pos}{seq_char}")

        # Increment ref_pos if the current character in reference is not an insertion
        if ref_char != '-':
            ref_pos += 1

    return ';'.join(differences)

    


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
df['Differences'] = df['Differences'].replace('', 'No mutations found')

# Save the final dataframe to a CSV file
df.to_csv(output_file, index=False)

output_file_report = output_file.replace('.csv', '_report.csv')

new_name = segment + ' ' + 'Differences' + ' ' + type
df.rename(columns={'Differences': new_name}, inplace=True)

# Save the final dataframe to a CSV file
df.to_csv(output_file_report, index=False)
