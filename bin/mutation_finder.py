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

print(type)


reference_file = os.path.join(reference_file, type + '/' + subtype + '/' + segment + '.fasta')
print("python reference: {}".format(reference_file))

print("python segment: {}".format(segment))


def find_differences(reference, seq):
    # Align the sequences
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    alignments = aligner.align(reference, seq)
    
    best_alignment = max(alignments)
    print(best_alignment)
    
    differences = []
    ref_pos = 1  # Keep track of the position in the reference sequence
    for i in range(len(best_alignment[0])):
        ref_char = best_alignment[0][i]
        seq_char = best_alignment[1][i]

        if ref_char != seq_char:
            if ref_char == '-':  # Insertion in sequence compared to reference
                if differences and differences[-1].startswith(f"ins{ref_pos-1}"):
                    # If the last difference was an insertion at the same position, append this char to it
                    differences[-1] += seq_char
                else:
                    differences.append(f"ins{ref_pos-1}{seq_char}")  # Note the insertion point and the inserted char
            elif seq_char == '-':  # Deletion in sequence compared to reference
                differences.append(f"{ref_char}{ref_pos}del")  # Note the deletion
            else:  # Mismatch
                differences.append(f"{ref_char}{ref_pos}{seq_char}")

        if ref_char != '-':  # Don't advance reference position on insertions to reference
            ref_pos += 1

    return ';'.join(differences)
    


def check_frameshift(seq):
    if 'X' in seq:
        return seq.index('X') + 1
    else:
        return 'N.A.'

def process_differences(row):
        return row['Differences']
    
    #xs_count = sum([1 for x in row['Differences'].split(';') if 'X' in x])
    #if xs_count > 20:
    #    return 'Pos. FS ' + str(row['Frameshift/Poor Seq'])
    #else:
    #    return row['Differences']


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
df[['sample', 'Ref_Name']] = df['ID'].str.split('|', n=1, expand=True)

# Set the "Ref_Name" column to the 'segment' variable value
#df['Ref_Name'] = segment

# Combine the last character of the "sample" column with the first character of the "Ref_Name" column
#df['sample'] = df['sample'] + '_' + df['Ref_Name'].str[:1]

# Remove the first character of the "Ref_Name" column
#df['Ref_Name'] = df['Ref_Name'].str[2:]
#df['Ref_Name'] = segment

# Drop the original "ID" column
df.drop(columns=['Ref_Name'], inplace=True)
df.drop(columns=['ID'], inplace=True)


# Reorder the columns
df = df[['sample', 'Differences']]

# Save the final dataframe to a CSV file
df.to_csv(output_file, index=False)
