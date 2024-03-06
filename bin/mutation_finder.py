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


reference_file = os.path.join(reference_file, subtype + '/' + segment + '.fasta')

print("python segment: {}".format(segment))
print(sequence_file)
print(reference_file)

def find_differences(reference, seq):
    # Align the sequences
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.open_gap_score = -10
    alignments = aligner.align(reference, seq)
    
    # Get the alignment with the highest score
    #best_alignment = max(alignments, key=lambda x: x.score)
    best_alignment = max(alignments)
    print(best_alignment)

    
    # Find the differences between the sequences
    differences = []
    #for i in range(len(best_alignment)-1):
    for i in range(len(best_alignment[1])):
        if best_alignment[0][i] != best_alignment[1][i]:
            differences.append(f"{best_alignment[0][i]}{i+1}{best_alignment[1][i]}")
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
