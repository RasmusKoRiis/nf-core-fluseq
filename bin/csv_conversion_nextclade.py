import pandas as pd
import sys

# Get command-line arguments: first is csv file, second is an identifier for output files.
csv_file = sys.argv[1]
id = sys.argv[2]

# Load the CSV file into a DataFrame using ';' as the separator.
df = pd.read_csv(csv_file, sep=';')

# --- Step 0: Extract the segment information from the original seqName ---
# Example seqName: "000001-INFA|01-HA-H5N1"
# We assume that every row in the file uses the same format.
first_seq = df['seqName'].iloc[0]
# The segment is the letter(s) between the first and second '-' in the part after '|'
segment = first_seq.split('|')[1].split('-')[1]  # For the example, segment becomes "HA"

# --- Step 1: Create a working DataFrame with selected columns ---
# List of columns that we want to process.
target_cols = ["qc.overallStatus", "coverage", "aaDeletions", 
               "aaInsertions", "frameShifts", "aaSubstitutions", "totalNonACGTNs"]

# Extract the relevant columns; we keep the original 'seqName' for now.
nextclade_clean = df[['seqName'] + target_cols].copy()

# --- Step 2: Rename and update the seqName column ---
# Rename 'seqName' to 'Sample'
nextclade_clean.rename(columns={'seqName': 'Sample'}, inplace=True)
# Update the content of 'Sample' to only contain the part before the '|' character.
nextclade_clean['Sample'] = nextclade_clean['Sample'].apply(lambda x: x.split('|')[0] if isinstance(x, str) and '|' in x else x)

# --- Step 3: Define a function to remove prefix text up to and including colon (:) ---
def remove_prefix(val):
    """
    If the cell is a string, split on commas and for each part remove the text before (and including) the colon.
    If no colon is found in a part, it is returned unchanged.
    """
    if isinstance(val, str):
        parts = val.split(',')
        new_parts = [p.split(':', 1)[-1] if ':' in p else p for p in parts]
        return ",".join(new_parts)
    return val

# Apply the prefix removal function to each of the target columns.
for col in target_cols:
    nextclade_clean[col] = nextclade_clean[col].apply(remove_prefix)

# --- Step 4: Split the aaSubstitutions column into four columns according to custom rules ---
def split_aa_substitutions(substitutions):
    """
    Splits the aaSubstitutions string into up to 4 parts.
    If the value is NaN or "No mutation", returns a list with four entries of "No mutation".
    Otherwise, splits by comma and groups the changes in sets of 22.
    """
    if pd.isna(substitutions) or substitutions == "No mutation":
        return ["No mutation"] * 4
    changes = substitutions.split(',')
    if len(changes) <= 22:
        return [",".join(changes), "", "", ""]
    elif len(changes) <= 44:
        return [",".join(changes[:22]), ",".join(changes[22:]), "", ""]
    elif len(changes) <= 66:
        return [",".join(changes[:22]), ",".join(changes[22:44]), ",".join(changes[44:]), ""]
    else:
        return [",".join(changes[:22]), ",".join(changes[22:44]), ",".join(changes[44:66]), ",".join(changes[66:])]

# Create new columns by splitting the processed "aaSubstitutions" column.
# The splitting is done after the prefix removal.
split_cols = ['aaSubstitutions_1', 'aaSubstitutions_2', 'aaSubstitutions_3', 'aaSubstitutions_4']
nextclade_clean[split_cols] = pd.DataFrame(
    nextclade_clean['aaSubstitutions'].apply(split_aa_substitutions).tolist(),
    index=nextclade_clean.index
)

# --- Step 5: Rename the target columns using the segment variable ---
# Build a renaming dictionary for the original target columns.
renaming_dict = {col: f"NC_{segment}_{col}" for col in target_cols}
nextclade_clean.rename(columns=renaming_dict, inplace=True)

# Also, rename the split columns for substitutions using the same scheme.
split_renaming = {col: f"NC_{segment}_{col}" for col in split_cols}
nextclade_clean.rename(columns=split_renaming, inplace=True)

# --- Step 6: Replace empty mutation cells with defaults ---
# We apply custom defaults: for deletions use "No deletions", for insertions use "No insertions",
# and for all other mutation-related columns use "No mutations".
mutation_columns = [f"NC_{segment}_{col}" for col in target_cols + split_cols]
for col in mutation_columns:
    if "aaDeletions" in col:
        default_val = "No deletions"
    elif "aaInsertions" in col:
        default_val = "No insertions"
    else:
        default_val = "No mutations"
    nextclade_clean[col] = nextclade_clean[col].apply(
        lambda x: default_val if isinstance(x, str) and x.strip() == "" else x
    )

# --- (Optional) Step 7: Replace any remaining commas in cell values with semicolons ---
# If needed, replace all commas with semicolons to avoid delimiter conflicts.
nextclade_clean = nextclade_clean.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

# Make df for mutatins lookup

mutations_lookup = nextclade_clean[['Sample', f"NC_{segment}_aaSubstitutions"]].copy()
 

# --- Step 8: Save the results to CSV files ---
# Define output file names. Here, we save the same DataFrame into two files, as in your original script.
mutation_name = id + '_' + segment + '_nextclade_mutations.csv'
mutation_name_lookup = id + '_' + segment + '_nextclade_lookup_mutations.csv'

nextclade_clean.to_csv(mutation_name, index=False)
mutations_lookup.to_csv(mutation_name_lookup, index=False)
