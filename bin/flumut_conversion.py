import pandas as pd
import re
import sys

# Check if the necessary arguments were provided
if len(sys.argv) < 3:
    print("Usage: python script.py <input_tsv_file> <id>")
    sys.exit(1)

# Load your data
input_file = sys.argv[1]
id = sys.argv[2]
df = pd.read_csv(input_file, sep='\t')

# Define columns for output
desired_columns = [
    'Sample',
    'Flumut_HA1', 'Effect_HA1',
    'Flumut_HA2', 'Effect_HA2',
    'Flumut_NA', 'Effect_NA',
    'Flumut_PB1', 'Effect_PB1',
    'Flumut_PB2', 'Effect_PB2',
    'Flumut_PA', 'Effect_PA',
    'Flumut_M1', 'Effect_M1',
    'Flumut_M2', 'Effect_M2',
    'Flumut_NS1', 'Effect_NS1',
    'Flumut_NP', 'Effect_NP'
]

# Initialize an empty DataFrame with desired columns
output_df = pd.DataFrame(columns=desired_columns)

# Function to get the correct column name based on the marker
def get_column_names(marker):
    marker = str(marker)
    if marker.startswith('HA1'):
        return 'Flumut_HA1', 'Effect_HA1'
    elif marker.startswith('HA2'):
        return 'Flumut_HA2', 'Effect_HA2'
    elif marker.startswith('NA'):
        return 'Flumut_NA', 'Effect_NA'
    elif marker.startswith('PB1'):
        return 'Flumut_PB1', 'Effect_PB1'
    elif marker.startswith('PB2'):
        return 'Flumut_PB2', 'Effect_PB2'
    elif marker.startswith('PA'):
        return 'Flumut_PA', 'Effect_PA'
    elif marker.startswith('M1'):
        return 'Flumut_M1', 'Effect_M1'
    elif marker.startswith('M2'):
        return 'Flumut_M2', 'Effect_M2'
    elif marker.startswith('NS-1'):
        return 'Flumut_NS1', 'Effect_NS1'
    elif marker.startswith('NP'):
        return 'Flumut_NP', 'Effect_NP'
    else:
        return None, None

# Iterate through each sample
grouped = df.groupby('Sample')
rows = []

# Identify mutation/effect columns
mutation_cols = [c for c in desired_columns if c.startswith("Flumut_")]
effect_cols = [c for c in desired_columns if c.startswith("Effect_")]

for sample, group in grouped:
    # Build lists then join with ';' at the end (dedup + stable order)
    mut_lists = {c: [] for c in mutation_cols}
    eff_lists = {c: [] for c in effect_cols}

    # Track seen mutations per segment/column
    seen = {c: set() for c in mutation_cols}

    for _, row in group.iterrows():
        mutation_column, effect_column = get_column_names(row.get('Marker', ''))
        if not mutation_column or not effect_column:
            continue

        # Remove segment and colon from the mutation
        mutation_raw = str(row.get('Mutations in your sample', '')).strip()
        mutation = re.sub(r'^[^:]+:', '', mutation_raw).strip()
        effect = str(row.get('Effect', '')).strip()

        if not mutation:
            continue

        # Deduplicate mutations within each mutation column
        if mutation not in seen[mutation_column]:
            seen[mutation_column].add(mutation)
            mut_lists[mutation_column].append(mutation)
            eff_lists[effect_column].append(effect)

    # Create final row
    row_data = {col: '' for col in desired_columns}
    row_data['Sample'] = sample

    # Join lists into ';' separated strings
    for c in mutation_cols:
        row_data[c] = ";".join(mut_lists[c])
    for c in effect_cols:
        row_data[c] = ";".join(eff_lists[c])

    rows.append(row_data)

# Create the output DataFrame from the list of rows
output_df = pd.DataFrame(rows, columns=desired_columns)

# Save the output DataFrame to a new CSV file
output_file = f"{id}_flumut_report.csv"
output_df.to_csv(output_file, index=False)

print(f"Data has been successfully consolidated into a single row per sample and saved to '{output_file}'")
