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
for sample, group in grouped:
    # Create a dictionary to hold the consolidated row data
    row_data = {col: '' for col in desired_columns}
    row_data['Sample'] = sample
    
    # Iterate through the mutations in the group
    for _, row in group.iterrows():
        mutation_column, effect_column = get_column_names(row['Marker'])
        if mutation_column and effect_column:
            # Remove segment and colon from the mutation
            mutation = re.sub(r'^[^:]+:', '', row['Mutations in your sample'])
            
            # Append mutation and effect information to the corresponding columns
            if row_data[mutation_column]:
                row_data[mutation_column] += f";{mutation}"
                row_data[effect_column] += f";{row['Effect']}"
            else:
                row_data[mutation_column] = mutation
                row_data[effect_column] = row['Effect']
    
    # Append the row data to the list of rows
    rows.append(row_data)

# Create the output DataFrame from the list of rows
output_df = pd.DataFrame(rows, columns=desired_columns)

# Save the output DataFrame to a new CSV file
output_file = f"{id}_flumut_report.csv"
output_df.to_csv(output_file, index=False)

print(f"Data has been successfully consolidated into a single row per sample and saved to '{output_file}'")
