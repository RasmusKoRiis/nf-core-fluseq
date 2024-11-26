import pandas as pd
import sys

# Arguments for file paths and criteria
mutations_file = sys.argv[1]
output_file = sys.argv[2]
xlsx_file = sys.argv[3]
segment = sys.argv[4]
subtype = sys.argv[5]
sample_id = sys.argv[6]
mutation_type = sys.argv[7]  # Used for naming the mutations column dynamically

# If segment is NA change to NA1 because of Excel formatting of NA
segment_look = "NA1" if segment == "NA" else segment

# Read mutations data from a CSV file
mutations_df = pd.read_csv(mutations_file)

# Check if the value in the mutations column is NaN
if pd.notna(mutations_df.iloc[0, 1]):
    # Ensure that mutations are consistently stripped of spaces and converted to uppercase
    sample_mutations = [mut.strip().upper() for mut in mutations_df.iloc[0, 1].split(';') if pd.notna(mut)]
else:
    sample_mutations = []

# Convert sample mutations to a set for fast lookup
list_mutations_set = set(sample_mutations)

# Read the Excel file into a DataFrame
df = pd.read_excel(xlsx_file)

# Filter DataFrame based on segment
filtered_df = df[(df['segment'] == segment_look)]

# Initialize a dictionary to hold results
results_dict = {}

# Initialize a list to track samples with no matching mutations for specific conditions
no_matching_mutations_samples = []

# Loop through each row in the filtered DataFrame
for index, row in filtered_df.iterrows():
    mutations = row['mutation']
    # Ensure that mutations are consistently formatted (uppercase and stripped)
    if pd.notna(mutations):
        row_mutations_set = set(mut.strip().upper() for mut in str(mutations).split(';'))
    else:
        row_mutations_set = set()

    # Find the intersection of row mutations and predefined sample mutations
    matching_mutations = row_mutations_set.intersection(list_mutations_set)

    if matching_mutations:
        if sample_id in results_dict:
            # Directly update the set with matching mutations
            results_dict[sample_id].update(matching_mutations)
        else:
            # Start a new set of mutations for the sample
            results_dict[sample_id] = set(matching_mutations)
    else:
        # If no matching mutations and the segment is M2 and mutation type is inhibition
        if segment == "M2" and mutation_type.lower() == "inhibition":
            no_matching_mutations_samples.append(sample_id)

# Print only for cases where the segment is M2 and no inhibition mutations were found
if no_matching_mutations_samples:
    print(f"No inhibition mutations found for M2 segment in samples: {no_matching_mutations_samples}")

# Prepare results for the output DataFrame
results = [{'Sample': k, f"{segment} {mutation_type} mutations": ';'.join(v)} for k, v in results_dict.items()]

# Handle cases where no results are generated for the mutations
if not results:
    df_output = pd.DataFrame([{
        f"{segment} {mutation_type} mutations": 'No matching mutations found',
        'Sample': mutations_df.iloc[0, 0]
    }])
else:
    df_output = pd.DataFrame(results)

# Drop duplicates and write the output to a CSV file
df_output = df_output.drop_duplicates()
df_output.to_csv(output_file, index=False)

print(f"Results written to {output_file}")