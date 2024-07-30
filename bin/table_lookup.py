import pandas as pd
import sys

# Arguments for file paths and criteria
mutations_file = sys.argv[1]
output_file = sys.argv[2]
xlsx_file = sys.argv[3]
segment = sys.argv[4]
subtype = sys.argv[5]
id = sys.argv[6]
type = sys.argv[7]  # Used for naming the mutations column dynamically

# Read mutations data from a CSV file
mutations_df = pd.read_csv(mutations_file)
print(mutations_df)
print(segment)

# Check if the value in the mutations column is NaN
if pd.notna(mutations_df.iloc[0, 1]):
    # Split and strip spaces from mutations
    sample_mutations = [mut.strip() for mut in mutations_df.iloc[0, 1].split(';')]
else:
    sample_mutations = []

list_mutations_set = set(sample_mutations)

print("Predefined list of mutations:", list_mutations_set)

# Read in the DataFrame from the Excel file
df = pd.read_excel(xlsx_file)

# Filter DataFrame based on segment and subtype
filtered_df = df[(df['segment'] == segment) & (df['subtype'] == subtype)]
#filtered_df = df[(df['segment'] == segment)]

print(f"Number of rows in filtered DataFrame: {len(filtered_df)}")

# Initialize a dictionary to hold results
results_dict = {}

# Loop through each row in the filtered DataFrame
for index, row in filtered_df.iterrows():
    mutations = row['mutation']
    # Check if mutations is not NaN and then split; otherwise, set to an empty set
    if pd.notna(mutations):
        row_mutations_set = set(mut.strip() for mut in str(mutations).split(';'))
    else:
        row_mutations_set = set()

    print(f"Processing row index {index} with mutations: {row_mutations_set}")

    # Find the intersection of row mutations and predefined mutations
    matching_mutations = row_mutations_set.intersection(list_mutations_set)
    print(f"Matching mutations in row {index}: {matching_mutations}")

    if matching_mutations:
        # Convert matching mutations to string
        common_mutations_str = ';'.join(matching_mutations)
        if id in results_dict:
            results_dict[id].add(common_mutations_str)  # Add new mutations to the set
        else:
            results_dict[id] = set([common_mutations_str])  # Start a new set of mutations for the sample

# Prepare results for DataFrame
results = [{'Sample': k, f"{segment} {type} mutations": ';'.join(v)} for k, v in results_dict.items()]

# Handle cases where no results are generated
if not results:
    df_output = pd.DataFrame([{
        f"{segment} {type} mutations": 'No matching mutations found',
        'Sample': mutations_df.iloc[0, 0]
    }])
else:
    df_output = pd.DataFrame(results)

df_output = df_output.drop_duplicates()
df_output.to_csv(output_file, index=False)

print("Final results captured:")
print(df_output)
