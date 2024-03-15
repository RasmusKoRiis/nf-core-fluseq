import pandas as pd
import sys

mutation_file = sys.argv[1]
mutation_list = sys.argv[2]
output_file = sys.argv[3]
segment = sys.argv[4]
id = sys.argv[5]
type = sys.argv[6]

print(segment)

# Read in the dataframes
df1 = pd.read_csv(mutation_file)
df2 = pd.read_csv(mutation_list)

# We convert them into sets for easier comparison
sample_mutations = set(';'.join(df1["Differences"].dropna()).split(';'))
print(sample_mutations)
list_mutations = set(';'.join(df2[type].dropna()).split(';'))
print(list_mutations)

# Find common mutations
common_mutations = sample_mutations.intersection(list_mutations)
print(common_mutations)

# Convert the set of common mutations back to a string separated by ';' to store in a single DataFrame cell
common_mutations_str = ';'.join(common_mutations)

mutations = f'{segment + type}  adaptation mutations'

# Create a DataFrame with 'id' as the index and the common mutations string as the value in the specified column
df_output = pd.DataFrame({mutations: [common_mutations_str]}, index=[id])

# Save the DataFrame to a new CSV
df_output.to_csv(output_file)


print(segment)
print(mutations)
print(common_mutations_str)




