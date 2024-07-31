import pandas as pd
import sys

# Arguments for file paths and criteria
irma_stat = sys.argv[1]
output_file = sys.argv[2]
id = sys.argv[3]

# Read the file into a DataFrame
df = pd.read_csv(irma_stat, sep='\t')

# Discard Patterns and PairsAndWidows columns
df = df[['Record', 'Reads']]

# Extract IRMA columns with '4-'
irma_cols = df[df['Record'].str.contains('4-')].copy()
irma_cols['Record'] = irma_cols['Record'].str.split('-').str[-1]
irma_cols.set_index('Record', inplace=True)

# Transpose the IRMA columns to get them as rows
irma_cols = irma_cols.T

# Rename columns to remove any prefix before the final underscore
irma_cols.columns = ['DEPTH_' + col.split('_')[1] for col in irma_cols.columns]

# Add sample column
meta_id = id
irma_cols.insert(0, 'Sample', meta_id)

# Create the final dataframe
final_df = pd.DataFrame({
    'Sample': [meta_id],
    'IRMA_initial': df.loc[df['Record'] == '1-initial', 'Reads'].values[0],
    'IRMA_failQC': df.loc[df['Record'] == '2-failQC', 'Reads'].values[0],
    'IRMA_passQC': df.loc[df['Record'] == '2-passQC', 'Reads'].values[0],
    'IRMA_chimeric': df.loc[df['Record'] == '3-chimeric', 'Reads'].values[0],
    'IRMA_nomatch': df.loc[df['Record'] == '3-nomatch', 'Reads'].values[0],
    'IRMA_match': df.loc[df['Record'] == '3-match', 'Reads'].values[0],
    'IRMA_altmatch': df.loc[df['Record'] == '3-altmatch', 'Reads'].values[0],
})

# Append DEPTH columns
for col in irma_cols.columns[1:]:
    final_df[col] = irma_cols[col].values[0]

#FILL UNANLYZED COLUMNS WITH NA
keywords = {'DEPTH_HA', 'DEPTH_NA', 'DEPTH_PB2', 'DEPTH_PB1', 'DEPTH_PA', 'DEPTH_NP', 'DEPTH_NS', 'DEPTH_MP' }
columns_to_fill = [col for col in merged_data.columns if any(keyword in col for keyword in keywords)]
final_df[columns_to_fill] = final_df[columns_to_fill].fillna('NA')

# Save to CSV
final_df.to_csv(output_file, index=False)
