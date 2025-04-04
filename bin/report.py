import pandas as pd
import glob
import sys

#MERGE ALL DATA
# Get a list of all CSV files in the current directory
csv_files = glob.glob('*.csv')
samplesheet = sys.argv[1]

#LOAD SAMPLE SHEET
samplesheet_df = pd.read_csv(samplesheet, sep='\t')

# Rename the column
samplesheet_df.rename(columns={'SequenceID': 'Sample'}, inplace=True)
# Remove Barcode columns
samplesheet_df = samplesheet_df.drop(columns=['Barcode'])

# Initialize an empty DataFrame to store the merged data
merged_data = None

# Iterate over each CSV file
for i, file in enumerate(csv_files):
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file)

    # If this is the first iteration, initialize merged_data with the first DataFrame
    if merged_data is None:
        merged_data = df
    else:
        # Merge the DataFrame with the merged_data DataFrame based on the "Sample" column
        merged_data = pd.concat([merged_data, df], axis=0, ignore_index=True)

# Group by 'Sample' and combine the rows
merged_data = merged_data.groupby('Sample', as_index=False).first()

#NIPH spesific adustment
#replace ! with - in the Sample column 
merged_data['Sample'] = merged_data['Sample'].str.replace('!', '-')

#ADD CALCULATED FILEDS

#FILL COLUMNS NOT ANALYZED WITH NA
keywords = {'mutation', 'differences', 'frameshift', 'aadeletions', 'aainsertions', 'subtype', 'nextclade qc', 'clade', 'subclade', 'clade na',
            'coverage-ha', 'coverage-na', 'coverage-mp', 'coverage-ns', 'coverage-pa', 'coverage-np', 'coverage-pb1', 'coverage-pb2', 'glycosylation', 
            'depth-ha', 'depth-na', 'depth-mp', 'depth-ns', 'depth-pa', 'depth-np', 'depth-pb1', 'depth-pb2'}

# Convert the column names to lowercase and check if any keyword is in each column name
columns_to_fill = [col for col in merged_data.columns if any(keyword in col.lower() for keyword in keywords)]

# Fill NA values in the identified columns
merged_data[columns_to_fill] = merged_data[columns_to_fill].fillna('NA')

# Function to check if a column exists, and if not, create it with 'NA' values
def ensure_column(df, column_name):
    if column_name not in df.columns:
        df[column_name] = 'NA'

# Ensure required columns exist
ensure_column(merged_data, 'M2 inhibtion mutations')
ensure_column(merged_data, 'NA inhibtion mutations')
ensure_column(merged_data, 'PA inhibtion mutations')

# RESISTANCE COLUMN
ensure_column(merged_data, 'DR_Res_Adamantine')
merged_data['DR_Res_Adamantine'] = merged_data['M2 inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANI' if 'No matching mutations' in x else 'Review'))

ensure_column(merged_data, 'DR_Res_Oseltamivir')
merged_data['DR_Res_Oseltamivir'] = merged_data['NA inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANI' if 'No matching mutations' in x else 'Review'))

ensure_column(merged_data, 'DR_Res_Zanamivir')
merged_data['DR_Res_Zanamivir'] = merged_data['NA inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANI' if 'No matching mutations' in x else 'Review'))

ensure_column(merged_data, 'DR_Res_Peramivir')
merged_data['DR_Res_Peramivir'] = merged_data['NA inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANI' if 'No matching mutations' in x else 'Review'))

ensure_column(merged_data, 'DR_Res_Laninamivir')
merged_data['DR_Res_Laninamivir'] = merged_data['NA inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANI' if 'No matching mutations' in x else 'Review'))

ensure_column(merged_data, 'DR_Res_Baloxavir')
merged_data['DR_Res_Baloxavir'] = merged_data['PA inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANS' if 'No matching mutations' in x else 'Review'))

# Create DR_M2_Mut column
ensure_column(merged_data, 'DR_M2_Mut')
merged_data['DR_M2_Mut'] = merged_data.apply(lambda x: 'NA' if x['M2 inhibtion mutations'] == 'NA' else (x['M2 inhibtion mutations'] if x['DR_Res_Adamantine'] == 'Review' else 'No Mutations'), axis=1)

# Create DR_PA_Mut column
ensure_column(merged_data, 'DR_PA_Mut')
merged_data['DR_PA_Mut'] = merged_data.apply(lambda x: 'NA' if x['PA inhibtion mutations'] == 'NA' else (x['PA inhibtion mutations'] if x['DR_Res_Baloxavir'] == 'Review' else 'No Mutations'), axis=1)

# Create DR_NA_Mut column
ensure_column(merged_data, 'DR_NA_Mut')
merged_data['DR_NA_Mut'] = merged_data.apply(
    lambda x: 'NA' if x['NA inhibtion mutations'] == 'NA' else (
        x['NA inhibtion mutations'] if x['DR_Res_Oseltamivir'] == 'Review' or x['DR_Res_Zanamivir'] == 'Review' or x['DR_Res_Peramivir'] == 'Review' or x['DR_Res_Laninamivir'] == 'Review' else 
        'No Mutations'
    ), 
    axis=1
)


# SUBTYPE COLUMN
merged_data['Sekvens_Resultat'] = merged_data['Subtype'].apply(lambda x: 'A/H3N2' if x == 'H3N2' else 'A/H1N1' if x == 'H1N1' else 'B/Victoria' if x == 'VICVIC' else 'B/Victoria' if x == 'VIC' else 'B/Yamagata' if x == 'YAMYAM' else 'B/Yamagata' if x == 'YAM' else x)


#REMOVE UNESSESARY COLUMNS
# Drop all columns that contain the word "mammalian"
merged_data = merged_data[merged_data.columns.drop(list(merged_data.filter(regex='mammalian')))]

required_columns = ['Sample', 'Sekvens_Resultat', 'Coverage-HA', 'Coverage-M', 'Coverage-NA', 'Coverage-NP', 'Coverage-NS', 'Coverage-PA', 
                    'Coverage-PB1', 'Coverage-PB2', 'DEPTH_HA', 'DEPTH_MP', 'DEPTH_NA', 'DEPTH_NP', 'DEPTH_NS', 'DEPTH_PA', 'DEPTH_PB1', 
                    'DEPTH_PB2', 'DR_M2_Mut', 'DR_NA_Mut', 'DR_PA_Mut', 'DR_Res_Adamantine', 'DR_Res_Baloxavir', 'DR_Res_Oseltamivir', 
                    'DR_Res_Peramivir', 'DR_Res_Zanamivir', 'HA1 Differences human', 'HA1 Differences human_vaccine', 'HA2 Differences human', 
                    'HA2 Differences human_vaccine', 'IRMA_altmatch', 'IRMA_chimeric', 'IRMA_failQC', 'IRMA_initial', 'IRMA_match', 
                    'IRMA_nomatch', 'IRMA_passQC', 'M1 Differences human', 'M2 Differences human', 'M2 Differences inhibition_human', 
                    'M2 inhibtion mutations', 'NA Differences human', 'NA Differences human_vaccine', 'NA Differences inhibition_human', 
                    'NA inhibtion mutations', 'NP Differences human', 'NS1 Differences human', 'Nextclade QC HA1', 'Nextclade QC M1', 
                    'Nextclade QC NA', 'Nextclade QC NP', 'Nextclade QC NS', 'Nextclade QC PA', 'Nextclade QC PB1', 'Nextclade QC PB2', 
                    'PA Differences human', 'PA Differences inhibition_human', 'PA inhibtion mutations', 'PB1 Differences human', 
                    'PB2 Differences human', 'SigPep Differences human', 'Subtype', 'aaDeletions HA1', 'aaDeletions M1', 'aaDeletions M2', 
                    'aaDeletions NA', 'aaDeletions NP', 'aaDeletions NS', 'aaDeletions PA', 'aaDeletions PB1', 'aaDeletions PB2', 
                    'aaInsertions HA1', 'aaInsertions M1', 'aaInsertions M2', 'aaInsertions NA', 'aaInsertions NP', 'aaInsertions NS', 
                    'aaInsertions PA', 'aaInsertions PB1', 'aaInsertions PB2', 'clade', 'clade NA', 'frameShifts HA1', 'frameShifts M1', 
                    'frameShifts M2', 'frameShifts NA', 'frameShifts NP', 'frameShifts NS', 'frameShifts PA', 'frameShifts PB1', 
                    'frameShifts PB2', 'glycosylation', 'subclade']

for column in required_columns:
    if column not in merged_data.columns:
        merged_data[column] = 'NA'

#SORT DF BY COLUMNS

merged_data = merged_data.sort_index(axis=1)
columns = ['Sample', 'Sekvens_Resultat'] + [col for col in merged_data.columns if col not in ['Sample', 'Sekvens_Resultat']]
merged_data = merged_data[columns]

# Find samples in df_new that are not in merged_data
new_samples = samplesheet_df[~samplesheet_df['Sample'].isin(merged_data['Sample'])]

# Append the new samples to merged_data without filling with NaN
merged_data = pd.concat([merged_data, new_samples], ignore_index=True)

# Ensure all columns from both dataframes are present and correctly ordered
all_columns = columns + [col for col in new_samples.columns if col not in columns]
merged_data = merged_data.reindex(columns=all_columns)

# Fill all empty values in merged_data with NaN
merged_data = merged_data.applymap(lambda x: 'NA' if pd.isna(x) else x)

# Ensure alle numeric columns have onlue 5 decimals

numeric_cols = merged_data.select_dtypes(include='number').columns
merged_data[numeric_cols] = merged_data[numeric_cols].round(5)
merged_data["IRMA_noise"] = pd.to_numeric(merged_data["IRMA_noise"], errors='coerce')
merged_data["IRMA_noise"] = merged_data["IRMA_noise"].round(5)

# Identify all columns that contain "Coverage" in their name
coverage_columns = [col for col in merged_data.columns if "Coverage" in col]

# Convert these columns to numeric and round to 2 decimal places
for col in coverage_columns:
    merged_data[col] = pd.to_numeric(merged_data[col], errors='coerce')
    merged_data[col] = merged_data[col].round(2)

merged_data[coverage_columns] = merged_data[coverage_columns].fillna('NA')


# Write the merged data to a new CSV file
merged_data.to_csv('merged_report.csv', index=False)

