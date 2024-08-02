import pandas as pd
import glob

#MERGE ALL DATA

# Get a list of all CSV files in the current directory
csv_files = glob.glob('*.csv')

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
keywords = {'mutation', 'Differences', 'frameShift', 'aaDeletions', 'aaInsertions', 'Subtype', 'Nextclade QC','clade','subclade', 'clade NA',
            'Coverage-HA', 'Coverage-NA', 'Coverage-MP', 'Coverage-NS', 'Coverage-PA', 'Coverage-NP', 'Coverage-PB1', 'Coverage-PB2', 'glycosylation', 
            'DEPTH-HA', 'DEPTH-NA', 'DEPTH-MP', 'DEPTH-NS', 'DEPTH-PA', 'DEPTH-NP', 'DEPTH-PB1', 'DEPTH-PB2' }

columns_to_fill = [col for col in merged_data.columns if any(keyword in col for keyword in keywords)]
merged_data[columns_to_fill] = merged_data[columns_to_fill].fillna('NA')


# RESISTANCE COLUMN
merged_data['DR_Res_Adamantine'] = merged_data['M2 inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANI' if 'No matching mutations' in x else 'Review'))
merged_data['DR_Res_Oseltamivir'] = merged_data['NA inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANS' if 'No matching mutations' in x else 'Review'))
merged_data['DR_Res_Zanamivir'] = merged_data['NA inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANS' if 'No matching mutations' in x else 'Review'))
merged_data['DR_Res_Peramivir'] = merged_data['NA inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANS' if 'No matching mutations' in x else 'Review'))
merged_data['DR_Res_Baloxavir'] = merged_data['PA inhibtion mutations'].apply(lambda x: 'NA' if x == 'NA' else ('AANS' if 'No matching mutations' in x else 'Review'))

# Create DR_M2_Mut column
merged_data['DR_M2_Mut'] = merged_data.apply(lambda x: 'NA' if x['M2 inhibtion mutations'] == 'NA' else (x['M2 inhibtion mutations'] if x['DR_Res_Adamantine'] == 'Review' else 'L26;V27;A30;S31;G34;L38'), axis=1)

# Create DR_PA_Mut column
merged_data['DR_PA_Mut'] = merged_data.apply(lambda x: 'NA' if x['PA inhibtion mutations'] == 'NA' else (x['PA inhibtion mutations'] if x['DR_Res_Baloxavir'] == 'Review' else 'E23;L28;K34;A36;A37;I38;E119;E198;E199'), axis=1)

# Create DR_NA_Mut column
merged_data['DR_NA_Mut'] = merged_data.apply(lambda x: 'NA' if x['NA inhibtion mutations'] == 'NA' else (x['NA inhibtion mutations'] if (x['DR_Res_Adamantine'] == 'Review' or x['DR_Res_Oseltamivir'] == 'Review' or x['DR_Res_Zanamivir'] == 'Review' or x['DR_Res_Peramivir'] == 'Review') else 'E119;Q136;T148;D151;I222;R224;N245;N245-;A246-;T247-;G248-;K249-;A250-;K249;E27'), axis=1)

# SUBTYPE COLUMN
merged_data['Sekvens_Resultat'] = merged_data['Subtype'].apply(lambda x: 'A/H3N2' if x == 'H3N2' else 'A/H1N1' if x == 'H1N1' else 'B/Victoria' if x == 'VICVIC' else x)


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
                    'Nextclade QC NA1', 'Nextclade QC NP1', 'Nextclade QC NS1', 'Nextclade QC PA1', 'Nextclade QC PB11', 'Nextclade QC PB21', 
                    'PA Differences human', 'PA Differences inhibition_human', 'PA inhibtion mutations', 'PB1 Differences human', 
                    'PB2 Differences human', 'SigPep Differences human', 'Subtype', 'aaDeletions HA1', 'aaDeletions M1', 'aaDeletions M2', 
                    'aaDeletions NA', 'aaDeletions NP', 'aaDeletions NS', 'aaDeletions PA', 'aaDeletions PB1', 'aaDeletions PB2', 
                    'aaInsertions HA1', 'aaInsertions M1', 'aaInsertions M2', 'aaInsertions NA', 'aaInsertions NP', 'aaInsertions NS', 
                    'aaInsertions PA', 'aaInsertions PB1', 'aaInsertions PB2', 'clade', 'clade NA', 'frameShifts HA1', 'frameShifts M1', 
                    'frameShifts M2', 'frameShifts NA', 'frameShifts NP', 'frameShifts NS', 'frameShifts PA', 'frameShifts PB1', 
                    'frameShifts PB2', 'glycosylation', 'subclade', 'RunID', 'Instrument ID']

for column in required_columns:
    if column not in merged_data.columns:
        merged_data[column] = 'NA'

#SORT DF BY COLUMNS
merged_data = merged_data.sort_index(axis=1)
columns = ['Sample', 'Sekvens_Resultat'] + [col for col in merged_data.columns if col not in ['Sample', 'Sekvens_Resultat']]
merged_data = merged_data[columns]

# Write the merged data to a new CSV file
merged_data.to_csv('merged_report.csv', index=False)