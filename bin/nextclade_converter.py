import pandas as pd
import sys

def process_file(input_file, meta_id, segment):
    # Read the CSV file with proper delimiter
    df = pd.read_csv(input_file, delimiter=';')
    # Remove the part of the string after '|' in the seqName column
    df['seqName'] = df['seqName'].apply(lambda x: x.split('|')[0] if isinstance(x, str) else x)
    # Rename columns in the DataFrame
    df.rename(columns={'seqName': 'Sample', 'aaSubstitutions': 'Differences'}, inplace=True)


    # Filter the columns
    if 'HA' in segment:
        columns_to_keep = ['Sample', 'clade', 'subclade', 'glycosylation', 'coverage', 'frameShifts', 'Differences', 'aaDeletions', 'aaInsertions']
    elif 'NA' in segment:
        columns_to_keep = ['Sample', 'clade', 'glycosylation', 'coverage', 'frameShifts', 'Differences', 'aaDeletions', 'aaInsertions']
    else:
        columns_to_keep = ['Sample', 'clade', 'coverage', 'frameShifts', 'Differences', 'aaDeletions', 'aaInsertions']
                           

    filtered_df = df[columns_to_keep]
   
    if 'HA' in segment:

        # Ensure no null values in empty columns
        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['glycosylation'] = filtered_df['glycosylation'].fillna('No glycosylation')
        filtered_df['frameShifts'] = filtered_df['glycosylation'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')
        summary_df = filtered_df.copy()



        # Split the Differences into separate columns
        filtered_df['HA1'] = filtered_df['Differences'].apply(lambda x: ','.join([i.replace('HA1:', '') for i in x.split(',') if 'HA1' in i]))
        filtered_df['HA2'] = filtered_df['Differences'].apply(lambda x: ','.join([i.replace('HA2:', '') for i in x.split(',') if 'HA2' in i]))            
        filtered_df['SigPep'] = filtered_df['Differences'].apply(lambda x: ','.join([i.replace('SigPep:', '') for i in x.split(',') if 'SigPep' in i]))

        # Filter and save the HA1 mutations
        df_ha1 = filtered_df[['Sample', 'clade', 'subclade', 'glycosylation', 'coverage', 'frameShifts', 'HA1', 'aaDeletions', 'aaInsertions']].copy()
        df_ha1.rename(columns={'HA1': 'Differences'}, inplace=True)
        columns_to_remove = ['clade', 'subclade', 'glycosylation', 'coverage']
        df_ha1.drop(columns=columns_to_remove, inplace=True)
        print(type)
        df_ha1.to_csv(f'./{meta_id}_HA1_nextclade_{type}_mutation.csv', index=False)
        print(f"Filtered HA1 file saved as: ./{meta_id}_HA1_nextclade_mutation.csv")
            
        # Filter and save the HA2 mutations
        df_ha2 = filtered_df[['Sample', 'clade', 'subclade', 'glycosylation', 'coverage', 'frameShifts', 'HA2', 'aaDeletions', 'aaInsertions']].copy()
        df_ha2.rename(columns={'HA2': 'Differences'}, inplace=True)
        columns_to_remove = ['clade', 'subclade', 'glycosylation', 'coverage']
        df_ha2.drop(columns=columns_to_remove, inplace=True)
        df_ha2.to_csv(f'./{meta_id}_HA2_nextclade_{type}_mutation.csv', index=False)
        print(f"Filtered HA2 file saved as: ./{meta_id}_HA2_nextclade_mutation.csv")
            
        # Filter and save the SigPep mutations
        df_sigpep = filtered_df[['Sample', 'clade', 'subclade', 'glycosylation', 'coverage', 'frameShifts', 'SigPep', 'aaDeletions', 'aaInsertions']].copy()
        df_sigpep.rename(columns={'SigPep': 'Differences'}, inplace=True)
        columns_to_remove = ['clade', 'subclade', 'glycosylation', 'coverage']
        df_sigpep.drop(columns=columns_to_remove, inplace=True)
        df_sigpep.to_csv(f'./{meta_id}_SigPep_nextclade_{type}_mutation.csv', index=False)
        print(f"Filtered SigPep file saved as: ./{meta_id}_SigPep_nextclade_mutation.csv")


        # Make Sample Nextclade summary file
        columns_to_remove = ['frameShifts', 'Differences', 'aaDeletions', 'aaInsertions']
        summary_df.drop(columns=columns_to_remove, inplace=True)
        summary_df.to_csv(f'./{meta_id}_nextclade_summary.csv', index=False)      

    elif 'M' in segment:
         # Ensure no null values in empty columns
        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')

        # Split the Differences into separate columns
        filtered_df['M1'] = filtered_df['Differences'].apply(lambda x: ','.join([i.replace('M1:', '') for i in x.split(',') if 'M1' in i]))
        filtered_df['M2'] = filtered_df['Differences'].apply(lambda x: ','.join([i.replace('M2:', '') for i in x.split(',') if 'M2' in i]))            

        # Filter and save the M2 mutations
        df_m1 = filtered_df[['Sample', 'clade', 'coverage', 'frameShifts', 'M1', 'aaDeletions', 'aaInsertions']].copy()
        df_m1.rename(columns={'M1': 'Differences'}, inplace=True)
        columns_to_remove = ['clade', 'coverage']
        df_m1.drop(columns=columns_to_remove, inplace=True)
        df_m1.to_csv(f'./{meta_id}_M1_nextclade_mutation.csv', index=False)
        print(f"Filtered M1 file saved as: ./{meta_id}_M1_nextclade_{type}_mutation.csv")
            
        # Filter and save the M2 mutations
        df_m2 = filtered_df[['Sample', 'clade', 'coverage', 'frameShifts', 'M2', 'aaDeletions', 'aaInsertions']].copy()
        df_m2.rename(columns={'M2': 'Differences'}, inplace=True)
        columns_to_remove = ['clade', 'coverage']
        df_m2.drop(columns=columns_to_remove, inplace=True)
        df_m2.to_csv(f'./{meta_id}_M2_nextclade_{type}_mutation.csv', index=False)
        print(f"Filtered M2 file saved as: ./{meta_id}_M2_nextclade_mutation.csv")
    else:

        # Ensure no null values in ewmpty columns
        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')

        # Remove everything before ':' in the Differences column
        filtered_df['Differences'] = df['Differences'].str.replace(r'\b\w+:', '', regex=True)

        columns_to_remove = ['clade', 'coverage']
        filtered_df.drop(columns=columns_to_remove, inplace=True)


        # Create the new filename
        new_file_name = f"{meta_id}_{segment}_nextclade_{type}_mutation.csv"

        # Save the filtered dataframe to a new CSV file
        new_file_path = f"./{new_file_name}"
        filtered_df.to_csv(new_file_path, index=False)

        print(f"Filtered file saved as: {new_file_path}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_file> <meta_id> <segment>")
        sys.exit(1)

    input_file = sys.argv[1]
    meta_id = sys.argv[2]
    segment = sys.argv[3]
    type = sys.argv[4]

    process_file(input_file, meta_id, segment)
