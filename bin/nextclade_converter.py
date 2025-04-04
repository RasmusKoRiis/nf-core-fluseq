import pandas as pd
import sys

def transform_string(s):
    # Split the string by '-' and return the last part
    return s.split('-')[-1]

def process_file(input_file, meta_id, segment, type_):
    # Read the CSV file with proper delimiter
    df = pd.read_csv(input_file, delimiter=';')
    # Remove the part of the string after '|' in the seqName column
    df['seqName'] = df['seqName'].apply(lambda x: x.split('|')[0] if isinstance(x, str) else x)
    # Rename columns in the DataFrame
    df.rename(columns={'seqName': 'Sample', 'aaSubstitutions': 'Differences'}, inplace=True)

    # Filter the columns
    if 'HA' in segment:
        columns_to_keep = [
            'Sample', 'clade', 'subclade', 'glycosylation', 'coverage',
            'frameShifts', 'Differences', 'aaDeletions', 'aaInsertions',
            'qc.overallStatus', 'qc.mixedSites.totalMixedSites'
        ]
    elif 'NA' in segment:
        columns_to_keep = [
            'Sample', 'clade', 'coverage', 'frameShifts', 'Differences',
            'aaDeletions', 'aaInsertions', 'qc.overallStatus',
            'qc.mixedSites.totalMixedSites'
        ]
    else:
        columns_to_keep = [
            'Sample', 'clade', 'coverage', 'frameShifts', 'Differences',
            'aaDeletions', 'aaInsertions', 'qc.overallStatus',
            'qc.mixedSites.totalMixedSites'
        ]

    filtered_df = df[columns_to_keep]

    # -----------------------
    # CASE 1: Segment = HA
    # -----------------------
    if 'HA' in segment:
        # Ensure no null values in empty columns
        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['glycosylation'] = filtered_df['glycosylation'].fillna('No glycosylation')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')
        filtered_df['qc.mixedSites.totalMixedSites'] = filtered_df['qc.mixedSites.totalMixedSites'].fillna('0')
        summary_df = filtered_df.copy()

        # Split the Differences into separate columns
        filtered_df['HA1'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('HA1:', '') for i in x.split(',') if 'HA1' in i])
        )
        filtered_df['HA2'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('HA2:', '') for i in x.split(',') if 'HA2' in i])
        )
        filtered_df['SigPep'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('SigPep:', '') for i in x.split(',') if 'SigPep' in i])
        )

        # -----------------------
        # Save HA1 mutations
        # -----------------------
        df_ha1 = filtered_df[
            [
                'Sample', 'clade', 'subclade', 'glycosylation',
                'coverage', 'frameShifts', 'HA1', 'aaDeletions',
                'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()

        df_ha1.rename(columns={'HA1': 'Differences'}, inplace=True)
        columns_to_remove = ['clade', 'subclade', 'glycosylation', 'coverage', 'Differences']
        df_ha1.drop(columns=columns_to_remove, inplace=True)

        # Rename columns to segment-specific names
        df_ha1.rename(columns={'frameShifts': f'frameShifts {segment}1'}, inplace=True)
        df_ha1.rename(columns={'aaDeletions': f'aaDeletions {segment}1'}, inplace=True)
        df_ha1.rename(columns={'aaInsertions': f'aaInsertions {segment}1'}, inplace=True)
        df_ha1.rename(columns={'qc.overallStatus': f'Nextclade QC {segment}1'}, inplace=True)
        df_ha1.rename(columns={'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}1'}, inplace=True)

        # Replace all commas within cell values
        df_ha1 = df_ha1.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

        # Save
        df_ha1.to_csv(f'./{meta_id}_HA1_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered HA1 file saved as: ./{meta_id}_HA1_nextclade_{type_}_mutation.csv")

        # -----------------------
        # Save HA2 mutations
        # -----------------------
        df_ha2 = filtered_df[
            [
                'Sample', 'clade', 'subclade', 'glycosylation',
                'coverage', 'frameShifts', 'HA2', 'aaDeletions',
                'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()

        df_ha2.rename(columns={'HA2': 'Differences'}, inplace=True)
        columns_to_remove = ['clade', 'subclade', 'glycosylation', 'coverage', 'Differences']
        df_ha2.drop(columns=columns_to_remove, inplace=True)

        # Rename columns to segment-specific names
        df_ha2.rename(columns={'frameShifts': f'frameShifts {segment}2'}, inplace=True)
        df_ha2.rename(columns={'aaDeletions': f'aaDeletions {segment}2'}, inplace=True)
        df_ha2.rename(columns={'aaInsertions': f'aaInsertions {segment}2'}, inplace=True)
        df_ha2.rename(columns={'qc.overallStatus': f'Nextclade QC {segment}2'}, inplace=True)
        df_ha2.rename(columns={'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}2'}, inplace=True)

        # Replace commas within cell values
        df_ha2 = df_ha2.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

        # Save
        df_ha2.to_csv(f'./{meta_id}_HA2_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered HA2 file saved as: ./{meta_id}_HA2_nextclade_{type_}_mutation.csv")

        # -----------------------
        # Make HA summary file
        # -----------------------
        columns_to_remove = [
            'Differences', 'coverage', 'frameShifts', 'aaDeletions', 'aaInsertions',
            'Differences', 'qc.overallStatus', 'qc.mixedSites.totalMixedSites'
        ]
        summary_df.drop(columns=columns_to_remove, inplace=True)

        # Function for splitting glycosylation column
        def split_glycosylation(substitutions):
            if pd.isna(substitutions) or substitutions == "No glycosylation":
                return ["No glycosylation"] * 3
            changes = substitutions.split(',')
            changes = [c.strip() for c in changes]  # remove extra spaces
            if len(changes) <= 22:
                return [",".join(changes), "", ""]
            elif len(changes) <= 44:
                return [",".join(changes[:22]), ",".join(changes[22:]), ""]
            else:
                return [
                    ",".join(changes[:22]),
                    ",".join(changes[22:44]),
                    ",".join(changes[44:])
                ]

        # Split the glycosylation into separate columns
        summary_df[["HA_glycosylation_1", "HA_glycosylation_2", "HA_glycosylation_3"]] = pd.DataFrame(
            summary_df["glycosylation"].apply(split_glycosylation).tolist(),
            index=summary_df.index
        )

        # Replace commas within cell values
        summary_df = summary_df.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)
        summary_df.to_csv(f'./{meta_id}_nextclade_summary.csv', index=False)
        print(f"HA summary file saved as: ./{meta_id}_nextclade_summary.csv")

    # -----------------------
    # CASE 2: Segment = M
    # -----------------------
    elif 'M' in segment:
        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')
        filtered_df['qc.mixedSites.totalMixedSites'] = filtered_df['qc.mixedSites.totalMixedSites'].fillna('0')

        # Split the Differences
        filtered_df['M1'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('M1:', '') for i in x.split(',') if 'M1' in i])
        )
        filtered_df['M2'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('M2:', '') for i in x.split(',') if 'M2' in i])
        )

        # -----------------------
        # Save M1
        # -----------------------
        df_m1 = filtered_df[
            [
                'Sample', 'clade', 'coverage', 'frameShifts', 'M1',
                'aaDeletions', 'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()

        df_m1.rename(columns={'M1': 'Differences'}, inplace=True)
        columns_to_remove = ['clade', 'coverage', 'Differences']
        df_m1.drop(columns=columns_to_remove, inplace=True)

        # Rename columns
        df_m1.rename(columns={'frameShifts': f'frameShifts {segment}1'}, inplace=True)
        df_m1.rename(columns={'aaDeletions': f'aaDeletions {segment}1'}, inplace=True)
        df_m1.rename(columns={'aaInsertions': f'aaInsertions {segment}1'}, inplace=True)
        df_m1.rename(columns={'qc.overallStatus': f'Nextclade QC {segment}1'}, inplace=True)
        df_m1.rename(columns={'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}1'}, inplace=True)

        # Replace commas within cell values
        df_m1 = df_m1.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

        # Save
        df_m1.to_csv(f'./{meta_id}_M1_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered M1 file saved as: ./{meta_id}_M1_nextclade_{type_}_mutation.csv")

        # -----------------------
        # Save M2
        # -----------------------
        df_m2 = filtered_df[
            [
                'Sample', 'clade', 'coverage', 'frameShifts', 'M2',
                'aaDeletions', 'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()

        df_m2.rename(columns={'M2': 'Differences'}, inplace=True)
        columns_to_remove = ['clade', 'coverage', 'Differences']
        df_m2.drop(columns=columns_to_remove, inplace=True)

        # Rename columns
        df_m2.rename(columns={'frameShifts': f'frameShifts {segment}2'}, inplace=True)
        df_m2.rename(columns={'aaDeletions': f'aaDeletions {segment}2'}, inplace=True)
        df_m2.rename(columns={'aaInsertions': f'aaInsertions {segment}2'}, inplace=True)
        df_m2.rename(columns={'qc.overallStatus': f'Nextclade QC {segment}2'}, inplace=True)
        df_m2.rename(columns={'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}2'}, inplace=True)

        # Replace commas within cell values
        df_m2 = df_m2.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

        # Save
        df_m2.to_csv(f'./{meta_id}_M2_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered M2 file saved as: ./{meta_id}_M2_nextclade_{type_}_mutation.csv")

    # -----------------------
    # CASE 3: Segment = NA
    # -----------------------
    elif 'NA' in segment:
        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')

        # Save NA
        df_na = filtered_df[
            [
                'Sample', 'clade', 'coverage', 'frameShifts', 'Differences',
                'aaDeletions', 'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()

        columns_to_remove = ['coverage', 'Differences']
        df_na.drop(columns=columns_to_remove, inplace=True)

        # Rename columns
        df_na.rename(columns={'frameShifts': f'frameShifts {segment}'}, inplace=True)
        df_na.rename(columns={'aaDeletions': f'aaDeletions {segment}'}, inplace=True)
        df_na.rename(columns={'aaInsertions': f'aaInsertions {segment}'}, inplace=True)
        df_na.rename(columns={'clade': f'clade {segment}'}, inplace=True)
        df_na.rename(columns={'qc.overallStatus': f'Nextclade QC {segment}'}, inplace=True)
        df_na.rename(columns={'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}'}, inplace=True)

        # Replace commas within cell values
        df_na = df_na.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

        # Save
        df_na.to_csv(f'./{meta_id}_NA_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered NA file saved as: ./{meta_id}_NA_nextclade_{type_}_mutation.csv")

    # -----------------------
    # CASE 4: All other segments
    # -----------------------
    else:
        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')

        # If you want to remove everything before ":" in Differences:
        # but only if the row actually has any real mutations
        # (so we don't break the "No mutations" string)
        if not all(filtered_df['Differences'].str.contains("No mutations", na=False)):
            filtered_df['Differences'] = filtered_df['Differences'].str.replace(r'\b\w+:', '', regex=True)

        columns_to_remove = ['clade', 'coverage', 'Differences']
        filtered_df.drop(columns=columns_to_remove, inplace=True)

        # Rename columns
        filtered_df.rename(columns={'frameShifts': f'frameShifts {segment}'}, inplace=True)
        filtered_df.rename(columns={'aaDeletions': f'aaDeletions {segment}'}, inplace=True)
        filtered_df.rename(columns={'aaInsertions': f'aaInsertions {segment}'}, inplace=True)
        filtered_df.rename(columns={'qc.overallStatus': f'Nextclade QC {segment}'}, inplace=True)
        filtered_df.rename(columns={'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}'}, inplace=True)

        # Replace commas within cell values
        filtered_df = filtered_df.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

        # Create the new filename
        new_file_name = f"{meta_id}_{segment}_nextclade_{type_}_mutation.csv"
        new_file_path = f"./{new_file_name}"

        # Save
        filtered_df.to_csv(new_file_path, index=False)
        print(f"Filtered file saved as: {new_file_path}")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_file> <meta_id> <segment> <type_>")
        sys.exit(1)

    input_file = sys.argv[1]
    meta_id = sys.argv[2]
    segment = sys.argv[3]
    type_ = sys.argv[4]

    segment = transform_string(segment)
    process_file(input_file, meta_id, segment, type_)
