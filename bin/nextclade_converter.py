import pandas as pd
import sys

def transform_string(s):
    # Split the string by '-' and return the last part
    return s.split('-')[-1]

def ensure_columns(df: pd.DataFrame, cols, fill_value=pd.NA):
    """Add any missing columns with a constant fill value."""
    missing = [c for c in cols if c not in df.columns]
    if missing:
        print(f"[guard] Adding missing columns with NA: {', '.join(missing)}")
        # Use .loc to avoid chained-assignment warnings
        for c in missing:
            df.loc[:, c] = fill_value
    return df

def process_file(input_file, meta_id, segment, type_):
    # Read the CSV file with proper delimiter
    df = pd.read_csv(input_file, delimiter=';')

    # Remove the part of the string after '|' in the seqName column
    df['seqName'] = df['seqName'].apply(lambda x: x.split('|')[0] if isinstance(x, str) else x)

    # Rename columns in the DataFrame
    df.rename(columns={'seqName': 'Sample', 'aaSubstitutions': 'Differences'}, inplace=True)

    # -----------------------
    # CASE 1: Segment = HA
    # -----------------------
    if 'HA' in segment:
        columns_to_keep = [
            'Sample', 'legacy-clade', 'subclade', 'glycosylation', 'coverage',
            'frameShifts', 'Differences', 'aaDeletions', 'aaInsertions',
            'qc.overallStatus', 'qc.mixedSites.totalMixedSites'
        ]

        # Ensure all needed columns exist BEFORE slicing
        df = ensure_columns(df, columns_to_keep, fill_value=pd.NA)
        filtered_df = df.loc[:, columns_to_keep].copy()

        # Ensure no null values in empty columns
        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['glycosylation'] = filtered_df['glycosylation'].fillna('No glycosylation')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')
        filtered_df['qc.mixedSites.totalMixedSites'] = filtered_df['qc.mixedSites.totalMixedSites'].fillna('0')
        summary_df = filtered_df.copy()

        # Split the Differences into HA1 / HA2 / SigPep
        filtered_df['HA1'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('HA1:', '') for i in x.split(',') if isinstance(x, str) and 'HA1' in i])
        )
        filtered_df['HA2'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('HA2:', '') for i in x.split(',') if isinstance(x, str) and 'HA2' in i])
        )
        filtered_df['SigPep'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('SigPep:', '') for i in x.split(',') if isinstance(x, str) and 'SigPep' in i])
        )

        # -----------------------
        # Save HA1 mutations
        # -----------------------
        df_ha1 = filtered_df[
            [
                'Sample', 'legacy-clade', 'subclade', 'glycosylation',
                'coverage', 'frameShifts', 'HA1', 'aaDeletions',
                'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()

        df_ha1.rename(columns={'HA1': 'Differences'}, inplace=True)
        df_ha1.drop(columns=['legacy-clade', 'subclade', 'glycosylation', 'coverage', 'Differences'], inplace=True)

        # Rename columns to segment-specific names
        df_ha1.rename(columns={
            'frameShifts': f'frameShifts {segment}1',
            'aaDeletions': f'aaDeletions {segment}1',
            'aaInsertions': f'aaInsertions {segment}1',
            'qc.overallStatus': f'Nextclade QC {segment}1',
            'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}1'
        }, inplace=True)

        # Replace commas within cell values
        df_ha1 = df_ha1.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

        # Save
        df_ha1.to_csv(f'./{meta_id}_HA1_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered HA1 file saved as: ./{meta_id}_HA1_nextclade_{type_}_mutation.csv")

        # -----------------------
        # Save HA2 mutations
        # -----------------------
        df_ha2 = filtered_df[
            [
                'Sample', 'legacy-clade', 'subclade', 'glycosylation',
                'coverage', 'frameShifts', 'HA2', 'aaDeletions',
                'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()

        df_ha2.rename(columns={'HA2': 'Differences'}, inplace=True)
        df_ha2.drop(columns=['legacy-clade', 'subclade', 'glycosylation', 'coverage', 'Differences'], inplace=True)

        df_ha2.rename(columns={
            'frameShifts': f'frameShifts {segment}2',
            'aaDeletions': f'aaDeletions {segment}2',
            'aaInsertions': f'aaInsertions {segment}2',
            'qc.overallStatus': f'Nextclade QC {segment}2',
            'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}2'
        }, inplace=True)

        df_ha2 = df_ha2.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

        df_ha2.to_csv(f'./{meta_id}_HA2_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered HA2 file saved as: ./{meta_id}_HA2_nextclade_{type_}_mutation.csv")

        # -----------------------
        # HA summary file
        # -----------------------
        summary_df.drop(columns=[
            'Differences', 'coverage', 'frameShifts', 'aaDeletions', 'aaInsertions',
            'qc.overallStatus', 'qc.mixedSites.totalMixedSites'
        ], inplace=True, errors='ignore')

        def split_glycosylation(substitutions):
            if pd.isna(substitutions) or substitutions == "No glycosylation":
                return ["No glycosylation"] * 3
            changes = [c.strip() for c in str(substitutions).split(',') if c.strip()]
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

        summary_df[["HA_glycosylation_1", "HA_glycosylation_2", "HA_glycosylation_3"]] = pd.DataFrame(
            summary_df["glycosylation"].apply(split_glycosylation).tolist(),
            index=summary_df.index
        )

        summary_df = summary_df.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)
        summary_df.to_csv(f'./{meta_id}_nextclade_summary.csv', index=False)
        print(f"HA summary file saved as: ./{meta_id}_nextclade_summary.csv")

    # -----------------------
    # CASE 2: Segment = M
    # -----------------------
    elif 'M' in segment:
        columns_to_keep = [
            'Sample', 'clade', 'coverage', 'frameShifts', 'Differences',
            'aaDeletions', 'aaInsertions', 'qc.overallStatus',
            'qc.mixedSites.totalMixedSites'
        ]
        df = ensure_columns(df, columns_to_keep, fill_value=pd.NA)
        filtered_df = df.loc[:, columns_to_keep].copy()

        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')
        filtered_df['qc.mixedSites.totalMixedSites'] = filtered_df['qc.mixedSites.totalMixedSites'].fillna('0')

        filtered_df['M1'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('M1:', '') for i in str(x).split(',') if 'M1' in i])
        )
        filtered_df['M2'] = filtered_df['Differences'].apply(
            lambda x: ','.join([i.replace('M2:', '') for i in str(x).split(',') if 'M2' in i])
        )

        df_m1 = filtered_df[
            [
                'Sample', 'clade', 'coverage', 'frameShifts', 'M1',
                'aaDeletions', 'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()
        df_m1.rename(columns={'M1': 'Differences'}, inplace=True)
        df_m1.drop(columns=['clade', 'coverage', 'Differences'], inplace=True)
        df_m1.rename(columns={
            'frameShifts': f'frameShifts {segment}1',
            'aaDeletions': f'aaDeletions {segment}1',
            'aaInsertions': f'aaInsertions {segment}1',
            'qc.overallStatus': f'Nextclade QC {segment}1',
            'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}1'
        }, inplace=True)
        df_m1 = df_m1.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)
        df_m1.to_csv(f'./{meta_id}_M1_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered M1 file saved as: ./{meta_id}_M1_nextclade_{type_}_mutation.csv")

        df_m2 = filtered_df[
            [
                'Sample', 'clade', 'coverage', 'frameShifts', 'M2',
                'aaDeletions', 'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()
        df_m2.rename(columns={'M2': 'Differences'}, inplace=True)
        df_m2.drop(columns=['clade', 'coverage', 'Differences'], inplace=True)
        df_m2.rename(columns={
            'frameShifts': f'frameShifts {segment}2',
            'aaDeletions': f'aaDeletions {segment}2',
            'aaInsertions': f'aaInsertions {segment}2',
            'qc.overallStatus': f'Nextclade QC {segment}2',
            'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}2'
        }, inplace=True)
        df_m2 = df_m2.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)
        df_m2.to_csv(f'./{meta_id}_M2_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered M2 file saved as: ./{meta_id}_M2_nextclade_{type_}_mutation.csv")

    # -----------------------
    # CASE 3: Segment = NA
    # -----------------------
    elif 'NA' in segment:
        columns_to_keep = [
            'Sample', 'clade', 'coverage', 'frameShifts', 'Differences',
            'aaDeletions', 'aaInsertions', 'qc.overallStatus',
            'qc.mixedSites.totalMixedSites'
        ]
        df = ensure_columns(df, columns_to_keep, fill_value=pd.NA)
        filtered_df = df.loc[:, columns_to_keep].copy()

        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')

        df_na = filtered_df[
            [
                'Sample', 'clade', 'coverage', 'frameShifts', 'Differences',
                'aaDeletions', 'aaInsertions', 'qc.overallStatus',
                'qc.mixedSites.totalMixedSites'
            ]
        ].copy()

        df_na.drop(columns=['coverage', 'Differences'], inplace=True)
        df_na.rename(columns={
            'frameShifts': f'frameShifts {segment}',
            'aaDeletions': f'aaDeletions {segment}',
            'aaInsertions': f'aaInsertions {segment}',
            'clade': f'clade {segment}',
            'qc.overallStatus': f'Nextclade QC {segment}',
            'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}'
        }, inplace=True)

        df_na = df_na.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)
        df_na.to_csv(f'./{meta_id}_NA_nextclade_{type_}_mutation.csv', index=False)
        print(f"Filtered NA file saved as: ./{meta_id}_NA_nextclade_{type_}_mutation.csv")

    # -----------------------
    # CASE 4: All other segments
    # -----------------------
    else:
        columns_to_keep = [
            'Sample', 'clade', 'coverage', 'frameShifts', 'Differences',
            'aaDeletions', 'aaInsertions', 'qc.overallStatus',
            'qc.mixedSites.totalMixedSites'
        ]
        df = ensure_columns(df, columns_to_keep, fill_value=pd.NA)
        filtered_df = df.loc[:, columns_to_keep].copy()

        filtered_df['Differences'] = filtered_df['Differences'].fillna('No mutations')
        filtered_df['frameShifts'] = filtered_df['frameShifts'].fillna('No frameShifts')
        filtered_df['aaDeletions'] = filtered_df['aaDeletions'].fillna('No aaDeletions')
        filtered_df['aaInsertions'] = filtered_df['aaInsertions'].fillna('No aaInsertions')

        # Remove label prefixes like "PB2:" if there are real mutations
        if not all(filtered_df['Differences'].str.contains("No mutations", na=False)):
            filtered_df['Differences'] = filtered_df['Differences'].str.replace(r'\b\w+:', '', regex=True)

        filtered_df.drop(columns=['clade', 'coverage', 'Differences'], inplace=True)

        filtered_df.rename(columns={
            'frameShifts': f'frameShifts {segment}',
            'aaDeletions': f'aaDeletions {segment}',
            'aaInsertions': f'aaInsertions {segment}',
            'qc.overallStatus': f'Nextclade QC {segment}',
            'qc.mixedSites.totalMixedSites': f'Nextclade Mixed Sites {segment}'
        }, inplace=True)

        filtered_df = filtered_df.applymap(lambda x: x.replace(',', ';') if isinstance(x, str) else x)

        new_file_name = f"{meta_id}_{segment}_nextclade_{type_}_mutation.csv"
        new_file_path = f"./{new_file_name}"
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
