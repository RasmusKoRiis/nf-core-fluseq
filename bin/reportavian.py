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


# Write the merged data to a new CSV file
merged_data.to_csv('merged_report.csv', index=False)