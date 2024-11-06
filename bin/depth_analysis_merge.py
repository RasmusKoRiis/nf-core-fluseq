import pandas as pd
import glob
import sys

# MERGE ALL DATA
# Get a list of all CSV files in the current directory
csv_files = glob.glob('*.csv')

# Initialize an empty DataFrame to store the merged data
merged_data = pd.DataFrame()

# Check if there are any CSV files found
if not csv_files:
    print("No CSV files found in the current directory.")
    sys.exit(1)

# Iterate over each CSV file and concatenate them into one DataFrame
for file in csv_files:
    try:
        # Read the CSV file into a DataFrame
        df = pd.read_csv(file)

        # Concatenate with the merged_data DataFrame
        merged_data = pd.concat([merged_data, df], axis=0, ignore_index=True)
    except Exception as e:
        print(f"Error reading {file}: {e}")
        continue

# Write the merged data to a new CSV file if there is any data
if not merged_data.empty:
    merged_data.to_csv('depth_report.csv', index=False)
    print("Merged CSV file created: depth_report.csv")
else:
    print("No data found to merge.")
