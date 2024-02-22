import pandas as pd
import sys

# Load the data
input_csv = sys.argv[1]
output_tsv = sys.argv[2]

#Read csv file and add column names
df=pd.read_table(input_csv,names=["query","hit","pident","length","mismatch","gapopen" ,"qstart" ,"qend" ,"sstart" ,"send" ,"evalue" ,"bitscore"])

# Convert the 'pident' column to numeric for sorting
df['pident'] = pd.to_numeric(df['pident'])

# Sort the DataFrame by 'pident' column in descending order
df = df.sort_values(by='pident', ascending=False)

# Drop duplicates, keeping only the first occurrence (highest value)
df = df.drop_duplicates(subset='query')

#sorts by query to get  correct order of segments PB1, PB2, PA ....
df = df.sort_values(by='query')

# Reset the index of the DataFrame
df = df.reset_index(drop=True)

# save csv file
df.to_csv(output_tsv,sep="\t", header= False, index = False)




