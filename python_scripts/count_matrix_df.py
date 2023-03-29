import os
import pandas as pd
import numpy as np

np.seterr(divide = "ignore")

# create a new directory to save the extracted data files
if not os.path.exists('extracted_data'):
    os.makedirs('extracted_data')

# read in the count matrix file
count_matrix = pd.read_csv('orig_matrix.csv', index_col=0)

# extract the gene information from the first column
gene = count_matrix.index.values

# loop over each i-th column
for i in range(1, count_matrix.shape[1]):
    # extract the i-th column and column heading
    column = count_matrix.iloc[:, i]
    column_heading = count_matrix.columns[i]
    summation=sum(column[1:count_matrix.shape[1]])

    # extract the source of the brain tissue from the last row of the i-th column
    tissue_source = count_matrix.iloc[-1, i]

    # create a new dataframe with the gene, i-th column, and tissue source
    new_df = pd.DataFrame({'gene': gene, column_heading: column, 'tissue_source': tissue_source})
    new_df.drop(new_df.tail(36).index, inplace = True)
    new_df[column_heading] = pd.to_numeric(new_df[column_heading])
    new_df[column_heading] = np.where(new_df[column_heading] > 0,  np.log10((new_df[column_heading] / summation) * 1000000), new_df[column_heading])

    # save the new dataframe as a separate CSV file in the extracted_data directory with the i-th column heading as the file name
    file_name = f'extracted_data/{column_heading}.csv'
    new_df.to_csv(file_name, index=False)
