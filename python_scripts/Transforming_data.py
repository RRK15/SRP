import pandas as pd
import numpy as np

# read in the count matrix
df=pd.read_csv("orig_matrix.csv")

#dropping the last rows
df.drop(df.tail(36).index, inplace= True)

#loop over each column
for column in df:
    if column == "Gene_Id": #to avoid the gene columns
        continue
    else:
        df[column]=pd.to_numeric(df[column])
        summation=sum(df[column])
        df[column]= np.where(df[column] > 0,  np.log10((df[column] / summation) * 1000000), df[column]) #normalising
        df[column]=df[column].round()
        df[column]=df[column].astype(int)
        df[column]= np.where(df[column] < 0, 0, df[column])

#save new dataframe
file_name="../RNAseq/data/transformed_data.csv"
df.to_csv(file_name, index=False)
