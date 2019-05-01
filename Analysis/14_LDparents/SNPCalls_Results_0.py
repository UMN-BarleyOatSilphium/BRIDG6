# Import all required packages
import pandas as pd
import numpy as np
import multiprocessing as mp
import csv

# read in the dataframe
df = pd.read_csv('SNPCalls_Imputed.csv', index_col=0)

# Generate all unique combinations of columns and store as a list in pair
pair = [[i, j] for i in range(len(df.columns)-1) for j in range(i+1, len(df.columns))]

# Split the list subarrays
pair_split = np.array_split(pair, 10)

# Define the function that compares elements of the two columns and calculates percentage difference in the columns
# Fix the dataframe argument to take in the df of interest and let combination remain a variable to iterate over
def Perc_diff(combn, df = df):
    pair = df.iloc[:, combn]
    pairNoNA = pair.dropna()
    TotalSites = pairNoNA.shape[0]
    DiffSites = np.sum(pairNoNA[pairNoNA.columns[0]] != pairNoNA[pairNoNA.columns[1]])
    PercDiff = round((DiffSites * 100)/ TotalSites, 2)
    # Generate a list as output
    result = [pair.columns[0], pair.columns[1], TotalSites, DiffSites, PercDiff]
    
    return result

#Using multiprocessing for parallel computation, instantiate pool and assign the total number of cores available using mp.cpu_count()
pool = mp.Pool(mp.cpu_count())
# Apply the funtion over one subarray of pair_split at a time usine pool.map as follows
results_0 = pool.map(Perc_diff, [combn for combn in pair_split[0]])
# Do not forget to close the pool
pool.close()
# Convert the results list into a dataframe
Results_0 = pd.DataFrame(results_0, columns = ['Sample1', 'Sample2', 'Total_sites', 'Diff_Sites', 'Percent_Diff'])
# Save the results dataframe as a csv file
Results_0.to_csv('Results_0.csv')

