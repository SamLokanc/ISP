#----- imports -----
import pandas as pd
import os
import re
import sys

#----- functions -----
def generate_mutation_pattern(mutation):
    mutation_patterns = [mutation]

    #add patterns for mutation at different positions
    for i in range(1, len(mutation)):
        mutation_patterns.append(mutation[:i] + '-' + mutation[i:])

    mutation_patterns.append(mutation[1:])
    mutation_patterns.append(mutation[:len(mutation) - 1])
    
    return '|'.join(map(re.escape, mutation_patterns))

#----- create results directory -----
results_dir = 'results'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

#----- global variables from command line -----
if len(sys.argv) != 6:
    print("Usage: python3 MHC_script.py MHC_I_out_file_path MHC_II_out_file_path output_file.csv number_of_top_binders mutation_with_flanking_residues\n")
    sys.exit(1)

MHC_I_out_file_path = sys.argv[1]
MHC_II_out_file_path = sys.argv[2]
output_file_path = os.path.join(results_dir, sys.argv[3])
n =int(sys.argv[4])
mutation = sys.argv[5]

#----- generate mutation pattern -----
mutation_pattern = generate_mutation_pattern(mutation)

#----- overwrite out file -----
with open(output_file_path, 'w'):
    pass

#----- reading MHC out files into dataframes -----
#get HLA types, store in list in order of appearance
MHC_I = pd.read_csv(MHC_I_out_file_path, delimiter = '\t')
MHC_II = pd.read_csv(MHC_II_out_file_path, delimiter = '\t')
HLA_types_I = [col for col in MHC_I.columns if not col.startswith('Unnamed')]
HLA_types_II = [col for col in MHC_II.columns if not col.startswith('Unnamed')]
HLA_types = HLA_types_I + HLA_types_II

#read MHC out files into dataframes, skip first row for formatting
MHC_I = pd.read_csv(MHC_I_out_file_path, delimiter = '\t', skiprows = 1, header = None)
MHC_II = pd.read_csv(MHC_II_out_file_path, delimiter = '\t', skiprows = 1, header = None)

#remove target column from MHC_II dataframe
MHC_II.drop(MHC_II.columns[3], axis = 1, inplace = True)

#create dataframes array for later
dataframes = [MHC_I, MHC_II]

#----- subsetting dataframes into HLA types -----
#create array to store subsets
subsets = []

#loop thru MHCI and II dfs 
for df in dataframes:
    #remove last two columns
    df = df.iloc[:, :-2]

    #specify total columns + additional columns for each subset
    total_cols = df.shape[1]
    additional_cols_per_subset = 6

    #loop thru MHC_I df and create subsets for each HLA type
    for i in range(0, total_cols - 3, additional_cols_per_subset):
        #calculate end index for each subset
        subset_end = min(i + 3 + additional_cols_per_subset, total_cols)

        #create a subset containing the first 3 columns and additional columns
        subset = df.iloc[:, :3].join(df.iloc[:, i + 3:subset_end])

        #append the subset to the list
        subsets.append(subset)

#----- retreive best predicted binding peptides that contain mutation of interest -----
#loop thru subsets
for i, subset in enumerate(subsets):
    #ensure subsets from MHCI and II dfs have consistent format
    subset = subset.copy()
    subset.columns = subset.iloc[0].copy()
    subset = subset.iloc[1:].copy()
    subset.rename(columns = {'Rank' : 'EL_Rank', 'Rank_BA' : 'BA_Rank', 'Core' : 'core'}, inplace = True)
    subset.loc[:, 'EL_Rank'] = pd.to_numeric(subset['EL_Rank'], errors='coerce')
    subset.loc[:, 'BA_Rank'] = pd.to_numeric(subset['BA_Rank'], errors='coerce')

    #remove rows corresponding to peptides that do not contain the mutation of interest
    subset = subset[subset['core'].str.contains(mutation_pattern)]

    #create column that averages EL and BA ranks
    subset.loc[:, 'AVE_Rank'] = (subset['EL_Rank'] + subset['BA_Rank']) / 2

    #retrieve best binders by looking for lowest average rank in each subset
    subset['AVE_Rank'] = pd.to_numeric(subset['AVE_Rank'], errors='coerce')
    best_binders = subset.nsmallest(n, 'AVE_Rank')

    #write information to file
    with open(output_file_path, 'a') as file:
        file.write(f"{HLA_types[i]}:\n")
    
    best_binders[['Peptide', 'core','EL_Rank','BA_Rank','AVE_Rank']].to_csv(output_file_path, mode = 'a', header = not i, index = False)