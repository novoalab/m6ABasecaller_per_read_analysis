#########################################################################
### Script to filter out non-replicable m6A sites from per-read data. ###
#########################################################################

#Import the required libraries: 
import pandas as pd
import os
import argparse

#Start input parser:
parser = argparse.ArgumentParser()

parser.add_argument("-pr", "--per_read", help="Path to the input file containing per read modification data.")
parser.add_argument("-rs", "--replicable_sites", help="Path to the input file containing replicable sites (*_SummaryData_ReplicableSites.tsv).")
parser.add_argument("-o", "--output", help="Output name.")

#Read arguments from the command line
args = parser.parse_args()

#Import and parse input data:
m6a_sites = pd.read_table(args.per_read, sep=" ")
replicable_sites = pd.read_table(args.replicable_sites, sep="\t", usecols=[1,2])

#Loop over all chromosomes:
for chr_ref in m6a_sites['ref_name'].unique():
    
    #Subset m6a sites and replicable sites:
    subset_m6a = m6a_sites.loc[m6a_sites['ref_name'] == chr_ref]
    subset_replicable = replicable_sites.loc[replicable_sites['chr'] == chr_ref].iloc[:,1].unique()
    
    #Loop over all m6a sites:
    for i in range(0,len(subset_m6a)):
        value = subset_m6a.iloc[i,9]
        
        if pd.isnull(value):
            continue
        
        else:
            #Check if m6A sites from individual reads are in the list of replicable sites:
            value = value.split(",")
            filtered = list()
            
            for idx,j in enumerate(value):
                if int(j) in subset_replicable:
                    filtered.append(j)
                else:
                    continue
            
            #Update table with the filtered sites:
            if len(filtered)==0:
                subset_m6a.iloc[i,9] = float("nan")
                subset_m6a.iloc[i,1] = 0
            else:   
                subset_m6a.iloc[i,9] = ','.join(filtered)
                subset_m6a.iloc[i,1] = len(filtered)
    
    #Save to output file:
    output_file = 'Per_read_' + args.output + '_ReplicableSites.tsv'
    subset_m6a.to_csv(output_file, sep = '\t', index = False, mode = 'a', header = not os.path.exists(output_file))
