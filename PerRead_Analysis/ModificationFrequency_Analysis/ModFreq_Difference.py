##############################################################################
### Script to calculate modFreq differences at per gene and isoform level. ###
##############################################################################

#Import required libraries:
import pandas as pd
import os
import seaborn as sns
import numpy as np
import scipy as sp
import argparse

#Additional functions:
def calculate_modFreq_differences(filtered_data, plotting_sites):
    
    initial_df = True
    for site in plotting_sites:
        subset = filtered_data.loc[filtered_data['Index'] == site]
        
        #Get sample, gene and m6A site information:
        sample = subset['Sample'].unique()[0]
        gene = subset['GeneID'].unique()[0]
        m6A_site = int(float(subset['m6A_Site'].unique()[0]))
        
        #Loop through the modFreq from the different isoforms:
        for i in range(0,subset.shape[0]):
            Isoform1 = subset.iloc[i,2]
            modFreq_Isoform1 = subset.iloc[i,4]
            
            for j in range(i+1,subset.shape[0]):
                Isoform2 = subset.iloc[j,2]
                modFreq_Isoform2 = subset.iloc[j,4]
                
                #Calculate diffModFreq(I1-I2):
                diff = modFreq_Isoform1 - modFreq_Isoform2
                
                #Data to export into a file: 
                line = pd.DataFrame([[sample, gene, Isoform1, Isoform2, m6A_site, modFreq_Isoform1, modFreq_Isoform2, diff]], 
                                    columns=['Sample', 'GeneID', 'Isoform_1', 'Isoform_2', 'Site', 'ModFreq_I1', 'ModFreq_I2', 'ΔModFreq(I1-I2)'])
                
                #Save data to a text file:
                output_path='DiffModFreq_PerIsoform.tsv'
                line.to_csv(output_path, sep='\t', mode='a', index=False, header=not os.path.exists(output_path))
                
                #Also, create a dataframe for plotting purposes:
                if initial_df:
                    plotting = line
                    initial_df = False
                    
                else:
                    plotting = pd.concat([plotting, line], ignore_index=True, axis=0)
    
    return(plotting)

#Start input parser:
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Path to the input file containing modification frequencies at per isoform or gene level.")
parser.add_argument("-cov", "--coverage", default=40, type=int, help="Minimum coverage to include an isoform/gene in the analysis. Default = 40")
parser.add_argument("-o", "--output", help="Output name.")

#Read arguments from the command line
args = parser.parse_args()

#Read input data: 
data = pd.read_table(args.input, sep='\t')
filtered_data = data.loc[data['Coverage']>=args.coverage, ['Sample', 'GeneID', 'IsoformID', 'm6A_Site', 'ModFreq']]

#Create index (GeneID, m6A_Site):
filtered_data['Index'] = filtered_data['GeneID'] + '_' + filtered_data['m6A_Site'].astype(str)

#Extract sites present in more than one isoform:
unique_sites = filtered_data['Index'].value_counts().to_dict()
plotting_sites = list()

for key in unique_sites.keys():
    if unique_sites[key]>=2:
        plotting_sites.append(key)

#Calculate modFreq differences at per gene and isoform level:
plotting = calculate_modFreq_differences(filtered_data, plotting_sites)

#Plotting scatter plot with density:
#Density:
if plotting.shape[0]>=5:
    values = np.vstack([plotting.iloc[:,5], plotting.iloc[:,6]])
    kernel = sp.stats.gaussian_kde(values)(values)

    #Scatterplot:
    xy = sns.scatterplot(data=plotting, x=plotting.columns[5],
                        y=plotting.columns[6],
                        c=kernel,
                        cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True))
else:
    #Scatterplot:
    xy = sns.scatterplot(data=data, x=data.columns[0],
                        y=data.columns[1],
                        cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True))

#Plot differentially modified sites in red:
diff_sites = plotting.loc[abs(plotting.iloc[:,7])>=0.1]
xy = sns.scatterplot(data=diff_sites, x=diff_sites.columns[5],
                        y=diff_sites.columns[6], color='r')

#Add headers and labels:
xy.plot([0,1],[0,1], 'black', linewidth=2, linestyle="dashed")
xy.set(xlabel = 'Isoform 1', ylabel = 'Isoform 2')
xy.set_ylim(0,0.8)
xy.set_xlim(0,0.8)
xy.set_title('ModFreq difference between isoforms')

#Save figure into pdf:
outname = "ModFreq_Differences_PerIsoform.pdf"
xy.figure.savefig(outname, bbox_inches='tight')      