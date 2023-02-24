#############################################################################
### Script to analyse the co-occurance of m6A sites from read level data. ###
#############################################################################

#Import the required libraries: 
import pandas as pd
import numpy as np
import itertools
import math
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict

#Additional functions:
def parse_data (input_data,min_coverage):
    #Only include reads that were uniquely assigned:
    uniquely_assigned = input_data.loc[input_data.iloc[:,4]=="unique"]
    
    #Extract transcript IDs with coverage>=args.coverage:
    transcript_counts = uniquely_assigned['isoform_id'].value_counts()
    transcript_coverage = transcript_counts[transcript_counts>=min_coverage].index.tolist()

    #Print summary to stout:
    print('-- Transcript selection summary --')
    print('Total number of transcripts with uniquely assigned reads: {}'.format(len(transcript_counts)))
    print('Total number of transcripts with uniquely assigned reads>={}: {}'.format(min_coverage, len(transcript_coverage)))

    return(transcript_coverage)

def calculate_frequencies(transcript_data,transcript_id,min_expected,total_counts):
    
    #Extract all individual m6A positions: 
    sites = transcript_data.str.split(',',expand=True)
    
    m6A_sites = list()
    for col in range(0, sites.shape[1]):
        m6A_sites.append(sites.iloc[:,col].dropna().unique().tolist())

    m6A_sites = set(list(itertools.chain.from_iterable(m6A_sites))) 

    #Calculate the frequency (stoichiometries) for each individual site: 
    counts = defaultdict(lambda:1)

    #Loop over and get the site counts:
    for site in m6A_sites:
        for line in range(0,sites.shape[0]):
            sites_per_read = sites.iloc[line,:].tolist()
            
            #Update dictionary for individual frequencies (pA and pB):
            if site in sites_per_read:
                counts[site] +=1

    '''
    #Get the frequency:
    for key in list(counts.keys()):
        freq = counts.get(key)/sites.shape[0]

        #If lower than the threshold remove site from dictionary:
        if freq<min_freq:
            del counts[key]
        else:
            continue
    
    #Update the m6A sites list - only including those with enough stoichiometry:
    m6A_sites = list(counts.keys())
    '''

    #Only continue if there are 2 or more m6A sites in the transcript:
    if len(m6A_sites)>1:
        #Generate all possible pair-wise comparisons:
        pairwise_comparisons = list(itertools.combinations(m6A_sites, 2))
        #print('{} pairwise comparisons found:'.format(len(pairwise_comparisons)))

        #Analyse each pair-wise comparison:
        final_data = ""
        for comp in pairwise_comparisons:
            comparison = '{}-{}'.format(comp[0], comp[1])
            
            #Calculate frequency for individual sites:
            freq_A = counts[comp[0]]/total_counts
            freq_B = counts[comp[1]]/total_counts

            #Apply min_expected counts threshold:
            counts_expected = freq_A*freq_B*total_counts
            
            if counts_expected >= min_expected:
                print('Processing: '+ comparison)

                #Calculate observed frequency (pAB):
                for line in range(0,sites.shape[0]):
                    sites_per_read = sites.iloc[line,:].tolist()

                    #Update dictionary with observed frequency for pAB:
                    if comp[0] in sites_per_read and comp[1] in sites_per_read and comparison not in counts:
                        counts[comparison] = 1
                    elif comp[0] in sites_per_read and comp[1] in sites_per_read and comparison in counts:
                        counts[comparison] +=1

                #If there are no AB counts, update the dictionary anyway: 
                if comparison not in counts.keys():
                    counts[comparison] = 0
                
                #Calculate standard deviations from expected
                counts_expected = freq_A*freq_B*sites.shape[0]
                sd_from_expected = (counts[comparison] - counts_expected)/math.sqrt(counts_expected*(1-(freq_A*freq_B)))

                #Update final dataframe:
                to_update = pd.DataFrame([[transcript_id, comp[0], comp[1], sites.shape[0], freq_A, freq_B, counts[comparison]/sites.shape[0], freq_A*freq_B, sd_from_expected]])
                if type(final_data)==str:
                    final_data = to_update
                else: 
                    final_data = pd.concat([final_data, to_update])

            else:
                continue
    else:
        final_data = ""
        print('Not enough m6A sites to perform pairwise comparisons.')

    return(final_data)

#Start input parser:
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Path to the input file containing per read label data.")
parser.add_argument("-o", "--output", help="Output name.")
parser.add_argument("-cov", "--coverage", default=200, type=int, help="Minimum coverage to include a transcript in the analysis. Default = 200")
parser.add_argument("-min_exp", "--min_expected", default=2, type=int, help="Minimum expected counts to include a pairwise comparison in the analysis. Default = 2")

#Read arguments from the command line
args = parser.parse_args()

#Read the data:
input_data = pd.read_table(args.input, usecols=['read_id', 'chr', 'pos_m6A_sites', 'isoform_id', 'assignment_type'])
transcripts_to_analyse = parse_data(input_data, args.coverage)

#Extract reads from transcripts with coverage>=args.coverage:
coocurance_data = ""
for ind_transcript in transcripts_to_analyse:
    
    #Perform calculations on reads belonging to the same transcript:
    print('-- Analysing transcript: {} ---'.format(ind_transcript))
    
    #Subset data per transcript:
    transcript_subset = input_data.loc[input_data['isoform_id']==ind_transcript,'pos_m6A_sites']
    #print(transcript_subset.dropna(axis = 0, how = 'all'))
    if transcript_subset is not None:
        results = calculate_frequencies(transcript_subset.dropna(axis = 0, how = 'all'),ind_transcript, args.min_expected, transcript_subset.shape[0])

        if type(results)==str:
            continue
        elif type(coocurance_data)==str and type(results)!=str:
            coocurance_data = results
        elif type(coocurance_data)!=str and type(results)!=str:
            coocurance_data = pd.concat([coocurance_data,results])
    else:
        print('Transcript without m6A sites.')

coocurance_data.columns = ['Transcript','Site A','Site B', 'Transcript_counts','FreqA','FreqB', 'obs_FreqAB', 'exp_FreqAB','SdFromExpected']

#Output results in a tsv file:
coocurance_data.to_csv('{}_CoocuranceAnalysis.tsv'.format(args.output), sep='\t', index=False)

