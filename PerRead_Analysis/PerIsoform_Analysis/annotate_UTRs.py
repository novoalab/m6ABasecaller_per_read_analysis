#######################################################################
### Script to annotate the 5' and 3' UTRs from per read level data. ###
#######################################################################

#Import the required libraries: 
import pandas as pd
import itertools
import argparse
import statistics
import os
import warnings
warnings.filterwarnings("ignore")


#Additional functions:
def parse_data (input_data,min_coverage):
    #Only include reads that were uniquely assigned and full-length:
    uniquely_assigned = input_data.loc[(input_data.iloc[:,14]=="unique") & ((input_data.iloc[:,15]=="fsm") | (input_data.iloc[:,15]=="mono_exon_match"))]

    #Extract transcript IDs with coverage>=args.coverage:
    transcript_counts = uniquely_assigned['isoform_id'].value_counts()
    transcript_coverage = transcript_counts[transcript_counts>=min_coverage].index.tolist()

    #Print summary to stout:
    print('-- Transcript selection summary --')
    print('Total number of transcripts with uniquely assigned and full-length reads: {}'.format(len(transcript_counts)))
    print('Total number of transcripts with uniquely assigned and full-length reads>={}: {}'.format(min_coverage, len(transcript_coverage)))

    return(transcript_coverage, uniquely_assigned)

def extract_starts_ends(transcript_subset, start):

    #Subset reads based on their start/end coordinates:
    if start:
        coordinates = transcript_subset.loc[:,'start_aln'].value_counts()
    else:
        coordinates = transcript_subset.loc[:,'end_aln'].value_counts()

    #Bin all starts/ends represented by a significant number of reads (>=10% of the total number of reads belonging to that transcript):
    threshold_coordinates = transcript_subset.shape[0]*0.085
    valid = sorted(coordinates[coordinates>=threshold_coordinates].index.tolist())
    
    #If none of the coordinates meets the threshold requirement, take the most frequent position: 
    if not valid:
        valid = [coordinates.index.tolist()[0]]

    #Binning:
    bins = dict()
    bin = 0

    #Initial site:
    bins[bin] = [str(valid[0])]

    #Looping through the remaining sites:
    for i in range(1,len(valid)):
        diff = valid[i]-valid[i-1]

        if diff>=50:
            bin += 1
            bins[bin] = [str(valid[i])]
        else:
            bins[bin].append(str(valid[i]))
    
    #Calculate the median position per bin:
    median_coordinates = list()
    for key in bins.keys():
        s = statistics.median(map(int, bins.get(key)))
        median_coordinates.append(round(s))
    
    return(median_coordinates)

def start_end_permutations(starts,ends):
    #Permutation of all possible starts/ends:
    if len(starts)<len(ends):
        permut = itertools.permutations(ends, len(starts))
    else:
        permut = itertools.permutations(starts, len(ends))
        
    
    unique_combinations = list()
    
    for comb in permut:
        if len(starts)<len(ends):
            zipped = zip(starts, comb)
        else:
            zipped = zip(comb, ends)
        unique_combinations.append(list(zipped))

    #Reformat list with unique combinations:
    unique_combinations = list(itertools.chain.from_iterable(unique_combinations))

    #Return unique_combination list:
    return(unique_combinations)

#Start input parser:
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Path to the input file containing per read label data.")
parser.add_argument("-o", "--output", help="Output name.")
parser.add_argument("-cov", "--coverage", default=40, type=int, help="Minimum coverage to include a transcript in the analysis. Default = 40")

#Read arguments from the command line
args = parser.parse_args()

#Read the data:
input_data = pd.read_table(args.input)
transcripts_to_analyse, parsed_data = parse_data(input_data, args.coverage)

#Extract reads from transcripts with coverage>=args.coverage:
coocurance_data = ""
for ind_transcript in transcripts_to_analyse:
    #ind_transcript = "ENST00000395699"
    #Perform calculations on reads belonging to the same transcript:
    print('-- Analysing transcript: {} ---'.format(ind_transcript))
    
    #Subset data per transcript:
    transcript_subset = parsed_data.loc[parsed_data['isoform_id']==ind_transcript, ['read_id', 'start_aln', 'end_aln', 'assignment_type', 'assignment_data']]
    transcript_subset['start_aln'] = pd.to_numeric(transcript_subset['start_aln'])
    transcript_subset['end_aln'] = pd.to_numeric(transcript_subset['end_aln'])
    
    #Extract significant start/end coordinates:
    starts = extract_starts_ends(transcript_subset, True)
    ends = extract_starts_ends(transcript_subset, False)

    #Define all posible start/end combinations:
    combinations = start_end_permutations(starts,ends)
    
    #Extract reads with specific start/end of alignment:
    for i in range(0,len(combinations)):
        ind_start = combinations[i][0]
        ind_end = combinations[i][1]

        annotation_subset = transcript_subset.loc[(transcript_subset['start_aln']<=(ind_start+25)) & (transcript_subset['start_aln']>=(ind_start-25)) & (transcript_subset['end_aln']<=(ind_end+25)) & (transcript_subset['end_aln']>=(ind_end-25)),'read_id']
        
        #If none of the reads fill in the coordinates, skip them:
        if annotation_subset.empty:
            print('No reads match this coordinate set.')
            continue

        #Define updated transcript id:
        updated_transcript_id = ind_transcript+'_'+str(i)

        #Extract reads from the initial table and update their annotation:
        to_export = parsed_data[parsed_data['read_id'].isin(annotation_subset)]
        to_export.loc[:,'upd_transcript_id'] = updated_transcript_id
        to_export.loc[:,'upd_transcript_start'] = ind_start
        to_export.loc[:,'upd_transcript_end'] = ind_end

        #Save re-annotated reads into a table:
        output_path = 'per_read_reannotated_unique_fsm_monoex_'+args.output+'.tsv' 
        to_export.to_csv(output_path, sep="\t", index=False, mode='a', header=not os.path.exists(output_path))

    


    
    
    
