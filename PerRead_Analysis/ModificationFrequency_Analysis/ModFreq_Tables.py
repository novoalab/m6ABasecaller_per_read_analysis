###########################################################################################################
### Script to generate tables containing the modification frequency data at per gene and isoform level. ###
###########################################################################################################

#Import the required libraries: 
import pandas as pd
import os
import argparse
import pybedtools
import re

#Additional functions:
def create_dict_m6A(transcript, given_sites, m6a_sites_coverage):
    #Initialize dictionary to calculate frequencies:
    isoform_dict = dict()

    #Parse data per isoform:
    for k in range(0, transcript.shape[0]):
        if pd.isnull(transcript.iloc[k,1]):
            continue

        else:
            sites = str(transcript.iloc[k,1]).split(',')
    
            #Update the dictionary:
            for site in sites:
                site = int(float(site))
                if site not in isoform_dict:
                    isoform_dict[site] = 1
                else:
                    isoform_dict[site] +=1
                    
    #If there are given m6A sites, update dict in case of any missing site with coverage but with no m6a:
    if given_sites and m6a_sites_coverage is not None:
        for m6a in m6a_sites_coverage:
            if m6a not in isoform_dict.keys():
                isoform_dict[m6a] = 0
        print(isoform_dict)
    
    return(isoform_dict)

#Start input parser:
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Path to the input file containing reannotated per read data (output from annotate_UTRs.py).")
parser.add_argument("-gtf", "--gtf", help="Path to the annotation file - ideally only containing exon data (*.gtf).")
parser.add_argument("-o", "--output", help="Output name.")

#Read arguments from the command line
args = parser.parse_args()

#Read input data: 
data = pd.read_table(args.input, usecols=[2,9,12,13,16,17,18,19])

#Read input gtf file:
gtf_file = args.gtf
a = pybedtools.BedTool(gtf_file)

#Parse data at per gene level:
genes = data['gene_id'].unique()

#Get sample information:
samples = data['sample'].unique()

for sample in samples:
    data_sample = data.loc[data['sample'] == sample]

    for gene in genes:
        #Subset by gene:
        subset = data_sample.loc[data['gene_id'] == gene]
        
        #Extract all isoforms ID and coordinates present in the gene: 
        isoforms = subset['upd_transcript_id'].unique()
        isoforms_data = subset.iloc[:,[0,2,5,6,7]].drop_duplicates()
        
        #Create dictionary to calculate frequencies at per gene level:
        gene_dict = create_dict_m6A(subset, False, '')
        
        #Calculate modFreq per gene:
        total_reads_gene = subset.shape[0]
        for key in gene_dict.keys():
            modFreq = gene_dict.get(key)/total_reads_gene
            #Create line to export:
            line_gene = pd.DataFrame([[sample, gene, total_reads_gene, key, modFreq]], 
                                columns=['Sample', 'GeneID', 'Coverage', 'm6A_Site', 'ModFreq'])
            #Save data:
            output_path = args.output + '_ModFreq_PerGene.tsv'
            line_gene.to_csv(output_path, sep='\t', mode='a', index=False, header=not os.path.exists(output_path))

        #If the gene does not contain any m6A site - skip:
        if len(gene_dict) == 0:
            continue
        
        #For every isoform, check if it has coverage in the m6A sites:
        isoform_sites = dict()
        
        for isoform in isoforms:
            
            #Extract start/end coordinates of the isoform:
            start_isoform = int(float(isoforms_data.loc[isoforms_data['upd_transcript_id'] == isoform].iloc[:,3]))
            end_isoform = int(float(isoforms_data.loc[isoforms_data['upd_transcript_id'] == isoform].iloc[:,4]))
            chr_isoform = isoforms_data.loc[isoforms_data['upd_transcript_id'] == isoform].iloc[:,0].to_list()[0]
            isoform_isoquant = isoforms_data.loc[isoforms_data['upd_transcript_id'] == isoform].iloc[:,1].to_list()[0]

            #Analyse each m6A site:
            for key in gene_dict.keys():   
                
                #Check 1: is the site within start/end coordinates of the isoform?
                key = int(float(key))
                if (key>=start_isoform and key<=end_isoform) or (key>=end_isoform and key<=start_isoform):
                    
                    #Check 2: is the site in a region covered by this isoform?
                    #Create bed file containing the site:
                    bed_site = pybedtools.BedTool(chr_isoform + " " + str(key) + " " + str(key+1), from_string=True)

                    #Create intersection:
                    intersection = bed_site.intersect(a, wb = True)
                    
                    #Analyse the intersection:
                    for i in intersection:
                        regex_result = re.findall(isoform_isoquant, str(i).split('\t')[11], flags=re.IGNORECASE)

                        if len(regex_result)>0:
                            #Add to dict:
                            if isoform not in isoform_sites.keys():
                                isoform_sites[isoform] = [key]
                            else:
                                isoform_sites.get(isoform).append(key)
                        
                            #Once that we found a match, go to the next site:
                            break
                        
                        else:
                            #Site is discarded for this isoform as it doesnt have coverage
                            continue

                else:
                    continue
                    
            #Create dict to calculate frequencies in every site in the isoform:
            transcript = subset.loc[subset['upd_transcript_id'] == isoform]
            ind_isoform_dict = create_dict_m6A(transcript, True, isoform_sites.get(isoform))
            
            #Calculate modFreq per isoform:
            total_reads_isoform = transcript.shape[0]
            
            for key in ind_isoform_dict.keys():
                modFreq = ind_isoform_dict.get(key)/total_reads_isoform
                
                #Create line to export:
                line = pd.DataFrame([[sample, gene, isoform, total_reads_isoform, key, modFreq]], 
                                    columns=['Sample', 'GeneID', 'IsoformID', 'Coverage', 'm6A_Site', 'ModFreq'])
                
                #Save data:
                output_path = args.output + '_ModFreq_PerIsoform.tsv'
                line.to_csv(output_path, sep='\t', mode='a', index=False, header=not os.path.exists(output_path))


'''
### OLD VERSION
for gene in genes:
    #Subset by gene:
    subset = data.loc[data['gene_id'] == gene]
    
    #Create list to store dictionaries per gene:
    gene_dict = list()
    
    #Subset gene data based on isoforms:
    isoforms = subset['upd_transcript_id'].unique()
    
    for isoform in isoforms:
        transcript = subset.loc[subset['upd_transcript_id'] == isoform]
          
        #Create dictionary to calculate frequencies:
        isoform_dict = create_dict_m6A(transcript)
        
        #Only proceed if there is any m6A site present along the feature:
        if len(isoform_dict)==0:
            continue
            
        #Calculate modFreq per isoform:
        total_reads_isoform = transcript.shape[0]
        
        for key in isoform_dict.keys():
            modFreq = isoform_dict.get(key)/total_reads_isoform
            
            #Create line to export:
            line = pd.DataFrame([[sample, gene, isoform, total_reads_isoform, key, modFreq]], 
                                columns=['Sample', 'GeneID', 'IsoformID', 'Coverage', 'm6A_Site', 'ModFreq'])
            
            #Save data:
            output_path = args.output + '_ModFreq_PerIsoform.tsv'
            line.to_csv(output_path, sep='\t', mode='a', index=False, header=not os.path.exists(output_path))
        
        #Update dictionary list:
        gene_dict.append(isoform_dict)
        
    #Merge dictionaries in list:
    initial_merge = True
    for i in range(1,len(gene_dict)):
        
        if initial_merge:
            merged_dict = {k: gene_dict[0].get(k, 0) + gene_dict[i].get(k, 0) for k in set(gene_dict[0]) | set(gene_dict[i])}
            initial_merge = False
        else:
            merged_dict = {k: merged_dict.get(k, 0) + gene_dict[i].get(k, 0) for k in set(merged_dict) | set(gene_dict[i])}
    
    #Calculate modFreq per gene:
    total_reads_gene = subset.shape[0]
    
    if len(gene_dict) != 0:
        for key in merged_dict.keys():
                modFreq = merged_dict.get(key)/total_reads_gene
                
                #Create line to export:
                line_gene = pd.DataFrame([[sample, gene, total_reads_gene, key, modFreq]], 
                                    columns=['Sample', 'GeneID', 'Coverage', 'm6A_Site', 'ModFreq'])
                
                #Save data:
                output_path = args.output + '_ModFreq_PerGene.tsv'
                line_gene.to_csv(output_path, sep='\t', mode='a', index=False, header=not os.path.exists(output_path))

'''
    