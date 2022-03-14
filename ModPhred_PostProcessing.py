#Import required libraries:
import sys
import os
import pandas as pd
from functools import reduce
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
from venn import venn
import numpy as np
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
import pybedtools
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'

#Start input parser:
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Path to the input mod.gz file.")
parser.add_argument("-o", "--output", help="Output name.")
parser.add_argument("-r", "--reference", help="Reference file (*.fa).")
parser.add_argument("-c", "--conditions", nargs="+", help="Conditions included in the analysis. ie: WT, KO.")
parser.add_argument("-s", "--samples", nargs="+", help="Samples included in the analysis. ie: WT1, WT2, KO1, KO2.")
parser.add_argument("-l", "--labels", nargs="+", help="Group in which the sample is included. ie: 1, 1, 2, 2.")
parser.add_argument("-bed", "--bed_file", default=None, help="Bed file with genes to annotate the replicable m6A sites.")


#Read arguments from the command line
args = parser.parse_args()

##Define accesory funtions:
def InputCheck_Preprocessing(data, labels, samples, conditions):
    num_samples = int((len(data.columns)-5)/4)

    if num_samples==len(labels) and num_samples==len(samples) and len(labels)==len(samples) and len(set(labels))==len(conditions):

        #Rename the columns of the dataframe:
        features = ["_Coverage", "_Accuracy", "_ModFreq", "_MedianProbMod"]
        colnames = list(data.columns[0:5])

        for i in samples:
            for j in features:
                colnames.append(i+j)

        data.columns = colnames

    else:
        sys.exit('Inputs provided are not consistent. Please recheck them and try again.')
        
    return(data)
    
def VennDiagrams(dict_ids, output_file):
    fig = venn(dict_ids, cmap="plasma", fontsize=14, legend_loc="lower right", figsize=(12,12))
    plt.savefig(output_file+".pdf", dpi=300, bbox_inches='tight')   
    plt.close() 
    
def ScatterPlots_withR2(twoSamples_df, output, graph_title, log_scale, output_path):
    sns.set(rc={'figure.figsize':(7,6)})
    sns.set_theme(style="whitegrid")
       
    p = sns.scatterplot(x=twoSamples_df.iloc[:,0], 
                        y=twoSamples_df.iloc[:,1], 
                        #color=(0.6509803921568628, 0.807843137254902, 0.8901960784313725))
                        color="#3262b5")
    
    #Calculate R²:
    r, pv = sp.stats.pearsonr(twoSamples_df.iloc[:,0], twoSamples_df.iloc[:,1])
    ax = plt.gca()
    ax.text(.025, .9, 'r={:.3f}, p={:.2g}'.format(r, pv),
            transform=ax.transAxes)

    #Plot aesthetics:
    p.set(xlabel = twoSamples_df.columns[0], 
          ylabel = twoSamples_df.columns[1])
    
    p.set_title(graph_title)
    
    if log_scale:
        p.set(xscale="log", yscale="log")

        p.xaxis.set_major_formatter(mticker.ScalarFormatter())
        p.yaxis.set_major_formatter(mticker.ScalarFormatter())

        p.set_xlim(0.005,1)
        p.set_ylim(0.005,1)
        
    else:
        p.set_xlim(-0.1,1)
        p.set_ylim(-0.1,1)
    
    p.figure.savefig(output_path+"_Output/Plots/"+twoSamples_df.columns[0]+"_"+twoSamples_df.columns[1]+output, dpi=300)
    plt.close(p.figure)

def ScatterPlot_ChangingSites(data, conditions, output):
    sns.set(rc={'figure.figsize':(7,7)})
    sns.set_theme(style="whitegrid")
    
    #Define the color palette to use:
    if len(set(data["Status"]))==1:
        col_palette = "grey"
    elif len(set(data["Status"]))==2:
        col_palette = ["blue", "grey"]
    else:
        col_palette = ["grey", "blue", "red"]

    #Plotting:    
    xy = sns.scatterplot(data=data, x=data.columns[6], hue = data["Status"],
                        y=data.columns[7],
                        palette=col_palette)
    xy.plot([0,1],[0,1], 'black', linewidth=2, linestyle="dashed")
    xy.set(xlabel = conditions[0], ylabel = conditions[1])
    #xy.set(xscale="log", yscale="log")
    xy.set_ylim(0,0.6)
    xy.set_xlim(0,0.6)
    xy.set_title("Median (% Mod) - Replicable sites with Coverage>=50 in ALL samples")
    xy.xaxis.set_major_formatter(mticker.ScalarFormatter())
    xy.yaxis.set_major_formatter(mticker.ScalarFormatter())
    
    xy.figure.savefig(output+"_Output/Plots/"+conditions[0]+"_"+conditions[1]+"_ChangingStatus_ReplicablePeaks_AllCov_PeakBased.pdf", dpi=300)
    plt.close(xy.figure)
    
def DensityPlots(ModFreq_data, output, file, labels, samples):
    sns.set(rc={'figure.figsize':(11,8)})
    sns.set_theme(style="whitegrid")

    for j in range(0,ModFreq_data.shape[1]):
        c = sns.color_palette("tab10")[int(labels[j])-1]

        if "1_" in ModFreq_data.columns[j]:
            ls = "solid"
        else:
            ls = "dashed"

        #Create label:
        median_sample = np.median(ModFreq_data.iloc[:,j])
        ind_label = samples[j]+" (median = "+"{0:.3f}".format(median_sample)+")"
        ax = sns.kdeplot(data=ModFreq_data, x=(ModFreq_data.iloc[:,j].replace(0, np.nan)), label=ind_label, color = c, 
                         linestyle=ls, log_scale=True, linewidth = 3)

        ax.axvline(median_sample, linestyle=ls, color = c, linewidth = 2)

    ax.legend(loc='upper left', fontsize=15, frameon=False)
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.set_xlabel('m6A modification Frequency', fontsize=16)
    ax.set_ylabel('Density', fontsize=16)
    ax.tick_params(labelsize=14)
    ax.figure.savefig(output+file, dpi=300)
    plt.close()

def Barplots_ReplicableSites(total_sites, mod_AL1, mod_Both, samples_names, output):
    
    #Barplots of the total sites, sites modified in at least one rep, sites modified in all reps:
    x = ['Total sites', 'Modified in at least 1', 'Modified in all reps']
    y = [total_sites, mod_AL1, mod_Both]
    
    a = sns.barplot(x, y)
    
    #Annotate the barplot:
    for bar in a.patches:
        a.annotate(format(bar.get_height(), '.0f'),
                       (bar.get_x() + bar.get_width() / 2,
                        bar.get_height()-2), ha='center', va='center',
                       size=15, xytext=(0, 8),
                       textcoords='offset points')
        
    a.set_title(samples_names[0]+" and "+samples_names[1]+ " - Replicable sites (Cov>=50 and ModFreq>=0.05)")
    a.figure.savefig(output+"_Output/Plots/"+samples_names[0]+"_"+samples_names[1]+"_BarplotsReplicablePeaks_PeakBased.pdf", dpi=300)
    plt.close(a.figure)
    #plt.show()

def mergeLists(list_toMerge, type_merge):
    mergedList = reduce(lambda left, right:
                            pd.merge(left , right,
                            on = ["chr", "pos", "ref_base", "strand", "mod", "Site_ID"],
                            how = type_merge),
                            list_toMerge)
    
    return(mergedList)


def extractColumnData(data, initial_column, non_processed):
    column_list = list()
    n_reps = int((data.shape[1]-5)/4)
    initial = True
    
    if n_reps>1:
        for i in range(0,n_reps):
            if initial: 
                column_list.append(initial_column)
                initial = False
            elif non_processed:
                column_list.append(initial_column+(4*i))
            else:
                column_list.append(initial_column+(4*i)+1)
    else:
        column_list.append(initial_column)
    
    return(column_list)

def motifAnalysisMEME(replicable_sites, reference, output):
    bed_replicable_sites=pd.concat([replicable_sites.iloc[:,1],replicable_sites.iloc[:,2]-10,
                                replicable_sites.iloc[:,2]+10, replicable_sites.iloc[:,0],
                               replicable_sites.iloc[:,5], replicable_sites.iloc[:,4]], axis=1)

    bed_replicable_sites.to_csv("test.bed", sep="\t", index=False, header=False)
    
    print('Running motif analysis with MEME')
    cmd = "bedtools getfasta -fi %s -bed test.bed -s > Kmers_data.fa"%(reference)
    os.system(cmd)
    
    os.system("fasta-get-markov %s %s.ooc"%(reference, reference))
    output_path = output+"_Output"
    
    cmd = "meme -nostatus -dna -mod zoops -nmotifs 5 -minw 2 -maxw 10 -bfile %s.ooc -o %s/MEME %s"%(reference, output_path, "Kmers_data.fa")
    os.system(cmd)
    
    #Remove intermediate files:
    os.system("rm Kmers_data.fa test.bed")

## Main function:
def main():
    ##Load input data: 
    input_data = pd.read_table(args.input, skiprows=16)
    labels = args.labels
    samples = args.samples
    conditions = args.conditions
    output = args.output
    reference = args.reference
    
    data = InputCheck_Preprocessing(input_data, labels, samples, conditions)
    
    ##Create output directories if needed:
    #Main directory:
    if not os.path.exists(args.output+"_Output"):
        os.makedirs(args.output+"_Output")
    
    #Plots directory:
    if not os.path.exists(args.output+"_Output/Plots"):
        os.makedirs(args.output+"_Output/Plots")
    
    #Text files directory:
    if not os.path.exists(args.output+"_Output/Text_files"):
        os.makedirs(args.output+"_Output/Text_files") 
    
    ##COVERAGE BASED ANALYSIS:
    #Extract coverage and mod_freq information:
    coverage = extractColumnData(data, 5, True)
    mod_freq = extractColumnData(data, 7, True)

    #Determine the number of sites with coverage>50 in ALL samples:
    peaks_cov_AllSamples = data.loc[(data.iloc[:,coverage]>=50).all(axis=1)]
    peaks_cov_AllSamples["Site_ID"] = peaks_cov_AllSamples["chr"]+"_"+peaks_cov_AllSamples["pos"].astype("str")+"_"+peaks_cov_AllSamples["strand"]
    print("Total sites with coverage>50 in ALL samples: " + str(peaks_cov_AllSamples.shape[0]))

    #OUTPUT 1: Density plots with data from sites with cov>=50 in ALL samples:
    #ModFreq density plot:
    ModFreq_data = pd.concat([peaks_cov_AllSamples.iloc[:,mod_freq], peaks_cov_AllSamples.iloc[:,-1]], axis=1)
    DensityPlots(ModFreq_data.iloc[:,:-1], args.output+"_Output/Plots/"+output, '_DensityPlots_ModFreq_CoverageBased.pdf', labels, samples)

    #OUTPUT 2: VennDiagrams from modified sites identified in each individual sample + coverage>50 in ALL samples: 
    sites_indSamples = dict()

    for s in range(0,ModFreq_data.shape[1]-1):
        sites = set(ModFreq_data.loc[ModFreq_data.iloc[:,s] >= 0.05,"Site_ID"])
        ind_label = ModFreq_data.columns[s].replace("_ModFreq", "") + " - Modified sites: " + str(len(sites))
        sites_indSamples[ind_label] = sites

    #VennDiagram (input, dataframe with IDs, labels, samples): 
    VennDiagrams(sites_indSamples, args.output+"_Output/Plots/"+output+"_VennDiagram_ModifiedSitesPerSample_AllCoverage_CoverageBased")

    #Print the data to the user:
    for key in sites_indSamples.keys():
        print(key)
        
    #OUTPUT 3: Scatter plots within replicates including data from sites with cov>=50 in ALL samples:
    groups = list()
    for ind_group in labels:
        if ind_group not in groups:
            groups.append(ind_group)

    for n in groups: 

        #Extract replicates data:
        modFreq_reps = list()
        for count,label in enumerate(labels): 
            if label==n:
                modFreq_reps.append(count)

        #Scatter plot ModFreq comparing replicates within the same condition (non-log scale):
        ScatterPlots_withR2(ModFreq_data.iloc[:,modFreq_reps], "_ScatterPlot_ModFreq_CoverageBased.pdf", 
                            "% Mod - Sites with coverage>50 in all samples", False, output)

        #Scatter plot ModFreq comparing replicates within the same condition (log scale):
        ScatterPlots_withR2(ModFreq_data.iloc[:,modFreq_reps], "_ScatterPlot_LogScale_ModFreq_CoverageBased.pdf", 
                            "% Mod - Sites with coverage>50 in all samples", True, output)

    ##PEAK BASED ANALYSIS:
    #OUTPUT 4: VennDiagrams of replicable peaks (cov>=25 + freq>=0.05) - including sites with cov>25 in all reps too
    site_id_idx = peaks_cov_AllSamples.shape[1]-1
    columns = list(range(0,5))
    replicable_per_condition = list()

    for n in groups:
        replicates_data = list()
        samples_names = list()

        #Extract replicates data:
        for count,label in enumerate(labels):  
            if n==label:

                #Extract the columns from a specific sample:
                initial = 5+(count*4)
                final = initial + 4
                columns.extend(range(initial,final))

                #Append to the data from the other replicates (if any):
                columns.append(site_id_idx)
                replicates_data.append(peaks_cov_AllSamples.iloc[:,columns])
                columns = list(range(0,5))

                #Extract sample names:
                samples_names.append(samples[count])

            else:
                continue


        ##Analyse replicates data:
        reps = mergeLists(replicates_data, "outer")

        #Extract coverage columns:
        coverage_rep = extractColumnData(reps, 5, False)
        mod_freq_rep = extractColumnData(reps, 7, False)

        peaks_modfreq_AllSamples_condition = reps.loc[(reps.iloc[:,mod_freq_rep]>=0.05).all(axis=1)]
        replicable_per_condition.append(peaks_modfreq_AllSamples_condition)

        n_cov = reps.shape[0]

        #Identification of sites reported as modified (modfreq>=0.05) in at least one rep:
        n_cov_singleSamples = list()

        for ss in range(0,len(mod_freq_rep)):
            n_cov_singleSamples.append(reps.loc[(reps.iloc[:,mod_freq_rep[ss]]>=0.05)])

        inner_ncov = mergeLists(n_cov_singleSamples, "outer")
        n_freq = peaks_modfreq_AllSamples_condition.shape[0]

        #Barplots of the total sites, sites modified in at least one rep, sites modified in all reps:
        Barplots_ReplicableSites(n_cov, inner_ncov.shape[0], n_freq, samples_names, output)
        
    #OUTPUT 5: VennDiagrams of replicable peaks (cov>=50 + freq>=0.05) across conditions:
    #Plot the VennDiagram of replicable peaks across conditions: 
    replicable_sites_dict = dict()
    for count,condition in enumerate(conditions):
        replicable_sites_dict[condition] = set(replicable_per_condition[count].loc[:,"Site_ID"])

    VennDiagrams(replicable_sites_dict, output+"_VennDiagram_ReplicableSites_AcrossConditions_AllCoverage_PeakBased")

    #OUTPUT 6: Define changing sites across conditions and optional overlap with an annotation file provided by the user:
    replicable_sites_coordinates = mergeLists(replicable_per_condition, "outer").loc[:,"Site_ID"]
    replicable_sites = pd.merge(replicable_sites_coordinates, peaks_cov_AllSamples,
                                on = ["Site_ID"],
                                how = "inner")

    #Calculate ModFreq median:
    for j in range(0, len(conditions)):           
        replicable_sites[conditions[j]+"_medianModFreq"] = replicable_sites.filter(like=conditions[j]).filter(like='_ModFreq').median(axis=1)

    #Calculate ModFreq differences and ratios between conditions:
    replicable_sites["ModFreq(Cond1-Cond2)"] = replicable_sites.filter(like="_medianModFreq").iloc[:,0] - replicable_sites.filter(like="_medianModFreq").iloc[:,1]
    replicable_sites["Ratio(Cond1/Cond2)"] = replicable_sites.filter(like="_medianModFreq").iloc[:,0]/replicable_sites.filter(like="_medianModFreq").iloc[:,1]

    #Determine directionalities:
    status = list()
    for index,row in replicable_sites.iterrows():
        diff_ModFreq = row["ModFreq(Cond1-Cond2)"]
        ratio_ModFreq = row["Ratio(Cond1/Cond2)"]

        if (diff_ModFreq>=0.2 or ratio_ModFreq>=2):
            status.append("Decreased upon "+conditions[1])
        elif (diff_ModFreq<=(-0.2) or ratio_ModFreq<=0.5):
            status.append("Increased upon "+conditions[1])
        else:
            status.append("No changes")

    replicable_sites["Status"] = status

    #Optional annotation of the replicable m6A sites with bed file provided by the user:
    
    if args.bed_file is None:
        #Report replicable sites - all data:
        replicable_sites.to_csv(output+"_Output/Text_files/"+output+"_RawData_ReplicablePeaks_AllCoverage_PeakBased.tsv", sep="\t", index=False)

        #Report replicable sites - summary of analysis:
        replicable_sites_processed = pd.concat([replicable_sites.iloc[:,0:6], replicable_sites.iloc[:,-5:]], axis=1)
        replicable_sites_processed.to_csv(output+"_Output/Text_files/"+output+"_SummaryData_ReplicablePeaks_AllCoverage_PeakBased.tsv", sep="\t", index=False)
        
    else:
        bed_format = pd.concat([replicable_sites.iloc[:,1], replicable_sites.iloc[:,2], replicable_sites.iloc[:,2]+1, 
                       replicable_sites.iloc[:,0], pd.Series([0] * len(replicable_sites)),
                       replicable_sites.iloc[:,4]], axis=1)
        
        a = pybedtools.BedTool.from_dataframe(bed_format)
        b = pybedtools.BedTool(args.bed_file)

        #Intersection: 
        intersection = a.intersect(b, wb=True).to_dataframe().iloc[:,[3,9]]
        intersection.rename(columns={'name': 'Site_ID', 'blockCount': 'Gene'}, inplace=True)

        #Merge:
        replicable_sites = pd.merge(replicable_sites , intersection,
                                    on = ["Site_ID"],
                                    how = "outer")    
    
        #Report replicable sites - all data:
        replicable_sites.to_csv(output+"_Output/Text_files/"+output+"_RawData_ReplicablePeaks_AllCoverage_PeakBased.tsv", sep="\t", index=False)

        #Report replicable sites - summary of analysis:
        replicable_sites_processed = pd.concat([replicable_sites.iloc[:,0:6], replicable_sites.iloc[:,-6:]], axis=1)
        replicable_sites_processed.to_csv(output+"_Output/Text_files/"+output+"_SummaryData_ReplicablePeaks_AllCoverage_PeakBased.tsv", sep="\t", index=False)

    
    # OUTPUT 7: Scatter plot with changing sites (only sites with enough cov in all samples):
    ScatterPlot_ChangingSites(replicable_sites_processed, conditions, output)

    ##OUTPUT 8: Motiff analysis with MEME:
    motifAnalysisMEME(replicable_sites, reference, output)
    
if __name__ == "__main__":
    main()



