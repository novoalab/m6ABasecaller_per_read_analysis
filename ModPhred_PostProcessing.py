#Import required libraries:
import sys
import pandas as pd
from functools import reduce
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib_venn import venn2
import numpy as np
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
import argparse

#Start input parser:
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Path to the input mod.gz file.")
parser.add_argument("-o", "--output", help="Output name.")
parser.add_argument("-c", "--conditions", nargs="+", help="Conditions included in the analysis. ie: WT, KO.")
parser.add_argument("-s", "--samples", nargs="+", help="Samples included in the analysis. ie: WT1, WT2, KO1, KO2.")
parser.add_argument("-l", "--labels", nargs="+", help="Group in which the sample is included. ie: 1, 1, 2, 2.")

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
    
def VennDiagram_2groups(n_cond1, n_cond2, n_intersection, output, samples_names, graph_title):
    vd = venn2(subsets = (n_cond1-n_intersection, n_cond2-n_intersection, n_intersection), 
               set_labels = (samples_names[0], samples_names[1]))
    plt.title(graph_title)
    plt.savefig(samples_names[0]+"_"+samples_names[1]+output, dpi=300)
    plt.close()

def ScatterPlots_withR2(twoSamples_df, output, graph_title, log_scale):
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
    
    p.figure.savefig(twoSamples_df.columns[0]+"_"+twoSamples_df.columns[1]+output, dpi=300)
    plt.close(p.figure)

def ScatterPlot_ChangingSites(data, conditions):
    sns.set(rc={'figure.figsize':(7,7)})
    sns.set_theme(style="whitegrid")

    xy = sns.scatterplot(data=data, x=data.columns[6], hue = data["Status"],
                        y=data.columns[7], style = "Group",
                        palette=["grey", "blue", "red"])
    xy.plot([0,1],[0,1], 'black', linewidth=2, linestyle="dashed")
    xy.set(xlabel = conditions[0], ylabel = conditions[1])
    #xy.set(xscale="log", yscale="log")
    xy.set_ylim(0,0.6)
    xy.set_xlim(0,0.6)
    xy.set_title("Median (% Mod) - Replicable sites with Coverage>=50 in ALL samples")
    xy.xaxis.set_major_formatter(mticker.ScalarFormatter())
    xy.yaxis.set_major_formatter(mticker.ScalarFormatter())
    
    xy.figure.savefig(conditions[0]+"_"+conditions[1]+"_ChangingStatus_ReplicablePeaks_AllCov_PeakBased.pdf", dpi=300)
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

def Barplots_ReplicableSites(total_sites, mod_AL1, mod_Both, samples_names):
    
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
    a.figure.savefig(samples_names[0]+"_"+samples_names[1]+"_BarplotsReplicablePeaks_PeakBased.pdf", dpi=300)
    plt.close(a.figure)
    #plt.show()

def mergeLists(list_toMerge, type_merge):
    mergedList = reduce(lambda left, right:
                            pd.merge(left , right,
                            on = ["chr", "pos", "ref_base", "strand", "mod"],
                            how = type_merge),
                            list_toMerge)
    
    return(mergedList)

def extractColumnData(data, initial_column):
    column_list = list()
    n_reps = int((data.shape[1]-5)/4)

    if n_reps>1:
        for i in range(0,n_reps):
            column_list.append(initial_column+(4*i))

    else:
        column_list.append(initial_column)
    
    return(column_list)

## Main function:
def main():
    ##Load input data: 
    input_data = pd.read_table(args.input, skiprows=16)
    labels = args.labels
    samples = args.samples
    conditions = args.conditions
    output = args.output
    
    data = InputCheck_Preprocessing(input_data, labels, samples, conditions)
    
    ##COVERAGE BASED ANALYSIS:
    #Extract coverage and mod_freq information:
    coverage = extractColumnData(data, 5)
    mod_freq = extractColumnData(data, 7)

    peaks_cov_AllSamples = data.loc[(data.iloc[:,coverage]>=50).all(axis=1)]

    #OUTPUT 1: Density plots with data from sites with cov>=50 in ALL samples:
    #ModFreq density plot:
    ModFreq_data = peaks_cov_AllSamples.iloc[:,mod_freq]
    DensityPlots(ModFreq_data, output, '_DensityPlots_ModFreq_CoverageBased.pdf', labels, samples)

    #OUTPUT 2: Scatter plots within replicates including data from sites with cov>=50 in ALL samples:
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
                            "% Mod - Sites with coverage>50 in all samples", False)

        #Scatter plot ModFreq comparing replicates within the same condition (log scale):
        ScatterPlots_withR2(ModFreq_data.iloc[:,modFreq_reps], "_ScatterPlot_LogScale_ModFreq_CoverageBased.pdf", 
                            "% Mod - Sites with coverage>50 in all samples", True)
        
    ##PEAK BASED ANALYSIS:
    #OUTPUT 3: VennDiagrams of replicable peaks (cov>=25 + freq>=0.05) - including sites with cov>25 in all reps too
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
                replicates_data.append(data.iloc[:,columns])
                columns = list(range(0,5))

                #Extract sample names:
                samples_names.append(samples[count])

            else:
                continue

        ##Analyse replicates data:
        reps = mergeLists(replicates_data, "outer")

        #Extract coverage columns:
        coverage_rep = extractColumnData(reps, 5)
        mod_freq_rep = extractColumnData(reps, 7)

        peaks_cov_AllSamples_condition = reps.loc[(reps.iloc[:,coverage_rep]>=50).all(axis=1)]
        peaks_modfreq_AllSamples_condition = peaks_cov_AllSamples_condition.loc[(peaks_cov_AllSamples_condition.iloc[:,mod_freq_rep]>=0.05).all(axis=1)]
        replicable_per_condition.append(peaks_modfreq_AllSamples_condition)

        n_cov = peaks_cov_AllSamples_condition.shape[0]

        ##Assuming only two replicates - NEED TO BE IMPROVED:
        n_cov1 = peaks_cov_AllSamples_condition.loc[(peaks_cov_AllSamples_condition.iloc[:,mod_freq_rep[0]]>=0.05)]
        n_cov2 = peaks_cov_AllSamples_condition.loc[(peaks_cov_AllSamples_condition.iloc[:,mod_freq_rep[1]]>=0.05)]
        inner_ncov = mergeLists(list([n_cov1,n_cov2]), "outer")   
        n_freq = peaks_modfreq_AllSamples_condition.shape[0]

        #Barplots of the total sites, sites modified in at least one rep, sites modified in all reps:
        Barplots_ReplicableSites(n_cov, inner_ncov.shape[0], n_freq, samples_names)

    #OUTPUT 4 and 5: VennDiagrams of replicable peaks (cov>=50 + freq>=0.05) across conditions:
    reps_conditions_inner = mergeLists(replicable_per_condition, "inner")
    reps_conditions_outer_coordinates = mergeLists(replicable_per_condition, "outer").iloc[:,0:5]

    #Venn Diagram of replicate peaks across conditions:
    n_intersection = reps_conditions_inner.shape[0]
    n_cond1 = replicable_per_condition[0].shape[0]
    n_cond2 = replicable_per_condition[1].shape[0]

    VennDiagram_2groups(n_cond1, n_cond2, n_intersection, "_VennDiagram_ReplicableSites_AcrossConditions_PeakBased.pdf", 
                            conditions, "Replicable sites (Cov>=50 and ModFreq>=0.05) across conditions") 

    #Venn Diagram of replicate peaks across conditions with coverage>=25 in ALL samples:
    reps_conditions_outer = pd.merge(reps_conditions_outer_coordinates , data,
                                on = ["chr", "pos", "ref_base", "strand", "mod"],
                                how = "inner")

    outer_allCov = reps_conditions_outer.loc[(reps_conditions_outer.iloc[:,coverage]>=50).all(axis=1)]
    n_intersection = outer_allCov.loc[(outer_allCov.iloc[:,mod_freq]>=0.05).all(axis=1)]
    n_cond1 = pd.merge(outer_allCov,replicable_per_condition[0],
                                on = ["chr", "pos", "ref_base", "strand", "mod"],
                                how = "inner")

    n_cond2 = pd.merge(outer_allCov,replicable_per_condition[1],
                                on = ["chr", "pos", "ref_base", "strand", "mod"],
                                how = "inner")

    VennDiagram_2groups(n_cond1.shape[0], n_cond2.shape[0], n_intersection.shape[0], "_VennDiagram_ReplicableSites_AcrossConditions_CovAllSamples_PeakBased.pdf", 
                            conditions, "Replicable sites (Cov>=50 and ModFreq>=0.05) with Cov>=50 in all samples") 
    
    #OUTPUT 6 and 7: Data from replicable peaks across conditions that have enough coverage in ALL samples: 
    #Classify the results in 2 groups:
    group = list()
    for element in (outer_allCov.iloc[:,mod_freq]>=0.05).all(axis=1):
        if element:
            group.append(2)
        else:
            group.append(1)

    outer_allCov["Group"] = group

    #Calculate medians:
    for j in range(0, len(conditions)):           
        outer_allCov[conditions[j]+"_medianModFreq"] = outer_allCov.filter(like=conditions[j]).filter(like='_ModFreq').median(axis=1)

    #Calculate difference and ratios between conditions:
    medians_cond1 = outer_allCov.iloc[:,outer_allCov.shape[1]-2]
    medians_cond2 = outer_allCov.iloc[:,outer_allCov.shape[1]-1]
    outer_allCov["ModFreq(Cond1-Cond2)"] = medians_cond1 - medians_cond2
    outer_allCov["Ratio(Cond1/Cond2)"] = medians_cond1/medians_cond2

    #Determine directionalities:
    status = list()
    for index,row in outer_allCov.iterrows():
        diff_ModFreq = row["ModFreq(Cond1-Cond2)"]
        ratio_ModFreq = row["Ratio(Cond1/Cond2)"]
        g = row["Group"]

        #Update status:
        if g==1:
            if diff_ModFreq>0:
                status.append("Decreased upon "+conditions[1])
            else:
                status.append("Increased upon "+conditions[1])
        else:
            #Group 2: decrease/increase
            if (diff_ModFreq>=0.1 or ratio_ModFreq>=2):
                status.append("Decreased upon "+conditions[1])
            elif (diff_ModFreq<=(-0.1) or ratio_ModFreq<=0.5):
                status.append("Increased upon "+conditions[1])
            else:
                status.append("No changes")

    outer_allCov["Status"] = status
    outer_allCov.to_csv(output+"_RawData_ReplicablePeaks_AllCov_PeakBased.tsv", sep="\t", index=False)
    outer_allCov_Processed = pd.concat([outer_allCov.iloc[:,0:5], outer_allCov.iloc[:,-6:]], axis=1)
    outer_allCov_Processed.to_csv(output+"_ProcessedData_ReplicablePeaks_AllCov_PeakBased.tsv", sep="\t", index=False)
    
    ## OUTPUT 8: Scatter plot with changing sites (only sites with enough cov in all samples):
    ScatterPlot_ChangingSites(outer_allCov_Processed, conditions)

if __name__ == "__main__":
    main()



