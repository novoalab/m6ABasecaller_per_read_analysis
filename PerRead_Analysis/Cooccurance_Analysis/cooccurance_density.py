###################################################################
### Script to plot the density plot from SdFromExpected values. ###
###################################################################

#Import the required libraries: 
import pandas as pd
import numpy as np
import math
import statistics
from scipy.stats import mannwhitneyu
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

#Additional functions:
def DensityPlots(sd_data, output, coverage, min_expected):
    sns.set(rc={'figure.figsize':(11,8)})
    sns.set_theme(style="whitegrid")
    c = sns.color_palette("tab10")[0]

    #Create random distribution to compare:
    sd_random = statistics.stdev(sd_data.iloc[:,11])
    random_distribution = np.random.normal(0, sd_random, len(sd_data.iloc[:,11]))
    random_df = pd.DataFrame({'Values': random_distribution, 'Distribution': 'Random'})

    #Create label from co-occurance data:
    ind_label = '{} (n={})'.format(output, sd_data.shape[0])

    #Prepare plotting data: 
    cooc_df = pd.DataFrame({'Values': sd_data.iloc[:,11], 'Distribution': ind_label})
    plotting_data = pd.concat([random_df, cooc_df]).reset_index()
    print(plotting_data)
    #Perform wilcoxon test to compare both distributions: 
    res = mannwhitneyu(random_df.iloc[:,0], cooc_df.iloc[:,0])
    print(res)

    #Plotting:
    plt_title = 'Co-ocurance analysis - Coverage threshold = {} - Min.expected counts threshold = {}'.format(coverage, min_expected)
    ax = sns.kdeplot(data=plotting_data, x='Values', hue='Distribution')

    #ax = sns.kdeplot(data=plotting_data, x=Values, color = c, label=ind_label)

    #ax.legend(loc='upper right', fontsize=15, frameon=False)
    ax.set_title(plt_title)
    ax.set_xlabel('Standard deviation from expected', fontsize=16)
    ax.set_ylabel('Density', fontsize=16)
    ax.tick_params(labelsize=14)
    ax.figure.savefig('{}_DensityPlot_CoocuranceAnalysis.pdf'.format(args.output), dpi=300)
    plt.close()

#Start input parser:
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Path to the input file containing co-occurance data.")
parser.add_argument("-o", "--output", help="Output name.")
parser.add_argument("-cov", "--coverage", default=200, type=int, help="Minimum coverage to include a transcript in the analysis. Default = 200")
parser.add_argument("-min_exp", "--min_expected", default=2, type=int, help="Minimum expected counts to include a pairwise comparison in the analysis. Default = 2")

#Read arguments from the command line
args = parser.parse_args()

#Open input file:
input_data = pd.read_table(args.input)

#Plot the distribution of the sd_from_expected values:
DensityPlots(input_data, args.output, args.coverage, args.min_expected)
