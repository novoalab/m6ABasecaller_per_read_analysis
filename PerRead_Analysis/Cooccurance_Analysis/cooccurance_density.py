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
import scipy as sp

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
    
    #Perform wilcoxon test to compare both distributions: 
    res = mannwhitneyu(random_df.iloc[:,0], cooc_df.iloc[:,0])
    print(res)

    #Plotting:
    plt_title = 'Co-ocurance analysis - Coverage threshold = {} - Min.expected counts threshold = {}'.format(coverage, min_expected)
    ax = sns.kdeplot(data=plotting_data, x='Values', hue='Distribution')

    ax.set_title(plt_title)
    ax.set_xlabel('Standard deviation from expected', fontsize=16)
    ax.set_ylabel('Density', fontsize=16)
    ax.tick_params(labelsize=14)
    ax.figure.savefig('{}_DensityPlot_CoocuranceAnalysis.pdf'.format(args.output), dpi=300)
    plt.close()

def ScatterPlots(sd_data, output):
    distance = abs(sd_data.iloc[:,1]-sd_data.iloc[:,2])
    plotting = pd.DataFrame([distance, sd_data.iloc[:,11]]).transpose()
    plotting.columns = ['Distance', 'SdFromExpected']

    if plotting.shape[0]>=5:
        values = np.vstack([plotting.iloc[:,0].astype('float'), plotting.iloc[:,1].astype('float')])
        kernel = sp.stats.gaussian_kde(values)(values)

        #Scatterplot:
        xy = sns.scatterplot(data=plotting, x=plotting.iloc[:,0],
                            y=plotting.iloc[:,1],
                            c=kernel,
                            cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True))
    else:
        #Scatterplot:
        xy = sns.scatterplot(data=plotting, x='Distance',
                            y='SdFromExpected',
                            cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True))

    #Plot aesthetics:
    xy.set(xlabel = 'log(Genomic distance)', ylabel = 'Sd from expected')
    xy.set(xscale="log")

    #Define y-scale:
    maximum_sd = max(abs(plotting.iloc[:,1])) + 0.1
    xy.set_ylim((maximum_sd*(-1)),maximum_sd)

    #Title:
    xy.set_title("Sd from expected VS genomic distance between sites")

    #Calculate spearman:
    res = sp.stats.spearmanr(plotting.iloc[:,0], plotting.iloc[:,1])
    ax = plt.gca()
    ax.text(.85, .95, '{}={:.3f}'.format(r"$\rho$", res.correlation),
            transform=ax.transAxes)
    ax.figure.savefig('{}_ScatterPlot_CoocuranceAnalysis.pdf'.format(output), dpi=300)
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

#Plot the correlation between sd from expected and distance between sites:
ScatterPlots(input_data, args.output)

'''
#Scatter plot - under testing: 
distance = abs(input_data.iloc[:,1]-input_data.iloc[:,2])
plotting = pd.DataFrame([distance, input_data.iloc[:,11]]).transpose()
plotting.columns = ['Distance', 'SdFromExpected']


if plotting.shape[0]>=5:
    values = np.vstack([plotting.iloc[:,0].astype('float'), plotting.iloc[:,1].astype('float')])
    kernel = sp.stats.gaussian_kde(values)(values)

    #Scatterplot:
    xy = sns.scatterplot(data=plotting, x=plotting.iloc[:,0],
                        y=plotting.iloc[:,1],
                        c=kernel,
                        cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True))
else:
    #Scatterplot:
    xy = sns.scatterplot(data=plotting, x='Distance',
                        y='SdFromExpected',
                        cmap=sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True))

#Plot aesthetics:
xy.set(xlabel = 'log(Genomic distance)', ylabel = 'Sd from expected')
xy.set(xscale="log")

#Define y-scale:
maximum_sd = max(abs(plotting.iloc[:,1])) + 0.1
xy.set_ylim((maximum_sd*(-1)),maximum_sd)

#Title:
xy.set_title("Sd from expected VS genomic distance between sites")

#Calculate spearman:
res = sp.stats.spearmanr(plotting.iloc[:,0], plotting.iloc[:,1])
ax = plt.gca()
ax.text(.85, .95, '{}={:.3f}'.format(r"$\rho$", res.correlation),
        transform=ax.transAxes)
ax.figure.savefig('{}_Scatter_Plot_CoocuranceAnalysis.pdf'.format(args.output), dpi=300)
plt.close()

'''