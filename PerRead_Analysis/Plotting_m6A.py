######################################################
### Script to plot the m6A data at per read level. ###
######################################################

#Import the required libraries:
import pandas as pd
import seaborn as sns
import itertools
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse

#Additional functions:
def extract_shorter_UTR(transcript):
    coordinates = transcript.iloc[:,[5,6,7]].drop_duplicates()
    shorter_start = sorted(coordinates['upd_transcript_start'])[coordinates.shape[0]-1]
    shorter_end = sorted(coordinates['upd_transcript_end'])[0]
    
    return shorter_start, shorter_end

def extract_m6A_sites(transcript, only_common): 
    #Extract unique m6A sites within the transcript: 
    m6A_sites = list()
    for i in range(0, transcript.shape[0]):
        m6A_sites.append(transcript.iloc[i,1])

    unique_m6A_sites = set([str(i).split(',')[0] for i in set(m6A_sites)])
    unique_m6A_sites = [item for item in unique_m6A_sites if item != 'nan']

    if only_common:
        
        #Extract the shorter UTR from all the isoforms included:
        shorter_start, shorter_end = extract_shorter_UTR(transcript)
        
        #Remove those sites that fall outside the shorter UTRs from all isoforms included in the plot:
        ordered_sites = list()

        for i in unique_m6A_sites:
            i_int = int(float(i))

            if i_int >= shorter_start and i_int <= shorter_end:
                ordered_sites.append(i_int)

    else:
        ordered_sites = [int(float(i)) for i in unique_m6A_sites]
    
    #Order the list for plotting purposes:
    ordered_sites.sort()
    
    return(ordered_sites)

def create_m6A_matrix(transcript, ordered_sites):
    #Create presence/absence data for heatmap for every transcript:
    initial_df = True

    for ind_transcript in transcript['upd_transcript_id'].unique():
        #Subset per isoform:
        subset_isoform = transcript.loc[transcript['upd_transcript_id'] == ind_transcript]

        #Heatmap creation:
        initial_site = True
        for j in ordered_sites:
            column = list()
            j = str(j)
            
            for k in range(0, subset_isoform.shape[0]):
                if pd.isnull(subset_isoform.iloc[k,1]):
                    column.append(0)
                elif j in str(subset_isoform.iloc[k,1]).split(','): 
                    column.append(1)
                else: 
                    column.append(0)

            if initial_site:
                plotting = pd.DataFrame(column, columns=[j])
                initial_site = False

            else:
                plotting[j] = column

        #Add isoform information:
        plotting['Isoform'] = subset_isoform.iloc[:,5].tolist()
        plotting['Samples'] = subset_isoform.iloc[:,4].tolist()
        
        #Order dataframe by rows:
        if initial_df:
            plot_heatmap = plotting.sort_values(by=[str(i) for i in ordered_sites], ascending=False)
            initial_df = False

        else:
            plot_heatmap = pd.concat([plot_heatmap, plotting.sort_values(by=[str(i) for i in ordered_sites], ascending=False)])
        
    return(plot_heatmap)

#Start input parser:
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Path to the input file containing reannotated per read label data (output of annotate_UTRs.py).")
parser.add_argument("-o", "--output", help="Output name.")
parser.add_argument("-common", "--common_m6A", help="Only plot m6A shared by all the isoforms plotted together.", action="store_true")
parser.add_argument("-g_list", "--gene_list", nargs="+", help="List of genes to plot.", default="")
parser.add_argument("-i_list", "--isoform_list", nargs="+", help="List of isoforms to plot.", default="")
parser.add_argument("-per_isoform", "--per_isoform", help="Generate plots at per isoform level.", action="store_true")
parser.add_argument("-s_list", "--sites_list", nargs="+", help="List of sites to plot.", default="")

#parser.add_argument("-min_exp", "--min_expected", default=2, type=int, help="Minimum expected counts to include a pairwise comparison in the analysis. Default = 2")
#ENSG00000213741
#Read arguments from the command line
args = parser.parse_args()

#Read input data: 
data = pd.read_table(args.input, usecols=[0,9,12,13,16,17,18,19])

#Analyse the modification pattern from the same feature:
#Extract feature list:
if args.per_isoform:
    if len(args.isoform_list)==0:
        features = data['isoform_id'].unique()
    else:
        features = args.isoform_list
else:
    if len(args.gene_list)==0:
        features = data['gene_id'].unique()
    else:
        features = args.gene_list

#Start analysis:
for feature_id in features:
    #Subset data:
    if args.per_isoform:
        transcript = data.loc[data['isoform_id'] == feature_id]
    else:
        transcript = data.loc[data['gene_id'] == feature_id]
    
    #Only proceed if the feature has more than one UTR combination:
    if len(transcript['upd_transcript_id'].unique())>1:
        
        #Extract m6A sites to analyse:
        if len(args.sites_list) == 0:
            ordered_sites = extract_m6A_sites(transcript, args.common_m6A)
        else:
            ordered_sites = args.sites_list
        
        #If the gene is unmodified, proceed with the next one:
        if len(ordered_sites) == 0:
            continue
            
        #Create presence/absence matrix for heatmap:
        plotting_data = create_m6A_matrix(transcript, ordered_sites)
        
        #Add color information for different annotations: 
        ##By isoform:
        isoforms = plotting_data.pop("Isoform")
        pal = sns.cubehelix_palette(isoforms.unique().size, light=.9, dark=.1, reverse=False, start=1, rot=-2)
        lut = dict(zip(sorted(isoforms.unique(), key=lambda x: x.split("_")[-1]), pal))

        ##By sample:
        samples = plotting_data.pop("Samples")
        lut2 = dict(zip(samples.unique(), sns.color_palette("dark")))

        ##Provide the two color scales:
        row_colors = [samples.map(lut2),isoforms.map(lut)]
        
        #Plotting the heatmap:
        sns.set(rc={"figure.dpi":300, 'savefig.dpi':300})
        g = sns.clustermap(plotting_data, cmap="Blues", figsize=(16, 10), row_colors=row_colors,
                           row_cluster=None, col_cluster=None, yticklabels=False,
                           cbar_kws={'label': 'm6A', 'orientation': 'horizontal', "ticks":[0,1]})

        #Add legends:
        l1 = g.fig.legend(loc='lower left',bbox_to_anchor=(0.01, 0.47), frameon=True, title="Isoforms", 
                                  handles=[mpatches.Patch(color=c, label=l) for l, c in lut.items()], fontsize='14')
        l1.set_title(title='Isoforms',prop={'size':16})

        l2 = g.fig.legend(loc='lower left',bbox_to_anchor=(0.01, 0.7), frameon=True, title="Samples", 
                                  handles=[mpatches.Patch(color=c, label=l) for l, c in lut2.items()], fontsize='14')
        l2.set_title(title='Samples',prop={'size':16})

        #Adjust legend horizontal and top left corner:
        g.cax.set_position([.01, .05, .15, .03])#[.8, .96, .2, .03])
        # set title and labels for axes
        g.fig.suptitle("Modification pattern by feature: %s "%(feature_id), fontsize='16', y=0.85, x=0.62)
        g.ax_heatmap.set_xlabel("Modified positions", fontsize=14)
        g.ax_heatmap.set_ylabel("Reads",  fontsize=14)
        g.ax_heatmap.tick_params(right=False)

        #Save plot into figure:
        outname = args.output + "ModifiedPositions_" + feature_id + ".pdf"
        g.savefig(outname, bbox_inches='tight') 
        plt.close(g.fig)
        
    else:
        continue
                