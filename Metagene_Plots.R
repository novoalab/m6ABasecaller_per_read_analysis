#########################################
### Script to plot m6a metagene plots ###
#########################################

#Load required libraries: 
suppressMessages(library("Guitar"))
library("argparse")

##Argument parser:
#Create parser object
parser <- ArgumentParser()

#Define arguments:
parser$add_argument("-input", "--input_beds", nargs='+', type="character", help="Path to input files (bed).")
parser$add_argument("-gtf", "--gtf_file", nargs=1, type="character", help="Path to annotation file (gtf).")
parser$add_argument("-l", "--labels", nargs='+', type="character", help="Sample names labels to be used in the metagene plot.")
parser$add_argument("-o", "--output", type="character", help="Output filename.")

#Get command line options, if help option encountered - print help and exit:
args <- parser$parse_args()

##Main functionalities:
#GTF file:
print('Importing GTF file.')
gtfFile <- args$gtf_file
txdb2 <- makeTxDbFromGFF(file=gtfFile,dataSource="ensemblgenomes")

#BED files:
print('Start processing bed files.')
if (length(args$input)!=length(args$labels)){
  stop("Number of samples must be equal to number of labels. Try again. Bye!")
} 
stBedFiles <- args$input_beds

# Guitar Plot
pdf(paste(args$output, '_metagene_plot.pdf', sep=""))
GuitarPlot(txTxdb = txdb2,
           stBedFiles = stBedFiles,
           headOrtail = TRUE,
           enableCI = FALSE,
           mapFilterTranscript = TRUE,
           pltTxType = c("mrna"),
           stGroupName = args$labels)

dev.off()