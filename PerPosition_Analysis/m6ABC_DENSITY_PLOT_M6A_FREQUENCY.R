#### DENSITY PLOT FOR M6A FREQUENCY DISTRIBUTION
library(tidyr)
library(tidyverse) 
library(dplyr)
library(ggplot2)
library(ggrepel) 
library(reshape2)
library(ggpubr)

### set working directory. i.e.:
setwd("/no_backup_isis/enovoa/users/scruciani/ModPhred/XPORE/compare_ivt_reps/")

#### in this folder, where you had the mod.gz output file, you have generated a mod_gz.txt file that
### contains the table and is uncompressed with: zcat mod.gz | tail -n +17 > mod_gz.txt

### import the raw output of mmodphred in mod_gz.txt

mod_gz <- as.data.frame(read.table(file="/no_backup_isis/enovoa/users/scruciani/ModPhred/XPORE/compare_ivt_reps/mod_gz.txt", header=T, stringsAsFactors = F, sep="\t"))

## rename columns for easier manipulation. 
colnames(mod_gz) <- c("chr", "pos", "ref_base", "strand", "mod", 
                      "IVT1_depth", "IVT1_basecall_accuracy", "IVT1_mod_freq", "IVT1_mod_prob","IVT2_depth", "IVT2_basecall_accuracy", "IVT2_mod_freq", "IVT2_mod_prob", 
                      "KO2_depth", "KO2_basecall_accuracy", "KO2_mod_freq", "KO2_mod_prob","KO3_depth", "KO3_basecall_accuracy", "KO3_mod_freq", "KO3_mod_prob", 
                      "WT2_depth", "WT2_basecall_accuracy", "WT2_mod_freq", "WT2_mod_prob","WT3_depth", "WT3_basecall_accuracy", "WT3_mod_freq", "WT3_mod_prob")

#### filter by >=5% frequency in all WT/CTR samples you have and by >=25 reads
#### coverage in all the samples you want to compare

mod_gz <-  mod_gz[mod_gz$WT2_depth>=25 & mod_gz$WT3_depth>=25 & mod_gz$IVT1_depth>=25 & mod_gz$IVT2_depth>=25 & mod_gz$KO2_depth>=25 & mod_gz$KO3_depth>=25,]

mod_gz <-  mod_gz[ mod_gz$WT2_mod_freq >=0.05 & mod_gz$WT3_mod_freq >=0.05,]


### after filtering, only keep the columns with frequency information 

boxplot_table <- mod_gz[,c(grep("mod_freq", colnames(mod_gz)))]  

### rename columns based on sample and replicate - separated by a dash 

colnames(boxplot_table) <- c("IVT-1", "IVT-2", "KO-1", "KO-2", "WT-1", "WT-2")

### convert the dataframe in long format for boxplots and density plots

boxplot_table <-  melt(boxplot_table)

### generate the "sample" column 

boxplot_table <-  boxplot_table %>% separate(variable, c("sample", "replicate"), sep="-", remove=F)

### generate the final table for plotting
table_boxplot <- boxplot_table

## convert the sample column from character to factor
table_boxplot$sample <-  factor(table_boxplot$sample, levels = c( "WT", "KO", "IVT"))

### VERY IMPORTANT FOR LOG SCALED PLOT: add a small pseudocount (0.1% frequency) to all frequency values
### if you have frequencies that equal to 0, so you avoid errors when plotting in log scaled axis
table_boxplot$value <-  table_boxplot$value+0.001

# Load the data.table library
library(data.table)
# Convert the data.frame to a data.table
df <- setDT(table_boxplot)
# Calculate the median value for each variable level
median_table <- df[, .(median_value = median(value, na.rm = TRUE)), by = variable]

#### plot the density with median frequency per sample 
pdf(file="m6A_frequency_HEK293_log_pseudocount001.pdf", height=10,width=12,onefile=FALSE)
print(ggplot(table_boxplot, aes(x=value, fill=variable, label=sample, adjust=2)) +
        geom_density(alpha=0.3, adjust=2)+
        theme_bw()+
        scale_x_continuous(trans="log10", limits=c(0.001, 1))+
        ggtitle("m6A frequency")+
        xlab("m6A frequency")+
        ylab("Density") +
        geom_vline(data = median_table, aes(xintercept = median_value,color = variable), size=1, linetype="dashed")+
        theme(axis.text=element_text(size=15),strip.text = element_text(size=13),
              axis.title=element_text(size=17,face="bold"),legend.title = element_text(size = 20),
              plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
              legend.text = element_text(color = "black", size=15)))
dev.off()
