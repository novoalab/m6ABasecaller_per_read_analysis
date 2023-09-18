##### SCATTERPLOT M6A FREQUENCY
library(tidyr)
library(tidyverse) 
library(dplyr)
library(ggplot2)
library(ggrepel) 
library(reshape2)
library(ggpubr)

#### 1. if you want to plot REPLICABILITY then filter by >=5% frequency only if it is WT/CTR
#### and by >=25 reads coverage in the two replicates you want to compare

#### 2. if you want to plot CHANGES then filter by >=5% frequency in WT/CTR
#### and by >=25 reads coverage in the samples you want to compare


#### 1. if you want to plot REPLICABILITY
### set working directory. i.e.:
setwd("/no_backup_isis/enovoa/users/scruciani/ModPhred/XPORE/compare_ivt_reps/")

#### in this folder, where you had the mod.gz output file, you have generated a mod_gz.txt file that
### contains the table and is uncompressed with: zcat mod.gz | tail -n +17 > mod_gz.txt

### import the raw output of mmodphred in mod_gz.txt

mod_gz <- as.data.frame(read.table(file="/no_backup_isis/enovoa/users/scruciani/m6A_basecaller_benchmark/figure2/mouse_ELIGOS_WT/mod_WT", header=T, stringsAsFactors = F, sep="\t"))

## rename columns for easier manipulation. 
colnames(mod_gz) <- c("chr", "pos", "ref_base", "strand", "mod", "WT1_depth", "WT1_basecall_accuracy", "WT1_mod_freq", "WT1_mod_prob","WT2_depth", "WT2_basecall_accuracy", "WT2_mod_freq", "WT2_mod_prob")


mod_gz <-  mod_gz[mod_gz$WT1_depth>=25 & mod_gz$WT2_depth>=25 
                  & mod_gz$WT1_mod_freq >=0.05 & mod_gz$WT2_mod_freq >=0.05,]

###create the column "coordinate"

mod_gz$coord <- paste(mod_gz$chr, mod_gz$pos, sep="_")

#### keep only coordinate and modfreq of the two replicates
my_data <- mod_gz[c("coord",8,12)]

### if you are showing replicability of samples that are not WT/CTR, you might have
### some frequency values that are 0. in this case you should add a 0.01 pseudocount to all
### frequency values.

#### scatterplot to show replicability (squared) - log10 scaled
#### (dashed line shows 1:1 ratio x:y, "R" shows spearman's rho, 
#### dashed lines indicate the 5% threshold, a gradient of light blue to dark blue shows increasing 
#### data point density)

pdf(file="replicability_eligos_WT_threshold25_5freq-xaxis001.pdf", height=5,width=5,onefile=FALSE)
dcols<-densCols(my_data[,2],my_data[,3], colramp=colorRampPalette(blues9[-(1:3)]))
print(ggplot(my_data) +
        geom_point(aes(x=my_data[,2], y=my_data[,3], col=dcols), size=1 ) +
        scale_color_identity() +
        geom_abline(slope=1, intercept=0,linetype="dashed")+
        geom_vline(xintercept=0.05,linetype="dashed")+
        geom_hline(yintercept=0.05,linetype="dashed")+
        ggtitle("replicability ")+ 
        xlab("rep1")+
        ylab("rep2") +
        stat_cor(aes(x=my_data[,2], y=my_data[,3]),method = "spearman")+ #, label.x = 0.01, label.y = 1.2)+ # position of the R corr on the plot - it'll depend on your data and its limits, can make it a variable
        theme_bw()+
        scale_x_continuous(trans='log10', limits = c(0.01,1)) +
        scale_y_continuous(trans='log10', limits = c(0.01,1))+
        theme(axis.text.x = element_text(face="bold", color="black",size=11),
              axis.text.y = element_text(face="bold", color="black", size=11),
              plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
              axis.title.x = element_text(color="black", size=15, face="bold"),
              axis.title.y = element_text(color="black", size=15, face="bold"),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black", size=0.5), 	  
              legend.title = element_text(color = "black", size = 15,face="bold"),
              legend.text = element_text(color = "black", size=15),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              aspect.ratio = 1)
)
dev.off()

#### 2. if you want to plot CHANGES 

### set working directory. i.e.:
setwd("/no_backup/enovoa/users/scruciani/ModPhred/XPORE/compare_ivt_all_pooled2/")

#### in this folder, where you had the mod.gz output file, you have generated a mod_gz.txt file that
### contains the table and is uncompressed with: zcat mod.gz | tail -n +17 > mod_gz.txt

### import the raw output of mmodphred in mod_gz.txt

mod_gz <- as.data.frame(read.table(file="mod_gz.txt", header=T, stringsAsFactors = F, sep="\t"))

## rename columns for easier manipulation. 
colnames(mod_gz) <- c("chr", "pos", "ref_base", "strand", "mod", "IVT_depth", "IVT_basecall_accuracy", "IVT_mod_freq", "IVT_mod_prob",
                      "KO_depth", "KO_basecall_accuracy", "KO_mod_freq", "KO_mod_prob", 
                      "WT_depth", "WT_basecall_accuracy", "WT_mod_freq", "WT_mod_prob")

#### filter by >=5% frequency in WT/CTR
#### and by >=25 reads coverage in both samples you want to compare

mod_gz <- mod_gz[mod_gz$WT_mod_freq>=0.05 & mod_gz$WT_depth>=25 & mod_gz$IVT_depth>=25,]

###create the column "coordinate"

mod_gz$coord <-  paste(mod_gz$chr, mod_gz$pos, sep="_")

#### keep only coordinate and modfreq of the two samples to compare

my_data <- mod_gz[c("coord",8,12)]

##### add pseudocount (1%) for plotting data in log scaled axis
my_data$WT_mod_freq <-  my_data$WT_mod_freq+0.01
my_data$IVT_mod_freq <-  my_data$IVT_mod_freq+0.01

#### scatterplot to show CHANGES (squared) - log10 scaled
#### (dashed line shows 1:1 ratio x:y, dashed lines indicate the 5% threshold, 
#### a gradient of light blue to dark blue shows increasing data point density)

pdf(file="XPORE_WT_IVT_log-pseudocount01.pdf", height=5,width=5,onefile=FALSE)
dcols<-densCols(my_data[,2],my_data[,3], colramp=colorRampPalette(blues9[-(1:3)]))
print(ggplot(my_data) +
        geom_point(aes(x=my_data[,2], y=my_data[,3], col=dcols), size=1 ) +
        scale_color_identity() +
        geom_abline(slope=1, intercept=0,linetype="dashed")+
        geom_vline(xintercept=0.06,linetype="dashed")+ ### this is the 5% threshold corrected because of the pseudocount
        geom_hline(yintercept=0.06,linetype="dashed")+ ### this is the 5% threshold corrected because of the pseudocount
        ggtitle("m6A frequency")+ 
        xlab("WT")+
        ylab("IVT") +
        theme_bw()+
        scale_x_continuous(trans="log10", limits = c(0.01,1)) +
        scale_y_continuous(trans="log10",limits = c(0.01,1))+
        theme(axis.text.x = element_text(face="bold", color="black",size=11),
              axis.text.y = element_text(face="bold", color="black", size=11),
              plot.title = element_text(color="black", size=15, face="bold.italic",hjust = 0.5),
              axis.title.x = element_text(color="black", size=15, face="bold"),
              axis.title.y = element_text(color="black", size=15, face="bold"),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black", size=0.5), 	  
              legend.title = element_text(color = "black", size = 15,face="bold"),
              legend.text = element_text(color = "black", size=15),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              aspect.ratio = 1)
)
dev.off() 
