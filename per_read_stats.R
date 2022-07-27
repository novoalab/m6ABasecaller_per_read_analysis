#Per read stats table: 

#Required libraries:
library(dplyr)
library(argparse)
library(stringr)
library(ggplot2)

##Parse input:
#Create parser object
parser <- ArgumentParser()

parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("-i", "--input", nargs="+", type="character",
                    help="Path to the per_read table stats to be analysed globally.")
parser$add_argument("-o", "--output", nargs="+", type="character",
                    help="Path to the per_read table stats to be analysed globally.")
parser$add_argument("--per_gene", default=30, type="double",
                    help="Minimum counts per gene to be analyzed in future analysis [default %(default)s]")
parser$add_argument("--per_isoform", default=10, type="double",
                    help="Minimum counts per isoform to be analyzed in future analysis [default %(default)s]")

#Get command line options:
args <- parser$parse_args()

##Load input: 
initial <- TRUE
files <- args$input
#files <- c("per_read_WT_2hpf_rep1.tsv", "per_read_WT_2hpf_rep2.tsv")

#Parse one file at a time: 
for (file in files){
  
  #Extract sample name - based on file name:
  sample_name <- str_match(file, "per_read_\\s*(.*?)\\s*.tsv")[2]
  
  #Load data from individual sample:
  data <- read.table(file, header=T)
  
  ##Calculate stats: 
  #Number of m6A sites:
  mean_num_m6A_sites <- mean(data$num_m6a_sites)
  median_num_m6A_sites <- median(data$num_m6a_sites)
  sample_m6A_sites <- cbind(sample_name, data$num_m6a_sites)

  #Median pA tail length: 
  median_nanpolish <- median(data$nanopolish, na.rm = T)
  median_tailfindr <- median(data$tailfindr, na.rm = T)
  sample_nanopolish <- cbind(sample_name, data$nanopolish)
  sample_tailfindr <- cbind(sample_name, data$tailfindr)
  
  #Type of gene/isoform assignment: 
  test <- group_by(data, assignment_type) %>% summarise(Count = length(assignment_type))
  test$total <- sum(test$Count)
  test$percentage <- 100*(test$Count/test$total)
  
  #Number of reads per gene:
  data$gene_id <- replace(data$gene_id, data$gene_id==".", NA)
  genes <- data %>% count(gene_id)
  median_reads_per_gene <- median(genes$n)
  total_num_genes <- length(genes$gene_id)
  cov_num_genes <- nrow(subset(genes, n>=args$per_gene))
  per_cov_num_genes <- 100*(cov_num_genes/total_num_genes)
  sample_reads_per_gene <- cbind(sample_name, genes$n)

  #Number of reads per isoform:
  data$isoform_id <- replace(data$isoform_id, data$isoform_id==".", NA)
  isoforms <- data %>% count(isoform_id)
  median_reads_per_isoform <- median(isoforms$n)
  total_num_isoforms <- length(isoforms$isoform_id)
  cov_num_isoforms <- nrow(subset(isoforms, n>=args$per_isoform))
  per_cov_num_isoforms <- 100*(cov_num_isoforms/total_num_isoforms)
  sample_reads_per_isoform <- cbind(sample_name, isoforms$n)


  ##Generate final table:
  sample_stats <- c(sample_name, round(mean_num_m6A_sites,2), median_num_m6A_sites, median_nanpolish, median_tailfindr, 
                    total_num_genes, median_reads_per_gene, round(per_cov_num_genes,2), 
                    total_num_isoforms, median_reads_per_isoform, round(per_cov_num_isoforms,2))
  
  if (initial){
    #Stats:
    final_stats <- as.data.frame(rbind(sample_stats))
    colnames(final_stats) <- c("Sample", "Mean m6A sites per read", "Median m6A sites per read",
                               "Median pA length by Nanopolish (bp)", "Median pA length by Tailfindr (bp)",
                               "Total number of genes", "Median number of reads per gene", 
                               paste("Number of genes with counts >= ", args$per_gene, " (%)",sep=""),
                               "Total number of isoforms", "Median number of reads per isoform",
                               paste("Number of isoforms with counts >= ", args$per_isoform, " (%)",sep="")) 
    initial <- FALSE
    
    #m6A sites per read:
    data_m6A_sites <- as.data.frame(sample_m6A_sites)

    #Nanopolish data:
    data_nanopolish <- as.data.frame(sample_nanopolish)
    
    #Tailfindr data:
    data_tailfindr <- as.data.frame(sample_tailfindr)

    #Reads per gene:
    data_reads_per_gene <- as.data.frame(sample_reads_per_gene)

    #Reads per isoform:
    data_reads_per_isoform <- as.data.frame(sample_reads_per_isoform)
    
  } else {
    #Stats:
    final_stats <- rbind(final_stats, sample_stats)
    
    #m6A sites per read:
    data_m6A_sites <- rbind(data_m6A_sites, sample_m6A_sites)

    #Nanopolish data:
    data_nanopolish <- rbind(data_nanopolish, sample_nanopolish)
    
    #Tailfindr data:
    data_tailfindr <- rbind(data_tailfindr, sample_tailfindr)

    #Reads per gene:
    data_reads_per_gene <- rbind(data_reads_per_gene, sample_reads_per_gene)

    #Reads per isoform:
    data_reads_per_isoform <-rbind(data_reads_per_isoform, sample_reads_per_isoform)

  }
}

##Export final stats table: 
write.table(final_stats, paste(args$output,"_per_read_stats.tsv",sep=""), sep="\t", row.names=F, quote = F)

##Plotting violin plots for the different features:
data_m6A_sites$V2 <- factor(data_m6A_sites$V2, levels=c("0","1","2","3","4"))
data_m6A_sites$V2[data_m6A_sites$V2=="NA"] <- "5+"
pdf(paste(args$output,'_m6A_sites_per_read.pdf',sep=""),height=5,width=14,onefile=FALSE)
plt_violin_base <- ggplot(data_m6A_sites, aes(x=V2)) +
  geom_histogram(binwidth=1, stat="count") +
  xlab('m6A sites per read') +
  ylab('Counts') +
  theme_classic(base_size=13)+
  facet_grid(~ sample_name)

plt_violin_base
dev.off()

data_nanopolish$V2 <- as.numeric(data_nanopolish$V2)
pdf(paste(args$output,'_pA_Nanopolish.pdf',sep=""),height=5,width=14,onefile=FALSE)
plt_violin_base <- ggplot(data_nanopolish, aes(x=sample_name, y=V2)) +
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill="white")+
  ylab('PolyA tail length (bp)') + xlab("")+
  theme_classic(base_size=13)

plt_violin_base
dev.off()

data_tailfindr$V2 <- as.numeric(data_tailfindr$V2)
pdf(paste(args$output,'_pA_Tailfindr.pdf',sep=""),height=5,width=14,onefile=FALSE)
plt_violin_base <- ggplot(data_tailfindr, aes(x=sample_name, y=V2)) +
  geom_violin(trim=TRUE)+
  geom_boxplot(width=0.1, fill="white")+
  ylab('PolyA tail length (bp)') + xlab("")+
  theme_classic(base_size=13)

plt_violin_base
dev.off()

data_reads_per_gene$V2 <- as.numeric(data_reads_per_gene$V2)
pdf(paste(args$output,'_reads_per_gene.pdf',sep=""),height=5,width=14,onefile=FALSE)
plt_violin_base <- ggplot(data_reads_per_gene, aes(x=sample_name, y=V2)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  ylab('Reads per gene') + xlab("")+
  theme_classic(base_size=13)

plt_violin_base
dev.off()

#data_reads_per_isoform$V2 <- as.numeric(data_reads_per_isoform$V2)
#pdf(paste(args$output,'_reads_per_isoform.pdf',sep=""),height=5,width=14,onefile=FALSE)
#plt_violin_base <- ggplot(data_reads_per_isoform, aes(x=sample_name, y=V2)) +
#  geom_violin(trim=FALSE)+
#  geom_boxplot(width=0.1, fill="white")+
#  ylab('Reads per isoform') + xlab("")+
#  theme_classic(base_size=13)

#plt_violin_base
#dev.off()
