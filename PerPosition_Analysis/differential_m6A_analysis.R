#### 17/1_2/24 differential m6A methylation analysis 

setwd()
modgz <- "name_of_your_input_table"

## vector with the names of your conditions
conds <- c("Cond1", "Cond2")

#### get the list of sites from the modgz
m6Adf <- as.data.frame(read.table(modgz, header=T, stringsAsFactors = F, sep="\t"))

### create a vector to rename the columns of the modgz based on conditions and reps

conditions <- c("Cond1", "Cond2")
reps <- c("1", "2", "3")
fields <- c("Coverage", "basecall_acc", "ModFreq", "med_mod_prob")
            
# Generate all combinations
combinations <- expand.grid(Condition = conditions, Rep = reps, Fields = fields)
combinations <- combinations[order(combinations$Condition, combinations$Rep, combinations$Fields), ]

# Combine into a single vector of strings
conds_cols <- paste(combinations$Condition, combinations$Rep, combinations$Fields, sep = "_")

## rename the columns            
colnames(m6Adf) <- c("chr", "pos", ".", "strand", "mod", conds_cols)

## create a unique coordinate column
m6Adf$coord <- paste(m6Adf$chr, m6Adf$pos, sep="_")

### when comparing between 2 conditions, keep only sites that have >25 coverage in all reps of both 

m6Adf <- m6Adf[m6Adf$Cond1_1_Coverage>=25 & m6Adf$Cond1_2_Coverage>=25 & m6Adf$Cond1_3_Coverage>=25
               & m6Adf$Cond21_Coverage>=25 & m6Adf$Cond2_2_Coverage>=25 & m6Adf$Cond2_3_Coverage>=25,]

## and >5% frequency in all reps of at least one condition (replicable in that cond)

m6Adf <- m6Adf[(m6Adf$Cond1_1_ModFreq>=0.05 & m6Adf$Cond1_2_ModFreq>=0.05  & m6Adf$Cond1_3_ModFreq>=0.05)
               | (m6Adf$Cond21_ModFreq>=0.05 & m6Adf$Cond2_2_ModFreq>=0.05  & m6Adf$Cond2_3_ModFreq>=0.05),]


###  %diff 

### Cond1 vs Cond2s

DM  <- m6Adf

DM <- DM[, c("coord", "Cond1_1_ModFreq", "Cond1_2_ModFreq", "Cond1_3_ModFreq","Cond2_1_ModFreq", "Cond2_2_ModFreq", "Cond2_3_ModFreq")]


library(matrixTests)

DM_test <- row_t_welch(DM[,c("Cond1_1_ModFreq", "Cond1_2_ModFreq", "Cond1_3_ModFreq")], DM[,c("Cond21_ModFreq", "Cond2_2_ModFreq", "Cond2_3_ModFreq")])

compare_conds <-  cbind(DM[,c("coord","Cond1_1_ModFreq", "Cond1_2_ModFreq", "Cond1_3_ModFreq")], DM[,c("Cond21_ModFreq", "Cond2_2_ModFreq", "Cond2_3_ModFreq")], 
                          DM_test)

compare_conds = mutate(compare_conds, sig=ifelse(compare_conds$pvalue<0.05, "pval<0.05", "Not Sig"))

compare_conds_bed <- compare_conds %>% separate(coord, into=c("chr", "end"), sep="_")
compare_conds_bed$start <-  as.numeric(compare_conds_bed$end)-1
compare_conds_bed <-  compare_conds_bed[,c("chr", "start", "end")]

write.table(compare_conds_bed, "Cond2_vs_Cond1.bed", quote=F, col.names=F, row.names=F,sep="\t")
#### starting from this file, generate the bed and intersect it with gene names 
# bedtools intersect -a Cond2_vs_Cond1.bed -b /no_backup/enovoa/reference_fasta/mouse_genome/gencode.vM24.coord_gene_names.bed -wb > Cond2_vs_Cond1_gene_name.bed

gene_names <- read.table("Cond2_vs_Cond1_gene_name.bed", sep="\t", header=F, stringsAsFactors = F)
gene_names$coord <- paste(gene_names$V1,gene_names$V3, sep="_" )


compare_conds <- merge(unique(gene_names[,c("coord", "V7")]),compare_conds, by="coord", all.y=T)
compare_conds$mean.diff <-  -compare_conds$mean.diff



universe <- as.data.frame(unique(compare_conds$V7))
sig_changed <- as.data.frame(unique(compare_conds[compare_conds$pvalue<=0.05 & abs(compare_conds$mean.diff)>0.2,"V7"]))

write.table(universe, "Cond1_Cond2_m6A_universe.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sig_changed, "Cond1_Cond2_m6A_sig_changed.txt", quote=F, col.names=F, row.names=F, sep="\t")


pdf("Cond2s_vs_Cond1_enhanced_volcano.pdf")
EnhancedVolcano(compare_conds,
                title = 'Cond1 vs Cond2s',
                lab = compare_conds$V7,
                x = 'mean.diff',
                xlim = c(-1, 1),
                ylim = c(0,4),
                pCutoff = 0.05,
                FCcutoff = 0.2,
                labSize = 4.0,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                y = 'pvalue',
                max.overlaps = 40)
dev.off()




universe <- as.data.frame(unique(compare_conds$V7))
sig_changed_up <- as.data.frame(unique(compare_conds[compare_conds$pvalue<=0.05 & compare_conds$mean.diff>=0.1,"V7"]))
sig_changed_down <- as.data.frame(unique(compare_conds[compare_conds$pvalue<=0.05 & compare_conds$mean.diff<=(-0.1),"V7"]))


write.table(universe, "ESC_NPC_m6A_universe.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sig_changed_up, "ESC_NPC_m6A_sig_changed_up.txt", quote=F, col.names=F, row.names=F, sep="\t")
write.table(sig_changed_down, "ESC_NPC_m6A_sig_changed_down.txt", quote=F, col.names=F, row.names=F, sep="\t")

