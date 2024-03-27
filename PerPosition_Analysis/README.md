# Analysis of m6A sites at per position level

ModPhred outputs a file containing all sites that have a coverage of at least 25 reads and a modification frequency of 5% in at least one of the samples included in the analysis. It is called `mod.gz` and it is processed to obtain replicable m6A sites at per position level among other types of information. 

## Table of contents
- [1. Processing ModPhred output](#1-Processing-ModPhred-output)
- [2. Generation of metagene plots from m6ABasecaller results](#2-Generation-of-metagene-plots-from-m6ABasecaller-results)

## 1. Processing ModPhred output

To extract replicable sites from the `mod.gz` file as well as performing metagene and motif enrichment analysis, please use the script `ModPhred_PostProcessing.py`.

### 1.1. Usage: 

```
ModPhred_PostProcessing.py [-h] [-i INPUT] [-o OUTPUT] [-r REFERENCE]
                                  [-c CONDITIONS [CONDITIONS ...]]
                                  [-s SAMPLES [SAMPLES ...]]
                                  [-l LABELS [LABELS ...]] [-gtf GTF_FILE]
                                  [-decay] [-cov COVERAGE] [-modf MOD_FREQ]
                                  [-transcriptome]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input mod.gz file.
  -o OUTPUT, --output OUTPUT
                        Output name.
  -r REFERENCE, --reference REFERENCE
                        Reference file (*.fa).
  -c CONDITIONS [CONDITIONS ...], --conditions CONDITIONS [CONDITIONS ...]
                        Conditions included in the analysis. ie: WT, KO.
  -s SAMPLES [SAMPLES ...], --samples SAMPLES [SAMPLES ...]
                        Samples included in the analysis. 
                        They should go in the same order in which you
                        have them listed in the columns of the mod.gz
                        ie: WT1, WT2, KO1, KO2.
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...]
                        Group in which the sample is included, referring to
                        the order in which you have listed the conditions.
                        ie: 1, 1, 2, 2 for WT,WT,KO,KO
  -gtf GTF_FILE, --gtf_file GTF_FILE
                        Gtf file with genes to annotate the replicable m6A
                        sites (*.gtf).
  -decay, --decay       Include genes that have enough coverage in at least
                        one condition.
  -cov COVERAGE, --coverage COVERAGE
                        Coverage threshold to use in the analysis.
  -modf MOD_FREQ, --mod_freq MOD_FREQ
                        Modification frequency threshold to use in the
                        analysis.
  -transcriptome, --transcriptome
                        Transcriptome analysis - remove sites in the negative
                        strand.

```

Note: the program does not require a matching number of replicates per condition (i.e. you can have 2 reps for WT and 3 for KO)


### 1.2. Examples:

* Example 1: processing the demo data (2replicates, WT and KO conditions, default parameter settings)
```python
python ModPhred_PostProcessing.py -i mod.gz -s condition1_rep1 condition1_rep2 condition2_rep1 condition2_rep2 -c condition1 condition2 -l 1 1 2 2  -o Experiment_Name -gtf your_gene_annotation.gtf
```

*  Example 2: changing the coverage threshold (otherwise, default:50)
```python
python ModPhred_PostProcessing.py -i mod.gz -s condition1_rep1 condition1_rep2 condition2_rep1 condition2_rep2 -c condition1 condition2 -l 1 1 2 2  -o Experiment_Name -gtf your_gene_annotation.gtf -cov 30
```

* Example 3: run the program in *decay* mode: it will consider a site valid if it has enough coverage in all replicates from at least one of the conditions (ie: WT or KO). Default: sites are valid when there is enough coverage across all the samples. 
```python
python ModPhred_PostProcessing.py -i mod.gz -s condition1_rep1 condition1_rep2 condition2_rep1 condition2_rep2 -c condition1 condition2 -l 1 1 2 2  -o Experiment_Name -gtf your_gene_annotation.gtf -decay
```

### 1.3. Expected output:

#### Plots 

- Barplot of replicable sites per condition 

- Scatterplot of ModFreqs between replicates of the same condition

- Scatterplot of average ModFreqs of different conditions

- Venn diagram of sites found replicable in each condition 

- Venn diagram of sites in each replicate of each condition

#### MEME

-output of meme run on the sequence context of the replicable sites

#### Text files

- Raw Data table with raw information from mod.gz filtered according to user's parameters 

- Summary Data with modification frequencies, median frequencies per condition, ratios between conditions and changing status assigned

- Bedgraph with modification frequency difference between condition 1 (reference) and condition 2 per m6A position


## 2. Generation of metagene plots from m6ABasecaller results

First, you'll need to install dependencies: 
```bash
R
install.packages("argparse")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Guitar")
```

```bash
Rscript Metagene_Plots.R -i Sample1.bed Sample2.bed -o Test -gtf Annotation.gtf -l Sample-1 Sample-2
```
The bed file required as input has to contain 6 columns: chr, start, end, "m6A", a numeric value, strand

Expected output:

![image](https://user-images.githubusercontent.com/44866316/196667515-8d0993d9-249a-4f4f-bb5f-2233220a1bf4.png)

## 3. Generation of scatterplots from m6ABasecaller results (for replicability or for comparison) 

The script m6ABC_SCATTERPLOT_M6A_FREQUENCY.R contains the code to produce the scatteplots for checking replicability between replicates (with spearson's coefficient) and for comparing between samples. Please refer to the comments in the script for further detail.

## 4. Generation of density plots from m6ABasecaller results (for comparison) 

The script m6ABC_DENSITY_PLOT_M6A_FREQUENCY.R contains the code to produce the density plots for comparing between samples. Please refer to the comments in the script for further detail.
