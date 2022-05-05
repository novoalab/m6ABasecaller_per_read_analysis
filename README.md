# ModPhred_PostProcessing

# usage:

```
python ./ModPhred_PostProcessing/ModPhred_PostProcessing.py -i mod.gz -s condition1_rep1 condition1_rep2 condition2_rep1 condition2_rep2 -c condition1 condition2 -l 1 1 2 2  -o Experiment_Name -bed gene_coordinates_with_gene_names.bed
```
-i: this is the raw output file of ModPhred 

-s: names of the samples as you have input them to ModPhred, which correspond to the order of the columns in mod.gz. requires at least 2

-c: experimental conditions (i.e. WT KO). requires at least 1.

-l: ordinal number matching the sample names to the conditions (i.e. -s WT1 WT2 KO1 KO2 KO3 -c WT KO -l 1 1 2 2 2). 

-o: output folder name

-bed: if you want to know in which genes the m6A sites are found, provide a bed file with the gene start and end coordinates and with a gene identifier at your choice (gene name, gene ID... etc)

Note: the program does not require a matching number of replicates per condition (i.e. you can have 2 reps for WT and 3 for KO)
