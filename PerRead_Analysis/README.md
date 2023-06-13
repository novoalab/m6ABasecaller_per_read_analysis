# Analysis of m6A sites at per read level

The m6A basecaller allows the analysis of m6A data at per single molecule level. Before proceeding with this analysis, please first perform the analysis at per position level here (MISSING LINK). 

## 1. Generating per-read level tables with m6A, polyA tail and isoform information

### 1.1. Parse results at per read level from the m6A basecaller:
- Usage:

```bash
python ./PerRead_Table/get_mods.py -i input.bam -o output.tsv -q 15 -c all -s yes

awk '{sub(/True/, "-", $4)}1' output.tsv | awk '{sub(/False/, "+", $4)}1' - > PerRead.tsv
```

- Example:

```bash
python ./PerRead_Table/get_mods.py -i ./example_input/example_input.bam -o ./example_output/parse_m6A.tsv -q 15 -c all -s yes

awk '{sub(/True/, "-", $4)}1' ./example_output/parse_m6A.tsv | awk '{sub(/False/, "+", $4)}1' - > ./example_output/parse_strand_m6A.tsv
```

- Expected output: 
```bash
readname N_mod ref_name is_reverse start end readlen alnlen is_secondary mod_list
c297b2e7-baf1-4678-9a07-9eeee397e3ff 0 chr1 + 629651 630666 904 881 False NA
0db3c821-2930-4107-a4e5-4e3a3ad09d68 1 chr1 + 632766 633436 638 612 False 633101
de5af892-650c-4052-9d04-36cb9b49af76 0 chr1 + 632766 633443 622 593 False NA
9f40bed2-2428-4c22-8254-ff702e9d8e92 0 chr1 + 632766 633442 671 610 False NA
5ed78e6a-2f26-4346-82f6-98e3f95ad5b5 0 chr1 + 632766 633440 675 602 False NA
```

### 1.2. Filter out non-replicable sites from the m6A data at per read level:
- Usage:
```bash
python ./PerRead_Analysis/PerRead_Table/filter_mod_sites.py -pr PerRead.tsv -rs PerPosition_SummaryData_ReplicableSites.tsv -o Filtered

```

- Example:

```bash
python ./PerRead_Table/get_mods.py -i ./example_input/example_input.bam -o ./example_output/parse_m6A.tsv -q 15 -c all -s yes

awk '{sub(/True/, "-", $4)}1' ./example_output/parse_m6A.tsv | awk '{sub(/False/, "+", $4)}1' - > ./example_output/parse_strand_m6A.tsv
```

- Expected output:
```bash
readname N_mod ref_name strand start end readlen alnlen is_secondary mod_list
6c53fb4e-3deb-40dc-a23e-0321d969d2b5 1 chr1 - 4490932 4493594 2179 2175 False 4492362
d12b1362-fd55-4374-b46b-8ad13a3741f3 0 chr1 - 4490943 4492005 1093 1065 False NA
8613f34b-01cb-4c46-b419-f745147b28f1 0 chr1 - 4491371 4493583 1732 1705 False NA
6b124dcc-a013-4a62-b620-f40a47366480 2 chr1 - 4491381 4493181 1217 1205 False 4492245,4492362
```

### 1.3. Generate final table with m6A sites, isoform and polyA tail data for every read in the sample:

- Usage:
```bash
bash ./PerRead_Analysis/PerRead_Table/Generate_perread_tables.sh Per_read_Filtered_ReplicableSites.tsv pA_tail.txt isoquant.tsv sample_name output_name

```

- Expected output:
```bash
TO FILL IN
```

## 2. Re-annotating 5' and 3' UTR from reads belonging to the same isoform

First, if we are analysing multiple samples, we should first merge the per read tables from each one of them. 

- Usage:
```bash
python ./PerRead_Analysis/PerIsoform_Analysis/annotate_UTRs.py -i per_read_table.tsv -o Reannotated

```

- Command line example:
```bash
python ./PerRead_Analysis/PerIsoform_Analysis/annotate_UTRs.py -i XXX -o XXX

```

- Expected output:
```bash
TO FILL IN
```

## 3. Modification frequency analysis



## 4. Co-occurance analysis




