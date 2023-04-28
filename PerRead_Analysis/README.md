# Analysis of m6A sites at per read level

The m6A basecaller allows the analysis of m6A data at per single molecule level. Before proceeding with this analysis, please first perform the analysis at per position level here (MISSING LINK). 

## 1. Generating per-read level tables with m6A, polyA tail and isoform information

### 1.1. Parse results at per read level from the m6A basecaller:
To proceed with this step, please download this [GitHub](https://github.com/biocorecrg/nanomod_map) repository. 
- Usage:

```bash
python get_mods.py -i input.bam -o output.tsv -q 15 -c all -s yes

awk '{sub(/True/, "-", $4)}1' output.tsv | awk '{sub(/False/, "+", $4)}1' - > PerRead.tsv
```

- Expected output: 
```bash
readname N_mod ref_name strand start end readlen alnlen is_secondary mod_list
6c53fb4e-3deb-40dc-a23e-0321d969d2b5 1 chr1 - 4490932 4493594 2179 2175 False 4492362
d12b1362-fd55-4374-b46b-8ad13a3741f3 0 chr1 - 4490943 4492005 1093 1065 False NA
8613f34b-01cb-4c46-b419-f745147b28f1 1 chr1 - 4491371 4493583 1732 1705 False 4492352
6b124dcc-a013-4a62-b620-f40a47366480 4 chr1 - 4491381 4493181 1217 1205 False 4491528,4492245,4492352,4492362
```

### 1.2. Filter out non-replicable sites from the m6A data at per read level:
- Usage:
```bash
python ./PerRead_Analysis/PerRead_Table/filter_mod_sites.py -pr PerRead.tsv -rs PerPosition_SummaryData_ReplicableSites.tsv -o Filtered

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
bash ./PerRead_Analysis/PerRead_Table/Generate_perread_tables.sh Per_read_Filtered_ReplicableSites.tsv pA_tail.txt isoquant.tsv Filtered

```

- Expected output:
```bash
TO FILL IN
```
