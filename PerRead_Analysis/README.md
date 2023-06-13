# Analysis of m6A sites at per read level

The m6A basecaller allows the analysis of m6A data at per single molecule level. Before proceeding with this analysis, please first perform the analysis at per position level here (MISSING LINK). 

MISSING: how to run mop_tail and isoquant!!

## 1. Generating per-read level tables with m6A, polyA tail and isoform information

### 1.1. Parse results at per read level from the m6A basecaller:

To parse the bam file containing m6A probabilities encoded, use the script `./PerRead_Table/get_mods.py`: 

```bash
Usage: get_mods.py -i <input bam file> -o <OUTPUT PREFIX>

Options:
  -h, --help            show this help message and exit
  -i INPUT, --input=INPUT
                        Input bam file
  -q MINQUAL, --min_qual=MINQUAL
                        Mininum quality
  -o WOTUS, --output=WOTUS
                        Ouput prefix File
  -c CHR, --chromosome=CHR
                        Chromosome or all
  -s SIMPLE, --simple=SIMPLE
                        yes or no
```

- Command line example:

```bash
#Parse the bam file:
python ./PerRead_Table/get_mods.py -i ./example_input/example_input.bam -o ./example_output/parse_m6A.tsv -q 15 -c all -s yes

#Update strand information with awk:
awk '{sub(/True/, "-", $4)}1' ./example_output/parse_m6A.tsv | awk '{sub(/False/, "+", $4)}1' - > ./example_output/parse_strand_m6A.tsv
```

- Expected output: 
```bash
head ./example_output/parse_strand_m6A.tsv

readname N_mod ref_name is_reverse start end readlen alnlen is_secondary mod_list
c297b2e7-baf1-4678-9a07-9eeee397e3ff 0 chr1 + 629651 630666 904 881 False NA
0db3c821-2930-4107-a4e5-4e3a3ad09d68 1 chr1 + 632766 633436 638 612 False 633101
de5af892-650c-4052-9d04-36cb9b49af76 0 chr1 + 632766 633443 622 593 False NA
9f40bed2-2428-4c22-8254-ff702e9d8e92 0 chr1 + 632766 633442 671 610 False NA
5ed78e6a-2f26-4346-82f6-98e3f95ad5b5 0 chr1 + 632766 633440 675 602 False NA
```

### 1.2. OPTIONAL - Filter out non-replicable sites from the m6A data at per read level:

If you only want to keep at per read level those sites that are replicable at per position level, which is **recommended**, please use the script `./PerRead_Table/filter_mod_sites.py`: 

```bash
usage: filter_mod_sites.py [-h] [-pr PER_READ] [-rs REPLICABLE_SITES]
                           [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -pr PER_READ, --per_read PER_READ
                        Path to the input file containing per read
                        modification data.
  -rs REPLICABLE_SITES, --replicable_sites REPLICABLE_SITES
                        Path to the input file containing replicable sites
                        (*_SummaryData_ReplicableSites.tsv).
  -o OUTPUT, --output OUTPUT
                        Output name.
```

- Command line example:

```bash
python ./../PerRead_Table/filter_mod_sites.py -pr ./example_output/parse_strand_m6A.tsv -rs ./example_input/PerPosition_SummaryData_ReplicableSites.tsv -o m6A

```

- Expected output:
```bash
head ./example_output/Per_read_m6A_ReplicableSites.tsv

readname	N_mod	ref_name	is_reverse	start	end	readlen	alnlen	is_secondary	mod_list
c297b2e7-baf1-4678-9a07-9eeee397e3ff	0	chr1	+	629651	630666	904	881	False	
0db3c821-2930-4107-a4e5-4e3a3ad09d68	1	chr1	+	632766	633436	638	612	False	633101
de5af892-650c-4052-9d04-36cb9b49af76	0	chr1	+	632766	633443	622	593	False	
9f40bed2-2428-4c22-8254-ff702e9d8e92	0	chr1	+	632766	633442	671	610	False	
5ed78e6a-2f26-4346-82f6-98e3f95ad5b5	0	chr1	+	632766	633440	675	602	False
```

### 1.3. Generate final table with m6A sites, isoform and polyA tail data for every read in the sample:

To generate the complete per read table, containing m6A, polyA tail length and isoform features for every read, please use `./PerRead_Table/Generate_perread_tables.sh`: 

```bash
usage: Generate_perread_tables.sh m6a_data.tsv polyA_data.tsv isoquant.tsv sample_name output_name

```

- Command line example:

```bash
bash ./../PerRead_Table/Generate_perread_tables.sh ./example_output/Per_read_m6A_ReplicableSites.tsv ./example_input/polyA_tail.tsv./example_input/isoquant.tsv test_data test
```

- Expected output:
```bash
head per_read_test.tsv 

read_id	num_m6a_sites	chr	strand	start_aln	end_aln	read_length	bases_aligned	secondary_aln	pos_m6A_sites	tailfindr	nanopolish	isoform_id	gene_id	assignment_type	assignment_data	sample
000015ba-e6b9-4149-b5db-a2e1277b1a12	0	chr5	-	40832458	40835194	343	334	False	NA	NA	NA	ENST00000274242	ENSG00000145592	unique	fsm	test_data
00002b69-c6b2-417a-8445-d2884f2e5674	0	chr16	+	70346928	70373382	2916	2903	False	NA	70.72	70.72	ENST00000302243	ENSG00000168872	unique	fsm	test_data
0001c401-ef7a-411a-86b1-b31f528f5234	0	chr20	-	35302578	35411960	2186	2185	False	NA	47.35	47.35	ENST00000374385	ENSG00000101019	unique	fsm	test_data
000238d4-4cc9-4b75-aa74-288d69a1160e	0	chr8	-	98041726	98045532	470	467	False	NA	NA	NA	ENST00000287038	ENSG00000156482	unique	fsm	test_data
000277b1-3977-4ea2-82e1-0dd9a7167bb5	0	chr5	-	42799886	42811881	2031	1991	False	NA	111.12	111.12	ENST00000514985	ENSG00000250722	unique	fsm	test_data
```

## 2. Re-annotating 5' and 3' UTR from reads belonging to the same isoform

### 2.1. OPTIONAL - Merge per read level tables from multiple samples:

If we are analysing multiple samples, we should first merge the per read tables from each one of them so the annotation would be common in all of them. Thus, the results will be comparable across samples. Please use the code below to merge tables:

```bash
#Concatenate per read tables from multiple samples:
cat per_read_*.tsv > merged_test.tsv

#Only keep the first header:
sed -i -e '2,${/^read_id/d' -e '}' merged_test.tsv 
```

### 2.2. Re-annotation of 5' and 3' UTR:

The script `./PerIsoform_Analysis/annotate_UTRs.py` re-annotates all reads based on their 5' and 3' UTR coordinates and only reports full-length reads with their updated annotations. 

```bash
usage: annotate_UTRs.py [-h] [-i INPUT] [-o OUTPUT] [-cov COVERAGE]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input file containing per read label data.
  -o OUTPUT, --output OUTPUT
                        Output name.
  -cov COVERAGE, --coverage COVERAGE
                        Minimum coverage to include a transcript in the
                        analysis. Default = 40
```

- Command line example:
```bash
python ./PerIsoform_Analysis/annotate_UTRs.py -i ./example_output/per_read_test.tsv -o test

```

- Expected output:
```bash
head per_read_reannotated_unique_fsm_monoex_test.tsv

read_id	num_m6a_sites	chr	strand	start_aln	end_aln	read_length	bases_aligned	secondary_aln	pos_m6A_sites	tailfindr	nanopolish	isoform_id	gene_id	assignment_type	assignment_data	sample	upd_transcript_id	upd_transcript_start	upd_transcript_end
0008e625-06d8-4173-b112-63d8065ba1ee	0	chr15	+	69452831	69455545	493	475	False	NA	33.64	33.64	ENST00000260379	ENSG00000137818	unique	fsm	test_data	ENST00000260379_0	69452830	69455544
003414cf-180f-4dfa-8310-e331900d9b32	0	chr15	+	69452831	69455545	489	486	False	NA	38.63	38.63	ENST00000260379	ENSG00000137818	unique	fsm	test_data	ENST00000260379_0	69452830	69455544
003c1384-0f99-415b-9ce4-7bbad6e96638	0	chr15	+	69452829	69455545	481	477	False	NA	85.53	85.53	ENST00000260379	ENSG00000137818	unique	fsm	test_data	ENST00000260379_0	69452830	69455544
00433ec9-4acf-492f-96a0-324d5b64df5c	0	chr15	+	69452829	69455545	486	484	False	NA	91.06	91.06	ENST00000260379	ENSG00000137818	unique	fsm	test_data	ENST00000260379_0	69452830	69455544
005c0d0f-f923-4634-97fc-ed0110b0cfc6	0	chr15	+	69452830	69455545	483	480	False	NA	NA	NA	ENST00000260379	ENSG00000137818	unique	fsm	test_data	ENST00000260379_0	69452830	69455544
```

## 3. Modification frequency analysis

### 3.1. Calculate modification frequency at per isoform and gene level:

Before identifying sites whose modification frequency varies across isoforms, genes or samples, the modification frequency at per isoform and gene level should be calculated with the script `./ModificationFrequency_Analysis/ModFreq_Tables.py`:

```bash
usage: ModFreq_Tables.py [-h] [-i INPUT] [-gtf GTF] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input file containing reannotated per read
                        data (output from annotate_UTRs.py).
  -gtf GTF, --gtf GTF   Path to the annotation file - ideally only containing
                        exon data (*.gtf).
  -o OUTPUT, --output OUTPUT
                        Output name.
```

- Command line example:
```bash
python ./../ModificationFrequency_Analysis/ModFreq_Tables.py -i ./example_output/per_read_reannotated_unique_fsm_monoex_test.tsv -gtf ./example_input/Homo_sapiens.GRCh38.109.chr_annotation_chrM.sorted.exons.gtf -o test
```

- Expected output:
```bash
MISSING!!
```


## 4. Co-occurance analysis




