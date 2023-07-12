# Analysis of m6A sites at per read level

The m6ABasecaller allows the analysis of m6A data at per single molecule level. Before proceeding with this analysis, please first perform the analysis at per position level [here](https://github.com/novoalab/m6ABasecaller/tree/main/PerPosition_Analysis). 

Prior to running this section, please make sure that you have: 

- The bam files obtained as output of ModPhred pipeline in the minimap2 folder. In order to avoid artefacts due to mapping ambiguities, please filter these bams so that they only contain unique and primary alignments: 

```bash
samtools view -Sb -h  -F 3844 pre-filtered.bam > unique_primary_reads.bam
```

- A table with polyA length information, calculated with tailfindr and nanopolish. This table can be obtained with the mop_tail module of Master of Pores (see https://github.com/biocorecrg/MOP2) and will be in mop_tail/polya_common/fast5_joined.txt. 

- A table with Isoquant read assignment to isoforms. Isoquant (https://github.com/ablab/IsoQuant) should be run on the filtered bam files, we suggest to  use the latest annotation for your organism and to use the following command: 

```bash
isoquant.py --reference same_ref_as_the_alignment.fa --genedb annotation.gtf --stranded forward --complete_genedb --no_secondary --count_exons --bam unique_primary_reads.bam --data_type nanopore -o samplename
```
This command will output many tables, the one needed for the per read analysis will be in Isoquant output in 00_samplename/00_samplename.read_assignments.tsv


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
  -gtf GTF, --gtf GTF   Path to the annotation file - it should only contain
                        exon data (*.gtf).
  -o OUTPUT, --output OUTPUT
                        Output name.

```

- Command line example:
```bash
python ./ModificationFrequency_Analysis/ModFreq_Tables.py -i ./example_output/per_read_reannotated_unique_fsm_monoex_test.tsv -gtf ./example_input/Homo_sapiens.GRCh38.109.chr_annotation_chrM.sorted.exons.gtf -o test -rean
```

- Expected output:
```bash
head test_ModFreq_PerGene.tsv 

Sample	GeneID	Coverage	m6A_Site	ModFreq
test_data	ENSG00000133112	1904	45336884	0.12920168067226892
test_data	ENSG00000133112	1904	45337046	0.10504201680672269
test_data	ENSG00000133112	1904	45337294	0.3177521008403361
test_data	ENSG00000133112	1904	45337310	0.11922268907563026
test_data	ENSG00000196136	1019	94614605	0.20608439646712462


head test_ModFreq_PerIsoform.tsv 

Sample	GeneID	IsoformID	Coverage	m6A_Site	ModFreq
test_data	ENSG00000133112	ENST00000530705_0	434	45336884	0.5668202764976958
test_data	ENSG00000133112	ENST00000530705_0	434	45337046	0.4608294930875576
test_data	ENSG00000133112	ENST00000530705_0	434	45337294	0.35023041474654376
test_data	ENSG00000133112	ENST00000530705_0	434	45337310	0.15898617511520738
test_data	ENSG00000133112	ENST00000530705_1	1470	45337294	0.3081632653061224
```

### 3.2. Calculate difference in modification frequency (ΔModFreq):
To calculate ΔModFreq between isoforms, please use the following script: `./ModificationFrequency_Analysis/ModFreq_Difference.py`

```bash
usage: ModFreq_Difference.py [-h] [-i INPUT] [-cov COVERAGE] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input file containing modification
                        frequencies at per isoform or gene level.
  -cov COVERAGE, --coverage COVERAGE
                        Minimum coverage to include an isoform/gene in the
                        analysis. Default = 40
  -o OUTPUT, --output OUTPUT
                        Output name.
```

- Command line example:
```bash
python ./ModificationFrequency_Analysis/ModFreq_Difference.py -i ./example_output/test_ModFreq_PerIsoform.tsv -o test
```

- Expected output:
```bash
head DiffModFreq_PerIsoform.tsv 

Sample  GeneID  Isoform_1       Isoform_2       Site    ModFreq_I1      ModFreq_I2      ΔModFreq(I1-I2)
test_data       ENSG00000133112 ENST00000530705_0       ENST00000530705_1       45337310        0.15898617511520738     0.10748299319727893     0.05150318191792845
test_data       ENSG00000133112 ENST00000530705_0       ENST00000530705_1       45337294        0.35023041474654376     0.3081632653061224      0.042067149440421336
```

## 4. Co-occurance analysis

### 4.1. Calculate standard deviation from expected values:
To assess if a pair of sites co-occur or are mutually exclusive, it is needed to calculate the standard deviation from expected values using the script `./Cooccurance_Analysis/cooccurance_analysis.py`.

```bash
usage: cooccurance_analysis.py [-h] [-i INPUT] [-o OUTPUT] [-cov COVERAGE]
                               [-min_exp MIN_EXPECTED] [-rean]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input file containing per read label data.
  -o OUTPUT, --output OUTPUT
                        Output name.
  -cov COVERAGE, --coverage COVERAGE
                        Minimum coverage to include a transcript in the
                        analysis. Default = 200
  -min_exp MIN_EXPECTED, --min_expected MIN_EXPECTED
                        Minimum expected counts to include a pairwise
                        comparison in the analysis. Default = 2
  -rean, --reannotated  Input is per read data reannotated by annotate_UTRs.py
```

- Command line example:
```bash
python ./Cooccurance_Analysis/cooccurance_analysis.py -i ./example_output/per_read_reannotated_unique_fsm_monoex_test.tsv -o test -rean
```

- Expected output:
```bash
head test_CoocuranceAnalysis.tsv

Transcript	Site A	Site B	Transcript_counts	Counts_A	FreqA	Counts_B	FreqB	Counts_AB	obs_FreqAB	exp_FreqAB	SdFromExpected
ENST00000530705_1	45337310	45337294	1470	158	0.10748299319727891	453	0.3081632653061224	62	0.04217687074829932	0.03312231014854922	1.939902790252726
ENST00000393078_0	94614956	94614605	1019	429	0.4210009813542689	210	0.20608439646712462	93	0.09126594700686948	0.08676173315446169	0.5107988770871357
ENST00000393078_0	94614956	94614893	1019	429	0.4210009813542689	96	0.09421000981354269	45	0.04416094210009813	0.03966250658489678	0.7357776138717788
ENST00000393078_0	94614956	94614769	1019	429	0.4210009813542689	81	0.07948969578017664	31	0.03042198233562316	0.033465239931006655	-0.5401571091944367
ENST00000393078_0	94614956	94615034	1019	429	0.4210009813542689	71	0.06967615309126594	31	0.03042198233562316	0.029333728828413237	0.20587260578752073
```

### 4.2. Analyse the distribution of standard deviation from expected values:
Please use the script: `./Cooccurance_Analysis/cooccurance_density.py`.

```bash
usage: cooccurance_density.py [-h] [-i INPUT] [-o OUTPUT] [-cov COVERAGE]
                              [-min_exp MIN_EXPECTED]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Path to the input file containing co-occurance data.
  -o OUTPUT, --output OUTPUT
                        Output name.
  -cov COVERAGE, --coverage COVERAGE
                        Minimum coverage to include a transcript in the
                        analysis. Default = 200
  -min_exp MIN_EXPECTED, --min_expected MIN_EXPECTED
                        Minimum expected counts to include a pairwise
                        comparison in the analysis. Default = 2
```

- Command line example:
```bash
python ./Cooccurance_Analysis/cooccurance_density.py -i ./example_output/test_CoocuranceAnalysis.tsv -o test
```

- Expected output:

![Distribution_comparisons](https://github.com/novoalab/m6ABasecaller/blob/main/img/Density_Comparison.png)

```bash
MannwhitneyuResult(statistic=261.0, pvalue=0.2923883610418598)
```

![Sd_GenDistance](https://github.com/novoalab/m6ABasecaller/blob/main/img/GenomicDistance_Cooccurence.png)







