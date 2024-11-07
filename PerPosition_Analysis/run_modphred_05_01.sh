module load Singularity/3.2.1

## modphred default modprob threshold is 0.5
singularity exec modPhred/modphred-3.6.1.sif modPhred/run -f reference.fa -o ModPhred05 -i input/sample1 input/sample2 input/etc

### unzip and remove first rows of the output so you can further process it
zcat  ModPhred05/mod.gz | tail -n +17 >  ModPhred05/mod_gz.txt 

## from the bam files that are produced, run modphred with 0.1 modprob threshold
mkdir ModPhred01
mkdir ModPhred01/minimap2

cp ModPhred05/modPhred.pkl ModPhred01/modPhred.pkl
ln -s ModPhred05/minimap2/*.ba* ModPhred01/minimap2/

singularity exec modPhred/modphred-3.6.1.sif modPhred/src/mods_from_bams.py -f reference.fa --minModProb 0.1 -o ModPhred01 -i ModPhred01/minimap2/*.bam

### unzip and remove first rows of the output so you can further process it
zcat  ModPhred01/mod.gz | tail -n +17 >  ModPhred01/mod_gz.txt 
