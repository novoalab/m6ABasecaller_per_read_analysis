module load Singularity/3.2.1
singularity exec /users/enovoa/scruciani/soft/modPhred/modphred-3.6.1.sif /users/enovoa/scruciani/soft/modPhred/run -f /no_backup_isis/enovoa/reference_fasta/mouse_genome/GRCm38.p6.genome.fa -o /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/ModPhred -i /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/fast5_files/CTR1 /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/fast5_files/CTR2 /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/fast5_files/INH1 /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/fast5_files/INH2

## directly run also modprob01
mkdir /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/ModPhred/modprob01
mkdir /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/ModPhred/modprob01/minimap2

cp /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/ModPhred/modPhred.pkl /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/ModPhred/modprob01/modPhred.pkl

ln -s /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/ModPhred/minimap2/*.ba* /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/ModPhred/modprob01/minimap2/

singularity exec /users/enovoa/scruciani/soft/modPhred/modphred-3.6.1.sif /users/enovoa/scruciani/soft/modPhred/src/mods_from_bams.py -f /no_backup/enovoa/reference_fasta/GRCm38.p6.genome.fa --minModProb 0.1 -o /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/ModPhred/modprob01 -i /no_backup/enovoa/nextflow_outputs/METTL3-inhibitor/mESC_2reps_1stpooled/ModPhred/modprob01/minimap2/*.bam
