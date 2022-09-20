################################################
### Script to generate per read-level tables ###
################################################

#Import arguments:
m6a=$1
pA_tail=$2
isoforms=$3
output=$4

##Isoform data: 
#Extract transcript information: 
awk -v OFS="\t" '{print $1,$4}' $isoforms > test.isoforms.tsv

#Get assignment type: 
cut -f1,5,6 $isoforms | uniq > assignment_type.tsv

#Concatenate data for reads with multiple assignments (isoforms and/or genes): 
awk -F'\t' -v OFS='\t' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' test.isoforms.tsv | awk -v OFS="," '{first = $1; $1 = ""; print first, $0; }' - | sed 's/,,/\t/g' - > final.isoforms.txt
awk -F'\t' -v OFS='\t' '{x=$1;$1="";a[x]=a[x]$0}END{for(x in a)print x,a[x]}' <(cut -f1,2 assignment_type.tsv) | awk -v OFS="," '{first = $1; $1 = ""; print first, $0; }' - | sed 's/,,/\t/g' - > final.genes.txt

##Output name generation:
out_name="per_read_"$output".tsv"

## Merge all data:
#Generate the table: 
join -a 1 -e "NA" -o auto <(sort $m6a) <(sort $pA_tail | uniq) | join -a 1 -e "NA" -o auto - <(sort final.isoforms.txt) | join -a 1 -e "NA" -o auto -  <(sort final.genes.txt) | join -a 1 -e "NA" -o auto - <(cut -f 1,3 assignment_type.tsv | sort | uniq ) | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$14,$15,$16}' > $out_name

#Add header: 
sed -i '1i read_id\tnum_m6a_sites\tchr\tstrand\tstart_aln\tend_aln\tread_length\tbases_aligned\tsecondary_aln\tpos_m6A_sites\ttailfindr\tnanopolish\tisoform_id\tgene_id\tassignment_type' $out_name

##Clean intermediate files:
rm test.isoforms.tsv
rm assignment_type.tsv
rm final.genes.txt
rm final.isoforms.txt
