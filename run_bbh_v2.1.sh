#!/usr/bin/bash

# Author: Ksenia Krasileva
# E-mail: kseniak@berkeley.edu
# Prerequisites:

# KEGG database for species of interest
# BLAST
# HMMER
# Bedtools
# NLRparser (including meme)
# additional custom scripts

# Step 0. Prepare data files

# point to the scripts directory
scripts=/Users/Shared/eli/ID_pathways/scripts

# point to the KEGG directory
# Individual plant proteomes and their domain annotations are obtained from KEGG database.
# KEGG pathways gene assignment files are also from KEGG database

KEGG_files=/Users/Shared/eli/ID_pathways/datasets/kegg

# point to the NLR database file
# NLR database is obtained based on presence of NB-ARC domain in protein sequences (PF00931) using Biomart and hmmsearch.

NLRdb=/Users/Shared/eli/ID_pathways/datasets/NLR_database/SFile1_Phytozome12.1.6_PlantEnsembl43_refseq94.PF00931.protein.c1aS1.fa
NLRdbdir=/Users/Shared/eli/ID_pathways/datasets/NLR_database

# point to the location of the meme file containing NLR motifs
meme=~/Box/Krasileva_Lab/Research/ksenia/Projects/current/NLR-ID/KEGG/scripts/NLR-Parser-master/meme.xml

# point to directory where you want to store intermediate data analyses
temp=/Users/Shared/eli/ID_pathways/datasets/NLR_database/temp

# Create database indexes for blast searches
makeblastdb -in $NLRdb -dbtype prot

# Annotate domains in NLR database using pfamscan

perl /biotools/PfamScan/pfam_scan.pl -e_seq 1 -e_dom 1 -as -outfile NLRdb_pfamscan.out -cpu 8 -fasta $NLRdb -dir /biotools/biodb/ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0

#nb: double the spaces to match multiple spaces

egrep '(TIR)|(NB-ARC)|(RPW8)|(LRR)' NLRdb_pfamscan.out | sed -e 's/  */<TAB>/g' | cut -f 1,2,3 > NLRdb_pfamscan.bed

#run NLRannotator on NLR database

#unwrap fasta

awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $NLRdb > NLRdb_v2.0_unwrapped.protein.fa

split -l 400 $NLRdb $temp/NLRdb.200NLRs.

#run run_mast.sh to define all motifs in NLR database on the split files (because of sequence limit in mast)

for file in $temp/NLRdb.200NLRs.* 
do
mast -hit_list $meme $file > $temp/$file"_mastNLR.tsv"
done

cat ./temp/*_mastNLR.tsv > $NLRdbdir/NLRdb_v2.0_mast.tsv

cut -f 1,3,4 -d " " Pv12_pe_r85.v1.NLR.mast.tsv | grep -v "#" | sed 's/ \+/     /g' > Pv12_pe_r85.v1.NLR.mast.bed

# Step 1. Subtract genes with canonical NLR domains from plant proteomes

# point to the results directory

results=/Users/Shared/eli/ID_pathways/datasets/v2.0

for i in ath osa bdi nta sly

do

input_fa=$(find $KEGG_files/$i -name "*pep")
echo $input_fa
input_ids=$KEGG_files/$i/${i}_ncbi-proteinid.list
echo $input_ids
pfam_list=$KEGG_files/$i/${i}_pfam.list
echo $pfam_list
kegg_list=$KEGG_files/$i/${i}_pathway.list
echo $kegg_list

grep -v -f <(egrep 'NB-ARC|RPW8|TIR' $pfam_list | cut -f 1) $input_ids | sort | uniq > $results/${i}.nonNLR.ids

perl $scripts/K-get_fasta_from_ids.pl -i $temp/$i.nonNLR.ids -f $input_fa > $temp/$i.nonNLR.protein.fa

# Step 2.run blast of proteome against NLRs

makeblastdb -in $results/$i.nonNLR.protein.fa -dbtype prot -parse_seqids

#run blast using following parameters: `blastp -db $db -query $infasta -out $outfile -outfmt 6 -evalue 1e-5`

blastp -db $NLRdb -query $results/$i.nonNLR.protein.fa -out $results/${i}_blastp_NLRdb.e-5.tsv -outfmt 6 -evalue 1e-5 -num_threads 8

cut -f 2,9,10 $results/${i}_blastp_NLRdb.e-5.tsv > $results/${i}.hits.bed

# Step 3. Intersect canonical NLR domains with blast results to filter them out

bedtools intersect -v -a $results/${i}.hits.bed -b $results/NLRdb_pfamscan.bed > $results/${i}.hits.filtered.bed

#filter out additional canonical NLR motifs (CC, LRR)

bedtools intersect -v -a $results/${i}.hits.filtered.bed -b $NLRdbdir/NLRdb_v2.0_mast.bed > $results/${i}.hits.filtered2.bed

#merge blast hsps

bedtools merge -i <(sort -k1,1 -k2,2n -k3,3n $results/${i}.hits.filtered2.bed ) > $results/${i}.hits.filtered2.merged.bed

bedtools getfasta -fi $NLRdb -bed $results/${i}.hits.filtered2.merged.bed -fo $results/${i}.hits.filtered2.merged.fa

# Step 4. Run the reciprocal blast

blastp -db $results/$i.nonNLR.protein.fa -query $results/${i}.hits.filtered2.merged.fa -out $results/${i}.hits_blastp_$i.nonNLR.e-5.tsv -outfmt 6 -evalue 1e-5

cut -f 1,2 $results/${i}.hits_blastp_$i.nonNLR.e-5.tsv | sed 's/:.*<TAB>/<TAB>/' > $results/$i.bbp.tsv

# Step 5. Map reciprocal blast hits to KEGG pathways

join -1 1 -2 2 -o '0,1.2,2.1' <(sort -k 1,1 $kegg_list) <(sort -k2,2 $results/$i.bbp.tsv) > $results/$i.bbp.pathways.tsv

#convert into a human readable list

perl $scripts/K-reformat-bbp-table.pl -i <(awk -F' ' '{print $2" "$1}' $results/$i.bbp.tsv) > $results/$i.bbp.list.tsv

perl $scripts/K-reformat-bbp-table.pl -i $results/$i.bbp.pathways.tsv > $results/$i.bbp.pathways.list.tsv

done
