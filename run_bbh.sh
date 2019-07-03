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

# point to the scripts directory
scripts=/Users/kseniak/Box/Krasileva_Lab/Research/ksenia/Projects/NLR-ID/KEGG/scripts
meme=$scripts/NLR-Parser-master/meme.xml

# Step 0. Prepare data files

# Individual plant proteomes and their domain annotations are obtained from KEGG database.

# point to the KEGG directory
KEGG_files=/Users/kseniak/Box/Krasileva_Lab/Research/ksenia/Projects/NLR-ID/KEGG/genes/organisms

# NLR database is obtained based on presence of NB-ARC domain in protein sequences (PF00931) using Biomart and hmmsearch.

# point to the NLR database file
NLRdb=/Users/kseniak/Dropbox/work/manuscipts/ID_pathways/datasets/NLR_database/SFile1_Phytozome12.1.6_PlantEnsembl43_refseq94.PF00931.protein.c1aS1.fa

#makeblastdb -in $NLRdb -dbtype prot

# Annotate domains in NLR database using pfamscan

#perl /biotools/PfamScan/pfam_scan.pl -e_seq 1 -e_dom 1 -as -outfile $temp/NLRdb_pfamscan.out -cpu 8 -fasta $NLRdb -dir /biotools/biodb/ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0

#nb: double the spaces to match multiple spaces
#egrep '(TIR)|(NB-ARC)|(RPW8)|(LRR)' $temp/NLRdb_pfamscan.out | sed -e 's/  */<TAB>/g' | cut -f 1,2,3 > $temp/NLRdb_pfamscan.bed

#run NLRannotator on NLR database

#split -l 200 $NLRdb $temp/NLRdb.100NLRs.

#run run_mast.sh to define all motifs in NLR database on the split files (because of sequence limit in mast)

#mast -hit_list $meme $NLRdb > $temp/NLRdb_nlrannotator.mast

#cat *_mastNLR.tsv > Pv12_pe_r85.v1.NLR.mast.tsv

#cut -f 1,3,4 -d " " Pv12_pe_r85.v1.NLR.mast.tsv | grep -v "#" | sed 's/ \+/     /g' > Pv12_pe_r85.v1.NLR.mast.bed


# Step 1.subtract genes with canonical NLR domains from plant proteomes

# point to the directory for intermediate files (temp) and results directory
temp=/Users/kseniak/Box/Krasileva_Lab/Research/ksenia/Projects/NLR-ID/KEGG/results/v2.0/temp
results=/Users/kseniak/Box/Krasileva_Lab/Research/ksenia/Projects/NLR-ID/KEGG/results/v2.0/

for i in ath osa bdi mtr nta sly

do

input_fa=$(find $KEGG_files/$i -name "*pep")
echo $input_fa
input_ids=$KEGG_files/$i/${i}_ncbi-proteinid.list
echo $input_ids
pfam_list=$KEGG_files/$i/${i}_pfam.list
echo $pfam_list
kegg_lsit=$KEGG_files/$i/${i}_pathway.list
echo $kegg_list


grep -v -f <(egrep 'NB-ARC|RPW8|TIR' $pfam_list | cut -f 1) $input_ids | sort | uniq > $temp/${i}.nonNLR.ids

#perl $scripts/K-get_fasta_from_ids.pl -i $temp/$i.nonNLR.ids -f $input_fa > $temp/$i.nonNLR.protein.fa


# Step 2.run blast of proteome against NLRs

makeblastdb -in $temp/$i.nonNLR.protein.fa -dbtype prot -parse_seqids

#run blast using following parameters: `blastp -db $db -query $infasta -out $outfile -outfmt 6 -evalue 1e-5`

blastp -db $NLRdb -query $temp/$i.nonNLR.protein.fa -out $temp/${i}_blastp_NLRdb.e-5.tsv -outfmt 6 -evalue 1e-5

cut -f 2,9,10 $temp/${i}_blastp_NLRdb.e-5.tsv > $temp/${i}.hits.bed

done

exit 0

# Step 3. Intersect canonical NLR domains with blast results to filter them out

bedtools intersect -v -a $temp/${i}.hits.bed -b $temp/NLRdb_pfamscan.bed > $temp/${i}.hits.filtered.bed

#filter out additional canonical NLR motifs (CC, LRR)

bedtools intersect -v -a $temp/bdi_hits.filtered.bed -b $temp/NLRdb_nlrannotator.mast.bed > $outdir/bdi_hits.filtered2.bed

#merge blast hsps

bedtools merge -i <(sort -k1,1 -k2,2n -k3,3n $temp/${i}.hits.filtered2.bed ) > $temp/${i}.hits.filtered2.merged.bed

bedtools getfasta -fi $NLRdb -bed $temp/${i}.hits.filtered2.merged.bed -fo $temp/${i}.hits.filtered2.merged.fa

# Step 4. Run the reciprocal blast

blastp -db $temp/$i.nonNLR.protein.fa -query $temp/${i}.hits.filtered2.merged.fa -out $temp/${i}.hits_blastp_$i.nonNLR.e-5.tsv -outfmt 6 -evalue 1e-5

cut -f 1,2 $temp/${i}.hits_blastp_$i.nonNLR.e-5.tsv | sed 's/:.*<TAB>/<TAB>/' > $results/$i.bbp.tsv

# Step 5. Map reciprocal blast hits to KEGG pathways

join -1 1 -2 2 -o '0,1.2,2.1' <(sort -k 1,1 $KEGG_files/${i}_pathway.list) <(sort -k2,2 $results/$i.bbp.tsv) > $results/$i.bbp.pathways.tsv

#convert into a human readable list

perl $scripts/K-reformat-bbp-table.pl -i <(awk -F' ' '{print $2" "$1}' $results/$i_bbp.tsv) > $results/$i_bbp.list.tsv

perl $scripts/K-reformat-bbp-table.pl -i $results/$i.bbp.pathways.tsv > $results/$i.bbp.pathways.list.tsv

