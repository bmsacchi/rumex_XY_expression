#!/bin/bash

# confirming that putatively lost genes (present on X, not Y) are 
#transcriptfastsa="/ohta2/Rumex/Dovetail_XY_2023_TX_male/final_scaffolded_assembly/TX_maternal.all.maker.transcripts.fasta"
newfasta="data/TX_maternal_parsed.fa"
outfasta="data/TX_mat_geneloss.fa"
genomefasta="/ohta2/Rumex/Dovetail_XY_2023_TX_male/final_scaffolded_assembly/TX_paternal.fa"
genelist="data/genelist_tx_mat.txt"



#awk -F"-" '/^>/{print $1; next}1' $transcriptfasta > $newfasta


#while read i; do 
#	samtools faidx $newfasta $i >> $outfasta
#done < $genelist

echo "creating genome db"

makeblastdb -in $genomefasta \
	-input_type fasta -dbtype nucl -out patgenome

wait

echo "blasting mat sequences to pat geno"

blastn -db patgenome -query $outfasta \
	-max_hsps 1 \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore score" \
	-out data/blast_pat_genome.tsv



