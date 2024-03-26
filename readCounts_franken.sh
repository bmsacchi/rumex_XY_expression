#!/bin/bash

featureCounts -v

#featureCounts -p --countReadPairs -t exon -g gene_id -a genome/merged_TX_noMatPAR_main.gtf -o readCounts/countsFrankenTX.txt RNAbams/merged_TX_noMatPAR/Rh_02_RNA_S122_L004_Aligned.sortedByCoord.out.bam

#featureCounts -p --countReadPairs -t transcript -f -g gene_id -a genome/merged_TX_noMatPAR_main.gtf -o readCounts/countsFrankenTX_genes_noMulti RNAbams/merged_TX_noMatPAR/Rh_02_RNA_S122_L004_Aligned.sortedByCoord.NoMulti.out.bam
  
#featureCounts -p -t transcript -g gene_id -a genome/merged_TX_noMatPAR_main.gtf -o readCounts/countsFrankenTX_fragments RNAbams/merged_TX_noMatPAR/Rh_02_RNA_S122_L004_Aligned.sortedByCoord.out.bam


## don't be silly! summarise by exon! the last two columns are what i want, ignore the other bits!

## removing multi prior doesn't change anything
# freaturcounts only count those if you tell it to

#featureCounts -T 10 -p --countReadPairs -t exon -g gene_id -a genome/merged_TX_noMatPARlarge_txanno.gtf -o readCounts/countsFrankenTXannoNoMulti RNAbams/merged_TX_noMatPAR/Rh_02_RNA_S122_L004_Aligned.sortedByCoord.NoMulti.out.bam 

featureCounts -T 10 -p --countReadPairs -t exon -g gene_id -a genome/merged_TX_noMatPARlarge_txanno.gtf -o readCounts/countsFrankenTXanno RNAbams/merged_TX_noMatPAR/Rh_02_RNA_S122_L004_Aligned.sortedByCoord.out.bam
