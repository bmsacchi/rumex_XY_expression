#!/bin/bash

featureCounts -v

#featureCounts -p --countReadPairs -t exon -g gene_id -a genome/merged_TX_noMatPAR_main.gtf -o readCounts/countsFrankenTX.txt RNAbams/merged_TX_noMatPAR/Rh_02_RNA_S122_L004_Aligned.sortedByCoord.out.bam

featureCounts -p --countReadPairs -t transcript -f -g gene_id -a genome/merged_TX_noMatPAR_main.gtf -o readCounts/countsFrankenTX_genes_noMulti RNAbams/merged_TX_noMatPAR/Rh_02_RNA_S122_L004_Aligned.sortedByCoord.NoMulti.out.bam
  
#featureCounts -p -t transcript -g gene_id -a genome/merged_TX_noMatPAR_main.gtf -o readCounts/countsFrankenTX_fragments RNAbams/merged_TX_noMatPAR/Rh_02_RNA_S122_L004_Aligned.sortedByCoord.out.bam

