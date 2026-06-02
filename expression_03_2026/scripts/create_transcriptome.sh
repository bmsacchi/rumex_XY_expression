#!/bin/bash




gffread merged_TX_noMatPARlarge_txanno.gtf  \
	-g ../genome/merged_TX_noMatPAR_main.fa \
	-w tx_nomatpar_transcripts.fa



