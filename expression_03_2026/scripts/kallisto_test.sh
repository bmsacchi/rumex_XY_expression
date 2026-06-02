#!/bin/bash

kallisto quant -i tx_nomatpar_transcripts.idx -o test_output -t 12 \
	/ohta1/felix.beaudry/rawSequence/RNAseq/josh/pop/10.TM4_R1_clean.fastq.gz /ohta1/felix.beaudry/rawSequence/RNAseq/josh/pop/10.TM4_R2_clean.fastq.gz
