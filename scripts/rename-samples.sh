## Copyright 2021, iso-seq team
## isoseq-pipeline
## version 2.0

#!/bin/bash
set -eo pipefail

#=========================#
#=== ISOSEQ PREPROCESS ===#
#=========================#

####################################
#### Rename Sample Names to IDs ####
####################################

CWD="/scratch/saberi/isoseq/RTOpt"

cd $CWD

### L.v1.r1 <- RTOpt/Lexogen_v1_rep1
### L.v1.r2 <- RTOpt/Lexogen_v1_rep2
### L.v2.r1 <- RTOpt/Lexogen_v2_rep1
### L.v2.r2 <- RTOpt/Lexogen_v2_rep2
### M.v1.r1 <- RTOpt/Marathon_v1_rep1
### M.v1.r2 <- RTOpt/Marathon_v1_rep2
### M.v2.r1 <- RTOpt/Marathon_v2_rep1
### M.v2.r2 <- RTOpt/Marathon_v2_rep2
### N.v1.r1 <- RTOpt/NEB_rep1
### N.v1.r2 <- RTOpt/NEB_rep2

while IFS=' ' read -r project sample id
do 
    mkdir -p samples/${id}
    ln -f raw/${project}/RToptv2.demux*.${sample}--${sample}.bam  samples/${id}/${id}.ccs.bam
    ln -f raw/${project}/RToptv2.demux*.${sample}--${sample}.bam.pbi  samples/${id}/${id}.ccs.bam.pbi
    ln -f raw/${project}/pcr_handle.fa samples/${id}/${id}.primers.fasta

done < scripts/RTOpt.sample-ids.tsv
