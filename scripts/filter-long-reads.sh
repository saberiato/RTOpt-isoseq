## Copyright 2022, iso-seq team
## isoseq-pipeline
## version 2.0

#!/bin/bash
set -eo pipefail

CWD="/scratch/saberi/isoseq/RTOpt"
REF="$CWD/references"
ENV="/scratch/saberi/envs"
SRC="/scratch/saberi/isoseq/scripts"

NUM_THREADS=2
GENCODE_VERSION=40

MIN_LEN=5000

# source $ENV/init-conda.sh
module load StdEnv/2020 gcc/9.3.0 samtools/1.13 java/11.0.2
PICARD="$SRC/picard/picard.jar"

cd $CWD

for sn in N.v1.r1 N.v1.r2
do

    cd ${CWD}/samples/${sn}

    awk \
        -v min_len='$MIN_LEN' \
        '$3 > min_len { print $1 }' \
        ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.tsv \
            > ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.long-reads.qname.list


    java -jar $PICARD \
        FilterSamReads \
            I=${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.bam \
            O=${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.long-reads.bam \
            READ_LIST_FILE=${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.long-reads.qname.list \
            FILTER=includeReadList

    samtools index \
        -@ ${NUM_THREADS} \
        ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.long-reads.bam

done
