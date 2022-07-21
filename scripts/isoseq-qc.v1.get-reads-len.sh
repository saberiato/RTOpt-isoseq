## Copyright 2022, iso-seq team
## isoseq-pipeline
## version 2.0

#!/bin/bash
set -eo pipefail

CWD="/scratch/saberi/isoseq/RTOpt"
REF="$CWD/references"
ENV="/scratch/saberi/envs"
SRC="/scratch/saberi/isoseq/scripts"

NUM_THREADS=48
GENCODE_VERSION=40

# source $ENV/init-conda.sh
module load StdEnv/2020 gcc/9.3.0 samtools/1.13

cd $CWD

for s in L M N
do
    for v in 1 2
    do
        for r in 1 2
        do

            sn=${s}.v${v}.r${r}
            echo ${sn}
            SAMPLE_DIR=${CWD}/samples/${sn}
            if [ ! -d ${SAMPLE_DIR} ]; then continue; fi

            cd ${SAMPLE_DIR}

            samtools view \
                ${sn}.fl.bam \
                | awk '{print $1,length($10)}' \
                    > ${sn}.fl.read-len.list

            samtools view \
                ${sn}.flnc.bam \
                | awk '{print $1,length($10)}' \
                    > ${sn}.flnc.read-len.list

            samtools view \
                ${sn}.clustered.hq.bam \
                | awk '{print $1,length($10)}' \
                    > ${sn}.clustered.hq.read-len.list

            samtools view \
                ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.bam \
                | awk '{print $1,length($10)}' \
                    > ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.read-len.list

        done
    done
done