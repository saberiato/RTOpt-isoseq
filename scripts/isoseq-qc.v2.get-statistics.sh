## Copyright 2022, iso-seq team
## isoseq-pipeline
## version 2.0

#!/bin/bash
set -eo pipefail

ISOSEQ_WD="/scratch/saberi/isoseq"
CWD="${ISOSEQ_WD}/RTOpt"
ENV="/scratch/saberi/envs"

REF="${ISOSEQ_WD}/references"
SRC="${ISOSEQ_WD}/scripts"

NUM_THREADS=48
GENCODE_VERSION=40

# source $ENV/init-conda.sh

cd $CWD

module load StdEnv/2020 gcc/9.3.0 samtools/1.13
conda activate $ENV/isoseq_env

for s in L M N
do
    for v in 1 2
    do
        for r in 1 2
        do
            for type in fl flnc clustered.hq
            do

                sn=${s}.v${v}.r${r}
                echo ${sn}.${type}

                SAMPLE_DIR=${CWD}/samples/${sn}
                if [ ! -d ${SAMPLE_DIR} ]; then continue; fi

                cd ${SAMPLE_DIR}

                ### RAW
                samtools view \
                    ${sn}.${type}.bam \
                    | awk '{print $1,$5,length($10)}' \
                        > ${sn}.${type}.tsv

                ### ALIGNED
                samtools view \
                    ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.bam \
                    | awk '{print $1,$5,length($10)}' \
                        > ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.tsv

                ### ALIGNED HQ
                samtools view \
                    ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.filtered.bam \
                    | awk '{print $1,$5,length($10)}' \
                        > ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.filtered.tsv

            done
        done
    done
done