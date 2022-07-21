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
MM2="$SRC/minimap2/minimap2-2.24_x64-linux/minimap2"

cd $CWD

ref_fasta="$REF/genome/gencode.v${GENCODE_VERSION}.primary_assembly.genome.fa"

MIN_MAPQ="50"

qc () {

    sn=$1
    type=$2
    
    ### BAM to FASTQ
    samtools bam2fq \
        -@ $NUM_THREADS \
        ${sn}.${type}.bam \
            > ${sn}.${type}.fastq

    ### Align to Genome
    $MM2 \
        -t $NUM_THREADS \
        -a \
        -x splice:hq \
        -uf \
        --secondary=no \
        ${ref_fasta} \
        ${sn}.${type}.fastq \
        > ${sn}.${type}.aln-gene.mm2.splice-hq.sam

        ## -x splice:hq : -xsplice -C5 -O6,24 -B4

    samtools view \
        -b \
        -T ${ref_fasta} \
        ${sn}.${type}.aln-gene.mm2.splice-hq.sam \
        > ${sn}.${type}.aln-gene.mm2.splice-hq.bam 

    samtools sort \
        -@ $NUM_THREADS \
        ${sn}.${type}.aln-gene.mm2.splice-hq.bam \
        > ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.bam

    samtools index \
        -@ $NUM_THREADS \
        ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.bam

    ### Filter out supplementary mappings + unmapped
    samtools view \
        -b \
        -F 2308 \
        -q ${MIN_MAPQ} \
        ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.bam \
        > ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.filtered.bam

        ## -F 2308 : read unmapped + non-primary (secondary) + supplementary alignment

    samtools index \
        -@ $NUM_THREADS \
        ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.filtered.bam

    ### BAM to TSV (Extracting Reads' ID, MAPQ, and length)
    samtools view \
        ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.filtered.bam \
        | awk '{print $1,$5,length($10)}' \
            > ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.filtered.tsv

}

for s in L M N
do
    for v in 1 2
    do
        for r in 1 2
        do
            for type in fl flnc clustered # clustered.hq
            do

                sn=${s}.v${v}.r${r}
                echo ${sn}.${type}

                SAMPLE_DIR=${CWD}/samples/${sn}
                if [ ! -d ${SAMPLE_DIR} ]; then continue; fi

                cd ${SAMPLE_DIR}

                qc $sn $type

            done
        done
    done
done
