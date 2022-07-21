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
module load StdEnv/2020 gcc/9.3.0 samtools/1.13 java/11.0.2 r/4.1.2
PICARD="$SRC/picard/picard.jar"
MM2="$SRC/minimap2/minimap2-2.24_x64-linux/minimap2"

cd $CWD

ref_fasta="$REF/genome/gencode.v${GENCODE_VERSION}.primary_assembly.genome.fa"
ref_tr_fasta="$REF/transcriptome/gencode.v${GENCODE_VERSION}.transcripts.fa"

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

            #### Align to Genome
            $MM2 \
                -t $NUM_THREADS \
                -a \
                -x splice:hq \
                -uf \
                --secondary=no \
                ${ref_fasta} \
                ${sn}.clustered.hq.fasta.gz \
                > ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sam

                ## -x splice:hq : -xsplice -C5 -O6,24 -B4

            samtools view \
                -b \
                -T ${ref_fasta} \
                ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sam \
                > ${sn}.clustered.hq.aln-gene.mm2.splice-hq.bam 

            samtools sort \
                -@ $NUM_THREADS \
                ${sn}.clustered.hq.aln-gene.mm2.splice-hq.bam \
                > ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.bam

            samtools index \
                -@ $NUM_THREADS \
                ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.bam

            ### Filter out supplementary mappings + unmapped
            samtools view \
                -b \
                -F 2308 \
            ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.bam \
                > ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.bam

                ## -F 2308 : read unmapped + non-primary (secondary) + supplementary alignment

            samtools index \
                -@ $NUM_THREADS \
                ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.bam

            #### Align to Transcriptome
            $MM2 \
                -t $NUM_THREADS \
                -a \
                -x map-hifi \
                --secondary=no \
                ${ref_tr_fasta} \
                ${sn}.clustered.hq.fasta.gz \
                > ${sn}.clustered.hq.aln-tr.mm2.map-hifi.sam

                ## -x map-hifi : -k19 -w19 -U50,500 -g10k -A1 -B4 -O6,26 -E2,1 -s200

            samtools view \
                -b \
                -T ${ref_tr_fasta} \
                ${sn}.clustered.hq.aln-tr.mm2.map-hifi.sam \
                > ${sn}.clustered.hq.aln-tr.mm2.map-hifi.bam 

            samtools sort \
                -@ $NUM_THREADS \
                ${sn}.clustered.hq.aln-tr.mm2.map-hifi.bam \
                > ${sn}.clustered.hq.aln-tr.mm2.map-hifi.sorted.bam

            samtools index \
                -@ $NUM_THREADS \
                ${sn}.clustered.hq.aln-tr.mm2.map-hifi.sorted.bam

            ### Filter out supplementary mappings + unmapped
            samtools view \
                -b \
                -F 2052 \
            ${sn}.clustered.hq.aln-tr.mm2.map-hifi.sorted.bam \
                > ${sn}.clustered.hq.aln-tr.mm2.map-hifi.sorted.filtered.bam

                ## -F 2052 : read unmapped + supplementary alignment

            samtools index \
                -@ $NUM_THREADS \
                ${sn}.clustered.hq.aln-tr.mm2.map-hifi.sorted.filtered.bam


            #### BAM to TSV (Extracting Reads' ID and MAPQ)
            ### aln-gene
            samtools view \
                ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.bam \
                | awk '{print $1,$5,length($10)}' \
                    > ${sn}.clustered.hq.aln-gene.mm2.splice-hq.sorted.filtered.tsv

            ### aln-tr
            samtools view \
                ${sn}.clustered.hq.aln-tr.mm2.map-hifi.sorted.filtered.bam \
                | cut -f 1,5 \
                    > ${sn}.clustered.hq.aln-tr.mm2.map-hifi.sorted.filtered.tsv
    
        done
    done
done
