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
SQANTI3_DIR="${SRC}/SQANTI3/SQANTI3-5.0"
CDNA_CUPCAKE_DIR="${SRC}/cDNA_Cupcake/cDNA_Cupcake-28.0.0"

NUM_THREADS=48
NUM_THREADS_ALIGN=36
NUM_THREADS_SORT=12
NUM_THREADS_SQANTI3=8
GENCODE_VERSION=40

source $ENV/init-conda.sh

cd $CWD

module load StdEnv/2020 gcc/9.3.0 samtools/1.13 java/11.0.2 r/4.1.2
# conda activate $ENV/isoseq_env

MM2="$SRC/minimap2/minimap2-2.24_x64-linux/minimap2"

ref_fasta="$REF/genome/gencode.v${GENCODE_VERSION}.primary_assembly.genome.fa"

for s in L M N
do
    for v in 1 2
    do
        for r in 1 2
        do
            for type in fl flnc
            do

                sn=${s}.v${v}.r${r}
                echo ${sn}.${type}

                SAMPLE_DIR=${CWD}/samples/${sn}
                if [ ! -d ${SAMPLE_DIR} ]; then continue; fi

                cd ${SAMPLE_DIR}

                if [ ${sn}.${type} == "L.v1.r1.fl" ]; then continue; fi

                isoseq3 cluster \
                    ${sn}.${type}.aln-gene.mm2.splice-hq.sorted.filtered.bam \
                    ${sn}.${type}.clustered.bam \
                    --use-qvs \
                    --verbose \
                    --num-threads $NUM_THREADS

                ### BAM to FASTQ
                samtools bam2fq \
                    -@ $NUM_THREADS \
                    ${sn}.${type}.clustered.hq.bam \
                        > ${sn}.${type}.clustered.hq.fastq

                ### Align to Genome
                $MM2 \
                    -t $NUM_THREADS \
                    -a \
                    -x splice:hq \
                    -uf \
                    --secondary=no \
                    ${ref_fasta} \
                    ${sn}.${type}.clustered.hq.fastq \
                    > ${sn}.${type}.clustered.hq.aln-gene.mm2.splice-hq.sam

                  ## -x splice:hq : -xsplice -C5 -O6,24 -B4

                samtools view \
                    -b \
                    -T ${ref_fasta} \
                    ${sn}.${type}.clustered.hq.aln-gene.mm2.splice-hq.sam \
                    > ${sn}.${type}.clustered.hq.aln-gene.mm2.splice-hq.bam 

                samtools sort \
                    -@ $NUM_THREADS \
                    ${sn}.${type}.clustered.hq.aln-gene.mm2.splice-hq.bam \
                    > ${sn}.${type}.clustered.hq.aln-gene.mm2.splice-hq.sorted.bam

                samtools index \
                    -@ $NUM_THREADS \
                    ${sn}.${type}.clustered.hq.aln-gene.mm2.splice-hq.sorted.bam

                pbmm2 align \
                    ${ref_fasta} \
                    ${sn}.${type}.clustered.hq.bam \
                    ${sn}.${type}.clustered.hq.aln-gene.pbmm2.isoseq.sorted.bam \
                    --preset ISOSEQ \
                    --sort \
                    --log-level INFO \
                    --num-threads $NUM_THREADS_ALIGN \
                    --sort-threads $NUM_THREADS_SORT

                ### Collapsing aligned transcripts to generate a transcriptome --
                isoseq3 collapse \
                    ${sn}.${type}.clustered.hq.aln-gene.pbmm2.isoseq.sorted.bam \
                    ${sn}.${type}.transcriptome.gff \
                    --num-threads $NUM_THREADS

            done
        done    
    done
done

conda activate $ENV/sqanti3_env
export PYTHONPATH=$PYTHONPATH:$CDNA_CUPCAKE_DIR/sequence/
export PYTHONPATH=$PYTHONPATH:$CDNA_CUPCAKE_DIR/

ref_fasta="${REF}/genome/gencode.v${GENCODE_VERSION}.primary_assembly.genome.fa"
ref_gtf="${REF}/annotation/gencode.v${GENCODE_VERSION}.primary_assembly.annotation.gtf"
cage_peaks_bed="${REF}/cage-peaks/refTSS_v3.3_human_coordinate.hg38.gencode.bed"
polyA_motif_list="${REF}/polyA-peaks/human.polyA.list.txt"
polyA_peaks_bed="${REF}/polyA-peaks/atlas.clusters.2.0.GRCh38.96.gencode.bed"
intropolis_bed="${REF}/intropolis/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified"
tappas_ref_gff="${REF}/tappas/tappas.Homo_sapiens_GRCh38_Ensembl_86.gff3"

for s in N # L M
do
    for v in 1 # 2
    do
        for r in 2 # 1
        do
            for type in fl flnc
            do

                sn=${s}.v${v}.r${r}
                echo ${sn}.${type}

                SAMPLE_DIR=${CWD}/samples/${sn}
                if [ ! -d ${SAMPLE_DIR} ]; then continue; fi
                if [ ${sn}.${type} == "L.v1.r1.fl" ]; then continue; fi

                cd ${SAMPLE_DIR}

                mkdir -p SQANTI3

                ### Generating SQANTI3 QC report for the transcriptomes -----------------
                python3 ${SQANTI3_DIR}/sqanti3_qc.py \
                    ${sn}.${type}.transcriptome.gff \
                    ${ref_gtf} \
                    ${ref_fasta} \
                    --dir SQANTI3/ \
                    --output ${sn}.${type}.transcriptome.sqanti3 \
                    --fl_count ${sn}.${type}.transcriptome.abundance.txt \
                    --CAGE_peak ${cage_peaks_bed} \
                    --polyA_motif_list ${polyA_motif_list} \
                    --polyA_peak ${polyA_peaks_bed} \
                    --coverage ${intropolis_bed} \
                    --chunks ${NUM_THREADS_SQANTI3} \
                    --report skip

                    # --isoAnnotLite \
                    # --gff3 ${tappas_ref_gff} \
                    # --short_reads ${sn}.sqanti3.short_reads.fofn

                ### Filtering low-quality transcripts based on SQANTI3 criteria ----------
                ## 	The current filtering rules are as follow:
                ##		If a transcript is FSM, then it is kept unless the 3' end is unreliable (intrapriming).
                ##		If a transcript is not FSM, then it is kept only if all of below are true:
                ##			(1) 3' end is reliable.
                ##			(2) does not have a junction that is labeled as RTSwitching.
                ##			(3) all junctions are either canonical or has short read coverage above -c threshold.

                # python3 ${SQANTI3_DIR}/sqanti3_filter.py \
                #     ${sn}.${type}.transcriptome.sqanti3_classification.txt \
                #     ${sn}.${type}.transcriptome.sqanti3_corrected.fasta \
                #     ${sn}.${type}.transcriptome.sqanti3_corrected.gtf
        
            done
        done    
    done
done


for s in L M N
do
    for v in 1 2
    do
        for r in 1 2
        do
            for type in fl flnc
            do

                sn=${s}.v${v}.r${r}
                echo ${sn}.${type}

                SAMPLE_DIR=${CWD}/samples/${sn}
                if [ ! -d ${SAMPLE_DIR} ]; then continue; fi

                cd ${SAMPLE_DIR}

                samtools view \
                    ${sn}.${type}.clustered.hq.bam \
                    | cut -f 1 \
                        > ${sn}.${type}.clustered.hq.qname.list

            done
        done
    done
done
