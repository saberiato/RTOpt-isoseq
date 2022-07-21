## Copyright 2022, iso-seq team
## isoseq-pipeline
## version 2.0

#!/bin/bash
set -eo pipefail

#======================#
#=== ISOSEQ PROCESS ===#
#======================#

################################
#### PacBio isoseq pipeline ####
################################

CWD="/scratch/saberi/isoseq/RTOpt"
ENV="/scratch/saberi/envs"
REF="$CWD/references"
NUM_THREADS=48
NUM_THREADS_ALIGN=36
NUM_THREADS_SORT=12
GENCODE_VERSION=40

source $ENV/init-conda.sh
conda activate $ENV/isoseq_env

cd $CWD

ref_mmi="$REF/genome/gencode.v${GENCODE_VERSION}.primary_assembly.genome.mmi"

for s in L M N
do
    for v in 1 2
    do
        for r in 1 2
        do

            sn=${s}.v${v}.r${r}
            echo ${sn}

            if [ ! -d ${CWD}/samples/${sn} ]; then continue; fi

            cd ${CWD}/samples/${sn}

            ref_primers=${sn}.primers.fasta

            ## Removing primers sequences ---------------------------------
            lima \
                ${sn}.ccs.bam \
                ${ref_primers} \
                ${sn}.fl.bam \
                --isoseq \
                --peek-guess \
                --log-level INFO \
                --num-threads $NUM_THREADS

            ## Renaming FL BAM files
            mv ${sn}.fl.*_5p--*_3p.bam ${sn}.fl.bam
            mv ${sn}.fl.*_5p--*_3p.bam.pbi ${sn}.fl.bam.pbi
            mv ${sn}.fl.*_5p--*_3p.consensusreadset.xml ${sn}.fl.primers.consensusreadset.xml

            ### Removing polyA tails and other concatemers ------------------
            isoseq3 refine \
                ${sn}.fl.bam \
                ${ref_primers} \
                ${sn}.flnc.bam \
                --log-level DEBUG \
                --num-threads $NUM_THREADS

            ### Clustering full-length transcripts --------------------------
            isoseq3 cluster \
                ${sn}.flnc.bam \
                ${sn}.clustered.bam \
                --use-qvs \
                --verbose \
                --num-threads $NUM_THREADS_ALIGN

            ### Aligning clusters' representative transcript ----------------
            pbmm2 align \
                ${ref_mmi} \
                ${sn}.clustered.hq.bam \
                ${sn}.clustered.hq.aligned.bam \
                --preset ISOSEQ \
                --sort \
                --log-level INFO \
                --num-threads $NUM_THREADS_ALIGN \
                --sort-threads $NUM_THREADS_SORT

            ### Collapsing aligned transcripts to generate a transcriptome --
            isoseq3 collapse \
                ${sn}.clustered.hq.aligned.bam \
                ${sn}.transcriptome.gff \
                --num-threads $NUM_THREADS

        done
    done   
done
