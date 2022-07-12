# RTOpt IsoSeq

**Data**: (`graham:/scratch/hanig/RTOpt-isoseq`)
* `raw/`: Raw files (output of demultiplexing step).
* `samples/`: Renamed and processed files (through IsoSeq pipeline).
* `references/`: Reference genome and transcriptome for alignment. 

**Scripts**:
* `rename-samples.sh`: Rename sample ids to more convenient names.
* `isoseq-pipeline.sh`: Apply PacBio IsoSeq pipeline on `raw` data.
* `isoseq-qc.sh`: QC of IsoSeq pipeline based on alignment to reference genome and transcriptome.
