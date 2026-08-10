#!/bin/bash
#BSUB -q normal
#BSUB -M 50GB
#BSUB -R "rusage[mem=50GB]"
#BSUB -J salsa_step5
#BSUB -e salsa_step5.err
#BSUB -o salsa_step5.out
#BSUB -n 20

sample=FK4
samplepath=/home/levinsj/MultiOme/Fetal/cellranger-arc/FK4/outs/atac_possorted_bam.bam
sif_path=/home/levinsj/MultiOme/Analysis/sif/salsa_latest.sif
cellranger_ref=/home/levinsj/Applications/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
gatk_bundle_path=/home/levinsj/Applications/gatk
reference=/home/levinsj/Applications
project=/home/levinsj/MultiOme/Analysis/SALSA
barcodes=/home/levinsj/MultiOme/Barcodes/final_annotated/barcodes_for_salsa.csv
modailty=atac
HOST_WORK_DIR=${project}/WASP_ATAC


apptainer exec --cleanenv --writable-tmpfs -B ${HOST_WORK_DIR}:/WASP_ATAC/,/home/levinsj,${project} ${sif_path} bash ${project}/step5_filterbam.sh \
        --library_id ${sample} \
        --inputbam ${samplepath} \
        --modality atac \
        --barcodes ${barcodes} \
        --outputdir ${HOST_WORK_DIR}/WASP_ATAC \
        --outputbam ${sample}.bcfilter.bam
