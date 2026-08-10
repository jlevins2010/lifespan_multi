#!/bin/bash
#BSUB -q normal
#BSUB -M 50GB
#BSUB -R "rusage[mem=50GB]"
#BSUB -J salsa_step6
#BSUB -e err_step6.err
#BSUB -o out_step6.out
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
HOST_WORK_DIR=${project}/WASP_ATAC/

apptainer exec --cleanenv --writable-tmpfs -B ${HOST_WORK_DIR}:/WASP_ATAC/ ${sif_path} bash ${project}/step6_wasp.sh --library_id ${sample} --inputvcf ${project}/phasing/${sample}.pass.merged.hcphase.vcf.gz --inputbam ${HOST_WORK_DIR}/${sample}.bcfilter.bam --outputdir ${HOST_WORK_DIR} --outputbam ${sample}.hcphase.wasp.bam --atacref ${cellranger_ref} --genotype joint --modality atac --isphased true
