#!/bin/bash
#BSUB -q normal
#BSUB -M 100GB
#BSUB -R "rusage[mem=100GB]"
#BSUB -J salsa
#BSUB -e errorFiles/salsa_step1.err
#BSUB -o outputFiles/salsa_step1.out
#BSUB -n 20

sample=FK4
samplepath=/home/levinsj/MultiOme/Fetal/cellranger-arc/FK4/outs/atac_possorted_bam.bam
sif_path=/home/levinsj/MultiOme/Analysis/sif/salsa_latest.sif
cellranger_ref=/home/levinsj/Applications/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
gatk_bundle_path=/home/levinsj/Applications/gatk
reference=/home/levinsj/Applications
project=/home/levinsj/MultiOme/Analysis/SALSA
modality=atac

HOST_WORK_DIR=${project}/work_dir/${modality}/${sample}
mkdir -p ${HOST_WORK_DIR}

apptainer exec --cleanenv -B ${HOST_WORK_DIR}:/gatk_genotype/atac/${sample} ${sif_path} bash ${project}/step1_gatk_genotype.sh \
    --inputbam ${samplepath} \
    --library_id ${sample} \
    --reference ${cellranger_ref} \
    --gatk_bundle ${gatk_bundle_path} \
    --outputdir ${project}/atac_genotype \
    --outputvcf ${sample}_ATAC.vcf.gz \
    --modality atac \
    --verbose false
