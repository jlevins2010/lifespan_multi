#!/bin/bash
#BSUB -q normal
#BSUB -M 100GB
#BSUB -R "rusage[mem=100GB]"
#BSUB -J salsa_step3
#BSUB -e salsa_step3.err
#BSUB -o salsa_step3.out
#BSUB -n 20

sample=FK4
samplepath=/home/levinsj/MultiOme/Fetal/cellranger-arc/FK4/outs/atac_possorted_bam.bam
sif_path=/home/levinsj/MultiOme/Analysis/sif/salsa_latest.sif
cellranger_ref=/home/levinsj/Applications/refdata-cellranger-arc-GRCh38-2020-A-2.0.0
gatk_bundle_path=/home/levinsj/Applications/gatk
reference=/home/levinsj/Applications
project=/home/levinsj/MultiOme/Analysis/SALSA

HOST_WORK_DIR=${project}/phasing/
mkdir -p ${HOST_WORK_DIR}/${sample}

apptainer exec --cleanenv --writable-tmpfs -B ${HOST_WORK_DIR}:/phasing ${sif_path} bash ${project}/step3_phase_vcf.sh \
		--library_id ${sample} \
		--inputvcf /home/levinsj/MultiOme/Analysis/SALSA/merged_genotype/${sample}.vcf.gz \
		--outputdir ${project}/phasing/ \
		--outputvcf ${sample}.pass.merged.hcphase.vcf.gz \
		--phasingref ${project}/biallelic_SNVs \
		--hcphase \
		--snvonly 
