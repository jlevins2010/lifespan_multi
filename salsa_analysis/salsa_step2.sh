#!/bin/bash
#BSUB -q normal
#BSUB -M 100GB
#BSUB -R "rusage[mem=100GB]"
#BSUB -J salsa
#BSUB -e ../errorFiles/salsa_step2.err
#BSUB -o ../outputFiles/salsa_step2.out
#BSUB -n 20

sample=FK4
sif_path=/home/levinsj/MultiOme/Analysis/sif/salsa_latest.sif
vcfone=/home/levinsj/MultiOme/Analysis/SALSA/atac_genotype/FK4_ATAC.vcf.gz
vcftwo=/home/levinsj/MultiOme/Analysis/SALSA/rna_genotype/FK4_rna.vcf.gz
project=/home/levinsj/MultiOme/Analysis/SALSA
modality=merged

HOST_WORK_DIR=${project}/work_dir/${modality}/${sample}
mkdir -p ${HOST_WORK_DIR}

if [ -f ${vcfone} ] && [ -f ${vcftwo} ]; then
	echo "Both files exist."
else
	echo "One or both files are missing."
	exit 1
fi

apptainer exec --cleanenv -B ${HOST_WORK_DIR}:/gatk_genotype/${modality}/${sample} ${sif_path} bash ${project}/step2_merge_geno.sh \
    --vcfone ${vcfone} \
	--vcftwo ${vcftwo} \
    --library_id ${sample} \
    --outputdir ${project}/${modality}/ \
    --outputvcf ${sample}.vcf.gz 
