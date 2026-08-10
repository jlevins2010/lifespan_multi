#!/bin/bash
#BSUB -q normal
#BSUB -M 100GB
#BSUB -R "rusage[mem=100GB]"
#BSUB -J salsa
#BSUB -e ../errorFiles/salsa_preparePhase.err
#BSUB -o ../outputFiles/salsa_preparePhase.out
#BSUB -n 20

# Stop script immediately on error
set -e

sample=FK4
sif_path=/home/levinsj/MultiOme/Analysis/sif/salsa_latest.sif
vcfone=/home/levinsj/MultiOme/Analysis/SALSA/atac_genotype/FK4_ATAC.vcf.gz
vcftwo=/home/levinsj/MultiOme/Analysis/SALSA/rna_genotype/FK4_rna.vcf.gz
project=/home/levinsj/MultiOme/Analysis/SALSA
modality=merged

# directory setup
HOST_WORK_DIR=${project}/work_dir/temp
DOWNLOAD_DIR=${project}/biallelic_SNVs
mkdir -p ${HOST_WORK_DIR}
mkdir -p ${DOWNLOAD_DIR}

# Create the rename mapping file once, inside the work dir (safer than /tmp on clusters)
RENAME_FILE=${HOST_WORK_DIR}/rename_chrm.txt
for k in {1..22} X; do echo "${k} chr${k}"; done > ${RENAME_FILE}

# --- MAIN LOOP ---
for i in {1..22} X; do
    echo "Processing Chromosome ${i}..."
    
    # Define filename and URL variables
    inputvcf="ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
    url="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/${inputvcf}"
    
    # 1. Download
    wget -P ${DOWNLOAD_DIR} ${url}

    # 2. Annotate (Rename Chrs) - outputting to temp dir
    apptainer exec --cleanenv ${sif_path} bcftools annotate \
        ${DOWNLOAD_DIR}/${inputvcf} \
        --threads 4 \
        --rename-chrs ${RENAME_FILE} \
        -Oz -o ${HOST_WORK_DIR}/${inputvcf}

    # 3. Move processed file back to main folder
    mv ${HOST_WORK_DIR}/${inputvcf} ${DOWNLOAD_DIR}/

    # 4. Index
    apptainer exec --cleanenv ${sif_path} bcftools index -f --threads 4 --tbi ${DOWNLOAD_DIR}/${inputvcf}
done

echo "Completed successfully"
