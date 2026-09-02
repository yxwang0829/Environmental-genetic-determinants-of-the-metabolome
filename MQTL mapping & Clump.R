#!/bin/bash
################################################################################
# mQTL Mapping Pipeline - REGENIE Step 1
# Step 1: Generate genome-wide predictions
################################################################################
cd /path/to/your/project

# Input: Population subset (e.g., England, Scotland_Wales)
set=$1

source /path/to/miniconda/etc/profile.d/conda.sh
conda activate regenie4.1

# Run REGENIE step 1
regenie \
  --step 1 \
  --bed /path/to/genotype_data/ukb_cal_allChrs_hg38 \
  --keep ./data/${set}_caucasian_metabolome_sample_ids.txt \
  --extract /path/to/qc_pass.snplist \
  --phenoFile ./data/caucasian_regenie_phenotypes.txt \
  --covarFile ./data/caucasian_regenie_covariates.txt \
  --catCovarList Sex,Ethnicity,Month,Center,Provider,Batch \
  --maxCatLevels 26 \
  --apply-rint \
  --force-qt \
  --bsize 1000 \
  --threads 32 \
  --lowmem \
  --lowmem-prefix ./results/step1/step1_${set}_tmp_preds \
  --out ./results/step1/step1_${set}

################################################################################
# Step 2: Single variant association testing (per chromosome)
################################################################################
chr=$1
set=$2
cd /path/to/your/project
source /path/to/miniconda/etc/profile.d/conda.sh
conda activate regenie4.1

# Run REGENIE step 2
regenie \
  --step 2 \
  --keep ./data/${set}_caucasian_metabolome_sample_ids.txt \
  --extract ./data/common_low_freq_variants.ids \
  --bed /path/to/WGS_data/caucasian_chr${chr} \
  --phenoFile ./data/caucasian_regenie_phenotypes.txt \
  --covarFile ./data/caucasian_regenie_covariates.txt \
  --catCovarList Sex,Ethnicity,Month,Center,Provider,Batch \
  --maxCatLevels 26 \
  --pred ./results/step1/step1_${set}_pred.list \
  --apply-rint \
  --force-qt \
  --bsize 1000 \
  --out ./results/step2_single/${set}_chr${chr}_single_variant_step2
# Note: Only variants with MAF >= 0.01% are included by default

################################################################################
# Clumping: LD-based variant pruning
################################################################################

cd /path/to/your/project/mqtl_mapping/clump

input_file="../results/significant_mqtl_results.tsv"
output_dir="./clump_results"
mkdir -p "${output_dir}"

# Extract unique metabolites from input file
metabolites=$(tail -n +2 "${input_file}" | cut -f1 | sort | uniq)

for metabolite in ${metabolites}; do
  echo "Processing metabolite: ${metabolite}"
  
  # Create temporary file for this metabolite
  temp_input="${metabolite}_temp.tsv"
  head -n 1 "${input_file}" > "${temp_input}"
  grep "^${metabolite}[[:space:]]" "${input_file}" >> "${temp_input}"
  
  # Run clumping per chromosome
  for chr in {1..22}; do
    echo "  Processing chromosome ${chr}"
    
    /path/to/plink2 \
      --bfile /path/to/WGS_data/caucasian_chr${chr} \
      --rm-dup force-first \
      --keep ./data/clump_keep_samples.txt \
      --clump "${temp_input}" \
      --clump-p1 1.6e-10 \
      --clump-p2 1.6e-10 \
      --clump-r2 0.1 \
      --clump-kb 1000 \
      --clump-field P \
      --clump-snp-field ID \
      --chr ${chr} \
      --out "${output_dir}/clumping_${metabolite}_chr${chr}"
  done
  
  rm "${temp_input}"
done

echo "Clumping completed for all metabolites and chromosomes."