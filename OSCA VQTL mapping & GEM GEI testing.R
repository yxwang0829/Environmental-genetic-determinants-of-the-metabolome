#!/bin/bash
################################################################################
# OSCA VQTL Mapping & GEM GEI Testing Pipeline
# Description: Variance quantitative trait locus mapping and GxE interaction testing
################################################################################

#==============================================================================
# PART 1: OSCA VQTL Mapping
#==============================================================================

#------------------------------------------------------------------------------
# Step 1: Run VQTL mapping per chromosome
#------------------------------------------------------------------------------

pheno=$1
chr=$2
set="England"  # Options: England, Scotland_Wales
cd /path/to/your/project

# Run OSCA VQTL
./osca \
  --vqtl \
  --vqtl-mtd 2 \
  --bfile ./data/genotype_chr${chr} \
  --pheno ./data/phenotypes_${pheno}.txt \
  --keep ./data/sample_ids_${set}.txt \
  --thread-num 8 \
  --out ./results/vqtl/${set}_${pheno}_chr${chr}
# Note: Only variants with MAF >= 1% are included
# Only England samples used initially; Scotland_Wales tested later for significant hits

#------------------------------------------------------------------------------
# Step 2: Clump VQTL results
#------------------------------------------------------------------------------
cd ./results/clump

input_file="../vqtl/significant_vqtl_results.tsv"
output_dir="./clump_results"
mkdir -p "${output_dir}"

# Extract unique metabolites
metabolites=$(tail -n +2 "${input_file}" | cut -f1 | sort | uniq)

for metabolite in ${metabolites}; do
  echo "Processing metabolite: ${metabolite}"
  
  # Create temporary input file with required columns (SNP, CHR, BP, P)
  temp_input="${metabolite}_temp.tsv"
  head -n 1 "${input_file}" > "${temp_input}"
  grep "^${metabolite}[[:space:]]" "${input_file}" >> "${temp_input}"
  awk -F'\t' '{print $3, $4, $5, $13}' OFS='\t' "${temp_input}" > temp_file && mv temp_file "${temp_input}"
  
  # Run clumping per chromosome
  for chr in {1..22}; do
    echo "  Chromosome ${chr}"
    
    plink2 \
      --bfile ./data/genotype_chr${chr} \
      --rm-dup force-first \
      --keep ./data/clump_keep_samples.txt \
      --clump "${temp_input}" \
      --clump-p1 1.6e-10 \
      --clump-p2 1.6e-10 \
      --clump-r2 0.1 \
      --clump-kb 1000 \
      --clump-field P \
      --clump-snp-field SNP \
      --chr ${chr} \
      --out "${output_dir}/clumping_${metabolite}_chr${chr}"
  done
  
  rm "${temp_input}"
done

echo "Clumping completed"

#------------------------------------------------------------------------------
# Step 3: Extract lead SNPs from clumping results
#------------------------------------------------------------------------------
# cd ./clump_results

input_file="../vqtl/significant_vqtl_results.tsv"
metabolites=$(tail -n +2 "${input_file}" | cut -f1 | sort | uniq)

for metabolite in ${metabolites}; do
  temp_file="temp_${metabolite}.txt"
  
  # Merge clumping results from all chromosomes
  for file in clumping_${metabolite}_chr*.clumps; do
    if [[ -f "$file" ]]; then
      if [[ ! -s "$temp_file" ]]; then
        cat "$file" > "$temp_file"  # First file: include header
      else
        tail -n +2 "$file" >> "$temp_file"  # Subsequent files: append without header
      fi
    fi
  done
  
  mv "$temp_file" "${metabolite}_independent_vqtls.txt"
done


#==============================================================================
# PART 2: GEM GxE Interaction Testing
#==============================================================================

#------------------------------------------------------------------------------
# Generate job scripts for GEM GEI testing
#------------------------------------------------------------------------------
# Reference: https://github.com/large-scale-gxe-methods/GEM

cd /path/to/project

# Extract phenotype list from phenotype file
pheno_file="./data/phenotypes.txt"
pheno_list=$(awk 'NR==1 {for(i=3; i<=NF; i++) print $i}' ${pheno_file})

# Template job script
template_file="./scripts/GEM_GEI_jobs.sh"
> "${template_file}"

first_line=true

# Read exposure list (format: exposure_name, exposure_type, ...)
while read line; do
  # Skip header
  if [ "$first_line" = true ]; then
    first_line=false
    continue
  fi
  
  interacting_exposure=$(echo "$line" | awk -F'\t' '{print $7}')
  interacting_exposure_type=$(echo "$line" | awk -F'\t' '{print $6}')
  
  for chr in {1..22}; do
    job_name="${interacting_exposure}_chr${chr}"
    job_script="./scripts/GEI_testing/${job_name}_GEM_GEI.sh"
    log_file="./logs/${job_name}_GEM_GEI.log"
    
    # Build GEM command
    cmd="GEM_2.2.1 \\
  --bfile ./data/genotype_chr${chr} \\
  --include-snp-file ./data/vqtl_snps_chr${chr}.txt \\
  --pheno-file ./data/GEI_phenotypes.txt \\
  --delim \\t \\
  --sampleid-name IID \\
  --pheno-name \${pheno} \\
  --exposure-names ${interacting_exposure} \\
  --covar-names Sex Age Fasting Month Center Provider Batch PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 PC17 PC18 PC19 PC20 PC21 PC22 PC23 PC24 PC25 PC26 PC27 PC28 PC29 PC30 PC31 PC32 PC33 PC34 PC35 PC36 PC37 PC38 PC39 PC40 \\
  --robust 1 \\
  --threads 7 \\
  --maf 0.0001 \\
  --output-style meta \\
  --out ./GEM_output/\${pheno}_chr${chr}_${interacting_exposure}_GEM"
    
    # Add categorical covariate flag if exposure is categorical
    if [[ "$interacting_exposure_type" == "Binary" || "$interacting_exposure_type" == "CatOrd" || "$interacting_exposure_type" == "CatUnord" ]]; then
      cmd="${cmd} --categorical-names ${interacting_exposure}"
    fi
    
    # Create job script
    cat > "${job_script}" <<EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --output=${log_file}
#SBATCH --error=${log_file}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10GB

cd /path/to/project/GEI_testing/

# Loop through all phenotypes
pheno_list=\$(awk 'NR==1 {for(i=3; i<=NF; i++) print \$i}' ${pheno_file})

for pheno in \${pheno_list}; do
  output_file="./GEM_output/\${pheno}_chr${chr}_${interacting_exposure}_GEM.out"
  
  # Skip if output already exists
  if [[ -f "\${output_file}" ]]; then
    continue
  fi
  
  ${cmd}
done
EOF
    
    # Add to job submission template
    echo "sbatch ${job_script}" >> "${template_file}"
    
  done
done < ./data/GEI_exposure_map.txt

echo "GEM GEI job scripts generated: ${template_file}"
