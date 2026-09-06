# ============================================
# VEP Annotation Pipeline for vQTL Variants
# ============================================

conda activate vep

# (Optional) Install VEP cache if not already set up
# vep_install -a cf -s homo_sapiens -y GRCh38 -c /path/to/vep_cache --CONVERT

# Run VEP annotation
cd /path/to/project/vqtl_mapping
vep \
  --cache \
  --dir_cache /path/to/vep_cache \
  --cache_version 115 \
  --force_overwrite \
  --input_file input_variants.vcf \
  --offline \
  --output_file replicated_lead_vqtls_annotated.vcf \
  --no_stats \
  --vcf \
  --assembly GRCh38 \
  --fork 16 \
  --buffer_size 1000 \
  --mane \
  --canonical \
  --appris \
  --biotype \
  --protein \
  --pick \
  --pick_order mane_select,mane_plus_clinical,rank,canonical,appris,biotype,ensembl