#!/usr/bin/env bash
# Mainland-only sensitivity prep, power-boost variant.
# Input VCF is already mainland-only (filtered by 06_filter_mainland.R, which
# also re-runs SNPfiltR with per-SNP missingness recomputed on mainland samples).
# This script: autosome subset -> PLINK -> MAF 0.03, no LD pruning.
#
# Run from project root in WSL with `admixture_env` activated:
#   conda activate admixture_env
#   bash 07_snmf_prep_mainland.sh

set -euo pipefail

IN_VCF="data/vcf/ns_rad_filtered_mainland.vcf.gz"
OUT_DIR="data/admixture"
AUT_VCF="${OUT_DIR}/ns_rad_autosomes_mainland.vcf.gz"
PLINK_BASE="${OUT_DIR}/ns_rad_mainland"
PRUNED_BASE="${OUT_DIR}/ns_rad_pruned_mainland"

mkdir -p "${OUT_DIR}"

AUTOSOMES="chrI,chrII,chrIII,chrIV,chrV,chrVI,chrVII,chrVIII,chrIX,chrX,chrXI,chrXII,chrXIII,chrXIV,chrXV,chrXVI,chrXVII,chrXVIII,chrXX,chrXXI"
bcftools view -r "${AUTOSOMES}" -Oz -o "${AUT_VCF}" "${IN_VCF}"
tabix -f -p vcf "${AUT_VCF}"

plink --vcf "${AUT_VCF}" \
      --double-id --allow-extra-chr \
      --set-missing-var-ids @:# \
      --maf 0.03 \
      --make-bed --out "${PLINK_BASE}"

awk 'BEGIN{
        OFS="\t";
        m["chrI"]=1;  m["chrII"]=2;  m["chrIII"]=3;  m["chrIV"]=4;  m["chrV"]=5;
        m["chrVI"]=6; m["chrVII"]=7; m["chrVIII"]=8; m["chrIX"]=9;  m["chrX"]=10;
        m["chrXI"]=11; m["chrXII"]=12; m["chrXIII"]=13; m["chrXIV"]=14; m["chrXV"]=15;
        m["chrXVI"]=16; m["chrXVII"]=17; m["chrXVIII"]=18; m["chrXX"]=20; m["chrXXI"]=21;
        m["23"]=10;
     }
     { if ($1 in m) { $1 = m[$1] } else { print "Unmapped: "$1 > "/dev/stderr"; exit 1 } print }' \
     "${PLINK_BASE}.bim" > "${PLINK_BASE}.bim.tmp"
mv "${PLINK_BASE}.bim.tmp" "${PLINK_BASE}.bim"

# No LD pruning — sNMF tolerates correlated SNPs. Just copy bed -> pruned name
# and write ped/map (kept the "pruned" name to avoid renaming all downstream).
plink --bfile "${PLINK_BASE}" --make-bed --out "${PRUNED_BASE}"
plink --bfile "${PRUNED_BASE}" --recode --out "${PRUNED_BASE}"

N_SAMPLES=$(wc -l < "${PRUNED_BASE}.fam")
N_SNPS=$(wc -l < "${PRUNED_BASE}.bim")
echo ""
echo "----- mainland prep done -----"
echo "Samples: ${N_SAMPLES}"
echo "SNPs (MAF >= 0.03, no LD prune): ${N_SNPS}"
