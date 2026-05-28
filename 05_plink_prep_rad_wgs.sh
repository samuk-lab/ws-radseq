#!/usr/bin/env bash
# VCF -> PLINK -> MAF/LD filter -> PLINK binary set for DAPC (09_dapc_rad_wgs.R).
# Produces ns_rad_pruned.{bed,bim,fam,ped} used by the DAPC and mainland sNMF pipelines.
#
# Run from project root in WSL with `admixture_env` activated:
#   conda activate admixture_env
#   bash 05_plink_prep_rad_wgs.sh

set -euo pipefail

IN_VCF="data/vcf/ns_rad_filtered_v2.vcf.gz"
OUT_DIR="data/admixture"
AUT_VCF="${OUT_DIR}/ns_rad_autosomes.vcf.gz"
PLINK_BASE="${OUT_DIR}/ns_rad"
PRUNED_BASE="${OUT_DIR}/ns_rad_pruned"

mkdir -p "${OUT_DIR}"

# Autosomes only (drop chrY, chrXIX, chrM, chrUn)
AUTOSOMES="chrI,chrII,chrIII,chrIV,chrV,chrVI,chrVII,chrVIII,chrIX,chrX,chrXI,chrXII,chrXIII,chrXIV,chrXV,chrXVI,chrXVII,chrXVIII,chrXX,chrXXI"
bcftools view -r "${AUTOSOMES}" -Oz -o "${AUT_VCF}" "${IN_VCF}"
tabix -f -p vcf "${AUT_VCF}"

# VCF -> PLINK with MAF >= 0.05 and unique chr:pos variant IDs
plink --vcf "${AUT_VCF}" \
      --double-id --allow-extra-chr \
      --set-missing-var-ids @:# \
      --maf 0.05 \
      --make-bed --out "${PLINK_BASE}"

# Rename chr to integers (PLINK auto-maps chrX -> 23; remap to 10)
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

# LD prune and write the pruned bed + ped/map for LEA
plink --bfile "${PLINK_BASE}" --indep-pairwise 50 10 0.2 --out "${PRUNED_BASE}"
plink --bfile "${PLINK_BASE}" --extract "${PRUNED_BASE}.prune.in" \
      --make-bed --out "${PRUNED_BASE}"
plink --bfile "${PRUNED_BASE}" --recode --out "${PRUNED_BASE}"

N_SAMPLES=$(wc -l < "${PRUNED_BASE}.fam")
N_SNPS=$(wc -l < "${PRUNED_BASE}.bim")
echo ""
echo "----- snmf prep done -----"
echo "Samples: ${N_SAMPLES}"
echo "Pruned SNPs (MAF >= 0.05, LD r2 < 0.2): ${N_SNPS}"
