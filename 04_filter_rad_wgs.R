# SNPfiltR filtering for the joint RAD+WGS dataset (DAPC / ancestry pipeline).
# Applies split max_depth thresholds: RAD samples capped at 40, WGS at 80,
# since pooled depth is diluted by RAD low coverage and lets WGS paralogs through.
# Output is used by 05_plink_prep_rad_wgs.sh for LD pruning and PLINK prep.
#
# Input:  data/vcf/ns_rad_prefiltered_no_M.vcf.gz  (autosomes + chrY, post-script 01)
# Output: data/vcf/ns_rad_filtered_v2.vcf.gz

library("SNPfiltR")
library("vcfR")

vcf <- read.vcfR("data/vcf/ns_rad_prefiltered_no_M.vcf.gz")

vcf <- filter_biallelic(vcf)
vcf <- hard_filter(vcf, depth = 5, gq = 20)
vcf <- missing_by_snp(vcf, cutoff = 0.85)
vcf <- missing_by_sample(vcf, cutoff = 0.9)
vcf <- filter_allele_balance(vcf)

# --- split max_depth by sequencing method --------------------------------
all_samples <- colnames(vcf@gt)[-1]
wgs_samples <- grep("_L002", all_samples, value = TRUE)
rad_samples <- setdiff(all_samples, wgs_samples)
message(sprintf("Split max_depth: %d RAD (cap 40), %d WGS (cap 80)",
                length(rad_samples), length(wgs_samples)))

site_id <- function(v) paste0(vcfR::getCHROM(v), "_", vcfR::getPOS(v))

vcf_rad <- vcf; vcf_rad@gt <- vcf@gt[, c("FORMAT", rad_samples), drop = FALSE]
vcf_wgs <- vcf; vcf_wgs@gt <- vcf@gt[, c("FORMAT", wgs_samples), drop = FALSE]

both_pass <- intersect(site_id(max_depth(vcf_rad, maxdepth = 40)),
                      site_id(max_depth(vcf_wgs, maxdepth = 80)))
vcf <- vcf[site_id(vcf) %in% both_pass, ]

# --- iterative cleanup ---------------------------------------------------
vcf <- min_mac(vcf, min.mac = 1)
vcf <- missing_by_sample(vcf, cutoff = 0.9)
vcf <- missing_by_snp(vcf, cutoff = 0.85)
vcf <- min_mac(vcf, min.mac = 1)

write.vcf(vcf, "data/vcf/ns_rad_filtered_v2.vcf.gz")
message("Wrote data/vcf/ns_rad_filtered_v2.vcf.gz")
message("NOTE: re-bgzip + tabix-index this file in WSL before running 05_plink_prep_rad_wgs.sh:")
message("  zcat data/vcf/ns_rad_filtered_v2.vcf.gz | bgzip -c > tmp.vcf.gz && \\")
message("  mv tmp.vcf.gz data/vcf/ns_rad_filtered_v2.vcf.gz && \\")
message("  tabix -f data/vcf/ns_rad_filtered_v2.vcf.gz")
