# SNPfiltR run on the mainland-only subset (Cape Breton dropped FIRST).
# Per-SNP missingness and MAF are recomputed on the 129 mainland samples only,
# rescuing loci that failed in the full set due to BC/MJ-specific missingness.
#
# Input:  data/vcf/ns_rad_prefiltered_no_M.vcf.gz
# Output: data/vcf/ns_rad_filtered_mainland.vcf.gz
#
# After running this script, re-bgzip + index the output in WSL:
#   conda activate admixture_env
#   zcat data/vcf/ns_rad_filtered_mainland.vcf.gz | bgzip -c > tmp.vcf.gz && \
#   mv tmp.vcf.gz data/vcf/ns_rad_filtered_mainland.vcf.gz && \
#   tabix -f data/vcf/ns_rad_filtered_mainland.vcf.gz

library("SNPfiltR")
library("vcfR")

vcf <- read.vcfR("data/vcf/ns_rad_prefiltered_no_M.vcf.gz")

# Drop Cape Breton samples before any per-SNP filtering
all_samples  <- colnames(vcf@gt)[-1]
keep_samples <- grep("^BC[0-9]|^MJ[0-9]|^MJ2_or", all_samples,
                     value = TRUE, invert = TRUE)
message(sprintf("Subset: %d -> %d mainland samples",
                length(all_samples), length(keep_samples)))
vcf@gt <- vcf@gt[, c("FORMAT", keep_samples), drop = FALSE]

vcf <- filter_biallelic(vcf)
vcf <- hard_filter(vcf, depth = 5, gq = 20)
vcf <- missing_by_snp(vcf, cutoff = 0.85)
vcf <- missing_by_sample(vcf, cutoff = 0.9)
vcf <- filter_allele_balance(vcf)

# Split max_depth by sequencing method on the mainland subset
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

vcf <- min_mac(vcf, min.mac = 1)
vcf <- missing_by_sample(vcf, cutoff = 0.9)
vcf <- missing_by_snp(vcf, cutoff = 0.85)
vcf <- min_mac(vcf, min.mac = 1)

write.vcf(vcf, "data/vcf/ns_rad_filtered_mainland.vcf.gz")
message("Wrote data/vcf/ns_rad_filtered_mainland.vcf.gz (re-bgzip in WSL)")
