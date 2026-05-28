# Extract values needed to fill in paragraph placeholders:
#   1. Mainland sNMF: best CE for K=1 and K=2 -> delta CE
#   2. Mainland sNMF K=2: Q values for the 6 hybrid candidates
#   3. RAD+WGS DAPC: posterior P(wht) for the 6 hybrid candidates
#      (uses known n.pca=30, skips xvalDapc)

suppressPackageStartupMessages({
  library(LEA)
  library(adegenet)
  library(tidyverse)
})

hybrid_ids <- c("CL30", "CL34", "FCL113_S14_L002",
                "sal_riv27", "sal_riv41", "sal_riv42")

# =========================================================================
# 1 + 2. Mainland sNMF
# =========================================================================
adm_dir   <- "data/admixture"
base_main <- file.path(adm_dir, "ns_rad_pruned_mainland")
ped_main  <- paste0(base_main, ".ped")
geno_main <- paste0(base_main, ".geno")
proj_main <- paste0(base_main, ".snmfProject")

# Run sNMF fresh with 30 reps (enough for stable CE; bypasses Windows
# snmfClass save bug that caused 13_snmf_mainland.R to crash at 100 reps)
if (file.exists(proj_main)) {
  message("Loading existing mainland snmfProject ...")
  obj <- tryCatch(
    { o <- load.snmfProject(proj_main)
      # verify at least K=1 has runs
      cross.entropy(o, K = 1)
      o },
    error = function(e) NULL
  )
} else {
  obj <- NULL
}

if (is.null(obj)) {
  message("Running mainland sNMF (K=1-6, 30 reps) ...")
  if (file.exists(proj_main)) remove.snmfProject(proj_main)
  obj <- snmf(geno_main, K = 1:6, entropy = TRUE, repetitions = 30,
              project = "new", CPU = 4, seed = 42)
}

# Cross-entropy per K
Ks <- 1:6
ce_df <- map_dfr(Ks, function(K) {
  ce <- cross.entropy(obj, K = K)
  tibble(K = K, CE = as.numeric(ce))
})
best_ce <- ce_df %>% group_by(K) %>% summarise(best_CE = min(CE))
ce_k1 <- best_ce$best_CE[best_ce$K == 1]
ce_k2 <- best_ce$best_CE[best_ce$K == 2]
delta_ce <- ce_k2 - ce_k1

cat("\n--- Mainland sNMF cross-entropy ---\n")
print(as.data.frame(best_ce))
cat(sprintf("\nBest K=1 CE: %.6f\nBest K=2 CE: %.6f\nDelta CE (K2-K1): %.6f\n",
            ce_k1, ce_k2, delta_ce))

# K=2 Q values for hybrid candidates
best_run_k2 <- ce_df %>% filter(K == 2) %>% slice_min(CE, n = 1) %>% pull(K)
# get the run index with min CE at K=2
ce_k2_all <- cross.entropy(obj, K = 2)
best_rep   <- which.min(ce_k2_all)

q_mat <- Q(obj, K = 2, run = best_rep)
sample_ids_main <- read.table(ped_main)[[2]]

q_df <- data.frame(sample_id = sample_ids_main,
                   anc1 = q_mat[, 1],
                   anc2 = q_mat[, 2])

popmap <- read.table("meta/popmap.txt", header = TRUE, stringsAsFactors = FALSE)
q_df <- left_join(q_df, popmap[, c("sample_id", "species", "population")],
                  by = "sample_id")

# Identify which ancestry component corresponds to "wht"
# (component with higher mean Q in wht samples)
wht_mean_anc1 <- mean(q_df$anc1[q_df$species == "wht"], na.rm = TRUE)
wht_mean_anc2 <- mean(q_df$anc2[q_df$species == "wht"], na.rm = TRUE)
wht_anc <- if (wht_mean_anc1 > wht_mean_anc2) "anc1" else "anc2"
cat(sprintf("\nWhite-form ancestry component: %s\n", wht_anc))

q_df$q_wht <- if (wht_anc == "anc1") q_df$anc1 else q_df$anc2

cat("\n--- K=2 Q(white) for hybrid candidates ---\n")
hybrids_snmf <- q_df %>%
  filter(sample_id %in% hybrid_ids) %>%
  select(sample_id, population, species, q_wht) %>%
  arrange(q_wht)
print(hybrids_snmf)

cat(sprintf("\nQ(white) range for hybrid candidates: %.3f - %.3f\n",
            min(hybrids_snmf$q_wht), max(hybrids_snmf$q_wht)))

# For context, show full range of Q(white) for pure-form samples
pure_range <- q_df %>%
  filter(!sample_id %in% hybrid_ids) %>%
  group_by(species) %>%
  summarise(min_q_wht = min(q_wht), max_q_wht = max(q_wht))
cat("\n--- Q(white) range for non-hybrid samples ---\n")
print(as.data.frame(pure_range))

# =========================================================================
# 3. RAD+WGS DAPC posteriors (n.pca=30 already cross-validated)
# =========================================================================
cat("\n\n--- RAD+WGS DAPC posteriors ---\n")
raw_file <- file.path(adm_dir, "ns_rad_pruned.raw")
if (!file.exists(raw_file)) stop(".raw file not found — run 18_dapc_rad_wgs.R first")

geno <- read.PLINK(raw_file)

sample_ids_full <- indNames(geno)
meta <- popmap[match(sample_ids_full, popmap$sample_id), ]
pop(geno) <- factor(meta$species, levels = c("cmn", "wht"))

# Impute NAs
X <- as.matrix(geno)
col_means <- colMeans(X, na.rm = TRUE)
na_idx    <- which(is.na(X), arr.ind = TRUE)
X[na_idx] <- col_means[na_idx[, 2]]

# Run DAPC with known optimal n.pca (skip xvalDapc)
set.seed(42)
dapc_res <- dapc(geno, pop = pop(geno), n.pca = 30, n.da = 1)

post_df <- as.data.frame(dapc_res$posterior) %>%
  rownames_to_column("sample_id")

cat("\n--- Posterior P(white) for hybrid candidates ---\n")
hybrids_dapc <- post_df %>%
  filter(sample_id %in% hybrid_ids) %>%
  select(sample_id, p_wht = wht) %>%
  left_join(popmap[, c("sample_id", "species", "population")], by = "sample_id") %>%
  arrange(p_wht)
print(hybrids_dapc)

cat(sprintf("\nP(white) range for hybrid candidates: %.3f - %.3f\n",
            min(hybrids_dapc$p_wht), max(hybrids_dapc$p_wht)))

# Context: pure-form range
pure_dapc <- post_df %>%
  filter(!sample_id %in% hybrid_ids) %>%
  left_join(popmap[, c("sample_id", "species")], by = "sample_id") %>%
  group_by(species) %>%
  summarise(min_p_wht = min(wht), max_p_wht = max(wht))
cat("\n--- P(white) range for non-hybrid samples ---\n")
print(as.data.frame(pure_dapc))

# =========================================================================
# Summary for paragraph
# =========================================================================
cat("\n\n========= PARAGRAPH VALUES =========\n")
cat(sprintf("Delta CE (mainland sNMF K2 - K1): %.4f\n", delta_ce))
cat(sprintf("Hybrid Q(white) range (sNMF K=2): %.2f - %.2f\n",
            min(hybrids_snmf$q_wht), max(hybrids_snmf$q_wht)))
cat(sprintf("Hybrid P(white) range (DAPC):     %.2f - %.2f\n",
            min(hybrids_dapc$p_wht), max(hybrids_dapc$p_wht)))
n_mainland <- nrow(read.table(ped_main))
cat(sprintf("Mainland n:                        %d\n", n_mainland))
cat(sprintf("Hybrid fraction:                   %d / %d (%.1f%%)\n",
            length(hybrid_ids), n_mainland,
            100 * length(hybrid_ids) / n_mainland))
