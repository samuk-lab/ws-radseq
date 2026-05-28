# DAPC on the full RAD+WGS intersected dataset (ns_rad_pruned, 4,678 SNPs,
# n=224). Groups: species (cmn vs wht) from meta/popmap.txt.
# The faceted assignment plot shows
# within-population structure at the sympatric sites (CL, SR).
#
# Inputs:  data/admixture/ns_rad_pruned.bed/.bim/.fam  (05_plink_prep_rad_wgs.sh)
#          meta/popmap.txt
# Outputs: figures/dapc_rad_wgs_scores.pdf
#          figures/dapc_rad_wgs_assignment.pdf

if (!requireNamespace("adegenet", quietly = TRUE)) install.packages("adegenet")

library(adegenet)
library(tidyverse)

adm_dir  <- "data/admixture"
fig_dir  <- "figures"
raw_file <- file.path(adm_dir, "ns_rad_pruned.raw")

# ---- Create .raw file from PLINK bed (via plink in WSL) -----------------
if (!file.exists(raw_file)) {
  message("Creating .raw file via plink in WSL...")
  rc <- system(paste0(
    'wsl -e bash -c "',
    'source /home/ksamuk/miniconda3/etc/profile.d/conda.sh && ',
    'conda activate admixture_env && ',
    'cd /mnt/f/Dropbox/02_Projects/ns_radseq && ',
    'plink --bfile data/admixture/ns_rad_pruned ',
    '--recode A --out data/admixture/ns_rad_pruned 2>&1"'
  ))
  if (rc != 0) stop("plink --recode A failed")
}

# ---- Read genotypes (genlight object) -----------------------------------
message("Reading genotypes...")
geno <- read.PLINK(raw_file)

# ---- Attach species + population labels from popmap --------------------
popmap <- read.table("meta/popmap.txt", header = TRUE, stringsAsFactors = FALSE)
sample_ids <- indNames(geno)
meta <- popmap[match(sample_ids, popmap$sample_id), ]

if (any(is.na(meta$species)))
  warning(sprintf("%d samples not found in popmap", sum(is.na(meta$species))))

pop(geno) <- factor(meta$species, levels = c("cmn", "wht"))
message(sprintf("Samples: %d  (%d cmn / %d wht)",
                nInd(geno), sum(meta$species == "cmn"), sum(meta$species == "wht")))

# ---- Cross-validate number of PCs (xvalDapc) ----------------------------
set.seed(42)
message("Cross-validating number of PCs (100 reps, n.pca.max = 50) ...")
X <- as.matrix(geno)
message(sprintf("Imputing %d missing genotypes (%.2f%%) with column means ...",
                sum(is.na(X)), 100 * mean(is.na(X))))
col_means <- colMeans(X, na.rm = TRUE)
na_idx    <- which(is.na(X), arr.ind = TRUE)
X[na_idx] <- col_means[na_idx[, 2]]

xval <- xvalDapc(X, pop(geno),
                 n.pca.max     = 50,
                 training.set  = 0.9,
                 result        = "groupMean",
                 center        = TRUE,
                 scale         = FALSE,
                 n.rep         = 100,
                 xval.plot     = FALSE)
n_pca_opt <- as.integer(xval$`Number of PCs Achieving Lowest MSE`)
message(sprintf("Optimal n.pca: %d", n_pca_opt))

# ---- Run DAPC -----------------------------------------------------------
dapc_res <- dapc(geno, pop = pop(geno), n.pca = n_pca_opt, n.da = 1)

assign_acc <- mean(as.character(dapc_res$assign) == as.character(pop(geno)))
message(sprintf("Assignment accuracy (training): %.1f%%", 100 * assign_acc))

# ---- Build score data frame with population labels ----------------------
pop_levels <- c("BC", "CB", "CL_RAD", "CL_WGS", "MJ", "SR", "MJ_or_")

score_df <- data.frame(
  sample_id  = sample_ids,
  species    = factor(as.character(pop(geno)), levels = c("cmn", "wht")),
  population = meta$population,
  p_wht      = as.numeric(dapc_res$posterior[, "wht"]),
  method     = ifelse(grepl("_L002", sample_ids), "WGS", "RAD"),
  stringsAsFactors = FALSE
) %>%
  mutate(population = ifelse(population == "CL",
                             paste0("CL_", method),
                             population),
         population = factor(population, levels = pop_levels))

# ---- P(white) density + rug plot ----------------------------------------
p_density <- ggplot(score_df, aes(x = p_wht, fill = species, colour = species)) +
  geom_density(alpha = 0.45, linewidth = 0.8) +
  geom_rug(alpha = 0.5, linewidth = 0.5) +
  scale_fill_brewer(palette   = "Set1",
                    labels    = c("Common", "White"),
                    direction = -1) +
  scale_colour_brewer(palette   = "Set1",
                      labels    = c("Common", "White"),
                      direction = -1) +
  labs(x      = "Posterior P(white form)",
       y      = "Density",
       fill   = NULL,
       colour = NULL,
       title  = sprintf("DAPC RAD+WGS  (n.pca = %d, accuracy = %.0f%%)",
                        n_pca_opt, 100 * assign_acc)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

ggsave(file.path(fig_dir, "dapc_rad_wgs_scores.pdf"),
       p_density, width = 7, height = 4)

# ---- Posterior assignment bar plot, faceted by population ---------------
post_df <- as.data.frame(dapc_res$posterior) %>%
  rownames_to_column("sample_id") %>%
  left_join(score_df[, c("sample_id", "species", "population", "p_wht")],
            by = "sample_id") %>%
  arrange(population, p_wht) %>%
  group_by(population) %>%
  mutate(within_pos = row_number()) %>%
  ungroup() %>%
  pivot_longer(cols = c("cmn", "wht"),
               names_to = "assigned_to", values_to = "prob")

p_assign <- ggplot(post_df, aes(x = within_pos, y = prob, fill = assigned_to)) +
  geom_col(width = 1) +
  facet_grid(. ~ population, scales = "free_x", space = "free_x",
             switch = "x") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette   = "Set1",
                    labels    = c("Common", "White"),
                    direction = -1) +
  labs(x    = NULL,
       y    = "Posterior probability",
       fill = "Assigned to") +
  theme_bw(base_size = 12) +
  theme(axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        panel.spacing.x  = unit(0.15, "lines"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        legend.position  = "top")

ggsave(file.path(fig_dir, "dapc_rad_wgs_assignment.pdf"),
       p_assign, width = 12, height = 4)

message("Figures written:")
message("  figures/dapc_rad_wgs_scores.pdf")
message("  figures/dapc_rad_wgs_assignment.pdf")
