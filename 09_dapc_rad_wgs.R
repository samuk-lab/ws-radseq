# DAPC on the RAD+WGS intersected dataset (ns_rad_pruned), excluding the
# allopatric BC and MJ populations. Training is on the 40 WGS samples
# (20 CB/cmn + 20 CL/wht), which carry a priori morphology labels.
# RAD samples from CL, SR (and any remaining) are held out and assigned via
# predict.dapc(), avoiding the circularity of training on labels that were
# themselves derived by PCA on the same genetic data.
#
# SNPs are further restricted to those with <= 5% missingness in the RAD
# subset to avoid column-mean imputation biasing held-out samples toward the
# midpoint of LD1.
#
# Inputs:  data/admixture/ns_rad_pruned.bed/.bim/.fam  (05_plink_prep_rad_wgs.sh)
#          meta/popmap.txt
# Outputs: figures/dapc_rad_wgs_scores.pdf
#          figures/dapc_rad_wgs_assignment.pdf
#          figures/dapc_rad_wgs_scatter.pdf
#          figures/dapc_rad_wgs_scatter_radonly.pdf
#          figures/dapc_rad_wgs_combined.pdf
#          figures/dapc_rad_wgs_combined_all.pdf

if (!requireNamespace("adegenet", quietly = TRUE)) install.packages("adegenet")

library(adegenet)
library(tidyverse)
library(patchwork)

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

# ---- Attach metadata; identify WGS (a priori labels) vs RAD ------------
popmap <- read.table("meta/popmap.txt", header = TRUE, stringsAsFactors = FALSE)
sample_ids_all <- indNames(geno)
meta_all <- popmap[match(sample_ids_all, popmap$sample_id), ]

if (any(is.na(meta_all$species)))
  warning(sprintf("%d samples not found in popmap",
                  sum(is.na(meta_all$species))))

# Drop allopatric populations (BC, MJ; MJ_or_ rolled in with MJ).
excluded_pops <- c("BC", "MJ", "MJ_or_")
keep <- !meta_all$population %in% excluded_pops
sample_ids <- sample_ids_all[keep]
meta       <- meta_all[keep, ]

# WGS samples are the ones with a priori (morphology) species labels.
is_wgs <- grepl("_L002", sample_ids)

message(sprintf("Excluded populations: %s  (dropped %d samples)",
                paste(excluded_pops, collapse = ", "), sum(!keep)))
message(sprintf("Kept: %d total  |  WGS (train) %d  |  RAD (predict) %d",
                length(sample_ids), sum(is_wgs), sum(!is_wgs)))
message(sprintf("  WGS labels: %d cmn / %d wht",
                sum(is_wgs & meta$species == "cmn"),
                sum(is_wgs & meta$species == "wht")))

# ---- Drop SNPs with high missingness in the RAD subset -----------------
# RAD samples (~7.5% missing/sample) are the bottleneck; WGS averages 3%.
# Removing SNPs with >max_rad_miss missing in RAD ensures held-out samples
# aren't dragged toward the column-mean midpoint at sites they don't cover.
max_rad_miss <- 0.05
X <- as.matrix(geno)[keep, , drop = FALSE]
miss_rad   <- colMeans(is.na(X[!is_wgs, , drop = FALSE]))
snps_keep  <- which(miss_rad <= max_rad_miss)
message(sprintf("SNP filter: %d/%d SNPs retained (RAD missing <= %.0f%%)",
                length(snps_keep), ncol(X), 100 * max_rad_miss))
X <- X[, snps_keep, drop = FALSE]

# ---- Impute remaining missing genotypes with column means --------------
message(sprintf("Imputing %d missing genotypes (%.2f%%) with column means ...",
                sum(is.na(X)), 100 * mean(is.na(X))))
col_means <- colMeans(X, na.rm = TRUE)
na_idx    <- which(is.na(X), arr.ind = TRUE)
X[na_idx] <- col_means[na_idx[, 2]]

X_train  <- X[is_wgs, , drop = FALSE]
y_train  <- factor(meta$species[is_wgs], levels = c("cmn", "wht"))
X_predict <- X[!is_wgs, , drop = FALSE]

# ---- Cross-validate n.pca within the WGS training set ------------------
set.seed(42)
# n.pca.max must be < n_train; keep some headroom for xval folds.
n_pca_max <- min(30, sum(is_wgs) - 2)
message(sprintf("Cross-validating n.pca on WGS samples (n=%d, n.pca.max=%d)...",
                sum(is_wgs), n_pca_max))
xval <- xvalDapc(X_train, y_train,
                 n.pca.max     = n_pca_max,
                 training.set  = 0.9,
                 result        = "groupMean",
                 center        = TRUE,
                 scale         = FALSE,
                 n.rep         = 100,
                 xval.plot     = FALSE)
n_pca_opt <- as.integer(xval$`Number of PCs Achieving Lowest MSE`)
message(sprintf("Optimal n.pca: %d", n_pca_opt))

# ---- Fit DAPC on WGS samples only --------------------------------------
# dapc.matrix takes the raw matrix; pass training labels directly.
dapc_res <- dapc(X_train, grp = y_train, n.pca = n_pca_opt, n.da = 1)

train_acc <- mean(as.character(dapc_res$assign) == as.character(y_train))
message(sprintf("Training assignment accuracy (WGS): %.1f%%", 100 * train_acc))

# ---- Predict RAD samples as held-out ------------------------------------
pred <- predict.dapc(dapc_res, newdata = X_predict)

# ---- Combine posteriors into one tidy table -----------------------------
post_train   <- as.data.frame(dapc_res$posterior)
post_predict <- as.data.frame(pred$posterior)
post_all     <- rbind(post_train, post_predict)
post_all     <- post_all[match(sample_ids, rownames(post_all)), , drop = FALSE]

pop_levels <- c("CB", "CL_RAD", "CL_WGS", "SR")

score_df <- data.frame(
  sample_id  = sample_ids,
  species    = factor(meta$species, levels = c("cmn", "wht")),
  population = meta$population,
  p_wht      = as.numeric(post_all[, "wht"]),
  method     = ifelse(is_wgs, "WGS", "RAD"),
  set        = ifelse(is_wgs, "train (a priori)", "predict (held-out)"),
  stringsAsFactors = FALSE
) %>%
  mutate(population = ifelse(population == "CL",
                             paste0("CL_", method),
                             population),
         population = factor(population, levels = pop_levels),
         set        = factor(set, levels = c("train (a priori)",
                                             "predict (held-out)")))

# ---- P(white) density + rug plot, split by train vs predict ------------
p_density <- ggplot(score_df, aes(x = p_wht, fill = species, colour = species)) +
  geom_density(alpha = 0.45, linewidth = 0.8) +
  geom_rug(alpha = 0.5, linewidth = 0.5) +
  facet_wrap(~ set, ncol = 1, scales = "free_y") +
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
       title  = sprintf(
         "DAPC trained on WGS only  (n.pca = %d, train acc = %.0f%%)",
         n_pca_opt, 100 * train_acc)) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

ggsave(file.path(fig_dir, "dapc_rad_wgs_scores.pdf"),
       p_density, width = 7, height = 6)

# ---- Posterior assignment bar plot, faceted by population --------------
post_df <- post_all %>%
  rownames_to_column("sample_id") %>%
  left_join(score_df[, c("sample_id", "species", "population", "p_wht", "set")],
            by = "sample_id") %>%
  arrange(population, p_wht) %>%
  group_by(population) %>%
  mutate(within_pos = row_number()) %>%
  ungroup() %>%
  pivot_longer(cols = c("cmn", "wht"),
               names_to = "assigned_to", values_to = "prob")

p_assign <- ggplot(post_df, aes(x = within_pos, y = prob, fill = assigned_to)) +
  geom_col(width = 1, colour = "grey50", linewidth = 0.05) +
  facet_grid(. ~ population, scales = "free_x", space = "free_x",
             switch = "x") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(cmn = "steelblue", wht = "white"),
                    labels = c(cmn = "Common", wht = "White"),
                    breaks = c("cmn", "wht")) +
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

# ---- DAPC scatter: LD1 vs PC1, ellipses by population ------------------
# Only one DA exists (two groups), so pair LD1 with PC1 of the genotype
# matrix to get a 2D PCA-style view of the discriminant space.
pc1_all <- prcomp(X, center = TRUE, scale. = FALSE)$x[, 1]

ld1_all <- numeric(length(sample_ids))
names(ld1_all) <- sample_ids
ld1_all[rownames(dapc_res$ind.coord)] <- dapc_res$ind.coord[, 1]
ld1_all[rownames(pred$ind.scores)]    <- pred$ind.scores[, 1]

scatter_df <- data.frame(
  sample_id  = sample_ids,
  LD1        = ld1_all[sample_ids],
  PC1        = pc1_all[sample_ids],
  p_wht      = as.numeric(post_all[, "wht"]),
  population = score_df$population,
  set        = score_df$set,
  stringsAsFactors = FALSE
)

# Fillable shapes per population (21-25 accept both fill and outline colour).
pop_shapes <- setNames(c(21, 22, 23, 24)[seq_along(pop_levels)], pop_levels)

p_scatter <- ggplot(scatter_df,
                    aes(x = LD1, y = PC1,
                        fill = p_wht, shape = population)) +
  geom_point(size = 2.8, alpha = 0.95, stroke = 0.4, colour = "grey20") +
  scale_fill_gradient(low      = "steelblue",
                      high     = "white",
                      limits   = c(0, 1),
                      breaks   = c(0, 0.25, 0.5, 0.75, 1),
                      name     = "P(white)") +
  scale_shape_manual(values = pop_shapes, name = "Population") +
  guides(fill  = guide_colourbar(order = 1, barheight = 6),
         shape = guide_legend(order = 2,
                              override.aes = list(fill = "grey60",
                                                  colour = "grey20"))) +
  labs(x     = "DAPC LD1",
       y     = "Genotype PC1",
       title = "DAPC discriminant axis vs genotype PC1") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "dapc_rad_wgs_scatter.pdf"),
       p_scatter, width = 8.5, height = 6)

# ---- RAD-only scatter: PC1 recomputed on RAD subset --------------------
# LD1 is held-out from the WGS-trained DAPC; PC1 is recomputed using only
# RAD samples so the y-axis reflects variation among RAD individuals alone.
X_rad   <- X[!is_wgs, , drop = FALSE]
pc1_rad <- prcomp(X_rad, center = TRUE, scale. = FALSE)$x[, 1]

scatter_rad_df <- scatter_df %>%
  filter(set == "predict (held-out)") %>%
  mutate(PC1 = pc1_rad[sample_id])

rad_pop_levels <- intersect(pop_levels, unique(as.character(scatter_rad_df$population)))
rad_shapes     <- pop_shapes[rad_pop_levels]

p_scatter_rad <- ggplot(scatter_rad_df,
                        aes(x = LD1, y = PC1,
                            fill = p_wht, shape = population)) +
  geom_point(size = 2.8, alpha = 0.95, stroke = 0.4, colour = "grey20") +
  scale_fill_gradient(low      = "steelblue",
                      high     = "white",
                      limits   = c(0, 1),
                      breaks   = c(0, 0.25, 0.5, 0.75, 1),
                      name     = "P(white)") +
  scale_shape_manual(values = rad_shapes, name = "Population",
                     breaks = rad_pop_levels) +
  guides(fill  = guide_colourbar(order = 1, barheight = 6),
         shape = guide_legend(order = 2,
                              override.aes = list(fill = "grey60",
                                                  colour = "grey20"))) +
  labs(x     = "DAPC LD1 (WGS-trained, projected)",
       y     = "PC1 (RAD samples only)",
       title = "RAD-only: DAPC LD1 vs RAD-subset PC1") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "dapc_rad_wgs_scatter_radonly.pdf"),
       p_scatter_rad, width = 8.5, height = 6)

# ---- Combined figure: scatter (top) + structure-like assignment (bottom)
# Top panel shows held-out RAD samples only (CB and CL_WGS = WGS training
# samples are dropped from the scatter but remain in the bottom-panel bars).
scatter_combined_df <- scatter_df %>%
  filter(!population %in% c("CB", "CL_WGS"))

combined_pop_levels <- intersect(pop_levels,
                                 as.character(unique(scatter_combined_df$population)))
combined_shapes     <- pop_shapes[combined_pop_levels]

p_scatter_combined <- ggplot(scatter_combined_df,
                             aes(x = LD1, y = PC1,
                                 fill = p_wht, shape = population)) +
  geom_point(size = 3.4, alpha = 0.95, stroke = 0.4, colour = "grey20") +
  scale_fill_gradient(low      = "steelblue",
                      high     = "white",
                      limits   = c(0, 1),
                      breaks   = c(0, 0.25, 0.5, 0.75, 1),
                      name     = "P(white)") +
  scale_shape_manual(values = combined_shapes, name = "Population",
                     breaks = combined_pop_levels) +
  guides(fill  = guide_colourbar(order = 1, barheight = 6),
         shape = guide_legend(order = 2,
                              override.aes = list(fill = "grey60",
                                                  colour = "grey20"))) +
  labs(x = "DAPC LD1", y = "Genotype PC1") +
  theme_bw(base_size = 13) +
  theme(legend.position       = c(0.985, 0.985),
        legend.justification   = c(1, 1),
        legend.background      = element_rect(colour = "black",
                                              fill = "white",
                                              linewidth = 0.3),
        legend.box.margin      = margin(2, 2, 2, 2),
        legend.margin          = margin(4, 6, 4, 6),
        legend.spacing.y       = unit(2, "pt"),
        legend.key             = element_blank())

p_combined <- p_scatter_combined /
              (p_assign + labs(title = NULL)) +
  plot_layout(heights = c(2.2, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(fig_dir, "dapc_rad_wgs_combined.pdf"),
       p_combined, width = 9, height = 9)

# ---- Combined figure (all samples): scatter includes WGS training ------
p_scatter_combined_all <- ggplot(scatter_df,
                                 aes(x = LD1, y = PC1,
                                     fill = p_wht, shape = population)) +
  geom_point(size = 3.4, alpha = 0.95, stroke = 0.4, colour = "grey20") +
  scale_fill_gradient(low      = "steelblue",
                      high     = "white",
                      limits   = c(0, 1),
                      breaks   = c(0, 0.25, 0.5, 0.75, 1),
                      name     = "P(white)") +
  scale_shape_manual(values = pop_shapes, name = "Population",
                     breaks = pop_levels) +
  guides(fill  = guide_colourbar(order = 1, barheight = 6),
         shape = guide_legend(order = 2,
                              override.aes = list(fill = "grey60",
                                                  colour = "grey20"))) +
  labs(x = "DAPC LD1", y = "Genotype PC1") +
  theme_bw(base_size = 13) +
  theme(legend.position       = c(0.985, 0.985),
        legend.justification   = c(1, 1),
        legend.background      = element_rect(colour = "black",
                                              fill = "white",
                                              linewidth = 0.3),
        legend.box.margin      = margin(2, 2, 2, 2),
        legend.margin          = margin(4, 6, 4, 6),
        legend.spacing.y       = unit(2, "pt"),
        legend.key             = element_blank())

p_combined_all <- p_scatter_combined_all /
                  (p_assign + labs(title = NULL)) +
  plot_layout(heights = c(2.2, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(fig_dir, "dapc_rad_wgs_combined_all.pdf"),
       p_combined_all, width = 9, height = 9)

message("Figures written:")
message("  figures/dapc_rad_wgs_scores.pdf")
message("  figures/dapc_rad_wgs_assignment.pdf")
message("  figures/dapc_rad_wgs_scatter.pdf")
message("  figures/dapc_rad_wgs_scatter_radonly.pdf")
message("  figures/dapc_rad_wgs_combined.pdf")
message("  figures/dapc_rad_wgs_combined_all.pdf")
