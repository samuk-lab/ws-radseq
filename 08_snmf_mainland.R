# sNMF on the mainland-only subset (Cape Breton allopatric commons dropped).
# Input: ns_rad_pruned_mainland.{ped,geno} produced by 07_snmf_prep_mainland.sh
# Output: figures/snmf_mainland_*.pdf

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("LEA")

library("LEA")
library("tidyverse")

adm_dir   <- "data/admixture"
fig_dir   <- "figures"
base      <- file.path(adm_dir, "ns_rad_pruned_mainland")
ped_file  <- paste0(base, ".ped")
geno_file <- paste0(base, ".geno")
Ks        <- 1:6

if (!file.exists(geno_file)) ped2geno(ped_file, output.file = geno_file)

project_file <- paste0(geno_file, ".snmfProject")
if (file.exists(project_file)) {
  obj <- load.snmfProject(project_file)
} else {
  obj <- snmf(geno_file, K = Ks, entropy = TRUE, repetitions = 100,
              project = "new", CPU = 6, seed = 42)
}

ce_df <- map_dfr(Ks, function(K) {
  ce <- cross.entropy(obj, K = K)
  tibble(K = K, run = seq_along(ce), CE = as.numeric(ce))
})
best_run <- ce_df %>% group_by(K) %>% slice_min(CE, n = 1) %>% ungroup()
best_K   <- best_run$K[which.min(best_run$CE)]

ce_plot <- ggplot(ce_df, aes(x = K, y = CE)) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  geom_line(data = best_run, color = "steelblue", linewidth = 1) +
  geom_point(data = best_run, color = "steelblue", size = 3) +
  geom_point(data = filter(best_run, K == best_K),
             color = "red", size = 4) +
  scale_x_continuous(breaks = Ks) +
  labs(title = paste0("sNMF mainland-only cross-entropy (best K = ", best_K, ")"),
       x = "K", y = "Cross-entropy") +
  theme_bw(base_size = 14)

ggsave(file.path(fig_dir, "snmf_mainland_cross_entropy.pdf"),
       ce_plot, width = 6, height = 4)

sample_ids <- read.table(ped_file)[[2]]
popmap <- read.table("meta/popmap.txt", header = TRUE,
                     stringsAsFactors = FALSE) %>%
  select(sample_id, population)
pop_levels <- c("CB", "CL_RAD", "CL_WGS", "SR")

q_wide <- map_dfr(Ks, function(K) {
  br <- best_run$run[best_run$K == K]
  q  <- as.data.frame(Q(obj, K = K, run = br))
  colnames(q) <- paste0("anc", seq_len(K))
  q$sample_id <- sample_ids
  q$K <- K
  q
}) %>%
  left_join(popmap, by = "sample_id") %>%
  mutate(method = ifelse(grepl("_L002", sample_id), "WGS", "RAD"),
         population = ifelse(population == "CL",
                             paste0("CL_", method),
                             as.character(population)),
         population = factor(population, levels = pop_levels))

sort_pos <- q_wide %>%
  group_by(K, population) %>%
  group_modify(~ {
    Kval <- .y$K
    anc_cols <- paste0("anc", seq_len(Kval))
    anc_means <- sapply(anc_cols, function(a) mean(.x[[a]]))
    anc_rank  <- names(sort(anc_means, decreasing = TRUE))
    .x %>% arrange(across(all_of(anc_rank), desc)) %>%
      mutate(within_pos = row_number())
  }) %>%
  ungroup() %>%
  select(K, population, sample_id, within_pos)

q_long <- q_wide %>%
  left_join(sort_pos, by = c("K", "population", "sample_id")) %>%
  pivot_longer(cols = starts_with("anc"),
               names_to = "ancestry", values_to = "prop") %>%
  filter(!is.na(prop)) %>%
  mutate(K_label = paste0("K = ", K))

snmf_plot <- ggplot(q_long, aes(x = within_pos, y = prop, fill = ancestry)) +
  geom_col(width = 1) +
  facet_grid(K_label ~ population, scales = "free_x", space = "free_x",
             switch = "y") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_bw(base_size = 12) +
  theme(axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        panel.spacing.x  = unit(0.15, "lines"),
        strip.background = element_blank(),
        strip.placement  = "outside",
        legend.position  = "none")

ggsave(file.path(fig_dir, "snmf_mainland_K1-6_sorted.pdf"),
       snmf_plot, width = 10, height = 8)

message(sprintf("mainland best K: %d", best_K))
