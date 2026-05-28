
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("gdsfmt")
#BiocManager::install("SNPRelate")

library("gdsfmt")
library("vcfR")
library("tidyverse")
library("SNPRelate")
library("plotly")
library("viridis")


if (!file.exists("data/gds/all_chr.gds")){
  
  dir.create("data/gds")
  snpgdsVCF2GDS("data/vcf/ns_rad_filtered.vcf.gz", "data/gds/all_chr.gds", method="biallelic.only")
  
}

# read in genotypes and assign snp ids from ID column (chr_pos)
genofile <- snpgdsOpen("data/gds/all_chr.gds")
geno <- snpgdsGetGeno(genofile, with.id = TRUE)
geno$snp.id <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))

geno_df <- as.data.frame(t(geno$genotype))
names(geno_df) <- geno$sample.id
geno_df <- data.frame(snp_id = geno$snp.id, geno_df) %>%
  gather(key = "sample_id", value = "genotype", -snp_id) %>%
  mutate(chr = gsub("_.*", "", snp_id)) %>%
  mutate(pos = gsub(".*_", "", snp_id) %>% as.numeric) %>%
  select(chr, pos, snp_id, sample_id, genotype)

# create metadata df
meta_df <- data.frame(sample_id = geno$sample.id)
meta_df <- meta_df %>%
  mutate(population = gsub("[0-9]+", "", sample_id)) %>%
  mutate(population = ifelse(population == "sal_riv", "SR", population)) %>%
  mutate(method = ifelse(grepl("^MC|^F", sample_id), "WGS", "RAD"))

sex_df <- read.table("meta/sex_info.txt", header = TRUE)
names(sex_df) <- c("sample_id","sex")

meta_df <- left_join(meta_df, sex_df)

meta_df <- meta_df %>%
  mutate(sex = ifelse(grepl("^MC", sample_id), "M", sex)) %>%
  mutate(sex = ifelse(grepl("^FC", sample_id), "F", sex))

meta_df <- meta_df %>%
  mutate(population = ifelse(grepl("^.CL", sample_id), "CL", population)) %>%
  mutate(population = ifelse(grepl("^.CB|^FB", sample_id), "CB", population))

write.table(meta_df, "meta/meta_df.txt", row.names = FALSE, quote = FALSE)

#pca
autosome_snps <- snpgdsSNPList(genofile, sample.id=NULL) %>%
  filter(!(chromosome %in% c("Y", "XIX"))) %>%
  pull(snp.id)

sex_snps <- snpgdsSNPList(genofile, sample.id=NULL) %>%
  filter(chromosome %in% c("Y", "XIX")) %>%
  pull(snp.id)

all_snps <- snpgdsSNPList(genofile, sample.id=NULL) %>%
  pull(snp.id)

#pca_samples <- meta_df %>%
#  filter(!(sample_id %in% outliers)) %>%
#  filter(grepl("sal|SR|CL", sample_id)) %>%
#  pull(sample_id)

pca_samples <- meta_df %>%
  #filter(!grepl("BC|MJ", population)) %>%
  pull(sample_id)

pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=FALSE, 
                 sample.id = pca_samples, snp.id = all_snps)

tab <- data.frame(sample_id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],  
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)

# putative species assignment

meta_df %>%
  left_join(tab) %>%
  mutate(species = ifelse(EV1 > 0.04, "wht", "cmn")) %>%
  mutate(species = ifelse(is.na(species), "cb", species)) %>%
  select(sample_id, population, species) %>%
  write.table("meta/popmap.txt", row.names = FALSE)

fill_col <- RColorBrewer::brewer.pal(6, "Set1")

pca_dat <- left_join(tab, meta_df) %>%
  filter(!is.na(sex)) %>%
  mutate(population = ifelse(grepl(".*MJ.*", population), "MJ", population)) %>%
  mutate(population = ifelse(grepl("_L002", sample_id), paste0(population,"_WGS"), paste0(population,"_RAD"))) 

# rename EV to PC
pca_dat <- pca_dat %>%
  rename(PC1 = EV1, PC2 = EV2)

### PLOTS

pca_plot_S1 <- pca_dat %>%   
  ggplot(aes(x = PC1, y = PC2, fill = population, shape = sex, label = sample_id))+
  geom_point(size = 4, color = "black")+
  theme_bw(base_size = 20)+
  scale_fill_brewer(palette = "Set1")+ 
  #scale_fill_manual(values = fill_col)+
  scale_shape_manual(values = c(21, 24))+
  guides(fill = guide_legend(override.aes = list(color = fill_col)),
         color = guide_legend(override.aes = list(shape = 21, color = "black")))

# Supp Mat PCA
ggsave("figures/FigureS1.pdf", pca_plot_S1, width = 10, height = 8)


# only random sympatric sites
fill_col <- c("cornflowerblue", "skyblue", "white")
pca_plot_fig5 <- pca_dat %>%   
  filter(grepl("CL_RAD|SR_RAD|MJ_RAD|BC_RAD", population )) %>%
  mutate(species = ifelse(PC1 > 0.04, "White", "Common")) %>%
  mutate(species = ifelse(population == "CL_RAD" & (PC1 < 0.05) & (PC1 > 0.025), "Putative Hybrids", species)) %>%
  #mutate(population = ifelse(population == "CL_RAD", "Canal Lake", "Salmon River")) %>%
  ggplot(aes(x = PC1, y = PC2, fill = species, shape = sex, label = sample_id))+
  geom_point(size = 4, color = "black")+
  facet_wrap(~population)+
  theme_bw(base_size = 20)+
  #scale_fill_manual(values = fill_col)+
  scale_shape_manual(values = c(21, 24))+
  theme(strip.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25),
        legend.position = "inside",
        legend.position.inside = c(0.50, 0.4))+
  guides(fill = guide_legend(override.aes = list(color = fill_col),),
         color = guide_legend(override.aes = list(shape = 19, color = "black")))


ggsave("figures/Figure5.pdf", pca_plot_fig5, width = 8, height = 5)


#### other exploratory checks/figures here
#### not included in manuscript

# CANAL LAKE PCA
pca_samples_CL <- meta_df %>%
  filter(grepl("CL", sample_id)) %>%
  pull(sample_id)


pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=FALSE, 
                 sample.id = pca_samples_CL, snp.id = autosome_snps)

tab <- data.frame(sample_id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],  
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)

fill_col <- RColorBrewer::brewer.pal(3, "Set1")[1:2]

pca_plot <- left_join(tab, meta_df) %>%
  filter(!is.na(sex)) %>%
  mutate(population = ifelse(grepl("_L002", sample_id), "WGS Whites", "Unknown")) %>%
  ggplot(aes(x = EV1, y = EV2, fill = population, shape = sex, label = sample_id))+
  geom_point(size = 4, color = "black")+
  theme_bw(base_size = 20)+
  scale_fill_brewer(palette = "Set1")+
  scale_shape_manual(values = c(21, 24))+
  guides(fill = guide_legend(override.aes = list(color = fill_col)),
         color = guide_legend(override.aes = list(shape = 21, color = "black")))

ggplotly(pca_plot)

left_join(tab, meta_df) %>%
  ggplot(aes(x = EV3, y = EV4, color = population, shape = sex))+
  geom_point(size = 3)

pca_loadings <- snpgdsPCASNPLoading(pca, genofile, num.thread=1L, verbose=TRUE)

loadings <- t(pca_loadings$snploading)
loadings <- left_join(data.frame(snp.id = pca_loadings$snp.id, v2_loading = loadings[,2]), snpgdsSNPList(genofile, sample.id=NULL))
loadings %>%
  select(chromosome, position, v2_loading) %>%
  ggplot(aes(x = position, y = abs(v2_loading)))+
  facet_wrap(~chromosome)+
  geom_point(size = 0.5)


# check sex assignment
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=FALSE, snp.id = sex_snps, sample.id = pca_samples)

tab <- data.frame(sample_id = pca$sample.id,
                  EV1 = pca$eigenvect[,
                                      1],
                  EV2 = pca$eigenvect[,2],  
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)

pca_plot <- left_join(tab, meta_df) %>%
  ggplot(aes(x = EV1, y = EV2, color = population, shape = sex, label = sample_id))+
  geom_point(size = 3)
ggplotly(pca_plot)

#CL 47 is an extreme outlier -- removing


