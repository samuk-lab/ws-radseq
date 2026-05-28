#install.packages("devtools")
#library("devtools") 
#devtools::install_github("DevonDeRaad/SNPfiltR") 
library("SNPfiltR")
library("vcfR")

dir.create("data/vcf/filtered")

vcf_file <- "ns_rad_prefiltered_no_M.vcf.gz"

output_file <- "data/vcf/ns_rad_filtered.vcf.gz"

# apply filters
# this can actually be done using a pipe, but doing it this way allows the 
# object to reduce itself in size as you go (preventing out of memory issues)
# some steps are repeated to take advantage of progressively reducing memory
vcf_filtered <- read.vcfR(paste0("data/vcf/", vcf_file))
vcf_filtered <- filter_biallelic(vcf_filtered)
vcf_filtered <- hard_filter(vcf_filtered, depth = 5, gq = 20) 
vcf_filtered <- missing_by_snp(vcf_filtered, cutoff = 0.85) 
vcf_filtered <- missing_by_sample(vcf_filtered, cutoff = .9)
vcf_filtered <- filter_allele_balance(vcf_filtered)
vcf_filtered <- max_depth(vcf_filtered, maxdepth = 40)
vcf_filtered <- min_mac(vcf_filtered, min.mac = 1) 
vcf_filtered <- missing_by_sample(vcf_filtered, cutoff = .9)
vcf_filtered <- missing_by_snp(vcf_filtered, cutoff = 0.85) 
#vcf_filtered <- distance_thin(vcf_filtered, min.distance = 500)
vcf_filtered <- min_mac(vcf_filtered, min.mac = 1) 

write.vcf(vcf_filtered, output_file)

system("gunzip data/vcf/ns_rad_filtered.vcf.gz")
system("bgzip data/vcf/ns_rad_filtered.vcf")
system("tabix data/vcf/ns_rad_filtered.vcf.gz")
system("bcftools view -O z -o data/vcf/all_chr_filtered_only_Y.vcf.gz --regions chrY data/vcf/ns_rad_filtered.vcf.gz")


# mito genome
# not included

# vcf_file <- "chrM_genotyped.vcf"
# output_file <- gsub("_genotyped", "_filtered", vcf_file) %>% paste0("data/vcf/filtered/", ., ".gz")
# 
# vcf_filtered <- read.vcfR(paste0("data/vcf/", vcf_file)) %>%
#   filter_biallelic() %>% 
#   hard_filter(depth = 5, gq = 20) %>%
#   filter_allele_balance() %>%
#   missing_by_sample(cutoff = 0.90) %>%
#   min_mac() %>%
#   missing_by_snp() %>%
#   distance_thin(min.distance = 500)
# 
# output_file <- gsub("_genotyped", "_filtered", vcf_file) %>% paste0("data/vcf/filtered/", ., ".gz")




