# prefilter VCFs prior to input in SNPfiltR
# KS July 2023

# recompress with bgzip
cd data/vcf
mkdir prefiltered

vcftools --vcf raw/merged_rad_wgs.vcf --max-missing 0.7 --mac 2 --min-alleles 2 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out SNPs_only --stdout > prefiltered/ns_rad_prefiltered.vcf

# recompress with bgzip
cd prefiltered
ls *.vcf | xargs -n1 bgzip
ls *.vcf.gz | xargs -n1 tabix

# list and merge vcfs
#ls *.gz  > vcf_list.txt
#bcftools concat -f vcf_list.txt -Oz -o all_chr_prefiltered.vcf.gz

# apply GATK bp filters 
#bcftools view -O z -o ns_rad_prefiltered_bp.vcf.gz -e "MQ < 40.0 || MQRankSum < -12.5" ns_rad_prefiltered.vcf.gz

# remove INFO (only used for BP filters) and unused format fields
bcftools annotate --remove INFO,^FORMAT/DP,^FORMAT/AD,^FORMAT/GQ,^FORMAT/GT -O z -o ns_rad_prefiltered_no_info.vcf.gz ns_rad_prefiltered.vcf.gz

cd ..

mv prefiltered/ns_rad_prefiltered_no_info.vcf.gz .
tabix ns_rad_prefiltered_no_info.vcf.gz

# subset out M & Y chromosome
bcftools view -O z -o ns_rad_prefiltered_no_M.vcf.gz --regions chrIII,chrII,chrIV,chrIX,chrI,chrUn,chrVIII,chrVII,chrVI,chrV,chrXIII,chrXII,chrXIV,chrXIX,chrXI,chrXVIII,chrXVII,chrXVI,chrXV,chrXXI,chrXX,chrX,chrY ns_rad_prefiltered_no_info.vcf.gz
bcftools view -O z -o ns_rad_prefiltered_only_M.vcf.gz --regions chrM ns_rad_prefiltered_no_info.vcf.gz
bcftools view -O z -o ns_rad_prefiltered_only_Y.vcf.gz --regions chrY ns_rad_prefiltered_no_info.vcf.gz

ls *.vcf.gz | xargs -n1 tabix
#rm -r prefiltered

