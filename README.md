# ws_radseq

Population genomic analysis of white and common threespine sticklebacks from Nova Scotia, combining a labelled set of whole-genome sequenced (WGS) individuals with a broader random sample sequenced by RADseq. The analysis characterises population structure, identifies putative admixed individuals at sympatric sites, and provides genomic evidence discussed in the manuscript:

> **"Evidence of parental care as a novel reproductive isolating mechanism"**  
> Behrens et al. (2026, full citation forthcoming.)

---

## Overview

Two forms — referred to here as **white** (`wht`) and **common** (`cmn`) — co-occur at sympatric sites on the Nova Scotia. The WGS data from Canal Lake (CL) and Cherry Burton Road (CB) (from Sumarli et al., full citation forthcoming) contains individuals with confirmed white/common identities, whereas the RAD-seq data are random samples of individuals of unknown identity at Canal Lake and Salmon River (two sites where whites/commons are sympatric).

 A joint SNP dataset was built from the merged RAD+WGS call set and used for:

- **Principal components analysis (PCA)** — primary evidence of form-level structure
- **Ancestry estimation (sNMF)** — mainland-only supplemental analysis identifying individuals with mixed ancestry
- **Discriminant analysis of principal components (DAPC)** — posterior assignment probabilities used to quantify the degree of admixture in candidate hybrid individuals

---

## Repository structure

```
ns_radseq/
├── data/
│   ├── vcf/
│   │   ├── raw/                          # Original merged RAD+WGS joint VCF
│   │   ├── ns_rad_prefiltered_no_M.vcf.gz  # Post-script 01: biallelic SNPs, autosomes + chrY
│   │   ├── ns_rad_filtered.vcf.gz          # Post-script 02: SNPfiltR-filtered, used for PCA
│   │   ├── ns_rad_filtered_v2.vcf.gz       # Post-script 04: split-depth filtered, used for DAPC
│   │   └── ns_rad_filtered_mainland.vcf.gz # Post-script 06: mainland-only, used for sNMF
│   ├── admixture/
│   │   ├── mainland_samples.txt            # Sample list for mainland subsetting
│   │   ├── ns_rad_pruned.{bed,bim,fam,raw} # LD-pruned full dataset (DAPC input)
│   │   ├── ns_rad_pruned_mainland.{bed,bim,fam,ped,geno}  # LD-pruned mainland dataset (sNMF input)
│   │   └── ns_rad_pruned_mainland.snmfProject + .snmf/    # sNMF results (K=1-6, 100 reps)
│   └── gds/                              # GDS format files for SNPRelate PCA
├── meta/
│   ├── popmap.txt     # Sample metadata: sample_id, population, species
│   ├── meta_df.txt    # Extended metadata
│   └── sex_info.txt   # Sex assignment data
├── envs/
│   └── admixture_env.yaml  # Conda environment (plink, bcftools, tabix)
├── figures/           # All output figures (PDF/PNG)
└── *.R / *.sh         # Analysis scripts (see below)
```

---

## Scripts

Scripts are numbered in approximate execution order. Shell scripts (`.sh`) require WSL with the `admixture_env` conda environment; R scripts run natively on Windows.

| Script | Purpose |
|--------|---------|
| `01_prefilter_vcfs.sh` | Initial joint VCF filtering: biallelic SNPs, max-missing 0.7, MAC ≥ 2; removes INFO fields for downstream compatibility |
| `02_filter_vcfs_snpfiltR.R` | SNPfiltR quality filters on the full joint dataset (depth, GQ, allele balance, missingness); produces VCF for PCA |
| `03_pca.R` | PCA via SNPRelate on the filtered VCF; produces the main-text PCA figure |
| `04_filter_rad_wgs.R` | SNPfiltR filtering with split depth thresholds (RAD cap 40×, WGS cap 80×); produces the joint-filtered VCF used by the DAPC pipeline |
| `05_plink_prep_rad_wgs.sh` | Converts filtered VCF to PLINK format, applies MAF ≥ 0.05 and LD pruning (*r*² < 0.2); outputs `ns_rad_pruned.*` used by DAPC |
| `06_filter_mainland.R` | As `04`, but restricted to mainland samples (Cape Breton allopatric commons dropped); used for the supplemental sNMF analysis |
| `07_snmf_prep_mainland.sh` | As `05`, but for mainland-only VCF with relaxed MAF (≥ 0.03) and no LD pruning to maximise power |
| `08_snmf_mainland.R` | sNMF ancestry estimation (K = 1–6, 100 replicates) on the mainland subset; produces supplemental bar plots and cross-entropy figure |
| `09_dapc_rad_wgs.R` | DAPC on the full RAD+WGS dataset (4,678 SNPs, *n* = 224); cross-validates optimal PC retention, plots discriminant scores and posterior assignment probabilities faceted by population |
| `10_extract_paper_values.R` | Utility script: extracts cross-entropy values, sNMF Q scores, and DAPC posteriors for the six hybrid candidates and formats them for in-text reporting |

---

## Dependencies

For running the shell scripts used in the analysis, install the conda environment described in  envs/admixture_env.yaml:

### Conda (WSL)
Install via `conda env create -f envs/admixture_env.yaml`:
- `plink` 1.90b6.21
- `bcftools`
- `tabix`

### R packages
- **Bioconductor:** `LEA`, `SNPRelate`, `gdsfmt`
- **CRAN:** `SNPfiltR`, `vcfR`, `adegenet`, `tidyverse`, `viridis`, `plotly`
