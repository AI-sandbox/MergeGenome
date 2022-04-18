# MergeGenome - Toolkit for Merging VCF files

This repository includes a Python implementation of the MergeGenome toolkit, which underlies the importance of cleaning genomic sequences prior to analysis. The MergeGenome toolkit is designed to integrate DNA sequences from a query and a reference datasets in variant call format (VCF) while targeting data quality. MergeGenome is a robust pipeline of comprehensive steps to merge both datasets, including chromosome nomenclature standardization, SNP ambiguities removal, SNP flips detection, SNP mismatches elimination, and query/reference mismatches detection and, optionally, fixing.

This repository also includes the implementation of other common tasks related to merging genomic sequences, such as identifying the common markers (i.e. SNPs with identical CHROM, POS, REF, and ALT fields) between two datasets and subsetting the available data to those common markers.

# Steps & usage

1. **[Partition data into separate VCF files (one per chromosome)](readmes/README_1_partition_into_separate_files.md)**

2. **[Change chromosome notation](readmes/README_2_change_chrom_notation.md)**

3. **[Preprocess .vcf files](readmes/README_3_preprocess_vcf_files.md)**

4. **[Impute SNPs (from short to large dataset)](readmes/README_4_impute.md)**

5. **[Subset SNPs to another dataset](readmes/README_5_snps_subsetting.md)**

6. **[Concatenate and merge .vcf files](readmes/README_6_concatenate_and_merge_vcf_files.md)**

# Verification steps & usage

1. **[Plot of PCA/UMAP/tSNE](readmes/README_1_make_pca_plots.md)**

3. **[Plot of SNPs means](readmes/README_2_plot_snps_means.md)**
