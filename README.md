# MergeGenome - Toolkit for Merging VCF files

This repository includes a Python implementation of the MergeGenome toolkit, which underlies the importance of cleaning genomic sequences prior to analysis. The MergeGenome toolkit is designed to integrate DNA sequences from a query and a reference datasets in variant call format (VCF) while targeting data quality. MergeGenome is a robust pipeline of comprehensive steps to merge both datasets, including chromosome nomenclature standardization, SNP ambiguities removal, SNP flips detection, SNP mismatches elimination, and query/reference mismatches detection and/or fixing. MergeGenome works with any organism’s DNA sequences, which brings a broad solution to having access to more than one source of data but only being able to exploit one for statistical analysis.

This repository also includes the implementation of other common tasks related to merging genomic sequences, such as identifying the common markers (i.e. SNPs with identical CHROM, POS, REF, and ALT fields) between two datasets and subsetting the available data to those common markers.

# Preprocessing Steps

1. **[Partition data into a separate VCF file per chromosome](readmes/README_1_partition_into_separate_files.md)**

2. **[Rename chromosome notation](readmes/README_2_rename_chrom_notation.md)**

3. **[Clean VCF files](readmes/README_3_clean_vcf_files.md)**

4. **[Impute](readmes/README_4_impute.md)**

5. **[Subset SNPs to common markers with another dataset](readmes/README_5_snps_subsetting.md)**

6. **[Data Homogenization: Detect and fix mismatches](readmes/README_6_detect_and_fix_mismatches.md)**

7. **[Concat VCF files](readmes/README_7_concat_vcf_files.md)**

8. **[Merge VCF files](readmes/README_8_merge_vcf_files.md)**

# Merging Evaluation

9. **[Data Homogenization: Discriminator Performance](readmes/README_9_evaluate_with_discriminator.md)**

10. **[Plot SNP means (comparison)](readmes/README_10_evaluate_with_snp_means_plot.md)**

11. **[Plot Principal Component Analysis (PCA)](readmes/README_11_evaluate_with_pca_plots.md)**

# Other util commands

12. **[Remove SNPs with different means](readmes/README_12_remove_snps_different_means.md)**

13. **[Store indexes common markers](readmes/README_13_store_indexes_common_markers.md)**

14. **[Store allele data (from VCF to NPY or H5)](readmes/README_14_store_allele_data_npy_h5.md)**

15. **[Store allele data (from NPY to VCF)](readmes/README_15_store_allele_data_vcf.md)**
