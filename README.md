# MergeGenome - Toolkit for Merging VCF files

This repository includes a Python implementation of the MergeGenome toolkit, which underlies the importance of cleaning genomic sequences prior to analysis. The MergeGenome toolkit is designed to integrate DNA sequences from a query and a reference datasets in variant call format (VCF) while targeting data quality. MergeGenome is a robust pipeline of comprehensive steps to merge both datasets, including chromosome nomenclature standardization, SNP ambiguities removal, SNP flips detection, SNP mismatches elimination, and query/reference mismatches detection and, optionally, fixing.

This repository also includes the implementation of other common tasks related to merging genomic sequences, such as identifying the common markers (i.e. SNPs with identical CHROM, POS, REF, and ALT fields) between two datasets and subsetting the available data to those common markers.

# Steps

1. **[Partition data into a separate VCF file per chromosome](readmes/README_1_partition_into_separate_files.md)**

2. **[Rename chromosome notation](readmes/README_2_rename_chrom_notation.md)**

3. **[Clean VCF files](readmes/README_3_clean_vcf_files.md)**

4. **[Impute](readmes/README_4_impute.md)**

5. **[Subset SNPs to another dataset](readmes/README_5_snps_subsetting.md)**

6. **[Detect and fix mismatches](readmes/README_6_detect_and_fix_mismatches.md)**

7. **[Concat VCF files](readmes/README_7_concat_vcf_files.md)**

8. **[Merge VCF files](readmes/README_8_merge_vcf_files.md)**

# Other util preprocessing commands

9. **[Store_allele_data in NPY or H5 file](readmes/README_9_store_allele_data.md)**

<!---
--- 1. **[Plot of PCA/UMAP/tSNE](readmes/README_1_make_pca_plots.md)**

3. **[Plot of SNPs means](readmes/README_2_plot_snps_means.md)**

--->