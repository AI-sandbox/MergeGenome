# Introduction to merging genomic datasets

The small number of samples present in some datasets may be insufficient for analysis. This is especially the case of deep neural networks, which in general require large amounts of data for training. There are some ways to deal with the paucity of data without the need to obtain more real-world data. An example of this is generating synthetic data. However, if you have access to multiple datasets that share some SNPs, the simplest approach is to merge them by the common markers (i.e. all SNPs with identical CHROM, POS, REF, and ALT fields). The main drawback of merging datasets is that if a dataset contains much shorter DNA sequences than the others, a big amount of SNPs are gonna be lost. A way of avoiding this is through **imputation**, for example, of missing variables of the shorter dataset to have the same SNPs as in the larger dataset -the latter called **reference**-.

The limitations in the imputation algorithms and a lack of consensus in the way DNA sequences are read -among others-, make the merging process hard in some cases. This repository underlies the importance of cleaning genomic sequences prior to analysis, and explains the sequence of processing steps for merging a shorter dataset with a larger reference dataset.

# Preprocessing steps

1. **[Partition data into separate .vcf files (one for each chromosome)]***(readmes/README_1_partition_into_separate_files.md)

2. **[Standardize chromosome notation]**(readmes/README_2_change_chrom_notation.md)
