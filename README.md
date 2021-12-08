# Introduction to merging genomic datasets

The small number of samples present in some datasets may be insufficient for analysis. This is especially the case of deep neural networks, which in general require large amounts of data for training. There are some ways to deal with the paucity of data without the need to obtain more real-world data. An example of this is generating synthetic data. However, if we have access to multiple datasets that share some SNPs, the simplest approach is to merge them by the common SNPs. The main drawback of merging datasets is that, if a dataset contains much shorter DNA sequences than the others, a big amount of SNPs will be lost. A way of avoiding this is through **imputation**, for example, of missing variables of the shorter dataset to have the same SNPs as in the larger dataset, the latter called **reference**.

The limitations in the imputation algorithms and a lack of consensus in the way DNA sequences are read (among others), make the merging process hard in some cases. This repository underlies the importance of cleaning genomic sequences prior to analysis, and explains the steps for merging a shorter dataset with a larger dataset.

# Preprocessing steps

## 1. Partition data into separate .vcf files (one for each chromosome)

The computational time of preprocessing big .vcf files can be substantial. That is why, working with one .vcf file for each chromosome is in general better than treating the whole dataset at once. Moreover, splitting the dataset into separate files allows processing them in parallel.

In order to split the dataset in a .vcf file per chromosome, you can use the script in  `utils/partition_and_rename_chr.py`. 

## 2. Rename chromosome nomenclature (from "1" to "chr1")


TODO
