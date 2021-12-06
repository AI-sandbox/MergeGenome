# Preprocessing .vcf files

The small number of samples present in some datasets may be insufficient for analysis. This is especially the case of deep neural networks, which in general require large amounts of data for training. There are some ways to deal with the paucity of data without the need to obtain more real-world data. An example of this is generating synthetic data. However, if we have access to multiple datasets that share some SNPs, the simplest approach is to merge them by the common SNPs. In the case of having a dataset with much shorter genomic sequences in length, prior to merging, we can impute the missing values to have the same SNPs as in the larger reference dataset.

However, the merging process can be hard in some cases due to limitations in the imputation algorithms and a lack of consensus in the way DNA sequences are read, among others. 

In this repository, we explain a possible pipeline of (in process...)

### 1. Partition data into separate .vcf files (one for each chromosome)

First of all, it is important to keep in mind that the computational time of preprocessing .vcf files can be substantial. The more SNPs to preprocess, the larger the computational time. For that reason, working with one .vcf file for each chromosome seems better than treating the whole dataset at once. Moreover, a further benefit of splitting the data into separate files is that it allows you to preprocess them in parallel.

For the above reasons, in case your data is not already saved in different files, the first step will be to split it so that each new .vcf file contains the data of a particular chromosome.

### 2. Rename chromosome nomenclature (from "1" to "chr1")


In process...
