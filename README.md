# Preprocessing .vcf files

In process...

##### 1. Partitioning data into separate .vcf files and renaming chromosome nomenclature

First of all, it is important to keep in mind that the computational time of preprocessing .vcf files can be substantial. The more SNPs to preprocess, the larger the computational time. For that reason, working with one .vcf file for each chromosome seems better than treating the whole dataset at once. Moreover, a further benefit of splitting the data into separate files is that it allows you to preprocess them in parallel.

For the above reasons, in case your data is not already saved in different files, the first step will be to split it so that each new .vcf file contains the data of a particular chromosome. 

In process...
