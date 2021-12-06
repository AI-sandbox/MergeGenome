# preprocessing-vcf-files

First of all, it is important to keep in mind that the computational time of preprocessing .vcf files can be substantial. The more SNPs to preprocess, the larger the computational time. For that reason, working with one .vcf file for each chromosome is better than treating the whole dataset at once. Moreover, a further benefit of splitting the data in separate files is that it allows us to preprocess them in parallel.

For the above reasons, in case your data is not already saved in different .vcf files, the first step will be to split it so that each new .vcf file contains the data of a particular chromosome. The script that performs this step is **/scripts/partition_and_rename_chr.py**.
