## 1. Partition data into separate .vcf files (one for each chromosome)

The computational time of preprocessing big .vcf files can be substantial. That is why, working with one .vcf file for each chromosome is in general better than treating the whole dataset at once. Moreover, splitting the dataset into separate files allows processing them in parallel.

In order to split the dataset in a .vcf file per chromosome, you can use the script in `utils/partition_and_rename_chr.py`. 
