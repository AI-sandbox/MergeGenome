## 1. Partition data into separate .vcf files (one for each chromosome) and change chromosome nomenclature

The computational time of preprocessing big .vcf files can be substantial. In general, working with one .vcf file for each chromosome is better than treating the whole dataset at once because:

1. Most of the computations are faster when treating with smaller datasets.

2. Having separate files allows processing the data for each chromosome in parallel.

3. It prevents from having trouble with memory.

### Script

In order to split the dataset in a .vcf file per chromosome, you can use the script in `utils/partition_and_rename_chr.py`. 
