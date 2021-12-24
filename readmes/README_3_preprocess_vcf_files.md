## 3. Preprocess .vcf files

Following, are explained some important preprocessing steps to clean genomic data:

* **Remove undesired samples by sample ID (breed name)**: there are some breeds or species that are very different from the other samples in the dataset. These can be identified doing a 2-D PCA plot as outliers. It is better to remove these samples before imputation.

* **Remove ambiguous SNPs**: A/T, T/A, C/G, and G/C pairs are strand ambiguous since they are complementary, and it can be challenging to determine the DNA strand given the lack of consensus on how DNA sequences are read.

* **Correct SNP flips**: flipped SNPs between two datasets are positions in which the reference and the alternate are swapped. That is, the reference in one dataset is the alternate in the other dataset and vice-versa. They do not need to be removed because they can be easily corrected. Correcting a SNP flip between two datasets consists in swapping the reference by the alternate and vice-versa in one of the datasets. Thus, the 0's also must be changed by 1's and vice-versa. 

* **Remove SNP mismatches**: a SNP mismatch between two datasets takes place when the reference or alternate of a SNP at the same position differs among datasets. If SNP flips are not corrected previous to removing SNP mismatches, they will be lost. Hence, it is recommended to correct all SNP flips before removing SNP mismatches.

* **Rename missings**: some widely used software packages such as Beagle can give the following error: *"Caused by: java.lang.IllegalArgumentException: ERROR: invalid allele [-1]"* if missing allele values are encoded as -1. Changing the missing nomenclature to a dot "." should solve this issue.

In order to run all previous preprocessing steps, you can use the script in `scripts/vcf_preprocessing.py`.

