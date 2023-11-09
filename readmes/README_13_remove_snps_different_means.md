## Remove SNPs with different means

The MergeGenome remove-snps-different-means command removes all the common markers (i.e., SNPs at the same CHROM, POS, REF, and ALT) between a query and reference datasets that have a mean absolute difference higher than a threshold (by default, 0.1). It is a fast solution to identifying and removing outlier SNPs, which might harm any analysis built on top of corrupt features. The query and reference with the filtered SNPs are stored in a new VCF file.

The provided .vcf files must contain data for the same single chromosome, and the files need to appear in the same chromosomical order.

## Usage

```
$ python3 MergeGenome.py remove-snps-different-means -q <query_file_1>...<query_file_n> -r <reference_file_1>...<reference_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Path to query .vcf files with data for for a single chromosome each (required).
* -r, --reference LIST, Paths to reference .vcf files with data for for a single chromosome each (required).
* -o, --output-folder PATH, Path to output folder (required). Note: make sure a '/' appears at the end of the output folder.
* -t, --threshold FLOAT, All common SNPs with a mean absolute difference higher than this threshold will be removed (optional). Accepted values are between 0.0 and 1.0. Default=0.1.
* -d, --debug PATH, Path to .log/.txt file to store info/debug messages (optional).

**Output**

* One .vcf file with filtered data for each input query and reference files. Each new .vcf file will receive the same base name as the input file, but ending with '_removed.vcf'.
* If --debug, a .log or .txt file with information regarding the dimensions of the data (number of samples and number of SNPs), the chromosome in each file, and the amount SNPs removed.

`Example`

1. Remove all the common SNPs with a mean absolute difference higher than 0.2:

```
$ python3 MergeGenome.py remove-snps-different-means -q query_chr1.vcf -r reference_chr1.vcf -o ./output/ -t 0.2
```