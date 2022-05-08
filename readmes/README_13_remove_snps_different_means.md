## Remove SNPs with different means

Another approach to overcome data heterogeneity is to remove the SNPs with very different means between the query and the reference.

The remove-snps-different-means from the MergeGenome toolkit removes all the common SNPs with a mean absolute difference higher than specified threshold, and stres the result in a new VCF file for each chromosome.

## Usage

```
$ python3 MergeGenome.py remove-snps-different-means -r <reference_file_1>...<reference_file_n> -q <query_file_1>...<query_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Path to query .vcf files with data for each chromosome (required).
* -r, --reference LIST, Paths to reference .vcf files with data for each chromosome (required).
* -o, --output-folder PATH, Path to output folder to store the modified VCF files (required). Note: make sure a '/' appears at the end of the output folder.
* -t, --threshold FLOAT, All common SNPs with a mean absolute difference higher than the threshold will be removed (optional).
* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Example**

1. Remove all the common SNPs with a mean absolute difference higher than 0.1:

```
$ python3 MergeGenome.py remove-snps-different-means -q query_chr*.vcf -r reference_chr*.vcf -o ./output/ -t 0.1
```