## Subset SNPs to another dataset

There are multiple statistical challenges of high dimensional data. The fact that not all variants are equally important sometimes leads to select the most relevant features for further analysis and discard the often uninformative features. In case of having a third "subsetted" reference dataset with a pre-selection of important features or SNPs, it is usually beneficiat to subset the high-dimensional data to those relevant SNPs.

MergeGenome subset command finds the common markers (i.e. SNPs with identical CHROM, POS, REF, and ALT fields). between a query and a reference dataset containing data from a particular chromosome and subsets both datasets to only contain those common markers.

## Usage

```
$ python3 MergeGenome.py subset -q <query_file_1>...<query_file_n> -r <reference_file_1>...<reference_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Path to query .vcf files with data for for a single chromosome each (required).
* -r, --reference LIST, Paths to reference .vcf files with data for for a single chromosome each (optional).
* -o, --output-folder PATH, Path to output folder (required). Note: make sure a '/' appears at the end of the output folder.
* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Output**

* One .vcf file with subsetted SNPs for each input query and reference files. Each new .vcf file will receive the same base name as the input file, but ending with '_subsetted.vcf'.
* If --debug, a .log or .txt file with information regarding the dimensions of the data (number of samples and number of SNPs), the chromosome in each file, and the amount of samples and SNPs after filtering to the common markers.

**Example**

1. Find and subset to the common markers between the reference and the query files:

```
$ python3 MergeGenome.py subset -q query_chr1.vcf query_chr2.vcf -r reference_chr1.vcf reference_chr2.vcf -o ./output/
```
