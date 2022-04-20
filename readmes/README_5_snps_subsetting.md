## Subset

There are multiple statistical challenges of high dimensional data. The fact that not all variants are equally important sometimes lead to select the most relevant features for further analysis. In case of having a third "subsetted" dataset with a pre-selection of the important features, you can subset your your data to the common features.

MergeGenome subset command filters the common markers between a query and a reference dataset with data from a specific chromosome.

## Usage

```
$ python3 MergeGenome.py subset -r <reference_file_1>...<reference_file_n> -q <query_file_1>...<query_file_n> -o <output_folder>
```

Input flags include:

* -r, --reference LIST, Paths to reference .vcf files with data for each chromosome (required).
* -q, --query LIST, Path to query .vcf files with data for each chromosome (required).
* -o, --output-folder PATH, Path to output folder to store the modified VCF files (required).
* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Examples**

1. Keep the common markers between the reference and the query:

```
python3 MergeGenome.py subset -r reference_chr1.vcf -q query_chr1.vcf -o ./output/
```
