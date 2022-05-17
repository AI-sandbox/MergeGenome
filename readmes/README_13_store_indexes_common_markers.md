## Store indexes common markers

In some occasions, it is useful to identify which SNPs are or not common markers between a query and a reference datasets. For instance, when merging two datasets, after inducing shorter DNA sequences to have all the variants from the longer DNA sequences, it is desirable to study the homogeneity of the resulting dataset. For such analysis, it might be interesting to distinguish between the SNPs that were initially common markers (i.e., SNPs at the same CHROM, POS, REF, and ALT) and the SNPs that became common markers after imputation. 

The MergeGenome store-common-indexes command searches the indexes of the SNPs of the query dataset that are also present in the reference dataset and viceversa. The indexes of the common markers in each VCF file are saved in a NPY or H5 file.

The provided query and reference .vcf files must contain data for the same chromosome, in the same chromosomical order.

## Usage

```
$ python3 MergeGenome.py store-common-indexes -q <query_file> -r <reference_file> -o <output_folder>
```

Input flags include:

* -q, --query Path, Path to query .vcf file with data for a single or multiple chromosomes (required).
* -r, --reference Path, Path to reference .vcf file with data for a single or multiple chromosomes (required).
* -o, --output-folder PATH, Path to output folder (required). Note: make sure a '/' appears at the end of the output folder.
* -f, --file-format STR, Format of the output file (required). Options: '.npy', to store data in numpy format; '.h5' to store data in h5py format. Default='.npy'.
* -d, --debug PATH, Path to .log/.txt file to store info/debug messages (optional).

**Output**

* Two .npy or .h5 files, with the indexes of the common markers in the query and reference .vcf files. Each new file will receive the same base name as the input files, but ending with '_indexes.npy' or '_indexes.h5'.
* If --debug, a .log or .txt file with information regarding the dimensions of the data (number of samples and number of SNPs), the chromosome in particular, and the amount of common markers.

`Examples`

1. Store indexes of common markers between the query and the reference in .npy format:

```
$ python3 MergeGenome.py store-common-indexes -q query.vcf -r reference.vcf -o ./output/ -f .npy
```

2. Store indexes of common markers between the query and the reference in .h5 format:

```
$ python3 MergeGenome.py store-common-indexes -q query.vcf -r reference.vcf -o ./output/ -f .h5
```
