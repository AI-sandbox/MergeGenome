## Store indexes common markers

In some occasions, it may be useful to know which SNPs are common markers (or not) between two datasets. For instance, when merging two datasets after inducing shorter DNA sequences (query) to have all the variants from the longer DNA sequences (reference), it is desirable to study the homogeneity of the resulting dataset. For such analysis, it might be interesting to distinguish between the SNPs that were initially common markers and the SNPs that became common markers after imputing the query with respect to the reference, as the insights might be different. Knowing the indexes of the merged dataset that were also in the query allows plotting the common markers of different color than the imputed SNPs.

The store-common-indexes command from MergeGenome searches for the indexes of the SNPs of the query dataset that are also present in the reference dataset and viceversa. The two lists with the indexes of the common markers are saved in a NPY or H5 file.

## Usage

```
$ python3 MergeGenome.py store-common-indexes -q <query_file> -r <reference_file> -o <output_folder>
```

Input flags include:

* -q, --query PATH, Path to query .vcf file with data for all chromosome (required).
* -r, --reference PATH, Path to reference .vcf file with data for all chromosome (required).
* -o, --output-folder PATH, Path to output folder to store the modified VCF files (required). Note: make sure a '/' appears at the end of the output folder.
* -f, --file-format, Format of the output file (required).  **.npy:** store data in numpy format. **.h5:** store data in h5py format.
* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Examples**

1. Store indexes of common markers between the query and the reference in .npy format:

```
$ python3 MergeGenome.py store-common-indexes -q query.vcf -r reference.vcf -o ./output/ -f .npy
```

2. Store indexes of common markers between the query and the reference in .h5 format:

```
$ python3 MergeGenome.py store-common-indexes -q query.vcf -r reference.vcf -o ./output/ -f .h5
```
