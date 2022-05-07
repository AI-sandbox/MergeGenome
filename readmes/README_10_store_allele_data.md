## Store allele data (NPY or H5)

Allele data in calldata/GT is usually stored as a three-dimensional matrix, where the first dimension corresponds to the number of SNPs, the second dimension to the number of samples and the third dimension to the number of strands. As for the latter dimension, most genomic datasets store allele data for the two strands that constitute DNA sequences: the maternal strand and the paternal strand.

For most data analysis, it is often necessary to get rid of the third dimension. To flatten the third dimension, a possibility is to split the maternal and paternal strands into separate samples. This way, the number of samples in the dataset is doubled and both strands can be treated as two independent data samples hereafter. An alternative is to combine the maternal and paternal strands by averaging them. In this case, both strands can be treated as one data sample.

The store-allele command from MergeGenome supports both splitting and averaging maternal and paternal strands of any given VCF file, and it stores the resulting allele data in a NPY or H5 file.

## Usage

```
$ python3 MergeGenome.py store-allele -q <query_file> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Path to query .vcf files with data for each chromosome (required).
* -o, --output-folder PATH, Path to output folder to store the modified VCF files (required). Note: make sure a '/' appears at the end of the output folder.
* -s, --data-format STRING, Separate or average maternal and paternal strands (required).
    
&emsp;&ensp; separated: split maternal and paternal strands into separate samples. 
    
&emsp;&ensp; averaged: combine maternal and  paternal strands by averaging them.

* -f, --file-format, Format of the output file (required).
    
&emsp;&ensp; .npy: store data in numpy format.
    
&emsp;&ensp; .h5: store data in h5py format.

* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Examples**

1. Store allele data in .npy format by splitting the maternal and paternal strands:

```
$ python3 MergeGenome.py store-allele -q query.vcf -o ./output/ -s separated -f .npy
```

2. Store allele data in .h5 format by averaging the maternal and paternal strands:

```
$ python3 MergeGenome.py store-allele -q query.vcf -o ./output/ -s averaged -f .h5
```
