## Store allele data

Allele data in calldata/GT is usually stored as a three-dimensional matrix, where the first dimension corresponds to the number of SNPs, the second to the number of samples and the third to the number of strands. As for the latter dimension, most genomic datasets store allele data for the two strands available in DNA sequences: the maternal and the paternal strands.

For most data analysis, it is often necessary to get rid of the third dimension. To flatten the third dimension, a possibility is to split the maternal and paternal strands into separate samples. Thus, the number of samples in the dataset is doubled. Both strands can be treated as two independent data samples hereafter. An alternative is to combine the maternal and paternal strands by averaging them. Conversely, both strands can be treated as one data sample in this case.

## Usage

```
$ python3 MergeGenome.py store-allele -q <query_file> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Path to query .vcf files with data for each chromosome (required).
* -o, --output-folder PATH, Path to output folder to store the modified VCF files (required).
* -s, --data-format STRING, Separate or average maternal and paternal strands (required).
    separated: split maternal and paternal strands into separate samples. 
    averaged: combine maternal and  paternal strands by averaging them.
* -f, --file-format, Format of the output file (required).
    .npy: store data in numpy format.
    .h5: store data in h5py format.
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
