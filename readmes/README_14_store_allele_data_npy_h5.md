## Store allele data (from VCF to NPY or H5)

Allele data in calldata/GT is usually stored as a three-dimensional matrix, where the first dimension corresponds to the number of SNPs, the second dimension to the number of samples and the third dimension to the number of strands. As for the latter dimension, most genomic datasets store allele data for the two strands that constitute DNA sequences, that is, the maternal strand and the paternal strands.

For most data analysis, it is often necessary to get rid of the third dimension. In order to flatten the third dimension, a possibility is to split the maternal and paternal strands into separate samples. This way, the number of samples in the dataset is doubled, and both strands can be treated as two independent data samples hereafter. An alternative is to combine the maternal and paternal strands by averaging them. In this case, both strands can be treated as one data sample.

The MergeGenome store-allele command supports both splitting and averaging maternal and paternal strands of any given VCF file, and it stores the resulting allele data in a NPY or H5 file. By default, the strands are separated and the result is stores in NPY format.

## Usage

```
$ python3 MergeGenome.py store-allele -q <query_file> -o <output_folder>
```

Input flags include:

* -q, --query PATH, Path to query .vcf file (required).
* -o, --output-folder PATH, Path to output folder (required). Note: make sure a '/' appears at the end of the output folder.
* -s, --data-format STR, Separate or average maternal and paternal strands (required). Accepted values: 'separated' to split maternal and paternal strands into separate samples; 'averaged' to combine maternal and  paternal strands by averaging them.
* -f, --file-format STR, Format of the output file (required). Options: '.npy', to store data in numpy format; '.h5' to store data in h5py format. Default='.npy'.
* -d, --debug PATH, Path to .log/.txt file to store info/debug messages (optional).

**Output**

* A .npy or .h5 file with the separated or averaged maternal and paternal strands. The new file will receive the same base name as the input file, but ending with '_{data_format}.{file-format}'.
* If --debug, a .log or .txt file with information regarding the dimensions of the data (number of samples and number of SNPs), and the dimensions after splitting or averaging the maternal and paternal strands.

`Examples`

1. Store allele data in .npy format by splitting the maternal and paternal strands:

```
$ python3 MergeGenome.py store-allele -q query.vcf -o ./output/ -s separated -f .npy
```

2. Store allele data in .h5 format by averaging the maternal and paternal strands:

```
$ python3 MergeGenome.py store-allele -q query.vcf -o ./output/ -s averaged -f .h5
```
