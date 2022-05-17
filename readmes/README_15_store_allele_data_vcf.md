## Store allele data (from NPY to VCF)

The MergeGenome store-vcf command reads from a provided NPY file with a two-dimensional matrix containing splitted maternal and paternal strands and stores this data in the three-dimensional calldata/GT field of the provided VCF file. Additionaly, `--binary-indexes` flag allows removing all incorrect SNPs with a 0 in the the NPY file with the binary indexes. Therefore, the MergeGenome store-vcf command can be used to integrate the output of other preprocessing softwares (e.g., [DataFix](https://github.com/AI-sandbox/Datafix)) back into a VCF file.

## Usage

```
$ python3 MergeGenome.py store-vcf -q <query_vcf_file> --n <query_npy_file> o <output_folder>
```

Input flags include:

* -q, --vcf-query Path, Path to query .vcf file (required).
* -n, --npy-query Path, Path to query .npy file (required).
* -o, --output-folder PATH, Path to output folder (required). Note: make sure a '/' appears at the end of the output folder.
* -b, --binary-indexes PATH, Path to binary indexes, where 1 means the SNP is correct and 0 that it is to be removed (optional).
* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Output**

* A .vcf file with calldata/GT from the provided .npy file and, optionally, remove undesired SNPs. The new file will receive the same base name as the input .vcf file, but ending with '_from_npy.vcf'.
* If --debug, a .log or .txt file with information regarding the dimensions of the data (number of samples and number of SNPs).

`Examples`

1. Store allele data in .npy format by splitting the maternal and paternal strands:

```
$ python3 MergeGenome.py store-allele -q query.vcf -o ./output/ -s separated -f .npy
```

2. Store allele data in .h5 format by averaging the maternal and paternal strands:

```
$ python3 MergeGenome.py store-allele -q query.vcf -o ./output/ -s averaged -f .h5
```
