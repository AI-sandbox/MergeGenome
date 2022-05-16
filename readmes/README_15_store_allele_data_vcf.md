## Store allele data (from NPY to VCF)

TODO

## Usage

```
$ python3 MergeGenome.py store-vcf -q <query_vcf_file> --n <query_npy_file> o <output_folder>
```

Input flags include:

* -q, --vcf-query Path, Path to query .vcf file (required).
* -n, --npy-query Path, Path to query .npy file (required).
* -o, --output-folder PATH, Path to output folder to store the modified VCF files (required). Note: make sure a '/' appears at the end of the output folder.
* -b, --binary-indexes PATH, Path to binary indexes, where 1 means the SNP is correct and 0 that it is to be removed (optional).
* -d, --debug PATH, Path to file to store info/debug messages (optional).
