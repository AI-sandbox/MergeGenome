## Rename chromosome notation

A prerequisite to preprocess and merge VCF files with data coming from the same chromosome is to ensure the chromosome nomenclature in the variants/CHROM field is in the exact same format between files. Usually, the chromosome nomenclature is in one of two forms: *<chrom_number>* or *chr<chrom_number>*, depending on the data source.

By default, the rename command from MergeGenome changes the variants/CHROM field nomenclature from *<chrom_number>* to *chr<chrom_number>* or vice-versa. If the notation is in neither of these formats, the rename command will do nothing, except if a specific notation mapping is specified through `--rename-map`. The previous flag can be used to specify a dictionary with the notation mapping, where the keys are the old notations and the values are the new notations. The dictionary should contain a key-value pair for each chromosome notation to be modified. Note that, if the input file contains data for more than one chromosome, the number of key-value pairs in the dictionary should be equal to or smaller than the number of chromosomes in the dataset. In case a mapping is specified for a chromosome not present in the data file, this will be ignored.

## Usage

```
$ python3 MergeGenome.py rename -q <query_file> -o <output_folder>
```

Input flags include:

* -q, --query PATH, Path to input .vcf file with data for a single or multiple chromosomes (required).
* -o, --output-folder PATH, Path to output folder to store partitioned .vcf files (required). Note: make sure a '/' appears at the end of the output folder.
* -m, --rename-map DICT, Mapping from actual to new chromosome notation (optional).
* -d, --debug PATH, Path to .log/.txt file to store info/debug messages (optional).

**Output**

* A .vcf file with renamed chromosome notation. The new .vcf file will receive the same base name as the input file, but ending with '_renamed.vcf'.
* If --debug, a .log or .txt file with information regarding the dimensions of the data (number of samples and number of SNPs), the amount of chromosomes available and the changes in the chromosome notation.

**Examples**

1. Change chromosome notation from *<chrom_number>* to *chr<chrom_number>* (or vice-versa):

```
$ python3 MergeGenome.py rename -q query_allchr.vcf -o ./output/
```

2. Change chromosome notation from "1" to "chr_1", and from "2" to "chr_2". Also save debug info in log file:

```
$ python3 MergeGenome.py rename -q query.vcf -o ./output/ -m '{"1":"chr_1", "2":"chr_2"}' -d rename.log
```