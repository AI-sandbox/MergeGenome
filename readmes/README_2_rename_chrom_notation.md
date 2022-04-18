## 2. Rename chromosome notation

A necessary step prior to analyzing and preprocessing two VCF files with data from the same chromosome is to ensure the chromosome nomenclature in the variants/CHROM field is in the same format. Usually, the chromosome nomenclature is in one of two forms: *<chrom_number>* or *chr<chrom_number>*. By default, the rename command changes the format from *<chrom_number>* to *chr<chrom_number>* or (vice-versa). If the notation is in neither of these formats, it does nothing. To specify a different notation, use the `--rename-chr` flag to specify a notation mapping for all the chromosomes to modify, where the key is the old notation for the chromosome in particular and the value is the new notation.

## Usage

```
$ python3 MergeGenome.py rename -f <file_path> -o <output_folder>
```

Input flags include:

* -f, --file PATH, Path to .vcf file with data for a particular chromosome (required).
* -o, --output-folder PATH, Path to output folder to store the modified VCF file (required).
* -m, --rename-map DICT, Dictionary with mapping from actual to new chromosome notation (optional).
* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Examples**

1. Change chromosome notation from *<chrom_number>* to *chr<chrom_number>* (or vice-versa):

```
$ python3 MergeGenome.py rename -f array.vcf -o ./array_partitioned/
```

2. Change chromosome notation from "1" to "chr_1", and from "2" to "chr_2". Also save debug info in log file:

```
$ python3 MergeGenome.py rename -f array.vcf -o ./array_partitioned/ -m '{"1":"chr_1", "2":"chr_2"}' -d partition.log
```