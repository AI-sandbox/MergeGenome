## Partition data into separate VCF files (one per chromosome)

The first preprocessing step is to partition the data into manageable independent VCF files (one file per chromosome) in an effort to fasten data retrieval, smooth the management of raw data, and exploit parallel execution. Otherwise, accessing and manipulating large-scale multidimensional biological data might be excessively costly, both in terms of execution time and memory consumption.

Before writing the resulting VCF files, the chromosome notation in the variants/CHROM field can optionally be modified. The chromosome notation of the reference and the query datasets should be the same prior to merging. Standardizing the chromosome notation at this step is preferable to doing it in an added step for efficiency purposes. To rename the chromosome notation, use the `--rename-chr` flag. If the `--rename-chr` flag is called, the variants/CHROM format will change from *<chrom_number>* to *chr<chrom_number>* or vice-versa. If the notation of a chromosome is in neither of these formats, its notation will remain the same. To specify a different notation, use the `--rename-chr` flag to specify a notation mapping for all the chromosomes to modify, where the keys are the old notations and the values are the new notations.

## Usage

```
$ python3 MergeGenome.py partition -f <file_path> -o <output_folder>
```

Input flags include:

* -f, --file PATH, Path to input VCF file (required).
* -o, --output-folder PATH, Path to output folder to store the separate VCF files (required).
* -r, --rename-chr, Rename chromosome notation (optional).
* -m, --rename-map, Dictionary with mapping from actual to new chromosome notation. (optional)'
* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Examples**

1. Partition without renaming chromosome notation:

```
$ python3 MergeGenome.py partition -f array.vcf -o ./array_partitioned/
```

2. Partition and change chromosome notation from *<chrom_number>* to *chr<chrom_number>* (or vice-versa):

```
$ python3 MergeGenome.py partition -f array.vcf -o ./array_partitioned/ -r
```

2. Partition and change chromosome notation from "1" to "chr_1", and from "2" to "chr_2". Also save debug info in log file:

```
$ python3 MergeGenome.py partition -f array.vcf -o ./array_partitioned/ -r -m '{"1":"chr_1", "2":"chr_2"}' -d partition.log
```