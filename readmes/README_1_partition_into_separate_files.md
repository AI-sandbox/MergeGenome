## Partition data into a separate VCF file per chromosome

The first preprocessing step is to partition the data into manageable independent VCF files (one file per chromosome), in an effort to fasten data retrieval, smooth the management of raw data, and exploit parallel execution. Otherwise, accessing and manipulating large-scale multidimensional biological data might be excessively costly, both in terms of execution time and memory consumption.

Before writing the partitioned VCF files, the chromosome notation in the variants/CHROM field can optionally be modified. The chromosome notation of the reference and the query datasets should be the same prior to merging. Standardizing the chromosome notation at this step is preferable to doing it in an added step for efficiency purposes. To rename the chromosome notation, you can use the `--rename-chr` flag. If the `--rename-chr` flag is called, the variants/CHROM format will change from *<chrom_number>* to *chr<chrom_number>* or vice-versa by default. If the notation of a chromosome is in neither of these formats, its notation will remain the same. To specify a different notation mapping, you can use the `--rename-map` flag and define a dictionary where the keys are the old chromosome names and the values are the new chromosome names.

## Usage

```
$ python3 MergeGenome.py partition -q <query_file> -o <output_folder>
```

Input flags include:

* -q, --query PATH, Path to input .vcf file with data for multiple chromosomes (required).
* -o, --output-folder PATH, Path to output folder to store partitioned .vcf files (required). Note: make sure a '/' appears at the end of the output folder.
* -r, --rename-chr, To rename chromosome notation (optional).
* -m, --rename-map DICT, Mapping from actual to new chromosome notation (optional).
* -d, --debug PATH, Path to .log/.txt file to store info/debug messages (optional).

**Output**

* One .vcf file for each chromosome in <query_file>.

* If --debug, a .log or .txt file will be created with information regarding the dimensions of the data (number of samples and number of SNPs), the amount of chromosomes available and their dimension and, when applicable, the changes in the chromosome notation.

**Examples**

1. Partition .vcf data in a separate .vcf file per chromosome:

```
$ python3 MergeGenome.py partition -q query.vcf -o ./output/
```

2. Partition .vcf data in a separate .vcf file per chromosome and change chromosome notation from *<chrom_number>* to *chr<chrom_number>*:

```
$ python3 MergeGenome.py partition -q query.vcf -o ./output/ -r
```

3. Partition .vcf data in a separate .vcf file per chromosome and change chromosome notation from *chr<chrom_number>* to *<chrom_number>*:

```
$ python3 MergeGenome.py partition -q query.vcf -o ./output/ -r
```

4. Partition .vcf data in a separate .vcf file per chromosome and change chromosome notation from "1" to "chr_1", and from "2" to "chr_2". Also save debug info in .log file:

```
$ python3 MergeGenome.py partition -q query.vcf -o ./output/ -r -m '{"1":"chr_1", "2":"chr_2"}' -d partition.log
```