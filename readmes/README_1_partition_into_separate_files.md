## Partition data into separate VCF files (one per chromosome)

The first preprocessing step is to partition data into manageable independent VCF files (one file per chromosome) in an effort to fasten data retrieval, smooth the management of raw data, and exploit parallel execution. Otherwise, accessing and manipulating large-scale multidimensional biological data might be excessively costly both in terms of execution time and memory consumption.

Before writing the resulting VCF files, the chromosome notation in the variants/CHROM field can optionally be modified. The chromosome notation of the reference and the query datasets should be the same prior to merging. Standardizing the chromosome notation at this step is preferable to doing it in another step so as to fasten the preprocessing.

## Usage

```
$ python3 MergeGenome.py partition -f <file_path> -o <output_folder>
```

Input flags include:

* -f, --file PATH, Path to input VCF file (required).
* -o, --output-folder PATH, Path to output folder to store the separate VCF files (required).
* -r, --rename-chr, Rename chromosome notation (optional).
* -d, --debug PATH, Path to file to store info/debug messages (optional).