## Evaluate with SNP means plot

One way to evaluate the quality of the merged dataset is to compare the SNP means between the samples coming from the query and the samples coming from the reference. It is common that there are some small inherent differences in the SNP value distributions between two datasets due to them having non-identical populational or breed representations. However, in order to consider the data is homogeneous, the distribution of zeros and ones between all common markers should be similar. Accordingly, if most SNPs have similar means between datasets, this indicates that the samples were properly merged. Contrareously, if multiple SNPs diverge from the line y = x (where x is the mean of a SNP in the reference and y is the mean of the same SNP in the query), this indicates that the genomic data was not properly phased, there are sequencing inconsistencies or the imputation went wrong. It is common that when imputing SNP the most frequent allele is overrepresented.

The plot-means command from MergeGenome plots a scatter plot with the SNP means between all common markers between the given datasets, as long as they contain data for the same chromosomes.

## Usage

```
$ python3 MergeGenome.py plot-snp-means -r <reference_file_1>...<reference_file_n> -q <query_file_1>...<query_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Path to query .vcf files with data for each chromosome (required).
* -r, --reference LIST, Paths to reference .vcf files with data for each chromosome (required).
* -o, --output-folder PATH, Path to output folder to store the modified VCF files (required).
* -x, --x-axis-name STR, Name given to the query dataset that will appear in the x-axis (optional). Default=query.
* -y, --y-axis-name STR, Name given to the query dataset that will appear in the y-axis (optional). Default=reference.
* -f, --fontsize INT, Fontsize in the plot (optional). Default=25.
* -w, --figure-width INT, Figure width of the plot (optional). Default=26.
* -h, --figure-height INT, Figure height of the plot (optional). Default=15.
* -s, --size-points INT, Size of the points in the plot (optional). Default=0.1.
* -c, --color-points INT, Color of the points in the plot (optional). Default=#306998.
* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Examples**

1. Plot the SNP means between the common markers in the first chromosome in the query and reference datasets:

```
$ python3 MergeGenome.py plot-snp-means -q query_chr1.vcf -r reference_chr1.vcf -o ./output/
```

2. Plot the SNP means between the common markers of all the chromosomes in the query and reference datasets, setting the fontsize to 15 and the figure size to (16, 9). Also save debug info in log file:

```
$ python3 MergeGenome.py plot-snp-means -q query_chr*.vcf -r reference_chr*.vcf -o ./output/ -f 15 -w 16 -h 9 -d
```
