## Plot SNP means (comparison)

There are several straightforward techniques to qualitatively evaluate the merged dataset. One of them is to compare the SNP means between samples coming from one source or the other. If the data is homogeneous, the SNP means should be similar. Namely, the SNP means from the query should fall within a small bandwidth from the SNP means from the reference. Under any other circumstances, the data might contain inconsistent SNPs. There are several possible causes of diverging SNPs, such as a bad imputation or a sequencing error.

MergeGenome plot-snp-means command produces a scatter plot showing the SNP means for all the common markers between the provided query and reference .vcf files. The provided .vcf files can contain data for a single or multiple chromosomes. Either way, the chromosomes need to be appear in the same order.

Figure 1 shows an output example of MergeGenome plot-snp-means command. We can observe numerous diverging SNPs close to 0.0 and 1.0 in the vertical axis. There are also a few diverging SNPs on the top-right. Probably, the former come from a bad imputation in which the most frequent allele is overrepresented, while the latter come from a sequencing error.

![Figure 1. SNP means comparison](https://github.com/AI-sandbox/merge-vcf-files/blob/main/figures/snp_means_reference_and_query.png)

## Usage

```
$ python3 MergeGenome.py plot-snp-means -q <query_file_1>...<query_file_n> -r <reference_file_1>...<reference_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Paths to query .vcf files with data for a single or multiple chromosomes each (required).
* -r, --reference LIST, Paths to reference .vcf files with data for a single or multiple chromosomes each (required).
* -o, --output-folder PATH, Path to output folder. (required). Note: make sure a '/' appears at the end of the output folder.
* -x, --x-axis-name STR, Name given to reference dataset (x-axis) (optional). Default=reference.
* -y, --y-axis-name STR, Name given to query dataset (y-axis) (optional). Default=query.
* -f, --fontsize INT, Fontsize of all text in plot (optional). Default=25.
* -w, --figure-width INT, Figure width of plot (optional). Default=26.
* -i, --figure-height INT, Figure height of plot (optional). Default=15.
* -s, --size-points INT, Size of the points in plot (optional). Default=0.1.
* -c, --color-points STR, Color of points in plot (optional). Default=#306998.
* -d, --debug PATH, Path to .log/.txt file to store info/debug messages (optional).

**Output**

* A .png image with the SNP means of all the common markers between the reference and the query datasets. The .png image will have the name snp_means_{x_axis_name}_{y_axis_name}.pn
* If --debug, a .log or .txt file with information regarding the dimensions of the data (number of samples and number of SNPs), the chromosomes in each file, and the amount of common markers found.

`Examples`

1. Plot the SNP means between the common markers in the first chromosome in the query and reference datasets:

```
$ python3 MergeGenome.py plot-snp-means -q query_chr1.vcf -r reference_chr1.vcf -o ./output/
```

2. Plot the SNP means between the common markers of all the chromosomes in the query and reference datasets. Change plot configurations and save debug info in log file:

```
$ python3 MergeGenome.py plot-snp-means -q query_chr*.vcf -r reference_chr*.vcf -o ./output/ -q "Imputed Query" -r "Reference" -f 15 -w 16 -i 9 -d
```
