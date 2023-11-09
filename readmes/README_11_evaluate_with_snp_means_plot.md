## Plot SNP means (comparison)

There are several straightforward techniques to qualitatively evaluate the merged dataset. One of them is to compare the SNP means between samples coming from one source or the other. If the data is homogeneous, the SNP means should be similar. Namely, the SNP means from the query should fall within a small bandwidth from the SNP means from the reference. Under any other circumstances, the data might contain inconsistent SNPs. There are several possible causes of diverging SNPs, such as a bad imputation or a sequencing error.

MergeGenome plot-snp-means command produces a scatter plot showing the SNP means for all the common markers (i.e., SNPs at the same CHROM, POS, REF, and ALT) between the provided query and reference .vcf files. The provided .vcf files can contain data for a single or multiple chromosomes. Either way, the chromosomes need to appear in the same order. In case of observing different SNP means between the query and the reference, it is recommennded to run DataFix `detect_and_fix_mismatch` or MergeGenome `remove-snps-different-means` command. While the computationl cost of the former is much higher, it has the benefit of preserving the features by correcting them through imputation. For a more detailed description of the `detect_and_fix_mismatch` task, read the [documentation](readmes/README_6_ml_source_classifiers_for_snp_filtering.md). Likewise, for a more detailed description of the MergeGenome `remove-snps-different-means` command, read the [documentation](https://github.com/AI-sandbox/merge-vcf-files/blob/main/readmes/README_13_remove_snps_different_means.md).

Figure 1 shows an output example of MergeGenome plot-snp-means command. We can observe numerous diverging SNPs close to 0.0 and 1.0 in the vertical axis. There are also a few diverging SNPs on the top-right. Probably, the former come from a bad imputation since the most frequent allele is overrepresented, while the latter come from a sequencing error.

![Figure 1. SNP means comparison](https://github.com/AI-sandbox/merge-vcf-files/blob/main/figures/snp_means_reference_and_query.png)

*Figure 1*. SNP means comparison for the common markers between a query and a reference datasets.

## Usage

```
$ python3 MergeGenome.py plot-snp-means -q <query_file_1>...<query_file_n> -r <reference_file_1>...<reference_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Paths to query .vcf files with data for a single or multiple chromosomes each (required).
* -r, --reference LIST, Paths to reference .vcf files with data for a single or multiple chromosomes each (required).
* -o, --output-folder PATH, Path to output folder (required). Note: make sure a '/' appears at the end of the output folder.
* -idx, --indexes PATH, Path to SNP indexes with different color.
* -x, --x-axis-name STR, Name given to reference dataset (x-axis) (optional). Default="reference".
* -y, --y-axis-name STR, Name given to query dataset (y-axis) (optional). Default="query".
* -f, --fontsize FLOAT, Fontsize of all text in plot (optional). Default=35.
* -w, --figure-width FLOAT, Figure width of plot (optional). Default=20.
* -i, --figure-height FLOAT, Figure height of plot (optional). Default=20.
* -s, --size-points FLOAT, Size of points in plot (optional). Default=0.7.
* -a, --alpha FLOAT, Transparency of points in plot. Default=0.7.
* -c, --color-points STR, Color of all points not in indexes (optional). Default="#FA8460".
* -ci, --color-points-indexes STR, Color of all points in indexes (optional). Default="#7793F5".
* -d, --debug PATH, Path to .log/.txt file to store info/debug messages (optional).
* -l, --legend-points STR, Name given to all points not in indexes (optional). Default="SNPs not in indexes".
* -li, --legend-points-indexes STR, Name given to all points in indexes (optional). Default="SNPs in indexes".

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
