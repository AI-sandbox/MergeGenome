## Plot Principal Component Analysis (PCA)

The Principal Component Analysis (PCA) dimensionality reduction technique can be used in bioinformatics to identify outlier samples, where a sample is a DNA squence. For instance, when dealing with dog DNA sequences, samples from village or mixed dogs might be consistently different than from purebred dogs. Also, a sample might be very unsimilar to the rest because it sequencing went wrong. Either way, outlier samples might mislead any future analysis, so it is important to identify them and treat them properly.

Two-dimensional PCA plots can also be used in genomics to to see how the samples are distributed, and evaluate the quality of a merged dataset. If the the query and reference points fall into the same space, it suggests that the merged dataset is homogeneous. At least, two components do not capture any significant variation between samples from one source or the other.

The plot-pca command from MergeGenome plots a scatter plot after PCA. If only data from the query is provided, the PCA is trained and projected on all the SNPs data provided. If both data from the query and the reference are provided, the PCA is trained and projected on the common markers between both datasets. When train_both = False, then the data is trained on the query and projected on both the query and the reference.

## Usage

```
$ python3 MergeGenome.py plot-pca -q -q <query_file_1>...<query_file_n> -r <reference_file_1>...<reference_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Path to query .vcf files with data for each chromosome (required).
* -r, --reference LIST, Paths to reference .vcf files with data for each chromosome (required).
* -o, --output-folder PATH, Path to output folder to store the modified VCF files (required). Note: make sure a '/' appears at the end of the output folder.
* -t, --train-both, Train on both the query and the reference (optional). Default=False.
* -f, --fontsize INT, Fontsize in the plot (optional). Default=25.
* -w, --figure-width INT, Figure width of the plot (optional). Default=26.
* -i, --figure-height INT, Figure height of the plot (optional). Default=15.
* -s, --size-points INT, Size of the points in the plot (optional). Default=0.1.
* -cq, --color-points-query STR, Color of query points in the plot (optional). Default=#259988.
* -cr, --color-points-reference STR, Color of reference points in the plot (optional). Default=#EBD0A1.
* -d, --debug PATH, Path to file to store info/debug messages (optional).

`Examples`

1. Plot PCA on query data:

```
$ python3 MergeGenome.py plot-pca -q query_chr1.vcf -o ./output/
```

2. Plot PCA on the common markers between the query and the reference, trained on the query and projected on both datasets:

```
$ python3 MergeGenome.py plot-pca -q query_chr1.vcf -r reference_chr1.vcf -o ./output/
```

3. Plot PCA on the common markers between the query and the reference, trained and projected on both datasets:

```
$ python3 MergeGenome.py plot-pca -q query_chr1.vcf -r reference_chr1.vcf -o ./output/ -t
```
