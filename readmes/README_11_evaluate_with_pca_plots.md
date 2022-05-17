## Plot Principal Component Analysis (PCA)

The Principal Component Analysis (PCA) dimensionality reduction technique can be used in bioinformatics to identify outlier DNA sequences (samples) in any given dataset. Outlier samples can exist for several reasons, including natural genetic differences. For instance, samples from a village or mixed dog might be consistently different from purebred dogs because of their different biological origins. Also, it could be that improbable (not plausible) data forms an outlier. Either way, outlier samples can mislead any analysis conducted on top of it. The most practical solution is to remove them.

PCA can also be used as a tool to visualize how samples are distributed. Usually, clusters for the different ancestries or species present in the data are generated. For humans, the population structure is commonly correlated with geography. Since the PCA captures the variance in the data, this technique is extremely useful to qualitatively determine if joined samples from different sources constitute a homogeneous dataset. Of course, it does not prove the quality of the merged dataset. It only confirms that the first two principal components are not enough to capture any significant variation between samples from one source or the other.

MergeGenome plot-pca command produces a scatter plot with the two first PCA components. The provided .vcf files can contain data for a single or multiple chromosomes. Either way, the chromosomes need to appear in the same order. If only data from the query is provided, the PCA is trained and projected on the query SNPs. If both data from the query and the reference are provided, the PCA is trained and projected on the common markers (i.e., SNPs at the same CHROM, POS, REF, and ALT) between both datasets, except when `--train_query` is called, then the data is only trained on the query and projected on both datasets.

Figure 2.a) shows an output example of MergeGenome plot-pca command, in which the PCA has been trained on the common markers from the query and projected on both the query and the reference datasets. We can observe that the PCA points from the query do not fall into the same space as the PCA points from the reference, suggesting the data have different distributions between sources. However, after preprocessing with MergeGenome, the differences in the DNA sequences between both sources is no longer noticaeble, as shown in Figure 3. This is just an example of why it is necessary to properly preprocess genomic sequences prior to merging.

a) No preprocessing           |  b) MergeGenome Preprocessing
:-------------------------:|:-------------------------:
![](https://github.com/AI-sandbox/merge-vcf-files/blob/main/figures/trained_both_projected_both.png)  |  ![](https://github.com/AI-sandbox/merge-vcf-files/blob/main/figures/trained_both_projected_both_after_preprocessing.png)

*Figure 2*. PCA w/wo MergeGenome preprocessing.

## Usage

```
$ python3 MergeGenome.py plot-pca -q <query_file_1>...<query_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Paths to query .vcf files with data for a single or multiple chromosomes each (optional).
* -r, --reference LIST, Paths to reference .vcf files with data for a single or multiple chromosomes each (optional).
* -o, --output-folder PATH, Path to output folder. (required). Note: make sure a '/' appears at the end of the output folder.
* -t, --train-query, To train the PCA the PCA only on the query instead of on both datasets (optional). Default=False.
* -f, --fontsize INT, Fontsize of all text in plot (optional). Default=25.
* -w, --figure-width INT, Figure width of plot (optional). Default=26.
* -i, --figure-height INT, Figure height of plot (optional). Default=15.
* -s, --size-points INT, Size of points in plot (optional). Default=15.
* -a, --alpha, FLOAT, Transparency of points in plot (optional). Default=0.7.
* -cq, --color-points-query STR, Color of query points in the plot (optional). Default=#EBD0A1.
* -cr, --color-points-reference STR, Color of reference points in the plot (optional). Default=#259988.
* -d, --debug PATH, Path to .log/.txt file to store info/debug messages (optional).

**Output**

* A .png image with the first PCA components. The .png image will have the name 'trained_query' if only the query is provided, 'trained_both_projected_both' if both the query and the reference are provided and 'trained_query_projected_both.png' if both datasets are provided and `--train_query` flag is True.
* If --debug, a .log or .txt file with information regarding the dimensions of the data (number of samples and number of SNPs), the chromosomes in each file, and the amount of common markers found.

`Examples`

1. Plot PCA on query data:

```
$ python3 MergeGenome.py plot-pca -q query_chr1.vcf -o ./output/
```

2. Plot PCA on the common markers between the query and the reference, training the PCA only on the query:

```
$ python3 MergeGenome.py plot-pca -q query_chr1.vcf -r reference_chr1.vcf -o ./output/ -t
```

3. Plot PCA on the common markers between the query and the reference, traing the PCA on both datasets:

```
$ python3 MergeGenome.py plot-pca -q query_chr1.vcf -r reference_chr1.vcf -o ./output/
```
