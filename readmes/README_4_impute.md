## Impute

One way to avoid missing all non-matching SNPs when merging genomic datasets is through imputation (e.g. of low-resolution genotype arrays to induce them to have all SNPs from a higher-resolution panel). The potential applications of genotype imputation are several, from boosting the power of phenotype prediction models and genomewide association studies to allowing comparisons between studies.

There are several genotype imputation softwares to induce shorter DNA sequences to have all the variants from the longer DNA sequences. This way, all relevant information stored in the high-dimensional data is preserved. MergeGenome describes how to use Beagle v5.1 to impute the low-resolution dataset (query or target panel) to the high-resolution dataset (reference panel). This way, all sequences are forced to be of the same length, specifically, of the genotype arrays in the reference after cleaning.

To improve the imputation performance, it is recommended to provide Beagle with a genetic map.

## Usage

To use Beagle:

```
java -Xmx50g -jar <beagle.jar> gt=<query_file> ref=<reference_file> out=<output_fike> impute=True chrom=<chrom_number>  map=<beagle_map_file>
```

For more information, you can read the [documentation](https://faculty.washington.edu/browning/beagle/beagle_5.3_07Feb22.pdf).

**Example**

1. Impute target panel given the reference panel, for chromosome files from 1 to 38:

```
for chr in {1..38};
do

java -Xmx50g -jar /home/packages/beagle.jar gt=/home/query_chr${chr}_cleaned.vcf
ref=/home/reference_chr${chr}_cleaned.vcf out=/home/query_chr${chr}_cleaned_imputed impute=True chrom=chr${chr}
map=/home/genetic_map/chr${chr}_beagle_gmap.txt

done
```
