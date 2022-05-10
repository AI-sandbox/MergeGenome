## Impute

One way to avoid missing all non-matching SNPs when merging genomic datasets is through imputation, usually, of low-resolution genotype arrays to induce them to have all SNPs from a higher-resolution panel. The potential applications of genotype imputation are several, from boosting the power of phenotype prediction models and genomewide association studies to allowing comparisons between studies.

There are several genotype imputation softwares to induce shorter DNA sequences to have all the variants from the longer DNA sequences. This way, all relevant information in the high-dimensional data is preserved. MergeGenome describes how to use Beagle v5.1 to impute the low-resolution dataset - query or target panel - to the high-resolution dataset - reference panel -. This way, all sequences are forced to be of the same length, specifically, of the genotype arrays in the reference.

To improve the imputation performance, it is recommended to provide the imputation software with the genetic map.

## Usage

To use Beagle:

```
java -Xmx50g -jar <beagle.jar> gt=<query_file> ref=<reference_file> out=<output_fike> impute=True chrom=<chrom_number>  map=<beagle_map_file>
```

Input flags include:

* gt: input VCF file with the genotype data to be imputed.

* ref: reference panel in bref3 or VCF format. Each genotype must have two phased, non-missing alleles. If a VCF file is specified, the phased allele separator must be used ‘|’.

* out: output filename prefix. The prefix may be an absolute or relative filename, but it cannot be a directory name.

* impute: specifies whether markers that are present in the reference panel but absent in that target will be imputed (default: impute=true). This option has no effect if no reference (ref) panel is specified.

* chrom: argument [chrom]:[start]-[end] specifies a chromosome interval:[chrom] is the CHROM field in the input VCF file and [start] and [end] are the starting and ending positions.

* map: specifies a PLINK format genetic map with cM units.

**Output**

* logfile: gives a summary of the analysis that includes the Beagle version, the command line arguments, and compute time.

* vcf.gz: bgzip-compressed VCF file that contains phased, non-missing genotypes for all non-reference samples. The output vcf.gz file can be uncompressed with the unix gunzip utility.

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
