## Impute

One way to avoid missing all non-matching SNPs when merging genomic datasets is through imputation (e.g. of low-resolution genotype arrays to induce them to have all SNPs from a higher-resolution panel). The potential appplications of genotype imputation are several, from boosting the power of phenotype prediction models and genomewide association studies [1] to allowing comparisons between studies.

There are several genotype imputation softwares to induce shorter DNA sequences to have all the variants from the longer DNA sequences. This way, all relevant information stored in the high-dimensional data is preserved. 

MergeGenome describes how to use Beagle v5.1 [2] to impute the low-resolution dataset (query or target panel) to the high-resolution dataset (reference panel). This way, all sequences are forced to be of the same length, specifically, of the genotype arrays in the reference after cleaning. For this task, it is recommended to provide Beagle with the genetic map, as it usually improves imputation performance.

## Usage - Beagle

```
java -Xmx50g -jar <beagle.jar> gt=<query_file> ref=<reference_file> out=<output_fike> impute=True chrom=<chrom_number>  map=<beagle_map_file>
```

* gt=FILE, Input VCF file with the genotype data to be imputed.
* ref=FILE, Reference panel in bref3 or VCF format. Each genotype must have two phased, non-missing alleles. If a VCF file is specified, the phased allele separator must be used ‘|’.
* out=FILE, Output filename prefix. The prefix may be an absolute or relative filename, but it cannot be a directory name.
* impute=bool, Specifies whether markers that are present in the reference panel but absent in that target will be imputed (default: impute=true). This option has no effect if no reference (ref) panel is specified.
* chrom=[chrom]:[start]-[end], Specifies a chromosome interval:[chrom] is the CHROM field in the input VCF file and [start] and [end] are the starting and ending positions.
* map, Specifies a PLINK format genetic map with cM units.

**Output files**

* **logfile:** gives a summary of the analysis that includes the Beagle version, the command line arguments, and compute time.
 
* **vcf.gz:** bgzip-compressed VCF file that contains phased, non-missing genotypes for all non-reference samples. The output vcf.gz file can be uncompressed with the unix gunzip utility.

For more information about Beagle, read the [[documentation]](https://faculty.washington.edu/browning/beagle/beagle_5.3_07Feb22.pdf).

**Examples**

1. Impute target panel given the reference panel, for chromosome files from 1 too 38:

```
for chr in {1..38};
do

java -Xmx50g -jar /home/packages/beagle.jar gt=/home/query_chr${chr}_cleaned.vcf
ref=/home/reference_chr${chr}_cleaned.vcf out=/home/query_chr${chr}_cleaned_imputed impute=True chrom=chr${chr}
map=/home/miriambt/data/genetic_map/chr${chr}_beagle_gmap.txt

done
```

[1] William YS Wang, Bryan J Barratt, David G Clayton, and John A Todd. Genome-wide association studies: theoretical and practical concerns. Nature Reviews Genetics, 6(2):109–118, 2005.

[2] Brian L Browning and Sharon R Browning. Geno- type imputation with millions of reference samples. The American Journal of Human Genetics, 98(1):116– 126, 2016.
