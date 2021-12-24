## 4. Impute SNPs (from short to large dataset)

You can impute the missing variables of the shorter dataset to have the same SNPs as the larger dataset as follows:

```
for chr in {1..38};
do

java -Xmx50g -jar /home/packages/beagle.jar gt=/home/miriambt/data/array/cleaned/All_Pure_150k_chr${chr}_cleaned.vcf ref=/home/miriambt/data/whole_genome/cleaned/chr${chr}_unfiltered_phased_cleaned.vcf out=/home/miriambt/data/array/imputed/All_Pure_150k_chr${chr}_cleaned_imputed impute=True chrom=chr${chr}  map=/home/miriambt/data/genetic_map/chr${chr}_beagle_gmap.txt

done
```

* **gt:** input VCF file with the genotype data to be imputed.

* **ref:**  reference panel in bref3 or VCF format. Each genotype must have two phased, non-missing alleles. If a VCF file is specified, the phased allele separator must be used ‘|’.

* **out:** output filename prefix. The prefix may be an absolute or relative filename, but it cannot be a directory name.

* **impute:** specifies whether markers that are present in the reference panel but absent in that target will be imputed (default: impute=true). This option has no effect if no reference (ref) panel is specified.

* **chrom:** argument [chrom]:[start]-[end] specifies a chromosome interval:[chrom] is the CHROM field in the input VCF file and [start] and [end] are the starting and ending positions.

* **map:** specifies a PLINK format genetic map with cM units. You can use the script in scripts/format_genetic_map_beagle.py to obtain the genetic map in the correct format. Note: it is very important to specify the genetic map.

#### Output files

There are two output files: 

* **logfile:** gives a summary of the analysis that includes the Beagle version, the command line arguments, and compute time. 
* **vcf.gz:** bgzip-compressed VCF file that contains phased, non-missing genotypes for all non-reference samples. The output vcf.gz file can be uncompressed with the unix gunzip utility.
