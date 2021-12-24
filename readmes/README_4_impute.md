## 4. Impute SNPs (from short to large dataset)

You can impute the missing variables of the shorter dataset to have the same SNPs as the larger dataset as follows:

````{python}
for chr in {1..38};
do

java -Xmx50g -jar /home/packages/beagle.jar gt=/home/miriambt/data/array/cleaned/All_Pure_150k_chr${chr}_cleaned.vcf ref=/home/miriambt/data/whole_genome/cleaned/chr${chr}_unfiltered_phased_cleaned.vcf out=/home/miriambt/data/array/imputed/All_Pure_150k_chr${chr}_cleaned_imputed impute=True chrom=chr${chr}  map=/home/miriambt/data/genetic_map/chr${chr}_beagle_gmap.txt

done
````

* **gt:** VCF file containing the genomic sequences to be imputed.

* **ref:**  reference panel in bref3 or VCF format. Each genotype must have two phased, non-missing alleles.

* **out:** output filename prefix. The prefix may be an absolute or relative filename, but it cannot be a directory name.

* **impute:** specifies whether markers that are present in the reference panel but absent in that target will be imputed (default: impute=true). This option has no effect if no reference (ref) panel is specified.

* **chrom:** specifies a chromosome interval.

* **map:** specifies a PLINK format genetic map with cM units. You can use the script in scripts/format_genetic_map_beagle.py to obtain the genetic map in the correct format.
