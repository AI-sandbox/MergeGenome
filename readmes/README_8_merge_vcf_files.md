## Merge VCF files

If the .vcf files contain information about different samples (i.e., there is no overlap between sample IDs in all of the .vcf files), to create one single .vcf file with information for a particular chromosome, the merge function from bcftools should be used.

The result of this process is a single VCF file with the merged samples from the query and the reference for all chromosomes.

## Usage

To use merge from bcftools:

```
bgzip -ci <query_allchr.vcf> > <query_allchr.vcf.gz>
bgzip -ci <reference_allchr.vcf> > <reference_allchr.vcf.gz>

rm -f file.list
for F in <query_allchr.vcf>.gz <reference_allchr.vcf.gz>
do
      bcftools sort -O b -o ${F}.bcf $F
      bcftools index ${F}.bcf
      echo "${F}.bcf" >> file.list
done

bcftools merge --file-list file.list -O z -o <merged_file>
```

**Example**

```
bgzip -ci query_allchr.vcf > query_allchr.vcf.gz
bgzip -ci reference_allchr.vcf > reference_allchr.vcf.gz

rm -f file.list
for F in query_allchr.vcf.gz reference_allchr.vcf.gz
do
      bcftools sort -O b -o ${F}.bcf $F
      bcftools index ${F}.bcf
      echo "${F}.bcf" >> file.list
done

bcftools merge --file-list file.list -O z -o merged.vcf.gz
```
