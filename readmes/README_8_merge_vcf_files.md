## Merge VCF files

The horizontal combination through BCFTools `merge` consists in combining multiple VCF files containing the same information about different samples (i.e. there is no overlap between sample IDs in all of the VCF files). Thus, the `merge` command can be used to obtain a single VCF file with the merged samples from the query and the reference.

It is the user's responsibility to ensure that the sample names are consistent across all files. Unless the option `—-force-samples` is specified, the application will exit with an error. It is worth noting that only records from different files can be merged, never from the same file. Look into BCFtools `concat` for vertical merging.

## Usage

To use `merge` from bcftools:

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

For more information, you can read the [documentation](https://samtools.github.io/bcftools/bcftools.html#merge).

`Example`

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
