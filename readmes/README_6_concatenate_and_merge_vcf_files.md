## Concat/Merge VCF files

### Concat: Concatenate VCF files into single VCF file (same sample IDs)

If multiple .vcf files contain different information about the same samples (i.e., the same sample IDs can be found in all .vcf files) and the desire is to create one .vcf file with information from all the files, the concat function from bcftools should be used. Another way to think about this function is that it combines files vertically.

## Usage

To use concat from bcftools:

```
bcftools concat <file_1>...<file_n>
```

**Example**

```
bcftools concat file_chr1.vcf file_chr2.vcf file_chr3.vcf file_chr4.vcf file_chr5.vcf file_chr6.vcf file_chr7.vcf file_chr8.vcf file_chr9.vcf file_chr10.vcf file_chr11.vcf file_chr12.vcf file_chr13.vcf file_chr14.vcf file_chr15.vcf file_chr16.vcf file_chr17.vcf file_chr18.vcf file_chr19.vcf file_chr20.vcf file_chr21.vcf file_chr22.vcf file_chr23.vcf file_chr24.vcf file_chr25.vcf file_chr26.vcf file_chr27.vcf file_chr28.vcf file_chr29.vcf file_chr30.vcf file_chr31.vcf file_chr32.vcf file_chr33.vcf file_chr34.vcf file_chr35.vcf file_chr36.vcf file_chr37.vcf file_chr38.vcf > file_allchr.vcf
```

**Header error/solution**

It is pretty common that bcftools returns the error: *contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)* You can add the missing line in the header from the terminal before using the concat functinon to solve it:

```
for chr in {1..38}; do  sed -i "3s/^/##contig=<ID=chr${chr}>\n/" file_chr${chr}.vcf; done
```

### Merge: Merging VCF files (different sample IDs)

If multiple .vcf files contain information about different samples (i.e., there is no overlap between sample IDs in all of the .vcf files) and the desire is to create one .vcf file with information from all the files, the merge function from bcftools should be used.

## Usage

To use merge from bcftools:

```
bgzip -ci <file_all_chr1> > <file_all_chr1>.gz
bgzip -ci <file_allchr2> > <file_allchr2>.gz

rm -f file.list
for F in  <file_all_chr1>.gz <file_allchr2>.gz
do
      bcftools sort -O b -o ${F}.bcf $F
      bcftools index ${F}.bcf
      echo "${F}.bcf" >> file.list
done

bcftools merge  --file-list file.list  -O z -o <merged_file>
```

**Example**

```
bgzip -ci file_allchr.vcf > file_allchr.vcf.gz
bgzip -ci imputed_array_subset_allchr.vcf > imputed_array_subset_allchr.vcf.gz

rm -f file.list
for F in  file_allchr.vcf.gz imputed_array_subset_allchr.vcf.gz
do
      bcftools sort -O b -o ${F}.bcf $F
      bcftools index ${F}.bcf
      echo "${F}.bcf" >> file.list
done

bcftools merge  --file-list file.list  -O z -o merged.vcf.gz
```
