## 6. Concatenate and merge .vcf files

### Concat: Concatenate .vcf files into single .vcf file (same sample IDs)

If multiple .vcf files contain different information about the same samples (i.e., the same sample IDs can be found in all .vcf files) and the desire is to create one .vcf file with information from all the files, the concat function should be used.

Another way to think about this function is that it combines files vertically. 

In the example below, embark_chr{i}.vcf contains the data for chromosome i, where i goes from 1 to 38. The output is embark_allchr.vcf file, which contains the data for all the chromosomes from 1 to 38.

```
bcftools concat embark_chr1.vcf embark_chr2.vcf embark_chr3.vcf embark_chr4.vcf embark_chr5.vcf embark_chr6.vcf embark_chr7.vcf embark_chr8.vcf embark_chr9.vcf embark_chr10.vcf embark_chr11.vcf embark_chr12.vcf embark_chr13.vcf embark_chr14.vcf embark_chr15.vcf embark_chr16.vcf embark_chr17.vcf embark_chr18.vcf embark_chr19.vcf embark_chr20.vcf embark_chr21.vcf embark_chr22.vcf embark_chr23.vcf embark_chr24.vcf embark_chr25.vcf embark_chr26.vcf embark_chr27.vcf embark_chr28.vcf embark_chr29.vcf embark_chr30.vcf embark_chr31.vcf embark_chr32.vcf embark_chr33.vcf embark_chr34.vcf embark_chr35.vcf embark_chr36.vcf embark_chr37.vcf embark_chr38.vcf > embark_allchr.vcf
```

It is pretty common that bcftools returns this type of error:

> contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)

To solve it, you might want to add the missing line in the header from the terminal before using the concat functinon:

```
for chr in {1..38}; do  sed -i "3s/^/##contig=<ID=chr${chr}>\n/" embark_chr${chr}.vcf; done
```

### Merge: Merging .vcf files (different sample IDs)

If multiple .vcf files contain information about different samples (i.e., there is no overlap between sample IDs in all of the .vcf files) and the desire is to create one .vcf file with information from all the files, the merge function should be used.

```
bgzip -ci embark_allchr.vcf > embark_allchr.vcf.gz
```

```
rm -f file.list
for F in  embark_allchr.vcf.gz imputed_array_subset_allchr.vcf.gz
do
      bcftools sort -O b -o ${F}.bcf $F
      bcftools index ${F}.bcf
      echo "${F}.bcf" >> file.list
done

bcftools merge  --file-list file.list  -O z -o merged.vcf.gz
```
