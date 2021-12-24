## 1. Concatenate and merge .vcf files

#### Concatenate .vcf files into single file

You can use this command to concatenate all .vcf files containing data for different chromosomes. In the example below, embark_chr{i}.vcf contains the data for chromosome i, where i goes from 1 to 38. The output in embark_allchr.vcf contains the data for all the chromosomes.

```python
bcftools concat embark_chr1.vcf embark_chr2.vcf embark_chr3.vcf embark_chr4.vcf embark_chr5.vcf embark_chr6.vcf embark_chr7.vcf embark_chr8.vcf embark_chr9.vcf embark_chr10.vcf embark_chr11.vcf embark_chr12.vcf embark_chr13.vcf embark_chr14.vcf embark_chr15.vcf embark_chr16.vcf embark_chr17.vcf embark_chr18.vcf embark_chr19.vcf embark_chr20.vcf embark_chr21.vcf embark_chr22.vcf embark_chr23.vcf embark_chr24.vcf embark_chr25.vcf embark_chr26.vcf embark_chr27.vcf embark_chr28.vcf embark_chr29.vcf embark_chr30.vcf embark_chr31.vcf embark_chr32.vcf embark_chr33.vcf embark_chr34.vcf embark_chr35.vcf embark_chr36.vcf embark_chr37.vcf embark_chr38.vcf > embark_allchr.vcf
```

If bcftools returns this error:

contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)


for chr in {1..38}; do  sed -i "3s/^/##contig=<ID=chr${chr}>\n/" embark_chr${chr}.vcf; done

The computational time of preprocessing big .vcf files can be substantial. That is why, working with one .vcf file for each chromosome is in general better than 
