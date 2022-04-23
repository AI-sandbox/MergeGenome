## Concat VCF files

Finally, the last steps include usiing concat and merge commands from the software BCFtools to combine the files vertically and horizontally, respectively. The result of this process is a single VCF file with the merged samples from the query and the reference for all chromosomes.

If multiple .vcf files contain different information about the same samples (i.e., the same sample IDs can be found in all .vcf files) and the desire is to create one .vcf file with information from all the files, the concat function from bcftools can be used. Another way to think about this function is that it combines files vertically.

## Usage

To use concat from bcftools:

```
$ bcftools concat <file_1>...<file_n>
```

**Example**

```
$ bcftools concat query_chr1.vcf query_chr2.vcf query_chr3.vcf query_chr4.vcf query_chr5.vcf query_chr6.vcf query_chr7.vcf query_chr8.vcf query_chr9.vcf query_chr10.vcf query_chr11.vcf query_chr12.vcf query_chr13.vcf query_chr14.vcf query_chr15.vcf query_chr16.vcf query_chr17.vcf query_chr18.vcf query_chr19.vcf query_chr20.vcf query_chr21.vcf query_chr22.vcf query_chr23.vcf query_chr24.vcf query_chr25.vcf query_chr26.vcf query_chr27.vcf query_chr28.vcf query_chr29.vcf query_chr30.vcf query_chr31.vcf query_chr32.vcf query_chr33.vcf query_chr34.vcf query_chr35.vcf query_chr36.vcf query_chr37.vcf query_chr38.vcf > file_allchr.vcf
```

**Header error/solution**

It is quite common that bcftools returns the error: *contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)*. Thus, before using the concat functinon, you can solve this issue by adding the missing line in the header from the CLI:

```
$ sed -i "3s/^/##contig=<ID=<<chrom_number>>\n/" <file_1>; done
```

To easily apply this command to all chromosome files, you can iterative with a loop as follows:

```
$ for chr in {1..38}; do  sed -i "3s/^/##contig=<ID=chr${chr}>\n/" query_chr${chr}.vcf; done
```