## Concat VCF files

The `concat` and `merge` commands from the BCFtools software can be used to combine VCF files vertically and horizontally, respectively. The vertical combination through `concat` consists in combining multiple VCF files containing different information about the same samples (i.e., the same sample IDs can be found in all VCF files). All source files must have the same sample columns appearing in the same order. The input files must be sorted by chromosome and position. While the `partition` command from MergeGenome partitions a single VCF file into a VCF file per chromosome, the BCFtools `concat` command can be used to concatenate the chromosome VCF files back into one VCF.

## Usage

To use `concat` from BCFtools:

```
$ bcftools concat <query_1>...<query_n>
```

For more information, read the [documentation](https://samtools.github.io/bcftools/bcftools.html#concat).

`Example`

```
$ bcftools concat query_chr1.vcf query_chr2.vcf query_chr3.vcf query_chr4.vcf query_chr5.vcf query_chr6.vcf query_chr7.vcf query_chr8.vcf query_chr9.vcf query_chr10.vcf query_chr11.vcf query_chr12.vcf query_chr13.vcf query_chr14.vcf query_chr15.vcf query_chr16.vcf query_chr17.vcf query_chr18.vcf query_chr19.vcf query_chr20.vcf query_chr21.vcf query_chr22.vcf query_chr23.vcf query_chr24.vcf query_chr25.vcf query_chr26.vcf query_chr27.vcf query_chr28.vcf query_chr29.vcf query_chr30.vcf query_chr31.vcf query_chr32.vcf query_chr33.vcf query_chr34.vcf query_chr35.vcf query_chr36.vcf query_chr37.vcf query_chr38.vcf > file_allchr.vcf
```

## Header issues and solution

If the VCF files are missing some mandatory header, BCFtools `concat` command will throw an error. For instance, if a VCF file contains data for chromosome 1, but this information is not present in the header of the file, BCFtools will throw the error:

*contig 'chr1' is not defined in the header. (Quick workaround: index the file with tabix.)*.

To write the missing header at the top of the VCF file, `sed` can be used. For example, to specify that the VCF file contains data for chromosome 1:

```
sed -i '3s/^/##contig=<ID=chr1>\n/' query_chr1.vcf
```

If the header is missing in all VCF files, each file containing data for a particular chromosome, the header line can be added to each file with a loop, as follows:

```
$ for chr in {1..38}; do  sed -i "3s/^/##contig=<ID=chr${chr}>\n/" query_chr${chr}.vcf; done
```

Finally, to check all the header lines present in a VCF file, `head` is a quick workaround:

```
$ head -<lies_to_be_displayed> <query_n>
```