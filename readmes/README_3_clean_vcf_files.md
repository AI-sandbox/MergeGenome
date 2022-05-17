## Clean VCF files

It is highly recommended to clean genomic sequences to ensure data robustness and consistency between sources. MergeGenome includes options for the most relevant and common .vcf transformations, including the removal of undesired samples, the detection and removal of ambiguous SNPs, the detection and correction of SNP flips, the detection and removal of SNP mismatches, and finally, the renaming of missing values to make them compatible with other softwares. The removal of undesired samples, the detection and removal of ambiguous SNPs and the renaming of missing values apply to a single .vcf - query or reference - file, while the detection and correction of SNP flips and the detection and removal of SNP mismatches operate on a query and a reference file with data from a particular chromosome.

Below, is a short description and motivation for each of the cleaning steps supported by MergeGenome:

**Remove samples**

In many cases, it is desirable to remove a subset of DNA sequences before proceeding any further. For instance, DNA sequences from unusual breeds or species may add bias to any model built on top of that. MergeGenome `--remove-sample-ID-query` and `--remove-sample-ID-reference` flags remove a list of DNA samples from the respective query and reference .vcf files.

**Remove ambiguous SNPs**

A/T, T/A, C/G, and G/C pairs are strand-ambiguous because their components are complementary, being a challenge to determine the DNA strand. MergeGenome `--remove-ambiguous-snps-query` and `--remove-ambiguous-snps-reference` flags detect and remove ambiguous SNPs from the respective query and reference .vcf files.

**Correct SNP flips** 

An SNP flip happens when the REF (reference) and ALT (alternate) fields are swapped between SNPs at the same chromosome and position of the query and reference .vcf files. MergeGenome `--correct-snp-flips corrects` SNP flips by swapping the variants/REF and the variants/ALT of detected  SNP flips in the query, taking as true values the reference. Consequently, the zeros and ones in the calldata/GT field in the query are also swapped. An example of an SNP flip could be an SNP that in the query has variants/REF='A' and variants/ALT='C', but in the reference has variants/REF='C' and variants/ALT='A'. In this case, the nomenclature in the query would be changed to be equal to that of the reference.

**Remove SNP mismatches**

There is a mismatch between SNPs at the same chromosome and the position of the query and reference .vcf files whenever the variants/REF or variants/ALT fields do not coincide. Such SNPs are erroneous features in, at least, one of the two datasets. MergeGenome `--remove-mismatching-snps` flag detects and removed mismatching SNPs from the query and reference .vcf files. Not correcting the SNP flips before removing the SNP mismatches implies losing all SNPs with a flip. That is why, if `--remove-mismatching-snps` and `--correct-snp-flips` are True, MergeGenome will always correct the SNP flips before removing the SNP mismatches.

**Rename missing values**

If missing allele values in the calldata/GT field are badly encoded (e.g., with a -1), some widely used software packages such as Beagle can give the following error: *"Caused by: java.lang.IllegalArgumentException: ERROR: invalid allele [-1]"*. Changing the missing nomenclature to a dot "." usually solves this issue. MergeGenome `--rename-map-query` and `--rename-map-reference` rename missing values notation through a key-value pair, where the key is the actual missing value nomenclature and new is the new name adopted by any missing value in the calldata/GT field.

## Usage

```
$ python3 MergeGenome.py clean -q <query_file_1>...<query_file_n> -r <reference_file_1>...<reference_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Path to query .vcf files with data for for a single chromosome each (required).
* -r, --reference LIST, Paths to reference .vcf files with data for for a single chromosome each (optional).
* -o, --output-folder PATH, Path to output folder (required). Note: make sure a '/' appears at the end of the output folder.
* -s, --remove-sample-ID-query LIST, Sample IDs or substring of sample IDs to remove from the query (optional).
* -t, --remove-sample-ID-reference LIST, Sample IDs or substring of sample IDs to remove from the reference (optional).
* -a, --remove-ambiguous-snps-query, To remove (or not) ambiguous SNPs from the query (optional).
* -b, --remove-ambiguous-snps-reference, To remove (or not) ambiguous SNPs from the reference (optional).
* -f, --correct-snp-flips, To correct (or not) SNP flips in the query with respect to the reference (optional).
* -m, --remove-mismatching-snps, To remove (or not) mismatching SNPs (optional).
* -v, --rename-map-query, Mapping from old to new missing values notation for the query (optional).
* -w, --rename-map-reference, Mapping from old to new missing values notation for the reference (optional).
* -d, --debug PATH, Path to .log/.txt file to store info/debug messages (optional).

**Output**

* One .vcf file with transformed data for each input query and reference files. Each new .vcf file will receive the same base name as the input file, but ending with '_cleaned.vcf'.
* If --debug, a .log or .txt file with information regarding the dimensions of the data (number of samples and number of SNPs), the chromosome in each file, and the amount of samples or SNPs affected in each of the transformations.

`Examples`

1. Remove genomic sequences by sample ID from the query:

```
$ python3 MergeGenome.py clean -q query_chr1.vcf -o ./output/ -s wolf fox coyote dhole
```

2. Remove ambiguous SNPs from both the query and reference:

```
$ python3 MergeGenome.py clean -q query_chr1.vcf -r reference_chr1.vcf -o ./output/ -a -b
```

3. Correct SNP flips from the query with respect to the reference:

```
$ python3 MergeGenome.py clean -q query_chr1.vcf -r reference_chr1.vcf -o ./output/ -f
```

4. Remove SNP mismatches between the query and the reference:

```
$ python3 MergeGenome.py clean -q query_chr1.vcf -r reference_chr1.vcf -o ./output/ -m
```

5. Rename missing values from -1 to '.' both in the query and the reference:

```
$ python3 MergeGenome.py clean -q query_chr1.vcf -r reference_chr1.vcf -o ./output/ -v '{-1:"."}' -w '{-1:"."}'
```

6. Apply all previous cleaning steps for the query and reference, and for each chromosome:

```
$ python3 MergeGenome.py clean -q query_chr*.vcf -r reference_chr*.vcf -o ./output/ -s wolf fox coyote dhole -a -b -f -m -v '{-1:"."}' -w '{-1:"."}'
```



