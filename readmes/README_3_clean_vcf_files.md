## Clean VCF files

There are multiple preprocessing steps which are important to clean genomic data prior to merging. MergeGenome includes options for solving the most important ones: remove genomic sequences by sample ID, remove ambiguous SNPs, correct SNP flips, remove SNP mismatches, and rename missing values:

**Remove genomic sequences by sample ID**

In many cases, it is desirable to remove some genomic sequences before proceeding any further (e.g. because genomic sequences from unusual breeds or species may add bias when using them for training an imputation model for infering missing allele data). In some other cases, there is a list of specific genomics sequences that just want to be excluded from the analysis.

**Remove ambiguous SNPs**

A/T, T/A, C/G, and G/C pairs are strand-ambiguous because they are components are complementary, and it can be challenging to determine the DNA strand given the lack of consensus on how DNA sequences are read.

**Correct SNP flips** 

An SNP flip happens when the REF (reference) and ALT (alternate) fields are swapped between SNPs at the same position. MergeGenome corrects SNP flips by swapping the REF by the ALT (and vice-versa) in the query, and also the 0's by the 1's (and vice-versa).

**Remove SNP mismatches**

There is a mismatch between SNPs at the same position whenever the REF or ALT fields do not coincide. Such SNPs are erroneous features in, at least, one dataset. Note that if the SNP flips are not corrected before removing SNP mismatches, the former will be lost. Hence, MergeGenome corrects the SNP flips before removing the SNP mismatches.

**Rename missing values**

Some widely used software packages such as Beagle can give the following error: *"Caused by: java.lang.IllegalArgumentException: ERROR: invalid allele [-1]"* if missing allele values are encoded as -1. Changing the missing nomenclature to a dot "." should solve this issue. MergeGenome allows renaming missing values notation from the actual value (key) to the new value (value) with a dictionary.

## Usage

```
$ python3 MergeGenome.py clean -r <reference_file_1>...<reference_file_n> -q <query_file_1>...<query_file_n> -o <output_folder>
```

Input flags include:

* -q, --query LIST, Path to query .vcf files with data for each chromosome (required).
* -r, --reference LIST, Paths to reference .vcf files with data for each chromosome (optional).
* -o, --output-folder PATH, Path to output folder to store the modified VCF files (required).
* -s, --remove-sample-ID-query LIST, Sample IDs or substring of sample IDs to remove from the query (optional).
* -s, --remove-sample-ID-reference LIST, Sample IDs or substring of sample IDs to remove from the reference (optional).
* -a, --remove-ambiguous-snps-query, Remove (or not) ambiguous SNPs from the query (optional).
* -b, --remove-ambiguous-snps-reference, Remove (or not) ambiguous SNPs from the reference (optional).
* -f, --correct-snp-flips, Correct (or not) SNP in the query with respect to the reference (optional).
* -m, --remove-mismatching-snps, Remove (or not) mismatching SNPs (optional).
* -v, --rename-map-query, Dictionary with mapping from old to new missing notation in the reference (optional).
* -w, --rename-map-reference, Dictionary with mapping from old to new missing notation in the query (optional).
* -d, --debug PATH, Path to file to store info/debug messages (optional).

**Examples**

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



