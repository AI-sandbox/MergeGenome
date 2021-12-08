## 2. Standardize chromosome notation

The standard chromosome nomenclature in CHROM field can be in two forms, **chr{chromosome_number}** or **{chromosome_number}**. Although both are correct, the notation adopted in each dataset has to be the same for merging purposes.

In order to change the chromosome nomenclature from {chromosome_number} to "chr{chromosome_number}" of all files per chromosome, you can use the script in `utils/partition_and_rename_chr.py`.

````{python}
def rename_chromosome(vcf_data, before, after):
    '''
    Objective: rename variants/CHROM in vcf_data.
    Input:
        - vcf_data: allel.read_vcf output.
        - before: chromosome name before renaming.
        - after: chromosome name after renaming.
    Output:
        - vcf_data: vcf_data with renamed variants/CHROM.
    '''
    
    ## Rename variants/CHROM
    vcf_data['variants/CHROM'] = np.where(vcf_data['variants/CHROM'] == before, after, vcf_data['variants/CHROM'])
    
    return vcf_data
````