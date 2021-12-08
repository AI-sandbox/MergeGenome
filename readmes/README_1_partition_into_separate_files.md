## 1. Partition data into separate .vcf files (one for each chromosome)

The computational time of preprocessing big .vcf files can be substantial. That is why, working with one .vcf file for each chromosome is in general better than treating the whole dataset at once. Moreover, splitting the dataset into separate files allows processing them in parallel.

In order to split the dataset in a .vcf file per chromosome, you can use the script in `utils/partition_and_rename_chr.py`. 

````{python}
def filter_by_chromosome(vcf_data, chrom):
    '''
    Objective:
        - Filter vcf_data to keep the SNPs of the specified chromosome.
    Input:
        - vcf_data: allel.read_vcf output.
        - chrom: name/number of the chromosome for which we want to select the SNPs.
    Output:
        - vcf_data: vcf_data containing only the selected SNPs of the desired chromosome.
    '''
    
    ## Find the indexes of the SNPs of the desired chromosomme
    indexes_chr = vcf_data['variants/CHROM'] == str(chrom)
    
    ## Select data for the desired SNPs
    for key in vcf_data.keys():
        if key != 'samples':
            vcf_data[key] = vcf_data[key][indexes_chr]
                
    return vcf_data
````
