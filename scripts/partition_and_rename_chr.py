################################################################################
# Reads a single .vcf file
# Filters the data by chromosome
# Renames variants/CHROM nomenclature from '1' to 'chr1', ..., '38' to 'chr38'
# Writes the result in separate .vcf files, one for each chromosome
################################################################################

import sys
sys.path.append('/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing')

from utils.vcf_utils import read_vcf_file, filter_by_chromosome, rename_chromosome, write_vcf_file

################################################################################

## USER INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define base path to .vcf file with data for all chromosomes
PATH = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/All_Pure_150k.vcf'

## Define base path to output .vcf file with data for a single chromosome
OUTPUT_PATH = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/All_Pure_150k_chr*.vcf'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read .vcf file using skikit-allel
data = read_vcf_file(PATH)

for i in range (1, 39):
    ## Define path to output .vcf file with data for chromosome i
    output_path = OUTPUT_PATH.replace('*', str(i))
    
    ## Filter the data to include the SNPs for chromosome i
    data_chr_i = filter_by_chromosome(data.copy(), str(i))
    
    ## Rename variants/CHROM from '{number of chromosome}' to 'chr{number of chromosome}'
    data_chr_i = rename_chromosome(data_chr_i, str(i), 'chr'+str(i))
    
    ## Write data in output .vcf file
    write_vcf_file(data_chr_i, output_path)
