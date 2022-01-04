################################################################################
# Reads a single .vcf file with genomic data
# Filters the data by chromosome
# If the variants/CHROM nomenclature is in the form '1', ..., '38', 
# renames from '1' to 'chr1', ..., '38' to 'chr38'
# Writes the result in separate .vcf files, one for each chromosome
################################################################################

import sys
sys.path.append('/home/users/miriambt/my_work/dog-gen-to-phen/preprocessing')

from utils.vcf_utils import read_vcf_file, filter_by_chromosome, rename_chromosome, write_vcf_file

################################################################################

## USER INPUTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define path to single .vcf file with all chromosomes data
PATH = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/All_Pure_150k.vcf'

## Define path to output .vcf file with data for a single chromosome {}
OUTPUT_PATH = '/scratch/users/miriambt/data/dogs/formatted_data/4th_dataset/array/All_Pure_150k_chr{}.vcf'

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read .vcf file using skikit-allel
data = read_vcf_file(PATH)

for i in range (1, 39):
    ## For each chromosome number...
    
    ## Define path to output .vcf file with genomic data for chromosome i
    output_path = OUTPUT_PATH.format(i)
    
    ## Filter the data to include the SNPs for chromosome i
    data_chr_i = filter_by_chromosome(data.copy(), str(i))
    
    ## Rename variants/CHROM from {i} to chr{i}
    # If the variants/CHROM nomenclature is already in the chr{i} form, do nothing
    data_chr_i = rename_chromosome(data_chr_i, str(i), 'chr{}'.format(i))
    
    ## Write data for chromosome i in output .vcf file
    write_vcf_file(data_chr_i, output_path)
