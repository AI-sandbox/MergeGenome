################################################################################
# Function to partition the input file in a separate file per chromosome.
# Optionally, change the chromosome notation.
################################################################################

import os
import logging

from utils.io import read_vcf_file, write_vcf_file
from utils.vcf_utils import obtain_chromosomes, obtain_renamed_chrom, filter_by_chromosome, rename_chromosome


def partition_by_chromosome(input_path: str, output_folder: str, rename_chr: bool, rename_map: dict, logger: logging.Logger) -> None:
    """
    Partitions the input file in a separate .vcf file per chromosome.
    Optionally, change the chromosome notation. If it is in the form '<chr_number>',
    change it to 'chr<chr_number>', and viceversa.
    
    Args:
        input_path (str): path to .vcf file.
        output_folder (str): folder to save the output.
        rename_chr (bool): rename (or not) the chromosome nomenclature.
        rename_map (str): dictionary with mapping from actual to new chromosome notation.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        (None)
    
    """

    # Read the .vcf file using scikit-allel
    logger.debug('Reading the provided genomic file.')
    data = read_vcf_file(input_path)
    
    # Obtain the dimensions of the data
    logger.info(f'There are {len(data["variants/ID"])} SNPs and {len(data["samples"])} samples in total.')
    
    # Obtain the name of the chromosomes with available data
    chroms = obtain_chromosomes(data)
    logger.info(f'There are {len(chroms)} chromosomes in total.')
    
    for chrom in chroms:
        
        # Filter the data to include the data for the chromosome in particular
        logger.debug(f'Filtering data for chromosome {chrom}.')
        data_chrom = filter_by_chromosome(data.copy(), chrom)
        logger.info(f'There are {len(data_chrom["variants/ID"])} SNPs and {len(data_chrom["samples"])} samples in chromosome {chrom}.')
        
        # Obtain the new notation the chromosome
        new_chrom = obtain_renamed_chrom(rename_chr, chrom, rename_map)
            
        if new_chrom != chrom:
                
            # Renames variants/CHROM
            logger.debug(f'Renaming chromosome notation from {chrom} to {new_chrom}.')
            data_chrom = rename_chromosome(data_chrom, chrom, new_chrom)
        
        # Define output name to .vcf file
        output_name = f'{os.path.basename(input_path)[:-4]}_{new_chrom}.vcf'
        
        # Write data for chromosome in output .vcf file
        logger.debug(f'Writing VCF data for chromosome {new_chrom}.')
        # write_vcf_file(data_chrom, output_folder, output_name)
