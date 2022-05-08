################################################################################
# Partitions the input file in a separate file per chromosome.
# Optionally, changes the chromosome notation. By default, it renames the 
# chromosome notation from '<chr_number>' to 'chr<chr_number>' or from
# 'chr<chr_number>' to '<chr_number>'.
################################################################################

import os
import logging

from utils.io import read_vcf_file, write_vcf_file
from utils.vcf_utils import obtain_chromosomes, obtain_renamed_chrom, filter_by_chromosome, rename_chrom_field


def partition_by_chromosome(query_path: str, output_folder: str, rename_chr: bool, 
                            rename_map: dict, logger: logging.Logger) -> None:
    """    
    Partitions the input .vcf query file in a separate .vcf file per chromosome.
    If rename_chr=True, the chromosome notation is renamed. By default, the 
    chromosome notation is renamed from '<chr_number>' to 'chr<chr_number>' or 
    from 'chr<chr_number>' to '<chr_number>'. If rename_map is provided, the 
    change in the notation is from keys (actual notation) to values (new 
    notation).
    
    Args:
        query_path (str): path to query .vcf file with data for multiple chromosomes.
        output_folder (str): path to output folder to store separate .vcf files.
        rename_chr (bool): to rename (or not) chromosome notation.
        rename_map (str): mapping from actual to new chromosome notation.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        (None)
    
    """

    # Read the query .vcf file using scikit-allel
    # with data from multiple chromosomes
    logger.debug(f'Reading file {query_path}.')
    query = read_vcf_file(query_path)
    
    # Obtain the dimensions of the query
    logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples in total.')
    
    # Obtain the name of the chromosomes with available data
    chroms = obtain_chromosomes(query)
    logger.info(f'There are {len(chroms)} chromosomes in total.')
    
    for chrom in chroms:
        # For each chromosome in the query...
        # Filter the query to include the data for the chromosome in particular
        logger.debug(f'Filtering query for chromosome {chrom}.')
        query_chrom = filter_by_chromosome(query.copy(), chrom)
        logger.info(f'There are {len(query_chrom["variants/ID"])} SNPs and '\
                    f'{len(query_chrom["samples"])} samples in chromosome {chrom}.')
        
        # Obtain new chromosome notation
        # If rename_chr = False, the chromosome notation will remain the same
        # If rename_chr = True and the chromosome notation is in neither
        # '<chr_number>' or 'chr<chr_number>' formats, the chromosome notation 
        # will remain the same, except if rename_map is provided
        new_chrom = obtain_renamed_chrom(rename_chr, chrom, rename_map)
        
        if new_chrom != chrom:
            # If the new chromosome notation is different from the actual
            # chromosome notation...
            # Rename variants/CHROM field from chrom to new_chrom
            logger.debug(f'Renaming chromosome notation from {chrom} to {new_chrom}.')
            query_chrom = rename_chrom_field(query_chrom, chrom, new_chrom)
        
        # Define output name to .vcf file for the chromosome in particular
        # The .vcf file has the same base name as the input file, but 
        # ending with the new chromosome notation '<new_chrom>.vcf'
        # Example: query_chr1.vcf
        output_name = f'{os.path.basename(query_path)[:-4]}_{new_chrom}.vcf'
        
        # Write query data for the chromosome in particular in .vcf file
        # with name output_name inside folder output_folder
        logger.debug(f'Writing .vcf data for chromosome {new_chrom} in {output_folder}{output_name}.')
        #write_vcf_file(query_chrom, output_folder, output_name)
