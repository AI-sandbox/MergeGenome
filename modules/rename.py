################################################################################
# Changes the chromosome notation of a .vcf file with data from a single 
# or multiple chromosomes. By default, it renames the chromosome notation 
# from '<chr_number>' to 'chr<chr_number>' or from 'chr<chr_number>' to 
# '<chr_number>'.
################################################################################

import os
import logging

from utils.io import read_vcf_file, write_vcf_file
from utils.vcf_utils import obtain_chromosomes, obtain_renamed_chrom, rename_chrom_field


def rename_chromosome(query_path: str, output_folder: str, rename_map: dict, 
                      logger: logging.Logger) -> None:
    """
    Renames the chromosome notation of a query .vcf file with data from 
    a single or multiple chromosomes. By default, the chromosome notation 
    is renamed from '<chr_number>' to 'chr<chr_number>' or from 'chr<chr_number>' 
    to '<chr_number>'. If rename_map is provided, the change in the notation is 
    from keys (actual notation) to values (new notation).
    
    Args:
        query_path (str): path to query .vcf file with data for a single or 
        multiple chromosomes.
        output_folder (str): path to output folder to store separate .vcf files.
        rename_map (str): mapping from actual to new chromosome notation.
        logger (logging.Logger): debug/information tracker.
    
    Returns:
        (None)
    
    """
    
    # Read the query .vcf file using scikit-allel
    # with data from a single or multiple chromosomes
    logger.debug(f'Reading file {query_path}.')
    query = read_vcf_file(query_path)
    
    # Obtain the dimensions of the query
    logger.info(f'There are {len(query["variants/ID"])} SNPs and {len(query["samples"])} samples in total.')
    
    # Obtain the name of the chromosomes with available data
    chroms = obtain_chromosomes(query)
    logger.info(f'There are {len(chroms)} chromosomes in total.')
    
    # Define counter of modified chromosome notations
    counter = 0
    
    for chrom in chroms:
        # For each chromosome in the query...
        # Obtain new chromosome notation
        # If the chromosome notation is in neither '<chr_number>' or 'chr<chr_number>' formats,
        # the chromosome notation will remain the same, except if rename_map is provided
        new_chrom = obtain_renamed_chrom(True, chrom, rename_map)

        if new_chrom != chrom:
            # If the new chromosome notation is different from the actual
            # chromosome notation...
            # Rename variants/CHROM field from chrom to new_chrom
            logger.debug(f'Renaming chromosome notation from {chrom} to {new_chrom}.')
            query = rename_chrom_field(query, chrom, new_chrom)
            counter += 1
    
    if counter == 0:
        # If no modification was made...
        logger.debug('No change applied to variants/CHROM field.'\
                     'Consider specifying a mapping with --rename-map.')
    else:
        # Define output name to .vcf file with renamed data
        # The .vcf file has the same base name as the input file, but 
        # ending with '_renamed.vcf'
        output_name = f'{os.path.basename(query_path)[:-4]}_renamed.vcf'

        # Write query data with renamed variants/CHROM field in .vcf file
        # with name output_name inside folder output_folder
        logger.debug(f'Writing .vcf data with renamed variants/CHROM field in {output_folder}{output_name}.')
        # write_vcf_file(query, output_folder, output_name)
