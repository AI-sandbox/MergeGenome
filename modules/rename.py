################################################################################
# Change the chromosome notation for a particular chromosome data.
################################################################################

import os
import logging

from utils.io import read_vcf_file, write_vcf_file
from utils.vcf_utils import obtain_chromosomes, obtain_renamed_chrom, filter_by_chromosome, rename_chrom_field

def rename_chromosome(input_path: str, output_folder: str, rename_map: dict, logger: logging.Logger) -> None:
    """
    Changes the chromosome notation. If the chromosome notation is in the form '<chr_number>', it
    renames it to 'chr<chr_number>' (or vice-versa). If rename_map dictionary is provided,
    it renames the keys by the values.
    
    Args:
        input_path (str): path to .vcf file.
        output_folder (str): folder to save the output.
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
    
    # Ensure the file contains data of a particular chromosome
    chroms = obtain_chromosomes(data)
    assert len(chroms) == 1, 'The provided genomic file contains data for more than one chromosome.\
    Partition data into a separate VCF (partition command).'
    
    # Obtain chromosome
    chrom = chroms[0]
    
    # Obtain the new notation the chromosome
    new_chrom = obtain_renamed_chrom(True, chrom, rename_map)

    if new_chrom != chrom:
    
        # Renames variants/CHROM
        logger.debug(f'Renaming chromosome notation from {chrom} to {new_chrom}.')
        data = rename_chrom_field(data, chrom, new_chrom)
    
    # Define output name to .vcf file
    output_name = f'{os.path.basename(input_path)[:-4]}_{new_chrom}.vcf'
    
    # Write data for chromosome in output .vcf file
    logger.debug(f'Writing VCF data for chromosome {new_chrom}.')
    write_vcf_file(data, output_folder, output_name)
