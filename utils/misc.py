################################################################################
# Functions to check the validity of arguments and configurations
################################################################################

import os
import argparse
from typing import List

def check_arguments(paths: List, rename_chr: bool = False, rename_map: dict = None) -> None:
    """
    Checks that the file paths and formats are correct.
    Args:
        paths (List): paths to files.
    Returns:
        (None)
    """
    
    # Check all the paths exist
    for path in paths:
        assert os.path.isfile(path), f'{path} not found. Is it a file?'
        assert path.endswith('.vcf'), f'{path} is not in a supported format. Convert it to .vcf file.'
    
    # Check rename_map is only available if rename_chr is set to True
    if (rename_map is not None) & (rename_chr == False):
        raise argparse.ArgumentTypeError(f'--rename_map is only supported for --rename_chr.')
    

def check_configurations() -> None:
    """
    
    Returns:
        (None)
    """

    # TODO