################################################################################
# Functions to check the validity of arguments and configurations
################################################################################

import os
from typing import List

def check_arguments(paths: List) -> None:
    """
    Checks that the file paths and formats are correct.
    Args:
        paths (List): paths to files.
    Returns:
        (None)
    """
    
    ## Check all the paths exist
    for path in paths:
        assert os.path.isfile(path), f'{path} not found. Is it a file?'
        assert path.endswith('.vcf'), f'{path} is not in a supported format. Convert it to .vcf file.'
    

def check_configurations() -> None:
    """
    
    Returns:
        (None)
    """

    # TODO