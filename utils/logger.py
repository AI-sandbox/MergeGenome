################################################################################
# Functions to assist with the logger
################################################################################

import logging
import os

def get_logger(name: str, log_path: str) -> logging.Logger:
    """
    Creates a logger to record the events that occur as the software runs.
    
    Args:
        name (str): class, file, function of the logger.
        log_path (str): path to log file.
        
    Returns:
        logger (logging.Logger): debug/information tracker.
    
    """

    # Instantiate logger and set log level to DEBUG
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # Specify debug outputs to stream
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    # Specify the layout of log records in the final output
    formatter = logging.Formatter('[%(asctime)s]-[%(levelname)s]: %(message)s')

    # Store log messages in log file
    if log_path is not None:
        
        if os.path.exists(log_path):
            # Remove content of log file
            os.remove(log_path)
        
        logging.basicConfig(filename=log_path, encoding='utf-8', format='[%(asctime)s]-[%(levelname)s]: %(message)s')
    
    # Add both formatter and debug chanel to logger
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    logging.basicConfig(filename='aa.log', encoding='utf-8', 
                        format='[%(asctime)s]-[%(name)s]-[%(levelname)s]: %(message)s', level=logging.DEBUG) #format="%(message)s"

    return logger


def parser_msg() -> str:
    """
    Specifies the custom usage of the parser.
    
    Returns
        (str)
    
    """
    return '''MergeGenome.py
        partition
            [-f, --file: path to input .vcf file)
            [-o, --output-folder: path to output folder]
            [-r, --rename-chr: rename (or not) chromosome notation]
            [-d, --debug: path to file storing info/debug messages]
        '''