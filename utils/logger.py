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
        log_path (str): path to .log or .txt file.
        
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

    if log_path is not None:
        
        # Remove content from log file
        if os.path.exists(log_path):
            os.remove(log_path)
        
        # Store log messages in log file
        logging.basicConfig(filename=log_path, encoding='utf-8', format='[%(asctime)s]-[%(levelname)s]: %(message)s')
    
    # Add both formatter and debug chanel to logger
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger