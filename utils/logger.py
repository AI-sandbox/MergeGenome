################################################################################
# Functions to assist with the logger
################################################################################

import logging

def get_logger(name: str) -> logging.Logger:
    """
    Creates a logger to record the events that occur as the software runs.
    Args:
        name (str): class, file, function of the logger.
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
    formatter = logging.Formatter('[%(asctime)s]-[%(name)s]-[%(levelname)s]: %(message)s')

    # Add both formatter and debug chanel to logger
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    return logger


def parser_msg() -> str:
    """
    Specifies the custom usage of the parser.
    Returns
        (str)
    """
    return '''MergeGenome.py
        evaluate
            [-f: path to the input .vcf file]
            [-o: path to output folder]
        '''