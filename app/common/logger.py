""" Configurations for our logger """
import logging
import sys
from logging.handlers import RotatingFileHandler

FORMATTER = logging.Formatter(
    "%(asctime)s - %(levelname)s - %(name)s - %(messagePrefix)s - %(message)s"
)

def get_console_handler():
    """ Handler to display logs in the console """
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(FORMATTER)
    return console_handler

def get_log_file_name(module):
    """ Build the log file name """
    suffix = f'-{module}' if module else ''
    return f'./data/log/HLA_3D_{suffix}.log'

def get_file_handler(module=None):
    """ Handler to write logs into the log_file """
    file_handler = RotatingFileHandler(get_log_file_name(module), maxBytes=1*10**6, backupCount=2)
    file_handler.setFormatter(FORMATTER)
    return file_handler

def get_logger(logger_name, log_level=logging.INFO, module=None):
    """ Builds the logger for this application """
    logger = logging.getLogger(logger_name)
    logger.setLevel(log_level)
    logger.addHandler(get_console_handler())
    logger.addHandler(get_file_handler(module))
    logger.propagate = False
    return logger

# ######### logging ########

# def setup_logger(name, log_file, level=logging.DEBUG):
#     """To setup as many loggers as you want"""
#     formatter = logging.Formatter('%(name)s - %(asctime)s - %(levelname)s - %(message)s')
#     handler = logging.FileHandler(log_file)    
#     handler.setFormatter(formatter)

#     # Create & configure logger
#     logger = logging.getLogger(name)
#     logger.setLevel(level)
#     logger.addHandler(handler)

#     return logger

# logger = setup_logger('main_loger', 'data/log.log')