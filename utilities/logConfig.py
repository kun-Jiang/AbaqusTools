import os
import logging
from logging.handlers import RotatingFileHandler




class logconfig():
    def __init__(self, FileName,working_dir):
        # Check if the logger has been configured to avoid duplicate log entries
        if not logging.getLogger().handlers:
            # Configure the log file
            log_file_path = working_dir + '\\' + FileName
            if os.path.exists(log_file_path):
                os.remove(log_file_path)
            logging.basicConfig(filename=log_file_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
            # Configure the console output
            console_handler = logging.StreamHandler()
            console_handler.setLevel(logging.INFO)
            console_handler.setFormatter(logging.Formatter('%(message)s'))
            logging.getLogger().addHandler(console_handler)