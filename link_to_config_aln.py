# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 00:41:52 2022

@author: Anton
"""

import configparser
# Method to read config file settings
def read_config():
    config = configparser.ConfigParser()
    config.read('configurations_aln_file.ini')
    return config