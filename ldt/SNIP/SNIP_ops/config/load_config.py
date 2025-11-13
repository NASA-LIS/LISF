#!/usr/bin/env python3

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------

"""
SCRIPT: settings.py

Script for create JSON config file for generating AMSR2 based snow depth
retrievals with AI/ML.

REVISION HISTORY:
15 Aug 2025: Kehan Yang. Initial Specification.
18 Aug 2025: Eric Kemp, Code cleanup.
"""

# Standard modules
from datetime import datetime
import json
from pathlib import Path
import logging
import sys

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('load_config')

# Pylint flags catching general exceptions, which is unreasonable.
# We disable that warning here.
# pylint: disable=W0718

class Config:
    """Class for processing JSON config file."""
    def __init__(self, config_path=None):
        self.proj = "EPSG:4326"
        self.flag_cold = False
        self.flag_rain = False
        self.flag_rfi = False
        self.target_resolution = 0.1 # degrees
        self.time_window_hours = 6
        self.land_frac_th = 90
        if config_path:
            try:
                # Try to load from file
                with open(config_path, 'r', encoding="ascii") as f:
                    config_data = json.load(f)

                # Set attributes from file
                for key, value in config_data.items():
                    if key == 'target_datetime':
                        setattr(self, key, datetime.strptime(value, "%Y%m%d%H%M"))
                    else:
                        setattr(self, key, value)

                # Fix data types
                if hasattr(self, 'project_path'):
                    self.project_path = Path(self.project_path)

                return  # Successfully loaded, exit early

            except (FileNotFoundError, json.JSONDecodeError, Exception) as e:
                logger.error("Error loading config file %s: %s! "
                             "Please provide a valid configuration file",
                             config_path, e)
                sys.exit(1)  # Exit with error code
