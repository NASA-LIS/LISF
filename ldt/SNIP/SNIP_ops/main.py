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
SCRIPT: main.py

Driver script for generating snow depth retrievals from AMSR2 L1 files
using AI/ML.

REVISION HISTORY:
15 Aug 2025: Kehan Yang, Initial specification
18 Aug 2025: Eric Kemp, Code cleanup.
"""

# Standard modules
import argparse
import logging
import os
import sys

# Local modules
from config.load_config import Config
from data_processing.amsr2_reader import AMSR2DataProcessor
from ml_prediction.run_prediction import SnowDepthPredictor

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('main')

# Pylint flags catching general exceptions, which is excessive.  We
# disable the warnings.
# pylint: disable=W0718

class AMSR2SnowWorkflow:
    """Class for automating workflow for AMSR2 snow depth retrievals"""
    def __init__(self, config=None):
        self.config = config
        self.data_processor = AMSR2DataProcessor(config=config)
        self.sd_predictor = SnowDepthPredictor(config=config)

    def run_workflow(self):
        """Main workflow execution"""
        try:
            # Step 1: Read and process AMSR2 data
            target_datetime = self.config.target_datetime
            logging.info("Processing AMSR2 data for %s" ,target_datetime)

            datestr = target_datetime.strftime("%Y%m%d%H%M")
            # check if the file is already exist
            pmw_file = (f'{self.config.project_path}/{self.config.amsr2_merge_path}'
                        f'/AMSR2_L1R_combined_{datestr}.nc')
            if not os.path.exists(pmw_file):
                # if passive microwave input data is not merged, run pre-processing
                # to read AMSR2 L1R data and merge channels to one file.
                self.data_processor.process_l1r_data(target_datetime)

            # Step 2: ML SD prediction
            logging.info("Predicting snow depth with ML model")
            self.sd_predictor.run_pipeline(pmw_file)

        except Exception as e:
            logging.error("Workflow failed: %s", e)
            raise

    def _generate_output_path(self, dt):
        """Generate output filename"""
        filename = f"amsr2_snoice_0p1deg.{dt.strftime('%Y%m%d%H')}.nc"
        output_path = self.config.project_path / self.config.output_dir / filename
        return output_path

def process_single_config(config):
    """Process a single config file using AMSR2SnowWorkflow"""
    try:
        logging.info("Processing config for datetime %s", config.target_datetime)

        # Create workflow instance with config path
        amsr2workflow = AMSR2SnowWorkflow(config)
        amsr2workflow.run_workflow()

        logging.info("Successfully processed %s", config.target_datetime)
        return True, config

    except Exception as e:
        logging.error("Error processing %s : %s", config, e)
        return False, config



def main():
    """Main driver function"""
    parser = argparse.ArgumentParser(
        description='Generate SNIP config files for a period at ' + \
        '6-hour intervals (00, 06, 12, 18 UTC)')
    parser.add_argument('config_file',
                        help='Path to JSON configuration file')
    args = parser.parse_args()
    try:
        # Load configuration
        config = Config(args.config_file)

        # Process the single config
        process_single_config(config)

    except Exception as e:
        logger.error("Failed to process config %s", e)
        sys.exit(1)

if __name__ == "__main__":
    main()
