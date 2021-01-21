#!/usr/bin/env python3
"""
SCRIPT: customize_procoba_nwp.py

Customizes config file for procOBA_NWP program.

REVISION HISTORY:
03 Nov 2020:  Eric Kemp.  Initial specification.
"""

import configparser
import datetime
import os
import subprocess
import sys

import autotune

if __name__ == "__main__":

    # Construct automator object, which reads the command line and config
    # file.
    AUTOMATOR = autotune.AutomateTuning()

    # Create new blacklist.
    AUTOMATOR.create_blacklist()

    # Customize the procOBA_NWP run.
    AUTOMATOR.customize_procoba_nwp()

