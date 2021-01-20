#!/usr/bin/env python3
"""
SCRIPT: customize_procoba_sat.py

Customizes config file for procOBA_Sat program.

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

    # Check if blacklist file exists.  If so, we will use it.
    AUTOMATOR.check_gage_blacklist()

    # Customize the procOBA_Sat run.
    AUTOMATOR.customize_procoba_sat()
