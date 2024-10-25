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
SCRIPT: customize_procoba_nwp.py

Customizes config file for procOBA_NWP program.

REVISION HISTORY:
03 Nov 2020:  Eric Kemp.  Initial specification.
"""

import autotune

if __name__ == "__main__":

    # Construct automator object, which reads the command line and config
    # file.
    AUTOMATOR = autotune.AutomateTuning()

    # Create new blacklist.
    AUTOMATOR.create_blacklist()

    # Customize the procOBA_NWP run.
    AUTOMATOR.customize_procoba_nwp()
