#!/usr/bin/env python3
"""
SCRIPT: customize_lis_config.py

Customizes lis.config file with latest Bratseth error settings.

REVISION HISTORY:
03 Nov 2020:  Eric Kemp.  Initial specification.
"""

import autotune

if __name__ == "__main__":

    # Construct automator object, which reads the command line and config
    # file.
    AUTOMATOR = autotune.AutomateTuning()

    # Get the covariance settings
    for var in AUTOMATOR.varlist:
        AUTOMATOR.get_bratseth_err_settings(var)

    # Assemble new lines for lis.config file
    AUTOMATOR.assemble_new_lines()

    # Customize the lis.config file.
    AUTOMATOR.customize_lis_config()


