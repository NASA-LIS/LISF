#!/usr/bin/env python
"""
SCRIPT: autotune.py

Automates tuning of LIS Bratseth scheme based on prior runs, using
OBA (Observation, Background, Analaysis) files.

REVISION HISTORY:
27 Oct 2020:  Eric Kemp.  Initial specification.
"""

import configparser
import datetime
import os
import subprocess
import sys

class AutomateTuning:
    """Class to automate tuning of Bratseth covariances."""

    def _process_cmd_line(self):
        """Processes command line arguments."""
        if len(sys.argv) != 5:
            self.usage()
            sys.exit(1)

        self._process_cfg_file()

        varname = sys.argv[2]
        if varname not in ["rh2m", "t2m", "spd10m"]:
            print("[ERR] Invalid varname %s provided!" %(varname))
            sys.exit(1)
        self.varname = varname

        self.yyyymmddhh = sys.argv[3]
        year = int(self.yyyymmddhh[0:4])
        month = int(self.yyyymmddhh[4:6])
        day = int(self.yyyymmddhh[6:8])
        hour = int(self.yyyymmddhh[8:])
        self.enddt = \
            datetime.datetime(year=year, month=month, day=day, hour=hour)

        self.dd = sys.argv[4]
        days = int(self.dd)
        delta = datetime.timedelta(days=days)
        self.startdt = self.enddt - delta

    def _process_cfg_file(self):
        """Processes config file for this script."""
        if not os.path.exists(sys.argv[1]):
            print("[ERR] Cannot find config file %s" %(sys.argv[1]))
            print("Cannot continue....")
            sys.exit(1)

        cfgfile = sys.argv[1]
        config = configparser.ConfigParser()
        config.read(cfgfile)

        self.workdir = config.get('Input', 'workdir')
        self.scriptdir = config.get('Input', 'scriptdir')
        self.cfgdir = config.get('Input', 'cfgdir')
        self.bindir = config.get('Input', 'bindir')

    def __init__(self):
        """Initializes object"""
        self._process_cmd_line()
        self.use_blacklist = False

    def usage(self):
        """Print usage message for this script."""
        print("Usage: %s CFGFILE VARNAME YYYYMMDDHH DD" %(sys.argv[0]))
        print("   CFGFILE is name of config file")
        print("   VARNAME is name of variable to tune")
        print("   YYYYMMDDHH is end of training period")
        print("   DD is number of days in training period.")

    def create_blacklist(self):
        """High-level driver for creating blacklist for selected variable."""

        self.use_blacklist = True

        scriptfile = "%s/create_blacklist.py" %(self.scriptdir)
        if not os.path.exists(scriptfile):
            print("[WARN] Cannot find %s" %(scriptfile))
            print("Will skip blacklist for %s" %(self.varname))
            self.use_blacklist = False
            return

        cfgfile = "%s/blacklist_%s.cfg" %(self.cfgdir, self.varname)
        if not os.path.exists(cfgfile):
            print("[WARN] Cannot find %s" %(cfgfile))
            print("Will skip blacklist for %s" %(self.varname))
            self.use_blacklist = False
            return

        cmd = "%s %s blacklist_%s.txt" %(scriptfile, cfgfile, self.varname)
        cmd += " %s %s" %(self.yyyymmddhh, self.dd)
        print(cmd)
        try:
            subprocess.check_call(cmd, shell=True)
        except:
            print("[WARN] Problem running create_blacklist.py")
            print("Will skip blacklist for %s" %(varname))
            self.use_blacklist  = False
        return

    def customize_procoba_nwp(self):
        """Customizes config file for procOBA_NWP program."""
        cfgfile = "%s/procOBA_NWP.%s.config" %(self.cfgdir, self.varname)
        if not os.path.exists(cfgfile):
            print("[ERR] Cannot find %s" %(cfgfile))
            sys.exit(1)
            return

        # Copy and customize cfg file
        lines = open("%s" %(cfgfile), "r").readlines()
        cfgfile = "%s/procOBA_NWP.%s.config" %(self.workdir, self.varname)
        fd = open(cfgfile, "w")
        for line in lines:
            if "startyear:" in line:
                line = "startyear: %4.4d\n" %(self.startdt.year)
            elif "startmonth: " in line:
                line = "startmonth: %2.2d\n" %(self.startdt.month)
            elif "startday: " in line:
                line = "startday: %2.2d\n" %(self.startdt.day)
            elif "starthour: " in line:
                line = "starthour: %2.2d\n" %(self.startdt.hour)
            elif "endyear:" in line:
                line = "endyear: %4.4d\n" %(self.enddt.year)
            elif "endmonth: " in line:
                line = "endmonth: %2.2d\n" %(self.enddt.month)
            elif "endday: " in line:
                line = "endday: %2.2d\n" %(self.enddt.day)
            elif "endhour: " in line:
                line = "endhour: %2.2d\n" %(self.enddt.hour)
            elif "use_blacklist:" in line:
                option = "false\n"
                if self.use_blacklist:
                    option = "true\n"
                line = "use_blacklist: " + option
            elif "blacklist_file:" in line:
                line = "blacklist_file: blacklist_%s.txt" %(self.varname)
            fd.write(line)
        fd.close()

    # def run_procoba_nwp(self, varname):
    #     """Automate procOBA_NWP"""

    #     # Sanity check varname
    #     if varname not in ["precip", "rh2m", "spd10m", "t2m"]:
    #         print("[ERR] Unrecognized procOBA_NWP variable %" %(varname))
    #         print("Cannot continue")
    #         sys.exit(1)

    #     self.tune[varname] = True

    #     binfile = "%s/procOBA_NWP" %(self.bindir)
    #     if not os.path.exists(binfile):
    #         print("[WARN] Cannot find %s" %(binfile))
    #         print("Will skip tuning for %s" %(varname))
    #         self.tune[varname] = False
    #         return

    #     cfgfile = "%s/procOBA_NWP.%s.config" %(self.cfgdir, varname)
    #     if not os.path.exists(cfgfile):
    #         print("[WARN] Cannot find %s" %(cfgfile))
    #         print("Will skip tuning for %s" %(varname))
    #         self.tune[varname] = False
    #         return

    #     # Copy and customize cfg file
    #     lines = open("%s" %(cfgfile), "r").readlines()
    #     cfgfile = "%s/procOBA_NWP.%s.config" %(self.workdir, varname)
    #     fd = open(cfgfile, "w")
    #     for line in lines:
    #         if "startyear:" in line:
    #             line = "startyear: %4.4d\n" %(self.startdt.year)
    #         elif "startmonth: " in line:
    #             line = "startmonth: %2.2d\n" %(self.startdt.month)
    #         elif "startday: " in line:
    #             line = "startday: %2.2d\n" %(self.startdt.day)
    #         elif "starthour: " in line:
    #             line = "starthour: %2.2d\n" %(self.startdt.hour)
    #         elif "endyear:" in line:
    #             line = "endyear: %4.4d\n" %(self.enddt.year)
    #         elif "endmonth: " in line:
    #             line = "endmonth: %2.2d\n" %(self.enddt.month)
    #         elif "endday: " in line:
    #             line = "endday: %2.2d\n" %(self.enddt.day)
    #         elif "endhour: " in line:
    #             line = "endhour: %2.2d\n" %(self.enddt.hour)
    #         elif "use_blacklist:" in line:
    #             option = "false\n"
    #             if self.use_blacklist[varname]:
    #                 option = "true\n"
    #             line = "use_blacklist: " + option
    #         elif "blacklist_file:" in line:
    #             line = "blacklist_file: blacklist_%s.txt" %(varname)
    #         fd.write(line)
    #     fd.close()

    #     # Now run the program
    #     cmd = "mpirun -np 1 %s %s" %(binfile, cfgfile)
    #     print(cmd)
    #     try:
    #         subprocess.check_call(cmd, shell=True)
    #     except:
    #         print("[WARN] Problem running procOBA_NWP")
    #         print("Will skip tuning for %s" %(varname))
    #         self.tune[varname] = False
    #     return

if __name__ == "__main__":
    AUTOMATOR = AutomateTuning()

    # Process precip
    #AUTOMATOR.create_blacklist('precip')

    # Process 2-meter relative humidity
    AUTOMATOR.create_blacklist()
    AUTOMATOR.run_procoba_nwp()

    # Process 10-meter wind speed
    #AUTOMATOR.create_blacklist('spd10m')

    # Process 2-meter air temperature
    #AUTOMATOR.create_blacklist('t2m')

