#!/usr/bin/env python3
"""
SCRIPT: autotune.py

Automates tuning of LIS Bratseth scheme based on prior runs, using
OBA (Observation, Background, Analaysis) files.

REVISION HISTORY:
03 Nov 2020:  Eric Kemp.  Initial specification.
16 Dec 2020:  Eric Kemp.  Now calls create_blacklist via function call.
"""

import configparser
import datetime
import os
import subprocess
import sys

import create_blacklist

class AutomateTuning:
    """Class to automate tuning of Bratseth covariances."""

    def __init__(self):
        """Initializes object"""
        self._process_cmd_line()
        self.use_blacklist = False
        self.sigma2b = {}
        self.Lb = {}
        self.sigma2o = {}
        self.Lo = {}
        self.newlines = []

    def _process_cmd_line(self):
        """Processes command line arguments."""
        if len(sys.argv) not in [4, 5]:
            self.usage()
            sys.exit(1)

        self._process_cfg_file()

        self.yyyymmddhh = sys.argv[2]
        year = int(self.yyyymmddhh[0:4])
        month = int(self.yyyymmddhh[4:6])
        day = int(self.yyyymmddhh[6:8])
        hour = int(self.yyyymmddhh[8:])
        self.enddt = \
            datetime.datetime(year=year, month=month, day=day, hour=hour)

        self.dd = sys.argv[3]
        days = int(self.dd)
        delta = datetime.timedelta(days=days)
        self.startdt = self.enddt - delta


        if len(sys.argv) == 5:
            varname = sys.argv[4]
            if varname not in ["gage", "rh2m", "t2m", "spd10m", \
                               "cmorph", "geoprecip", "imerg", "ssmi"]:
                print("[ERR] Invalid varname %s provided!" %(varname))
                sys.exit(1)
            self.varname = varname
        else:
            self.varname = None

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
        self.lis_cfg_template = config.get('Input', 'lis_cfg_template')
        self.varlist =  config.get('Input', 'varlist').split(",")

    def usage(self):
        """Print usage message for this script."""
        print("Usage: %s CFGFILE YYYYMMDDHH DD [VARNAME]" %(sys.argv[0]))
        print("   CFGFILE is name of config file")
        print("   YYYYMMDDHH is end of training period")
        print("   DD is number of days in training period")
        print("   VARNAME is name of variable to tune (if applicable)")

    def create_blacklist(self):
        """High-level driver for creating blacklist for selected variable."""

        if self.varname is None:
            print("[ERR] Variable not specified on command line!")
            sys.exit(1)

        self.use_blacklist = True

        cfgfile = "%s/blacklist_%s.cfg" %(self.cfgdir, self.varname)
        if not os.path.exists(cfgfile):
            print("[WARN] Cannot find %s" %(cfgfile))
            print("Will skip blacklist for %s" %(self.varname))
            self.use_blacklist = False
            return

        blacklistfilename = "blacklist_%s.txt" %(self.varname)
        create_blacklist.create_blacklist(cfgfile, blacklistfilename, \
                                          self.yyyymmddhh, self.dd)
        return

    def check_gage_blacklist(self):
        """Checks to see if gage blacklist file exists."""
        filename = "blacklist_gage.txt"
        if os.path.exists(filename):
            self.use_blacklist = True

    def customize_procoba_nwp(self):
        """Customizes config file for procOBA_NWP program."""

        if self.varname is None:
            print("[ERR] Variable not specified on command line!")
            sys.exit(1)

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

    def customize_procoba_sat(self):
        """Customizes config file for procOBA_Sat."""

        if self.varname is None:
            print("[ERR] Variable not specified on command line!")
            sys.exit(1)

        cfgfile = "%s/procOBA_Sat.%s.config" %(self.cfgdir, self.varname)
        if not os.path.exists(cfgfile):
            print("[ERR] Cannot find %s" %(cfgfile))
            sys.exit(1)
            return

        # Copy and customize cfg file
        lines = open("%s" %(cfgfile), "r").readlines()
        cfgfile = "%s/procOBA_Sat.%s.config" %(self.workdir, self.varname)
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
                line = "blacklist_file: blacklist_gage.txt"
            fd.write(line)
        fd.close()

    def get_bratseth_err_settings(self, varname):
        """Fetches Bratseth error settings from files."""
        if varname in ["gage", "rh2m", "spd10m", "t2m"]:
            if varname == "gage":
                paramfile = "%s_nwp.param" %(varname)
            else:
                paramfile = "%s.param" %(varname)
            if not os.path.exists(paramfile):
                print("[WARN] Cannot find param file %s" %(paramfile))
                self.sigma2o[varname] = -9999
                self.sigma2b[varname] = -9999
                self.Lb[varname] = -9999
                return
            lines = open(paramfile, "r").readlines()
            for line in lines:
                # Here sfc reports are "obs", NWP is "background"
                if "SIGMA2_obs:" in line:
                    sigma2 = float(line.split()[-1])
                    self.sigma2o[varname] = sigma2
                    continue
                if "SIGMA2_back:" in line:
                    sigma2 = float(line.split()[-1])
                    self.sigma2b[varname] = sigma2
                    continue
                if "L_back:" in line:
                    Lb = float(line.split()[-1])*1000 # km to m
                    self.Lb[varname] = Lb
                    continue
        elif varname in ["cmorph", "geoprecip", "imerg", "ssmi"]:
            paramfile = "gage_%s_rescaled.param" %(varname)
            if not os.path.exists(paramfile):
                print("[WARN] Cannot find param file %s" %(paramfile))
                self.sigma2o[varname] = -9999
                self.Lo[varname] = -9999
                return
            lines = open(paramfile, "r").readlines()
            for line in lines:
                # NOTE:  Here the satellite obs data are the "background"
                # We now store them as obs since Bratseth will use NWP as
                # the background.
                if "SIGMA2_back:" in line:
                    sigma2 = float(line.split()[-1])
                    self.sigma2o[varname] = sigma2
                    continue
                if "L_back:" in line:
                    Lo = float(line.split()[-1])*1000 # km to m
                    self.Lo[varname] = Lo
                    continue

    def assemble_new_lines(self):
        """Assemble new lines for lis.config file."""
        lines = []

        if self.Lb["gage"] > 0:
            line = "AGRMET GALWEM Precip background error scale length (m):"
            line += " %s\n" %(self.Lb["gage"])
            lines.append(line)

        if self.sigma2b["gage"] > 0:
            line = "AGRMET GALWEM Precip background error variance:"
            line += " %s\n" %(self.sigma2b["gage"])
            lines.append(line)

        if self.sigma2o["gage"] > 0:
            line = "AGRMET GALWEM Precip Gauge observation error variance:"
            line += " %s\n" %(self.sigma2o["gage"])
            lines.append(line)

        if self.Lo["geoprecip"] > 0:
            line = \
        "AGRMET GALWEM Precip GEOPRECIP observation error scale length (m):"
            line += " %s\n" %(self.Lo["geoprecip"])
            lines.append(line)

        if self.sigma2o["geoprecip"] > 0:
            line = "AGRMET GALWEM Precip GEOPRECIP observation error variance:"
            line += " %s\n" %(self.sigma2o["geoprecip"])
            lines.append(line)

        if self.Lo["ssmi"] > 0:
            line = \
                "AGRMET GALWEM Precip SSMI observation error scale length (m):"
            line += " %s\n" %(self.Lo["ssmi"])
            lines.append(line)

        if self.sigma2o["ssmi"] > 0:
            line = "AGRMET GALWEM Precip SSMI observation error variance:"
            line += " %s\n" %(self.sigma2o["ssmi"])
            lines.append(line)

        if self.Lo["cmorph"] > 0:
            line = \
            "AGRMET GALWEM Precip CMORPH observation error scale length (m):"
            line += " %s\n" %(self.Lo["cmorph"])
            lines.append(line)

        if self.sigma2o["cmorph"] > 0:
            line = "AGRMET GALWEM Precip CMORPH observation error variance:"
            line += " %s\n" %(self.sigma2o["cmorph"])
            lines.append(line)

        if self.Lo["imerg"] > 0:
            line = \
                "AGRMET GALWEM Precip IMERG observation error scale length (m)"
            line += " %s\n" %(self.Lo["imerg"])
            lines.append(line)

        if self.sigma2o["imerg"] > 0:
            line = "AGRMET GALWEM Precip IMERG observation error variance:"
            line += " %s\n" %(self.sigma2o["imerg"])
            lines.append(line)

        if self.Lb["t2m"] > 0:
            line = "AGRMET GALWEM T2M background error scale length (m):"
            line += " %s\n" %(self.Lb["t2m"])
            lines.append(line)

        if self.sigma2b["t2m"] > 0:
            line = "AGRMET GALWEM T2M background error variance:"
            line += " %s\n" %(self.sigma2b["t2m"])
            lines.append(line)

        if self.sigma2o["t2m"] > 0:
            line = "AGRMET GALWEM T2M station observation error variance:"
            line += " %s\n" %(self.sigma2o["t2m"])
            lines.append(line)

        if self.Lb["rh2m"] > 0:
            line = "AGRMET GALWEM RH2M background error scale length (m):"
            line += " %s\n" %(self.Lb["rh2m"])
            lines.append(line)

        if self.sigma2b["rh2m"] > 0:
            line = "AGRMET GALWEM RH2M background error variance:"
            line += " %s\n" %(self.sigma2b["rh2m"])
            lines.append(line)

        if self.sigma2o["rh2m"] > 0:
            line = "AGRMET GALWEM RH2M station observation error variance:"
            line += " %s\n" %(self.sigma2o["rh2m"])
            lines.append(line)

        if self.Lb["spd10m"] > 0:
            line = "AGRMET GALWEM SPD10M background error scale length (m):"
            line += " %s\n" %(self.Lb["spd10m"])
            lines.append(line)

        if self.sigma2b["spd10m"] > 0:
            line = "AGRMET GALWEM SPD10M background error variance:"
            line += " %s\n" %(self.sigma2b["spd10m"])
            lines.append(line)

        if self.sigma2o["spd10m"] > 0:
            line = "AGRMET GALWEM SPD10M station observation error variance:"
            line += " %s\n" %(self.sigma2o["spd10m"])
            lines.append(line)

        self.newlines = lines

    def customize_lis_config(self):
        """Customize a lis.config file with latest error settings."""
        tmpl = self.lis_cfg_template
        if not os.path.exists(tmpl):
            print("[ERR] Cannot find LIS config template file %s" %(tmpl))
            sys.exit(1)
        lines = open(tmpl, "r").readlines()
        newfile = "%s/lis.config" %(self.workdir)
        fd = open(newfile, "w")
        for line in lines:
            # Pass through lines that don't specify error covariance settings
            if "AGRMET GALWEM" not in line:
                fd.write(line)
                continue
            if "error" not in line:
                fd.write(line)
                continue
            if "background" not in line and "observation" not in line:
                fd.write(line)
                continue
            if "variance" not in line and "scale length" not in line:
                fd.write(line)
                continue
            # At this point we have a Bratseth error covariance setting.
            # See if we replaced it.  If not, pass it through.
            string = line.split(":")[0]
            found = False
            for newline in self.newlines:
                if string in newline:
                    found = True
                    break
            if not found:
                fd.write(line)

        # Now append the new Bratseth settings at the end.
        fd.write("# NEW AUTOTUNED ERROR COVARIANCE SETTINGS\n")
        for line in self.newlines:
            fd.write(line)
        fd.close()

