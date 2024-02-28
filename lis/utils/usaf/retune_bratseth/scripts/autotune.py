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
SCRIPT: autotune.py

Automates tuning of LIS Bratseth scheme based on prior runs, using
OBA (Observation, Background, Analaysis) files.

REVISION HISTORY:
03 Nov 2020:  Eric Kemp.  Initial specification.
16 Dec 2020:  Eric Kemp.  Now calls create_blacklist via function call.
13 Dec 2021:  Eric Kemp.  Added checks for missing satellite data.
19 Oct 2022:  Eric Kemp.  Added support for GFS settings.
"""

import configparser
import datetime
import os
import sys

import create_blacklist

class AutomateTuning:
    """Class to automate tuning of Bratseth covariances."""

    def __init__(self):
        """Initializes object"""
        self.workdir = "NULL"
        self.scriptdir = "NULL"
        self.cfgdir = "NULL"
        self.bindir = "NULL"
        self.lis_cfg_template = "NULL"
        self.varlist = "NULL"

        self.use_blacklist = False
        self.sigma2b = {}
        self.l_b = {}
        self.sigma2o = {}
        self.l_o = {}
        self.newlines = []
        self.for_gfs = False
        self._process_cmd_line()

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

        self.day = sys.argv[3]
        days = int(self.day)
        delta = datetime.timedelta(days=days)
        self.startdt = self.enddt - delta


        if len(sys.argv) == 5:
            varname = sys.argv[4]
            if varname not in ["gage", "rh2m", "t2m", "spd10m", \
                               "cmorph", "geoprecip", "imerg", "ssmi"]:
                print(f"[ERR] Invalid varname {varname} provided!")
                sys.exit(1)
            self.varname = varname
        else:
            self.varname = None

    def _process_cfg_file(self):
        """Processes config file for this script."""
        if not os.path.exists(sys.argv[1]):
            print(f"[ERR] Cannot find config file {sys.argv[1]}")
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
        self.varlist = config.get('Input', 'varlist').split(",")
        self.for_gfs = config.getboolean('Input', 'gfs')

    def usage(self):
        """Print usage message for this script."""
        print(f"Usage: {sys.argv[0]} CFGFILE YYYYMMDDHH DD [VARNAME]")
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

        cfgfile = f"{self.cfgdir}/blacklist_{self.varname}.cfg"
        if not os.path.exists(cfgfile):
            print(f"[WARN] Cannot find {cfgfile}")
            print(f"Will skip blacklist for {self.varname}")
            self.use_blacklist = False
            return

        blacklistfilename = f"blacklist_{self.varname}.txt"
        create_blacklist.create_blacklist(cfgfile, blacklistfilename, \
                                          self.yyyymmddhh, self.day)
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

        cfgfile = f"{self.cfgdir}/procOBA_NWP.{self.varname}.config"
        if not os.path.exists(cfgfile):
            print(f"[ERR] Cannot find {cfgfile}")
            sys.exit(1)
            return

        # Copy and customize cfg file
        with open(f"{cfgfile}", "r", encoding="ascii") as file:
            lines = file.readlines()
        cfgfile = f"{self.workdir}/procOBA_NWP.{self.varname}.config"
        entries = {
            "startyear: " : f"startyear: {self.startdt.year:04d}\n",
            "startmonth: " : f"startmonth: {self.startdt.month:02d}\n",
            "startday: " : f"startday: {self.startdt.day:02d}\n",
            "starthour: " : f"starthour: {self.startdt.hour:02d}\n",
            "endyear: " : f"endyear: {self.enddt.year:04d}\n",
            "endmonth: " : f"endmonth: {self.enddt.month:02d}\n",
            "endday: " : f"endday: {self.enddt.day:02d}\n",
            "endhour: " : f"endhour: {self.enddt.hour:02d}\n",
            "blacklist_file:" : f"blacklist_file: blacklist_{self.varname}.txt"
        }
        with open(cfgfile, "w", encoding="ascii") as file:
            for line in lines:
                if "use_blacklist:" in line:
                    option = "false\n"
                    if self.use_blacklist:
                        option = "true\n"
                    line = "use_blacklist: " + option
                else:
                    for item in entries.items():
                        if item[0] in line:
                            line = item[1]
                            break
                file.write(line)

    def customize_procoba_sat(self):
        """Customizes config file for procOBA_Sat."""

        if self.varname is None:
            print("[ERR] Variable not specified on command line!")
            sys.exit(1)

        cfgfile = f"{self.cfgdir}/procOBA_Sat.{self.varname}.config"
        if not os.path.exists(cfgfile):
            print(f"[ERR] Cannot find {cfgfile}")
            sys.exit(1)
            return

        # Copy and customize cfg file
        with open(f"{cfgfile}", "r", encoding="ascii") as file:
            lines = file.readlines()
        cfgfile = f"{self.workdir}/procOBA_Sat.{self.varname}.config"
        entries = {
            "startyear: " : f"startyear: {self.startdt.year:04d}\n",
            "startmonth: " : f"startmonth: {self.startdt.month:02d}\n",
            "startday: " : f"startday: {self.startdt.day:02d}\n",
            "starthour: " : f"starthour: {self.startdt.hour:02d}\n",
            "endyear: " : f"endyear: {self.enddt.year:04d}\n",
            "endmonth: " : f"endmonth: {self.enddt.month:04d}\n",
            "endday: " : f"endday: {self.enddt.day:04d}\n",
            "endhour: " : f"endhour: {self.enddt.hour:04d}\n",
            "blacklist_file:" : "blacklist_file: blacklist_gage.txt"
        }
        with open(cfgfile, "w", encoding="ascii") as file:
            for line in lines:
                if "use_blacklist:" in line:
                    option = "false\n"
                    if self.use_blacklist:
                        option = "true\n"
                    line = "use_blacklist: " + option
                else:
                    for item in entries.items():
                        if item[0] in line:
                            line = item[1]
                            break
                file.write(line)

    def get_bratseth_err_settings(self, varname):
        """Fetches Bratseth error settings from files."""
        if varname in ["gage", "rh2m", "spd10m", "t2m"]:
            paramfile = f"{varname}.param"
            if varname == "gage":
                paramfile = f"{varname}_nwp.param"
            if not os.path.exists(paramfile):
                print(f"[WARN] Cannot find param file {paramfile}")
                self.sigma2o[varname] = -9999
                self.sigma2b[varname] = -9999
                self.l_b[varname] = -9999
                return
            with open(paramfile, "r", encoding="ascii") as file:
                lines = file.readlines()
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
                    l_b = float(line.split()[-1])*1000 # km to m
                    self.l_b[varname] = l_b
                    continue
        elif varname in ["cmorph", "geoprecip", "imerg", "ssmi"]:
            paramfile = f"gage_{varname}_rescaled.param"
            if not os.path.exists(paramfile):
                print(f"[WARN] Cannot find param file {paramfile}")
                self.sigma2o[varname] = -9999
                self.l_o[varname] = -9999
                return
            with open(paramfile, "r", encoding="ascii") as file:
                lines = file.readlines()
            for line in lines:
                # NOTE:  Here the satellite obs data are the "background"
                # We now store them as obs since Bratseth will use NWP as
                # the background.
                if "SIGMA2_back:" in line:
                    sigma2 = float(line.split()[-1])
                    self.sigma2o[varname] = sigma2
                    continue
                if "L_back:" in line:
                    l_o = float(line.split()[-1])*1000 # km to m
                    self.l_o[varname] = l_o
                    continue

    def _assemble_new_lines_lb_gage(self, nwp, lines):
        if self.l_b["gage"] > 0:
            line = f"AGRMET {nwp} Precip background error scale length (m):"
            line += f" {self.l_b['gage']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_sigma2b_gage(self, nwp, lines):
        if self.sigma2b["gage"] > 0:
            line = f"AGRMET {nwp} Precip background error variance:"
            line += f" {self.sigma2b['gage']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_sigma2o_gage(self, nwp, lines):
        if self.sigma2o["gage"] > 0:
            line = f"AGRMET {nwp} Precip Gauge observation error variance:"
            line += f" {self.sigma2o['gage']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_lo_geoprecip(self, nwp, lines):
        if "geoprecip" in self.l_o:
            if self.l_o["geoprecip"] > 0:
                line = \
                   f"AGRMET {nwp} Precip GEOPRECIP observation error scale length (m):"
                line += f" {self.l_o['geoprecip']}\n"
                lines.append(line)
        return lines

    def _assemble_new_lines_sigma2o_geoprecip(self, nwp, lines):
        if "geoprecip" in self.sigma2o:
            if self.sigma2o["geoprecip"] > 0:
                line = f"AGRMET {nwp} Precip GEOPRECIP observation error variance:"
                line += f" {self.sigma2o['geoprecip']}\n"
                lines.append(line)
        return lines

    def _assemble_new_lines_lo_ssmi(self, nwp, lines):
        if "ssmi" in self.l_o:
            if self.l_o["ssmi"] > 0:
                line = \
                    f"AGRMET {nwp} Precip SSMI observation error scale length (m):"
                line += f" {self.l_o['ssmi']}%s\n"
                lines.append(line)
        return lines

    def _assemble_new_lines_sigma2o_ssmi(self, nwp, lines):
        if "ssmi" in self.sigma2o:
            if self.sigma2o["ssmi"] > 0:
                line = f"AGRMET {nwp} Precip SSMI observation error variance:"
                line += f" {self.sigma2o['ssmi']}\n"
                lines.append(line)
        return lines

    def _assemble_new_lines_lo_cmorph(self, nwp, lines):
        if "cmorph" in self.l_o:
            if self.l_o["cmorph"] > 0:
                line = \
                       f"AGRMET {nwp} Precip CMORPH observation error scale length (m):"
                line += f" {self.l_o['cmorph']}\n"
                lines.append(line)
        return lines

    def _assemble_new_lines_sigma2o_cmorph(self, nwp, lines):
        if "cmorph" in self.sigma2o:
            if self.sigma2o["cmorph"] > 0:
                line = f"AGRMET {nwp} Precip CMORPH observation error variance:"
                line += f" {self.sigma2o['cmorph']}\n"
                lines.append(line)
        return lines

    def _assemble_new_lines_lo_imerg(self, nwp, lines):
        if "imerg" in self.l_o:
            if self.l_o["imerg"] > 0:
                line = \
                    f"AGRMET {nwp} Precip IMERG observation error scale length (m):"
                line += f" {self.l_o['imerg']}\n"
                lines.append(line)
        return lines

    def _assemble_new_lines_sigma2o_imerg(self, nwp, lines):
        if "imerg" in self.sigma2o:
            if self.sigma2o["imerg"] > 0:
                line = f"AGRMET {nwp} Precip IMERG observation error variance:"
                line += f" {self.sigma2o['imerg']}\n"
                lines.append(line)
        return lines

    def _assemble_new_lines_lb_t2m(self, nwp, lines):
        if self.l_b["t2m"] > 0:
            line = f"AGRMET {nwp} T2M background error scale length (m):"
            line += f" {self.l_b['t2m']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_sigma2b_t2m(self, nwp, lines):
        if self.sigma2b["t2m"] > 0:
            line = f"AGRMET {nwp} T2M background error variance:"
            line += f" {self.sigma2b['t2m']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_sigma2o_t2m(self, nwp, lines):
        if self.sigma2o["t2m"] > 0:
            line = f"AGRMET {nwp} T2M station observation error variance:"
            line += f" {self.sigma2o['t2m']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_lb_rh2m(self, nwp, lines):
        if self.l_b["rh2m"] > 0:
            line = f"AGRMET {nwp} RH2M background error scale length (m):"
            line += f" {self.l_b['rh2m']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_sigma2b_rh2m(self, nwp, lines):
        if self.sigma2b["rh2m"] > 0:
            line = f"AGRMET {nwp} RH2M background error variance:"
            line += f" {self.sigma2b['rh2m']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_sigma2o_rh2m(self, nwp, lines):
        if self.sigma2o["rh2m"] > 0:
            line = f"AGRMET {nwp} RH2M station observation error variance:"
            line += f" {self.sigma2o['rh2m']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_lb_spd10m(self, nwp, lines):
        if self.l_b["spd10m"] > 0:
            line = f"AGRMET {nwp} SPD10M background error scale length (m):"
            line += f" {self.l_b['spd10m']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_sigma2b_spd10m(self, nwp, lines):
        if self.sigma2b["spd10m"] > 0:
            line = f"AGRMET {nwp} SPD10M background error variance:"
            line += f" {self.sigma2b['spd10m']}\n"
            lines.append(line)
        return lines

    def _assemble_new_lines_sigma2o_spd10m(self, nwp, lines):
        if self.sigma2o["spd10m"] > 0:
            line = f"AGRMET {nwp} SPD10M station observation error variance:"
            line += f" {self.sigma2o['spd10m']}\n"
            lines.append(line)
        return lines

    def assemble_new_lines(self):
        """Assemble new lines for lis.config file."""
        lines = []

        nwp = "GALWEM"
        if self.for_gfs:
            nwp = "GFS"

        lines = self._assemble_new_lines_lb_gage(nwp, lines)
        lines = self._assemble_new_lines_sigma2b_gage(nwp, lines)
        lines = self._assemble_new_lines_sigma2o_gage(nwp, lines)

        lines = self._assemble_new_lines_lo_geoprecip(nwp, lines)
        lines = self._assemble_new_lines_sigma2o_geoprecip(nwp, lines)

        lines = self._assemble_new_lines_lo_ssmi(nwp, lines)
        lines = self._assemble_new_lines_sigma2o_ssmi(nwp, lines)

        lines = self._assemble_new_lines_lo_cmorph(nwp, lines)
        lines = self._assemble_new_lines_sigma2o_cmorph(nwp, lines)

        lines = self._assemble_new_lines_lo_imerg(nwp, lines)
        lines = self._assemble_new_lines_sigma2o_imerg(nwp, lines)

        lines = self._assemble_new_lines_lb_t2m(nwp, lines)
        lines = self._assemble_new_lines_sigma2b_t2m(nwp, lines)
        lines = self._assemble_new_lines_sigma2o_t2m(nwp, lines)

        lines = self._assemble_new_lines_lb_rh2m(nwp, lines)
        lines = self._assemble_new_lines_sigma2b_rh2m(nwp, lines)
        lines = self._assemble_new_lines_sigma2o_rh2m(nwp, lines)

        lines = self._assemble_new_lines_lb_spd10m(nwp, lines)
        lines = self._assemble_new_lines_sigma2b_spd10m(nwp, lines)
        lines = self._assemble_new_lines_sigma2o_spd10m(nwp, lines)

        self.newlines = lines

    def customize_lis_config(self):
        """Customize a lis.config file with latest error settings."""
        tmpl = self.lis_cfg_template
        if not os.path.exists(tmpl):
            print(f"[ERR] Cannot find LIS config template file {tmpl}")
            sys.exit(1)
        with open(tmpl, "r", encoding="ascii") as file:
            lines = file.readlines()
        newfile = f"{self.workdir}/lis.config"
        with open(newfile, "w", encoding="ascii") as file:
            for line in lines:
                # Pass through lines that don't specify error covariance
                if "AGRMET GALWEM" not in line and "AGRMET GFS" not in line:
                    file.write(line)
                    continue
                if "error" not in line:
                    file.write(line)
                    continue
                if "background" not in line and "observation" not in line:
                    file.write(line)
                    continue
                if "variance" not in line and "scale length" not in line:
                    file.write(line)
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
                    file.write(line)

            # Now append the new Bratseth settings at the end.
            file.write("# NEW AUTOTUNED ERROR COVARIANCE SETTINGS\n")
            for line in self.newlines:
                file.write(line)
