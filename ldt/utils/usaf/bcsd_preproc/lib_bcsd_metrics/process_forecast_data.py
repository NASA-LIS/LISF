#!/usr/bin/env python3
"""
#------------------------------------------------------------------------------
#
# SCRIPT:  process_forecast_data.py
#
# PURPOSE: Convert NNME GRIB2 files to netCDF.  Based on
# process_forecast_data.scr by Ryan Zamora.
#
# REQUIREMENTS as of 04 Nov 2021:
# * Python 3.9
# * Climate Data Operators (CDO) software
#
# REVISION HISTORY:
# 04 Nov 2021: Eric Kemp/SSAI, first version.
#
#------------------------------------------------------------------------------
"""

# Standard modules
import glob
import os
import subprocess
import sys

# Internal functions
def _usage():
    """Print command line usage."""
    txt = \
        f"[INFO] Usage: {sys.argv[0]} syr eyr fcst_init_monthday srcdir outdir"
    txt += " forcedir grid_description patchdir ic1 ic2 ic3"
    print(txt)

def _read_cmd_args():
    """Read command line arguments."""

    if len(sys.argv) != 12:
        print("[ERR] Invalid number of command line arguments!")
        _usage()
        sys.exit(1)

    args = {
        "syr" : sys.argv[1],
        "eyr" : sys.argv[2],
        "fcst_init_monthday" : sys.argv[3],
        "srcdir" : sys.argv[4],
        "outdir" : sys.argv[5],
        "forcedir" : sys.argv[6],
        "grid_description" : sys.argv[7],
        "patchdir" : sys.argv[8],
        "ic1" : sys.argv[9],
        "ic2" : sys.argv[10],
        "ic3" : sys.argv[11],
    }
    ic1 = args['ic1']
    ic2 = args['ic2']
    ic3 = args['ic3']
    args['all_ensmembers'] = ["00", "06", "12", "18", \
                              "00", "06", "12", "18", \
                              "00", "06", "12", "18"]
    args['all_monthdays'] = [ic1, ic1, ic1, ic1, \
                             ic2, ic2, ic2, ic2, \
                             ic3, ic3, ic3, ic3]

    return args

def _set_input_file_info(input_fcst_year, input_fcst_month, input_fcst_var):
    """Set input file information"""
    cutoff_refor_yyyymm = 201103
    cutoff_oper_yyyymm = 202112
    current_yyyymm = (100 * input_fcst_year) + input_fcst_month

    # Up to apr1 2011 - Refor_HPS (dlwsfc, dswsfc, q2m, wnd10m), Refor_FL
    # (prate, pressfc, tmp2m)
    if current_yyyymm <= cutoff_refor_yyyymm:
        if input_fcst_var in ["dlwsfc", "dswsfc", "q2m", "wnd10m"]:
            subdir = "Refor_HPS"
        elif input_fcst_var in ["prate", "pressfc", "tmp2m"]:
            subdir = "Refor_FL"
        file_pfx = input_fcst_var
        file_sfx = "time.grb2"
    #From may01 2011 to jan01 2021 - Oper_TS (all fields)
    elif current_yyyymm <= cutoff_oper_yyyymm:
        subdir = "Oper_TS"
        file_pfx = f"{input_fcst_var}.01"
        file_sfx = "daily.grb2"
    # From feb01 2021 - OperRT_TS (all fields)
    elif current_yyyymm > cutoff_oper_yyyymm:
        subdir = "OperRT_TS"
        file_pfx = f"{input_fcst_var}.01"
        file_sfx = "daily.grb2"
    return subdir, file_pfx, file_sfx

def _migrate_to_monthly_files(outdirs, temp_name, wanted_months,
                              fcst_init, args):
    """Migrate variables into monthly netCDF files"""

    outdir_6hourly = outdirs["outdir_6hourly"]
    outdir_monthly = outdirs["outdir_monthly"]
    final_name_pfx = f"{fcst_init['monthday']}.cfsv2."
    reftime = \
         f"{fcst_init['year']}-{fcst_init['month']}-{fcst_init['day']}"
    reftime += f",{fcst_init['hour']}:00:00,1hour"

    # Merge all variables into a single file
    cmd = "cdo --no_history merge "
    cmd += f"{outdir_6hourly}/junk1_*_{temp_name}"
    cmd += f" {outdir_6hourly}/junk2_{temp_name}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Subset data to only include the 9 forecast months
    cmd = "cdo --no_history selmon"
    for month in wanted_months:
        cmd += f",{month}"
    cmd += f" {outdir_6hourly}/junk2_{temp_name}" + \
        f" {outdir_6hourly}/junk3_{temp_name}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Convert all variable names and set missing value to -9999.f
    cmd = "cdo --no_history -setmissval,-9999. " + \
        "-chname,DSWRF_surface,SLRSF,DLWRF_surface,LWS," + \
        "PRATE_surface,PRECTOT,PRES_surface,PS," + \
        "SPFH_2maboveground,Q2M,TMP_2maboveground,T2M," + \
        "UGRD_10maboveground,U10M,VGRD_10maboveground,V10M" + \
        f" {outdir_6hourly}/junk3_{temp_name}" + \
        f" {outdir_6hourly}/junk4_{temp_name}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Add in windspeed variable
    cmd = "cdo --no_history "
    cmd += "aexpr,'WIND10M=sqrt(U10M*U10M + V10M*V10M)' "
    cmd += f"{outdir_6hourly}/junk4_{temp_name} "
    cmd += f"{outdir_6hourly}/junk5_{temp_name}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

    cmd = "cdo --no_history "
    cmd += "setattribute,WIND10M@long_name='Wind Speed',"
    cmd += "WIND10M@units='m/s',WIND10M@short_name='wnd10m',"
    cmd += "WIND10M@level='10 m above ground' "
    cmd += f"{outdir_6hourly}/junk5_{temp_name} "
    cmd += f"{outdir_6hourly}/junk6_{temp_name}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

    cmd = "cdo --no_history "
    cmd += f"-remapbil,{args['grid_description']} "
    cmd += f"-setreftime,{reftime} "
    cmd += f"{outdir_6hourly}/junk6_{temp_name} "
    cmd += f"{outdir_6hourly}/junk7_{temp_name}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Output 6-hourly data into monthly files
    cmd = "cdo --no_history -f nc4c -z zip_1 -splityearmon "
    cmd += f"{outdir_6hourly}/junk7_{temp_name} "
    cmd += f"{outdir_6hourly}/{final_name_pfx}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Output monthly data into monthly files
    cmd = \
        "cdo --no_history -L -f nc4c -z zip_1 -splityearmon -monmean "
    cmd += f"{outdir_6hourly}/junk7_{temp_name} "
    cmd += f"{outdir_monthly}/{final_name_pfx}"
    print(cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Cleanup intermediate files
    files = glob.glob(f"{outdir_6hourly}/junk*{temp_name}")
    for tmpfile in files:
        os.unlink(tmpfile)

def _print_reftime(fcst_init, ens_num):
    """Print reftime to standard out"""
    reftime = \
        f"{fcst_init['year']}-{fcst_init['month']}-{fcst_init['day']}"
    reftime += f",{fcst_init['hour']}:00:00,1hour"
    txt = f"[INFO] ENS{ens_num}: " + \
        f"{fcst_init['year']}-{fcst_init['monthday']}:" + \
        f"{fcst_init['hour']}"
    print(txt)

def _driver():
    """Main driver."""
    args = _read_cmd_args()
    fcst_init = {}
    fcst_init["monthday"] = args['fcst_init_monthday']
    outdirs = {}
    for year in range(int(args['syr']), (int(args['eyr']) + 1)):
        print(f"[INFO] {fcst_init['monthday']} {year}")

        if fcst_init['monthday'] == "jan01":
            fcst_init["year"] = year - 1
        else:
            fcst_init["year"] = year

        for ens_num in range(1, (len(args['all_ensmembers']) + 1)):
            monthday = args['all_monthdays'][ens_num - 1]
            temp_name = f"cfsv2.{fcst_init['year']}{monthday}.nc"
            fcst_init['date'] = f"{fcst_init['year']}{monthday}"
            fcst_init['month'] = int(monthday[0:2])
            fcst_init['day'] = int(monthday[2:4])
            fcst_init['hour'] = args['all_ensmembers'][ens_num - 1]
            fcst_init['timestring'] = f"{fcst_init['date']}{fcst_init['hour']}"
            wanted_months = []
            for i in range(int(fcst_init['month']), 13):
                wanted_months.append(i)
            for i in range(1, int(fcst_init['month'])):
                wanted_months.append(i)
            wanted_months = wanted_months[1:10]

            _print_reftime(fcst_init, ens_num)

            outdirs['outdir_6hourly'] = \
                f"{args['outdir']}/6-Hourly/" + \
                f"{fcst_init['monthday']}/{year}/ens{ens_num}"
            if not os.path.exists(outdirs['outdir_6hourly']):
                os.makedirs(outdirs['outdir_6hourly'])
            outdirs['outdir_monthly'] = \
                f"{args['outdir']}/Monthly/" + \
                f"{fcst_init['monthday']}/{year}/ens{ens_num}"
            if not os.path.exists(outdirs['outdir_monthly']):
                os.makedirs(outdirs['outdir_monthly'])

            for varname in ["prate", "pressfc", "tmp2m", "dlwsfc", "dswsfc",
                            "q2m", "wnd10m"]:
                print(f"[INFO] {varname}")
                subdir, file_pfx, file_sfx = \
                    _set_input_file_info(fcst_init['year'],
                                         fcst_init['month'],
                                         varname)
                indir = f"{args['forcedir']}/{subdir}/"
                indir += f"{fcst_init['year']}/{fcst_init['date']}"

                # Convert GRIB file to netCDF and handle missing/corrupted data
                cmd = f"{args['srcdir']}/convert_forecast_data_to_netcdf.py"
                cmd += f" {indir} {file_pfx} {fcst_init['timestring']}"
                cmd += f" {file_sfx} {outdirs['outdir_6hourly']}"
                cmd += f" {temp_name} {varname} {args['patchdir']}"
                print(cmd)
                subprocess.run(cmd, shell=True, check=True)

            _migrate_to_monthly_files(outdirs, temp_name, wanted_months,
                                      fcst_init, args)

    print("[INFO] Done processing CFSv2 forecast files")

if __name__ == "__main__":
    _driver()
