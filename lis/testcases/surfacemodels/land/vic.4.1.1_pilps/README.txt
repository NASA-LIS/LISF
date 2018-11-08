README.txt - VIC data for PILPS-2E
Author: Ted Bohn

This directory contains forcing data and land surface parameters used in
PILPS-2E to model the Torne basin.  Directory contents include:

README.txt	this file
global.pilps-2e	global control file containing model settings
forcing/	directory containing standard-format atmospheric forcing data
forcing.orig/	directory containing original atmospheric forcing data
landcover/	directory containing vegetation and lake parameter files:
	world_lib_pilps_LAI.txt	vegetation class attributes
	veg.param.txt		vegetation parameters, lakes not considered
	veg.param.wlakes	vegatation parameters, lakes considered
	lake.uncal.init		lake parameters
soil/		directory containing soil parameter files:
	soillist.gen1		contains names of other soil files
	soillist.gen1.4.1.x	contains names of other soil files, for
				use with VIC 4.1.0 and later
	(other files)		other soil files - arcinfo format
snowband/	directory containing snow (elevation) band parameter files
result/		directory for storing model results


To run a simulation for this basin, do the following:

1. Replace all occurrences of <base_path> in global.pilps-2e and
   soil/soillist.gen1 with the full path of the directory containing
   global.pilps-2e.

2. If the files in the forcing directory have been compressed via gzip, you
   can still run vic if you set the COMPRESS option in global.pilps-2e to
   TRUE.  Otherwise, you must uncompress them (via gunzip) before running vic.

3. The ARC_SOIL option in global.pilps-2e must be set to TRUE (this
   is not the default).

4. If you are running VIC 4.1.0 or later, you will need to replace
	SOIL            <base_path>/soil/soillist.gen1
    with
	SOIL            <base_path>/soil/soillist.gen1.4.1.x
    in global.pilps-2e.

5. If you want to take snow (elevation) bands into account, un-comment
   (remove any '#' at the beginning of the line) the following line
   in global.pilps-2e:

	SNOW_BAND	5	<base_path>/snowband/torne.elevband.txt

   and comment (insert a '#' at the beginning of the line) the following
   line:

	SNOW_BAND	1

6. If you are running VIC 4.1.0 or later, you may opt to model lakes.  If so,
   un-comment the following lines in global.pilps-2e (remove any '#' at the
   beginning of the line):

	VEGPARAM	<base_path>/landcover/veg.param.wlakes
	LAKES		<base_path>/landcover/lake.uncal.init
	LAKE_PROFILE	TRUE

   In addition, you must make sure that the following line is commented
   out (insert a '#' at the beginning of the line):

	VEGPARAM	<base_path>/landcover/veg.param.txt

7. To run VIC using these input files, type
	vicNl -g <base_path>/global.pilps-2e
   where <base_path> is the full path of the directory containing
   global.pilps-2e.


NOTE ABOUT THE FORCING DATA:

The forcing files have been modified from their original form to work
with the official releases of vic 4.0.x.  The original forcing files can
be found in the forcing.orig directory.  The original forcing files
contain the following variables:

	Precip	Precipitation Rate (mm/s)
	Tair	Near Surface Air Temperature (C)
	Psurf	Near Surface Air Pressure (Pa)
	Wind_n	Near Surface Northward Wind Component (m/s)
	Wind_e	Near Surface Eastward Wind Component (m/s)
	Qair	Near Surface Specific Humidity (kg/kg)
	SWdown	Surface Incident Shortwave Radiation (W/m2)
	LWdown  Surface Incident Longwave Radiation (W/m2)

The new forcing files contain the following variables:

	PREC		Total Precipitation (mm)
	AIR_TEMP	Near Surface Air Temperature (C)
	PRESSURE	Near Surface Air Pressure (kPa)
	WIND		Near Surface Wind Speed (m/s)
	VP		Near Surface Vapor Pressure (kPa)
	SHORTWAVE	Surface Incident Shortwave Radiation (W/m2)
	LONGWAVE	Surface Incident Longwave Radiation (W/m2)

Variables were converted as follows:

	PREC	 = Precip * 3600 seconds/timestep
	AIR_TEMP = Tair
	PRESSURE = Psurf / 1000
	WIND	 = sqrt(Wind_n*Wind_n + Wind_e*Wind_e)
	VP	 = Qair * PRESSURE
	SHORTWAVE= SWdown
	LONGWAVE = LWdown

