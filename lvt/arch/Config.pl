#!/usr/bin/perl

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.4
#
# Copyright (c) 2022 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------
# 6 Jan 2012: Sujay Kumar, Initial Specification

#Find the architecture

if(defined($ENV{LVT_ARCH})){
   $sys_arch = $ENV{LVT_ARCH};
   # The Cray/Intel environment is almost identical to the Linux/Intel
   # environment.  There are two modifications that must be made to the
   # Linux/Intel configuration settings to make them work on the Cray.
   # So reset the sys_arch variable to "intel_ifc" and set a flag to 
   # enable the cray modifications.
   if($sys_arch eq "cray_ifc"){
      $sys_arch = "linux_ifc";
      $cray_modifications = 1;
   }
   else{
      $cray_modifications = 0;
   }
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the LVT architecture option \n";
   print "and linker using the LVT_ARCH variable.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}


if(defined($ENV{LVT_FC})){
   $sys_fc = $ENV{LVT_FC};
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the Fortran90 compiler \n";
   print "and linker using the LVT_FC variable.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}

$sys_cc = $ENV{LVT_CC};
if(defined($ENV{LVT_CC})){
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the C compiler \n";
   print "using the LVT_CC variable.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}

#print "Parallelism (0-serial, 1-dmpar, default=1): ";
#$par_lev=<stdin>;
#if($par_lev eq "\n"){
#   $par_lev=1;
#}
$par_lev=0;
if($par_lev == 1) {
   $sys_par = " -DSPMD ";
}

print "Optimization level (-3=strict checks with warnings, -2=strict checks, -1=debug, 0,1,2,3, default=2): ";
$opt_lev=<stdin>;
if($opt_lev eq "\n"){
   $opt_lev=2;
}

if($opt_lev == -3) {
   # Default flags for C.
    $sys_c_opt = "-g";
    if($sys_arch eq "linux_ifc"){
	$sys_opt = "-g -warn";
	$sys_opt .= 
           " -check bounds,format,output_conversion,pointers,stack,uninit";
	$sys_opt .= " -fp-stack-check -ftrapuv ";

	$sys_c_opt = "-g -Wall -Wcast-qual -Wcheck -Wdeprecated";
	$sys_c_opt .= " -Wextra-tokens -Wformat";
	$sys_c_opt .= " -Wformat-security -Wmissing-declarations";
	$sys_c_opt .= " -Wmissing-prototypes -Wpointer-arith -Wremarks";
	$sys_c_opt .= " -Wreturn-type -Wshadow -Wsign-compare";
	$sys_c_opt .= " -Wstrict-prototypes -Wtrigraphs -Wuninitialized";
	$sys_c_opt .= " -Wunused-function -Wunused-parameter";
	$sys_c_opt .= " -Wunused-variable -Wwrite-strings";
	# Run-time flags
	$sys_c_opt .= " -check=conversions,stack,uninit";
	$sys_c_opt .= " -fp-stack-check -fp-trap=common -fp-trap-all=common";
	$sys_c_opt .= " -ftrapuv";
    }
    elsif($sys_arch eq "linux_pgi") {
	print "Optimization level $opt_lev is not defined for $sys_arch.\n";
	print "Using '-g'\n";
	$sys_opt = "-g ";
    }
    elsif($sys_arch eq "linux_absoft") {
	print "Optimization level $opt_lev is not defined for $sys_arch.\n";
	print "Using '-g'\n";
	$sys_opt = "-g ";
    }
    elsif($sys_arch eq "linux_lf95") {
	print "Optimization level $opt_lev is not defined for $sys_arch.\n";
	print "Using '-g'\n";
    }
    elsif($sys_arch eq "Darwin_gfortran" || $sys_arch eq "linux_gfortran") {
	$sys_opt = "-g -Wall -Wcharacter-truncation";
	$sys_opt .= " -Wconversion-extra -Wextra -Wrealloc-lhs";
	$sys_opt .= " -Wrealloc-lhs-all";
	# Run-time options
	$sys_opt .= " -ffpe-trap=invalid,zero,overflow";
	$sys_opt .= " -fcheck=all,no-array-temps ";

	$sys_c_opt = "-g -Wall -Wextra -Wpedantic -Wformat -Wtraditional";
	$sys_c_opt .= " -Wconversion";
	# Run-time flags
	$sys_c_opt .= " -fstack-protector-all -fstack-check -ftrapv";
    }
    elsif($sys_arch eq "AIX") {
	print "Optimization level $opt_lev is not defined for $sys_arch.\n";
	print "Using '-g'\n";
	$sys_opt = "-g ";
    }
    elsif($sys_arch eq "cray_cray") {
	print "Optimization level $opt_lev is not defined for $sys_arch.\n";
	print "Using '-g'\n";
	$sys_opt = "-g ";
	$sys_c_opt = "-g ";
    }
}

if($opt_lev == -2) {
   if($sys_arch eq "linux_ifc") {
      $sys_opt = "-g -check bounds,format,output_conversion,pointers,stack,uninit -fp-stack-check -ftrapuv ";
      $sys_c_opt = "-g ";
   }
   elsif($sys_arch eq "linux_gfortran") {
      $sys_opt = "-g -Wall -fbounds-check ";
      $sys_c_opt = "-g ";
   }
   elsif($sys_arch eq "cray_cray") {
	   print "Optimization level $opt_lev is not defined for $sys_arch.\n";
      print "Using '-g'\n";
      $sys_opt = "-g ";
      $sys_c_opt = "-g ";
   }
}
if($opt_lev == -1) {
   $sys_opt = "-g -O0";
   $sys_c_opt = "-g -O0";
}
elsif($opt_lev == 0) {
   $sys_opt = "-O0 ";
   $sys_c_opt = "";
}
elsif($opt_lev == 1) {
   $sys_opt = "-O1 ";
   $sys_c_opt = "";
}
elsif($opt_lev == 2) {
   if($sys_arch eq "cray_cray") {
      $sys_opt = "-O2 -h ipa2,scalar0,vector0 ";
      $sys_c_opt = "";
   }
   else {
   $sys_opt = "-O2 ";
   $sys_c_opt = "";
   }
}
elsif($opt_lev == 3) {
   $sys_opt = "-O3 ";
   $sys_c_opt = "";
}
print "Assume little/big_endian data format (1-little, 2-big, default=2): ";
$use_endian=<stdin>;
if($use_endian eq "\n") {
   $use_endian=2
}


$sys_esmfmod_path = $ENV{LVT_MODESMF};
$sys_esmflib_path = $ENV{LVT_LIBESMF};

if((defined($ENV{LVT_MODESMF})) && (defined($ENV{LVT_LIBESMF}))){
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the ESMF library paths using\n";
   print "the LVT_MODESMF and LVT_LIBESMF variables.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}

print "Use GRIBAPI/ECCODES? (1-gribapi, 2-eccodes, default=2): ";
$use_gribapi=<stdin>;
chomp($use_gribapi);
if($use_gribapi eq ""){
   $use_gribapi=2;
}

if($use_gribapi == 1) {
   if(defined($ENV{LVT_GRIBAPI})){
      $sys_gribapi_path = $ENV{LVT_GRIBAPI};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_gribapi=$sys_gribapi_path.$inc;
      $lib_gribapi=$sys_gribapi_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the GRIBAPI library path using\n";
      print "the LVT_GRIBAPI variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   if(defined($ENV{LVT_JASPER}) && defined($ENV{LVT_OPENJPEG})){
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LVT_JASPER variable\n";
      print "or the LVT_OPENJPEG variable, not both.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   elsif(defined($ENV{LVT_JASPER})){
      $sys_jpeg2000_path = $ENV{LVT_JASPER};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-ljasper";
   }
   elsif(defined($ENV{LVT_OPENJPEG})){
      $sys_jpeg2000_path = $ENV{LVT_OPENJPEG};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-lopenjp2";
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LVT_JASPER variable\n";
      print "or the LVT_OPENJPEG variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}
elsif($use_gribapi == 2) {
   if(defined($ENV{LVT_ECCODES})){
      $sys_gribapi_path = $ENV{LVT_ECCODES};
      $inc = "/include/";
      if ($sys_arch eq "cray_cray") {
         $lib = "/lib64/";
      }
      else {
         $lib = "/lib/";
      }
      $inc_gribapi=$sys_gribapi_path.$inc;
      $lib_gribapi=$sys_gribapi_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the ECCODES library path using\n";
      print "the LVT_ECCODES variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   if(defined($ENV{LVT_JASPER}) && defined($ENV{LVT_OPENJPEG})){
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LVT_JASPER variable\n";
      print "or the LVT_OPENJPEG variable, not both.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   elsif(defined($ENV{LVT_JASPER})){
      $sys_jpeg2000_path = $ENV{LVT_JASPER};
      $inc = "/include/";
      if ($cray_modifications == 1) {
         $lib = "/lib/";
      }
      else {
         $lib = "/lib64/";
      }
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-ljasper";
   }
   elsif(defined($ENV{LVT_OPENJPEG})){
      $sys_jpeg2000_path = $ENV{LVT_OPENJPEG};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-lopenjp2";
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LVT_JASPER variable\n";
      print "or the LVT_OPENJPEG variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}
else {
   print "--------------ERROR---------------------\n";
   print "Note that several forcing datasets (e.g.\n";
   print "NLDAS, GDAS, ...) and several remotely sensed\n";
   print "datasets cannot be used without GRIBAPI/ECCODES\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}

print "Use NETCDF? (1-yes, 0-no, default=1): ";
$use_netcdf=<stdin>;
if($use_netcdf eq "\n"){
   $use_netcdf=1;
}
if($use_netcdf == 1) {
   print "NETCDF version (3 or 4, default=4): ";
   $netcdf_v=<stdin>;
   if($netcdf_v eq "\n"){
      $netcdf_v=4;
   }
   if(defined($ENV{LVT_NETCDF})){
      $sys_netcdf_path = $ENV{LVT_NETCDF};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_netcdf=$sys_netcdf_path.$inc;
      $lib_netcdf=$sys_netcdf_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the NETCDF path using\n";
      print "the LVT_NETCDF variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   print "NETCDF use shuffle filter? (1-yes, 0-no, default = 1): ";
   $netcdf_shuffle=<stdin>;
   if($netcdf_shuffle eq "\n"){
      $netcdf_shuffle=1;
   }
   print "NETCDF use deflate filter? (1-yes, 0-no, default = 1): ";
   $netcdf_deflate=<stdin>;
   if($netcdf_deflate eq "\n"){
      $netcdf_deflate=1;
   }
   print "NETCDF use deflate level? (1 to 9-yes, 0-no, default = 9): ";
   $netcdf_deflate_level=<stdin>;
   if($netcdf_deflate_level eq "\n"){
      $netcdf_deflate_level=9;
   }
}


print "Use HDF4? (1-yes, 0-no, default=1): ";
$use_hdf4=<stdin>;
if($use_hdf4 eq "\n"){
   $use_hdf4=1;
}
if($use_hdf4 == 1) {
   if(defined($ENV{LVT_HDF4})){
      $sys_hdf4_path = $ENV{LVT_HDF4};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_hdf4=$sys_hdf4_path.$inc;
      $lib_hdf4=$sys_hdf4_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the HDF4 path using\n";
      print "the LVT_HDF4 variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}

print "Use HDF5? (1-yes, 0-no, default=1): ";
$use_hdf5=<stdin>;
if($use_hdf5 eq "\n"){
   $use_hdf5=1;
}
if($use_hdf5 == 1) {
   if(defined($ENV{LVT_HDF5})){
      $sys_hdf5_path = $ENV{LVT_HDF5};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_hdf5=$sys_hdf5_path.$inc;
      $lib_hdf5=$sys_hdf5_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the HDF5 path using\n";
      print "the LVT_HDF5 variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}


print "Use HDFEOS? (1-yes, 0-no, default=1): ";
$use_hdfeos=<stdin>;
if($use_hdfeos eq "\n"){
   $use_hdfeos=1;
}
if($use_hdfeos == 1) {
   if($use_hdf4 == 0) {
      print "--------------ERROR---------------------\n";
      print "HDF4 should be enabled to have HDFEOS ";
      print "support working properly.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   else{
      if(defined($ENV{LVT_HDFEOS})){
         $sys_hdfeos_path = $ENV{LVT_HDFEOS};
         $inc = "/include/";
         $lib = "/lib/";
         $inc_hdfeos=$sys_hdfeos_path.$inc;
         $lib_hdfeos=$sys_hdfeos_path.$lib;
      }
      else {
         print "--------------ERROR---------------------\n";
         print "Please specify the HDFEOS path using\n";
         print "the LVT_HDFEOS variable.\n";
         print "Configuration exiting ....\n";
         print "--------------ERROR---------------------\n";
         exit 1;
      }
   }
}

print "Enable AFWA-specific grib configuration settings? (1-yes, 0-no, default=0): ";
$use_afwagrib=<stdin>;
if($use_afwagrib eq "\n"){
   $use_afwagrib=0;
}

print "Enable GeoTIFF support? (1-yes, 0-no, default=1): ";
$enable_geotiff=<stdin>;
if($enable_geotiff eq "\n"){
   $enable_geotiff=1;
}
if($enable_geotiff == 1) {
    if(defined($ENV{LVT_GDAL})){
	$sys_gdal_path = $ENV{LVT_GDAL};
	$lib = "/lib/";
	$lib_gdal=$sys_gdal_path.$lib;
    }
    else {
	print "--------------ERROR---------------------\n";
	print "Please specify the GDAL path using\n";
	print "the LVT_GDAL variable.\n";
	print "GDAL can be obtained from www.gdal.org\n";
	print "Configuration exiting ....\n";
	print "--------------ERROR---------------------\n";
	exit 1;
    }
    if(defined($ENV{LVT_FORTRANGIS})){
	$sys_fortrangis_path = $ENV{LVT_FORTRANGIS};
	$inc1 = "/libfortrangis/";
	$inc2 = "/libfortranc/";
	$lib = "/lib/";
#EMK
#	$inc_fortrangis1=$sys_fortrangis_path.$inc1;
#	$inc_fortrangis2=$sys_fortrangis_path.$inc2;
	$include = "/include/";
	$inc_fortrangis1=$sys_fortrangis_path.$include;
	$inc_fortrangis2=$sys_fortrangis_path.$include;

	$lib_fortrangis=$sys_fortrangis_path.$lib;
    }
    else {
	print "--------------ERROR---------------------\n";
	print "Please specify the FORTRANGIS path using\n";
	print "the LVT_FORTRANGIS variable.\n";
	print "FORTRANGIS can be obtained from http://fortrangis.sourceforge.net\n";
	print "Configuration exiting ....\n";
	print "--------------ERROR---------------------\n";
	exit 1;
    }
}

print "Use MATLAB support? (1-yes, 0-no, default=0): ";
$use_matlab=<stdin>;
if($use_matlab eq "\n"){
   $use_matlab=0;
}


if(defined($ENV{LVT_JPEG})){
   $libjpeg = "-L".$ENV{LVT_JPEG}."/lib"." -ljpeg";
}
else{
   $libjpeg = "-ljpeg";
}

if(defined($ENV{LVT_RPC})){
   $librpc = $ENV{LVT_RPC};
}
else{
   $librpc = "";
}

if($sys_arch eq "linux_ifc") {
   if ($use_endian == 1 ) {
      $cflags = "-c ".$sys_c_opt." -traceback -DIFC";
      $fflags77= "-c ".$sys_opt." -traceback -nomixed-str-len-arg -names lowercase -convert little_endian -assume byterecl ".$sys_par." -DIFC -I\$(MOD_ESMF) ";
      $fflags =" -c ".$sys_opt." -u -traceback -fpe0  -nomixed-str-len-arg -names lowercase -convert little_endian -assume byterecl ".$sys_par."-DIFC -I\$(MOD_ESMF) ";
   }
   else {
      $cflags = "-c ".$sys_c_opt." -traceback -DIFC";
      $fflags77= "-c ".$sys_opt." -traceback -nomixed-str-len-arg -names lowercase -convert big_endian -assume byterecl ".$sys_par." -DIFC -I\$(MOD_ESMF) ";
      $fflags =" -c ".$sys_opt." -u -traceback -fpe0  -nomixed-str-len-arg -names lowercase -convert big_endian -assume byterecl ".$sys_par."-DIFC -I\$(MOD_ESMF) ";
   }
   $ldflags= " -L\$(LIB_ESMF) -lesmf -lstdc++ -limf -lm -lrt ";
}
elsif($sys_arch eq "linux_pgi") {
   $cflags = "-c -DLITTLE_ENDIAN -DPGI";
   $fflags77= "-c ".$sys_opt."-C -s -Rb -Rs -g -gopt -Mbounds -Minform=inform -Minfo=all -DPGI -Mbyteswapio -r4 -i4 -Mpreprocess ".$sys_par."-I\$(MOD_ESMF) ";
   $fflags ="-c ".$sys_opt."-C -s -Rb -Rs -g -gopt -Mbounds -Minform=inform -Minfo=all  -DPGI -Mbyteswapio -r4 -i4 -Mpreprocess ".$sys_par." -I\$(MOD_ESMF) ";
   $ldflags= " -L\$(LIB_ESMF) -lesmf -pgcpplibs -ldl -lrt -Mlfs ";
}
elsif($sys_arch eq "linux_absoft") {
}
elsif($sys_arch eq "linux_lf95") {
}
elsif($sys_arch eq "Darwin_gfortran" || $sys_arch eq "linux_gfortran") {
   if($use_endian == 1) {
      $endian = "";
   }
   else {
      $endian = "-fconvert=big-endian";
   }
   $cflags = "-c ".$sys_c_opt." -DGFORTRAN ";
   $fflags77= "-c -pass-exit-codes ".$sys_opt." -fbacktrace ".$sys_par." ".$endian." -DHIDE_SHR_MSG -DNO_SHR_VMATH -DGFORTRAN -DHIDE_MPI -I\$(MOD_ESMF) ";
   $fflags =" -c -pass-exit-codes -ffree-line-length-0 ".$sys_opt." -fbacktrace ".$sys_par." ".$endian." -DHIDE_SHR_MSG -DNO_SHR_VMATH -DGFORTRAN -DHIDE_MPI -I\$(MOD_ESMF) ";
   $ldflags= " -L\$(LIB_ESMF) -lesmf -lstdc++ ";
}
elsif($sys_arch eq "AIX") {
   $cflags = "-c -w -g -qfullpath -q64 -qcpluscmt";
   $fflags77= "-c ".$sys_opt."-c -g -qkeepparm -qsuffix=f=f:cpp=F90 -q64 -WF,-DAIX, ".$sys_par." -I\$(MOD_ESMF) ";
   $fflags ="-c ".$sys_opt."-c -g -qkeepparm -qsuffix=f=f:cpp=F90 -q64 -WF,-DAIX, ".$sys_par." -I\$(MOD_ESMF) ";
   $ldflags= "-q64 -bmap:map -bloadmap:lm -lmass  -L\$(LIB_ESMF) -lesmf -lstdc++ -limf -lm -lrt";
}
elsif($sys_arch eq "cray_cray") {
   if($use_endian == 1) {
      $fflags77= "-c ".$sys_opt." ".$sys_par." -DCRAYFTN -I\$(MOD_ESMF) ";
      $fflags =" -c ".$sys_opt."-ef -Ktrap=fp  ".$sys_par."-DCRAYFTN -I\$(MOD_ESMF) ";
      $ldflags= " -hdynamic -L\$(LIB_ESMF) -lesmf -lstdc++ -lrt";
   }
   else {
      $fflags77= "-c ".$sys_opt." ".$sys_par." -DCRAYFTN -I\$(MOD_ESMF) ";
      $fflags =" -c ".$sys_opt."-ef -Ktrap=fp  ".$sys_par."-DCRAYFTN -I\$(MOD_ESMF) ";
      $ldflags= " -hbyteswapio -hdynamic -L\$(LIB_ESMF) -lesmf -lstdc++ -lrt";
   }

   $cflags = "-c ".$sys_c_opt." -DCRAYFTN";

}



if($par_lev == 1) {
   $ldflags = $ldflags;
}

if($use_gribapi == 1) {
   $fflags77 = $fflags77." -I\$(INC_GRIBAPI) ";
   $fflags = $fflags." -I\$(INC_GRIBAPI) ";
   $ldflags = $ldflags." -L\$(LIB_GRIBAPI) -lgrib_api_f90 -lgrib_api -L\$(LIB_JPEG2000) ".$ljpeg2000;
}
elsif($use_gribapi == 2) {
   $fflags77 = $fflags77." -I\$(INC_ECCODES)";
   $fflags = $fflags." -I\$(INC_ECCODES)";
   $ldflags = $ldflags." -L\$(LIB_ECCODES) -leccodes_f90 -leccodes -L\$(LIB_JPEG2000) ".$ljpeg2000;
}
if($use_netcdf == 1) {
   $fflags77 = $fflags77." -I\$(INC_NETCDF) ";
   $fflags = $fflags." -I\$(INC_NETCDF) ";
   if($netcdf_v == 3) {
      $ldflags = $ldflags." -L\$(LIB_NETCDF) -lnetcdf";
   }
   else{
      $ldflags = $ldflags." -L\$(LIB_NETCDF) -lnetcdff -lnetcdf";
   }
}
if($use_hdfeos == 1){
   $fflags77 = $fflags77." -I\$(INC_HDFEOS) ";
   $fflags = $fflags." -I\$(INC_HDFEOS) ";
   $ldflags = $ldflags." -L\$(LIB_HDFEOS) -lhdfeos -lGctp";
}
if($use_hdf4 == 1){
   $fflags77 = $fflags77." -I\$(INC_HDF4) ";
   $fflags = $fflags." -I\$(INC_HDF4) ";
   $ldflags = $ldflags." -L\$(LIB_HDF4) -lmfhdf -ldf ".$libjpeg." -lz ".$librpc;
}
if($use_hdf5 == 1){
   $fflags77 = $fflags77." -I\$(INC_HDF5) ";
   $fflags = $fflags." -I\$(INC_HDF5) ";
   $ldflags = $ldflags." -L\$(LIB_HDF5) -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl";
   if ( $cray_modifications == 1 ) {
      $ldflags = $ldflags." -lz";
   }
}

if($enable_geotiff== 1){
   $fflags77 = $fflags77." -I\$(INC_FORTRANGIS1) -I\$(INC_FORTRANGIS2)";
   $fflags = $fflags." -I\$(INC_FORTRANGIS1) -I\$(INC_FORTRANGIS2)";
   $ldflags = $ldflags." -L\$(LIB_FORTRANGIS) -lfortrangis -lfortranc -L\$(LIB_GDAL) -lgdal";
   if ( $cray_modifications == 1 ) {
      $ldflags = $ldflags." ".$ljpeg2000." ".$libjpeg." -lz -lstdc++";
   }
}

open(conf_file,">configure.lvt");
printf conf_file "%s%s\n","FC              = $sys_fc";
printf conf_file "%s%s\n","FC77            = $sys_fc";
printf conf_file "%s%s\n","LD              = $sys_fc";
printf conf_file "%s%s\n","CC              = $sys_cc";
printf conf_file "%s%s\n","AR              = ar";
printf conf_file "%s%s\n","MOD_ESMF        = $sys_esmfmod_path";
printf conf_file "%s%s\n","LIB_ESMF        = $sys_esmflib_path";
printf conf_file "%s%s\n","INC_JPEG2000      = $inc_jpeg2000";
printf conf_file "%s%s\n","LIB_JPEG2000      = $lib_jpeg2000";
if($use_gribapi == 2) {
printf conf_file "%s%s\n","INC_ECCODES     = $inc_gribapi";
printf conf_file "%s%s\n","LIB_ECCODES     = $lib_gribapi";
}
else {
printf conf_file "%s%s\n","INC_GRIBAPI     = $inc_gribapi";
printf conf_file "%s%s\n","LIB_GRIBAPI     = $lib_gribapi";
}
printf conf_file "%s%s\n","INC_NETCDF      = $inc_netcdf";
printf conf_file "%s%s\n","LIB_NETCDF      = $lib_netcdf";
printf conf_file "%s%s\n","INC_HDF4        = $inc_hdf4";
printf conf_file "%s%s\n","LIB_HDF4        = $lib_hdf4";
printf conf_file "%s%s\n","INC_HDF5        = $inc_hdf5";
printf conf_file "%s%s\n","LIB_HDF5        = $lib_hdf5";
printf conf_file "%s%s\n","INC_HDFEOS      = $inc_hdfeos";
printf conf_file "%s%s\n","LIB_HDFEOS      = $lib_hdfeos";
printf conf_file "%s%s\n","INC_FORTRANGIS1 = $inc_fortrangis1";
printf conf_file "%s%s\n","INC_FORTRANGIS2 = $inc_fortrangis2";
printf conf_file "%s%s\n","LIB_FORTRANGIS  = $lib_fortrangis";
printf conf_file "%s%s\n","LIB_GDAL        = $lib_gdal";
printf conf_file "%s%s\n","CFLAGS          = $cflags";
printf conf_file "%s%s\n","FFLAGS77        = $fflags77";
printf conf_file "%s%s\n","FFLAGS          = $fflags";
printf conf_file "%s%s\n","LDFLAGS         = $ldflags";
close(conf_file);

print "-----------------------------------------------------\n";
print " configure.lvt file generated successfully\n";
print "-----------------------------------------------------\n";

open(misc_file,">LVT_misc.h");

if($use_netcdf == 1) {
   if($netcdf_v == 3) {
      printf misc_file "%s\n","#define USE_NETCDF3 ";
   }
   else{
      printf misc_file "%s\n","#define USE_NETCDF4 ";
   }
}
else{
   printf misc_file "%s\n","#undef USE_NETCDF3 ";
   printf misc_file "%s\n","#undef USE_NETCDF4 ";
}

if($use_gribapi == 1) {
   printf misc_file "%s\n","#define USE_GRIBAPI ";
}
elsif($use_gribapi == 2) {
   printf misc_file "%s\n","#define USE_GRIBAPI ";
   printf misc_file "%s\n","#define USE_ECCODES ";
}
else{
   printf misc_file "%s\n","#undef USE_GRIBAPI ";
}

if($use_hdf4 == 1) {
   printf misc_file "%s\n","#define USE_HDF4 ";
}
else{
   printf misc_file "%s\n","#undef USE_HDF4 ";
}

if($use_hdfeos == 1) {
   printf misc_file "%s\n","#define USE_HDFEOS2 ";
}
else{
   printf misc_file "%s\n","#undef USE_HDFEOS2 ";
}

if($use_hdf5 == 1) {
   printf misc_file "%s\n","#define USE_HDF5 ";
}
else{
   printf misc_file "%s\n","#undef USE_HDF5 ";
}

if($use_afwagrib == 1) {
   printf misc_file "%s\n","#define AFWA_GRIB_CONFIGS ";
}
else{
   printf misc_file "%s\n","#undef AFWA_GRIB_CONFIGS ";
}

if($use_matlab == 1) {
   printf misc_file "%s\n","#define USE_MATLAB_SUPPORT ";
}
else{
   printf misc_file "%s\n","#undef USE_MATLAB_SUPPORT ";
}

if($enable_geotiff == 1) {
   printf misc_file "%s\n","#define USE_GDAL ";
}
else{
   printf misc_file "%s\n","#undef USE_GDAL ";
}

printf misc_file "%s\n","#undef INC_WATER_PTS";
printf misc_file "%s\n","#undef USE_GRIB2";
printf misc_file "%s\n","#undef SPMD";
close(misc_file);

open(netcdf_file,">LVT_NetCDF_inc.h");
printf netcdf_file "%s %d \n","#define NETCDF_shuffle ", $netcdf_shuffle;
printf netcdf_file "%s %d \n","#define NETCDF_deflate ", $netcdf_deflate;
printf netcdf_file "%s %d \n","#define NETCDF_deflate_level ", $netcdf_deflate_level;
close(netcdf_file);


