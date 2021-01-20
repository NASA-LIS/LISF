#!/usr/bin/perl

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.3
#
# Copyright (c) 2020 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------
# 6 Jan 2012: Sujay Kumar, Initial Specification

#
# Process environment and configure options
#

if(defined($ENV{LIS_ARCH})){
   $sys_arch = $ENV{LIS_ARCH};
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
   print "Please specify the LIS architecture option \n";
   print "and linker using the LIS_ARCH variable.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}


if(defined($ENV{LIS_FC})){
   $sys_fc = $ENV{LIS_FC};
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the Fortran90 compiler \n";
   print "and linker using the LIS_FC variable.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}


if(defined($ENV{LIS_CC})){
   $sys_cc = $ENV{LIS_CC};
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the C compiler \n";
   print "using the LIS_CC variable.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}


print "Parallelism (0-serial, 1-dmpar, default=1): ";
$par_lev=<stdin>;
chomp($par_lev);
if($par_lev eq ""){
   $par_lev=1;
}

if($par_lev == 1) {
   $sys_par = "-DSPMD";
}
else {
   $sys_par = "";
}

my $par_table = ();
$par_table{"linux_ifc"}{"0"}       = "-DHIDE_MPI";
$par_table{"linux_ifc"}{"1"}       = "-DUSE_INCLUDE_MPI";
$par_table{"linux_pgi"}{"0"}       = "";
$par_table{"linux_pgi"}{"1"}       = "";
$par_table{"linux_absoft"}{"0"}    = "";
$par_table{"linux_absoft"}{"1"}    = "";
$par_table{"linux_lf95"}{"0"}      = "";
$par_table{"linux_lf95"}{"1"}      = "";
$par_table{"Darwin_gfortran"}{"0"} = "-DHIDE_MPI";
$par_table{"Darwin_gfortran"}{"1"} = "";
$par_table{"linux_gfortran"}{"0"}  = "-DHIDE_MPI";
$par_table{"linux_gfortran"}{"1"}  = "";
$par_table{"AIX"}{"0"}             = "";
$par_table{"AIX"}{"1"}             = "";

$sys_par_d = $par_table{$sys_arch}{$par_lev};

#print "Use openMP parallelism (1-yes, 0-no, default=0): ";
#$use_omp=<stdin>;
#chomp($use_omp);
#if($use_omp eq "") {
#   $use_omp=0
#}

#if ( $use_omp ) {
#   my %omp_table = ();
#
#   $omp_table{"linux_ifc"}       = "-openmp";
#   $omp_table{"linux_pgi"}       = "";
#   $omp_table{"linux_absoft"}    = "";
#   $omp_table{"linux_lf95"}      = "";
#   $omp_table{"Darwin_gfortran"} = "";
#   $omp_table{"linux_gfortran"}  = "";
#   $omp_table{"AIX"}             = "";
#
#   $sys_omp = $omp_table{$sys_arch};
#}
#else {
#   $sys_omp = "";
#}


print "Optimization level (-3=strict checks with warnings, -2=strict checks, -1=debug, 0,1,2,3, default=2): ";
$opt_lev=<stdin>;
chomp($opt_lev);
if($opt_lev eq ""){
   $opt_lev=2;
}

if($opt_lev == -3) {
   # Default flags for C.
   $sys_c_opt = "-g";
   if($sys_arch eq "linux_ifc"){
       $sys_opt = "-g -warn";
       $sys_opt .= 
	   " -check bounds,format,output_conversion,pointers,stack,uninit";
       $sys_opt .= " -fp-stack-check -ftrapuv";

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
      $sys_opt = "-g";
   }
   elsif($sys_arch eq "linux_absoft") {
      print "Optimization level $opt_lev is not defined for $sys_arch.\n";
      print "Using '-g'\n";
      $sys_opt = "-g";
   }
   elsif($sys_arch eq "linux_lf95") {
      print "Optimization level $opt_lev is not defined for $sys_arch.\n";
      print "Using '-g'\n";
   }
   elsif($sys_arch eq "Darwin_gfortran" || $sys_arch eq "linux_gfortran") {
      #print "Optimization level $opt_lev is not defined for $sys_arch.\n";
      #print "Using '-g'\n";
      #$sys_opt = "-g";
      $sys_opt = "-g -Wall -Wcharacter-truncation";
      $sys_opt .= " -Wconversion-extra -Wextra -Wpedantic -Wrealloc-lhs";
      $sys_opt .= " -Wrealloc-lhs-all";
      # Run-time options
      $sys_opt .= " -ffpe-trap=invalid,zero,overflow";
      $sys_opt .= " -fcheck=all,no-array-temps";

      $sys_c_opt = "-g -Wall -Wextra -Wpedantic -Wformat -Wtraditional";
      $sys_c_opt .= " -Wconversion";
      # Run-time flags
      $sys_c_opt .= " -fstack-protector-all -fstack-check -ftrapv";

   }
   elsif($sys_arch eq "AIX") {
      print "Optimization level $opt_lev is not defined for $sys_arch.\n";
      print "Using '-g'\n";
      $sys_opt = "-g";
   }
}
elsif($opt_lev == -2) {
   # Default flags for C.
   $sys_c_opt = "-g";
   if($sys_arch eq "linux_ifc"){
       $sys_opt = "-g ";
       $sys_opt .= 
	   " -check bounds,format,output_conversion,pointers,stack,uninit";
       $sys_opt .= " -fp-stack-check -ftrapuv";

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
      $sys_opt = "-g";
   }
   elsif($sys_arch eq "linux_absoft") {
      print "Optimization level $opt_lev is not defined for $sys_arch.\n";
      print "Using '-g'\n";
      $sys_opt = "-g";
   }
   elsif($sys_arch eq "linux_lf95") {
      print "Optimization level $opt_lev is not defined for $sys_arch.\n";
      print "Using '-g'\n";
      $sys_opt = "-g";
   }
   elsif($sys_arch eq "Darwin_gfortran" || $sys_arch eq "linux_gfortran") {
      print "Optimization level $opt_lev is not defined for $sys_arch.\n";
      print "Using '-g'\n";
      $sys_opt = "-g";
   }
   elsif($sys_arch eq "AIX") {
      print "Optimization level $opt_lev is not defined for $sys_arch.\n";
      print "Using '-g'\n";
      $sys_opt = "-g";
   }
}
elsif($opt_lev == -1) {
   $sys_opt = "-g -O0";
   $sys_c_opt = "-g -O0";
}
elsif($opt_lev == 0) {
   $sys_opt = "-O0";
   $sys_c_opt = "-O0";
}
elsif($opt_lev == 1) {
   $sys_opt = "-O1";
   $sys_c_opt = "";
}
elsif($opt_lev == 2) {
   $sys_opt = "-O2";
   $sys_c_opt = "";
}
elsif($opt_lev == 3) {
   $sys_opt = "-O3";
   $sys_c_opt = "";
}


print "Assume little/big_endian data format (1-little, 2-big, default=2): ";
$use_endian=<stdin>;
chomp($use_endian);
if($use_endian eq "") {
   $use_endian=2
}

my $end_table = ();

$end_table{"linux_ifc"}{"1"}       = "-convert little_endian";
$end_table{"linux_ifc"}{"2"}       = "-convert big_endian";
$end_table{"linux_pgi"}{"1"}       = "";
$end_table{"linux_pgi"}{"2"}       = "";
$end_table{"linux_absoft"}{"1"}    = "";
$end_table{"linux_absoft"}{"2"}    = "";
$end_table{"linux_lf95"}{"1"}      = "";
$end_table{"linux_lf95"}{"2"}      = "";
$end_table{"Darwin_gfortran"}{"1"} = "";
$end_table{"Darwin_gfortran"}{"2"} = "-fconvert=big-endian";
$end_table{"linux_gfortran"}{"1"}  = "";
$end_table{"linux_gfortran"}{"2"}  = "-fconvert=big-endian";
$end_table{"AIX"}{"1"}             = "";
$end_table{"AIX"}{"2"}             = "";

$sys_endian = $end_table{$sys_arch}{$use_endian};


if((defined($ENV{LIS_MODESMF})) && (defined($ENV{LIS_LIBESMF}))){
   $sys_esmfmod_path = $ENV{LIS_MODESMF};
   $sys_esmflib_path = $ENV{LIS_LIBESMF};
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the ESMF library paths using\n";
   print "the LIS_MODESMF and LIS_LIBESMF variables.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}


print "Use GRIBAPI/ECCODES? (0-neither, 1-gribapi, 2-eccodes, default=2): ";
$use_gribapi=<stdin>;
chomp($use_gribapi);
if($use_gribapi eq ""){
   $use_gribapi=2;
}

if($use_gribapi == 1) {
   if(defined($ENV{LIS_GRIBAPI})){
      $sys_gribapi_path = $ENV{LIS_GRIBAPI};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_gribapi=$sys_gribapi_path.$inc;
      $lib_gribapi=$sys_gribapi_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the GRIBAPI library path using\n";
      print "the LIS_GRIBAPI variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   if(defined($ENV{LIS_JASPER}) && defined($ENV{LIS_OPENJPEG})){
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LIS_JASPER variable\n";
      print "or the LIS_OPENJPEG variable, not both.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   elsif(defined($ENV{LIS_JASPER})){
      $sys_jpeg2000_path = $ENV{LIS_JASPER};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-ljasper";
   }
   elsif(defined($ENV{LIS_OPENJPEG})){
      $sys_jpeg2000_path = $ENV{LIS_OPENJPEG};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-lopenjp2";
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LIS_JASPER variable\n";
      print "or the LIS_OPENJPEG variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}
elsif($use_gribapi == 2) {
   if(defined($ENV{LIS_ECCODES})){
      $sys_gribapi_path = $ENV{LIS_ECCODES};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_gribapi=$sys_gribapi_path.$inc;
      $lib_gribapi=$sys_gribapi_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the ECCODES library path using\n";
      print "the LIS_ECCODES variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   if(defined($ENV{LIS_JASPER}) && defined($ENV{LIS_OPENJPEG})){
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LIS_JASPER variable\n";
      print "or the LIS_OPENJPEG variable, not both.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   if(defined($ENV{LIS_JASPER})){
      $sys_jpeg2000_path = $ENV{LIS_JASPER};
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
   elsif(defined($ENV{LIS_OPENJPEG})){
      $sys_jpeg2000_path = $ENV{LIS_OPENJPEG};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-lopenjp2";
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LIS_JASPER variable\n";
      print "or the LIS_OPENJPEG variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }

#    print "Grib Table Version (default = 128): ";
#    $grib_table_version=<stdin>;
#    if($grib_table_version eq "\n"){
#	$grib_table_version=128;
#    }
#    print "Grib Center Id (default = 57): ";
#    $grib_center_id=<stdin>;
#    if($grib_center_id eq "\n"){
#	$grib_center_id=57;
#    }
#    print "Grib Subcenter Id (default = 2): ";
#    $grib_subcenter_id=<stdin>;
#    if($grib_subcenter_id eq "\n"){
#	$grib_subcenter_id =2;
#    }
#    print "Grib Grid Id (default = 255): ";
#    $grib_grid_id=<stdin>;
#    if($grib_grid_id eq "\n"){
#	$grib_grid_id =255;
#    }
#    print "Grib Process Id (default = 88): ";
#    $grib_process_id=<stdin>;
#    if($grib_process_id eq "\n"){
#	$grib_process_id =88;
#    }
}
else {
   print "--------------WARNING---------------------\n";
   print "Note that several forcing datasets (e.g.\n";
   print "NLDAS, GDAS, ...) cannot be used without\n";
   print "GRIBAPI/ECCODES\n";
   print "--------------WARNING---------------------\n";
}


print "Enable AFWA-specific grib configuration settings? (1-yes, 0-no, default=0): ";
$use_afwagrib=<stdin>;
chomp($use_afwagrib);
if($use_afwagrib eq ""){
   $use_afwagrib=0;
}


print "Use NETCDF? (1-yes, 0-no, default=1): ";
$use_netcdf=<stdin>;
chomp($use_netcdf);
if($use_netcdf eq ""){
   $use_netcdf=1;
}

if($use_netcdf == 1) {
   print "NETCDF version (3 or 4, default=4): ";
   $netcdf_v=<stdin>;
   if($netcdf_v eq "\n"){
      $netcdf_v=4;
   }

   if(defined($ENV{LIS_NETCDF})){
      $sys_netcdf_path = $ENV{LIS_NETCDF};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_netcdf=$sys_netcdf_path.$inc;
      $lib_netcdf=$sys_netcdf_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the NETCDF path using\n";
      print "the LIS_NETCDF variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }

   print "NETCDF use shuffle filter? (1-yes, 0-no, default = 1): ";
   $netcdf_shuffle=<stdin>;
   chomp($netcdf_shuffle);
   if($netcdf_shuffle eq ""){
      $netcdf_shuffle=1;
   }

   print "NETCDF use deflate filter? (1-yes, 0-no, default = 1): ";
   $netcdf_deflate=<stdin>;
   chomp($netcdf_deflate);
   if($netcdf_deflate eq ""){
      $netcdf_deflate=1;
   }

   print "NETCDF use deflate level? (1 to 9-yes, 0-no, default = 9): ";
   $netcdf_deflate_level=<stdin>;
   chomp($netcdf_deflate_level);
   if($netcdf_deflate_level eq ""){
      $netcdf_deflate_level=9;
   }
}


print "Use HDF4? (1-yes, 0-no, default=1): ";
$use_hdf4=<stdin>;
chomp($use_hdf4);
if($use_hdf4 eq ""){
   $use_hdf4=1;
}

if($use_hdf4 == 1) {
   if(defined($ENV{LIS_HDF4})){
      $sys_hdf4_path = $ENV{LIS_HDF4};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_hdf4=$sys_hdf4_path.$inc;
      $lib_hdf4=$sys_hdf4_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the HDF4 path using\n";
      print "the LIS_HDF4 variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}


print "Use HDF5? (1-yes, 0-no, default=1): ";
$use_hdf5=<stdin>;
chomp($use_hdf5);
if($use_hdf5 eq ""){
   $use_hdf5=1;
}

if($use_hdf5 == 1) {
   if(defined($ENV{LIS_HDF5})){
      $sys_hdf5_path = $ENV{LIS_HDF5};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_hdf5=$sys_hdf5_path.$inc;
      $lib_hdf5=$sys_hdf5_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the HDF5 path using\n";
      print "the LIS_HDF5 variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}


print "Use HDFEOS? (1-yes, 0-no, default=1): ";
$use_hdfeos=<stdin>;
chomp($use_hdfeos);
if($use_hdfeos eq ""){
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
      if(defined($ENV{LIS_HDFEOS})){
         $sys_hdfeos_path = $ENV{LIS_HDFEOS};
         $inc = "/include/";
         $lib = "/lib/";
         $inc_hdfeos=$sys_hdfeos_path.$inc;
         $lib_hdfeos=$sys_hdfeos_path.$lib;
      }
      else {
         print "--------------ERROR---------------------\n";
         print "Please specify the HDFEOS path using\n";
         print "the LIS_HDFEOS variable.\n";
         print "Configuration exiting ....\n";
         print "--------------ERROR---------------------\n";
         exit 1;
      }
   }
}


print "Use MINPACK? (1-yes, 0-no, default=0): ";
$use_minpack=<stdin>;
chomp($use_minpack);
if($use_minpack eq ""){
   $use_minpack=0;
}

if($use_minpack == 1) {
   if(defined($ENV{LIS_MINPACK})){
      $sys_minpack_path = $ENV{LIS_MINPACK};
      $lib_minpack=$sys_minpack_path;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the MINPACK path using\n";
      print "the LIS_MINPACK variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}


print "Use LIS-CRTM? (1-yes, 0-no, default=0): ";
$use_crtm=<stdin>;
chomp($use_crtm);
if($use_crtm eq ""){
   $use_crtm=0;
}

if($use_crtm == 1) {
   if(defined($ENV{LIS_CRTM})){
      $sys_crtm_path = $ENV{LIS_CRTM};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_crtm=$sys_crtm_path.$inc;
      $lib_crtm=$sys_crtm_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the CRTM path using\n";
      print "the LIS_CRTM variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   if(defined($ENV{LIS_CRTM_PROF})){
      $sys_crtm_prof_path = $ENV{LIS_CRTM_PROF};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_crtm_prof=$sys_crtm_prof_path.$inc;
      $lib_crtm_prof=$sys_crtm_prof_path.$lib;
   }
   else {
      print "--------------ERROR----------------------------------\n";
      print "Please specify the CRTM profile utility path using\n";
      print "the LIS_CRTM_PROF variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR----------------------------------\n";
      exit 1;
   }
}


print "Use LIS-CMEM? (1-yes, 0-no, default=0): ";
$use_cmem=<stdin>;
chomp($use_cmem);
if($use_cmem eq ""){
   $use_cmem=0;
}

if($use_cmem == 1) {
   if(defined($ENV{LIS_CMEM})){
      $sys_cmem_path = $ENV{LIS_CMEM};
      $inc = "/src/";
      $lib = "/";
      $inc_cmem=$sys_cmem_path.$inc;
      $lib_cmem=$sys_cmem_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the CMEM path using\n";
      print "the LIS_CMEM variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}

print "Use LIS-LAPACK? (1-yes, 0-no, default=0): ";
$use_lapack=<stdin>;
chomp($use_lapack);
if($use_lapack eq ""){
   $use_lapck=0;
}

if($use_lapack == 1) {
   if(defined($ENV{LIS_LAPACK})){
      $sys_lapack_path = $ENV{LIS_LAPACK};
      $lib = "/";
      $lib_lapack=$sys_lapack_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the LAPACK path using\n";
      print "the LIS_LAPACK variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}

print "Use LIS-MKL-LAPACK? (1-yes, 0-no, default=0): ";
$use_mkllapack=<stdin>;
chomp($use_mkllapack);
if($use_mkllapack eq ""){
   $use_mkllapck=0;
}

if($use_mkllapack == 1) {
   if(defined($ENV{LIS_MKL_LAPACK})){
      $sys_lapack_path = $ENV{LIS_MKL_LAPACK};
      $lib = "/";
      $lib_lapack=$sys_lapack_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the MKL-LAPACK path using\n";
      print "the LIS_MKL_LAPACK variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}


if(defined($ENV{LIS_JPEG})){
   $libjpeg = "-L".$ENV{LIS_JPEG}."/lib"." -ljpeg";
}
else{
   $libjpeg = "-ljpeg";
}

if($ENV{ESMF_TRACE} > 0){
   $use_esmf_trace = 1 
}


#
# Compiler flags
#

if($sys_arch eq "linux_ifc") {
   $cflags = "-c ".$sys_omp." ".$sys_c_opt." -traceback -DIFC -DLINUX";
   $fflags77= "-c ".$sys_omp." ".$sys_opt." -traceback -nomixed_str_len_arg -names lowercase ".$sys_endian." -assume byterecl ".$sys_par." -DHIDE_SHR_MSG -DNO_SHR_VMATH -DIFC -DLINUX -I\$(MOD_ESMF) ".$sys_par_d;
   $fflags ="-c ".$sys_omp." ".$sys_opt." -u -traceback -fpe0 -nomixed_str_len_arg -names lowercase ".$sys_endian." -assume byterecl ".$sys_par." -DHIDE_SHR_MSG -DNO_SHR_VMATH -DIFC -DLINUX -I\$(MOD_ESMF) ".$sys_par_d;
   $ldflags= $sys_omp." -L\$(LIB_ESMF) -lesmf -lstdc++ -limf -lrt";
}
elsif($sys_arch eq "linux_pgi") {
   $cflags = "-c -DLITTLE_ENDIAN -DPGI";
   $fflags77= "-c ".$sys_opt." -C -s -Rb -Rs -g -gopt -Mbounds -Minform=inform -Minfo=all -DHIDE_SHR_MSG  -DHIDE_MPI -DUSE_INCLUDE_MPI -DNO_SHR_VMATH -DPGI -Mbyteswapio -r4 -i4 -Mpreprocess ".$sys_par." -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
   $fflags ="-c ".$sys_opt." -C -s -Rb -Rs -g -gopt -Mbounds -Minform=inform -Minfo=all -DHIDE_SHR_MSG  -DHIDE_MPI -DUSE_INCLUDE_MPI -DNO_SHR_VMATH -DPGI -Mbyteswapio -r4 -i4 -Mpreprocess" .$sys_par." -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
   $ldflags= " -L\$(LIB_ESMF) -lesmf -lesmf -lrt -lstd -lC -lnspgc -lpgc -lm -Mlfs";
}
elsif($sys_arch eq "linux_absoft") {
}
elsif($sys_arch eq "linux_lf95") {
}
elsif($sys_arch eq "linux_gfortran") {
   $cflags = "-c ".$sys_c_opt." -DGFORTRAN -DLINUX";
   $fflags77= "-c ".$sys_opt." -fbacktrace ".$sys_par." ".$sys_endian." -DHIDE_SHR_MSG -DNO_SHR_VMATH -DGFORTRAN -DLINUX ".$sys_par_d." -I\$(MOD_ESMF)";
   $fflags ="-c -ffree-line-length-0 ".$sys_opt." -fbacktrace ".$sys_par." ".$sys_endian." -DHIDE_SHR_MSG -DNO_SHR_VMATH -DGFORTRAN -DLINUX ".$sys_par_d." -I\$(MOD_ESMF)";
   $ldflags= " -L\$(LIB_ESMF) -lesmf -lstdc++ -lz";
}
elsif($sys_arch eq "Darwin_gfortran") {
   $cflags = "-c ".$sys_c_opt." -DGFORTRAN";
   $fflags77= "-c ".$sys_opt." -fbacktrace ".$sys_par." ".$sys_endian." -DHIDE_SHR_MSG -DNO_SHR_VMATH -DGFORTRAN ".$sys_par_d." -I\$(MOD_ESMF)";
   $fflags ="-c -ffree-line-length-0 ".$sys_opt." -fbacktrace ".$sys_par." ".$sys_endian." -DHIDE_SHR_MSG -DNO_SHR_VMATH -DGFORTRAN ".$sys_par_d." -I\$(MOD_ESMF)";
   $ldflags= " -L\$(LIB_ESMF) -lesmf -lstdc++";
}
elsif($sys_arch eq "AIX") {
   $cflags = "-c -w -g -qfullpath -q64 -qcpluscmt";
   $fflags77= "-c ".$sys_opt." -g -qkeepparm -qsuffix=f=f:cpp=F90 -q64 -WF,-DAIX,".$sys_par." -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
   $fflags ="-c ".$sys_opt." -g -qkeepparm -qsuffix=f=f:cpp=F90 -q64 -WF,-DAIX,".$sys_par." -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
   $ldflags= "-q64 -bmap:map -bloadmap:lm -lmass -L\$(LIB_ESMF) -lesmf -lstdc++ -limf -lm -lrt";
}

if($par_lev == 1) {
   if (index($sys_fc, "-mt_mpi") != -1) {
      $ldflags = $ldflags." -lmpi_mt";
   }
   elsif ($cray_modifications == 1) {
      $ldflags = $ldflags." -lmpich";
   }
   else{
      $ldflags = $ldflags." -lmpi";
   }
}

if($use_gribapi == 1) {
   $fflags77 = $fflags77." -I\$(INC_GRIBAPI)";
   $fflags = $fflags." -I\$(INC_GRIBAPI)";
   $ldflags = $ldflags." -L\$(LIB_GRIBAPI) -lgrib_api_f90 -lgrib_api -L\$(LIB_JPEG2000) ".$ljpeg2000;
}
elsif($use_gribapi == 2) {
   $fflags77 = $fflags77." -I\$(INC_ECCODES)";
   $fflags = $fflags." -I\$(INC_ECCODES)";
   $ldflags = $ldflags." -L\$(LIB_ECCODES) -leccodes_f90 -leccodes -L\$(LIB_JPEG2000) ".$ljpeg2000;
}
if($use_netcdf == 1) {
   $fflags77 = $fflags77." -I\$(INC_NETCDF)";
   $fflags = $fflags." -I\$(INC_NETCDF)";
   if($netcdf_v == 3) {
      $ldflags = $ldflags." -L\$(LIB_NETCDF) -lnetcdff -lnetcdf -lz"
   }
   else{
      $ldflags = $ldflags." -L\$(LIB_NETCDF) -lnetcdff -lnetcdf";
   }
}
if($use_hdfeos == 1){
   $fflags77 = $fflags77." -I\$(INC_HDFEOS)";
   $fflags = $fflags." -I\$(INC_HDFEOS)";
   $ldflags = $ldflags." -L\$(LIB_HDFEOS) -lhdfeos -lGctp";
}
if($use_hdf4 == 1){
   $fflags77 = $fflags77." -I\$(INC_HDF4)";
   $fflags = $fflags." -I\$(INC_HDF4)";
   $ldflags = $ldflags." -L\$(LIB_HDF4) -lmfhdf -ldf ".$libjpeg." -lz";
}
if($use_hdf5 == 1){
   $fflags77 = $fflags77." -I\$(INC_HDF5)";
   $fflags = $fflags." -I\$(INC_HDF5)";
   $ldflags = $ldflags." -L\$(LIB_HDF5) -lhdf5_fortran -lhdf5_hl -lhdf5";
   if ( $cray_modifications == 1 ) {
      $ldflags = $ldflags." -lz";
   }
}
if($use_crtm == 1){
   $fflags77 = $fflags77." -I\$(INC_CRTM) -I\$(INC_PROF_UTIL)";
   $fflags = $fflags." -I\$(INC_CRTM) -I\$(INC_PROF_UTIL)";
   $ldflags = $ldflags." -L\$(LIB_CRTM) -lCRTM -L\$(LIB_PROF_UTIL) -lProfile_Utility";
}

if($use_cmem == 1){
   $fflags77 = $fflags77." -I\$(INC_CMEM)";
   $fflags = $fflags." -I\$(INC_CMEM)";
   $ldflags = $ldflags." -L\$(LIB_CMEM) -llis_mem";
}

if($use_minpack == 1){
   $ldflags = $ldflags." -L\$(LIB_MINPACK) -lminpack";
}

if($use_lapack == 1){
   $ldflags = $ldflags." -L\$(LIB_LAPACK) -llapack -lblas";
}

if($use_mkllapack == 1){
   #Changed to be able to use mkl Wendy Sharples
   $ldflags = $ldflags." -L\$(LIB_LAPACK) -lmkl_rt";
}

if($use_esmf_trace == 1){
   $fflags = $fflags." -DESMF_TRACE";
}

#
# Kluge for LIS/Jules 5.0
#

$cflags = $cflags." -DLIS_JULES";
$fflags77 = $fflags77." -DLIS_JULES";
$fflags = $fflags." -DLIS_JULES";

#
# Write configure.lis and related files
#

open(conf_file,">configure.lis");
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
printf conf_file "%s%s\n","LIB_MINPACK     = $lib_minpack";
printf conf_file "%s%s\n","INC_CRTM        = $inc_crtm";
printf conf_file "%s%s\n","LIB_CRTM        = $lib_crtm";
printf conf_file "%s%s\n","INC_PROF_UTIL   = $inc_crtm_prof";
printf conf_file "%s%s\n","LIB_PROF_UTIL   = $lib_crtm_prof";
printf conf_file "%s%s\n","INC_CMEM        = $inc_cmem";
printf conf_file "%s%s\n","LIB_CMEM        = $lib_cmem";
printf conf_file "%s%s\n","LIB_LAPACK      = $lib_lapack";
printf conf_file "%s%s\n","CFLAGS          = $cflags";
printf conf_file "%s%s\n","FFLAGS77        = $fflags77";
printf conf_file "%s%s\n","FFLAGS          = $fflags";
printf conf_file "%s%s\n","LDFLAGS         = $ldflags";
close(conf_file);

print "-----------------------------------------------------\n";
print " configure.lis file generated successfully\n";
print "-----------------------------------------------------\n";


open(misc_file,">LIS_misc.h");
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

if($use_afwagrib == 1) {
   printf misc_file "%s\n","#define AFWA_GRIB_CONFIGS ";
}
else{
   printf misc_file "%s\n","#undef AFWA_GRIB_CONFIGS ";
}

if($use_hdf4 == 1) {
   printf misc_file "%s\n","#define USE_HDF4 ";
}
else{
   printf misc_file "%s\n","#undef USE_HDF4 ";
}

if($use_hdf5 == 1) {
   printf misc_file "%s\n","#define USE_HDF5 ";
}
else{
   printf misc_file "%s\n","#undef USE_HDF5 ";
}

if($use_minpack == 1) {
   printf misc_file "%s\n","#define USE_MINPACK ";
}
else{
   printf misc_file "%s\n","#undef USE_MINPACK ";
}

if($use_crtm == 1 || $use_cmem == 1) {
   printf misc_file "%s\n","#define RTMS ";
}
else{
   printf misc_file "%s\n","#undef RTMS ";
}

if($use_lapack == 1) {
   printf misc_file "%s\n","#define LAPACK ";
}
else{
   printf misc_file "%s\n","#undef LAPACK ";
}

if($use_mkllapack == 1) {
   printf misc_file "%s\n","#define MKL_LAPACK ";
}
else{
   printf misc_file "%s\n","#undef MKL_LAPACK ";
}


printf misc_file "%s\n","#undef INC_WATER_PTS";
printf misc_file "%s\n","#undef COUPLED";
close(misc_file);

open(netcdf_file,">LIS_NetCDF_inc.h");
printf netcdf_file "%s %d \n","#define NETCDF_shuffle ", $netcdf_shuffle;
printf netcdf_file "%s %d \n","#define NETCDF_deflate ", $netcdf_deflate;
printf netcdf_file "%s %d \n","#define NETCDF_deflate_level ", $netcdf_deflate_level;
close(netcdf_file);

#open(grib_file,">LIS_Grib_inc.h");
#printf grib_file "%s %d \n","#define GRIB_Table_Version ", $grib_table_version;
#printf grib_file "%s %d \n","#define GRIB_Center_Id ", $grib_center_id;
#printf grib_file "%s %d \n","#define GRIB_Subcenter_Id ", $grib_subcenter_id;
#printf grib_file "%s %d \n","#define GRIB_Grid_Id ", $grib_grid_id;
#printf grib_file "%s %d \n","#define GRIB_Process_Id ", $grib_process_id;
#close(grib_file);
