#!/usr/bin/perl

#-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
# NASA Goddard Space Flight Center
# Land Information System Framework (LISF)
# Version 7.5
#
# Copyright (c) 2024 United States Government as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
#-------------------------END NOTICE -- DO NOT EDIT-----------------------
# 6 Jan 2012: Sujay Kumar, Initial Specification

#Find the architecture

if(defined($ENV{LDT_ARCH})){
   $sys_arch = $ENV{LDT_ARCH};
   # The Cray/Intel environment is almost identical to the Linux/Intel
   # environment.  There are two modifications that must be made to the
   # Linux/Intel configuration settings to make them work on the Cray.
   # So reset the sys_arch variable to "linux_ifc" and set a flag to
   # enable the Cray modifications.
   if($sys_arch eq "cray_ifc"){
      $sys_arch = "linux_ifc";
      $cray_modifications = 1;
   }
   else{
      $cray_modifications = 0;
   }
   # The IBM/GNU environment is almost identical to the Linux/GNU
   # environment.  There is one modification that must be made to the
   # Linux/GNU configuration settings to make them work on the IBM Power9.
   # So reset the sys_arch variable to "linux_gfortran" and set a flag to
   # enable the IBM modifications.
   if($sys_arch eq "ibm_gfortran"){
      $sys_arch = "linux_gfortran";
      $ibm_modifications = 1;
   }
   else{
      $ibm_modifications = 0;
   }
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the LDT architecture option \n";
   print "and linker using the LDT_ARCH variable.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}


if(defined($ENV{LDT_FC})){
   $sys_fc = $ENV{LDT_FC};
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the Fortran90 compiler \n";
   print "and linker using the LDT_FC variable.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}


if(defined($ENV{LDT_CC})){
   $sys_cc = $ENV{LDT_CC};
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the C compiler \n";
   print "using the LDT_CC variable.\n";
   print "Configuration exiting ....\n";
   print "--------------ERROR---------------------\n";
   exit 1;
}


print "Parallelism (0-serial, 1-dmpar, default=0): ";
$par_lev=<stdin>;
chomp($par_lev);
if($par_lev eq ""){
   $par_lev=0;
}

if($par_lev == 1) {
   $sys_par = " -DSPMD ";
}
else {
   $sys_par = "";
}

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
      $sys_opt = "-g -check bounds,format,output_conversion,pointers,stack,uninit ";
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
   $sys_opt = "-g ";
   $sys_c_opt = "-g ";
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
      # EMK 14 Apr 2022...For best conformity to IEEE standard, use fp0
      $sys_opt = "-O2 -h ipa2,scalar0,vector0,fp0 ";
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
chomp($use_endian);
if($use_endian eq "") {
   $use_endian=2
}


$sys_esmfmod_path = $ENV{LDT_MODESMF};
$sys_esmflib_path = $ENV{LDT_LIBESMF};

if((defined($ENV{LDT_MODESMF})) && (defined($ENV{LDT_LIBESMF}))){
}
else{
   print "--------------ERROR---------------------\n";
   print "Please specify the ESMF paths using\n";
   print "the LDT_MODESMF and LDT_LIBESMF variables.\n";
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
   if(defined($ENV{LDT_GRIBAPI})){
      $sys_gribapi_path = $ENV{LDT_GRIBAPI};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_gribapi=$sys_gribapi_path.$inc;
      $lib_gribapi=$sys_gribapi_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the GRIBAPI library path using\n";
      print "the LDT_GRIBAPI variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   if(defined($ENV{LDT_JASPER}) && defined($ENV{LDT_OPENJPEG})){
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LDT_JASPER variable\n";
      print "or the LDT_OPENJPEG variable, not both.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   elsif(defined($ENV{LDT_JASPER})){
      $sys_jpeg2000_path = $ENV{LDT_JASPER};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-ljasper";
   }
   elsif(defined($ENV{LDT_OPENJPEG})){
      $sys_jpeg2000_path = $ENV{LDT_OPENJPEG};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-lopenjp2";
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LDT_JASPER variable\n";
      print "or the LDT_OPENJPEG variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}
elsif($use_gribapi == 2) {
   if(defined($ENV{LDT_ECCODES})){
      $sys_gribapi_path = $ENV{LDT_ECCODES};
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
      print "the LDT_ECCODES variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   if(defined($ENV{LDT_JASPER}) && defined($ENV{LDT_OPENJPEG})){
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LDT_JASPER variable\n";
      print "or the LDT_OPENJPEG variable, not both.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
   elsif(defined($ENV{LDT_JASPER})){
      $sys_jpeg2000_path = $ENV{LDT_JASPER};
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
   elsif(defined($ENV{LDT_OPENJPEG})){
      $sys_jpeg2000_path = $ENV{LDT_OPENJPEG};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_jpeg2000=$sys_jpeg2000_path.$inc;
      $lib_jpeg2000=$sys_jpeg2000_path.$lib;
      $ljpeg2000="-lopenjp2";
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the JPEG2000 library path using\n";
      print "either the LDT_JASPER variable\n";
      print "or the LDT_OPENJPEG variable.\n";
      print "Configuration exiting ....\n";
      print "--------------ERROR---------------------\n";
      exit 1;
   }
}
else {
   print "--------------WARNING---------------------\n";
   print "Note that several forcing datasets (e.g.\n";
   print "NLDAS, GDAS, ...) and several remotely sensed\n";
   print "datasets cannot be used without GRIBAPI/ECCODES\n";
   print "--------------WARNING---------------------\n";
}


$use_omp=0;
#print "Use openMP parallelism (1-yes, 0-no, default=0): ";
#$use_omp=<stdin>;
#if($use_omp eq "\n") {
#   $use_omp=0
#}

$use_netcdf=1;
if($use_netcdf == 1) {
   print "NETCDF version (3 or 4, default=4)?: ";
   $netcdf_v=<stdin>;
   if($netcdf_v eq "\n"){
      $netcdf_v=4;
   }
   if(defined($ENV{LDT_NETCDF})){
      $sys_netcdf_path = $ENV{LDT_NETCDF};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_netcdf=$sys_netcdf_path.$inc;
      $lib_netcdf=$sys_netcdf_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the NETCDF path using\n";
      print "the LDT_NETCDF variable.\n";
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
   if(defined($ENV{LDT_HDF4})){
      $sys_hdf4_path = $ENV{LDT_HDF4};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_hdf4=$sys_hdf4_path.$inc;
      $lib_hdf4=$sys_hdf4_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the HDF4 path using\n";
      print "the LDT_HDF4 variable.\n";
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
   if(defined($ENV{LDT_HDF5})){
      $sys_hdf5_path = $ENV{LDT_HDF5};
      $inc = "/include/";
      $lib = "/lib/";
      $inc_hdf5=$sys_hdf5_path.$inc;
      $lib_hdf5=$sys_hdf5_path.$lib;
   }
   else {
      print "--------------ERROR---------------------\n";
      print "Please specify the HDF5 path using\n";
      print "the LDT_HDF5 variable.\n";
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
      if(defined($ENV{LDT_HDFEOS})){
         $sys_hdfeos_path = $ENV{LDT_HDFEOS};
         $inc = "/include/";
         $lib = "/lib/";
         $inc_hdfeos=$sys_hdfeos_path.$inc;
         $lib_hdfeos=$sys_hdfeos_path.$lib;
      }
      else {
         print "--------------ERROR---------------------\n";
         print "Please specify the HDFEOS path using\n";
         print "the LDT_HDFEOS variable.\n";
         print "Configuration exiting ....\n";
         print "--------------ERROR---------------------\n";
         exit 1;
      }
   }
}

print "Enable GeoTIFF support? (1-yes, 0-no, default=1): ";
$enable_geotiff=<stdin>;
if($enable_geotiff eq "\n"){
   $enable_geotiff=1;
}
if($enable_geotiff == 1) {
    if(defined($ENV{LDT_GDAL})){
	$sys_gdal_path = $ENV{LDT_GDAL};
	$lib = "/lib/";
	$lib_gdal=$sys_gdal_path.$lib;
    }
    else {
	print "--------------ERROR---------------------\n";
	print "Please specify the GDAL path using\n";
	print "the LDT_GDAL variable.\n";
	print "GDAL can be obtained from www.gdal.org\n";
	print "Configuration exiting ....\n";
	print "--------------ERROR---------------------\n";
	exit 1;
    }
    if(defined($ENV{LDT_FORTRANGIS})){
	$sys_fortrangis_path = $ENV{LDT_FORTRANGIS};
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
	print "the LDT_FORTRANGIS variable.\n";
	print "FORTRANGIS can be obtained from http://fortrangis.sourceforge.net\n";
	print "Configuration exiting ....\n";
	print "--------------ERROR---------------------\n";
	exit 1;
    }
}

# EMK...Add LIBGEOTIFF support for Air Force
print "Enable LIBGEOTIFF support? (1-yes, 0-no, default=1): ";
$enable_libgeotiff=<stdin>;
if($enable_libgeotiff eq "\n"){
    $enable_libgeotiff=1;
}
if($enable_libgeotiff == 1) {
    if(defined($ENV{LDT_LIBGEOTIFF})){
	$sys_libgeotiff_path = $ENV{LDT_LIBGEOTIFF};
	$lib = "/lib/";
	$lib_libgeotiff=$sys_libgeotiff_path.$lib;
	$include = "/include/";
	$inc_libgeotiff=$sys_libgeotiff_path.$include;
    }
    else {
	print "--------------ERROR---------------------\n";
	print "Please specify the LIBGEOTIFF path using\n";
	print "the LDT_LIBGEOTIFF variable.\n";
	print "LIBGEOTIFF can be obtained from https://github.com/OSGeo/libgeotiff\n";
	print "Configuration exiting ....\n";
	print "--------------ERROR---------------------\n";
	exit 1;
    }
    if ($cray_modifications == 1) {
       if(defined($ENV{LDT_LIBTIFF})){
          $sys_libtiff_path = $ENV{LDT_LIBTIFF};
       }
       else {
             print "--------------ERROR---------------------\n";
             print "Please specify the LIBTIFF path using\n";
             print "the LDT_LIBTIFF variable.\n";
             print "Configuration exiting ....\n";
             print "--------------ERROR---------------------\n";
             exit 1;
       }
       if(defined($ENV{LDT_LIBJBIG})){
          $sys_libjbig_path = $ENV{LDT_LIBJBIG};
       }
       else {
             print "--------------ERROR---------------------\n";
             print "Please specify the LIBJBIG path using\n";
             print "the LDT_LIBJBIG variable.\n";
             print "Configuration exiting ....\n";
             print "--------------ERROR---------------------\n";
             exit 1;
       }
       if(defined($ENV{LDT_LIBLZMA})){
          $sys_liblzma_path = $ENV{LDT_LIBLZMA};
       }
       else {
             print "--------------ERROR---------------------\n";
             print "Please specify the LIBLZMA path using\n";
             print "the LDT_LIBLZMA variable.\n";
             print "Configuration exiting ....\n";
             print "--------------ERROR---------------------\n";
             exit 1;
       }
    }
}

print "Include date/time stamp history? (1-yes, 0-no, default=1): ";
$use_history=<stdin>;
if($use_history eq "\n"){
   $use_history=1;
}


if(defined($ENV{LDT_JPEG})){
   $libjpeg = "-L".$ENV{LDT_JPEG}."/lib"." -ljpeg";
}
else{
   $libjpeg = "-ljpeg";
}

if(defined($ENV{LDT_RPC})){
   $librpc = $ENV{LDT_RPC};
}
else{
   $librpc = "";
}

if($sys_arch eq "linux_ifc") {
   if($use_omp == 1) {
      if($use_endian == 1) {
         $fflags77= "-c -openmp ".$sys_opt."-nomixed-str-len-arg -names lowercase -convert little_endian -assume byterecl ".$sys_par." -DIFC -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
         $fflags =" -c -openmp ".$sys_opt."-u -traceback -fpe0  -nomixed-str-len-arg -names lowercase -convert little_endian -assume byterecl ".$sys_par."-DIFC -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
         $ldflags= " -openmp -L\$(LIB_ESMF) -lesmf -lstdc++ -limf -lm -lrt -lz";
      }
      else {
         $fflags77= "-c -openmp ".$sys_opt."-nomixed-str-len-arg -names lowercase -convert big_endian -assume byterecl ".$sys_par." -DIFC -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
         $fflags =" -c -openmp ".$sys_opt."-u -traceback -fpe0  -nomixed-str-len-arg -names lowercase -convert big_endian -assume byterecl ".$sys_par."-DIFC -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
         $ldflags= " -openmp -L\$(LIB_ESMF) -lesmf -lstdc++ -limf -lm -lrt -lz";
      }
   }
   else {
      if($use_endian == 1) {
         $fflags77= "-c ".$sys_opt."-nomixed-str-len-arg -names lowercase -convert little_endian -assume byterecl ".$sys_par." -DIFC -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
         $fflags =" -c ".$sys_opt."-u -traceback -fpe0  -nomixed-str-len-arg -names lowercase -convert little_endian -assume byterecl ".$sys_par."-DIFC -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
         $ldflags= " -L\$(LIB_ESMF) -lesmf -lstdc++ -limf -lm -lrt -lz";
      }
      else {
         $fflags77= "-c ".$sys_opt."-nomixed-str-len-arg -names lowercase -convert big_endian -assume byterecl ".$sys_par." -DIFC -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
         $fflags =" -c ".$sys_opt."-u -traceback -fpe0  -nomixed-str-len-arg -names lowercase -convert big_endian -assume byterecl ".$sys_par."-DIFC -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
         $ldflags= " -L\$(LIB_ESMF) -lesmf -lstdc++ -limf -lm -lrt -lz";
      }
   }

   $cflags = "-c ".$sys_c_opt." -DIFC";

}
elsif($sys_arch eq "linux_pgi") {
   $cflags = "-c -DLITTLE_ENDIAN -DPGI";
   $fflags77= "-c ".$sys_opt."-C -s -Rb -Rs -g -gopt -Mbounds -Minform=inform -Minfo=all -DPGI -Mbyteswapio -r4 -i4 -Mpreprocess ".$sys_par."-I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
   $fflags ="-c ".$sys_opt."-C -s -Rb -Rs -g -gopt -Mbounds -Minform=inform -Minfo=all  -DPGI -Mbyteswapio -r4 -i4 -Mpreprocess ".$sys_par." -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
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
   $cflags = "-c -DGFORTRAN ";
   $fflags77= "-c -pass-exit-codes ".$sys_opt." ".$sys_par." ".$endian." -DGFORTRAN -I\$(MOD_ESMF) ";
   $fflags =" -c -pass-exit-codes -ffree-line-length-0 ".$sys_opt." ".$sys_par." ".$endian." -DGFORTRAN -I\$(MOD_ESMF) ";
   $ldflags= " -L\$(LIB_ESMF) -lesmf -lstdc++ -lz ";
}
elsif($sys_arch eq "AIX") {
   $cflags = "-c -w -g -qfullpath -q64 -qcpluscmt";
   $fflags77= "-c ".$sys_opt."-c -g -qkeepparm -qsuffix=f=f:cpp=F90 -q64 -WF,-DAIX, ".$sys_par." -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
   $fflags ="-c ".$sys_opt."-c -g -qkeepparm -qsuffix=f=f:cpp=F90 -q64 -WF,-DAIX, ".$sys_par." -I\$(MOD_ESMF) -DUSE_INCLUDE_MPI";
   $ldflags= "-q64 -bmap:map -bloadmap:lm -lmass  -L\$(LIB_ESMF) -lesmf -lstdc++ -limf -lm -lrt -lz";
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
   if (index($sys_fc, "-mt_mpi") != -1) {
      $ldflags = $ldflags." -lmpi_mt";
   }
   elsif ($cray_modifications == 1 || $sys_arch eq "cray_cray") {
      $ldflags = $ldflags." -lmpich";
   }
   elsif ($ibm_modifications == 1) {
      $ldflags = $ldflags." -lmpi_ibm_mpifh";
   }
   else{
      $ldflags = $ldflags." -lmpi";
   }
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

# EMK...Added LIGBEOTIFF for Air Force
if($enable_libgeotiff== 1){
    # Kluge for koehr; locally installed the tiff library.
    if(defined($ENV{LDT_TIFF})){
       $INC_LIBTIFF = $ENV{LDT_TIFF}."/include";
       $cflags = $cflags." -I".$INC_LIBTIFF." -I\$(INC_LIBGEOTIFF)";
    }
    else {
       $cflags = $cflags." -I\$(INC_LIBGEOTIFF)";
    }
    $tiffpath = "";
    $tiffdeps = "";
    # Kluge for koehr; locally installed the tiff library.
    if(defined($ENV{LDT_TIFF})){
       $tiffpath = " -L".$ENV{LDT_TIFF}."/lib";
    }
    if ( $cray_modifications == 1 ) {
       $tiffpath = " -L".$sys_libtiff_path;
       $tiffdeps = " -L".$sys_libjbig_path." -ljbig "." -L".$sys_liblzma_path." -llzma ";
    }
    $ldflags = $ldflags." -L\$(LIB_LIBGEOTIFF) ".$tiffpath." -ltiff -lgeotiff -lm -lz ".$libjpeg." ".$tiffdeps;
}


open(conf_file,">configure.ldt");
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
# EMK...Added LIBGEOTIFF for Air Force
printf conf_file "%s%s\n","INC_LIBGEOTIFF  = $inc_libgeotiff";
printf conf_file "%s%s\n","LIB_LIBGEOTIFF  = $lib_libgeotiff";

printf conf_file "%s%s\n","LIB_GDAL        = $lib_gdal";
printf conf_file "%s%s\n","CFLAGS          = $cflags";
printf conf_file "%s%s\n","FFLAGS77        = $fflags77";
printf conf_file "%s%s\n","FFLAGS          = $fflags";
printf conf_file "%s%s\n","LDFLAGS         = $ldflags";
close(conf_file);

print "-----------------------------------------------------\n";
print " configure.ldt file generated successfully\n";
print "-----------------------------------------------------\n";

open(misc_file,">LDT_misc.h");

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

if($enable_geotiff == 1) {
   printf misc_file "%s\n","#define USE_GDAL ";
}
else{
   printf misc_file "%s\n","#undef USE_GDAL ";
}

# EMK Added LIBGEOTIFF support for Air Force
if($enable_libgeotiff == 1) {
    printf misc_file "%s\n","#define USE_LIBGEOTIFF ";
}
else{
    printf misc_file "%s\n","#undef USE_LIBGEOTIFF ";
}


if($use_history == 1) {
   printf misc_file "%s\n","#undef LDT_SKIP_HISTORY ";
}
else{
   printf misc_file "%s\n","#define LDT_SKIP_HISTORY ";
}

#printf misc_file "%s\n","#undef SPMD";
printf misc_file "%s\n","#define BILHEADER_FILE_READ_METHOD_1_";
printf misc_file "%s\n","#define BILREAD_ASSUME_BIG_ENDIAN_";

close(misc_file);

open(netcdf_file,">LDT_NetCDF_inc.h");
printf netcdf_file "%s %d \n","#define NETCDF_shuffle ", $netcdf_shuffle;
printf netcdf_file "%s %d \n","#define NETCDF_deflate ", $netcdf_deflate;
printf netcdf_file "%s %d \n","#define NETCDF_deflate_level ", $netcdf_deflate_level;
close(netcdf_file);

