!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
#include "LVT_NetCDF_inc.h"
!BOP
! 
! !MODULE: LVT_historyMod
! \label(LVT_historyMod)
!
! !INTERFACE:
module LVT_historyMod
! 
! !USES:
  use LVT_histDataMod
  use LVT_logMod
  use LVT_coreMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  use grib_api

#if (defined SPMD)
#if (defined USE_INCLUDE_MPI)
  implicit none
  include 'mpif.h' 
#else
  use mpi
  implicit none
#endif
#endif

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  The code in this file provides interfaces to manage LIS model output 
!  in different file formats. 
!  
! \subsubsection{Overview}
! The module provides generic interfaces to manage output from different
! land surface models. Currently LIS supports three data formats. 
! \begin{itemize}
! \item{binary: in binary, big-endian.}
! \item{grib1: WMO grib1 format}
! \item{netcdf: NCAR unidata netcdf format}
! \end{itemize}
! This module also provides generic interfaces for writing and reading
! model restart files. LIS uses the following convention on the 
! filename extensions: 
! \begin{itemize}
! \item{.1gd4r: binary output, 1 correspondes to one variable, g represents
! gridded data, d represnts direct access, 4r represents 4 byte real}
! \item{.grb: files in grib1 format}
! \item{.nc: files in netcdf format}
! \item{*rst: restart files. Currently all restart files are written as
! binary format}
! \end{itemize}
! This module also manages any data aggregatiion from computing nodes
! when parallel processing is employed. 
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  10 Feb 2004: Sujay Kumar; Initial Specification
!   4 Jul 2008: Sujay Kumar; Redesigned with generic routines that can handle
!                   model output in different formats for all LSMs
!  13 Jul 2008: Sujay Kumar; Added support for GRIB2
! 
!EOP
!BOP
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_writevar_bin               ! write a variable into a binary file
  public :: LVT_writevar_gridded           ! write a variable in grid space to a file
  public :: LVT_writevar_restart
  public :: LVT_readvar_restart
  public :: LVT_writevar_data_header       ! writes header information into the output file
  public :: LVT_close_data_header          ! closes entries for header information
  public :: LVT_tile2grid
!  public :: LVT_voidLISModelData

!EOP

!BOP
! 
! !ROUTINE: LVT_writevar_bin
! \label{LVT_writevar_bin}
!
! !INTERFACE:
  interface LVT_writevar_bin
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! This interface provides routines for writing variables (real)
! in a binary file. The aggregation from decomposed tile spaces in 
! each processor to their respective grid spaces and the aggregations 
! of these grid spaces from each processor are also performed by 
! this routine. The interface also provides options to write some 
! diagnostic statistics about the variable being written. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP 
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure writevar_bin_real     
!     module procedure writevar_bin_real_direct
!EOP 
  end interface



!BOP
! 
! !ROUTINE: LVT_writevar_gridded
! \label{LVT_writevar_gridded}
!
! !INTERFACE:
  interface LVT_writevar_gridded
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This interface provides routines for writing variables (real)
!  in grid space to a file. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP 
! !PRIVATE MEMBER FUNCTIONS:
     module procedure writevar_gridded_tile
     module procedure writevar_gridded_1dgrid
     module procedure writevar_gridded_1dgrid_ensem
     module procedure writevar_gridded_1dgrid_ensem_nfield
!EOP
  end interface

!BOP
! 
! !ROUTINE: LVT_readvar_restart
! \label{LVT_readvar_restart}
!
! !INTERFACE:
  interface LVT_readvar_restart
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This interface provides routines for writing variables (real)
!  in grid space to a file. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP 
! !PRIVATE MEMBER FUNCTIONS:
     module procedure readvar_restart_real
     module procedure readvar_restart_tile_real
     module procedure readvar_restart_int
!EOP
  end interface

!BOP
! 
! !ROUTINE: LVT_writevar_restart
! \label{LVT_writevar_restart}
!
! !INTERFACE:
  interface LVT_writevar_restart
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This interface provides routines for writing variables (real)
!  in grid space to a file. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP 
! !PRIVATE MEMBER FUNCTIONS:
     module procedure writevar_restart_real
     module procedure writevar_restart_tile_real
     module procedure writevar_restart_int
!EOP
  end interface

contains

!BOP
! 
! !ROUTINE: LVT_writevar_data_header
! \label{LVT_writevar_data_header}
!
! !INTERFACE: 
  subroutine LVT_writevar_data_header(ftn, ftn_meta_out, short_name, &
       standard_name, long_name, units,  &
       varid,vlevels, nMetricLevs, entryno)
! 
! !USES:     
    implicit none
!
! !INPUT PARAMETERS: 
    integer,           intent(in) :: ftn
    integer,           intent(in) :: ftn_meta_out
    character(len=*),  intent(in) :: short_name
    character(len=*),  intent(in) :: standard_name
    character(len=*),  intent(in) :: long_name
    character(len=*),  intent(in) :: units
    integer                       :: varid
    integer                       :: vlevels
    integer                       :: nMetricLevs
    integer                       :: entryno

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine writes the global header information and the 
!   variable-specific header information into the output file.   
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    !integer               :: dimID(5),dimID_t(4)
    integer,save          :: dimID(5) ! EMK preserve between calls
    integer               :: dimID_t(4)
    integer               :: tdimID
    character*8             :: xtime_begin_date
    character*6             :: xtime_begin_time
    character*50            :: xtime_units
    character*50            :: xtime_timeInc
! Note that the fix to add lat/lon to the NETCDF output will output
! undefined values for the water points. 
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    integer, dimension(8) :: values
    integer :: t
    integer :: shuffle, deflate, deflate_level

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(trim(LVT_rc%lvt_out_format).eq."netcdf") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 

       if(entryno.eq.1) then 
          if(allocated(LVT_histData%xlat%value)) then 
             deallocate(LVT_histData%xlat%value)
          endif
          if(allocated(LVT_histData%xlat%count)) then 
             deallocate(LVT_histData%xlat%count)
          endif
          if(allocated(LVT_histData%xlat%unittypes)) then 
             deallocate(LVT_histData%xlat%unittypes)
          endif
          if(allocated(LVT_histData%xlon%value)) then 
             deallocate(LVT_histData%xlon%value)
          endif
          if(allocated(LVT_histData%xlon%count)) then 
             deallocate(LVT_histData%xlon%count)
          endif
          if(allocated(LVT_histData%xlon%unittypes)) then 
             deallocate(LVT_histData%xlon%unittypes)
          endif
          call date_and_time(date,time,zone,values)
          LVT_histData%xlat%short_name = "latitude"
          LVT_histData%xlat%long_name = "latitude"
          LVT_histData%xlat%standard_name = "latitude"
          LVT_histData%xlat%units = "degree_north"
          LVT_histData%xlat%nunits = 1
          LVT_histData%xlat%format = 'F'
          LVT_histData%xlat%vlevels = 1
          LVT_histData%xlat%timeAvgOpt = 0 
          LVT_histData%xlat%startNlevs = 1
          LVT_histData%xlat%endNlevs = 1
          if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             allocate(LVT_histData%xlat%value(LVT_LIS_rc(1)%ntiles,&
                  1,LVT_histData%xlat%vlevels))
          else
             allocate(LVT_histData%xlat%value(LVT_rc%ngrid,&
                  1,LVT_histData%xlat%vlevels))
          endif
          allocate(LVT_histData%xlat%count(1,&
               1,LVT_histData%xlat%vlevels))
          LVT_histData%xlat%count = 1
          allocate(LVT_histData%xlat%unittypes(1))
          LVT_histData%xlat%unittypes(1) = "degree_north"

          LVT_histData%xlon%short_name = "longitude"
          LVT_histData%xlon%long_name = "longitude"
          LVT_histData%xlon%standard_name = "longitude"
          LVT_histData%xlon%units = "degree_east"
          LVT_histData%xlon%nunits = 1
          LVT_histData%xlon%format = 'F'
          LVT_histData%xlon%vlevels = 1
          LVT_histData%xlon%timeAvgOpt = 0 
          LVT_histData%xlon%startNlevs = 1
          LVT_histData%xlon%endNlevs = 1
          if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             allocate(LVT_histData%xlon%value(LVT_LIS_rc(1)%ntiles,&
                  1,LVT_histData%xlat%vlevels))
          else
             allocate(LVT_histData%xlon%value(LVT_rc%ngrid,&
                  1,LVT_histData%xlon%vlevels))
          endif
          allocate(LVT_histData%xlon%count(1,&
               1,LVT_histData%xlon%vlevels))
          LVT_histData%xlon%count = 1
          allocate(LVT_histData%xlon%unittypes(1))
          LVT_histData%xlon%unittypes(1) = "degree_north"
          
          if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
             call LVT_verify(nf90_def_dim(ftn,'ntiles',LVT_LIS_rc(1)%ntiles,&
                  dimID(1)))
          elseif(LVT_rc%lvt_wopt.eq."2d gridspace") then 
             call LVT_verify(nf90_def_dim(ftn,'east_west',LVT_rc%gnc,dimID(1)))
             call LVT_verify(nf90_def_dim(ftn,'north_south',LVT_rc%gnr,dimID(2)))
          elseif(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
             call LVT_verify(nf90_def_dim(ftn,'east_west',LVT_rc%gnc,dimID(1)))
             call LVT_verify(nf90_def_dim(ftn,'north_south',LVT_rc%gnr,dimID(2)))
             call LVT_verify(nf90_def_dim(ftn,'ensemble',LVT_rc%nensem,dimID(3)))
          endif
!LVT output is writing output for a single time
          call LVT_verify(nf90_def_dim(ftn,'time',1,tdimID))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value",&
               LVT_rc%udef))
!write lat/lon information
          if(trim(LVT_rc%lvt_wopt).eq."1d tilespace") then
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlat%short_name),&
                  nf90_float,&
                  dimids = dimID(1:1), varID=LVT_histData%xlat%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlat%varid_def,&
                  shuffle,deflate,deflate_level))
#endif
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlon%short_name),&
                  nf90_float,&
                  dimids = dimID(1:1), varID=LVT_histData%xlon%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlon%varid_def,&
                  shuffle,deflate,deflate_level))
#endif             

          elseif(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlat%short_name),&
                  nf90_float,&
                  dimids = dimID(1:2), varID=LVT_histData%xlat%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlat%varid_def,&
                  shuffle,deflate,deflate_level))
#endif
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlon%short_name),&
                  nf90_float,&
                  dimids = dimID(1:2), varID=LVT_histData%xlon%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlon%varid_def,&
                  shuffle,deflate,deflate_level))
#endif
          elseif(trim(LVT_rc%lvt_wopt).eq."2d ensemble gridspace") then 
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlat%short_name),&
                  nf90_float,&
                  dimids = dimID(1:2), varID=LVT_histData%xlat%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlat%varid_def,&
                  shuffle,deflate,deflate_level))
#endif
             call LVT_verify(nf90_def_var(ftn,&
                  trim(LVT_histData%xlon%short_name),&
                  nf90_float,&
                  dimids = dimID(1:2), varID=LVT_histData%xlon%varid_def))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,&
                  LVT_histData%xlon%varid_def,&
                  shuffle,deflate,deflate_level))
#endif          
          endif

          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "units",trim(LVT_histData%xlat%units)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "standard_name",trim(LVT_histData%xlat%standard_name)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "long_name",trim(LVT_histData%xlat%long_name)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "scale_factor",1.0))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "add_offset",0.0))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "missing_value",LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlat%varid_def,&
               "_FillValue",LVT_rc%udef))

          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "units",trim(LVT_histData%xlon%units)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "standard_name",trim(LVT_histData%xlon%standard_name)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "long_name",trim(LVT_histData%xlon%long_name)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "scale_factor",1.0))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "add_offset",0.0))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "missing_value",LVT_rc%udef))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xlon%varid_def,&
               "_FillValue",LVT_rc%udef))


!define time field
          call LVT_verify(nf90_def_var(ftn,'time',&
               nf90_float,dimids = tdimID,varID=LVT_histData%xtimeID))
          write(xtime_units,200) LVT_rc%yr, LVT_rc%mo, LVT_rc%da, &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
200       format ('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
               I2.2,':',I2.2)
          write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
               LVT_rc%yr, LVT_rc%mo, LVT_rc%da
          write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
               LVT_rc%hr, LVT_rc%mn, LVT_rc%ss
          write(xtime_timeInc, fmt='(I20)') &
               LVT_rc%ts
          
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "units",trim(xtime_units)))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "long_name","time"))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "time_increment",trim(adjustl(xtime_timeInc))))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "begin_date",xtime_begin_date))
          call LVT_verify(nf90_put_att(ftn,LVT_histData%xtimeID,&
               "begin_time",xtime_begin_time))

          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
               "LVT land surface analysis output"))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
               trim(LVT_rc%institution)))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
               "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
               date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
               "Kumar_etal_GMD_2012"))
          call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
               "website: http://lis.gsfc.nasa.gov/"))

          !grid information
          if(trim(LVT_rc%domain).eq."latlon") then !latlon
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "EQUIDISTANT CYLINDRICAL"))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(9)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(10)))       
          elseif(trim(LVT_rc%domain).eq."mercator") then 
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "MERCATOR"))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
          elseif(trim(LVT_rc%domain).eq."lambert") then !lambert conformal
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "LAMBERT CONFORMAL"))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
             
          elseif(trim(LVT_rc%domain).eq."polar") then ! polar stereographic
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
                  "POLAR STEREOGRAPHIC"))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
          endif
       
          if(nMetricLevs.gt.1) then 
             call LVT_verify(nf90_def_dim(ftn,trim(short_name)//'_levels',&
                  nMetricLevs, dimID(5)))
          endif

       endif

       if(LVT_rc%lvt_wopt.eq."1d tilespace") then
          if(vlevels.gt.1) then 
             call LVT_verify(nf90_def_dim(ftn,trim(short_name)//'_profiles',&
                  vlevels, dimID(4)))
          endif
          if(vlevels.eq.1) then 
             call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                  dimids = dimID(1:1), varID=varId))
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                  shuffle,deflate,deflate_level))
#endif
          else
             dimID_t(1) = dimID(1)
             dimID_t(2) = dimID(4)
             call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                  dimids = dimID_t(1:2), varID=varId),&
                  'nf90_def_var failed in LVT_historyMod')
#if(defined USE_NETCDF4) 
             call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                  shuffle,deflate,deflate_level),&
                  'nf90_def_var_deflate failed in LVT_historyMod')
#endif
          endif

       elseif(LVT_rc%lvt_wopt.eq."2d gridspace") then 
          ! The following are the 3-d fields
          if(vlevels.gt.1) then 
             call LVT_verify(nf90_def_dim(ftn,trim(short_name)//'_profiles',&
                  vlevels, dimID(4)))
          endif
          if(vlevels.eq.1) then 
             if(nMetricLevs.eq.1) then              
                call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                     dimids = dimID(1:2), varID=varId))
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level))
#endif
             else
                dimID_t(1) = dimID(1)
                dimID_t(2) = dimID(2)
                dimID_t(3) = dimID(5)

                call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                     dimids = dimID_t(1:3), varID=varId))
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level))
#endif
             endif
          else
             if(nMetricLevs.eq.1) then
                dimID_t(1) = dimID(1)
                dimID_t(2) = dimID(2)
                dimID_t(3) = dimID(4)
                
                call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                     dimids = dimID_t(1:3), varID=varId),&
                     'nf90_def_var failed in LVT_historyMod')
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level),&
                     'nf90_def_var_deflate failed in LVT_historyMod')
#endif
             else
                dimID_t(1) = dimID(1)
                dimID_t(2) = dimID(2)
                dimID_t(3) = dimID(4)
                dimID_t(4) = dimID(5)
                
                call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                     dimids = dimID_t(1:4), varID=varId),&
                     'nf90_def_var failed in LVT_historyMod')
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level),&
                     'nf90_def_var_deflate failed in LVT_historyMod')
#endif                
             endif
          endif
       elseif(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
          ! The following are the 3-d fields
          if(vlevels.gt.1) then 
             call LVT_verify(nf90_def_dim(ftn,trim(short_name)//'_profiles',&
                  vlevels, dimID(4)))
          endif
          if(vlevels.eq.1) then 
             if(nMetricLevs.eq.1) then 
                call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                     dimids = dimID(1:3), varID=varId))
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level))
#endif
             else
                dimID_t(1) = dimID(1)
                dimID_t(2) = dimID(2)
                dimID_t(3) = dimID(3)
                dimID_t(4) = dimID(5)
                call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                     dimids = dimID_t(1:4), varID=varId))
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level))
#endif
             endif
          else
             if(nMetricLevs.eq.1) then
                call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                     dimids = dimID(1:4), varID=varId),&
                     'nf90_def_var failed in LVT_historyMod')
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level),&
                     'nf90_def_var_deflate failed in LVT_historyMod')
#endif
             else
                call LVT_verify(nf90_def_var(ftn,trim(short_name),nf90_float,&
                     dimids = dimID(1:5), varID=varId),&
                     'nf90_def_var failed in LVT_historyMod')
#if(defined USE_NETCDF4) 
                call LVT_verify(nf90_def_var_deflate(ftn,varID,&
                     shuffle,deflate,deflate_level),&
                     'nf90_def_var_deflate failed in LVT_historyMod')
#endif
             endif
          endif
       endif

       call LVT_verify(nf90_put_att(ftn,varID,&
            "units",trim(units)))
       call LVT_verify(nf90_put_att(ftn,varID,&
            "standard_name",trim(standard_name)))
       call LVT_verify(nf90_put_att(ftn,varID,&
            "long_name",trim(long_name)))
       call LVT_verify(nf90_put_att(ftn,varID,&
            "scale_factor",1.0))
       call LVT_verify(nf90_put_att(ftn,varID,&
            "add_offset",0.0))
       call LVT_verify(nf90_put_att(ftn,varID,&
            "missing_value",LVT_rc%udef))
       call LVT_verify(nf90_put_att(ftn,varID,&
            "_FillValue",LVT_rc%udef))
#endif 
    elseif(trim(LVT_rc%lvt_out_format).eq."binary") then !binary
       if(entryno.eq.1) then 
          write(ftn_meta_out,*) 'DIMENSIONS '
          write(ftn_meta_out,*) 'east_west ',LVT_rc%gnc
          write(ftn_meta_out,*) 'north_south ',LVT_rc%gnr
          if(trim(short_name).eq."SoilMoist") then 
             write(ftn_meta_out,*) "soilm_profiles", vlevels
          endif
          if(trim(short_name).eq."SoilTemp") then 
             write(ftn_meta_out,*) "soilt_profiles", vlevels
          endif

          write(ftn_meta_out,*) ' '
          write(ftn_meta_out,*) 'Missing value ',LVT_rc%udef
          write(ftn_meta_out,*) ' '
          write(ftn_meta_out,*) 'GRID INFORMATION '
          if(trim(LVT_rc%domain).eq."latlon") then !latlon
             write(ftn_meta_out,*) &
                  "MAP_PROJECTION: EQUIDISTANT CYLINDRICAL"
             write(ftn_meta_out,*) &
                  "SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)
             write(ftn_meta_out,*) &
                  "SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)
             write(ftn_meta_out,*) &
                  "DX", &
                  LVT_rc%gridDesc(9)
             write(ftn_meta_out,*) &
                  "DY", &
                  LVT_rc%gridDesc(10)
          elseif(trim(LVT_rc%domain).eq."mercator") then 
             write(ftn_meta_out,*) &
                  "MAP_PROJECTION: MERCATOR"
             write(ftn_meta_out,*) &
                  "SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)
             write(ftn_meta_out,*) &
                  "SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)
             write(ftn_meta_out,*) &
                  "TRUELAT1", &
                  LVT_rc%gridDesc(10)
             write(ftn_meta_out,*) &
                  "STANDARD_LON", &
                  LVT_rc%gridDesc(11)
             write(ftn_meta_out,*) &
                  "DX", &
                  LVT_rc%gridDesc(8)
             write(ftn_meta_out,*) &
                  "DY", &
                  LVT_rc%gridDesc(9)
          elseif(trim(LVT_rc%domain).eq."lambert") then !lambert conformal
             write(ftn_meta_out,*) &
                  "MAP_PROJECTION: LAMBERT CONFORMAL"
             write(ftn_meta_out,*) &
                  "SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)
             write(ftn_meta_out,*) &
                  "SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)
             write(ftn_meta_out,*) &
                  "TRUELAT1", &
                  LVT_rc%gridDesc(10)
             write(ftn_meta_out,*) &
                  "TRUELAT2", &
                  LVT_rc%gridDesc(7)
             write(ftn_meta_out,*) &
                  "STANDARD_LON", &
                  LVT_rc%gridDesc(11)
             write(ftn_meta_out,*) &
                  "DX", &
                  LVT_rc%gridDesc(8)
             write(ftn_meta_out,*) &
                  "DY", &
                  LVT_rc%gridDesc(9)
             
          elseif(trim(LVT_rc%domain).eq."polar") then ! polar stereographic
             write(ftn_meta_out,*) &
                  "MAP_PROJECTION POLAR STEREOGRAPHIC"
             write(ftn_meta_out,*) &
                  "SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)
             write(ftn_meta_out,*) &
                  "SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)
             write(ftn_meta_out,*) &
                  "TRUELAT1", &
                  LVT_rc%gridDesc(10)
             write(ftn_meta_out,*) &
                  "ORIENT", &
                  LVT_rc%gridDesc(7)
             write(ftn_meta_out,*) &
                  "STANDARD_LON", &
                  LVT_rc%gridDesc(11)
             write(ftn_meta_out,*) &
                  "DX", &
                  LVT_rc%gridDesc(8)
             write(ftn_meta_out,*) &
                  "DY", &
                  LVT_rc%gridDesc(9)
          endif
       endif
       write(ftn_meta_out,*) 'VARIABLE: ',trim(short_name),vlevels
       
    endif
  end subroutine LVT_writevar_data_header

  subroutine LVT_close_data_header(ftn)


    integer, intent(in) :: ftn 
    integer             :: c,r,t,gid,iret,col,row
    real                :: gtmp(LVT_rc%gnc,LVT_rc%gnr)

    if(trim(LVT_rc%lvt_out_format).eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)     
      call LVT_verify(nf90_enddef(ftn))

      call LVT_verify(nf90_put_var(ftn,LVT_histData%xtimeID,0.0))

      if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
         do t=1,LVT_LIS_rc(1)%ntiles
            col = LVT_LIS_domain(1)%tile(t)%col
            row = LVT_LIS_domain(1)%tile(t)%row
            gid = LVT_domain%gindex(col,row)
            LVT_histData%xlat%value(t,1,1) = LVT_domain%grid(gid)%lat
            LVT_histData%xlon%value(t,1,1) = LVT_domain%grid(gid)%lon
         enddo
      else
      ! set lat/lons 
         do t=1,LVT_rc%ngrid
            LVT_histData%xlat%value(t,1,1) = LVT_domain%grid(t)%lat
            LVT_histData%xlon%value(t,1,1) = LVT_domain%grid(t)%lon
         enddo
      endif

      gtmp = LVT_rc%udef
      do r=1,LVT_rc%gnr
         do c=1,LVT_rc%gnc
            gid = LVT_domain%gindex(c,r)
            if(gid.ne.-1) then 
               gtmp(c,r) = LVT_histData%xlat%value(gid,1,1)
            endif
         enddo
      enddo
      
      if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
         iret = nf90_put_var(ftn,LVT_histData%xlat%varid_def, &
              LVT_histData%xlat%value(:,1,1), (/1/),&
              (/LVT_LIS_rc(1)%ntiles/))
      else
         iret = nf90_put_var(ftn,LVT_histData%xlat%varid_def, &
              gtmp, (/1,1/),&
              (/LVT_rc%gnc,LVT_rc%gnr/))
      endif

      gtmp = LVT_rc%udef
      do r=1,LVT_rc%gnr
         do c=1,LVT_rc%gnc
            gid = LVT_domain%gindex(c,r)
            if(gid.ne.-1) then 
               gtmp(c,r) = LVT_histData%xlon%value(gid,1,1)
            endif
         enddo
      enddo
      
      if(LVT_rc%lvt_wopt.eq."1d tilespace") then 
         iret = nf90_put_var(ftn,LVT_histData%xlon%varid_def, &
              LVT_histData%xlon%value(:,1,1), (/1/),&
              (/LVT_LIS_rc(1)%ntiles/))
      else
         iret = nf90_put_var(ftn,LVT_histData%xlon%varid_def, &
              gtmp, (/1,1/),&
              (/LVT_rc%gnc,LVT_rc%gnr/))
      endif
#endif
   endif
 end subroutine LVT_close_data_header

!BOP
! 
! !ROUTINE: writevar_bin_real
! \label{writevar_bin_real}
!
! !INTERFACE:
! Private name: call using LVT_writevar_bin
  subroutine writevar_bin_real(ftn, var)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in) :: ftn
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  Write a real variable to a binary output file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    real                :: var(LVT_rc%ngrid)
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)

    if(trim(LVT_rc%lvt_wopt).eq."1d tilespace") then 
       
       print*, '1d tilespace of output is not supported in LVT'
       print*, 'program stopping in writevar_bin_real'
       stop
!       call gather_tiled_vector_output(gtmp1,var)

!       if ( LVT_masterproc ) then
!          write(ftn) gtmp1
!          deallocate(gtmp1)
!       endif

    elseif(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
       call gather_gridded_output(gtmp, var)

       if ( LVT_masterproc ) then
          write(ftn) gtmp
          deallocate(gtmp)
       endif
    elseif(trim(LVT_rc%lvt_wopt).eq."1d gridspace") then 
       call gather_gridded_vector_output(gtmp1, var)

       if ( LVT_masterproc ) then
          write(ftn) gtmp1
          deallocate(gtmp1)
       endif
    endif
  end subroutine writevar_bin_real


!BOP
! 
! !ROUTINE: writevar_gridded_1dgrid
! \label{writevar_gridded_1dgrid}
!
! !INTERFACE:
! Private name: call LVT_writevar_gridded
  subroutine writevar_gridded_1dgrid(ftn,  var, varid, dim1)
! 
! !USES:
    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Writes a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
    real                            :: var(LVT_rc%ngrid)
    integer                         :: varid
    integer,   intent(in), optional :: dim1
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)
    real, allocatable :: gtmp2(:)
    integer       :: l
    integer       :: iret
    integer       :: c,r,gid,ntiles, count1, ierr


    if(trim(LVT_rc%lvt_out_format).eq."binary") then
       if(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
          if(LVT_masterproc) then
             allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
             allocate(gtmp1(LVT_rc%glbngrid))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LVT_gdeltas(LVT_localPet),MPI_REAL,gtmp1,&
               LVT_gdeltas(:),LVT_goffsets(:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1 = var
#endif
          if(LVT_masterproc) then 
             gtmp = LVT_rc%udef
             count1=1
             do l=1,LVT_npes
                do r=LVT_nss_halo_ind(l),LVT_nse_halo_ind(l)
                   do c=LVT_ews_halo_ind(l),LVT_ewe_halo_ind(l)
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         gtmp(c,r) = gtmp1(gid)
                      endif
                   enddo
                enddo
             enddo
             write(ftn) gtmp
             deallocate(gtmp)
          endif
          deallocate(gtmp1)
       elseif(trim(LVT_rc%lvt_wopt).eq."1d gridspace") then 
          if(LVT_masterproc) then 
             allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
             allocate(gtmp1(LVT_rc%glbngrid))
             allocate(gtmp2(LVT_rc%glbngrid_red))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LVT_gdeltas(LVT_localPet),MPI_REAL,gtmp1,&
               LVT_gdeltas(:),LVT_goffsets(:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1 = var
#endif
          gtmp2 = gtmp1
             
          write(ftn) gtmp2
          deallocate(gtmp)
          deallocate(gtmp2)
          deallocate(gtmp1)
       endif
    elseif(trim(LVT_rc%lvt_out_format).eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
       if(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
          if(LVT_masterproc) then 
             allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
             allocate(gtmp1(LVT_rc%glbngrid))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LVT_gdeltas(LVT_localPet),MPI_REAL,gtmp1,&
               LVT_gdeltas(:),LVT_goffsets(:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1 = var
#endif
          if(LVT_masterproc) then 
             gtmp = LVT_rc%udef
             count1=1
             do l=1,LVT_npes
                do r=LVT_nss_halo_ind(l),LVT_nse_halo_ind(l)
                   do c=LVT_ews_halo_ind(l),LVT_ewe_halo_ind(l)
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         gtmp(c,r) = gtmp1(gid)
                      endif
                   enddo
                enddo
             enddo
             if(PRESENT(dim1)) then 
                iret = nf90_put_var(ftn,varid, gtmp, (/1,1,dim1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr,1/))
             else
                iret = nf90_put_var(ftn,varid, gtmp, (/1,1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr/))
             endif
             deallocate(gtmp)
          endif
          deallocate(gtmp1)
       elseif(trim(LVT_rc%lvt_wopt).eq."1d gridspace") then 
          if(LVT_masterproc) then 
             allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
             allocate(gtmp1(LVT_rc%glbngrid))
             allocate(gtmp2(LVT_rc%glbngrid_red))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1))
             gtmp1 = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LVT_gdeltas(LVT_localPet),MPI_REAL,gtmp1,&
               LVT_gdeltas(:),LVT_goffsets(:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1 = var
#endif
          gtmp2 = gtmp1
             
          iret = nf90_put_var(ftn,varid,gtmp2,(/1/),(/LVT_rc%glbngrid_red/))
          deallocate(gtmp)
          deallocate(gtmp2)
          deallocate(gtmp1)
       endif
#endif
    endif

  end subroutine writevar_gridded_1dgrid


!BOP
! 
! !ROUTINE: writevar_gridded_1dgrid_ensem
! \label{writevar_gridded_1dgrid_ensem}
!
! !INTERFACE:
! Private name: call LVT_writevar_gridded
  subroutine writevar_gridded_1dgrid_ensem(ftn,  var, varid, dim1)
! 
! !USES:
    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Writes a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
    real                            :: var(LVT_rc%ngrid,LVT_rc%nensem)
    integer                         :: varid
    integer,   intent(in), optional :: dim1
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:,:)
    real, allocatable :: gtmp2(:)
    real, allocatable :: gtmp3(:,:,:)

    real, allocatable :: gtmp1_1d(:)


    real          :: tmpval
    integer       :: igrib
    integer       :: count_tmpval
    integer       :: l,m
    integer       :: iret
    integer       :: c,r,gid,ntiles, count1, ierr


    if(trim(LVT_rc%lvt_out_format).eq."binary") then

       if(LVT_rc%nensem.gt.1) then 
          write(*,*) 'binary output for ensembles is currently not supported..'
          write(*,*) 'program stopping...'
          stop       
       endif
       if(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
          if(LVT_masterproc) then
             allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
             allocate(gtmp1_1d(LVT_rc%glbngrid))
             gtmp = 0.0
             gtmp1_1d = 0.0
          else
             allocate(gtmp1_1d(1))
             gtmp1_1d = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LVT_gdeltas(LVT_localPet),MPI_REAL,gtmp1_1d,&
               LVT_gdeltas(:),LVT_goffsets(:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1_1d = var(:,1)
#endif
          if(LVT_masterproc) then 
             gtmp = LVT_rc%udef
             count1=1
             do l=1,LVT_npes
                do r=LVT_nss_halo_ind(l),LVT_nse_halo_ind(l)
                   do c=LVT_ews_halo_ind(l),LVT_ewe_halo_ind(l)
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         gtmp(c,r) = gtmp1_1d(gid)
                      endif
                   enddo
                enddo
             enddo
             write(ftn) gtmp
             deallocate(gtmp)
          endif
          deallocate(gtmp1_1d)
       elseif(trim(LVT_rc%lvt_wopt).eq."1d gridspace") then 
          if(LVT_masterproc) then 
             allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
             allocate(gtmp1_1d(LVT_rc%glbngrid))
             allocate(gtmp2(LVT_rc%glbngrid_red))
             gtmp = 0.0
             gtmp1_1d = 0.0
          else
             allocate(gtmp1_1d(1))
             gtmp1_1d = 0.0
          endif
#if (defined SPMD)      
          call MPI_GATHERV(var,LVT_gdeltas(LVT_localPet),MPI_REAL,gtmp1_1d,&
               LVT_gdeltas(:),LVT_goffsets(:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
          gtmp1_1d = var(:,1)
#endif
          gtmp2 = gtmp1_1d
             
          write(ftn) gtmp2
          deallocate(gtmp)
          deallocate(gtmp2)
          deallocate(gtmp1_1d)
       endif
    elseif(trim(LVT_rc%lvt_out_format).eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
       if(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
          if(LVT_masterproc) then 
             allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
             allocate(gtmp1(LVT_rc%glbngrid,LVT_rc%nensem))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1,1))
             gtmp1 = 0.0
          endif
          do m=1,LVT_rc%nensem
#if (defined SPMD)      
             call MPI_GATHERV(var(:,m),LVT_gdeltas(LVT_localPet),MPI_REAL,gtmp1(:,m),&
                  LVT_gdeltas(:),LVT_goffsets(:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
             gtmp1(:,m) = var(:,m)
#endif
          enddo
          if(LVT_masterproc) then 
             gtmp = LVT_rc%udef
             count1=1
             do l=1,LVT_npes
                do r=LVT_nss_halo_ind(l),LVT_nse_halo_ind(l)
                   do c=LVT_ews_halo_ind(l),LVT_ewe_halo_ind(l)
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         tmpval = 0.0
                         count_tmpval = 0 
                         do m=1,LVT_rc%nensem
                            if(gtmp1(gid,m).ne.LVT_rc%udef) then 
                               tmpval = tmpval +gtmp1(gid,m)
                               count_tmpval = count_tmpval + 1
                            endif
                         enddo
                         if(count_tmpval.gt.0) then 
                            gtmp(c,r) = tmpval/count_tmpval
                         else
                            gtmp(c,r) = LVT_rc%udef
                         endif
                      endif
                   enddo
                enddo
             enddo
             if(PRESENT(dim1)) then 
                iret = nf90_put_var(ftn,varid, gtmp, (/1,1,dim1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr,1/))
             else
                iret = nf90_put_var(ftn,varid, gtmp, (/1,1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr/))
             endif
             deallocate(gtmp)
          endif
          deallocate(gtmp1)
       elseif(trim(LVT_rc%lvt_wopt).eq."2d ensemble gridspace") then 
          if(LVT_masterproc) then 
             allocate(gtmp3(LVT_rc%gnc,LVT_rc%gnr,LVT_rc%nensem))
             allocate(gtmp1(LVT_rc%glbngrid,LVT_rc%nensem))
             gtmp3 = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp3(1,1,1))
             gtmp1 = 0.0
          endif
          do m=1,LVT_rc%nensem
#if (defined SPMD)      
             call MPI_GATHERV(var(:,m),LVT_gdeltas(LVT_localPet),MPI_REAL,gtmp1(:,m),&
                  LVT_gdeltas(:),LVT_goffsets(:),MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
             gtmp1(:,m) = var(:,m)
#endif
          enddo
          if(LVT_masterproc) then 
             gtmp3 = LVT_rc%udef
             count1=1
             do l=1,LVT_npes
                do r=LVT_nss_halo_ind(l),LVT_nse_halo_ind(l)
                   do c=LVT_ews_halo_ind(l),LVT_ewe_halo_ind(l)
                      gid = LVT_domain%gindex(c,r)
                      if(gid.ne.-1) then 
                         gtmp3(c,r,:) = gtmp1(gid,:)
                      endif
                   enddo
                enddo
             enddo
             if(PRESENT(dim1)) then 
                iret = nf90_put_var(ftn,varid, gtmp3, (/1,1,1,dim1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr,LVT_rc%nensem,1/))

             else
                iret = nf90_put_var(ftn,varid, gtmp3, (/1,1,1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr,LVT_rc%nensem/))
             endif
             deallocate(gtmp3)
          endif
          deallocate(gtmp1)
       
       elseif(trim(LVT_rc%lvt_wopt).eq."1d gridspace") then 
          write(*,*) '1d gridspace style currently not supported..'
          write(*,*) 'program stopping...'
          stop
       endif
#endif
    endif

  end subroutine writevar_gridded_1dgrid_ensem


!BOP
! 
! !ROUTINE: writevar_gridded_1dgrid_ensem_nfield
! \label{writevar_gridded_1dgrid_ensem_nfield}
!
! !INTERFACE:
! Private name: call LVT_writevar_gridded
  subroutine writevar_gridded_1dgrid_ensem_nfield(ftn,  nMetricLevs, var, varid, dim1)
! 
! !USES:
    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
    integer, intent(in)             :: nMetricLevs
    real                            :: var(LVT_rc%ngrid,nMetricLevs,LVT_rc%nensem)
    integer                         :: varid
    integer,   intent(in), optional :: dim1

! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Writes a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    real, allocatable :: gtmp(:,:,:)
    real, allocatable :: gtmp1(:,:,:)
    real, allocatable :: gtmp2(:)
    real, allocatable :: gtmp3(:,:,:,:)

    real, allocatable :: gtmp1_1d(:)


    real          :: tmpval
    integer       :: igrib
    integer       :: count_tmpval
    integer       :: l,m,k
    integer       :: iret
    integer       :: c,r,gid,ntiles, count1, ierr


    if(trim(LVT_rc%lvt_out_format).eq."binary") then

       write(*,*) 'binary output for ensembles is currently not supported..'
       write(*,*) 'program stopping...'
       stop       

    elseif(trim(LVT_rc%lvt_out_format).eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
       if(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
          if(LVT_masterproc) then 
             allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr,nMetricLevs))
             allocate(gtmp1(LVT_rc%glbngrid,nMetricLevs,LVT_rc%nensem))
             gtmp = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp1(1,1,1))
             gtmp1 = 0.0
          endif
          do k=1,nMetricLevs
             do m=1,LVT_rc%nensem
#if (defined SPMD)      
                call MPI_GATHERV(var(:,k,m),LVT_gdeltas(LVT_localPet),&
                     MPI_REAL,gtmp1(:,k,m),&
                     LVT_gdeltas(:),LVT_goffsets(:),&
                     MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
                gtmp1(:,k,m) = var(:,k,m)
#endif
             enddo
          enddo
          if(LVT_masterproc) then 
             gtmp = LVT_rc%udef
             do k=1,nMetricLevs
                count1=1
                do l=1,LVT_npes
                   do r=LVT_nss_halo_ind(l),LVT_nse_halo_ind(l)
                      do c=LVT_ews_halo_ind(l),LVT_ewe_halo_ind(l)
                         gid = LVT_domain%gindex(c,r)
                         if(gid.ne.-1) then 
                            tmpval = 0.0
                            count_tmpval = 0 
                            do m=1,LVT_rc%nensem
                               if(gtmp1(gid,k,m).ne.LVT_rc%udef) then 
                                  tmpval = tmpval +gtmp1(gid,k,m)
                                  count_tmpval = count_tmpval + 1
                               endif
                            enddo
                            if(count_tmpval.gt.0) then 
                               gtmp(c,r,k) = tmpval/count_tmpval
                            else
                               gtmp(c,r,k) = LVT_rc%udef
                            endif
                         endif
                      enddo
                   enddo
                enddo
             enddo
             if(PRESENT(dim1)) then 
                iret = nf90_put_var(ftn,varid, gtmp, (/1,1,1,dim1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr,nMetricLevs,1/))
             else
                iret = nf90_put_var(ftn,varid, gtmp, (/1,1,1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr,nMetricLevs/))
             endif
             deallocate(gtmp)
          endif
          deallocate(gtmp1)
       elseif(trim(LVT_rc%lvt_wopt).eq."2d ensemble gridspace") then 
          if(LVT_masterproc) then 
             allocate(gtmp3(LVT_rc%gnc,LVT_rc%gnr,nMetricLevs,LVT_rc%nensem))
             allocate(gtmp1(LVT_rc%glbngrid,nMetricLevs,LVT_rc%nensem))
             gtmp3 = 0.0
             gtmp1 = 0.0
          else
             allocate(gtmp3(1,1,1,1))
             gtmp1 = 0.0
          endif
          do k=1,nMetricLevs
             do m=1,LVT_rc%nensem
#if (defined SPMD)      
                call MPI_GATHERV(var(:,k,m),LVT_gdeltas(LVT_localPet),&
                     MPI_REAL,gtmp1(:,k,m),&
                     LVT_gdeltas(:),LVT_goffsets(:),&
                     MPI_REAL,0,MPI_COMM_WORLD,ierr)
#else 
                gtmp1(:,k,m) = var(:,k,m)
#endif
             enddo
          enddo

          if(LVT_masterproc) then 
             gtmp3 = LVT_rc%udef
             do k=1,nMetricLevs
                count1=1
                do l=1,LVT_npes
                   do r=LVT_nss_halo_ind(l),LVT_nse_halo_ind(l)
                      do c=LVT_ews_halo_ind(l),LVT_ewe_halo_ind(l)
                         gid = LVT_domain%gindex(c,r)
                         if(gid.ne.-1) then 
                            gtmp3(c,r,k,:) = gtmp1(gid,k,:)
                         endif
                      enddo
                   enddo
                enddo
             enddo

             if(PRESENT(dim1)) then 
                iret = nf90_put_var(ftn,varid, gtmp3, (/1,1,1,1,dim1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr,nMetricLevs,LVT_rc%nensem,1/))

             else
                iret = nf90_put_var(ftn,varid, gtmp3, (/1,1,1,1/),&
                     (/LVT_rc%gnc,LVT_rc%gnr,nMetricLevs,LVT_rc%nensem/))
             endif
             deallocate(gtmp3)
          endif
          deallocate(gtmp1)
       
       elseif(trim(LVT_rc%lvt_wopt).eq."1d gridspace") then 
          write(*,*) '1d gridspace style currently not supported..'
          write(*,*) 'program stopping...'
          stop
       endif
#endif
    endif

  end subroutine writevar_gridded_1dgrid_ensem_nfield


!BOP
! 
! !ROUTINE: writevar_gridded_tile
! \label{writevar_gridded_tile}
!
! !INTERFACE:
! Private name: call LVT_writevar_gridded
  subroutine writevar_gridded_tile(ftn, var, varid, dummy, dim1)
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Writes a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 
    real                            :: var(LVT_LIS_rc(1)%ntiles)
    integer                         :: varid
    real                            :: dummy
    integer,   intent(in), optional :: dim1
!EOP
    real, allocatable :: gtmp(:,:)
    real, allocatable :: gtmp1(:)
    real, allocatable :: gtmp2(:)
    integer       :: l
    integer       :: iret
    integer       :: c,r,gid,ntiles, count1, ierr


    if(trim(LVT_rc%lvt_out_format).eq."binary") then
       if(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
          allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
          call LVT_tile2grid(gtmp,var)

          write(ftn) gtmp
          deallocate(gtmp)
       elseif(trim(LVT_rc%lvt_wopt).eq."1d tilespace") then 
          write(ftn) var          
       endif
    elseif(trim(LVT_rc%lvt_out_format).eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
       if(trim(LVT_rc%lvt_wopt).eq."2d gridspace") then 
          allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
          
          call LVT_tile2grid(gtmp,var)
          if(PRESENT(dim1)) then 
             iret = nf90_put_var(ftn,varid, gtmp, (/1,1,dim1/),&
                  (/LVT_rc%gnc,LVT_rc%gnr,1/))
          else
             iret = nf90_put_var(ftn,varid, gtmp, (/1,1/),&
                  (/LVT_rc%gnc,LVT_rc%gnr/))
          endif
          deallocate(gtmp)
       elseif(LVT_rc%lvt_wopt.eq."1d tilespace") then 
          if(PRESENT(dim1)) then 
             iret = nf90_put_var(ftn,varid, var, (/1,dim1/),&
                  (/LVT_LIS_rc(1)%ntiles,1/))
          else
             iret = nf90_put_var(ftn,varid, var, (/1/),&
                  (/LVT_LIS_rc(1)%ntiles/))
          endif
       elseif(trim(LVT_rc%lvt_wopt).eq."1d gridspace") then 

          print*, 'output style not supported ..'
          stop
       endif
#endif
    endif

  end subroutine writevar_gridded_tile

!BOP
! 
! !ROUTINE: gather_gridded_output
! \label{gather_gridded_output}
!
! !INTERFACE:
subroutine gather_gridded_output(gtmp, var)
! 
! !USES: 

   implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! This routine gathers the output data into a gridded array.
!
! This process aggregates tiles to their grid, accounting for ensemble runs.
!
! This process accounts for the halo.
!
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the gridded output data
!   \item [var]
!     output data to process
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  29 Oct 2008:  James Geiger; Initial code
! 
!EOP
!BOP
! !ARGUMENTS: 
   real, allocatable :: gtmp(:,:)
   real          :: var(LVT_rc%ngrid)

!EOP

   real :: var1(LVT_rc%ngrid)
   real, allocatable :: gtmp1(:)
   integer :: i,c,r,m,t,l
   integer :: count1,gid,ntiles,ierr

   if ( LVT_masterproc ) then 
      allocate(gtmp(LVT_rc%gnc,LVT_rc%gnr))
      allocate(gtmp1(LVT_rc%glbngrid))
      gtmp = LVT_rc%udef
      gtmp1 = 0.0
   else
      allocate(gtmp1(1))
      gtmp1 = 0.0
   endif

   
   gtmp1 = var
   count1=1

   do l=1,LVT_npes
      do r=LVT_nss_halo_ind(l),LVT_nse_halo_ind(l)
         do c=LVT_ews_halo_ind(l),LVT_ewe_halo_ind(l)
            gid = c+(r-1)*LVT_rc%gnc
            ntiles = LVT_LIS_domain(1)%ntiles_pergrid(gid)
            if(ntiles.ne.0) then                          
               if(r.ge.LVT_nss_ind(l).and.&
                    r.le.LVT_nse_ind(l).and.&
                    c.ge.LVT_ews_ind(l).and.&
                    c.le.LVT_ewe_ind(l))then !points not in halo
                  gtmp(c,r) = gtmp1(count1)
               endif
               count1 = count1 + 1
            endif
         enddo
      enddo
   enddo
   deallocate(gtmp1)
end subroutine gather_gridded_output


!BOP
! 
! !ROUTINE: gather_gridded_vector_output
! \label{gather_gridded_vector_output}
!
! !INTERFACE:
subroutine gather_gridded_vector_output( gtmp, var)
! 
! !USES: 

   implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! This routine gathers the output data into a gridded 1d array.
!
! This process aggregates tiles to their grid, accounting for ensemble runs.
!
! This process accounts for the halo.
!    
!  \begin{description}
!   \item [n]
!     index of the current nest
!   \item [gtmp]
!     return array for the gridded output data
!   \item [var]
!     output data to process
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  30 Jan 2009: Sujay Kumar, Initial Code
! 
!EOP
! !ARGUMENTS: 
   real, allocatable :: gtmp(:)
   real          :: var(LVT_rc%ngrid)
!EOP

   integer :: i,c,r,m,t,l
   integer :: ntiles, gid, count1
   integer :: ierr

   allocate(gtmp(LVT_rc%glbngrid))
   gtmp = var
  
 end subroutine gather_gridded_vector_output

!BOP
! 
! !ROUTINE: LVT_tile2grid
! \label{LVT_tile2grid}
!
! !INTERFACE:
  subroutine LVT_tile2grid(gvar,tvar)
! 
! !USES:


    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This routine converts a tile space variable to the corresponding
!  grid space. The aggregation involves weighted average of each tile
!  in a grid cell based on the vegetation distribution. 
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [tvar]
!     variable dimensioned in the tile space. 
!   \item [gvar]
!     variable after converstion to the grid space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!
! !ARGUMENTS:     
    real              :: gvar(LVT_rc%lnc,LVT_rc%lnr)
    real, intent(in)  :: tvar(LVT_LIS_rc(1)%ntiles)
!
!
!EOP
    integer           :: i,c,r,m,t
    integer           :: gcount(LVT_rc%lnc,LVT_rc%lnr)

    gvar = 0.0
    gcount = 0
    do i=1,LVT_LIS_rc(1)%ntiles,LVT_LIS_rc(1)%nensem
       do m=1,LVT_LIS_rc(1)%nensem
          t = i+m-1
          c = LVT_LIS_domain(1)%tile(t)%col
          r = LVT_LIS_domain(1)%tile(t)%row
          gvar(c,r) = gvar(c,r)+&
               tvar(t)*LVT_LIS_domain(1)%tile(t)%fgrd/&
               LVT_LIS_rc(1)%nensem
          gcount(c,r) = gcount(c,r) + 1
       enddo
    enddo

    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(gcount(c,r).eq.0) then 
             gvar(c,r) = LVT_rc%udef
          endif
       enddo
    enddo

  end subroutine LVT_tile2grid




!BOP
! 
! !ROUTINE: LVT_tavgLISModelData
! \label{LVT_tavgLISModelData}
!
! !INTERFACE: 
  subroutine LVT_tavgLISModelData
! 
! !USES:
    use LVT_statsDataMod,only : LVT_stats
    implicit none 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine invokes the calls to compute temporal averages of 
!   desired set of LIS output variables, based on the specified 
!   temporal averaging frequency
!  
!   The routines invoked are: 
!   \begin{description}
!    \item[tavgSingleVar](\ref{tavgSingleVar})
!     computes the temporal average for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer      :: index
    integer      :: gid, nl, k

    type(LVT_metadataEntry), pointer :: dataEntry
    type(LVT_metadataEntry), pointer :: ebal, swnet,lwnet,qle,qh,qg
    type(LVT_metadataEntry), pointer :: qf,qa,qv,delsurfheat
    type(LVT_metadataEntry), pointer :: wbal,rainf,snowf,qs,qsb,runoff
    type(LVT_metadataEntry), pointer :: delswe,delintercept,delsoilmoist
    type(LVT_metadataEntry), pointer :: delcoldcont, delsurfstor
    type(LVT_metadataEntry), pointer :: evap, potevap, wrsi
    type(LVT_metadataEntry), pointer :: br, ef,et
    type(LVT_metadataEntry), pointer :: sweoverp, swe, precip
    type(LVT_metadataEntry), pointer :: etoverp,qsoverp, qsboverp
    type(LVT_metadataEntry), pointer :: rootmoist, soilmoist
    type(LVT_metadataEntry), pointer :: roottemp, soiltemp
    type(LVT_metadataEntry), pointer :: watertabled, tws, gws, wt


  end subroutine LVT_tavgLISModelData


!BOP
! 
! !ROUTINE: writevar_restart_real
! \label{writevar_restart_real}
!
! !INTERFACE:
! Private name: call LVT_writevar_gridded
  subroutine writevar_restart_real(ftn,  var)
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
    real                            :: var(LVT_rc%ngrid)
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Writes a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 

    write(ftn) var

  end subroutine writevar_restart_real


!BOP
! 
! !ROUTINE: writevar_restart_tile_real
! \label{writevar_restart_tile_real}
!
! !INTERFACE:
! Private name: call LVT_writevar_gridded
  subroutine writevar_restart_tile_real(ftn,  var, tileflag)
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
    real                            :: var(LVT_LIS_rc(1)%ntiles)
    integer                         :: tileflag
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Writes a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 

    write(ftn) var

  end subroutine writevar_restart_tile_real

!BOP
! 
! !ROUTINE: writevar_restart_int
! \label{writevar_restart_int}
!
! !INTERFACE:
! Private name: call LVT_writevar_gridded
  subroutine writevar_restart_int(ftn,  var)
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
    integer                         :: var(LVT_rc%ngrid)
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Writes a int variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 

    write(ftn) var

  end subroutine writevar_restart_int


!BOP
! 
! !ROUTINE: readvar_restart_real
! \label{readvar_restart_real}
!
! !INTERFACE:
! Private name: call LVT_readvar_gridded
  subroutine readvar_restart_real(ftn,  var)
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
    real                            :: var(LVT_rc%ngrid)
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Reads a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 

    read(ftn) var

  end subroutine readvar_restart_real

!BOP
! 
! !ROUTINE: readvar_restart_tile_real
! \label{readvar_restart_tile_real}
!
! !INTERFACE:
! Private name: call LVT_readvar_gridded
  subroutine readvar_restart_tile_real(ftn,  var, tileflag)
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
    real                            :: var(LVT_LIS_rc(1)%ntiles)
    integer                         :: tileflag
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Reads a real variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 

    read(ftn) var

  end subroutine readvar_restart_tile_real


!BOP
! 
! !ROUTINE: readvar_restart_int
! \label{readvar_restart_int}
!
! !INTERFACE:
! Private name: call LVT_readvar_gridded
  subroutine readvar_restart_int(ftn,  var)
! 
! !USES:

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)             :: ftn
    integer                         :: var(LVT_rc%ngrid)
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  Reads a int variable to a binary 
!  sequential access, gridded file as either 2-dimensional/1-dimensional gridded field. 
!  (The 1-d gridded field consists of a vector of valid land points)
!  After reading the global data, the routine subroutine subsets
!  the data for each processor's domain.
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variables being written, dimensioned in the tile space
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! 
! !ARGUMENTS: 

    read(ftn) var

  end subroutine readvar_restart_int

end module LVT_historyMod
