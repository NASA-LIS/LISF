!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: SMOPSsm_obsMod
! \label(SMOPSsm_obsMod)
! !DESCRIPTION: adopted from LDT SMOPS reader
! This module handles the observation plugin for the 
!  NOAA/NESDIS-based ASCAT L2 Soil Moisture (SM) Files (aka, SMOPS)
!
! A. Background for SMOPS Operational SM Files 
! (from http://www.ospo.noaa.gov/Products/land/smops/index.html)
!
! The NOAA/NESDIS Office of Satellite and Prodcut Operations (OSPO)  
! produces the Soil Moisture Operational Products System (SMOPS), which 
! combines soil moisture retrievals from multi-satellites/sensors to provide 
! a global soil moisture map with more spatial and temporal coverage. The SMOPS 
! first retrieves soil moisture from AMSR-E* on board NASA EOS-Aqua, and then 
! combines its baseline retrievals with those from other available satellites/
! sensors, including ASCAT and SMOS, to improve the spatial and
! temporal coverage of the AMSR-E observations. (*SMOPS will get the soil
! moisture retrieval from AMSR-2 onboard GCOM-W in future).
! 
! The global soil moisture maps are generated in 6-hourly and daily intervals 
! with the latest 6 and 24 hours worth of soil moisture retrievals from 
! multi-satellites/algorithms, and mapped with a cylindrical projection on 
! 0.25x0.25 degree grids. For each grid point of the map, the output includes 
! soil moisture values (%vol/vol) of the surface (top 1-5 cm) soil layer with 
! associated quality information and metadata. The 6-hourly product is available
! in GRIB2 format at standard forecast times (00Z, 06Z, 12Z and 18Z), and daily 
! product is available in both GRIB2 and netCDF4 formats.
!
! The SMOPS ASCAT layer is the just the gridded original ASCAT soil moisture 
! values.  They only do CDF matching for the blended layer in the product.
!
! There are three versions of the SMOPS datasets.
!
! According to the use by the 557th Weather Wing:
!
!                          version_1.3 <  2016-10-31T12:00:00
!   2016-10-31T12:00:00 <= version_2.0 <  2017-08-24T12:00:00
!                          version_3.0 >= 2017-08-24T12:00:00
!
! Also, NESDIS has generated SMOPS version 3.0 datasets starting
! from 2012-08-01.
!
! For SMOPS versions 1.3 and 2.0, LDT processes only ASCAT A and B soil moisture
! observations.
!
! ____________________________________________________________________________
!
! 1. Real-time SMOPS Files from OSPO
!
!    There is no stand-alone, operational gridded ASCAT L2 soil moisture
! data product available, although there is a data layer in the operational
! SMOPS products file for the ASCAT soil moisture data.  Users can get the
! SMOPS data file, and extract the gridded ASCAT data.  For operations, 
! only about 7 days worth of data are kept available in the operation.
!
! The SMOPS data files carry the following name convention:
!   NPR_SMOPS_CMAP_Dyyyymmddhh.gr2 = 6 hour Product   (grib2 format)
!   NPR_SMOPS_CMAP_Dyyyymmdd.gr2   = Daily Product   (grib2 format)
!   NPR_SMOPS_CMAP_Dyyyymmdd.nc4   = Archive Product (netcdf4 format)
!
! The product's time latency is about 2.5 ~ 8.5 hours.  
!  (personal communication with Limin Zhao at NESDIS/OSPO)
!
! Here is webpage for the operational SMOSP products.  Test data 
! files available there, if you like to check it out.
!   http://www.ospo.noaa.gov/Products/land/smops/index.html
!
!
! 2. Derived Archive of SMOPS Files from NESDIS/STAR
!  - With help from Jicheng Liu (jicheng.liu@noaa.gov)
!  -  and Jerry Zhan (xiwu.zhan@noaa.gov) 
!
!  Files downloaded from:
!  ftp://ftp.star.nesdis.noaa.gov/pub/smcd/emb/jliu/SMOPS/
!
!  Years available:  2007-2012
!  Data format:  Grib-2
!
!   The data fields are not in the current official GRIB Table yet. To be 
! able to see meaningful variable/layer names, update your GRIB2 Table locally 
! and recompile the wgrib2 tool with the gribtab_SMOPS.dat file provided
! in this directory. 
!
! As you can see from this table, there are 7 soil moisture layers (layer 
! 1-7, with 2 of them as spare layers so far) and hour, minute and QA layers for 
! each soil moisture layer, which make the total number of layers 28.
! The parameter number in this table is sufficient enough for you to read 
! the individual layers. 
!
!   NESDIS/STAR researchers ran SMOPS on the historical ASCAT data back to 2007. 
! The only modification made to SMOPS was to make it take ASCAT SM only files.
! The original files are in BUFR format and downloaded from EUMETSAT (2007-2011) 
! and from KNMI (2011-present).  STAR researchers checked the files from the
! two sources and found they basically are the same. With that, they  produce 
! what's exactly coming out from OSPO SMOPS now in terms of data format.  At OSPO,  
! the grib2 format is the near real time format with the data latency of 5.5 hours 
! while the NetCDF format files are produced two days later to maximize the input 
! data availability for archiving purposes.   
!
!   For QC flags, the current ASCAT layer only carries two bytes of QA info 
! from the original data, which are "Estimated Error in Soil moisture" and 
! "Soil Moisture Quality". It's hard to squeeze anything else in unless we 
! change the output format from current SMOPS. If others are a must for 
! you, we may temporarily squeeze them in for the time this rerun but 
! later on, you'll have different ones from OSPO unless we make a code 
! updating at OSPO with output format change, which is always a major 
! issue here.
!
! !INTERFACE:
module SMOPSsm_obsMod
! 
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the Land Parameter
!  Retrieval Model (LPRM) AMSR-E soil moisture product
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  12 Dec 2014: Sujay Kumar, Initial Specification
!  30 May 2018 : Mahdi Navari, Updated to 
!                support reading ASCAT Metop A and B & real time data 
!                support reading different version of the SMOPS data 
!                support binary QC flag as implemented in the LIS and LDT 
! 
!EOP
! 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOPSsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOPSsmobs
!EOP
  type, public :: smopssmobsdec ! rtsmopsl2smobsdec
    integer                :: useRealtime
     character*100          :: odir
     integer                :: useASCAT
!     integer                :: useWindSat
     integer                :: useSMOS
     integer                :: useAMSR2
     integer                :: useSMAP
     integer                :: mo
     logical                :: startmode
     integer                :: smopsnc, smopsnr
     real*8                 :: version2_time, version3_time
!     real,    allocatable   :: smobs(:,:)
     character(len=17)      :: version
     type(proj_info)        :: smopsproj 
     integer, allocatable       :: n11(:)
     real,  allocatable         :: rlat(:)
     real,  allocatable         :: rlon(:)
     real,  allocatable         :: maxv(:,:)
     real,  allocatable         :: minv(:,:)
  end type smopssmobsdec 

  type(smopssmobsdec), allocatable:: SMOPSsmobs(:)
contains
  
!BOP
! 
! !ROUTINE: SMOPSsm_obsInit
! \label{SMOPSsm_obsInit}
!
! !INTERFACE: 
  subroutine SMOPSsm_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading LPRM AMSRE soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20) ! real                  :: gridDesci(50)
    integer                 :: n 
    integer                 :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real                    :: upgmt
    


    if(.not.allocated(SMOPSsmobs)) then 
       allocate(SMOPSsmobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%odir, &
         label='SMOPS soil moisture observation directory:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture observation directory: not defined')


    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useASCAT, &
         label='SMOPS soil moisture use ASCAT data:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture use ASCAT data: not defined')


    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useRealtime, &
         label='SMOPS ASCAT use realtime data:', rc=status)
    call LVT_verify(status, &
         'SMOPS ASCAT use realtime data: not defined')

!    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useWindSat, &
!         label='SMOPS soil moisture use WindSat data:', rc=status)
!    call LVT_verify(status, &
!        'SMOPS soil moisture use WindSat data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useSMOS, &
         label='SMOPS soil moisture use SMOS data:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture use SMOS data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useAMSR2, &
         label='SMOPS soil moisture use AMSR2 data:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture use AMSR2 data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, SMOPSsmobs(i)%useSMAP, &
         label='SMOPS soil moisture use SMAP data:', rc=status)
    call LVT_verify(status, &
         'SMOPS soil moisture use SMAP data: not defined')

    call ESMF_ConfigFindLabel(LVT_config, "SMOPS soil moisture version:", rc=status)
    call ESMF_ConfigGetAttribute(LVT_config,SMOPSsmobs(i)%version,&
         default='date-based', rc=status)

    !call LVT_update_timestep(LVT_rc, 86400)   !    ????????

    if(SMOPSsmobs(i)%useASCAT + SMOPSsmobs(i)%useAMSR2 + & !+ SMOPSsmobs(i)%useWindSat + &
         SMOPSsmobs(i)%useSMOS + SMOPSsmobs(i)%useSMAP.gt.1) then
       write(LVT_logunit,*) '[ERR] Please do not select multiple sensor sources'
       write(LVT_logunit,*) '[ERR] simultaneously for LVT preprocessing..'
       write(LVT_logunit,*) '[ERR] If concurrent use of these data sources are desired,'
       write(LVT_logunit,*) '[ERR] please generate the CDF for each source separately '
       write(LVT_logunit,*) '[ERR] (using LVT) and then supply them to LIS'
       call LVT_endrun()
    endif


       if ( SMOPSsmobs(i)%version == 'NESDIS V3.0 REGEN' ) then
          yr1 = 2012; mo1 = 8; da1 = 1; hr1 = 0; mn1 = 0; ss1 = 0
          call LVT_date2time(SMOPSsmobs(i)%version3_time,updoy,upgmt,&
             yr1,mo1,da1,hr1,mn1,ss1)
       else
          ! Actual version dates
          yr1 = 2016; mo1 = 10; da1 = 31; hr1 = 12; mn1 = 0; ss1 = 0
          call LVT_date2time(SMOPSsmobs(i)%version2_time,updoy,upgmt,&
             yr1,mo1,da1,hr1,mn1,ss1)

          yr1 = 2017; mo1 = 8; da1 = 24; hr1 = 12; mn1 = 0; ss1 = 0
          call LVT_date2time(SMOPSsmobs(i)%version3_time,updoy,upgmt,&
             yr1,mo1,da1,hr1,mn1,ss1)
       endif


    SMOPSsmobs(i)%startmode = .true. 

    SMOPSsmobs(i)%smopsnc = 1440 !   rtsmopsnc = 1440
    SMOPSsmobs(i)%smopsnr = 720  ! rtsmopsnr = 720
    
    call map_set(PROJ_LATLON, -89.875,-179.875,&
         0.0, 0.25,0.25, 0.0,&
         SMOPSsmobs(i)%smopsnc,SMOPSsmobs(i)%smopsnr,&
         SMOPSsmobs(i)%smopsproj)
    
    gridDesci = 0
    gridDesci(1) = 0
    gridDesci(2) = 1440
    gridDesci(3) = 720
    gridDesci(4) = -89.875
    gridDesci(5) = -179.875
    gridDesci(6) = 128
    gridDesci(7) = 89.875
    gridDesci(8) = 179.875
    gridDesci(9) = 0.25
    gridDesci(10) = 0.25
    gridDesci(20) = 64

    allocate(SMOPSsmobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(SMOPSsmobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(SMOPSsmobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(SMOPSsmobs(i)%maxv(LVT_rc%lnc,LVT_rc%lnr))
    allocate(SMOPSsmobs(i)%minv(LVT_rc%lnc,LVT_rc%lnr))

    SMOPSsmobs(i)%maxv  = -1000000.0
    SMOPSsmobs(i)%minv  = 1000000.0
    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,SMOPSsmobs(i)%rlat, &
         SMOPSsmobs(i)%rlon,SMOPSsmobs(i)%n11)

  end subroutine SMOPSsm_obsinit


end module SMOPSsm_obsMod
