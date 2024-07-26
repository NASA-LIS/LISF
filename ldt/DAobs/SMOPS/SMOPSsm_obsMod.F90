!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: SMOPSsm_obsMod
! 
! !DESCRIPTION: 
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
! !REVISION HISTORY: 
!  04 May 2013: Sujay Kumar, Initial Specification
!  05 July 2017 : Mahdi Navari, Modified to read AMSR2, SMAP
!      Dec 2017 : Mahdi Navari, Condition added to remove pixcls close to openwater
!
module SMOPSsm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOPSsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOPSsmobs
!EOP
  type, public :: smopssmobsdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     integer                :: useASCAT
     integer                :: useSMOS
     integer                :: useAMSR2
     integer                :: useSMAP
     real                   :: search_radius
     integer                :: mo
     integer                :: smopsnc, smopsnr
     real*8                 :: version2_time, version3_time
     real,    allocatable   :: smobs(:,:)
     integer, allocatable   :: n11(:)
     character(len=17)      :: version
     type(proj_info)        :: smopsproj
  end type smopssmobsdec

  type(smopssmobsdec), allocatable:: SMOPSsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMOPSsm_obsInit
! \label{SMOPSsm_obsInit}
! 
! !INTERFACE: 
  subroutine SMOPSsm_obsinit()
! !USES: 
   use LDT_coreMod,      only : LDT_rc, LDT_config
   use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
   use LDT_timeMgrMod,   only : LDT_clock, LDT_calendar, LDT_date2time
   use LDT_logMod

  ! MN : added to read the greenness and land cover
   use LDT_gfracMod,     only : LDT_gfrac_struc
   use LDT_LMLCMod

    implicit none
! !ARGUMENTS: 
!
! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading RT SMOPS soil moisture data. 
! 
!EOP
    integer                 :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 
    integer                 :: updoy,yr1,mo1,da1,hr1,mn1,ss1
    real                    :: upgmt
!______________________________________________________

    allocate(SMOPSsmobs(LDT_rc%nnest))

    write(LDT_logunit,*) '[INFO] NOAA SMOPS soil moisture field selected ...'

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOPS soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOPSsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'SMOPS soil moisture observation directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOPS soil moisture use ASCAT data:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOPSsmobs(n)%useASCAT, &
            rc=status)
       call LDT_verify(status, &
            'SMOPS soil moisture use ASCAT data: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOPS soil moisture use SMOS data:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOPSsmobs(n)%useSMOS, &
            rc=status)
       call LDT_verify(status, &
            'SMOPS soil moisture use SMOS data: not defined')
    enddo
 
    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOPS soil moisture use AMSR2 data:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOPSsmobs(n)%useAMSR2, &
            rc=status)
       call LDT_verify(status, &
            'SMOPS soil moisture use AMSR2 data: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOPS soil moisture use SMAP data:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOPSsmobs(n)%useSMAP, &
            rc=status)
       call LDT_verify(status, &
            'SMOPS soil moisture use SMAP data: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, "SMOPS soil moisture version:", rc=status)
    do n=1, LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,SMOPSsmobs(n)%version,&
          default='date-based', rc=status)
    enddo

! MN: get attribute of the search radious    
    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOPS search radius for openwater proximity detection:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOPSsmobs(n)%search_radius, &
            rc=status)
       call LDT_verify(status, &
            'SMOPS search radius for openwater proximity detection: not defined')
    enddo

 
    do n=1,LDT_rc%nnest

       if(SMOPSsmobs(n)%useASCAT + SMOPSsmobs(n)%useAMSR2 + & 
          SMOPSsmobs(n)%useSMAP + SMOPSsmobs(n)%useSMOS.gt.1) then
          write(LDT_logunit,*) '[ERR] Please do not select multiple sensor sources'
          write(LDT_logunit,*) '[ERR]  simultaneously for LDT preprocessing ...'
          write(LDT_logunit,*) '[ERR] If concurrent use of these data sources are desired,'
          write(LDT_logunit,*) '[ERR]  please generate the CDF for each source separately '
          write(LDT_logunit,*) '[ERR]  (using LDT) and then supply them to LIS.'
          call LDT_endrun()          
       endif

       if ( SMOPSsmobs(n)%version == 'NESDIS V3.0 REGEN' ) then
          yr1 = 2012; mo1 = 8; da1 = 1; hr1 = 0; mn1 = 0; ss1 = 0
          call LDT_date2time(SMOPSsmobs(n)%version3_time,updoy,upgmt,&
             yr1,mo1,da1,hr1,mn1,ss1)
       else
          ! Actual version dates
          yr1 = 2016; mo1 = 10; da1 = 31; hr1 = 12; mn1 = 0; ss1 = 0
          call LDT_date2time(SMOPSsmobs(n)%version2_time,updoy,upgmt,&
             yr1,mo1,da1,hr1,mn1,ss1)

          yr1 = 2017; mo1 = 8; da1 = 24; hr1 = 12; mn1 = 0; ss1 = 0
          call LDT_date2time(SMOPSsmobs(n)%version3_time,updoy,upgmt,&
             yr1,mo1,da1,hr1,mn1,ss1)
       endif

       allocate(SMOPSsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       SMOPSsmobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
       SMOPSsmobs(n)%smopsnc = 1440
       SMOPSsmobs(n)%smopsnr = 720

       call map_set(PROJ_LATLON, -89.875,-179.875,&
            0.0, 0.25,0.25, 0.0,&
            SMOPSsmobs(n)%smopsnc,SMOPSsmobs(n)%smopsnr,&
            SMOPSsmobs(n)%smopsproj)
       
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
       
       allocate(SMOPSsmobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call neighbor_interp_input(n, gridDesci, SMOPSsmobs(n)%n11)

       ! For ASCAT field, landcover and greenness fraction
       !  parameters can be used to screen for desert / city
       !  locations.  Check here if those parameters were written
       !  read in:
       if(LDT_gfrac_struc(n)%gfrac%selectOpt.gt.0) then
          write(LDT_logunit,*) '[INFO] Greenness fraction and landcover parameters '
          write(LDT_logunit,*) '[INFO]  have been selected, which can screen desert'
          write(LDT_logunit,*) '[INFO]  and urban based areas in ASCAT. '
       else
          if( SMOPSsmobs(n)%useASCAT > 0 ) then
             write(LDT_logunit,*) '[INFO] SMOPS-ASCAT field selected ... '
             write(LDT_logunit,*) '[WARN] ** Enabling greenness fraction '// &
                                  'will improve the SMOPS ASCAT filtering.'
          endif
       endif

    enddo
  end subroutine SMOPSsm_obsinit
     
end module SMOPSsm_obsMod
