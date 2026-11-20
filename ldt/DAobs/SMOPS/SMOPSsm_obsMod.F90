!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
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
! There are four versions of the SMOPS datasets.
!
! According to the use by the 557th Weather Wing:
!
!                          version_1.3 <  2016-10-31T12:00:00
!   2016-10-31T12:00:00 <= version_2.0 <  2017-08-24T12:00:00
!   2017-08-24T12:00:00 <= version_3.0 <  2024-07-04T18:00:00
!                          version_4.0 >= 2024-07-04T18:00:00
!
! SMOPS version 4.0 drops ASCAT-A (Metop-A) and SMOS and adds ASCAT-C
! (Metop-C). For ASCAT, LDT reads ASCAT-A and ASCAT-B up to version 3.0
! and ASCAT-B and ASCAT-C from version 4.0, as LIS does.
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
! The 6-hourly files archived by the 557th Weather Wing ("AGRMET"
! naming convention) changed names twice:
!
!   smops_dYYYYMMDD_sHH0000_cness.gr2                 before 2017-10-05
!   NPR_SMOPS_CMAP_DYYYYMMDDHH.gr2                    2017-10-05 to
!                                                     2024-07-04T18
!   NPR-SMOPS-CMAP-6hrly_vNrM_blend_sYYYYMMDDHH*.grib2  from 2024-07-04T18
!
! The NPR_SMOPS_CMAP_DYYYYMMDDHH.gr2 files hold the 6 hours BEFORE HH;
! the other two hold the 6 hours starting at HH.
!
! The NPR-SMOPS-CMAP-6hrly names carry the SMOPS version: v4r0 from
! 2024-07-04T18, and v5r0 (AMSR3 added) from 2026-12-08. LDT takes the
! version of these files from the name. The GRIB2 parameter numbers of
! v5r0 are not known yet, so LDT stops on a v5r0 file; they are to be
! added in read_SMOPS_native once confirmed from a v5 file. Files whose
! names carry no version are dated: from 2026-12-08 (version5_time) LDT
! stops unless "SMOPS soil moisture version:" is set.
!
! With "SMOPS soil moisture use realtime data: 1", LDT reads every
! 6-hourly window that starts within the current LDT time step. With a
! 1-day time step ("LIS output timestep: 1da") the four windows of the day
! are combined into one daily field, the latest valid retrieval winning
! (as in the NESDIS daily product). With a 6-hour time step each window
! is logged separately, and "Temporal averaging interval: 1da" then
! averages them per day.
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
!  28 Oct 2026: Ehsan Jalilvand, Read the 6-hourly SMOPS archive (three
!               naming conventions) and SMOPS v4 (ASCAT-B and ASCAT-C),
!               following the LIS SMOPS ASCAT reader
!  12 Nov 2026: Ehsan Jalilvand, SMOPS v5: the version of the
!               NPR-SMOPS-CMAP-6hrly files is read from their names, and
!               LDT stops on versions whose GRIB2 parameter numbers it
!               does not have; "SMOPS soil moisture version:" is checked
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
     integer                :: useRealtime
     character(len=10)      :: conv
     real                   :: search_radius
     integer                :: mo
     integer                :: smopsnc, smopsnr
     real*8                 :: version2_time, version3_time, version4_time
     real*8                 :: version5_time
     real*8                 :: naming3_time
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
       select case ( trim(SMOPSsmobs(n)%version) )
       case ( 'date-based', 'NESDIS V3.0 REGEN', '1.3', '2.0', '3.0', '4.0' )
       case default
          write(LDT_logunit,*) '[ERR] SMOPS soil moisture version: "', &
               trim(SMOPSsmobs(n)%version), '" is not supported.'
          write(LDT_logunit,*) '[ERR] Use "date-based", "NESDIS V3.0 REGEN", ', &
               '"1.3", "2.0", "3.0" or "4.0".'
          call LDT_endrun()
       end select
    enddo

    ! 6-hourly archive (as in LIS). Optional; default is the daily files.
    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOPS soil moisture use realtime data:', rc=status)
    if (status == 0) then
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config, SMOPSsmobs(n)%useRealtime, &
               rc=status)
          call LDT_verify(status, &
               'SMOPS soil moisture use realtime data: not defined')
       enddo
    else
       SMOPSsmobs(:)%useRealtime = 0
    endif

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOPS soil moisture naming convention:', rc=status)
    if (status == 0) then
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config, SMOPSsmobs(n)%conv, &
               rc=status)
          call LDT_verify(status, &
               'SMOPS soil moisture naming convention: not defined')
       enddo
    else
       SMOPSsmobs(:)%conv = 'AGRMET'
    endif

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
          ! version2_time was left unset here; data before the regen
          ! start are read as version 1.3.
          SMOPSsmobs(n)%version2_time = SMOPSsmobs(n)%version3_time
       else
          ! Actual version dates
          yr1 = 2016; mo1 = 10; da1 = 31; hr1 = 12; mn1 = 0; ss1 = 0
          call LDT_date2time(SMOPSsmobs(n)%version2_time,updoy,upgmt,&
             yr1,mo1,da1,hr1,mn1,ss1)

          yr1 = 2017; mo1 = 8; da1 = 24; hr1 = 12; mn1 = 0; ss1 = 0
          call LDT_date2time(SMOPSsmobs(n)%version3_time,updoy,upgmt,&
             yr1,mo1,da1,hr1,mn1,ss1)
       endif

       ! SMOPS v4 (ASCAT-B and ASCAT-C), same date as in LIS
       yr1 = 2024; mo1 = 7; da1 = 4; hr1 = 18; mn1 = 0; ss1 = 0
       call LDT_date2time(SMOPSsmobs(n)%version4_time,updoy,upgmt,&
          yr1,mo1,da1,hr1,mn1,ss1)

       ! SMOPS v5 (AMSR3 added) operational from 2026-12-08 (NESDIS notice
       ! to users). Used only for files whose names carry no version.
       yr1 = 2026; mo1 = 12; da1 = 8; hr1 = 0; mn1 = 0; ss1 = 0
       call LDT_date2time(SMOPSsmobs(n)%version5_time,updoy,upgmt,&
          yr1,mo1,da1,hr1,mn1,ss1)

       ! Switch from smops_d*_cness.gr2 to NPR_SMOPS_CMAP_D*.gr2 names
       yr1 = 2017; mo1 = 10; da1 = 5; hr1 = 0; mn1 = 0; ss1 = 0
       call LDT_date2time(SMOPSsmobs(n)%naming3_time,updoy,upgmt,&
          yr1,mo1,da1,hr1,mn1,ss1)

       if (SMOPSsmobs(n)%useRealtime == 1) then
          ! Each 6-hourly window must fall in exactly one time step.
          if (LDT_rc%ts < 21600.0 .or. mod(nint(LDT_rc%ts), 21600) /= 0) then
             write(LDT_logunit,*) '[ERR] SMOPS 6-hourly files need a time '// &
                  'step that is a multiple of 6 hours'
             write(LDT_logunit,*) '[ERR] (LIS output timestep: 6hr or 1da).'
             write(LDT_logunit,*) '[ERR] Current time step (s): ', LDT_rc%ts
             call LDT_endrun()
          endif
          write(LDT_logunit,*) '[INFO] Reading 6-hourly SMOPS files, '// &
               'naming convention: ', trim(SMOPSsmobs(n)%conv)
       else
          write(LDT_logunit,*) '[INFO] Reading daily SMOPS files '// &
               '(<dir>/YYYY/NPR_SMOPS_CMAP_DYYYYMMDD.gr2)'
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
