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
! !MODULE: ARM_obsMod
!  \label(ARM_obsMod)
!
! !INTERFACE:
module ARM_obsMod
! 
! !USES: 
  use ESMF

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the Atmospheric 
!  Radiation Measurement (ARM) in-situ data. The datasets handled
!  in this plugin include measurements of surface properties 
!  (latent, sensible and ground heat fluxes, soil moisture and 
!  soil temperature profiles). 
!  
!  web link: http://www.archive.arm.gov/armlogin/login.jsp
!
!  The flux datasets contain bulk 
!  aerodynamic estimates of fluxes at 30 min intervals from 
!  EBBR instrument (energy balance bowen ratio) and/or 
!  from the ECOR instrument (eddy correlation)  
!  (ref: http://www.arm.gov/instruments/ebbr, 
!        http://www.arm.gov/instruments/ecor)
!
!  The soil moisture and temperature profiles data is from the
!  SWATS instruments (soil water and temperature system). 
!  (ref: http://www.arm.gov/instruments/swats)
!
!  TODO: 
!  Rainfall measurements are from SWATSPCP instrument system 
!  (soil water and temperature profiling system rain gauges)
!
!  The data file naming convention is: 
!
!  (sss)(nn)(inst)(qqq)(Fn).(ln).YYYYMMDD.hhmmss.cdf  
!
!  where:
!  sss    is the site identifier (e.g. sgp, twp, nsa)
!  nn     is optionally the data integration period in minutes (e.g. 1, 5,1440)
! inst    is the instrument basename (e.g. mwr, wsi, mpl)
! qqq     is an optional qualifier that distinguishes these data from other data sets produced by the same instrument
! Fn      is the facility designation (e.g. C1, E13, B4)
! ln      is the data level (e.g. a0, a1, b1, c1)
!  
!  The data level flags represent: 
!  
!  00    - raw data
!  a0    - data converted to netcdf
!  b1    - basic quality control checks applied
!  b2-b9 - b1 data with additional QC's
!  c1    - derived data products
!  c2-c9 - c2 data with additional processing.  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  18 Nov 2010   Sujay Kumar  Initial Specification
! 
!EOP
! 
! 
!
  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ARM_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ARMobs
!EOP
  type, public :: armobsdec
     character*100           :: odir
     integer                 :: n_stns 
     character*30            :: site_id
     integer                 :: baebbr_select
     integer                 :: ebbr_select
     integer                 :: ecor_select
     integer                 :: swats_select
     integer                 :: smos_select

     real,           allocatable :: stnlat(:)
     real,           allocatable :: stnlon(:)
     character*100,  allocatable :: stn_name(:)
     real,           allocatable :: stnbd(:)

     logical                 :: startFlag
     integer                 :: nts
     real                    :: udef
     type(ESMF_Time)         :: starttime
     type(ESMF_TimeInterval) :: baebbr_ts
     type(ESMF_TimeInterval) :: ebbr_ts
     type(ESMF_TimeInterval) :: ecor_ts
     type(ESMF_TimeInterval) :: swats_ts
     type(ESMF_TimeInterval) :: smos_ts

     integer                 :: da

     real, allocatable           :: baebbr_qle(:,:)
     real, allocatable           :: baebbr_qh(:,:)
     real, allocatable           :: baebbr_qg(:,:)
     real, allocatable           :: baebbr_netrad(:,:)
     real, allocatable           :: baebbr_solarrad(:,:)
     integer, allocatable        :: baebbr_tindex(:,:)

     real, allocatable           :: ebbr_sfsm(:,:)
     real, allocatable           :: ebbr_sfst(:,:)
     integer, allocatable        :: ebbr_tindex(:,:)

     real, allocatable           :: ecor_qle(:,:)
     real, allocatable           :: ecor_qh(:,:)
     integer, allocatable        :: ecor_tindex(:,:)

     real, allocatable           :: swats_sfst(:,:)  
     real, allocatable           :: swats_sfsm(:,:) 
     real, allocatable           :: swats_rzst(:,:)  
     real, allocatable           :: swats_rzsm(:,:) 
     integer, allocatable        :: swats_tindex(:,:)

     real, allocatable           :: smos_pcp(:,:)
     real, allocatable           :: smos_wspd(:,:)
     real, allocatable           :: smos_prs(:,:)
     real, allocatable           :: smos_sh(:,:)
     real, allocatable           :: smos_temp(:,:)
     real, allocatable           :: smos_snow(:,:)
     integer, allocatable        :: smos_tindex(:,:)

  end type armobsdec

  type(armobsdec), allocatable:: armobs(:)

contains
  
!BOP
! 
! !ROUTINE: ARM_obsInit
! \label{ARM_obsInit}
!
! !INTERFACE: 
  subroutine ARM_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading ARM data. 
!
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
!! !ARGUMENTS: 
    integer,   intent(IN) :: i 
!!EOP
    integer                 :: status
    integer                 :: ftn
    character*100           :: stnlist_file
    character*100           :: currentLine
    integer                 :: k, iloc


    if(.not.allocated(armobs)) then 
       allocate(armobs(LVT_rc%nDataStreams))
    endif
!
    armobs(i)%startFlag = .true. 

    call ESMF_ConfigGetAttribute(LVT_Config, armobs(i)%odir, &
         label='ARM observation directory:', rc=status)
    call LVT_verify(status, 'ARM observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, armobs(i)%site_id, & 
         label='ARM site identifier name:',rc=status)
    call LVT_verify(status, 'ARM site identifier name: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, stnlist_file, &
         label='ARM station list file:', rc=status)
    call LVT_verify(status, 'ARM station list file: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, armobs(i)%baebbr_select, &
         label='ARM use BAEBBR data:', rc=status)
    call LVT_verify(status, 'ARM use BAEBBR data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, armobs(i)%ebbr_select, &
         label='ARM use EBBR data:', rc=status)
    call LVT_verify(status, 'ARM use EBBR data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, armobs(i)%ecor_select, &
         label='ARM use ECOR flux data:', rc=status)
    call LVT_verify(status, 'ARM use ECOR flux data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, armobs(i)%swats_select, &
         label='ARM use SWATS data:', rc=status)
    call LVT_verify(status, 'ARM use SWATS data: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, armobs(i)%smos_select, &
         label='ARM use SMOS data:', rc=status)
    call LVT_verify(status, 'ARM use SMOS data: not defined')

    write(LVT_logunit,*) '[INFO] Processing ARM data locations '

    ftn = LVT_getNextUnitNumber()
    open(ftn, file=trim(stnlist_file), form='formatted')
    read(ftn,*)
    read(ftn,*) armobs(i)%n_stns
    read(ftn,*) 

    allocate(armobs(i)%stn_name(armobs(i)%n_stns))
    allocate(armobs(i)%stnlat(armobs(i)%n_stns))
    allocate(armobs(i)%stnlon(armobs(i)%n_stns))
    allocate(armobs(i)%stnbd(armobs(i)%n_stns))

    do k=1,armobs(i)%n_stns
       read(ftn,'(a)') currentLine
       iloc = Index(currentLine, ";")
       read(currentLine(1:iloc -1), *) armobs(i)%stn_name(k)
       currentLine = currentLine(iloc+1:Len(currentLine))
       
       iloc = Index(currentLine, ";")
       read(currentLine(1:iloc -1), *) armobs(i)%stnlat(k)
       currentLine = currentLine(iloc+1:Len(currentLine))

       iloc = Index(currentLine, ";")
       read(currentLine(1:iloc -1), *) armobs(i)%stnlon(k)
       currentLine = currentLine(iloc+1:Len(currentLine))

       read(currentLine, *) armobs(i)%stnbd(k)
       
       write(LVT_logunit,*) '[INFO] ',armobs(i)%stn_name(k), &
            armobs(i)%stnlat(k), armobs(i)%stnlon(k),armobs(i)%stnbd(k)
    enddo
    
    call LVT_releaseUnitNumber(ftn)

    armobs(i)%udef   = 9999.9


    call ESMF_TimeIntervalSet(armobs(i)%baebbr_ts, s=1800,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: ARM_obsMod ')

    call ESMF_TimeIntervalSet(armobs(i)%ebbr_ts, s=1800,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: ARM_obsMod ')

    call ESMF_TimeIntervalSet(armobs(i)%ecor_ts, s=1800,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: ARM_obsMod ')

    call ESMF_TimeIntervalSet(armobs(i)%swats_ts, s=3600,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: ARM_obsMod ')

    call ESMF_TimeIntervalSet(armobs(i)%smos_ts, s=1800,rc=status)
    call LVT_verify(status, 'Error in timeintervalset: ARM_obsMod ')
    
    if(armobs(i)%baebbr_select.eq.1) then 
       allocate(armobs(i)%baebbr_qle(armobs(i)%n_stns, 48))
       allocate(armobs(i)%baebbr_qh(armobs(i)%n_stns, 48))
       allocate(armobs(i)%baebbr_qg(armobs(i)%n_stns, 48))        
       allocate(armobs(i)%baebbr_netrad(armobs(i)%n_stns, 48))        
       allocate(armobs(i)%baebbr_solarrad(armobs(i)%n_stns, 48))        
       allocate(armobs(i)%baebbr_tindex(armobs(i)%n_stns, 48))        

       armobs(i)%baebbr_qle = LVT_rc%udef
       armobs(i)%baebbr_qh  = LVT_rc%udef
       armobs(i)%baebbr_qg  = LVT_rc%udef
       armobs(i)%baebbr_netrad  = LVT_rc%udef
       armobs(i)%baebbr_solarrad  = LVT_rc%udef
    endif

    if(armobs(i)%ebbr_select.eq.1) then 
       allocate(armobs(i)%ebbr_sfsm(armobs(i)%n_stns, 48))        
       allocate(armobs(i)%ebbr_sfst(armobs(i)%n_stns, 48))        
       allocate(armobs(i)%ebbr_tindex(armobs(i)%n_stns,48))

       armobs(i)%ebbr_sfsm = LVT_rc%udef
       armobs(i)%ebbr_sfst = LVT_rc%udef
    endif

    if(armobs(i)%ecor_select.eq.1) then 
       allocate(armobs(i)%ecor_qle(armobs(i)%n_stns, 48))
       allocate(armobs(i)%ecor_qh(armobs(i)%n_stns, 48))
       allocate(armobs(i)%ecor_tindex(armobs(i)%n_stns,48))

       armobs(i)%ecor_qle = LVT_rc%udef
       armobs(i)%ecor_qh  = LVT_rc%udef

    endif
       
    if(armobs(i)%swats_select.eq.1) then 
       allocate(armobs(i)%swats_sfsm(armobs(i)%n_stns,24))
       allocate(armobs(i)%swats_sfst(armobs(i)%n_stns,24))
       
       allocate(armobs(i)%swats_rzsm(armobs(i)%n_stns,24))
       allocate(armobs(i)%swats_rzst(armobs(i)%n_stns,24))
       allocate(armobs(i)%swats_tindex(armobs(i)%n_stns, 24))

       armobs(i)%swats_sfsm       = LVT_rc%udef
       armobs(i)%swats_sfst       = LVT_rc%udef
       armobs(i)%swats_rzsm       = LVT_rc%udef
       armobs(i)%swats_rzst       = LVT_rc%udef
    endif
    
    if(armobs(i)%smos_select.eq.1) then 
       allocate(armobs(i)%smos_pcp(armobs(i)%n_stns,48))
       allocate(armobs(i)%smos_prs(armobs(i)%n_stns,48))
       allocate(armobs(i)%smos_wspd(armobs(i)%n_stns,48))
       allocate(armobs(i)%smos_sh(armobs(i)%n_stns,48))
       allocate(armobs(i)%smos_temp(armobs(i)%n_stns,48))
       allocate(armobs(i)%smos_snow(armobs(i)%n_stns,48))
       allocate(armobs(i)%smos_tindex(armobs(i)%n_stns,48))

       armobs(i)%smos_pcp  = LVT_rc%udef
       armobs(i)%smos_prs  = LVT_rc%udef
       armobs(i)%smos_wspd  = LVT_rc%udef
       armobs(i)%smos_sh  = LVT_rc%udef
       armobs(i)%smos_temp  = LVT_rc%udef
       armobs(i)%smos_snow  = LVT_rc%udef
    endif

    armobs(i)%da = -1

    call LVT_update_timestep(LVT_rc, 1800)
  end subroutine ARM_obsinit


end module ARM_obsMod
