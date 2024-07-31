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
! !MODULE: SNODAS_obsMod
! \label(SNODAS_obsMod)
!
! !INTERFACE:
module SNODAS_obsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the observation plugin for the NOAA National Weather 
!  Service's National Operational Hydrologic Remote Sensing Center (NOHRSC)
!  SNOw Data Assimilation (SNODAS) data. SNODAS is a modeling and data
!  assimilation system developed by NOHRSC to provide the best possible 
!  estimates of snow cover and associated parameters to support hydrologic
!  modeling and analysis. The aim of SNODAS is to provide a physically
!  consistent framework to integrate snow data from satellite, airborne 
!  platforms, and ground stations with model estimates of snow cover 
!  (Carroll et al. 2001). SNODAS includes procedures to ingest and downscale 
!  output from the Numerical Weather Prediction (NWP) models, and to simulate
!  snowcover using a physically based, spatially-distributed energy- and 
!  mass-balance snow model. SNODAS also includes procedures to assimilate
!  satellite-derived, airborne, and ground-based observations of snow covered 
!  area and Snow Water Equivalent (SWE).
! 
!   Website: http://nsidc.org/data/g02158.html
! 
! Citation: National Operational Hydrologic Remote Sensing Center. 2004. 
! Snow Data Assimilation System (SNODAS) Data Products at NSIDC. Boulder, 
! Colorado USA: National Snow and Ice Data Center. Digital media.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  6 May 2010   Sujay Kumar  Initial Specification
! 
!EOP
! !USES: 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SNODAS_obsinit !Initializes structures for reading SNODAS data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SNODASobs !Object to hold SNODAS observation attributes
!EOP
  type, public :: snodasobsdec
     character*100           :: odir
     integer                 :: nc
     integer                 :: nr
     real                    :: udef
     real,        allocatable    :: stnlat(:)
     real,        allocatable    :: stnlon(:)
     type(ESMF_Clock)        :: clock
     type(ESMF_Time)         :: startTime, startMonth
     type(ESMF_TimeInterval) :: timestep
     real,  allocatable          :: swe(:,:)
     real,  allocatable          :: snwd(:,:)
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)
     integer, allocatable        :: n11(:)
     integer, allocatable        :: n12(:)
     integer, allocatable        :: n21(:)
     integer, allocatable        :: n22(:)
     real, allocatable           :: w11(:), w12(:)
     real, allocatable           :: w21(:), w22(:)
  end type snodasobsdec

  type(snodasobsdec), allocatable :: snodasobs(:)

contains
  
!BOP
! 
! !ROUTINE: SNODAS_obsInit
! \label{SNODAS_obsInit}
!
! !INTERFACE:
  subroutine SNODAS_obsinit(i)
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
!  for reading SNODAS data. The SNODAS data is provides in the WGS 1984 grid
!  with 30 arc second resolution. The domain extents are from (24.9504, -124.7337) 
!  to (52.8754, -66.9421). The data entries are 16-bit signed integers. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer            :: status, rc
    integer            :: ftn, k
    real*8             :: tdur
    integer            :: syr, smo, sda, shr, smn, sss
    integer            :: eyr, emo, eda, ehr, emn, ess
    integer            :: ts
    character*100      :: coordfile
    character*100      :: mdata
    real               :: xi1,xj1,xmesh,orient,alat1,alon1
    integer            :: t
    real               :: gridDesci(50)

    if(.not.allocated(snodasobs)) then 
       allocate(snodasobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, snodasobs(i)%odir, &
         label='SNODAS observation directory:',rc=status)
    call LVT_verify(status, 'SNODAS observation directory: not defined')

    snodasobs(i)%udef = -9999.0
    syr =2007
    smo = 11
    sda = 1
    shr = 0
    smn = 0 
    sss = 0 

    ts = 86400
    call LVT_update_timestep(LVT_rc, 86400)

    call ESMF_TimeSet(snodasobs(i)%startTime,yy=syr, &
         mm=smo, dd=sda, h=shr, m=smn, &
         s = 0, calendar=LVT_calendar, rc=status)
    call LVT_verify(status, 'snodas: starttime set failed')


    call ESMF_TimeIntervalSet(snodasobs(i)%timestep, s=ts, rc=status)
    call LVT_verify(status, 'error in setting timestep (snodasobs(i))')

    snodasobs(i)%nc = 6935
    snodasobs(i)%nr = 3351

    gridDesci = 0 
    gridDesci(1) = 0
    gridDesci(2) = snodasobs(i)%nc
    gridDesci(3) = snodasobs(i)%nr
    griddesci(4) = 24.94958
    gridDesci(5) = -124.73375
    gridDesci(6) = 128
    gridDesci(7) = 52.8704166
    gridDesci(8) = -66.9420833
    gridDesci(9) = 0.00833
    gridDesci(10) =0.00833 
    gridDesci(20) = 64

    allocate(snodasobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(snodasobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    allocate(snodasobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(snodasobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(snodasobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(snodasobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    allocate(snodasobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(snodasobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(snodasobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(snodasobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr,snodasobs(i)%rlat,snodasobs(i)%rlon,&
         snodasobs(i)%n11,snodasobs(i)%n12,snodasobs(i)%n21,&
         snodasobs(i)%n22,snodasobs(i)%w11,snodasobs(i)%w12,&
         snodasobs(i)%w21,snodasobs(i)%w22)

    allocate(snodasobs(i)%swe(LVT_rc%lnc,LVT_rc%lnr))
    snodasobs(i)%swe = LVT_rc%udef
    allocate(snodasobs(i)%snwd(LVT_rc%lnc,LVT_rc%lnr))
    snodasobs(i)%snwd = LVT_rc%udef

  end subroutine SNODAS_obsinit


end module SNODAS_obsMod
