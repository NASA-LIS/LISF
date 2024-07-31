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
! !MODULE: CMCSNWD_obsMod
!  \label(CMCSNWD_obsMod)
!
! !INTERFACE:
module CMCSNWD_obsMod
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
!  This module handles the observation plugin for the Canadian Meteorological 
!  Center (CMC) Daily Snow Depth Analysis data. This data set consists of Northern 
! Hemisphere snow depth analysis data processed by the Canadian Meteorological 
! Centre (CMC). Snow depth data obtained from surface synoptic observations 
! (synops), meteorological aviation reports (metars), and special aviation 
! reports (SAs) were acquired from the World Meteorological Organization
!  (WMO) information system for use in the CMC analyses. This CMC data set 
! includes daily observations from 1998 through 2009 and will be updated annually. 
! Monthly averages and monthly climatologies of snow depth and estimated Snow Water
! Equivalent (SWE) are provided, where SWE was estimated using a density look-up table. 
! Note that LVT does not handle these climatologies. 
! 
! Website: http://nsidc.org/data/nsidc-0447.html
! 
! Citation: Brown, Ross D. and Bruce Brasnett. 2010. Canadian Meteorological 
! Centre (CMC) Daily Snow Depth Analysis Data. © Environment Canada, 2010.
! Boulder, Colorado USA: National Snow and Ice Data Center. Digital media.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  23 Apr 2010   Sujay Kumar  Initial Specification
! 
!EOP
!
! 
!
! 
!

  PRIVATE 

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: CMCSNWD_obsinit !Initializes structures for reading CMCSNWD data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: CMCSNWDobs !Object to hold CMCSNWD observation attributes
!EOP
  type, public :: cmcsnwdobsdec
     character*100               :: odir
     integer                     :: nc
     integer                     :: nr
     real                        :: udef
     real,        allocatable    :: stnlat(:)
     real,        allocatable    :: stnlon(:)
     type(ESMF_Clock)            :: clock
     type(ESMF_Time)             :: startTime, startMonth
     type(ESMF_TimeInterval)     :: timestep
     real,  allocatable          :: snwd(:,:)
     integer, allocatable        :: n11(:)
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)
  end type cmcsnwdobsdec

  type(cmcsnwdobsdec), allocatable :: cmcsnwdobs(:)

contains
  
!BOP
! 
! !ROUTINE: CMCSNWD_obsInit
! \label{CMCSNWD_obsInit}
!
! !INTERFACE: 
  subroutine CMCSNWD_obsinit(i)
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
!  for reading CMCSNWD data. 
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
    integer            :: ts
    character*100      :: coordfile
    character*100      :: mdata
    real               :: xi1,xj1,xmesh,orient,alat1,alon1
    integer            :: t
    real               :: gridDesci(50)

    if(.not.allocated(cmcsnwdobs)) then 
       allocate(cmcsnwdobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, cmcsnwdobs(i)%odir, &
         label='CMC snow depth observation directory:',rc=status)
    !call LVT_verify(status, 'CMCSNWD observation directory: not defined')
    ! Modified by Shugong Wang 11/20/2011
    call LVT_verify(status, 'CMC snow depth observation directory: not defined')

    cmcsnwdobs(i)%udef = -1
    syr = 2002
    smo = 1
    sda = 1
    shr = 0
    smn = 0 
    call ESMF_TimeSet(cmcsnwdobs(i)%startTime,yy=syr, &
         mm=smo, dd=sda, h=shr, m=smn, &
         s = 0, calendar=LVT_calendar, rc=status)
    call LVT_verify(status, 'cmcsnwd: starttime set failed')
    
    ts = 86400
    call ESMF_TimeIntervalSet(cmcsnwdobs(i)%timestep, s=ts, rc=status)
    call LVT_verify(status, 'error in setting timestep (cmcsnwdobs)')

    call LVT_update_timestep(LVT_rc, 86400)

    cmcsnwdobs(i)%nc = 706
    cmcsnwdobs(i)%nr = 706

    xi1 = 1.0 - 353
    xj1 = 1.0 - 353
    xmesh = 23.8125
    orient = 190.0

    call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)
    
    gridDesci = 0 
    gridDesci(1) = 5
    gridDesci(2) = cmcsnwdobs(i)%nc
    gridDesci(3) = cmcsnwdobs(i)%nr
    griddesci(4) = alat1
    gridDesci(5) = alon1
    gridDesci(6) = 8
    gridDesci(7) = orient
    gridDesci(8) = xmesh
    gridDesci(9) = xmesh
    gridDesci(10) = 0.0
    gridDesci(13) = 1
    gridDesci(20) = 0 

    allocate(cmcsnwdobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(cmcsnwdobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    allocate(cmcsnwdobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))

    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,cmcsnwdobs(i)%rlat,cmcsnwdobs(i)%rlon,&
         cmcsnwdobs(i)%n11)

!    call LVT_initializeLSMdataEntry(LVT_MOC_SNOWDEPTH,i,&
!         LVT_getdataEntryUnits(LVT_MOC_SNOWDEPTH),"-",1,1)
    allocate(cmcsnwdobs(i)%snwd(LVT_rc%lnc,LVT_rc%lnr))
    cmcsnwdobs(i)%snwd = LVT_rc%udef

  end subroutine CMCSNWD_obsinit


end module CMCSNWD_obsMod
