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
! !MODULE: IMDPRCP_obsMod
! \label(IMDPRCP_obsMod)
!
! !INTERFACE:
module IMDPRCP_obsMod
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
!  This module handles the observation plugin for the 
!  Indian Meteorological Department (IMD) 
!  unified gauge-based analysis of daily precipitation. 
!  This data is part of products suite from the IMD unified
!  precipitation project. 
!  
!  Coverage is from 1 Jan 1960 - 2013, daily at 
!  0.25 deg lat/lon resolution
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  09 Aug 2017   Sujay Kumar  Initial Specification
! 
!EOP

  PUBLIC :: IMDPRCP_obsinit
  PUBLIC :: IMDPRCPobs

  type, public :: imdprcpobsdec
     character*100        :: odir
     integer              :: nc, nr
     integer              :: yr
     type(ESMF_Time)         :: startTime
     type(ESMF_TimeInterval) :: timeStep
     
     real,    allocatable     :: rlat(:)
     real,    allocatable     :: rlon(:)
     integer, allocatable     :: n11(:)
     integer, allocatable     :: n12(:)
     integer, allocatable     :: n21(:)
     integer, allocatable     :: n22(:)
     real,    allocatable     :: w11(:)
     real,    allocatable     :: w12(:)
     real,    allocatable     :: w21(:)
     real,    allocatable     :: w22(:)
     real,    allocatable     :: rainf(:,:,:)
  end type imdprcpobsdec

  type(imdprcpobsdec), allocatable :: imdprcpobs(:)

contains
  
!BOP
! 
! !ROUTINE: IMDPRCP_obsInit
! \label{IMDPRCP_obsInit}
!
! !INTERFACE: 
  subroutine IMDPRCP_obsinit(i)
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
!  for reading IMD PCP data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    real               :: gridDesci(50)
    integer            :: status
    
    if(.not.allocated(imdprcpobs)) then 
       allocate(imdprcpobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, imdprcpobs(i)%odir, &
         label='IMD PCP observation directory:',rc=status)
    call LVT_verify(status, 'IMD PCP observation directory: not defined')

    gridDesci = 0 
    call LVT_update_timestep(LVT_rc, 86400)

    allocate(imdprcpobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(imdprcpobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(imdprcpobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(imdprcpobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(imdprcpobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(imdprcpobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(imdprcpobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(imdprcpobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(imdprcpobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(imdprcpobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    imdprcpobs(i)%nc = 135
    imdprcpobs(i)%nr = 129
    
    gridDesci(1) = 0 
    gridDesci(2) = 135
    gridDesci(3) = 129 
    gridDesci(4) = 6.50
    gridDesci(5) = 66.50
    gridDesci(6) = 128 
    gridDesci(7) = 38.50
    gridDesci(8) = 100.00
    gridDesci(9) = 0.25
    gridDesci(10) = 0.25
    gridDesci(20) = 64
    
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         imdprcpobs(i)%rlat, imdprcpobs(i)%rlon, &  
         imdprcpobs(i)%n11, imdprcpobs(i)%n12,   & 
         imdprcpobs(i)%n21, imdprcpobs(i)%n22,   & 
         imdprcpobs(i)%w11, imdprcpobs(i)%w12,   & 
         imdprcpobs(i)%w21, imdprcpobs(i)%w22)

    allocate(imdprcpobs(i)%rainf(imdprcpobs(i)%nc,&
         imdprcpobs(i)%nr,&
         366))
    imdprcpobs(i)%yr = -1

    call ESMF_TimeSet(imdprcpobs(i)%startTime,  yy=LVT_rc%yr, &
         mm = 1, &
         dd = 1, &
         h = 0, &
         m = 0, &
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status, 'error in setting IMD start time')

!    call ESMF_TimeSet(imdprcpobs(i)%stopTime, yy=eyr, &
!         mm = emo, &
!         dd = eda, &
!         h = ehr, &
!         m = emn, &
!         calendar = LVT_calendar, &
!         rc=status)
!    call LVT_verify(status, 'error in setting scan stop time')

    call ESMF_TimeIntervalSet(imdprcpobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (imdprcpobs)')
        
  end subroutine IMDPRCP_obsinit


end module IMDPRCP_obsMod
