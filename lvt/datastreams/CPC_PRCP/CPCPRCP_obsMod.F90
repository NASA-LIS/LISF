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
! !MODULE: CPCPRCP_obsMod
! \label(CPCPRCP_obsMod)
!
! !INTERFACE:
module CPCPRCP_obsMod
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
!  This module handles the observation plugin for the CPC 
!  unified gauge-based analysis of daily precipitation. 
!  This data is part of products suite from the CPC unified
!  precipitation project. 
!  
!  Coverage is from 1 Jan 1979 - present, daily at 
!  0.25 deg lat/lon resolution
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  11 May 2011   Sujay Kumar  Initial Specification
! 
!EOP

  PUBLIC :: CPCPRCP_obsinit
  PUBLIC :: CPCPRCPobs

  type, public :: cpcprcpobsdec
     character*100        :: odir
     integer              :: userealtime
     character*30         :: domainType
     integer              :: nc, nr
     
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
  end type cpcprcpobsdec

  type(cpcprcpobsdec), allocatable :: cpcprcpobs(:)

contains
  
!BOP
! 
! !ROUTINE: CPCPRCP_obsInit
! \label{CPCPRCP_obsInit}
!
! !INTERFACE: 
  subroutine CPCPRCP_obsinit(i)
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
!  for reading CPC PCP data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    real               :: gridDesci(50)
    integer            :: status
    
    if(.not.allocated(cpcprcpobs)) then 
       allocate(cpcprcpobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, cpcprcpobs(i)%odir, &
         label='CPC PCP observation directory:',rc=status)
    call LVT_verify(status, 'CPC PCP observation directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, cpcprcpobs(i)%domainType, &
         label='CPC PCP domain type (CONUS or GLOBAL):',rc=status)
    call LVT_verify(status, 'CPC PCP domain type (CONUS or GLOBAL): not defined')

    call ESMF_ConfigGetAttribute(LVT_config, cpcprcpobs(i)%userealtime, &
         label='CPC PCP use real time data:',rc=status)
    call LVT_verify(status, 'CPC PCP use real time data: not defined')

    gridDesci = 0 
    call LVT_update_timestep(LVT_rc, 86400)

    allocate(cpcprcpobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(cpcprcpobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(cpcprcpobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(cpcprcpobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(cpcprcpobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(cpcprcpobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(cpcprcpobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(cpcprcpobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(cpcprcpobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(cpcprcpobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    if(cpcprcpobs(i)%domainType.eq."CONUS") then 
       cpcprcpobs(i)%nc = 300
       cpcprcpobs(i)%nr = 120
       
       gridDesci(1) = 0 
       gridDesci(2) = 300
       gridDesci(3) = 120 
       gridDesci(4) = 20.125
       gridDesci(5) = -129.875
       gridDesci(6) = 128 
       gridDesci(7) = 49.875
       gridDesci(8) = -55.125
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64

    elseif(cpcprcpobs(i)%domainType.eq."GLOBAL") then 
       cpcprcpobs(i)%nc = 720
       cpcprcpobs(i)%nr = 360

       gridDesci(1) = 0 
       gridDesci(2) = 720
       gridDesci(3) = 360 
       gridDesci(4) = -89.75
       gridDesci(5) = -179.75
       gridDesci(6) = 128 
       gridDesci(7) = 89.75
       gridDesci(8) = 179.75
       gridDesci(9) = 0.5
       gridDesci(10) = 0.5
       gridDesci(20) = 64

    endif

    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         cpcprcpobs(i)%rlat, cpcprcpobs(i)%rlon, &  
         cpcprcpobs(i)%n11, cpcprcpobs(i)%n12,   & 
         cpcprcpobs(i)%n21, cpcprcpobs(i)%n22,   & 
         cpcprcpobs(i)%w11, cpcprcpobs(i)%w12,   & 
         cpcprcpobs(i)%w21, cpcprcpobs(i)%w22)

  end subroutine CPCPRCP_obsinit


end module CPCPRCP_obsMod
