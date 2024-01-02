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
! !MODULE: ERAinterimLandobsMod
! \label(ERAinterimLandobsMod)
!
! !INTERFACE:
module ERAinterimLandobsMod
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
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ERAinterimLandobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ERAIlandobs !Object to hold ERAinterimLand observation attributes
!EOP

  type, public :: erainterimlanddec
     character*100           :: odir
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     integer                 :: da
     real                    :: datares
     logical                 :: startFlag
     type(ESMF_Time)         :: starttime

  end type erainterimlanddec
     
  type(erainterimlanddec), allocatable :: ERAIlandObs(:)

contains
  
!BOP
! 
! !ROUTINE: ERAinterimLandobsInit
! \label{ERAinterimLandobsInit}
!
! !INTERFACE: 
  subroutine ERAinterimLandobsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod    
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the ERA Interim land data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status
    integer               :: k 

    if(.not.allocated(ERAIlandobs)) then 
       allocate(ERAIlandobs(LVT_rc%nDataStreams))
    endif
   
    call ESMF_ConfigGetAttribute(LVT_Config, ERAIlandobs(i)%odir, &
         label='ERA interim land data directory:', rc=status)
    call LVT_verify(status, 'ERA interim land data directory: not defined')

    ERAIlandobs(i)%gridDesc = 0
        
    ERAIlandobs(i)%nc = 480
    ERAIlandobs(i)%nr = 241

    !filling the items needed by the interpolation library
    ERAIlandobs(i)%gridDesc(1) = 0  
    ERAIlandobs(i)%gridDesc(2) = ERAIlandobs(i)%nc
    ERAIlandobs(i)%gridDesc(3) = ERAIlandobs(i)%nr
    ERAIlandobs(i)%gridDesc(4) = -90.000
    ERAIlandobs(i)%gridDesc(5) = -180.000
    ERAIlandobs(i)%gridDesc(7) = 90.000
    ERAIlandobs(i)%gridDesc(8) = 180.000
    ERAIlandobs(i)%gridDesc(6) = 128
    ERAIlandobs(i)%gridDesc(9) = 0.75
    ERAIlandobs(i)%gridDesc(10) = 0.75
    ERAIlandobs(i)%gridDesc(20) = 0

    ERAIlandobs(i)%datares  = 0.75

    if(LVT_isAtAfinerResolution(ERAIlandobs(i)%datares)) then
       
       allocate(ERAIlandobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(ERAIlandobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(ERAIlandobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input(ERAIlandobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            ERAIlandobs(i)%rlat, &
            ERAIlandobs(i)%rlon, &
            ERAIlandobs(i)%n11)
    else
       allocate(ERAIlandobs(i)%n11(ERAIlandobs(i)%nc*&
            ERAIlandobs(i)%nr))
       call upscaleByAveraging_input(ERAIlandobs(i)%gridDesc,&
            LVT_rc%gridDesc,ERAIlandobs(i)%nc*ERAIlandobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,ERAIlandobs(i)%n11)
    endif

    call ESMF_TimeSet(ERAIlandobs(i)%startTime,yy=1900, &
         mm = 1, &
         dd = 1, &
         h = 0, &
         m = 0, &
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status, 'Error in ESMF_TimeSet: ERAinterimLandobsinit')
    
    call LVT_update_timestep(LVT_rc, 3600)

    ERAIlandobs(i)%da = -1
    ERAIlandobs(i)%startFlag = .true.

  end subroutine ERAinterimLandobsinit

  
end module ERAinterimLandobsMod
