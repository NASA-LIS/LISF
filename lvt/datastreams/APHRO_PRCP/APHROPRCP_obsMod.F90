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
! !MODULE: APHROPRCP_obsMod
! \label(APHROPRCP_obsMod)
!
! !INTERFACE:
module APHROPRCP_obsMod
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
!  APHRODITE (Asian Precipitaton - Highly-Resolved 
!  Observational Data Integration Towards Evaluation)
!  daily gridded precipitation data. 
!
!  https://climatedataguide.ucar.edu/climate-data/aphrodite-asian-precipitation-highly-resolved-observational-data-integration-towards
!  
!  Coverage is from 1951 - 2007, daily at 
!  0.25 deg lat/lon resolution
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  09 Aug 2017   Sujay Kumar  Initial Specification
! 
!EOP

  PUBLIC :: APHROPRCP_obsinit
  PUBLIC :: APHROPRCPobs

  type, public :: aphroprcpobsdec
     character*100        :: odir
     character*100        :: loc
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
  end type aphroprcpobsdec

  type(aphroprcpobsdec), allocatable :: aphroprcpobs(:)

contains
  
!BOP
! 
! !ROUTINE: APHROPRCP_obsInit
! \label{APHROPRCP_obsInit}
!
! !INTERFACE: 
  subroutine APHROPRCP_obsinit(i)
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
!  for reading APHRO PCP data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    real               :: gridDesci(50)
    integer            :: status
    
    if(.not.allocated(aphroprcpobs)) then 
       allocate(aphroprcpobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, aphroprcpobs(i)%odir, &
         label='APHRO PCP data directory:',rc=status)
    call LVT_verify(status, 'APHRO PCP data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_config, aphroprcpobs(i)%loc, &
         label='APHRO PCP data region:',rc=status)
    if(status.ne.0) then 
       write(LVT_logunit,*) '[ERR] APHRO PCP data region: not defined'
       write(LVT_logunit,*) '[ERR] options are: '
       write(LVT_logunit,*) "[ERR] 'MA' (for Monsoon Asia)' "
       write(LVT_logunit,*) "[ERR] 'ME' (for Middle East)' "
       write(LVT_logunit,*) "[ERR] 'RU' (for Northern Eurasia)' "
       write(LVT_logunit,*) "[ERR] 'PR' (for Combined Eurasia)' "
       call LVT_endrun()
    endif

    gridDesci = 0 
    call LVT_update_timestep(LVT_rc, 86400)

    allocate(aphroprcpobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(aphroprcpobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(aphroprcpobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(aphroprcpobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(aphroprcpobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(aphroprcpobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(aphroprcpobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(aphroprcpobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(aphroprcpobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(aphroprcpobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    if(aphroprcpobs(i)%loc.eq."MA") then 
       aphroprcpobs(i)%nc = 360
       aphroprcpobs(i)%nr = 280
       
       gridDesci(1) = 0 
       gridDesci(2) = 360
       gridDesci(3) = 280 
       gridDesci(4) = -14.875
       gridDesci(5) = 60.125
       gridDesci(6) = 128 
       gridDesci(7) = 54.875
       gridDesci(8) = 149.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64
       
       call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
            LVT_rc%lnc*LVT_rc%lnr, &
            aphroprcpobs(i)%rlat, aphroprcpobs(i)%rlon, &  
            aphroprcpobs(i)%n11, aphroprcpobs(i)%n12,   & 
            aphroprcpobs(i)%n21, aphroprcpobs(i)%n22,   & 
            aphroprcpobs(i)%w11, aphroprcpobs(i)%w12,   & 
            aphroprcpobs(i)%w21, aphroprcpobs(i)%w22)
    else
       write(LVT_logunit,*) "[ERR] The Aphrodite plugin only supports"
       write(LVT_logunit,*) "[ERR] the MA region currently"
       call LVT_endrun()
    endif
    
    allocate(aphroprcpobs(i)%rainf(aphroprcpobs(i)%nc,&
         aphroprcpobs(i)%nr,&
         366))
    aphroprcpobs(i)%yr = -1
    
    call ESMF_TimeSet(aphroprcpobs(i)%startTime,  yy=LVT_rc%yr, &
         mm = 1, &
         dd = 1, &
         h = 0, &
         m = 0, &
         calendar = LVT_calendar, &
         rc=status)
    call LVT_verify(status, 'error in setting APHRO start time')

!    call ESMF_TimeSet(aphroprcpobs(i)%stopTime, yy=eyr, &
!         mm = emo, &
!         dd = eda, &
!         h = ehr, &
!         m = emn, &
!         calendar = LVT_calendar, &
!         rc=status)
!    call LVT_verify(status, 'error in setting scan stop time')

    call ESMF_TimeIntervalSet(aphroprcpobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (aphroprcpobs)')
        
  end subroutine APHROPRCP_obsinit


end module APHROPRCP_obsMod
