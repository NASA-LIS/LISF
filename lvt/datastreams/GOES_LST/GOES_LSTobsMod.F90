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
! !MODULE: GOES_LSTobsMod
! \label(GOES_LSTobsMod)
!
! !INTERFACE:
module GOES_LSTobsMod
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
  PUBLIC :: GOES_LSTobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GOESLSTobs !Object to hold GOESLST observation attributes
!EOP

  type, public :: goeslstdec
     character*100           :: odir
     integer                 :: nc, nr
     integer, allocatable        :: n11(:)
     real,    allocatable        :: lst(:)
     real                    :: gridDesc(50)
  end type goeslstdec
     
  type(goeslstdec), allocatable :: GOESLSTObs(:)

contains
  
!BOP
! 
! !ROUTINE: GOES_LSTobsInit
! \label{GOES_LSTobsInit}
!
! !INTERFACE: 
  subroutine GOES_LSTobsinit(i)
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
!   for reading the GOES LST data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    if(.not.allocated(GOESLSTobs)) then 
       allocate(GOESLSTobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, GOESLSTobs(i)%odir, &
         label='GOES LST data directory:',rc=status)
    call LVT_verify(status, 'GOES LST data directory: not defined')

    call LVT_update_timestep(LVT_rc, 10800)

    allocate(Goeslstobs(I)%lst(LVT_rc%lnc*LVT_rc%lnr))

    goeslstobs(i)%gridDesc = 0
    
  end subroutine GOES_LSTobsinit


end module GOES_LSTobsMod
