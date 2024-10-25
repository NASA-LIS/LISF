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
! !MODULE: UWET_obsMod
! \label(UWET_obsMod)
!
! !INTERFACE:
module UWET_obsMod
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
!  19 Jul 2013   Sujay Kumar  Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: UWET_obsinit !Initializes structures for reading UWET data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: UWETobs !Object to hold UWET observation attributes
!EOP

  type, public :: uwetdec
     character*100           :: odir
     integer                 :: nc, nr
     integer, allocatable        :: n11(:)
     real,    allocatable        :: qle(:)
     integer                 :: yr
     integer                 :: mo
     real                    :: gridDesc(50)
     logical                 :: startFlag
  end type uwetdec
     
  type(uwetdec), allocatable :: UWETObs(:)

contains
  
!BOP
! 
! !ROUTINE: UWET_obsInit
! \label{UWET_obsInit}
!
! !INTERFACE: 
  subroutine UWET_obsinit(i)
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
!   for reading the UWET data, including the computation of spatial 
!   interpolation weights.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status
    real                  :: cornerlat1, cornerlat2
    real                  :: cornerlon1, cornerlon2

    if(.not.allocated(UWETobs)) then 
       allocate(UWETobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, UWETObs(i)%odir, &
         label='UW ET data directory: ',rc=status)
    call LVT_verify(status, 'UW ET data directory: not defined')

    allocate(UWETobs(i)%qle(LVT_rc%lnc*LVT_rc%lnr))

    uwetobs(i)%gridDesc = 0
    
    uwetobs(i)%nr = 560
    uwetobs(i)%nc = 1160

    allocate(Uwetobs(I)%n11(uwetobs(i)%nc*uwetobs(i)%nr))

    !filling the items needed by the interpolation library
    uwetobs(i)%gridDesc(1) = 0  
    uwetobs(i)%gridDesc(2) = uwetobs(i)%nc
    uwetobs(i)%gridDesc(3) = uwetobs(i)%nr
    uwetobs(i)%gridDesc(4) = 25.025
    uwetobs(i)%gridDesc(5) = -124.975
    uwetobs(i)%gridDesc(7) = 52.975
    uwetobs(i)%gridDesc(8) = -67.025
    uwetobs(i)%gridDesc(6) = 128
    uwetobs(i)%gridDesc(9) = 0.05
    uwetobs(i)%gridDesc(10) = 0.05
    uwetobs(i)%gridDesc(20) = 64

    call upscaleByAveraging_input(uwetobs(i)%gridDesc,&
         LVT_rc%gridDesc,uwetobs(i)%nc*uwetobs(i)%nr,&
         LVT_rc%lnc*LVT_rc%lnr,uwetobs(i)%n11)

    UWETobs(i)%mo = LVT_rc%mo
    UWETobs(i)%yr = -1
    uwetobs(i)%startFlag = .true.

    if(LVT_rc%tavgInterval.lt.2592000) then 
       write(LVT_logunit,*) '[ERR] The time averaging interval must be greater than'
       write(LVT_logunit,*) '[ERR] equal to a month since the UWET data is monthly'
       call LVT_endrun()
    endif
    call LVT_update_timestep(LVT_rc, 2592000)
  end subroutine UWET_obsinit


end module UWET_obsMod
