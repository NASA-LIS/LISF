!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: UASWE_obsMod
! \label(UASWE_obsMod)
!
! !INTERFACE:
module UASWE_obsMod
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
!  University of Arizona (UA) SWE data.
!  The UA SWE data is provides in the NAD 1983 grid with 
!  4km resolution. The domain extents are from (24, -125) to (50, -66.5). 
!  The data entries are 16-bit signed integers.
!  
!  Temporal coverage is from 1 Jan 1982 - 2017 
! 
! !FILES USED:
! 
!EOP

  PUBLIC :: UASWE_obsinit
  PUBLIC :: uasweobs

  type, public :: uasweobsdec
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
     real,    allocatable     :: swe(:,:,:)
  end type uasweobsdec

  type(uasweobsdec), allocatable :: uasweobs(:)

contains
  
!BOP
! 
! !ROUTINE: UASWE_obsinit
! \label{UASWE_obsinit}
!
! !INTERFACE: 
  subroutine UASWE_obsinit(i)
! 
! !USES: 
    use ESMF
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
!  for reading UA SWE data. 
! 
! !FILES USED:
!
!EOP

    real               :: gridDesci(50)
    integer            :: status
    
    if(.not.allocated(uasweobs)) then 
       allocate(uasweobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, uasweobs(i)%odir, &
         label='UA SWE observation directory:',rc=status)
    call LVT_verify(status, 'UA SWE observation directory: not defined')

    gridDesci = 0 
    call LVT_update_timestep(LVT_rc, 86400)

    allocate(uasweobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasweobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(uasweobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasweobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasweobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasweobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(uasweobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasweobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasweobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasweobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    uasweobs(i)%nc = 1405
    uasweobs(i)%nr = 621
    
    gridDesci(1) = 0 
    gridDesci(2) = 1405
    gridDesci(3) = 621
    gridDesci(4) = 24.0833340 
    gridDesci(5) = -125.0000
    gridDesci(6) = 128
    gridDesci(7) = 49.9166679
    gridDesci(8) = -66.5000
    gridDesci(9) = 0.04166662697178698
    gridDesci(10) = 0.04166662697178698
    gridDesci(20) = 64
    
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         uasweobs(i)%rlat,  uasweobs(i)%rlon, &  
         uasweobs(i)%n11, uasweobs(i)%n12,   & 
         uasweobs(i)%n21, uasweobs(i)%n22,   & 
         uasweobs(i)%w11, uasweobs(i)%w12,   & 
         uasweobs(i)%w21, uasweobs(i)%w22)

    allocate(uasweobs(i)%swe(uasweobs(i)%nc,&
         uasweobs(i)%nr,366))
    
    call ESMF_TimeIntervalSet(uasweobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (uasweobs)')
        
  end subroutine UASWE_obsinit


end module UASWE_obsMod
