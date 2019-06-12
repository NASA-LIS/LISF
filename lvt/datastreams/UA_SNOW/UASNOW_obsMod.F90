!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: UASNOW_obsMod
! \label(UASNOW_obsMod)
!
! !INTERFACE:
module UASNOW_obsMod
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
!  University of Arizona (UA) SWE/Snow Depth data.
!  The UA Snow data is provides in the NAD 1983 grid with 
!  4km resolution. The domain extents are from (24, -125) to (50, -66.5). 
!  The data entries are 16-bit signed integers.
!  
!  Temporal coverage is from 1 Jan 1982 - 2017 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  28 May 2019: Rhae Sung Kim, Initial Specification!
!EOP

  PUBLIC :: UASNOW_obsinit
  PUBLIC :: uasnowobs

  type, public :: uasnowobsdec
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
     real,    allocatable     :: snwd(:,:,:)
  end type uasnowobsdec

  type(uasnowobsdec), allocatable :: uasnowobs(:)

contains
  
!BOP
! 
! !ROUTINE: UASNOW_obsinit
! \label{UASNOW_obsinit}
!
! !INTERFACE: 
  subroutine UASNOW_obsinit(i)
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
!  for reading UA Snow data. 
! 
! !FILES USED:
!
!EOP

    real               :: gridDesci(50)
    integer            :: status
    
    if(.not.allocated(uasnowobs)) then 
       allocate(uasnowobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, uasnowobs(i)%odir, &
         label='UA Snow observation directory:',rc=status)
    call LVT_verify(status, 'UA Snow observation directory: not defined')

    gridDesci = 0 
    call LVT_update_timestep(LVT_rc, 86400)

    allocate(uasnowobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnowobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(uasnowobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnowobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnowobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnowobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(uasnowobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnowobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnowobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnowobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    uasnowobs(i)%nc = 1405
    uasnowobs(i)%nr = 621
    
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
         uasnowobs(i)%rlat,  uasnowobs(i)%rlon, &  
         uasnowobs(i)%n11, uasnowobs(i)%n12,   & 
         uasnowobs(i)%n21, uasnowobs(i)%n22,   & 
         uasnowobs(i)%w11, uasnowobs(i)%w12,   & 
         uasnowobs(i)%w21, uasnowobs(i)%w22)

    allocate(uasnowobs(i)%swe(uasnowobs(i)%nc,&
         uasnowobs(i)%nr,366))
    allocate(uasnowobs(i)%snwd(uasnowobs(i)%nc,&
         uasnowobs(i)%nr,366))    
 
    call ESMF_TimeIntervalSet(uasnowobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (uasnowobs)')
        
  end subroutine UASNOW_obsinit


end module UASNOW_obsMod
