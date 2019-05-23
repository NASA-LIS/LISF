!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!BOP
! 
! !MODULE: UASNWD_obsMod
! \label(UASNWD_obsMod)
!
! !INTERFACE:
module UASNWD_obsMod
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
!  University of Arizona (UA) Snow Depth (SNWD) data.
!  The UA SNWD data is provides in the NAD 1983 grid with 
!  4km resolution. The domain extents are from (24, -125) to (50, -66.5). 
!  The data entries are 16-bit signed integers.
!  
!  Temporal coverage is from 1 Jan 1982 - 2017 
! 
! !FILES USED:
! 
!EOP

  PUBLIC :: UASNWD_obsinit
  PUBLIC :: uasnwdobs

  type, public :: uasnwdobsdec
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
     real,    allocatable     :: snwd(:,:,:)
  end type uasnwdobsdec

  type(uasnwdobsdec), allocatable :: uasnwdobs(:)

contains
  
!BOP
! 
! !ROUTINE: UASNWD_obsinit
! \label{UASNWD_obsinit}
!
! !INTERFACE: 
  subroutine UASNWD_obsinit(i)
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
!  for reading UA SNWD data. 
! 
! !FILES USED:
!
!EOP

    real               :: gridDesci(50)
    integer            :: status
    
    if(.not.allocated(uasnwdobs)) then 
       allocate(uasnwdobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, uasnwdobs(i)%odir, &
         label='UA SNWD observation directory:',rc=status)
    call LVT_verify(status, 'UA SNWD observation directory: not defined')

    gridDesci = 0 
    call LVT_update_timestep(LVT_rc, 86400)

    allocate(uasnwdobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnwdobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(uasnwdobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnwdobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnwdobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnwdobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(uasnwdobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnwdobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnwdobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(uasnwdobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    uasnwdobs(i)%nc = 1405
    uasnwdobs(i)%nr = 621
    
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
         uasnwdobs(i)%rlat,  uasnwdobs(i)%rlon, &  
         uasnwdobs(i)%n11, uasnwdobs(i)%n12,   & 
         uasnwdobs(i)%n21, uasnwdobs(i)%n22,   & 
         uasnwdobs(i)%w11, uasnwdobs(i)%w12,   & 
         uasnwdobs(i)%w21, uasnwdobs(i)%w22)

    allocate(uasnwdobs(i)%snwd(uasnwdobs(i)%nc,&
         uasnwdobs(i)%nr,366))
    
    call ESMF_TimeIntervalSet(uasnwdobs(i)%timestep, s=86400, rc=status)
    call LVT_verify(status, 'error in setting timestep (uasnwdobs)')
        
  end subroutine UASNWD_obsinit


end module UASNWD_obsMod
