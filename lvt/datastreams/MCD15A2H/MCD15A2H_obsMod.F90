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
! !MODULE: MCD15A2H_obsMod
! \label(MCD15A2H_obsMod)
!
! !INTERFACE:
module MCD15A2H_obsMod
! 
! !USES:   
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!   This subroutine provides the observation plugin for reading the 
!   MCD15A2H 500m LAI data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  20 Oct 2020   Sujay Kumar  Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MCD15A2H_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: Mcd15a2hobs
!EOP
  type, public :: mcd15a2hobsdec
     character*100           :: odir
     integer                 :: nc, nr
     real, allocatable           :: rlat(:)
     real, allocatable           :: rlon(:)
     integer, allocatable        :: n11(:)
     integer, allocatable        :: n12(:)
     integer, allocatable        :: n21(:)
     integer, allocatable        :: n22(:)     
     real,    allocatable        :: w11(:)
     real,    allocatable        :: w12(:)
     real,    allocatable        :: w21(:)
     real,    allocatable        :: w22(:)
     logical                 :: startFlag
     type(proj_info)         :: mod_proj
    real                 :: gridDesci(50)
  end type mcd15a2hobsdec

  type(mcd15a2hobsdec), allocatable :: mcd15a2hobs(:)

contains
  
!BOP
! 
! !ROUTINE: MCD15A2H_obsInit
! \label{MCD15A2H_obsInit}
!
! !INTERFACE: 
  subroutine MCD15A2H_obsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_logMod
    use LVT_histDataMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!
!   This subroutine initializes and sets up the data structures required
!   for reading the MCD15A2H data, including the setup of spatial interpolation
!   weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer              :: status
    integer              :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real                 :: upgmt
    real                 :: cornerlat1, cornerlat2, cornerlon1, cornerlon2

    if(.not.allocated(mcd15a2hobs)) then 
       allocate(mcd15a2hobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, mcd15a2hobs(i)%odir, &
         label='MCD15A2H data directory:', rc=status)
    call LVT_verify(status, 'MCD15A2H data directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    cornerlat1 = max(-59.9978927,nint((LVT_rc%gridDesc(4)+59.9978927)/0.00416667)*0.00416667-59.9978927-2*0.00416667)
    cornerlon1 = max(-179.9979167,nint((LVt_rc%gridDesc(5)+179.9979167)/0.00416667)*0.00416667-179.9979167-2*0.00416667)
    cornerlat2 = min(89.9979167,nint((LVT_rc%gridDesc(7)+59.9978927)/0.00416667)*0.01-59.9978927+2*0.00416667)
    cornerlon2 = min(179.9979167,nint((LVT_rc%gridDesc(8)+179.9979167)/0.00416667)*0.00416667-179.9979167+2*0.00416667)
    
    mcd15a2hobs(i)%nr = nint((cornerlat2-cornerlat1)/0.00416667)+1
    mcd15a2hobs(i)%nc = nint((cornerlon2-cornerlon1)/0.00416667)+1

    mcd15a2hobs(i)%gridDesci = 0 
    mcd15a2hobs(i)%gridDesci(1) = 0 
    mcd15a2hobs(i)%gridDesci(2) = mcd15a2hobs(i)%nc
    mcd15a2hobs(i)%gridDesci(3) = mcd15a2hobs(i)%nr
    mcd15a2hobs(i)%gridDesci(4) = cornerlat1
    mcd15a2hobs(i)%gridDesci(5) = cornerlon1
    mcd15a2hobs(i)%gridDesci(7) = cornerlat2
    mcd15a2hobs(i)%gridDesci(8) = cornerlon2
    mcd15a2hobs(i)%gridDesci(6) = 128
    mcd15a2hobs(i)%gridDesci(9) = 0.00416667
    mcd15a2hobs(i)%gridDesci(10) = 0.00416667
    mcd15a2hobs(i)%gridDesci(20) = 64

    if(LVT_isAtAfinerResolution(mcd15a2hobs(i)%gridDesci(9))) then
       allocate(mcd15a2hobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(mcd15a2hobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       
       allocate(mcd15a2hobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(mcd15a2hobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(mcd15a2hobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(mcd15a2hobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
       allocate(mcd15a2hobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(mcd15a2hobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(mcd15a2hobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(mcd15a2hobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
       

       call bilinear_interp_input(mcd15a2hobs(i)%gridDesci,LVT_rc%gridDesc,&
            LVT_rc%lnc*LVT_rc%lnr, &
            mcd15a2hobs(i)%rlat, mcd15a2hobs(i)%rlon,&
            mcd15a2hobs(i)%n11, mcd15a2hobs(i)%n12, &
            mcd15a2hobs(i)%n21, mcd15a2hobs(i)%n22, & 
            mcd15a2hobs(i)%w11, mcd15a2hobs(i)%w12, &
            mcd15a2hobs(i)%w21, mcd15a2hobs(i)%w22)
       
    else
       
       allocate(mcd15a2hobs(i)%n11(&
            mcd15a2hobs(i)%nc*mcd15a2hobs(i)%nr))
       call upscaleByAveraging_input(mcd15a2hobs(i)%gridDesci, LVT_rc%gridDesc,&
            mcd15a2hobs(i)%nc*mcd15a2hobs(i)%nr, &
            LVT_rc%lnc*LVT_rc%lnr, mcd15a2hobs(i)%n11)
    endif

    mcd15a2hobs(i)%startflag = .true. 

    call map_set(PROJ_LATLON, mcd15a2hobs(i)%gridDesci(4), mcd15a2hobs(i)%gridDesci(5), &
         0.0, mcd15a2hobs(i)%gridDesci(9), mcd15a2hobs(i)%gridDesci(10), 0.0, &
         mcd15a2hobs(i)%nc, mcd15a2hobs(i)%nr, mcd15a2hobs(i)%mod_proj)
  end subroutine MCD15A2H_obsinit


end module MCD15A2H_obsMod
