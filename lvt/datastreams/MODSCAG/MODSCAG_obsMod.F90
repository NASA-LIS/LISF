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
! !MODULE: MODSCAG_obsMod
! \label(MODSCAG_obsMod)
!
! !INTERFACE:
module MODSCAG_obsMod
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
!   MODSCAG fractional snow cover product. Note that 
!   this is a resampled data at 0.005deg
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  20 Feb 2018   Sujay Kumar  Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MODSCAG_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: Modscagobs
!EOP
  type, public :: modscagobsdec
     character*100           :: odir
     integer                 :: modis_nc, modis_nr
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
  end type modscagobsdec

  type(modscagobsdec), allocatable :: modscagobs(:)

contains
  
!BOP
! 
! !ROUTINE: MODSCAG_obsInit
! \label{MODSCAG_obsInit}
!
! !INTERFACE: 
  subroutine MODSCAG_obsinit(i)
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
!   for reading the MODSCAG data, including the setup of spatial interpolation
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

    if(.not.allocated(modscagobs)) then 
       allocate(modscagobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, modscagobs(i)%odir, &
         label='MODSCAG observation directory:', rc=status)
    call LVT_verify(status, 'MODSCAG observation directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    cornerlat1 = max(0.0025,nint((LVT_rc%gridDesc(4)-0.0025)/0.005)*0.005+0.0025-2*0.005)
    cornerlon1 = max(0.0025,nint((LVt_rc%gridDesc(5)-0.0025)/0.005)*0.005+0.0025-2*0.005)
    cornerlat2 = min(89.9975,nint((LVT_rc%gridDesc(7)-0.0025)/0.005)*0.005+0.0025+2*0.005)
    cornerlon2 = min(179.9975,nint((LVT_rc%gridDesc(8)-0.0025)/0.005)*0.005+0.0025+2*0.005)
    
    modscagobs(i)%modis_nr = nint((cornerlat2-cornerlat1)/0.005)+1
    modscagobs(i)%modis_nc = nint((cornerlon2-cornerlon1)/0.005)+1

    modscagobs(i)%gridDesci = 0 
    modscagobs(i)%gridDesci(1) = 0 
    modscagobs(i)%gridDesci(2) = modscagobs(i)%modis_nc
    modscagobs(i)%gridDesci(3) = modscagobs(i)%modis_nr
    modscagobs(i)%gridDesci(4) = cornerlat1
    modscagobs(i)%gridDesci(5) = cornerlon1
    modscagobs(i)%gridDesci(7) = cornerlat2
    modscagobs(i)%gridDesci(8) = cornerlon2
    modscagobs(i)%gridDesci(6) = 128
    modscagobs(i)%gridDesci(9) = 0.005
    modscagobs(i)%gridDesci(10) = 0.005
    modscagobs(i)%gridDesci(20) = 64

    if(LVT_isAtAfinerResolution(modscagobs(i)%gridDesci(9))) then
       allocate(modscagobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(modscagobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       
       allocate(modscagobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(modscagobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(modscagobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(modscagobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
       allocate(modscagobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(modscagobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(modscagobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(modscagobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
       

       call bilinear_interp_input(modscagobs(i)%gridDesci,LVT_rc%gridDesc,&
            LVT_rc%lnc*LVT_rc%lnr, &
            modscagobs(i)%rlat, modscagobs(i)%rlon,&
            modscagobs(i)%n11, modscagobs(i)%n12, &
            modscagobs(i)%n21, modscagobs(i)%n22, & 
            modscagobs(i)%w11, modscagobs(i)%w12, &
            modscagobs(i)%w21, modscagobs(i)%w22)
       
    else
       
       allocate(modscagobs(i)%n11(&
            modscagobs(i)%modis_nc*modscagobs(i)%modis_nr))
       call upscaleByAveraging_input(modscagobs(i)%gridDesci, LVT_rc%gridDesc,&
            modscagobs(i)%modis_nc*modscagobs(i)%modis_nr, &
            LVT_rc%lnc*LVT_rc%lnr, modscagobs(i)%n11)
    endif

    modscagobs(i)%startflag = .true. 

    call map_set(PROJ_LATLON, modscagobs(i)%gridDesci(4), modscagobs(i)%gridDesci(5), &
         0.0, modscagobs(i)%gridDesci(9), modscagobs(i)%gridDesci(10), 0.0, &
         modscagobs(i)%modis_nc, modscagobs(i)%modis_nr, modscagobs(i)%mod_proj)
  end subroutine MODSCAG_obsinit


end module MODSCAG_obsMod
