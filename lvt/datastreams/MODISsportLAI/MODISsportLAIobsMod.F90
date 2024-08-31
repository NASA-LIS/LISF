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
! !MODULE: MODISsportLAIobsMod
! \label(MODISsportLAIobsMod)
!
! !INTERFACE:
module MODISsportLAIobsMod
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
  PUBLIC :: MODISsportLAIobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: MODISsportLAIobs !Object to hold MODISsportLAI observation attributes
!EOP

  type, public :: MODISsportLAIdec
     character*100           :: odir
     integer                 :: nc, nr
     real                    :: gridDesc(50)
     logical                 :: startFlag
     real                    :: datares

     real,    allocatable :: rlat(:)
     real,    allocatable :: rlon(:)

     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)

  end type MODISsportLAIdec
     
  type(MODISsportLAIdec), save :: MODISsportLAIObs(2)

contains
  
!BOP
! 
! !ROUTINE: MODISsportLAIobsInit
! \label{MODISsportLAIobsInit}
!
! !INTERFACE: 
  subroutine MODISsportLAIobsinit(i)
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
!   for reading the GIMMSAVHRR NDVI data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    call ESMF_ConfigGetAttribute(LVT_Config, MODISsportLAIobs(i)%odir, &
         label='MODIS SPORT LAI data directory:',rc=status)
    call LVT_verify(status, 'MODIS SPORT LAI data directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    MODISsportLAIobs(i)%gridDesc = 0
        
    MODISsportLAIobs(i)%nc = 1456
    MODISsportLAIobs(i)%nr = 625

    !filling the items needed by the interpolation library
    MODISsportLAIobs(i)%gridDesc(1) = 0  
    MODISsportLAIobs(i)%gridDesc(2) = MODISsportLAIobs(i)%nc
    MODISsportLAIobs(i)%gridDesc(3) = MODISsportLAIobs(i)%nr
    MODISsportLAIobs(i)%gridDesc(4) = 24.82
    MODISsportLAIobs(i)%gridDesc(5) = -125.02
    MODISsportLAIobs(i)%gridDesc(7) = 49.78
    MODISsportLAIobs(i)%gridDesc(8) = -66.82
    MODISsportLAIobs(i)%gridDesc(6) = 128
    MODISsportLAIobs(i)%gridDesc(9) = 0.04
    MODISsportLAIobs(i)%gridDesc(10) = 0.04
    MODISsportLAIobs(i)%gridDesc(20) = 64

    MODISsportLAIobs(i)%datares  = 0.04

    if(LVT_isAtAfinerResolution(MODISsportLAIobs(i)%datares)) then
       
       allocate(MODISsportLAIobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(MODISsportLAIobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(MODISsportLAIobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(MODISsportLAIobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(MODISsportLAIobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(MODISsportLAIobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
       allocate(MODISsportLAIobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(MODISsportLAIobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(MODISsportLAIobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(MODISsportLAIobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
       
       call bilinear_interp_input(MODISsportLAIobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            MODISsportLAIobs(i)%rlat, &
            MODISsportLAIobs(i)%rlon, &
            MODISsportLAIobs(i)%n11,&
            MODISsportLAIobs(i)%n12,&
            MODISsportLAIobs(i)%n21,&
            MODISsportLAIobs(i)%n22,&
            MODISsportLAIobs(i)%w11,&
            MODISsportLAIobs(i)%w12,&
            MODISsportLAIobs(i)%w21,&
            MODISsportLAIobs(i)%w22)
    else
       allocate(MODISsportLAIobs(i)%n11(&
            MODISsportLAIobs(i)%nc*MODISsportLAIobs(i)%nr))
       call upscaleByAveraging_input(MODISsportLAIobs(i)%gridDesc,&
            LVT_rc%gridDesc,MODISsportLAIobs(i)%nc*MODISsportLAIobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,MODISsportLAIobs(i)%n11)
    endif

    MODISsportLAIobs(i)%startFlag = .false. 

  end subroutine MODISsportLAIobsinit


end module MODISsportLAIobsMod
