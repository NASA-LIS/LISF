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
! !MODULE: GLASSalbedoobsMod
! \label(GLASSalbedoobsMod)
!
! !INTERFACE:
module GLASSalbedoobsMod
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
! 14 Sep 2018   Abheera Hazra - Adapted GLASSlai to GLASSalbedo 
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GLASSalbedoobsinit !Initializes structures for reading MOD16A2 data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GLASSalbedoobs !Object to hold GLASSalbedo observation attributes
!EOP

  type, public :: glassalbedodec
     character*200           :: odir
     character*200           :: source
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

  end type glassalbedodec
     
  type(glassalbedodec), save :: GLASSalbedoObs(2)

contains
  
!BOP
! 
! !ROUTINE: GLASSalbedoobsInit
! \label{GLASSalbedoobsInit}
!
! !INTERFACE: 
  subroutine GLASSalbedoobsinit(i)
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

    call ESMF_ConfigGetAttribute(LVT_Config, GLASSalbedoobs(i)%odir, &
         label='GLASS ALBEDO data directory:',rc=status)
    call LVT_verify(status, 'GLASS ALBEDO data directory: not defined')

! source = "AVHRR" or "MODIS"
    call ESMF_ConfigGetAttribute(LVT_Config, GLASSalbedoobs(i)%source, &
         label='GLASS ALBEDO data source:',rc=status)
    call LVT_verify(status, 'GLASS ALBEDO data source: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    glassalbedoobs(i)%gridDesc = 0
        
    glassalbedoobs(i)%nc = 7200
    glassalbedoobs(i)%nr = 3600

    !filling the items needed by the interpolation library
    glassalbedoobs(i)%gridDesc(1) = 0  
    glassalbedoobs(i)%gridDesc(2) = glassalbedoobs(i)%nc
    glassalbedoobs(i)%gridDesc(3) = glassalbedoobs(i)%nr
    glassalbedoobs(i)%gridDesc(4) = -89.875
    glassalbedoobs(i)%gridDesc(5) = -179.875
    glassalbedoobs(i)%gridDesc(7) = 89.875
    glassalbedoobs(i)%gridDesc(8) = 179.875
    glassalbedoobs(i)%gridDesc(6) = 128
    glassalbedoobs(i)%gridDesc(9) = 0.05
    glassalbedoobs(i)%gridDesc(10) = 0.05
    glassalbedoobs(i)%gridDesc(20) = 64

    glassalbedoobs(i)%datares  = 0.05

    if(LVT_isAtAfinerResolution(glassalbedoobs(i)%datares)) then
       
       allocate(glassalbedoobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glassalbedoobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glassalbedoobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glassalbedoobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glassalbedoobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glassalbedoobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glassalbedoobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glassalbedoobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glassalbedoobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
       allocate(glassalbedoobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
       
       call bilinear_interp_input(glassalbedoobs(i)%gridDesc, &
            LVT_rc%gridDesc, &
            LVT_rc%lnc*LVT_rc%lnr, &
            glassalbedoobs(i)%rlat, &
            glassalbedoobs(i)%rlon, &
            glassalbedoobs(i)%n11,&
            glassalbedoobs(i)%n12,&
            glassalbedoobs(i)%n21,&
            glassalbedoobs(i)%n22,&
            glassalbedoobs(i)%w11,&
            glassalbedoobs(i)%w12,&
            glassalbedoobs(i)%w21,&
            glassalbedoobs(i)%w22)
    else
       allocate(glassalbedoobs(i)%n11(glassalbedoobs(i)%nc*glassalbedoobs(i)%nr))
       call upscaleByAveraging_input(glassalbedoobs(i)%gridDesc,&
            LVT_rc%gridDesc,glassalbedoobs(i)%nc*glassalbedoobs(i)%nr,&
            LVT_rc%lnc*LVT_rc%lnr,glassalbedoobs(i)%n11)
    endif

    glassalbedoobs(i)%startFlag = .false. 

  end subroutine GLASSalbedoobsinit


end module GLASSalbedoobsMod
