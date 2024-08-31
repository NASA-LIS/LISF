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
! !MODULE: ESACCIsm_obsMod
! \label(ESACCIsm_obsMod)
!
! !INTERFACE:
module ESACCIsm_obsMod
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
!  This module handles the observation plugin for the Land Parameter
!  Retrieval Model (LPRM) AMSR-E soil moisture product
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
! 
!EOP
! 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ESACCIsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ESACCIsmobs
!EOP
  type, public :: esaccismobsdec
     character*100              :: odir
     integer                    :: version
     real                       :: gridDesci(50)
     integer                    :: esaccinc, esaccinr
     type(proj_info)            :: esacciproj
     integer, allocatable       :: n11(:)
     integer, allocatable       :: n12(:)
     integer, allocatable       :: n21(:)
     integer, allocatable       :: n22(:)
     real,  allocatable         :: w11(:)
     real,  allocatable         :: w12(:)
     real,  allocatable         :: w21(:)
     real,  allocatable         :: w22(:)
     real,  allocatable         :: rlat(:)
     real,  allocatable         :: rlon(:)
     real,    allocatable       :: smobs(:,:)
     logical                    :: startmode     

  end type esaccismobsdec

  type(esaccismobsdec), allocatable:: ESACCIsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: ESACCIsm_obsInit
! \label{ESACCIsm_obsInit}
!
! !INTERFACE: 
  subroutine ESACCIsm_obsinit(i)
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
!  for reading ESACCI AMSRE soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: npts, status
    real                  :: gridDesci(50)

    if(.not.allocated(ESACCIsmobs)) then 
       allocate(ESACCIsmobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, ESACCIsmobs(i)%odir, &
         label='ESA CCI soil moisture data directory:', rc=status)
    call LVT_verify(status, 'ESA CCI soil moisture data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_Config, ESACCIsmobs(i)%version, &
         label='ESA CCI soil moisture data version:', rc=status)
    call LVT_verify(status, 'ESA CCI soil moisture data version: not defined')
    
    call LVT_update_timestep(LVT_rc, 86400)

    ESACCIsmobs(i)%startmode = .true. 

    allocate(ESACCIsmobs(i)%smobs(LVT_rc%lnc,LVT_rc%lnr))

    ESACCIsmobs(i)%esaccinc = 1440
    ESACCIsmobs(i)%esaccinr = 720
    
    call map_set(PROJ_LATLON, -89.875,-179.875,&
         0.0, 0.25,0.25, 0.0,&
         ESACCIsmobs(i)%esaccinc,ESACCIsmobs(i)%esaccinr,&
         ESACCIsmobs(i)%esacciproj)
    
    gridDesci = 0 
    gridDesci(1) = 0 
    gridDesci(2) = 1440
    gridDesci(3) = 720
    gridDesci(4) = -89.875
    gridDesci(5) = -179.875
    gridDesci(6) = 128
    gridDesci(7) = 89.875
    gridDesci(8) = 179.875
    gridDesci(9) = 0.25
    gridDesci(10) = 0.25
    gridDesci(20) = 64
    
    allocate(ESACCIsmobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ESACCIsmobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(ESACCIsmobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ESACCIsmobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ESACCIsmobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ESACCIsmobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(ESACCIsmobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ESACCIsmobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ESACCIsmobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(ESACCIsmobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,ESACCIsmobs(i)%rlat, &
         ESACCIsmobs(i)%rlon,ESACCIsmobs(i)%n11, &
         ESACCIsmobs(i)%n12, ESACCIsmobs(i)%n21, &
         ESACCIsmobs(i)%n22, ESACCIsmobs(i)%w11, &
         ESACCIsmobs(i)%w12, ESACCIsmobs(i)%w21, &
         ESACCIsmobs(i)%w22)

  end subroutine ESACCIsm_obsinit


end module ESACCIsm_obsMod
