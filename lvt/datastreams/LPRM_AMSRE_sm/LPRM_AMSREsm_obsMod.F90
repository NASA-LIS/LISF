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
! !MODULE: LPRM_AMSREsm_obsMod
! \label(LPRM_AMSREsm_obsMod)
!
! !INTERFACE:
module LPRM_AMSREsm_obsMod
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
  PUBLIC :: LPRM_AMSREsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: LPRM_AMSREsmobs
!EOP
  type, public :: lprmamsresmobsdec

     character*100           :: odir
     integer                 :: rawdata
     character*10            :: version
     character*10            :: channel
     real                    :: gridDesci(50)
     integer                :: lprmnc, lprmnr
     type(proj_info)        :: lprmproj
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
     real,    allocatable    :: smobs(:,:)
     logical             :: startmode     

  end type lprmamsresmobsdec

  type(lprmamsresmobsdec), allocatable:: LPRM_AMSREsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: LPRM_AMSREsm_obsInit
! \label{LPRM_AMSREsm_obsInit}
!
! !INTERFACE: 
  subroutine LPRM_AMSREsm_obsinit(i)
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
!  for reading LPRM AMSRE soil moisture data. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: npts, status
    real                  :: gridDesci(50)


    if(.not.allocated(LPRM_AMSREsmobs)) then 
       allocate(LPRM_AMSREsmobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_Config, LPRM_AMSREsmobs(i)%odir, &
         label='LPRM AMSR-E soil moisture observation directory:', rc=status)
    call LVT_verify(status, &
         'LPRM AMSR-E soil moisture observation directory: not defined')
    

    call ESMF_ConfigGetAttribute(LVT_config, LPRM_AMSREsmobs(i)%version, &
         label='LPRM AMSR-E soil moisture data version:', rc=status)
    if(status.ne.0) then 
       write(LVT_logunit,*) "[ERR] 'LPRM AMSR-E soil moisture data version:' not defined"
       write(LVT_logunit,*) "[ERR] supported options are: "
       write(LVT_logunit,*) "[ERR] 'V05' or 'GES-DISC'"
       call LVT_endrun()
    endif

    if(LPRM_AMSREsmobs(i)%version.eq."GES-DISC") then 
       call ESMF_ConfigGetAttribute(LVT_config, LPRM_AMSREsmobs(i)%channel, &
            label='LPRM AMSR-E soil moisture data channel:', rc=status)
       if(status.ne.0) then 
          write(LVT_logunit,*) "[ERR] 'LPRM AMSR-E soil moisture data channel:' not defined"
          write(LVT_logunit,*) "[ERR] supported options are: "
          write(LVT_logunit,*) "[ERR] 'X-band' or 'C-band'"
          call LVT_endrun()
       endif
    elseif(LPRM_AMSREsmobs(i)%version.eq."V05") then 
       call ESMF_ConfigGetAttribute(LVT_config, LPRM_AMSREsmobs(i)%rawdata, &
            label='LPRM AMSR-E soil moisture use raw data:', rc=status)
       call LVT_verify(status, 'LPRM AMSR-E soil moisture use raw data: not defined')
       
    endif

    LPRM_AMSREsmobs(i)%startmode = .true. 

    allocate(LPRM_AMSREsmobs(i)%smobs(LVT_rc%lnc,LVT_rc%lnr))

    LPRM_AMSREsmobs(i)%lprmnc = 1440
    LPRM_AMSREsmobs(i)%lprmnr = 720
    
    call map_set(PROJ_LATLON, -89.875,-179.875,&
         0.0, 0.25,0.25, 0.0,&
         LPRM_AMSREsmobs(i)%lprmnc,LPRM_AMSREsmobs(i)%lprmnr,&
         LPRM_AMSREsmobs(i)%lprmproj)
    
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
    
    allocate(LPRM_AMSREsmobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(LPRM_AMSREsmobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(LPRM_AMSREsmobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(LPRM_AMSREsmobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(LPRM_AMSREsmobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(LPRM_AMSREsmobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))
    
    allocate(LPRM_AMSREsmobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(LPRM_AMSREsmobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(LPRM_AMSREsmobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(LPRM_AMSREsmobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))
    
    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,LPRM_AMSREsmobs(i)%rlat, &
         LPRM_AMSREsmobs(i)%rlon,LPRM_AMSREsmobs(i)%n11, &
         LPRM_AMSREsmobs(i)%n12, LPRM_AMSREsmobs(i)%n21, &
         LPRM_AMSREsmobs(i)%n22, LPRM_AMSREsmobs(i)%w11, &
         LPRM_AMSREsmobs(i)%w12, LPRM_AMSREsmobs(i)%w21, &
         LPRM_AMSREsmobs(i)%w22)

    call LVT_update_timestep(LVT_rc, 86400)

  end subroutine LPRM_AMSREsm_obsinit


end module LPRM_AMSREsm_obsMod
