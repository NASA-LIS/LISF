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
! !ROUTINE: ISCCP_TskinobsMod
! \label(ISCCP_TskinobsMod)
!
! !INTERFACE:
module ISCCP_TskinobsMod
! 
! !USES:   
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
! 
!EOP

  PRIVATE 

  PUBLIC :: ISCCP_TskinobsInit

  PUBLIC :: ISCCP_Tskin_obs
  
  type, public :: isccptskinobsdec
     character*100  :: odir
     integer        :: mi
     real, allocatable  :: rlat(:)
     real, allocatable  :: rlon(:)
     integer, allocatable :: n11(:)
  end type isccptskinobsdec

  type(isccptskinobsdec), allocatable  :: isccp_tskin_obs(:)

contains

!BOP
! 
! !ROUTINE: ISCCP_TskinobsInit
! \label(ISCCP_TskinobsInit)
!
! !INTERFACE:
  subroutine ISCCP_TskinobsInit(i)
! 
! !USES:   
    use ESMF
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN)    :: i  !index of the observation type
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    integer                  :: rc
    real                     :: gridDesci(50)
    type(ESMF_Config)        :: modelSpecConfig

    if(.not.allocated(isccp_tskin_obs)) then 
       allocate(isccp_tskin_obs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, isccp_tskin_obs(i)%odir, &
         label='ISCCP Tskin data directory:',rc=rc)
    call LVT_verify(rc, 'ISCCP Tskin data directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)

    gridDesci = 0 
    
    gridDesci(1) = 0 
    gridDesci(2) = 360
    gridDesci(3) = 180
    gridDesci(4) = -89.50
    gridDesci(5) = -179.50
    gridDesci(6) = 128
    gridDesci(7) = 89.50
    gridDesci(8) = 179.50
    gridDesci(9) = 1.000
    gridDesci(10) = 1.000
    gridDesci(20) = 64.0
    
    isccp_tskin_obs(i)%mi = 360*180
    allocate(isccp_tskin_obs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(isccp_tskin_obs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))
    allocate(isccp_tskin_obs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    
    call neighbor_interp_input(gridDesci,LVT_rc%gridDesc(:),&
         LVT_rc%lnc*LVT_rc%lnr,&
         isccp_tskin_obs(i)%rlat,&
         isccp_tskin_obs(i)%rlon,&
         isccp_tskin_obs(i)%n11)

  end subroutine ISCCP_TskinobsInit


end module ISCCP_TskinobsMod
