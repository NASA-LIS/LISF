!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: LPRM_AMSREsm_obsMod
! 
! !DESCRIPTION: 
!  This module handles the observation plugin for the 
!  Land Parameter Retrieval Model (LPRM) AMSR-E soil moisture
!  retrievals
! 
!   Ref: Owe et al. 2008; Multi-sensor historical climatology of satellite-
!   derived global land surface moisture. 
!   
! !REVISION HISTORY: 
!  21 July 2010: Sujay Kumar, Initial Specification
!
module LPRM_AMSREsm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
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

     character(len=LDT_CONST_PATH_LEN)          :: odir
     integer                :: mo
     integer                :: rawdata
     real,    allocatable       :: smobs(:,:)
     integer                :: lprmnc, lprmnr
     type(proj_info)        :: lprmproj
     integer, allocatable       :: n11(:)
  end type lprmamsresmobsdec

  type(lprmamsresmobsdec), allocatable:: LPRM_AMSREsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: LPRM_AMSREsm_obsInit
! \label{LPRM_AMSREsm_obsInit}
! 
! !INTERFACE: 
  subroutine LPRM_AMSREsm_obsinit()
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_config
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod, only : LDT_clock, LDT_calendar
    use LDT_logMod,     only : LDT_verify, LDT_logunit

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading LPRM AMSRE soil moisture data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(LPRM_AMSREsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'AMSR-E(LPRM) soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, LPRM_AMSREsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'AMSR-E(LPRM) soil moisture observation directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'AMSR-E(LPRM) use raw data:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, LPRM_AMSREsmobs(n)%rawdata, &
            rc=status)
       call LDT_verify(status, &
            'AMSR-E(LPRM) use raw data: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(LPRM_AMSREsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       LPRM_AMSREsmobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
       LPRM_AMSREsmobs(n)%lprmnc = 1440
       LPRM_AMSREsmobs(n)%lprmnr = 720

       call map_set(PROJ_LATLON, -89.875,-179.875,&
            0.0, 0.25,0.25, 0.0,&
            LPRM_AMSREsmobs(n)%lprmnc,LPRM_AMSREsmobs(n)%lprmnr,&
            LPRM_AMSREsmobs(n)%lprmproj)
       
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
    
       
       allocate(LPRM_AMSREsmobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call neighbor_interp_input(n, gridDesci, &
            LPRM_AMSREsmobs(n)%n11)
    enddo
  end subroutine LPRM_AMSREsm_obsinit
     

end module LPRM_AMSREsm_obsMod
