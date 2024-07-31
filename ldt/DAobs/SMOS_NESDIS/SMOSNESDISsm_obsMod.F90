!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: SMOSNESDISsm_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! SMOS soil moisture retrievals from NOAA NESDIS. 
!
!
!   
! !REVISION HISTORY: 
!  21 Aug 2016: Sujay Kumar, Initial Specification
!
module SMOSNESDISsm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSNESDISsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMOSNESDISsmobs
!EOP
  type, public :: smosnesdissmobsdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     integer                :: mo
     real,    allocatable   :: smobs(:,:)
     integer                :: nc, nr
     type(proj_info)        :: proj
     integer, allocatable   :: n11(:)
  end type smosnesdissmobsdec

  type(smosnesdissmobsdec), allocatable:: SMOSNESDISsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: SMOSNESDISsm_obsInit
! \label{SMOSNESDISsm_obsInit}
! 
! !INTERFACE: 
  subroutine SMOSNESDISsm_obsinit()
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
!  for reading SMOSNESDIS soil moisture data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(SMOSNESDISsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMOS NESDIS soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMOSNESDISsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'SMOS NESDIS soil moisture observation directory: not defined')
    enddo


    do n=1,LDT_rc%nnest

       allocate(SMOSNESDISsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       SMOSNESDISsmobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
       SMOSNESDISsmobs(n)%nc = 1440
       SMOSNESDISsmobs(n)%nr = 720

       call map_set(PROJ_LATLON, -89.875,-179.875,&
            0.0, 0.25,0.25, 0.0,&
            SMOSNESDISsmobs(n)%nc,SMOSNESDISsmobs(n)%nr,&
            SMOSNESDISsmobs(n)%proj)
       
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
           
       allocate(SMOSNESDISsmobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       call neighbor_interp_input (n, gridDesci,&
            SMOSNESDISsmobs(n)%n11)
    enddo
  end subroutine SMOSNESDISsm_obsinit
     
end module SMOSNESDISsm_obsMod
