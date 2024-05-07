!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: ASCATTUWsm_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
!  NOAA/NESDIS-based ASCAT L2 Soil Moisture (SM) Files (aka, SMOPS)
!
!
!   
! !REVISION HISTORY: 
!  04 May 2013: Sujay Kumar, Initial Specification
!
module ASCATTUWsm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ASCATTUWsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ASCATTUWsmobs
!EOP
  type, public :: ascattuwsmobsdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     integer                :: mo
     real,    allocatable   :: smobs(:,:)
     logical                :: startmode 
     integer                :: nc, nr
     type(proj_info)        :: ascattuwproj
     integer, allocatable   :: n11(:)
  end type ascattuwsmobsdec

  type(ascattuwsmobsdec), allocatable:: ASCATTUWsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: ASCATTUWsm_obsInit
! \label{ASCATTUWsm_obsInit}
! 
! !INTERFACE: 
  subroutine ASCATTUWsm_obsinit()
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
!  for reading RT SMOPS soil moisture data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(ASCATTUWsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'ASCAT (TUW) soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, ASCATTUWsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'ASCAT (TUW) soil moisture observation directory: not defined')
    enddo

    do n=1,LDT_rc%nnest
       ASCATTUWsmobs(n)%startmode = .true. 

       allocate(ASCATTUWsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       ASCATTUWsmobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
       ASCATTUWsmobs(n)%nc = 1440
       ASCATTUWsmobs(n)%nr = 600

       call map_set(PROJ_LATLON, -59.875,-179.875,&
            0.0, 0.25,0.25, 0.0,&
            ASCATTUWsmobs(n)%nc,ASCATTUWsmobs(n)%nr,&
            ASCATTUWsmobs(n)%ascattuwproj)
       
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 1440
       gridDesci(3) = 720
       gridDesci(4) = -59.875
       gridDesci(5) = -179.875
       gridDesci(6) = 128
       gridDesci(7) = 89.875
       gridDesci(8) = 179.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64
    
       
       allocate(ASCATTUWsmobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call neighbor_interp_input(n, gridDesci, ASCATTUWsmobs(n)%n11)

    enddo
  end subroutine ASCATTUWsm_obsinit
     
end module ASCATTUWsm_obsMod
