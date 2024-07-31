!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: WindSatsm_obsMod
! 
! !DESCRIPTION: 
!  This module handles the observation plugin for the 
!  WindSat soil moisture retrievals.  This plugin processes
!  the retrievals based on 10, 18.7 and 37GHz channels. For more
!  details, please see: 
! 
!   Li et al, "WindSat global soil moisture retrieval and validation",
!   IEEE Transactions on Geoscience and Remote Sensing, 2009
!   
! !REVISION HISTORY: 
!  28 Jan 2013: Sujay Kumar, Initial Specification
!
module WindSatsm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: WindSatsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: WindSatsmobs
!EOP
  type, public :: windsatsmobsdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     integer                :: mi
     real,    allocatable   :: smobs(:,:)
     real,    allocatable   :: tmobs(:,:)
     logical                :: startmode 
     integer                :: nc, nr
     type(proj_info)        :: windsatproj
     integer, allocatable   :: n11(:)
  end type windsatsmobsdec

  type(windsatsmobsdec), allocatable:: WindSatsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: WindSatsm_obsInit
! \label{WindSatsm_obsInit}
! 
! !INTERFACE: 
  subroutine WindSatsm_obsinit()
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
!  for reading WindSat soil moisture data. 
! 
!EOP
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(WindSatsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'WindSat soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, WindSatsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'WindSat soil moisture observation directory: not defined')
    enddo

    do n=1,LDT_rc%nnest
       WindSatsmobs(n)%startmode = .true. 

       allocate(WindSatsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))
       allocate(WindSatsmobs(n)%tmobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       WindSatsmobs(n)%smobs = -9999.0
       WindSatsmobs(n)%tmobs = -9999.0

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
           
       gridDesci = 0 
       gridDesci(1) = 9 
       gridDesci(2) = 1383
       gridDesci(3) = 586
       gridDesci(4) = -90.0
       gridDesci(5) = -179.6096
       gridDesci(7) = 83.33788
       gridDesci(8) = 180.1301
       gridDesci(9) = 1 !Ml
    
       WindSatsmobs(n)%nc = 1383
       WindSatsmobs(n)%nr = 586
       WindSatsmobs(n)%mi = WindSatsmobs(n)%nc*&
            WindSatsmobs(n)%nr

       
       allocate(WindSatsmobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
       
       call neighbor_interp_input(n, gridDesci, WindSatsmobs(n)%n11)

    enddo
  end subroutine WindSatsm_obsinit
     

end module WindSatsm_obsMod
