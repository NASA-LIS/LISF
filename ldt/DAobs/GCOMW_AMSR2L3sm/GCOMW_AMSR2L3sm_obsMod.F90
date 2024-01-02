!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: GCOMW_AMSR2L3sm_obsMod
! 
! !DESCRIPTION: 
!
!  This plugin handles the processing of the L3 soil moisture retrievals from
!  AMSR2 instrument onboard GCOM-W1 satellite. 
!   
! !REVISION HISTORY: 
!  01 Oct 2012: Sujay Kumar, Initial Specification
!
module GCOMW_AMSR2L3sm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GCOMW_AMSR2L3sm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GCOMW_AMSR2L3smobs
!EOP
  type, public :: amsr2smobsdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     integer                :: mo
     real,    allocatable   :: smobs(:,:)
     integer                :: amsr2nc, amsr2nr
     type(proj_info)        :: amsr2proj
     real                   :: datares
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)
     real,  allocatable     :: w11(:)
     real,  allocatable     :: w12(:)
     real,  allocatable     :: w21(:)
     real,  allocatable     :: w22(:)
  end type amsr2smobsdec

  type(amsr2smobsdec), allocatable:: GCOMW_AMSR2L3smobs(:)

contains
  
!BOP
! 
! !ROUTINE: GCOMW_AMSR2L3sm_obsInit
! \label{GCOMW_AMSR2L3sm_obsInit}
! 
! !INTERFACE: 
  subroutine GCOMW_AMSR2L3sm_obsinit()
! !USES: 
    use LDT_coreMod
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod, only : LDT_clock, LDT_calendar
    use LDT_logMod,     only : LDT_verify, LDT_logunit

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading GCOMW_AMSR2 soil moisture data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(GCOMW_AMSR2L3smobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'GCOMW AMSR2 L3 soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, GCOMW_AMSR2L3smobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'GCOMW AMSR2 L3 soil moisture observation directory: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(GCOMW_AMSR2L3smobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       GCOMW_AMSR2L3smobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
       GCOMW_AMSR2L3smobs(n)%amsr2nc = 3600
       GCOMW_AMSR2L3smobs(n)%amsr2nr = 1800

       call map_set(PROJ_LATLON, -89.95,-179.95,&
            0.0, 0.10,0.10, 0.0,&
            GCOMW_AMSR2L3smobs(n)%amsr2nc,GCOMW_AMSR2L3smobs(n)%amsr2nr,&
            GCOMW_AMSR2L3smobs(n)%amsr2proj)
       
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 3600
       gridDesci(3) = 1800
       gridDesci(4) = -89.95
       gridDesci(5) = -179.95
       gridDesci(6) = 128
       gridDesci(7) = 89.95
       gridDesci(8) = 179.95
       gridDesci(9) = 0.10
       gridDesci(10) = 0.10
       gridDesci(20) = 64

       GCOMW_AMSR2L3smobs(n)%datares = 0.10
    
       if(LDT_isLDTatAfinerResolution(n,GCOMW_AMSR2L3smobs(n)%datares)) then 
          
          allocate(GCOMW_AMSR2L3smobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3smobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3smobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3smobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          
          allocate(GCOMW_AMSR2L3smobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3smobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3smobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GCOMW_AMSR2L3smobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          
          call bilinear_interp_input(n, gridDesci, &
               GCOMW_AMSR2L3smobs(n)%n11, &
               GCOMW_AMSR2L3smobs(n)%n12, GCOMW_AMSR2L3smobs(n)%n21, &
               GCOMW_AMSR2L3smobs(n)%n22, GCOMW_AMSR2L3smobs(n)%w11, &
               GCOMW_AMSR2L3smobs(n)%w12, GCOMW_AMSR2L3smobs(n)%w21, &
               GCOMW_AMSR2L3smobs(n)%w22)
       else
          allocate(GCOMW_AMSR2L3smobs(n)%n11(GCOMW_AMSR2L3smobs(n)%amsr2nc*&
               GCOMW_AMSR2L3smobs(n)%amsr2nr))
          call upscaleByAveraging_input(gridDesci,&
               LDT_rc%gridDesc(n,:),&
               GCOMW_AMSR2L3smobs(n)%amsr2nc*GCOMW_AMSR2L3smobs(n)%amsr2nr,&
               LDT_rc%lnc(n)*LDT_rc%lnr(n),&
               GCOMW_AMSR2L3smobs(n)%n11)
       endif

    enddo
  end subroutine GCOMW_AMSR2L3sm_obsinit
     
end module GCOMW_AMSR2L3sm_obsMod
