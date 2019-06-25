!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: NASASMAPsm_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! NASASMAP soil moisture retrievals
!
!
!   
! !REVISION HISTORY: 
!  21 Aug 2016: Sujay Kumar, Initial Specification
!  12 Feb 2018: Mahdi Navari, openwater proximity detection was added
! 			edited to read New version of the SPL3SMP_R14 (file structure
! 			 differs from the previous versions)
!  04 Jun 2019: Sujay Kumar, Updated to support SMAP L2 retrievals 

module NASASMAPsm_obsMod
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: NASASMAPsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: NASASMAPsmobs
!EOP
  type, public :: smapsmobsdec

     character*100          :: odir
     character*20           :: data_designation
     real                   :: search_radius
     integer                :: mo
     real,    allocatable   :: smobs(:,:)
     integer                :: nc, nr
     type(proj_info)        :: proj
     integer, allocatable   :: n11(:)
  end type smapsmobsdec

  type(smapsmobsdec), allocatable:: NASASMAPsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: NASASMAPsm_obsInit
! \label{NASASMAPsm_obsInit}
! 
! !INTERFACE: 
  subroutine NASASMAPsm_obsinit()
! !USES: 
    use ESMF
    use LDT_coreMod,    only : LDT_rc, LDT_config
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod, only : LDT_clock, LDT_calendar
    use LDT_logMod,     only : LDT_verify, LDT_logunit

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading NASASMAP soil moisture data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(NASASMAPsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'NASA SMAP soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, NASASMAPsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'NASA SMAP soil moisture observation directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'NASA SMAP soil moisture data designation:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, &
            NASASMAPsmobs(n)%data_designation, &
            rc=status)
       call LDT_verify(status, &
            'NASA SMAP soil moisture data designation: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'NASA SMAP search radius for openwater proximity detection:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, NASASMAPsmobs(n)%search_radius, &
            rc=status)
       call LDT_verify(status, &
            'NASA SMAP search radius for openwater proximity detection: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(NASASMAPsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       NASASMAPsmobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
       if(NASASMAPsmobs(n)%data_designation.eq."SPL3SMP") then 
          NASASMAPsmobs(n)%nc = 964
          NASASMAPsmobs(n)%nr = 406

          gridDesci = 0 
          gridDesci(1) = 9
          gridDesci(2) = 964
          gridDesci(3) = 406
          gridDesci(9) = 4 !M36 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.36 
          gridDesci(11) = 1 !for the global switch

       elseif(NASASMAPsmobs(n)%data_designation.eq."SPL3SMP_E") then 
          NASASMAPsmobs(n)%nc = 3856
          NASASMAPsmobs(n)%nr = 1624

          gridDesci = 0 
          gridDesci(1) = 9
          gridDesci(2) = 3856
          gridDesci(3) = 1624
          gridDesci(9) = 5 !M09 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.09 
          gridDesci(11) = 1 !for the global switch

       elseif(NASASMAPsmobs(n)%data_designation.eq."SPL2SMP") then 
          NASASMAPsmobs(n)%nc = 964
          NASASMAPsmobs(n)%nr = 406

          gridDesci = 0 
          gridDesci(1) = 9
          gridDesci(2) = 964
          gridDesci(3) = 406
          gridDesci(9) = 4 !M36 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.36 
          gridDesci(11) = 1 !for the global switch

       elseif(NASASMAPsmobs(n)%data_designation.eq."SPL2SMP_E") then 
          NASASMAPsmobs(n)%nc = 3856
          NASASMAPsmobs(n)%nr = 1624

          gridDesci = 0 
          gridDesci(1) = 9
          gridDesci(2) = 3856
          gridDesci(3) = 1624
          gridDesci(9) = 5 !M09 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.09 
          gridDesci(11) = 1 !for the global switch
       endif
       allocate(NASASMAPsmobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       call neighbor_interp_input (n, gridDesci,&
            NASASMAPsmobs(n)%n11)
       

    enddo
  end subroutine NASASMAPsm_obsinit
     
end module NASASMAPsm_obsMod
