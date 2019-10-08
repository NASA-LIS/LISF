!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: NASASMAPvod_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! NASASMAP vegetation optical depth (vod) retrievals
!
! !REVISION HISTORY: 
!  26 Mar 2019: Sujay Kumar, Initial Specification
!
module NASASMAPvod_obsMod
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: NASASMAPvod_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: NASASMAPvodobs
!EOP
  type, public :: smapvodobsdec

     character*100          :: odir
     character*20           :: data_designation
     integer                :: mo
     integer                :: nc, nr
     type(proj_info)        :: proj
     integer, allocatable   :: n11(:)
  end type smapvodobsdec

  type(smapvodobsdec), allocatable:: NASASMAPvodobs(:)

contains
  
!BOP
! 
! !ROUTINE: NASASMAPvod_obsInit
! \label{NASASMAPvod_obsInit}
! 
! !INTERFACE: 
  subroutine NASASMAPvod_obsinit()
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
!  for reading NASA SMAP vegetation optical depth data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(NASASMAPvodobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'NASA SMAP vegetation optical depth observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, NASASMAPvodobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'NASA SMAP vegetation optical depth observation directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'NASA SMAP vegetation optical depth data designation:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, &
            NASASMAPvodobs(n)%data_designation, &
            rc=status)
       call LDT_verify(status, &
            'NASA SMAP vegetation optical depth data designation: not defined')
    enddo

    do n=1,LDT_rc%nnest

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%vod_obs, &
            "-",1,1)
       LDT_DAobsData(n)%vod_obs%selectStats = 1
    
       if(NASASMAPvodobs(n)%data_designation.eq."SPL3SMP") then 
          NASASMAPvodobs(n)%nc = 964
          NASASMAPvodobs(n)%nr = 406

          gridDesci = 0 
          gridDesci(1) = 9
          gridDesci(2) = 964
          gridDesci(3) = 406
          gridDesci(9) = 4 !M36 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.36 
          gridDesci(11) = 1 !for the global switch
           
       elseif(NASASMAPvodobs(n)%data_designation.eq."SPL3SMP_E") then 
          NASASMAPvodobs(n)%nc = 3856
          NASASMAPvodobs(n)%nr = 1624

          gridDesci = 0 
          gridDesci(1) = 9
          gridDesci(2) = 3856
          gridDesci(3) = 1624
          gridDesci(9) = 5 !M09 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.09 
          gridDesci(11) = 1 !for the global switch
           
       elseif(NASASMAPvodobs(n)%data_designation.eq."SPL2SMP") then 
          NASASMAPvodobs(n)%nc = 964
          NASASMAPvodobs(n)%nr = 406

          gridDesci = 0 
          gridDesci(1) = 9
          gridDesci(2) = 964
          gridDesci(3) = 406
          gridDesci(9) = 4 !M36 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.36 
          gridDesci(11) = 1 !for the global switch

       elseif(NASASMAPvodobs(n)%data_designation.eq."SPL2SMP_E") then 
          NASASMAPvodobs(n)%nc = 3856
          NASASMAPvodobs(n)%nr = 1624

          gridDesci = 0 
          gridDesci(1) = 9
          gridDesci(2) = 3856
          gridDesci(3) = 1624
          gridDesci(9) = 5 !M09 grid
          gridDesci(20) = 64
          gridDesci(10) = 0.09 
          gridDesci(11) = 1 !for the global switch

       endif
       
       allocate(NASASMAPvodobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       call neighbor_interp_input (n, gridDesci,&
            NASASMAPvodobs(n)%n11)

    enddo
  end subroutine NASASMAPvod_obsinit
     
end module NASASMAPvod_obsMod
