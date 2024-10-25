!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: SMAPEOPLSMobsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! SMAP_E_OPL soil moisture retrievals
!
! !REVISION HISTORY: 
!  06 Jun 2022: Yonghwan Kwon, Initial Specification
!
module SMAPEOPLSMobsMod
! !USES:
  use ESMF
  use map_utils

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAPEOPLSMobsinit  !Initializes structures for reading SMAP_E_OPL SM data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SMAPEOPLsmobs !Object to hold SMAPEOPLsm observation attributes
!EOP
  type, public :: smapeoplsmdec

     character*100        :: odir
     integer              :: nc, nr
     real                 :: gridDesci(50)
     real,    allocatable :: smobs(:,:)
     real                 :: search_radius
     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)

  end type smapeoplsmdec

  type(smapeoplsmdec),allocatable :: SMAPEOPLsmobs(:)

contains

!BOP
! 
! !ROUTINE: SMAPEOPLSMobsinit
! \label{SMAPEOPLSMobsinit}
!
! !INTERFACE: 
  subroutine SMAPEOPLSMobsinit()
! 
! !USES:
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading the SMAP_E_OPL soil moisture data.
! 
!EOP

    integer                 :: status, rc
    integer                 :: n

    allocate(SMAPEOPLsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMAP_E_OPL soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMAPEOPLsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'SMAP_E_OPL soil moisture observation directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'SMAP_E_OPL search radius for openwater proximity detection:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, SMAPEOPLsmobs(n)%search_radius, &
            rc=status)
       call LDT_verify(status, &
            'SMAP_E_OPL search radius for openwater proximity detection: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(SMAPEOPLsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       SMAPEOPLsmobs(n)%smobs = -9999.0

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1

       SMAPEOPLsmobs(n)%nc = 2560
       SMAPEOPLsmobs(n)%nr = 1920

       SMAPEOPLsmobs(n)%gridDesci(1) = 0
       SMAPEOPLsmobs(n)%gridDesci(2) = SMAPEOPLsmobs(n)%nc
       SMAPEOPLsmobs(n)%gridDesci(3) = SMAPEOPLsmobs(n)%nr
       SMAPEOPLsmobs(n)%gridDesci(4) = -89.9531250
       SMAPEOPLsmobs(n)%gridDesci(5) = -179.9296875
       SMAPEOPLsmobs(n)%gridDesci(6) = 128
       SMAPEOPLsmobs(n)%gridDesci(7) = 89.9531250
       SMAPEOPLsmobs(n)%gridDesci(8) = 179.9296875
       SMAPEOPLsmobs(n)%gridDesci(9) = 0.1406250  !dlon
       SMAPEOPLsmobs(n)%gridDesci(10) = 0.0937500 !dlat
       SMAPEOPLsmobs(n)%gridDesci(20) = 64

       if(LDT_isLDTatAfinerResolution(n,0.0937500)) then

          allocate(SMAPEOPLsmobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(SMAPEOPLsmobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(SMAPEOPLsmobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(SMAPEOPLsmobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(SMAPEOPLsmobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(SMAPEOPLsmobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(SMAPEOPLsmobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(SMAPEOPLsmobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input (n, &
               SMAPEOPLsmobs(n)%gridDesci,&
               SMAPEOPLsmobs(n)%n11,&
               SMAPEOPLsmobs(n)%n12,&
               SMAPEOPLsmobs(n)%n21,&
               SMAPEOPLsmobs(n)%n22,&
               SMAPEOPLsmobs(n)%w11,&
               SMAPEOPLsmobs(n)%w12,&
               SMAPEOPLsmobs(n)%w21,&
               SMAPEOPLsmobs(n)%w22)

       else

          allocate(SMAPEOPLsmobs(n)%n11(SMAPEOPLsmobs(n)%nc*&
               SMAPEOPLsmobs(n)%nr))

          call upscaleByAveraging_input (&
               SMAPEOPLsmobs(n)%gridDesci,&
               LDT_rc%gridDesc(n,:),&
               SMAPEOPLsmobs(n)%nc*&
               SMAPEOPLsmobs(n)%nr,&
               LDT_rc%lnc(n)*LDT_rc%lnr(n),&
               SMAPEOPLsmobs(n)%n11)

       endif
    enddo

   end subroutine SMAPEOPLSMobsinit

end module SMAPEOPLSMobsMod
