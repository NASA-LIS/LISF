!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: VIIRSGVFobsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! VIIRS GVF retrievals
!
! !REVISION HISTORY: 
!  07 Oct 2021: Yonghwan Kwon, Initial Specification
!
module VIIRSGVFobsMod
! !USES:
  use ESMF
  use map_utils

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: VIIRSGVFobsinit !Initializes structures for reading VIIRS GVF data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: VIIRSgvfobs !Object to hold VIIRSgvf observation attributes
!EOP
  type, public :: viirsgvfdec

     character*100        :: odir
     integer              :: nc, nr
     real                 :: gridDesci(50)
     real,    allocatable :: gvfobs(:,:)
     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)

  end type viirsgvfdec

  type(viirsgvfdec),allocatable :: VIIRSgvfobs(:)

contains

!BOP
! 
! !ROUTINE: VIIRSGVFobsinit
! \label{VIIRSGVFobsinit}
!
! !INTERFACE: 
  subroutine VIIRSGVFobsinit()
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
!  for reading the VIIRS GVF data.
! 
!EOP

    !integer                 :: npts
    !type(ESMF_TimeInterval) :: alarmInterval
    !type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    integer                 :: n

    allocate(VIIRSgvfobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'VIIRS GVF data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, VIIRSgvfobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'VIIRS GVF data directory: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(VIIRSgvfobs(n)%gvfobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       VIIRSgvfobs(n)%gvfobs = -9999.0

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%gvf_obs, &
            "-",1,1)
       LDT_DAobsData(n)%gvf_obs%selectStats = 1

       VIIRSgvfobs(n)%nc = 10000
       VIIRSgvfobs(n)%nr = 5000 

       VIIRSgvfobs(n)%gridDesci(1) = 0
       VIIRSgvfobs(n)%gridDesci(2) = VIIRSgvfobs(n)%nc
       VIIRSgvfobs(n)%gridDesci(3) = VIIRSgvfobs(n)%nr
       VIIRSgvfobs(n)%gridDesci(4) = -89.982
       VIIRSgvfobs(n)%gridDesci(5) = -179.982
       VIIRSgvfobs(n)%gridDesci(6) = 128
       VIIRSgvfobs(n)%gridDesci(7) = 89.982
       VIIRSgvfobs(n)%gridDesci(8) = 179.982
       VIIRSgvfobs(n)%gridDesci(9) = 0.036
       VIIRSgvfobs(n)%gridDesci(10) = 0.036
       VIIRSgvfobs(n)%gridDesci(20) = 64

       if(LDT_isLDTatAfinerResolution(n,0.036)) then
    
          allocate(VIIRSgvfobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(VIIRSgvfobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(VIIRSgvfobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(VIIRSgvfobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(VIIRSgvfobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(VIIRSgvfobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(VIIRSgvfobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(VIIRSgvfobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input (n, &
               VIIRSgvfobs(n)%gridDesci,&
               VIIRSgvfobs(n)%n11,&
               VIIRSgvfobs(n)%n12,&
               VIIRSgvfobs(n)%n21,&
               VIIRSgvfobs(n)%n22,&
               VIIRSgvfobs(n)%w11,&
               VIIRSgvfobs(n)%w12,&
               VIIRSgvfobs(n)%w21,&
               VIIRSgvfobs(n)%w22)

       else

          allocate(VIIRSgvfobs(n)%n11(VIIRSgvfobs(n)%nc*&
               VIIRSgvfobs(n)%nr))

          call upscaleByAveraging_input (&
               VIIRSgvfobs(n)%gridDesci,&
               LDT_rc%gridDesc(n,:),&
               VIIRSgvfobs(n)%nc*&
               VIIRSgvfobs(n)%nr,&
               LDT_rc%lnc(n)*LDT_rc%lnr(n),&
               VIIRSgvfobs(n)%n11)

       endif
    enddo   
  end subroutine VIIRSGVFobsinit

end module VIIRSGVFobsMod
