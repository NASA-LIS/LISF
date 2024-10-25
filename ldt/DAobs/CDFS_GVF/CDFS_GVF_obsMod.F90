!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: CDFSGVFobsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! CDFS GVF retrievals
!
! !REVISION HISTORY: 
!  04 Mar 2022: Yonghwan Kwon, Initial Specification
!
module CDFSGVFobsMod
! !USES:
  use ESMF
  use map_utils

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: CDFSGVFobsinit !Initializes structures for reading CDFS GVF data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: CDFSgvfobs !Object to hold CDFSgvf observation attributes
!EOP
  type, public :: cdfsgvfdec

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

  end type cdfsgvfdec

  type(cdfsgvfdec),allocatable :: CDFSgvfobs(:)

contains

!BOP
! 
! !ROUTINE: CDFSGVFobsinit
! \label{CDFSGVFobsinit}
!
! !INTERFACE: 
  subroutine CDFSGVFobsinit()
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
!  for reading the CDFS GVF data.
! 
!EOP

    !integer                 :: npts
    !type(ESMF_TimeInterval) :: alarmInterval
    !type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    integer                 :: n

    allocate(CDFSgvfobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'CDFS GVF data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, CDFSgvfobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'CDFS GVF data directory: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(CDFSgvfobs(n)%gvfobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       CDFSgvfobs(n)%gvfobs = -9999.0

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%gvf_obs, &
            "-",1,1)
       LDT_DAobsData(n)%gvf_obs%selectStats = 1

       CDFSgvfobs(n)%nc = 7200
       CDFSgvfobs(n)%nr = 3600 

       CDFSgvfobs(n)%gridDesci(1) = 0
       CDFSgvfobs(n)%gridDesci(2) = CDFSgvfobs(n)%nc
       CDFSgvfobs(n)%gridDesci(3) = CDFSgvfobs(n)%nr
       CDFSgvfobs(n)%gridDesci(4) = -89.975
       CDFSgvfobs(n)%gridDesci(5) = -179.975
       CDFSgvfobs(n)%gridDesci(6) = 128
       CDFSgvfobs(n)%gridDesci(7) = 89.975
       CDFSgvfobs(n)%gridDesci(8) = 179.975
       CDFSgvfobs(n)%gridDesci(9) = 0.05
       CDFSgvfobs(n)%gridDesci(10) = 0.05
       CDFSgvfobs(n)%gridDesci(20) = 64

       if(LDT_isLDTatAfinerResolution(n,0.05)) then
    
          allocate(CDFSgvfobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(CDFSgvfobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(CDFSgvfobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(CDFSgvfobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(CDFSgvfobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(CDFSgvfobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(CDFSgvfobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(CDFSgvfobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input (n, &
               CDFSgvfobs(n)%gridDesci,&
               CDFSgvfobs(n)%n11,&
               CDFSgvfobs(n)%n12,&
               CDFSgvfobs(n)%n21,&
               CDFSgvfobs(n)%n22,&
               CDFSgvfobs(n)%w11,&
               CDFSgvfobs(n)%w12,&
               CDFSgvfobs(n)%w21,&
               CDFSgvfobs(n)%w22)

       else

          allocate(CDFSgvfobs(n)%n11(CDFSgvfobs(n)%nc*&
               CDFSgvfobs(n)%nr))

          call upscaleByAveraging_input (&
               CDFSgvfobs(n)%gridDesci,&
               LDT_rc%gridDesc(n,:),&
               CDFSgvfobs(n)%nc*&
               CDFSgvfobs(n)%nr,&
               LDT_rc%lnc(n)*LDT_rc%lnr(n),&
               CDFSgvfobs(n)%n11)

       endif
    enddo   
  end subroutine CDFSGVFobsinit

end module CDFSGVFobsMod
