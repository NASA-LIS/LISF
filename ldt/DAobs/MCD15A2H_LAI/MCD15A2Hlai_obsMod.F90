!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: MCD15A2Hlai_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! MODIS MCD15A2H LAI retrievals
!
! !REVISION HISTORY: 
!  12 Nov 2020: Wanshu Nie, Initial Specification
!
module MCD15A2Hlai_obsMod
! !USES:
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none
 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: MCD15A2Hlai_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: MCD15A2Hlaiobs
!EOP
  type, public :: MCD15A2Hlai_dec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     character*100          :: version
     real                   :: gridDesci(50)
     real,    allocatable   :: laiobs(:,:)
     integer                :: climofill
     integer                :: qcflag
     integer                :: nc, nr
     type(proj_info)        :: proj
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n22(:)
     real,    allocatable   :: w11(:)
     real,    allocatable   :: w12(:)
     real,    allocatable   :: w21(:)
     real,    allocatable   :: w22(:)

  end type MCD15A2Hlai_dec

  type(MCD15A2Hlai_dec),allocatable :: MCD15A2Hlaiobs(:)

contains

!BOP
! 
! !ROUTINE: MCD15A2Hlai_obsInit
! \label{MCD15A2Hlai_obsInit}
! 
! !INTERFACE: 
  subroutine MCD15A2Hlai_obsinit()
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
!  for reading MODIS MCD15A2H LAI data. 
! 
!EOP


    integer                 :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: cornerlat1,cornerlon1
    real                    :: cornerlat2,cornerlon2
    integer                 :: n

    allocate(MCD15A2Hlaiobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'MCD15A2H LAI data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, MCD15A2Hlaiobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'MCD15A2H LAI data directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'MCD15A2H LAI data version:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, MCD15A2Hlaiobs(n)%version, &
            rc=status)
       call LDT_verify(status, &
            'MCD15A2H LAI data version: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'MCD15A2H LAI apply climatological fill values:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, MCD15A2Hlaiobs(n)%climofill, &
            rc=status)
       call LDT_verify(status, &
            'MCD15A2H LAI apply climatological fill values: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'MCD15A2H LAI apply QC flags:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, MCD15A2Hlaiobs(n)%qcflag, &
            rc=status)
       call LDT_verify(status, &
            'MCD15A2H LAI apply QC flags: not defined')
    enddo

   do n=1,LDT_rc%nnest

       allocate(MCD15A2Hlaiobs(n)%laiobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       MCD15A2Hlaiobs(n)%laiobs = -9999.0

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%lai_obs, &
            "-",1,1)
       LDT_DAobsData(n)%lai_obs%selectStats = 1


       cornerlat1 = max(-59.9978927, nint((LDT_rc%gridDesc(n,4)+59.9978927)/0.00416667)*0.00416667-59.9978927-50*0.00416667)
       cornerlon1 = max(-179.9979167, nint((LDT_rc%gridDesc(n,5)+179.9979167)/0.00416667)*0.00416667-179.9979167-50*0.00416667)
       cornerlat2 = min(89.9979167, nint((LDT_rc%gridDesc(n,7)+59.9978927)/0.00416667)*0.00416667-59.9978927+50*0.00416667)
       cornerlon2 = min(179.9979167, nint((LDT_rc%gridDesc(n,8)+179.9979167)/0.00416667)*0.00416667-179.9979167+50*0.00416667)


       MCD15A2Hlaiobs(n)%nc = nint((cornerlon2-cornerlon1)/0.00416667)+1
       MCD15A2Hlaiobs(n)%nr = nint((cornerlat2-cornerlat1)/0.00416667)+1

       MCD15A2Hlaiobs(n)%gridDesci(1) = 0
       MCD15A2Hlaiobs(n)%gridDesci(2) = MCD15A2Hlaiobs(n)%nc
       MCD15A2Hlaiobs(n)%gridDesci(3) = MCD15A2Hlaiobs(n)%nr
       MCD15A2Hlaiobs(n)%gridDesci(4) = cornerlat1
       MCD15A2Hlaiobs(n)%gridDesci(5) = cornerlon1
       MCD15A2Hlaiobs(n)%gridDesci(6) = 128
       MCD15A2Hlaiobs(n)%gridDesci(7) = cornerlat2
       MCD15A2Hlaiobs(n)%gridDesci(8) = cornerlon2
       MCD15A2Hlaiobs(n)%gridDesci(9) = 0.00416667
       MCD15A2Hlaiobs(n)%gridDesci(10) = 0.00416667
       MCD15A2Hlaiobs(n)%gridDesci(20) = 64


       if(LDT_isLDTatAfinerResolution(n,0.00416667)) then

          allocate(MCD15A2Hlaiobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(MCD15A2Hlaiobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(MCD15A2Hlaiobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(MCD15A2Hlaiobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(MCD15A2Hlaiobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(MCD15A2Hlaiobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(MCD15A2Hlaiobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(MCD15A2Hlaiobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input (n, &
               MCD15A2Hlaiobs(n)%gridDesci,&
               MCD15A2Hlaiobs(n)%n11,&
               MCD15A2Hlaiobs(n)%n12,&
               MCD15A2Hlaiobs(n)%n21,&
               MCD15A2Hlaiobs(n)%n22,&
               MCD15A2Hlaiobs(n)%w11,&
               MCD15A2Hlaiobs(n)%w12,&
               MCD15A2Hlaiobs(n)%w21,&
               MCD15A2Hlaiobs(n)%w22)

       else
  
          allocate(MCD15A2Hlaiobs(n)%n11(MCD15A2Hlaiobs(n)%nc*&
               MCD15A2Hlaiobs(n)%nr))

          call upscaleByAveraging_input (&
               MCD15A2Hlaiobs(n)%gridDesci,&
               LDT_rc%gridDesc(n,:),&
               MCD15A2Hlaiobs(n)%nc*&
               MCD15A2Hlaiobs(n)%nr,&
               LDT_rc%lnc(n)*LDT_rc%lnr(n),&
               MCD15A2Hlaiobs(n)%n11)

       endif
    enddo
  end subroutine MCD15A2Hlai_obsinit

end module MCD15A2Hlai_obsMod

