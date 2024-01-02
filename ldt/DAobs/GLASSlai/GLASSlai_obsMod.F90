!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: GLASSlai_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! GLASS LAI retrievals
!
! !REVISION HISTORY: 
!  26 Mar 2019: Sujay Kumar, Initial Specification
!
module GLASSlai_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GLASSlai_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GLASSlaiobs
!EOP
  type, public :: glasslaiobsdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     character*100          :: source
     real                   :: gridDesci(50)
     integer                :: mo
     real,    allocatable   :: laiobs(:,:)
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
  end type glasslaiobsdec

  type(glasslaiobsdec), allocatable:: GLASSlaiobs(:)

contains
  
!BOP
! 
! !ROUTINE: GLASSlai_obsInit
! \label{GLASSlai_obsInit}
! 
! !INTERFACE: 
  subroutine GLASSlai_obsinit()
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
!  for reading NASA SMAP vegetation optical depth data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: cornerlat1,cornerlon1
    real                    :: cornerlat2,cornerlon2
    integer                 :: n 

    allocate(GLASSlaiobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'GLASS LAI data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, GLASSlaiobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'GLASS LAI data directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'GLASS LAI data source:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, GLASSlaiobs(n)%source, &
            rc=status)
       call LDT_verify(status, &
            'GLASS LAI source: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(GLASSlaiobs(n)%laiobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       GLASSlaiobs(n)%laiobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%lai_obs, &
            "-",1,1)
       LDT_DAobsData(n)%lai_obs%selectStats = 1
    

       cornerlat1 = max(-59.975, &
            nint((LDT_rc%gridDesc(n,4)+59.975)/0.05)*0.05-59.975-2*0.05)
       cornerlon1 = max(-179.975, &
            nint((LDT_rc%gridDesc(n,5)+179.975)/0.05)*0.05-179.975-2*0.05)
       cornerlat2 = min(89.975, &
            nint((LDT_rc%gridDesc(n,7)+59.975)/0.05)*0.05-59.975+2*0.05)
       cornerlon2 = min(179.975, &
            nint((LDT_rc%gridDesc(n,8)+179.975)/0.05)*0.05-179.975+2*0.05)

       GLASSlaiobs(n)%nr = nint((cornerlat2 - cornerlat1)/0.05)+1
       GLASSlaiobs(n)%nc = nint((cornerlon2 - cornerlon1)/0.05)+1

       GLASSlaiobs(n)%gridDesci = 0 
       GLASSlaiobs(n)%gridDesci(1) = 0
       GLASSlaiobs(n)%gridDesci(2) = GLASSlaiobs(n)%nc
       GLASSlaiobs(n)%gridDesci(3) = GLASSlaiobs(n)%nr
       GLASSlaiobs(n)%gridDesci(4) = cornerlat1
       GLASSlaiobs(n)%gridDesci(5) = cornerlon1
       GLASSlaiobs(n)%gridDesci(6) = 128
       GLASSlaiobs(n)%gridDesci(7) = cornerlat2
       GLASSlaiobs(n)%gridDesci(8) = cornerlon2
       GLASSlaiobs(n)%gridDesci(9) = 0.05
       GLASSlaiobs(n)%gridDesci(10) = 0.05 
       GLASSlaiobs(n)%gridDesci(20) = 64
       
       if(LDT_isLDTatAfinerResolution(n,0.05)) then

          allocate(GLASSlaiobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
          allocate(GLASSlaiobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
          allocate(GLASSlaiobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
          allocate(GLASSlaiobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
          allocate(GLASSlaiobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
          allocate(GLASSlaiobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
          allocate(GLASSlaiobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
          allocate(GLASSlaiobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       

          call bilinear_interp_input (n, &
               GLASSlaiobs(n)%gridDesci,&
               GLASSlaiobs(n)%n11,&
               GLASSlaiobs(n)%n12,&
               GLASSlaiobs(n)%n21,&
               GLASSlaiobs(n)%n22,&
               GLASSlaiobs(n)%w11,&
               GLASSlaiobs(n)%w12,&
               GLASSlaiobs(n)%w21,&
               GLASSlaiobs(n)%w22)

       else
          
          allocate(GLASSlaiobs(n)%n11(GLASSlaiobs(n)%nc*&
               GLASSlaiobs(n)%nr))

          call upscaleByAveraging_input (&
               GLASSlaiobs(n)%gridDesci,&
               LDT_rc%gridDesc(n,:),&
               GLASSlaiobs(n)%nc*&
               GLASSlaiobs(n)%nr,&
               LDT_rc%lnc(n)*LDT_rc%lnr(n),&
               GLASSlaiobs(n)%n11)

       endif
    enddo
  end subroutine GLASSlai_obsinit
     
end module GLASSlai_obsMod
