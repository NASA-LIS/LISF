!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: THySM_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! THySM soil moisture retrievals
!
! !REVISION HISTORY: 
!  29 Mar 2021: Sujay Kumar, Initial Specification
!
module THySM_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: THySM_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: THySMobs
!EOP
  type, public :: thysmobsdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     integer                :: mo
     integer                :: nc, nr
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)

     real, allocatable      :: w11(:)
     real, allocatable      :: w12(:)
     real, allocatable      :: w21(:)
     real, allocatable      :: w22(:)
  end type thysmobsdec

  type(thysmobsdec), allocatable:: THySMobs(:)

contains
  
!BOP
! 
! !ROUTINE: THySM_obsInit
! \label{THySM_obsInit}
! 
! !INTERFACE: 
  subroutine THySM_obsinit()
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
!  for reading the LPRM vegetation optical depth data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(THySMobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'THySM observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, THySMobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'THySM observation directory: not defined')
    enddo

    do n=1,LDT_rc%nnest

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
       THySMobs(n)%nc = 5800
       THySMobs(n)%nr = 2800

       gridDesci(1) = 0
       gridDesci(2) = 5800
       gridDesci(3) = 2800
       gridDesci(4) = 25.005
       gridDesci(5) = -124.995
       gridDesci(6) = 128 
       gridDesci(7) = 52.995
       gridDesci(8) = -67.005
       gridDesci(9) = 0.01
       gridDesci(10) = 0.01
       gridDesci(20) = 64
           
       allocate(THySMobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(THySMobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(THySMobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(THySMobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       

       allocate(THySMobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(THySMobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(THySMobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(THySMobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       

       call bilinear_interp_input (n, gridDesci,&
            THySMobs(n)%n11,&
            THySMobs(n)%n12,& 
            THySMobs(n)%n21,& 
            THySMobs(n)%n22,& 
            THySMobs(n)%w11,&
            THySMobs(n)%w12,& 
            THySMobs(n)%w21,& 
            THySMobs(n)%w22) 
    enddo
  end subroutine THySM_obsinit
     
end module THySM_obsMod
