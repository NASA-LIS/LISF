!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: LPRMvod_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! LPRM vegetation optical depth (vod) retrievals
!
! !REVISION HISTORY: 
!  28 May 2019: Sujay Kumar, Initial Specification
!
module LPRMvod_obsMod
! !USES: 
  use ESMF
  use map_utils

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: LPRMvod_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: LPRMvodobs
!EOP
  type, public :: lprmvodobsdec

     character*100          :: odir
     character*20           :: data_designation
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
  end type lprmvodobsdec

  type(lprmvodobsdec), allocatable:: LPRMvodobs(:)

contains
  
!BOP
! 
! !ROUTINE: LPRMvod_obsInit
! \label{LPRMvod_obsInit}
! 
! !INTERFACE: 
  subroutine LPRMvod_obsinit()
! !USES: 
    use ESMF
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

    allocate(LPRMvodobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'LPRM vegetation optical depth observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, LPRMvodobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'LPRM vegetation optical depth observation directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'LPRM vegetation optical depth data designation:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, &
            LPRMvodobs(n)%data_designation, &
            rc=status)
       if(status.ne.0) then 
          write(LDT_logunit,*) "[ERR] 'LPRM vegetation optical depth data designation:' not defined"
          write(LDT_logunit,*) "[ERR] supported options are 'C-band' or 'X-band'"
          call LDT_endrun()
       endif
    enddo

    do n=1,LDT_rc%nnest

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%vod_obs, &
            "-",1,1)
       LDT_DAobsData(n)%vod_obs%selectStats = 1
    
       LPRMvodobs(n)%nc = 1440
       LPRMvodobs(n)%nr =  720

       gridDesci(1) = 0
       gridDesci(2) = 1440
       gridDesci(3) = 720
       gridDesci(4) = -89.875
       gridDesci(5) = -179.875
       gridDesci(6) = 128 
       gridDesci(7) = 89.875
       gridDesci(8) = 179.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64
           
       allocate(LPRMvodobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(LPRMvodobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(LPRMvodobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(LPRMvodobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       

       allocate(LPRMvodobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(LPRMvodobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(LPRMvodobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       allocate(LPRMvodobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       

       call bilinear_interp_input (n, gridDesci,&
            LPRMvodobs(n)%n11,&
            LPRMvodobs(n)%n12,& 
            LPRMvodobs(n)%n21,& 
            LPRMvodobs(n)%n22,& 
            LPRMvodobs(n)%w11,&
            LPRMvodobs(n)%w12,& 
            LPRMvodobs(n)%w21,& 
            LPRMvodobs(n)%w22) 
    enddo
  end subroutine LPRMvod_obsinit
     
end module LPRMvod_obsMod
