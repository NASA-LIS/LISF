!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: GCOMWAMSR2TB_ANNdataMod
! 
! !DESCRIPTION: 
!
!   
! !REVISION HISTORY: 
!  20 Jan 2018: Sujay Kumar, Initial Specification
!
module GCOMWAMSR2TB_ANNdataMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GCOMWAMSR2TB_ANNdatainit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GCOMWAMSR2TBobs
!EOP
  type, public :: GCOMWAMSR2TBobsdec

     character(len=LDT_CONST_PATH_LEN) :: odir
     integer                :: nc
     integer                :: nr
     real,    allocatable   :: TB_A_H(:,:,:)
     real,    allocatable   :: TB_A_V(:,:,:)
     real,    allocatable   :: TB_D_H(:,:,:)
     real,    allocatable   :: TB_D_V(:,:,:)
     character*100, allocatable          :: tb_a_h_name(:)
     character*100, allocatable          :: tb_a_v_name(:)
     character*100, allocatable          :: tb_d_h_name(:)
     character*100, allocatable          :: tb_d_v_name(:)
     logical                :: startmode 
     integer, allocatable   :: n11(:)
  end type GCOMWAMSR2TBobsdec

  type(GCOMWAMSR2TBobsdec), allocatable:: GCOMWAMSR2TBobs(:)

contains
  
!BOP
! 
! !ROUTINE: GCOMWAMSR2TB_ANNdatainit
! \label{GCOMWAMSR2TB_ANNdatainit}
! 
! !INTERFACE: 
  subroutine GCOMWAMSR2TB_ANNdatainit()
! !USES: 
    use LDT_coreMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading SYNTHETIC soil moisture data. 
! 
!EOP
    integer                 :: npts
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 
    real                    :: ts
    character*3             :: fnest
    character*20            :: stime

    n = 1

    allocate(GCOMWAMSR2TBobs(LDT_rc%nnest))

    call LDT_update_timestep(LDT_rc, n, 86400.0)

    write(fnest,'(i3.3)') n
    call LDT_registerAlarm("GCOMW AMSR2 Tb alarm "//trim(fnest),&
         86400.0, 86400.0)

    call ESMF_ConfigFindLabel(LDT_config, &
         'GCOMW AMSR2 Tb data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, GCOMWAMSR2TBobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'GCOMW AMSR2 Tb data directory: not defined')
    enddo

    do n=1,LDT_rc%nnest
       GCOMWAMSR2TBobs(n)%startmode = .true. 

       GCOMWAMSR2TBobs(n)%nc = 1440
       GCOMWAMSR2TBobs(n)%nr = 720

       gridDesci = 0
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

       allocate(GCOMWAMSR2TBobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       call neighbor_interp_input (n, gridDesci,&
            GCOMWAMSR2TBobs(n)%n11)

       allocate(GCOMWAMSR2TBobs(n)%TB_A_H(LDT_rc%lnc(n),LDT_rc%lnr(n),2))
       allocate(GCOMWAMSR2TBobs(n)%TB_A_V(LDT_rc%lnc(n),LDT_rc%lnr(n),2))

       allocate(GCOMWAMSR2TBobs(n)%TB_D_H(LDT_rc%lnc(n),LDT_rc%lnr(n),2))
       allocate(GCOMWAMSR2TBobs(n)%TB_D_V(LDT_rc%lnc(n),LDT_rc%lnr(n),2))

       GCOMWAMSR2TBobs(n)%TB_A_H = -9999.0
       GCOMWAMSR2TBobs(n)%TB_A_V = -9999.0

       GCOMWAMSR2TBobs(n)%TB_D_H = -9999.0
       GCOMWAMSR2TBobs(n)%TB_D_V = -9999.0
    
       allocate(GCOMWAMSR2TBobs(n)%TB_A_H_name(2))
       allocate(GCOMWAMSR2TBobs(n)%TB_A_V_name(2))
       allocate(GCOMWAMSR2TBobs(n)%TB_D_H_name(2))
       allocate(GCOMWAMSR2TBobs(n)%TB_D_V_name(2))

       GCOMWAMSR2TBobs(n)%TB_A_H_name(1) = 'TB_A_H_18'
       GCOMWAMSR2TBobs(n)%TB_A_H_name(2) = 'TB_A_H_36'

       GCOMWAMSR2TBobs(n)%TB_D_H_name(1) = 'TB_D_H_18'
       GCOMWAMSR2TBobs(n)%TB_D_H_name(2) = 'TB_D_H_36'

       GCOMWAMSR2TBobs(n)%TB_A_V_name(1) = 'TB_A_V_18'
       GCOMWAMSR2TBobs(n)%TB_A_V_name(2) = 'TB_A_V_36'

       GCOMWAMSR2TBobs(n)%TB_D_V_name(1) = 'TB_D_V_18'
       GCOMWAMSR2TBobs(n)%TB_D_V_name(2) = 'TB_D_V_36'

    enddo
  end subroutine GCOMWAMSR2TB_ANNdatainit
     
end module GCOMWAMSR2TB_ANNdataMod
