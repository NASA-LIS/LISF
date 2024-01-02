!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module ERAILandrunoffdataMod
!BOP
! 
! !MODULE: ERAILandrunoffdataMod
! 
! !DESCRIPTION: 
!  This module contains the data structures and routines to handle 
!  runoff data from GLDAS 1.0. This implementation handles both
!  1 degree and 0.25 degree products from different LSMs.
!
! !REVISION HISTORY: 
! 8 Jan 2016: Sujay Kumar, initial implementation
! 
! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: ERAILandrunoffdata_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: ERAILandrunoffdata_struc
  
  type, public :: ERAILandrunoffdatadec
     
     real                    :: outInterval 
     character(len=LIS_CONST_PATH_LEN) :: odir 
     integer                 :: nc, nr
     integer, allocatable    :: n11(:)
     type(ESMF_Time)         :: startTime
  end type ERAILandrunoffdatadec

  type(ERAILandrunoffdatadec), allocatable :: ERAILandrunoffdata_struc(:)

contains
 
!BOP
!
! !ROUTINE: ERAILandrunoffdata_init
! \label{ERAILandrunoffdata_init}
! 
  subroutine ERAILandrunoffdata_init
    !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod

    integer              :: n 
    integer              :: status
    character*10         :: time
    real                 :: gridDesc(50)

    allocate(ERAILandrunoffdata_struc(LIS_rc%nnest))
       
    call ESMF_ConfigFindLabel(LIS_config,&
         "ERA interim land runoff data output directory:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,ERAILandrunoffdata_struc(n)%odir,rc=status)
       call LIS_verify(status,&
            "ERA interim land runoff data output directory: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "ERA interim land runoff data output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "ERA interim land runoff data output interval: not defined")

       call LIS_parseTimeString(time,ERAILandrunoffdata_struc(n)%outInterval)
    
       gridDesc = 0.0

       
       ERAILandrunoffdata_struc(n)%nc = 480
       ERAILandrunoffdata_struc(n)%nr = 241
       
       gridDesc(1) = 0  
       gridDesc(2) = ERAILandrunoffdata_struc(n)%nc
       gridDesc(3) = ERAILandrunoffdata_struc(n)%nr
       gridDesc(4) = -90.000
       gridDesc(5) = -180.000
       gridDesc(7) = 90.000
       gridDesc(8) = 180.000
       gridDesc(6) = 128
       gridDesc(9) = 0.75
       gridDesc(10) = 0.75
       gridDesc(20) = 0
       
       if(LIS_isAtAfinerResolution(n,0.75)) then
          
          allocate(ERAILandrunoffdata_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call neighbor_interp_input(n,gridDesc, &
               ERAILandrunoffdata_struc(n)%n11)
       else
          allocate(ERAILandrunoffdata_struc(n)%n11(&
               ERAILandrunoffdata_struc(n)%nc*ERAILandrunoffdata_struc(n)%nr))
          call upscaleByAveraging_input(gridDesc,&
               LIS_rc%gridDesc(n,:),&
               ERAILandrunoffdata_struc(n)%nc*ERAILandrunoffdata_struc(n)%nr,&
               LIS_rc%lnc(n)*LIS_rc%lnr(n),ERAILandrunoffdata_struc(n)%n11)
       endif

       call ESMF_TimeSet(ERAILandrunoffdata_struc(n)%startTime, yy=1900, &
            mm = 1, dd = 1, h = 0 , m = 0, calendar=LIS_calendar, &
            rc=status)
       
    enddo
    
  end subroutine ERAILandrunoffdata_init
end module ERAILandrunoffdataMod
