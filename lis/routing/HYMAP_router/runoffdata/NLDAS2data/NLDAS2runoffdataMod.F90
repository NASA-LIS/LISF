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
module NLDAS2runoffdataMod
!BOP
! 
! !MODULE: NLDAS2runoffdataMod
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
  public :: NLDAS2runoffdata_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: NLDAS2runoffdata_struc
  
  type, public :: NLDAS2runoffdatadec
     
     real                    :: outInterval 
     character(len=LIS_CONST_PATH_LEN) :: odir 
     character*20            :: model_name
     integer                 :: nc, nr
     integer, allocatable    :: n11(:)
  end type NLDAS2runoffdatadec

  type(NLDAS2runoffdatadec), allocatable :: NLDAS2runoffdata_struc(:)

contains
 
!BOP
!
! !ROUTINE: NLDAS2runoffdata_init
! \label{NLDAS2runoffdata_init}
! 
  subroutine NLDAS2runoffdata_init
    !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod

    integer              :: n 
    integer              :: status
    character*10         :: time
    real                 :: gridDesc(50)

    allocate(NLDAS2runoffdata_struc(LIS_rc%nnest))
       
    call ESMF_ConfigFindLabel(LIS_config,&
         "NLDAS2 runoff data output directory:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,NLDAS2runoffdata_struc(n)%odir,rc=status)
       call LIS_verify(status,&
            "NLDAS2 runoff data output directory: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "NLDAS2 runoff data model name:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            NLDAS2runoffdata_struc(n)%model_name,rc=status)
       call LIS_verify(status,&
            "NLDAS2 runoff data model name: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "NLDAS2 runoff data output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "NLDAS2 runoff data output interval: not defined")

       call LIS_parseTimeString(time,NLDAS2runoffdata_struc(n)%outInterval)
    
       gridDesc = 0.0

       
       NLDAS2runoffdata_struc(n)%nc = 464
       NLDAS2runoffdata_struc(n)%nr = 224
       
       gridDesc(1) = 0  
       gridDesc(2) = NLDAS2runoffdata_struc(n)%nc
       gridDesc(3) = NLDAS2runoffdata_struc(n)%nr
       gridDesc(4) = 25.0625
       gridDesc(5) = -124.9375
       gridDesc(7) = 52.9375
       gridDesc(8) = -67.0625
       gridDesc(6) = 128
       gridDesc(9) = 0.125
       gridDesc(10) = 0.125
       gridDesc(20) = 64
       
       if(LIS_isAtAfinerResolution(n,0.125)) then
          
          allocate(NLDAS2runoffdata_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call neighbor_interp_input(n,gridDesc, &
               NLDAS2runoffdata_struc(n)%n11)
       else
          allocate(NLDAS2runoffdata_struc(n)%n11(&
               NLDAS2runoffdata_struc(n)%nc*NLDAS2runoffdata_struc(n)%nr))
          call upscaleByAveraging_input(gridDesc,&
               LIS_rc%gridDesc(n,:),&
               NLDAS2runoffdata_struc(n)%nc*NLDAS2runoffdata_struc(n)%nr,&
               LIS_rc%lnc(n)*LIS_rc%lnr(n),NLDAS2runoffdata_struc(n)%n11)
       endif
    enddo
    
  end subroutine NLDAS2runoffdata_init
end module NLDAS2runoffdataMod
