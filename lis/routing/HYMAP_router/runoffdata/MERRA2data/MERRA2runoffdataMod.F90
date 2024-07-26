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
module MERRA2runoffdataMod
!BOP
! 
! !MODULE: MERRA2runoffdataMod
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
  public :: MERRA2runoffdata_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: MERRA2runoffdata_struc
  
  type, public :: MERRA2runoffdatadec
     
     real                    :: outInterval 
     character(len=LIS_CONST_PATH_LEN) :: odir 
     integer                 :: nc, nr
     integer, allocatable    :: n11(:)
  end type MERRA2runoffdatadec

  type(MERRA2runoffdatadec), allocatable :: MERRA2runoffdata_struc(:)

contains
 
!BOP
!
! !ROUTINE: MERRA2runoffdata_init
! \label{MERRA2runoffdata_init}
! 
  subroutine MERRA2runoffdata_init
    !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod

    integer              :: n 
    integer              :: status
    character*10         :: time
    real                 :: gridDesc(50)

    allocate(MERRA2runoffdata_struc(LIS_rc%nnest))
       
    call ESMF_ConfigFindLabel(LIS_config,&
         "MERRA2 runoff data output directory:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,MERRA2runoffdata_struc(n)%odir,rc=status)
       call LIS_verify(status,&
            "MERRA2 runoff data output directory: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "MERRA2 runoff data output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "MERRA2 runoff data output interval: not defined")

       call LIS_parseTimeString(time,MERRA2runoffdata_struc(n)%outInterval)
    
       gridDesc = 0.0

       
       MERRA2runoffdata_struc(n)%nc = 576
       MERRA2runoffdata_struc(n)%nr = 361
       
       gridDesc(1) = 0  
       gridDesc(2) = MERRA2runoffdata_struc(n)%nc
       gridDesc(3) = MERRA2runoffdata_struc(n)%nr
       gridDesc(4) = -90.000
       gridDesc(5) = -180.000
       gridDesc(7) = 90.000
       gridDesc(8) = 179.375
       gridDesc(6) = 128
       gridDesc(9) = 0.625
       gridDesc(10) = 0.5
       gridDesc(20) = 0
       
       if(LIS_isAtAfinerResolution(n,0.5)) then
          
          allocate(MERRA2runoffdata_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call neighbor_interp_input(n,gridDesc, &
               MERRA2runoffdata_struc(n)%n11)
       else
          allocate(MERRA2runoffdata_struc(n)%n11(&
               MERRA2runoffdata_struc(n)%nc*MERRA2runoffdata_struc(n)%nr))
          call upscaleByAveraging_input(gridDesc,&
               LIS_rc%gridDesc(n,:),&
               MERRA2runoffdata_struc(n)%nc*MERRA2runoffdata_struc(n)%nr,&
               LIS_rc%lnc(n)*LIS_rc%lnr(n),MERRA2runoffdata_struc(n)%n11)
       endif
    enddo
    
  end subroutine MERRA2runoffdata_init
end module MERRA2runoffdataMod
