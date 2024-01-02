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
module GWBMIPrunoffdataMod
!BOP
! 
! !MODULE: GWBMIPrunoffdataMod
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
  public :: GWBMIPrunoffdata_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: GWBMIPrunoffdata_struc
  
  type, public :: GWBMIPrunoffdatadec
     
     real                    :: outInterval 
     character(len=LIS_CONST_PATH_LEN) :: odir 
     character*50            :: model_prefix
     integer                 :: nc, nr
     integer, allocatable    :: n11(:)
     type(ESMF_Time)         :: startTime
  end type GWBMIPrunoffdatadec

  type(GWBMIPrunoffdatadec), allocatable :: GWBMIPrunoffdata_struc(:)

contains
 
!BOP
!
! !ROUTINE: GWBMIPrunoffdata_init
! \label{GWBMIPrunoffdata_init}
! 
  subroutine GWBMIPrunoffdata_init
    !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod

    integer              :: n 
    integer              :: status
    character*10         :: time
    real                 :: gridDesc(50)

    allocate(GWBMIPrunoffdata_struc(LIS_rc%nnest))
       
    call ESMF_ConfigFindLabel(LIS_config,&
         "GWBMIP runoff data output directory:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            GWBMIPrunoffdata_struc(n)%odir,rc=status)
       call LIS_verify(status,&
            "GWBMIP runoff data output directory: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "GWBMIP runoff data model name prefix:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            GWBMIPrunoffdata_struc(n)%model_prefix,rc=status)
       call LIS_verify(status,&
            "GWBMIP runoff data model name prefix: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "GWBMIP runoff data output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "GWBMIP runoff data output interval: not defined")

       call LIS_parseTimeString(time,GWBMIPrunoffdata_struc(n)%outInterval)
    
       gridDesc = 0.0

       
       GWBMIPrunoffdata_struc(n)%nc = 360
       GWBMIPrunoffdata_struc(n)%nr = 180
       
       gridDesc(1) = 0  
       gridDesc(2) = GWBMIPrunoffdata_struc(n)%nc
       gridDesc(3) = GWBMIPrunoffdata_struc(n)%nr
       gridDesc(4) = -89.5
       gridDesc(5) = -179.5
       gridDesc(7) = 89.5
       gridDesc(8) = 179.5
       gridDesc(6) = 128
       gridDesc(9) = 1.0
       gridDesc(10) = 1.0
       gridDesc(20) = 0
       
       if(LIS_isAtAfinerResolution(n,1.0)) then
          
          allocate(GWBMIPrunoffdata_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call neighbor_interp_input(n,gridDesc, &
               GWBMIPrunoffdata_struc(n)%n11)
       else
          allocate(GWBMIPrunoffdata_struc(n)%n11(&
               GWBMIPrunoffdata_struc(n)%nc*GWBMIPrunoffdata_struc(n)%nr))
          call upscaleByAveraging_input(gridDesc,&
               LIS_rc%gridDesc(n,:),&
               GWBMIPrunoffdata_struc(n)%nc*GWBMIPrunoffdata_struc(n)%nr,&
               LIS_rc%lnc(n)*LIS_rc%lnr(n),GWBMIPrunoffdata_struc(n)%n11)
       endif
       
       call ESMF_TimeSet(GWBMIPrunoffdata_struc(n)%startTime,yy=1979, &
            mm = 1, &
            dd = 1, &
            h = 0, &
            m = 0, &
            calendar = LIS_calendar, &
            rc=status)
       call LIS_verify(status, 'Error in ESMF_TimeSet: GWBMIPrunoffdata_init')
    enddo
    
  end subroutine GWBMIPrunoffdata_init
end module GWBMIPrunoffdataMod
