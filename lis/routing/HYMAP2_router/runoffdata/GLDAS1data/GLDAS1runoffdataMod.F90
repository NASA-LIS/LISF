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
module GLDAS1runoffdataMod
!BOP
! 
! !MODULE: GLDAS1runoffdataMod
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
  public :: GLDAS1runoffdata_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: GLDAS1runoffdata_struc
  
  type, public :: GLDAS1runoffdatadec
     
     real                    :: outInterval 
     character(len=LIS_CONST_PATH_LEN) :: odir 
     character*20            :: model_name
     real                    :: datares
     integer                 :: nc, nr
     integer, allocatable    :: n11(:)
  end type GLDAS1runoffdatadec

  type(GLDAS1runoffdatadec), allocatable :: GLDAS1runoffdata_struc(:)

contains
 
!BOP
!
! !ROUTINE: GLDAS1runoffdata_init
! \label{GLDAS1runoffdata_init}
! 
  subroutine GLDAS1runoffdata_init
    !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod

    integer              :: n 
    integer              :: status
    character*10         :: time
    real                 :: gridDesc(50)

    allocate(GLDAS1runoffdata_struc(LIS_rc%nnest))
       
    call ESMF_ConfigFindLabel(LIS_config,&
         "GLDAS1 runoff data output directory:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,GLDAS1runoffdata_struc(n)%odir,rc=status)
       call LIS_verify(status,&
            "GLDAS1 runoff data output directory: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "GLDAS1 runoff data model name:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            GLDAS1runoffdata_struc(n)%model_name,rc=status)
       call LIS_verify(status,&
            "GLDAS1 runoff data model name: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "GLDAS1 runoff data spatial resolution (degree):",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            GLDAS1runoffdata_struc(n)%datares,rc=status)
       call LIS_verify(status,&
            "GLDAS1 runoff data spatial resolution (degree): not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "GLDAS1 runoff data output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "GLDAS1 runoff data output interval: not defined")

       call LIS_parseTimeString(time,GLDAS1runoffdata_struc(n)%outInterval)
    
       gridDesc = 0.0

       if(GLDAS1runoffdata_struc(n)%datares .eq. 0.25) then 
          
          GLDAS1runoffdata_struc(n)%nc = 1440
          GLDAS1runoffdata_struc(n)%nr = 600
          
          gridDesc(1) = 0  
          gridDesc(2) = GLDAS1runoffdata_struc(n)%nc
          gridDesc(3) = GLDAS1runoffdata_struc(n)%nr
          gridDesc(4) = -59.875
          gridDesc(5) = -179.875
          gridDesc(7) = 89.875
          gridDesc(8) = 179.875
          gridDesc(6) = 128
          gridDesc(9) = 0.25
          gridDesc(10) = 0.25
          gridDesc(20) = 64
          
          
       elseif(GLDAS1runoffdata_struc(n)%datares.eq. 1.0) then 
          
          GLDAS1runoffdata_struc(n)%nc = 360
          GLDAS1runoffdata_struc(n)%nr = 150
          
          gridDesc(1) = 0  
          gridDesc(2) = GLDAS1runoffdata_struc(n)%nc
          gridDesc(3) = GLDAS1runoffdata_struc(n)%nr
          gridDesc(4) = -59.5
          gridDesc(5) = -179.5
          gridDesc(7) = 89.5
          gridDesc(8) = 179.5
          gridDesc(6) = 128
          gridDesc(9) = 1.0
          gridDesc(10) = 1.0
          gridDesc(20) = 64
          
       endif
       
       if(LIS_isAtAfinerResolution(n,GLDAS1runoffdata_struc(n)%datares)) then
          
          allocate(GLDAS1runoffdata_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
          call neighbor_interp_input(n,gridDesc, &
               GLDAS1runoffdata_struc(n)%n11)
       else
          allocate(GLDAS1runoffdata_struc(n)%n11(&
               GLDAS1runoffdata_struc(n)%nc*GLDAS1runoffdata_struc(n)%nr))
          call upscaleByAveraging_input(gridDesc,&
               LIS_rc%gridDesc(n,:),&
               GLDAS1runoffdata_struc(n)%nc*GLDAS1runoffdata_struc(n)%nr,&
               LIS_rc%lnc(n)*LIS_rc%lnr(n),GLDAS1runoffdata_struc(n)%n11)
       endif
    enddo
    
  end subroutine GLDAS1runoffdata_init
end module GLDAS1runoffdataMod
