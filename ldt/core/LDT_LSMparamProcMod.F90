!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"
module LDT_LSMparamProcMod
!BOP
!
! !MODULE: LDT_LSMparamProcMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage the processing of 
!   LSM parameters
!
! !REVISION HISTORY: 
!  14 Aug 2014:  Sujay Kumar;  Initial Specification
! 
  use ESMF
  use LDT_coreMod
  use LDT_logMod

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_LSMparams_init
  public :: LDT_LSMparams_writeHeader
  public :: LDT_LSMparams_writeData

!BOP 
! 
! !ROUTINE: LDT_LSMparams_init
! \label{LDT_LSMparams_init}
! 
! !INTERFACE:
  interface LDT_LSMparams_init
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure LSMparams_init_LIS
     module procedure LSMparams_init_LISHydro
! 
! !DESCRIPTION:
! This interface provides routines for writing LSM parameters in both 
! in the preprocessing mode for LIS as well as in the LISHydro(WRFHydro) 
! preprocessing mode. 
!EOP 
  end interface

contains  
!BOP
! !ROUTINE: LSMparams_init_LIS
! \label{LSMparams_init_LIS}
!
! !INTERFACE: 
  subroutine LSMparams_init_LIS()
    integer :: flag
    flag = 0

    if(LDT_rc%lsm.ne."none") then 
       call lsmparamprocinit(trim(LDT_rc%lsm)//char(0),flag)
    endif
  end subroutine LSMparams_init_LIS


!BOP
! !ROUTINE: LSMparams_init_LISHydro
! \label{LSMparams_init_LISHydro}
!
! !INTERFACE: 
  subroutine LSMparams_init_LISHydro(flag)
    
    integer   :: flag

    flag = 1

    if(LDT_rc%lsm.ne."none") then 
       if(LDT_rc%lsm.ne."Noah.2.7.1".or.&
            LDT_rc%lsm.ne."Noah.3.2".or.&
            LDT_rc%lsm.ne."Noah.3.3".or.&
            LDT_rc%lsm.ne."Noah.3.6".or.&
            LDT_rc%lsm.ne."Noah.3.9".or.&
            LDT_rc%lsm.ne."Noah-MP.3.6".or.&
            LDT_rc%lsm.ne."Noah-MP.4.0.1") then 

          call lsmparamprocinit(trim(LDT_rc%lsm)//char(0),flag)
       else
          write(LDT_logunit,*)'[ERR] Support for this LSM in the LISHydro preprocessing mode'
          write(LDT_logunit,*)'[ERR] is not supported'
          call LDT_endrun()
       endif
    endif
  end subroutine LSMparams_init_LISHydro



!BOP
! !ROUTINE: LDT_LSMparams_writeHeader
! \label{LDT_LSMparams_writeHeader}
!
! !INTERFACE: 
  subroutine LDT_LSMparams_writeHeader(n,ftn,dimID, monthID)
    integer     :: n
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: monthID

    if(LDT_rc%lsm.ne."none") then 
       call lsmparamprocwriteheader(trim(LDT_rc%lsm)//char(0),&
            n,ftn,dimID, monthID)
    endif

  end subroutine LDT_LSMparams_writeHeader


!BOP
! !ROUTINE: LDT_LSMparams_writeData
! \label{LDT_LSMparams_writeData}
!
! !INTERFACE: 
  subroutine LDT_LSMparams_writeData(n,ftn)

    integer     :: n
    integer     :: ftn

    if(LDT_rc%lsm.ne."none") then 
       call lsmparamprocwritedata(trim(LDT_rc%lsm)//char(0),&
            n,ftn)
    endif
  end subroutine LDT_LSMparams_writeData

  
end module LDT_LSMparamProcMod
