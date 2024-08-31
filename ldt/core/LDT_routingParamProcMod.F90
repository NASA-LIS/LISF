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
module LDT_routingParamProcMod
!BOP
!
! !MODULE: LDT_routingParamProcMod
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

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_routingParams_init
  public :: LDT_routingParams_writeHeader
  public :: LDT_routingParams_writeData

contains  
!BOP
! !ROUTINE: LDT_routingParams_init
! \label{LDT_routingParams_init}
!
! !INTERFACE: 
  subroutine LDT_routingParams_init()
    
    if(LDT_rc%routingmodel.ne."none") then 
       call routingparamprocinit(trim(LDT_rc%routingmodel)//char(0))
    endif
  end subroutine LDT_routingParams_init


!BOP
! !ROUTINE: LDT_routingParams_writeHeader
! \label{LDT_routingParams_writeHeader}
!
! !INTERFACE: 
  subroutine LDT_routingParams_writeHeader(n,ftn,dimID, monthID)
    integer     :: n
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: monthID

    if(LDT_rc%routingmodel.ne."none") then 
       call routingparamprocwriteheader(trim(LDT_rc%routingmodel)//char(0),&
            n,ftn,dimID, monthID)
    endif
  end subroutine LDT_routingParams_writeHeader


!BOP
! !ROUTINE: LDT_routingParams_writeData
! \label{LDT_routingParams_writeData}
!
! !INTERFACE: 
  subroutine LDT_routingParams_writeData(n,ftn)

    integer     :: n
    integer     :: ftn

    if(LDT_rc%routingmodel.ne."none") then 
       call routingparamprocwritedata(trim(LDT_rc%routingmodel)//char(0),&
            n,ftn)
    endif
  end subroutine LDT_routingParams_writeData

  
end module LDT_routingParamProcMod
