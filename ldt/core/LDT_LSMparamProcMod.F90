!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_LSMparams_init
  public :: LDT_LSMparams_writeHeader
  public :: LDT_LSMparams_writeData

contains  
!BOP
! !ROUTINE: LDT_LSMparams_init
! \label{LDT_LSMparams_init}
!
! !INTERFACE: 
  subroutine LDT_LSMparams_init()
    
    if(LDT_rc%lsm.ne."none") then 
       call lsmparamprocinit(trim(LDT_rc%lsm)//char(0))
    endif
  end subroutine LDT_LSMparams_init


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
