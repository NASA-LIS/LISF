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
module LDT_lakeParamProcMod
!BOP
!
! !MODULE: LDT_lakeParamProcMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage the processing of 
!   Lake model parameters
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
  public :: LDT_lakeParams_init
  public :: LDT_lakeParams_writeHeader
  public :: LDT_lakeParams_writeData

contains  
!BOP
! !ROUTINE: LDT_lakeParams_init
! \label{LDT_lakeParams_init}
!
! !INTERFACE: 
  subroutine LDT_lakeParams_init()
    
    if(LDT_rc%lakemodel.ne."none") then 
       call lakeparamprocinit(trim(LDT_rc%lakemodel)//char(0))
    endif

  end subroutine LDT_lakeParams_init


!BOP
! !ROUTINE: LDT_lakeParams_writeHeader
! \label{LDT_lakeParams_writeHeader}
!
! !INTERFACE: 
  subroutine LDT_lakeParams_writeHeader(n,ftn,dimID, monthID)
    integer     :: n
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: monthID

    if(LDT_rc%lakemodel.ne."none") then 
       call lakeparamprocwriteheader(trim(LDT_rc%lakemodel)//char(0),&
            n,ftn,dimID, monthID)
    endif

  end subroutine LDT_lakeParams_writeHeader


!BOP
! !ROUTINE: LDT_lakeParams_writeData
! \label{LDT_lakeParams_writeData}
!
! !INTERFACE: 
  subroutine LDT_lakeParams_writeData(n,ftn)

#if ( defined SPMD )
    use mpi
#endif

    integer   :: n
    integer   :: ftn
    integer   :: ierr

    if(LDT_rc%lakemodel.ne."none") then 

       print *, '(LDT_lakeParamProcMod): processed lake model data ', LDT_localPet
#if ( defined SPMD )
       call MPI_Barrier(MPI_COMM_WORLD, ierr)
       print *, ierr
#endif

       call lakeparamprocwritedata(trim(LDT_rc%lakemodel)//char(0),&
            n,ftn)

    endif

  end subroutine LDT_lakeParams_writeData
  
end module LDT_lakeParamProcMod
