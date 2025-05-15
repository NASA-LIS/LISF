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
module RAPID_parmsMod
!BOP
!
! !MODULE: RAPID_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to process RAPID 
!  parameter data. 
!
! !REVISION HISTORY:
!
!  23 Apr 2025: Yeosang Yoon: Initial implementation
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: RAPIDparms_init    !allocates memory for required structures
  public :: RAPIDparms_writeHeader
  public :: RAPIDparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: RAPID_struc

  type, public :: rapid_type_dec
     character*50  :: rapid_proj
  end type rapid_type_dec

  type(rapid_type_dec), allocatable :: RAPID_struc(:)

contains

!BOP
! 
! !ROUTINE: RAPIDparms_init
! \label{RAPIDparms_init}
! 
! !INTERFACE:
  subroutine RAPIDparms_init
! !USES:
!  none
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the RAPID datasets 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[RAPIDparmssetup](\ref{RAPIDparmssetup}) \newline
!    calls the registry to invoke the RAPIDparms setup methods. 
!  \end{description}
!
!EOP
  end subroutine RAPIDparms_init

  subroutine RAPIDparms_writeHeader(n,ftn,dimID,monthID)

    integer   :: n
    integer   :: ftn
    integer   :: dimID(3)
    integer   :: monthID

  end subroutine RAPIDparms_writeHeader
  
  subroutine RAPIDparms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

  end subroutine RAPIDparms_writeData

end module RAPID_parmsMod
