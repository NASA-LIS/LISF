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
module LISrunoffdataMod
!BOP
! 
! !MODULE: LISrunoffdataMod
! 
! !DESCRIPTION: 
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
  public :: LISrunoffdata_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: LISrunoffdata_struc
  
  type, public :: LISrunoffdatadec
     
     real             :: outInterval 
     character(len=LIS_CONST_PATH_LEN) :: odir
  end type LISrunoffdatadec

  type(LISrunoffdatadec), allocatable :: LISrunoffdata_struc(:)

contains
 
!BOP
!
! !ROUTINE: LISrunoffdata_init
! \label{LISrunoffdata_init}
! 
  subroutine LISrunoffdata_init
    !USES: 
    use LIS_coreMod
    use LIS_logMod
    use LIS_timeMgrMod

    integer              :: n 
    integer              :: status
    character*10         :: time

    allocate(LISrunoffdata_struc(LIS_rc%nnest))
       
    call ESMF_ConfigFindLabel(LIS_config,&
         "LIS runoff output directory:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,LISrunoffdata_struc(n)%odir,rc=status)
       call LIS_verify(status,&
            "LIS runoff output directory: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "LIS runoff output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "LIS runoff output interval: not defined")

       call LIS_parseTimeString(time,LISrunoffdata_struc(n)%outInterval)
    enddo
    
    
  end subroutine LISrunoffdata_init
end module LISrunoffdataMod
