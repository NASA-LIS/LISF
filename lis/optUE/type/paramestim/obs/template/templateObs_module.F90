!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: templateObs_module
! 
! !DESCRIPTION: 
!   This is a template for the definition of observations for optUE
!   
! !REVISION HISTORY: 
!  09 Jul 09    Sujay Kumar;   Initial Specification
! 
module templateObs_module
! !USES: 
  use ESMF
!EOP
  implicit none
  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: templateObs_setup
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------

contains
!BOP
! 
! !ROUTINE: templateObs_setup
! \label{templateObs_setup}
! 
! !INTERFACE: 
  subroutine templateObs_setup(Obs_State)
! !USES: 
    use LIS_coreMod,   only : LIS_rc
    use LIS_logMod,    only : LIS_verify

    implicit none 

! !ARGUMENTS: 
    type(ESMF_State)       ::  Obs_State(LIS_rc%nnest)
! 
! !DESCRIPTION: 
!   
!   This is a empty, template for the definition of 
!   observation (for optUE) data structures
!  
!   The arguments are: 
!   \begin{description}
!    \item[Obj\_Space]   observation/Objective space object 
!   \end{description}
!EOP
    integer :: n, status
    
    do n=1, LIS_rc%nnest
       call ESMF_AttributeSet(Obs_State(n),"Data Update Status",&
            .false., rc=status)
       call LIS_verify(status)
    enddo

  end subroutine templateObs_setup
  
end module TemplateObs_module
