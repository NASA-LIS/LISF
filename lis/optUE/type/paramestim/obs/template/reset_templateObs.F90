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
! !ROUTINE: reset_templateObs
! \label{reset_templateObs}
!
! !REVISION HISTORY:
!  09 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine reset_templateObs(Obj_Space)
! !USES: 
  use ESMF

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
! Empty template for reseting the optUE-observation data structures
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP


end subroutine reset_templateObs
