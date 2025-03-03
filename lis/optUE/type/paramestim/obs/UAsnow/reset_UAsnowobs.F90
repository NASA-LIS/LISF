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
! !ROUTINE: reset_UAsnowobs
! \label{reset_UAsnowobs}
!
! !REVISION HISTORY:
!  2 May 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine reset_UAsnowobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,        only : LIS_rc
  use UAsnow_obsMod,       only : UAsnow_obs_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  resets the UA snow data structures for parameter
!  optimization
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP

  integer            :: n 
  
end subroutine reset_UAsnowobs
