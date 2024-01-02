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
! !ROUTINE: reset_ARSsmobs
! \label{reset_ARSsmobs}
!
! !REVISION HISTORY:
!  02 Feb 2018: Soni Yatheendradas; Initial Specification
!
! !INTERFACE: 
subroutine reset_ARSsmobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,        only : LIS_rc
  use ARSsm_obsMod,       only : ARSsm_obs_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  resets the ARS soil moisture data structures
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_State] observations state
!  \end{description}
!
!EOP

  integer            :: n 
  
  do n=1,LIS_rc%nnest
     ARSsm_obs_struc(n)%yr = -1
  enddo

end subroutine reset_ARSsmobs
