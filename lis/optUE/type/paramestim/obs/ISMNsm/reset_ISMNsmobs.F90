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
! !ROUTINE: reset_ISMNsmobs
! \label{reset_ISMNsmobs}
!
! !REVISION HISTORY:
!  21 Sep 2018: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine reset_ISMNsmobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,        only : LIS_rc
  use ISMNsm_obsMod,       only : ISMNsm_obs_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  resets the ISMN soil moisture data structures
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
     ISMNsm_obs_struc(n)%yr = -1
  enddo

end subroutine reset_ISMNsmobs
