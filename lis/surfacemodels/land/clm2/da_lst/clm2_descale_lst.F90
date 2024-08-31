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
! !ROUTINE: clm2_descale_lst
! \label{clm2_descale_lst}
!
! !REVISION HISTORY:
! 1 Apr 2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine clm2_descale_lst(n, LSM_State)

! !USES:
  use ESMF

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  descales the land surface temperature related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[LSM\_State] ESMF State container for LSM state variables
!  \end{description}
!EOP

end subroutine clm2_descale_lst

