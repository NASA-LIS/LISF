!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp50_descale_usafsi
! \label{noahmp50_descale_usafsi}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; update to work with refactored NoahMP (v5.0 and newer)
!
! !INTERFACE:
subroutine noahmp50_descale_usafsi(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Returns the snow related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP

end subroutine noahmp50_descale_usafsi

