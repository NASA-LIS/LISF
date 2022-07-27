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
! !ROUTINE: noahmp401_setburnareavars
!  \label{noahmp401_setburnareavars}
!
! !REVISION HISTORY:
! 24 Jul 2022: Sujay Kumar; Initial Specification
! 

! !INTERFACE:
subroutine noahmp401_setburnareavars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP

end subroutine noahmp401_setburnareavars

