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
! !ROUTINE: noahmp401_updateburnareavars
!  \label{noahmp401_updateburnareavars}
!
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp401_updateburnareavars(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP

end subroutine noahmp401_updateburnareavars

