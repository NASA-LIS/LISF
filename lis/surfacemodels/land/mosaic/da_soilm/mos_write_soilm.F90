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
! !ROUTINE: mos_write_soilm
!  \label{mos_write_soilm}
!
! !REVISION HISTORY:
! 06 Oct 2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine mos_write_soilm(ftn, n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use mos_lsmMod
  use LIS_logMod,    only : LIS_logunit
  use LIS_historyMod, only : LIS_writevar_restart

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: ftn
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to mosaic's
!  model space. 
! 
!EOP
  call LIS_writevar_restart(ftn, n, mos_struc(n)%mos%water1/20.0)
  call LIS_writevar_restart(ftn, n, mos_struc(n)%mos%water2/1480.0)
  call LIS_writevar_restart(ftn, n, mos_struc(n)%mos%water3/2000.0)

end subroutine mos_write_soilm

