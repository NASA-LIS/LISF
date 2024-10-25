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
! !ROUTINE: cable_finalize
! \label{cable_finalize}
!
! !REVISION HISTORY:
!  26 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!
! !INTERFACE:
subroutine cable_finalize
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use cable_arrays
  use cable_lsmMod,       only : cable_struc
!
! !DESCRIPTION:
!  This subroutine cleans up the allocated memory structures in CABLE
!
!EOP
  implicit none
  
  integer :: n
  
  do n = 1,LIS_rc%nnest
     deallocate(cable_struc(n)%cable)
  enddo
  deallocate(cable_struc)
  
! Deallocate CABLE main variables:
  call dealloc_cbm_var(air, 1, model_structure)
  call dealloc_cbm_var(bgc, 1, model_structure)
  call dealloc_cbm_var(canopy, 1, model_structure)
  call dealloc_cbm_var(met, 1, model_structure)
  call dealloc_cbm_var(bal, 1, model_structure)
  call dealloc_cbm_var(rad, 1, model_structure)
  call dealloc_cbm_var(rough, 1, model_structure)
  call dealloc_cbm_var(soil, 1, model_structure)
  call dealloc_cbm_var(ssoil, 1, model_structure)
  call dealloc_cbm_var(sum_flux, 1, model_structure)
  call dealloc_cbm_var(veg, 1, model_structure)
  call dealloc_cbm_var(model_structure, 1)
  
end subroutine cable_finalize
