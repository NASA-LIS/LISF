!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: finalize_climatology
! \label{finalize_climatology}
! 
! 
! !INTERFACE:
subroutine finalize_climatology( findex )
! !USES:
  use LIS_coreMod, only : LIS_rc
  use climatology_VariablesMod, only : climvars_struc
!  use climatology_SpatialInterpMod, only : climatology_interp_finalize
!
!
!EOP  
  implicit none
  integer :: findex
  integer :: n

! Deallocate the forcing fields:
  deallocate( climvars_struc )

! Deallocate the spatial interpolation arrays:
!  call climatology_interp_finalize( findex )


end subroutine finalize_climatology
