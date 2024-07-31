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
! !MODULE: finalize_metForcGenerated
! \label{finalize_metForcGenerated}
! 
! 
! !INTERFACE:
subroutine finalize_metForcGenerated( findex )
! !USES:
  use LIS_coreMod, only : LIS_rc
  use metForcGen_SpatialInterpMod, only : metForcGen_interp_finalize
  use metForcGen_VariablesMod, only : forcvars_struc
!
!
!EOP  
  implicit none
  integer :: findex
  integer :: n

! Deallocate the forcing fields:
  deallocate( forcvars_struc )

! Deallocate the spatial interpolation arrays:
  call metForcGen_interp_finalize( findex )


end subroutine finalize_metForcGenerated
