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

module accumDGVM

#if (defined DGVM)

  use LIS_precisionMod
  implicit none

! public methods

  public :: accumDGVMini

!=======================================================================
CONTAINS
!=======================================================================

  subroutine accumDGVMini
    
    use LIS_precisionMod
    use clm2_lsmMod
    use accumulMod, only : accext
    use clm_varmap, only : begpatch,endpatch
    use LIS_timeMgrMod, only: get_nstep
    implicit none
    
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize time variant DGVM variables
! 
! Method: 
! This routine must be called after the restart file is read
! because the time manager is initialized in the restart file
! and that is needed to obtain the time step
! 
! Author: Mariana Vertenstein
! 
!-----------------------------------------------------------------------
! $Id: accumDGVM.F90,v 1.6 2004/11/24 22:57:03 jim Exp $
!-----------------------------------------------------------------------

! ------------------------ local variables ---------------------------
    integer nstep ! time step
!-----------------------------------------------------------------------

! Initialize variables to be time accumulated for various purposes 

    nstep = get_nstep()
    call accext ('T10     ', clm(begpatch:endpatch)%t10    , nstep)
    call accext ('TDA     ', clm(begpatch:endpatch)%t_mo   , nstep)
    call accext ('AGDD0   ', clm(begpatch:endpatch)%agdd0  , nstep)
    call accext ('AGDD5   ', clm(begpatch:endpatch)%agdd5  , nstep)
    call accext ('FNPSN10 ', clm(begpatch:endpatch)%fnpsn10, nstep)
    call accext ('PREC365 ', clm(begpatch:endpatch)%prec365, nstep)
    call accext ('AGDDTW  ', clm(begpatch:endpatch)%agddtw , nstep)
    call accext ('AGDD    ', clm(begpatch:endpatch)%agdd   , nstep)
    
    return
  end subroutine accumDGVMini

#endif

end module accumDGVM
