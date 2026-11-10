!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp52_getsnowpred
! \label{noahmp52_getsnowpred}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; Modified for refactored NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine noahmp52_getsnowpred(n, k, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noahmp52_lsmMod
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  real                   :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))
!EOP

! !DESCRIPTION:
!  This routine computes the obspred ('Hx') term for SNOW DA assimilation
!  instances.

  integer                :: t

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     snwd(t) = NoahMP52_struc(n)%noahmp52(t)%snowh ! Keep in meters
  enddo

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       snwd,&
       obs_pred)

end subroutine noahmp52_getsnowpred

