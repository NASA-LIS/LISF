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
! !ROUTINE: noahmp401_getsnippred
! \label{noahmp401_getsnippred}
!
! !REVISION HISTORY:
! 18 Jul 2025: Eric Kemp; Initial specification (copied from USAFSI version)
!
! !INTERFACE:
subroutine noahmp401_getsnippred(n, k, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_DAobservationsMod
  use noahmp401_lsmMod

  ! Defaults
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  real                   :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))

! !DESCRIPTION:
!  This routine computes the obspred ('Hx') term for SNIP DA assimilation
!  instances.
!EOP

  integer                :: t

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     snwd(t) = noahmp401_struc(n)%noahmp401(t)%snowh ! Keep in meters
  enddo

  call LIS_convertPatchSpaceToObsEnsSpace(n, k, &
       LIS_rc%lsm_index, &
       snwd, &
       obs_pred)

end subroutine noahmp401_getsnippred

