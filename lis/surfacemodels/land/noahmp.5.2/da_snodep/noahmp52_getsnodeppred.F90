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
! !ROUTINE: noahmp52_getsnodeppred
! \label{noahmp52_getsnodeppred}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine noahmp52_getsnodeppred(n, k, obs_pred)

! !USES:
  use LIS_coreMod, only : LIS_rc
  use noahmp52_lsmMod
  use SNODEPobs_Mod, only: SNODEP_obs_obj
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  real                   :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))
!EOP

  integer                :: t

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (SNODEP_obs_obj(n)%mesh .eq. 8) then
        snwd(t) = &
             NoahMP52_struc(n)%noahmp52(t)%snowh*39.37 !convert from m to inch
     elseif (SNODEP_obs_obj(n)%mesh .eq. 16 .or. &
          SNODEP_obs_obj(n)%mesh .eq. 25) then
        snwd(t) = NoahMP52_struc(n)%noahmp52(t)%snowh
     endif
  enddo

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       snwd,&
       obs_pred)

end subroutine noahmp52_getsnodeppred

