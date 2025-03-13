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
! !ROUTINE: noahmp50_getsnodeppred
! \label{noahmp50_getsnodeppred}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!
! !INTERFACE:
subroutine noahmp50_getsnodeppred(n, k, obs_pred)

! !USES:
  use LIS_coreMod, only : LIS_rc
  use noahmp50_lsmMod
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
        snwd(t) = NoahMP50_struc(n)%noahmp50(t)%snowh*39.37 !convert from meter to inch
     elseif (SNODEP_obs_obj(n)%mesh .eq. 16 .or. SNODEP_obs_obj(n)%mesh .eq. 25) then
        snwd(t) = NoahMP50_struc(n)%noahmp50(t)%snowh
     endif
  enddo

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       snwd,&
       obs_pred)
  
end subroutine noahmp50_getsnodeppred

