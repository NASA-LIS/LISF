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
! !ROUTINE: noahmp50_getusafsipred
! \label{noahmp50_getusafsipred}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; update to work with refactored NoahMP (v5.0 and newer)
!
! !INTERFACE:
subroutine noahmp50_getusafsipred(n, k, obs_pred)

! !USES:
  use LIS_coreMod, only : LIS_rc
  use NoahMP50_lsmMod
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  real                   :: snwd(LIS_rc%npatch(n,LIS_rc%lsm_index))
!EOP

! !DESCRIPTION:
!  This routine computes the obspred ('Hx') term for USAFSI DA assimilation
!  instances.

  integer                :: t

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     snwd(t) = Noahmp50_struc(n)%noahmp50(t)%snowh ! Keep in meters
  enddo

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       snwd,&
       obs_pred)
  
end subroutine noahmp50_getusafsipred

