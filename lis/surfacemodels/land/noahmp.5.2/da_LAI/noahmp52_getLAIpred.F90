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
! !ROUTINE: noahmp52_getLAIpred
! \label{noahmp52_getLAIpred}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine noahmp52_getLAIpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_DAobservationsMod
  use noahmp52_lsmMod
!EOP

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the Soil moisture obs pred (model's estimate of
!  observations) for data assimilation
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  integer                :: t
  real                   :: lai(LIS_rc%npatch(n,LIS_rc%lsm_index))

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     lai(t) = NoahMP52_struc(n)%noahmp52(t)%lai
  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       lai,&
       obs_pred)

end subroutine noahmp52_getLAIpred

