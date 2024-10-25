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
! !ROUTINE: noahmp36_gettwspred
! \label{noahmp36_gettwspred}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp36_gettwspred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_DAobservationsMod
  use noahmp36_tws_DAlogMod
  use noahmp36_lsmMod

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
  real                   :: tws(LIS_rc%npatch(n,LIS_rc%lsm_index))

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
     tws(t) = (NOAHMPpred_struc(n)%clmnwater(1,t) + &	
          NOAHMPpred_struc(n)%clmnwater(2,t)+&
          NOAHMPpred_struc(n)%clmnwater(3,t))/3.	
  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       tws,&
       obs_pred)
end subroutine noahmp36_gettwspred

