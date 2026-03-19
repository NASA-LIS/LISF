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
! !ROUTINE: HYMAP2_getSWOTpred
! \label{HYMAP2_getSWOTpred}
!
! !REVISION HISTORY:
! 15 Apr 24: Yeosang Yoon; Initial specification;
!                          copied from HYMAP2_getWLpred
!
! !INTERFACE:
subroutine HYMAP2_getSWOTpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use HYMAP2_routingMod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the SWOT obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  integer                :: i,t,m
  real                   :: wl(HYMAP2_routing_struc(n)%nseqall*LIS_rc%nensem(n))


  do i=1,HYMAP2_routing_struc(n)%nseqall
     do m=1,LIS_rc%nensem(n)
        t=(i-1)*LIS_rc%nensem(n)+m
        wl(t) = HYMAP2_routing_struc(n)%sfcelv(i,m)
     enddo
  enddo

  call HYMAP2_convertPatchSpaceToObsEnsSpace(n,k,&
       wl,&
       obs_pred)
  
end subroutine HYMAP2_getSWOTpred
