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
! !ROUTINE: clm2_getsynlstpred
! \label{clm2_getsynlstpred}
!
! !REVISION HISTORY:
! 1 Apr 2007: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clm2_getsynlstpred(n, obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use clm2_lsmMod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the LST obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  integer                :: i,t,m,gid

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
        obs_pred(gid,m)= clm2_struc(n)%clm(t)%t_veg
     enddo
  enddo
  
end subroutine clm2_getsynlstpred
