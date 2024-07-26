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
! !ROUTINE: noah271_getTskinPred
! \label{noah271_getTskinPred}
!
! !REVISION HISTORY:
! 13Nov2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noah271_getTskinPred(n, obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use noah271_lsmMod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
!
! !DESCRIPTION:
!
!  Returns the Tskin obs pred (model's estimate of 
!  observations) for data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[obs\_pred] model's estimate of observations \newline
!  \end{description}
!EOP
  integer                :: i,t,m,gid

  do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_domain(n)%tile(t)%index
        obs_pred(gid,m)= noah271_struc(n)%noah(t)%t1
!        obs_pred(gid,m)= noah271_struc(n)%noah(t)%stc(1)
     enddo
  enddo
  
end subroutine noah271_getTskinPred

