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
! !ROUTINE: noah271_getmultismpred
! \label{noah271_getmultismpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah271_getmultismpred(n, obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_domain
  use noah271_lsmMod
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  real                   :: obs_pred(LIS_rc%ngrid(n)*4,LIS_rc%nensem(n))
  integer                :: count1(LIS_rc%ngrid(n)*4,LIS_rc%nensem(n))
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
  integer                :: i,t,m,gid
  integer                :: tid

  obs_pred = 0.0
  count1 = 0.0 

  do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_domain(n)%tile(t)%index   
        tid = gid
        obs_pred(tid,m)= obs_pred(tid,m)+ noah271_struc(n)%noah(t)%smc(1)
        count1(tid,m) = count1(tid,m) +1
     enddo
  enddo

  do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_domain(n)%tile(t)%index   
        tid = LIS_rc%ngrid(n)*1 + gid
        obs_pred(tid,m)= obs_pred(tid,m)+ noah271_struc(n)%noah(t)%smc(2)
        count1(tid,m) = count1(tid,m) +1
     enddo
  enddo

  do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_domain(n)%tile(t)%index   
        tid = LIS_rc%ngrid(n)*2 + gid
        obs_pred(tid,m)= obs_pred(tid,m)+ noah271_struc(n)%noah(t)%smc(3)
        count1(tid,m) = count1(tid,m) +1
     enddo
  enddo

  do i=1,LIS_rc%ntiles(n),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_domain(n)%tile(t)%index   
        tid = LIS_rc%ngrid(n)*3 + gid
        obs_pred(tid,m)= obs_pred(tid,m)+ noah271_struc(n)%noah(t)%smc(4)
        count1(tid,m) = count1(tid,m) +1
     enddo
  enddo
  
  do i=1,LIS_rc%ngrid(n)*4
     do m=1,LIS_rc%nensem(n)
        obs_pred(i,m) = obs_pred(i,m)/(count1(i,m))
     enddo
  enddo

end subroutine noah271_getmultismpred

