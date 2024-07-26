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
! !ROUTINE: clsmf25_getsnwdpred
! \label{clsmf25_getsnwdpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!
! !INTERFACE:
subroutine clsmf25_getsnwdpred(n, obs_pred)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use clsmf25_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
!EOP

  integer                :: i,t,m,gid

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_surface(n,1)%tile(t)%index
        obs_pred(gid,m)= &
             ((clsmf25_struc(n)%cat_progn(t)%sndz(1)*1000.0) + &
             (clsmf25_struc(n)%cat_progn(t)%sndz(2)*1000.0) + &
             (clsmf25_struc(n)%cat_progn(t)%sndz(3)*1000.0))
     enddo
  enddo
  
end subroutine clsmf25_getsnwdpred

