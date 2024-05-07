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
! !ROUTINE: clm2_getsynsmpred
! \label{clm2_getsynsmpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clm2_getsynsmpred(n, obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use clm2_lsmMod
  use clm2_varcon,    only : denice, denh2o
!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
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
  integer                :: count1(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  integer                :: i,t,m,gid
  real                   :: vol_ice, vol_liq, eff_porosity

  count1 = 0 
  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = i+m-1
        gid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
        vol_ice = min(clm2_struc(n)%clm(t)%watsat(1), &
             clm2_struc(n)%clm(t)%h2osoi_ice(1)/&
             (clm2_struc(n)%clm(t)%dz(1)*denice))
        eff_porosity = clm2_struc(n)%clm(t)%watsat(1)-vol_ice
        vol_liq = min(eff_porosity, &
             clm2_struc(n)%clm(t)%h2osoi_liq(1)/&
             (clm2_struc(n)%clm(t)%dz(1)*denh2o))
        
        obs_pred(gid,m)= obs_pred(gid,m) + vol_liq + vol_ice
        count1(gid,m) = count1(gid,m) + 1
     enddo
  enddo
  
  do i=1,LIS_rc%ngrid(n)
     do m=1,LIS_rc%nensem(n)
        obs_pred(i,m) = obs_pred(i,m)/(count1(i,m))
     enddo
  enddo
end subroutine clm2_getsynsmpred

