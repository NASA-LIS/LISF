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
! !ROUTINE: clsmf25_getsmpred
! \label{clsmf25_getsmpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine clsmf25_getsmpred(n, k, obs_pred)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc,LIS_surface
  use clsmf25_lsmMod
  use clsmf25_model
  use LIS_DAobservationsMod

!EOP

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  integer, intent(in)    :: k
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))
  integer                :: count1(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n))
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
  real                   :: sfmc(LIS_rc%npatch(n,LIS_rc%lsm_index))

  CALL CALC_SOIL_MOIST (                       &
       LIS_rc%npatch(n,LIS_rc%lsm_index),      &
       clsmf25_struc(n)%cat_param%vegcls,  &
       clsmf25_struc(n)%cat_param%dzsf,    &
       clsmf25_struc(n)%cat_param%vgwmax,  &
       clsmf25_struc(n)%cat_param%cdcr1,   &
       clsmf25_struc(n)%cat_param%cdcr2,   &
       clsmf25_struc(n)%cat_param%wpwet,   &
       clsmf25_struc(n)%cat_param%poros,   &
       clsmf25_struc(n)%cat_param%psis,    &
       clsmf25_struc(n)%cat_param%bee,     &
       clsmf25_struc(n)%cat_param%ars1,    & 
       clsmf25_struc(n)%cat_param%ars2,    & 
       clsmf25_struc(n)%cat_param%ars3,    & 
       clsmf25_struc(n)%cat_param%ara1,    & 
       clsmf25_struc(n)%cat_param%ara2,    & 
       clsmf25_struc(n)%cat_param%ara3,    & 
       clsmf25_struc(n)%cat_param%ara4,    & 
       clsmf25_struc(n)%cat_param%arw1,    & 
       clsmf25_struc(n)%cat_param%arw2,    & 
       clsmf25_struc(n)%cat_param%arw3,    & 
       clsmf25_struc(n)%cat_param%arw4,    &
       clsmf25_struc(n)%cat_progn%srfexc,  &
       clsmf25_struc(n)%cat_progn%rzexc,   &
       clsmf25_struc(n)%cat_progn%catdef,  &
       sfmc,    & 
       clsmf25_struc(n)%cat_diagn%rzmc,    & 
       clsmf25_struc(n)%cat_diagn%prmc)

  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       sfmc, &
       obs_pred)

end subroutine clsmf25_getsmpred

