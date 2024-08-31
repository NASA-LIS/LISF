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
! !ROUTINE: jules50_getsmpred
! \label{jules50_getsmpred}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 23Apr2018: Mahdi Navari: Modified for JULES 5.0 
!
! !INTERFACE:
subroutine jules50_getsmpred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use jules50_lsmMod
  use jules50_dasoilm_Mod
  !use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use jules_soil_mod, only:  dzsoil 
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
  real                   :: smc1(LIS_rc%npatch(n,LIS_rc%lsm_index))

  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)
! apply the unit conversion.
! In the EnKF equation X(+) = X(-) + K(y-Hx), Units of the terms y and
! Hx must be the same. we just need to change the unit of Hx. The Unit of Hx
! from JULES is [kg/m2] therefore we need to change the unit to
! be consistent with that of y [m3/m3].
! for the unit conversion FROM [kg/m2] to [M3/M3] we need layer thickness.
! Here the thickness of the first layer. (which is 10 cm in both PILDAS and JULES) has been used
!     smc1(t) = jules50_struc(n)%jules50(t)%smcl_soilt(1) * 1/LIS_sfmodel_struc(n)%lyrthk(1)*1/1000  ! [kg/m2]*1/m*1/(kg/m3) --> [m3/m3]

!MN: Eric changed the unit of LIS_sfmodel_struc%lyrthk so LIS_sfmodel_struc%lyrthk is in cm, we
!use dzsoil instead.

 smc1(t) = jules50_struc(n)%jules50(t)%smcl_soilt(1) * 1/dzsoil(1)*1/1000  ! [kg/m2]*1/m*1/(kg/m3) --> [m3/m3] 

  enddo
  call LIS_convertPatchSpaceToObsEnsSpace(n,k,&
       LIS_rc%lsm_index, &
       smc1,&
       obs_pred)

end subroutine jules50_getsmpred

