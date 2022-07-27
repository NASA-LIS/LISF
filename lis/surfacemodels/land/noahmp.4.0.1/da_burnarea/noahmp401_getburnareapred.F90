!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_getburnareapred
! \label{noahmp401_getburnareapred}
!
! !REVISION HISTORY:
! 24 Jul 2022: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp401_getburnareapred(n, k,obs_pred)
! !USES:
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_dataAssimMod
  use LIS_DAobservationsMod
  use noahmp401_lsmMod
  use noahmp401_dasoilm_Mod
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

  
end subroutine noahmp401_getburnareapred

