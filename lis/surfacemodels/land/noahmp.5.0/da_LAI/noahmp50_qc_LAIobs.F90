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
! !ROUTINE: noahmp50_qc_LAIobs
! \label{noahmp50_qc_LAIobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 1 Aug 2016: Mahdi Navari; Modified for Noahmp401 
!
! !INTERFACE:
subroutine noahmp50_qc_LAIobs(n,k,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_DAobservationsMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the soil moisture observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  soil is frozen or (3) ground is fully or partially covered 
!  with snow MN:(4) ground is covered with vegatation (more than 50%). 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP


end subroutine noahmp50_qc_LAIobs

