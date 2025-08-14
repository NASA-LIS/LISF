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
! !ROUTINE: noah39_qc_snipobs
! \label{noah39_qc_snipobs}
!
! !REVISION HISTORY:
! 17 Jul 2025: Eric Kemp; Initial specification (copied from USAFSI version)
!
! !INTERFACE:
subroutine noah39_qc_snipobs(n, OBS_State)
  ! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use noah39_lsmMod

  ! Defaults
  implicit none

  ! !ARGUMENTS:
  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
  !
  ! !DESCRIPTION:
  !
  !  This subroutine performs any model-based QC of the observation
  !  prior to data assimilation. Here the snow observations
  !  are flagged when LSM indicates that (1) rain is falling (2)
  !  ground is fully or partially covered with snow.
  !
  !  The arguments are:
  !  \begin{description}
  !  \item[n] index of the nest \newline
  !  \item[OBS\_State] ESMF state container for observations \newline
  !  \end{description}
  !
  !EOP
  type(ESMF_Field)         :: obs_snow_field

  real, pointer            :: snipobs(:)
  integer                  :: status


  call ESMF_StateGet(OBS_State, "Observation01", obs_snow_field, &
       rc=status)
  call LIS_verify(status, &
       "ESMF_StateGet failed in noah39_qc_snipobs")
  call ESMF_FieldGet(obs_snow_field, localDE=0, farrayPtr=snipobs, &
       rc=status)
  call LIS_verify(status, &
       "ESMF_FieldGet failed in noah39_qc_snipobs")

end subroutine noah39_qc_snipobs

