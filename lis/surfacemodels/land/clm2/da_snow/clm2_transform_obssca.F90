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
! !ROUTINE: clm2_transform_obssca
! \label{clm2_transform_obssca}
!
! !REVISION HISTORY:
! 25Jun2006: Sujay Kumar: Initial Specification
! 03May2007: K. Arseanult -- Added MODIS SCA (MOD10A1) Obs Option
! 19Feb2009: K. Arsenault -- Updated CLM2 DI Code for SCF; based on G. de Lannoy's code.
!
! !INTERFACE:
subroutine clm2_transform_obssca(n,OBS_State)

! !USES:
  use ESMF
  use LIS_coreMod,  only : LIS_rc
  use LIS_logMod,   only : LIS_verify
!EOP
  implicit none

  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the MODIS SCA observation state 
!  to the LSM SWE state.
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP
  type(ESMF_Field)         :: obs_sca_field

  real, pointer            :: scaobs(:)
  integer                  :: t
  integer                  :: N_obs_size
  integer                  :: status

! -------------------------------------------------------------------------------

  call ESMF_AttributeGet( OBS_State, name="Number Of Observations",&
       value=N_obs_size, rc=status)
  call LIS_verify(status, 'attributeget error in clm2_transform_obssca')
  
  call ESMF_StateGet( OBS_State, "Observation01", obs_sca_field,&
       rc=status)
  call LIS_verify(status, 'stateget error in clm2_transform_obssca')

  call ESMF_FieldGet( obs_sca_field, localDE=0, farrayPtr=scaobs,rc=status)
  call LIS_verify(status,'fieldget error in clm2_transform_obssca')
  
  do t = 1, N_obs_size
     if ( scaobs(t) .ne.LIS_rc%udef) then 
        scaobs(t) = scaobs(t)/1000.0    ! IF Obs are originally in mm --> m
     endif
  enddo
  

end subroutine clm2_transform_obssca

