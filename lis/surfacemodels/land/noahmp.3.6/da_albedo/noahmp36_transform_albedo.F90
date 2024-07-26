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
! !ROUTINE: noahmp36_transform_albedo
! \label{noahmp36_transform_albedo}
!
! !REVISION HISTORY:
!  22 Dec 2017: Sujay Kumar, Initial specification 
!
! !INTERFACE:
subroutine noahmp36_transform_albedo(n,OBS_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use noahmp36_lsmMod
!EOP
  implicit none

  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine transforms the albedo state
!  (mm) to the lsm state
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!EOP

!nothing to be done

#if 0 
  type(ESMF_Field)         :: obs_lai_field
  real, pointer            :: laiobs(:)
  integer                  :: t
  integer                  :: N_obs_size
  integer                  :: status
  
  call ESMF_AttributeGet(OBS_State,name="Number Of Observations",&
       value=N_obs_size,rc=status)
  call LIS_verify(status, 'attributeget error in noahmp36_transform_LAI')
  call ESMF_StateGet(OBS_State,"Observation01",obs_lai_field,&
       rc=status)
  call LIS_verify(status,'stateget error in noahmp36_transform_LAI')
  call ESMF_FieldGet(obs_lai_field,localDE=0,farrayPtr=laiobs,rc=status)
  call LIS_verify(status,'fieldget error in noahmp36_transform_LAI')

  ! If using 8th mesh LAI data, convert it from inches to meters.
  ! 16th mesh LAI data are already in meters.
  do t=1,N_obs_size
     if(laiobs(t).ne.LIS_rc%udef) then 
        laiobs(t) = laiobs(t)
     endif
  enddo
#endif
end subroutine noahmp36_transform_albedo
