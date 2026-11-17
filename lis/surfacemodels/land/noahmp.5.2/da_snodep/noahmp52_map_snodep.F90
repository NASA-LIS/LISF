!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp52_map_snodep
! \label{noahmp52_map_snodep}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; Modified for NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine noahmp52_map_snodep(n,k,OBS_State,LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,   only  : LIS_verify
  use LIS_lsmMod
  use noahmp52_lsmMod

  implicit none
! !ARGUMENTS:
  integer, intent(in)      :: n
  integer, intent(in)      :: k
  type(ESMF_State)         :: OBS_State
  type(ESMF_State)         :: LSM_Incr_State
! !DESCRIPTION:
!
!  This subroutine directly maps the observation state to the corresponding
!  variables in the LSM state for SNODEP data assimilation.
!
!  The arguments are:
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF State for observations \newline
!  \item[LSM\_State] ESMF State for LSM state variables \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: sweIncrField
  type(ESMF_Field)         :: obs_snodep_field
  real, pointer            :: sweincr(:)
  type(ESMF_Field)         :: snodIncrField
  real, pointer            :: snodincr(:)
  real                     :: tmpsneqv
  real, pointer            :: snodepobs(:)
  integer                  :: t
  integer                  :: status
  integer                  :: obs_state_count
  integer                  :: st_id, en_id
  character*100,allocatable    :: obs_state_objs(:)
  real, allocatable            :: noahmp52_swe(:)
  real, allocatable            :: noahmp52_snod(:)
  real, allocatable            :: snod(:)

  allocate(noahmp52_swe(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(noahmp52_snod(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(snod(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"Snowdepth",snodincrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(snodincrField,localDE=0,farrayPtr=snodincr,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(OBS_State,itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))

  call ESMF_StateGet(OBS_State,itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_snodep_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_snodep_field,localDE=0,farrayPtr=snodepobs,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     noahmp52_swe(t)  = NoahMP52_struc(n)%noahmp52(t)%sneqv
     noahmp52_snod(t) = NoahMP52_struc(n)%noahmp52(t)%snowh
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     call LIS_lsm_DAmapTileSpaceToObsSpace(n,k,t,st_id,en_id)

! Assume here that st_id and en_id are the same and that we are
! working with an model grid finer than the observation grid

     if(snodepobs(st_id).ge.0) then
        if(noahmp52_snod(t).gt.1e-6) then
           tmpsneqv = noahmp52_swe(t)/noahmp52_snod(t)
        else
           tmpsneqv = 0.0
        endif

        snod(t) = snodepobs(st_id)

! Based on SNODEP, we manually update SWE
        if(snod(t).lt.2.54E-3) tmpsneqv = 0.0
        if(snod(t).ge.2.54E-3.and.tmpsneqv.lt.0.001) then
           tmpsneqv = 0.20
        endif
        sweincr(t) = tmpsneqv*snod(t) - noahmp52_swe(t)
        snodincr(t) = snod(t) - noahmp52_snod(t)
     else
        sweincr(t) = 0
        snodincr(t) = 0
     endif
  enddo

  deallocate(obs_state_objs)
  deallocate(noahmp52_swe)
  deallocate(noahmp52_snod)
  deallocate(snod)

end subroutine noahmp52_map_snodep

