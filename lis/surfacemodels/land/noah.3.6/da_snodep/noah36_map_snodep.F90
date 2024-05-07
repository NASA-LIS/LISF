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
! !ROUTINE: noah36_map_snodep
! \label{noah36_map_snodep}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!  21 Jul 2011: James Geiger; Modified for Noah 3.2
!
! !INTERFACE:
subroutine noah36_map_snodep(n,k,OBS_State,LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use LIS_constantsMod, only  : LIS_CONST_TKFRZ
  use LIS_logMod,   only  : LIS_logunit, LIS_verify
  use LIS_lsmMod
  use noah36_lsmMod

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
  real, allocatable            :: noah36_swe(:)
  real, allocatable            :: noah36_snod(:)
  real, allocatable            :: snod(:)

  allocate(noah36_swe(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(noah36_snod(LIS_rc%npatch(n,LIS_rc%lsm_index)))
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
     noah36_swe(t)  = noah36_struc(n)%noah(t)%sneqv
     noah36_snod(t) = noah36_struc(n)%noah(t)%snowh
  enddo

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     call LIS_lsm_DAmapTileSpaceToObsSpace(n,k,t,st_id,en_id)

! Assume here that st_id and en_id are the same and that we are
! working with an model grid finer than the observation grid

     if(snodepobs(st_id).ge.0) then 
        if(noah36_snod(t).gt.1e-6) then 
           tmpsneqv = noah36_swe(t)/noah36_snod(t)
        else
           tmpsneqv = 0.0
        endif

        snod(t) = snodepobs(st_id)

! Based on SNODEP, we manually update SWE
        if(snod(t).lt.2.54E-3) tmpsneqv = 0.0
        if(snod(t).ge.2.54E-3.and.tmpsneqv.lt.0.001) then 
           tmpsneqv = 0.20
        endif        
        sweincr(t) = tmpsneqv*snod(t) - noah36_swe(t)
        snodincr(t) = snod(t) - noah36_snod(t)
!        if(sweincr(t).gt.20) then 
!           print*, 'sweincr ',t, sweincr(t), snod(t), tmpsneqv, noah36_swe(t)
!        endif
     else
        sweincr(t) = 0 
        snodincr(t) = 0 
     endif
  enddo
!  stop
  deallocate(obs_state_objs)
  deallocate(noah36_swe)
  deallocate(noah36_snod)
  deallocate(snod)

end subroutine noah36_map_snodep
   
