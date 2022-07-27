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
! !ROUTINE: noahmp401_map_burnarea
! \label{noahmp401_map_burnarea}
!
! !REVISION HISTORY:
! 24 Jul 2022: Sujay Kumar; Initial Specification

!
! !INTERFACE:
subroutine noahmp401_map_burnarea(n,k,OBS_State,LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_constantsMod, only  : LIS_CONST_TKFRZ
  use LIS_logMod,   only  : LIS_logunit, LIS_verify
  use LIS_lsmMod
  use LIS_timeMgrMod
  use noahmp401_lsmMod
  use NOAHMP_TABLES_401, ONLY : DKSAT_TABLE,REFDK_TABLE, REFKDT_TABLE,FSATMX_TABLE
 
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
  type(ESMF_Field)         :: obs_burnarea_field
  integer                  :: t
  integer                  :: status
  integer                  :: obs_state_count
  real                     :: dksat_in, dksat_out, delta
  real                     :: kdt_in, kdt_out
  real                     :: fsatmx_in, fsatmx_out
  real*8                   :: time
  integer                  :: cdoy
  real                     :: cgmt
  integer                  :: SOILTYP           ! soil type index [-]
  real                     :: burn_doy
  real                     :: REFKDT, REFDK, FSATMX
   integer                  :: st_id, en_id
  real, pointer            :: burnareaobs(:)
  character*100,allocatable    :: obs_state_objs(:)

  call LIS_tick(time,cdoy,cgmt,LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
       LIS_rc%hr,LIS_rc%mn,LIS_rc%ss,0.0)
  
  call ESMF_StateGet(OBS_State,itemCount=obs_state_count,rc=status)
  call LIS_verify(status)
  allocate(obs_state_objs(obs_state_count))

  call ESMF_StateGet(OBS_State,itemNameList=obs_state_objs,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(OBS_State,obs_state_objs(1),obs_burnarea_field,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(obs_burnarea_field,localDE=0,farrayPtr=burnareaobs,rc=status)
  call LIS_verify(status)
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     
     call LIS_lsm_DAmapTileSpaceToObsSpace(n,k,t,st_id,en_id)
! Assume here that st_id and en_id are the same and that we are
! working with an model grid finer than the observation grid

     !top soil hydraulic conductivity
     SOILTYP = noahmp401_struc(n)%noahmp401(t)%soiltype        
     DKSAT_IN = DKSAT_TABLE(SOILTYP)   ! MAXSMC (SOILTYP)
     REFKDT = REFKDT_TABLE
     REFDK = REFDK_TABLE
     FSATMX_IN = FSATMX_TABLE
     !from documentation : REFKDT is a tunable parameter that significantly
     !impacts surface infiltration and hence partitioning of total runoff
     !into surface and subsurface runoff. Increasing REFKDT decreases
     !surface runoff
     
     KDT_IN  = REFKDT * DKSAT_IN / REFDK
     dksat_out = dksat_in
     kdt_out = kdt_in
     
     if(st_id.gt.0.and.en_id.gt.0) then 
        if(burnareaobs(st_id).ge.0) then 
           
           burn_doy = burnareaobs(st_id)
           if(burn_doy.le.cdoy) then
              delta = (cdoy-burn_doy)/32
              if(delta.eq.0) then
                 delta = 0.01
              endif
              !slow recovery (in 120days) 2delta provides a 2 month recovery. 
              dksat_out = dksat_in * (1-exp(-delta))
              kdt_out = kdt_in/(1-exp(-delta))
              fsatmx_out = fsatmx_in/(1-exp(-delta))
           endif
           NOAHMP401_struc(n)%noahmp401(t)%param%dksat(1:4) = dksat_out
           !           NOAHMP401_struc(n)%noahmp401(t)%param%kdt = kdt_out
           NOAHMP401_struc(n)%noahmp401(t)%param%fsatmx = fsatmx_out
        endif
     endif
  enddo

  deallocate(obs_state_objs)

  
end subroutine noahmp401_map_burnarea
   
