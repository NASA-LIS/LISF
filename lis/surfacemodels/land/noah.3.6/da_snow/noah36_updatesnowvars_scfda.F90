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
! !ROUTINE: noah36_updatesnowvars_scfda
! \label{noah36_updatesnowvars_scfda}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! may 30 2014: Yuqiong Liu; adpated for SCF EnKF assimilation
! !INTERFACE:
subroutine noah36_updatesnowvars_scfda(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_surface
  use noah36_lsmMod
  use LIS_logMod,   only : LIS_logunit, LIS_verify
  !use enkf_Mod,  only: enkf_struc
  use LIS_DAobservationsMod, only : LIS_OBS_State
  use ANSASCFsnow_Mod, only: ANSASCFsnow_struc

  implicit none

! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Returns the snow related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \item[LSM\_Incr\_State] ESMF State container for LSM state increments \newline
!  \end{description}
!
!EOP


  external noah36_getscfpred
  external noah36_transform_snow
  external noah36_map_snow_DI

  type(ESMF_Field)       :: sweField, sweIncrField
  type(ESMF_Field)       :: snodField, snodIncrField

  integer                :: t, m, gid
  integer                :: status
  real, pointer          :: swe(:), sweincr(:)
  real, pointer          :: snod(:), snodincr(:)
  integer                :: i,k1,n0,n1,n2
  real, pointer          :: sweincr1(:), snodincr1(:)
  integer, allocatable       :: idx1(:)
  real                   :: obs_pred(LIS_rc%ngrid(n),LIS_rc%nensem(n))

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Snowdepth",snodIncrField,rc=status)
  call LIS_verify(status)

 
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodIncrField,localDE=0,farrayPtr=snodincr,rc=status)
  call LIS_verify(status)

  allocate(sweincr1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(snodincr1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     sweincr1(t) = sweincr(t)
     snodincr1(t) = snodincr(t)
  enddo

  if (ANSASCFsnow_struc(n)%useEnKFwithDI .eq. 1) then
     
     allocate(idx1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
     idx1 = 0
     call noah36_getscfpred(n,obs_pred)

     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index),LIS_rc%nensem(n)
        n1=0
        n2=0
        do m=1,LIS_rc%nensem(n)
           t = i+m-1
           gid = LIS_surface(n,1)%tile(t)%index
           if (obs_pred(gid,m) .le. 1) n1=n1+1
           if (obs_pred(gid,m) .ge. 99) n2=n2+1
        enddo
        if (n1.eq.LIS_rc%nensem(n) .or. n2.eq.LIS_rc%nensem(n)) then
           idx1(i:i+LIS_rc%nensem(n)-1)=1
        endif
     enddo
     
     n0=0
     do i=1, LIS_rc%ndas 
        if (LIS_rc%daset(i).eq."ANSA SCF" .or. LIS_rc%daset(i).eq."MODIS SCF") then
           k1=i
           n0=n0+1
        endif
     enddo
     if (n0.ne.1) then
        call LIS_verify(1,'Must assimilate one and only one SCF dataset')
     endif

     call noah36_transform_snow(n, LIS_OBS_State(n,k1))
     call noah36_map_snow_DI(n,LIS_OBS_State(n,k1),LSM_Incr_State)

     call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
     call LIS_verify(status, 'failed to retrieve SWE state in noah36_updatesnowvars')
     call ESMF_StateGet(LSM_Incr_State,"Snowdepth",snodIncrField,rc=status)
     call LIS_verify(status, 'failed to retrieve Snowdepth state in noah36_updatesnowvars')  
     call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
     call LIS_verify(status, 'failed to retrieve SWE field in noah36_updatesnowvars')
     call ESMF_FieldGet(snodIncrField,localDE=0,farrayPtr=snodincr,rc=status)
     call LIS_verify(status,'failed to retrieve Snowdepth field in noah36_updatesnowvars')
 
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        if (idx1(t).eq.1)   then
           sweincr1(t) = sweincr(t)
           snodincr(t) = snodincr(t)
        endif   
     enddo
  endif 

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)       
     swe(t) = swe(t) + sweincr1(t)
     snod(t) = snod(t) + snodincr1(t)
  enddo
end subroutine noah36_updatesnowvars_scfda

