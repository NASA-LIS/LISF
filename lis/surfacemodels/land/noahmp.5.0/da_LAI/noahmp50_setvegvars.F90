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
! !ROUTINE: noahmp50_setvegvars
!  \label{noahmp50_setvegvars}
!
! !REVISION HISTORY:
! 13 Feb 2020: Sujay Kumar; Initial Specification
! 
! Apply the update if it met the update conditions
! Update conditions: 
!                  1- Prior SM(sh2o) + increment > MIN_THRESHOLD 
!                  2- Prior SM(sh2o) + increment < sm_threshold
! There are 3 cases 
! 1- If all the ensemble members met the update conditions --> apply the update
! 2- If more than 50% of the ensemble members met the update condition --> 
!    apply the update for that members and set the other member to the mean 
!    value of the ensemble (i.e. mean of the members that met the conditions)
! 3- If less then 50% of the ensemble members met the update conditions --> 
!    adjust the states    


! !INTERFACE:
subroutine noahmp50_setvegvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp50_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP
  real, parameter        :: MIN_THRESHOLD = 0.02 
  type(ESMF_Field)       :: laiField,sm1Field
  real                   :: MAX_threshold
  real                   :: sm_threshold
  real                   :: delta1
  integer                :: SOILTYP           ! soil type index [-]
  integer                :: t,i,gid,m,t_unpert
  integer                :: status
  real, pointer          :: lai(:),soilm1(:)
  real                   :: lfmass
  logical                :: flag_tmp(LIS_rc%nensem(n))
  logical                :: update_flag(LIS_rc%ngrid(n))
  logical                :: ens_flag(LIS_rc%nensem(n))
  logical                :: update_flag_tile(LIS_rc%npatch(n,LIS_rc%lsm_index))
  logical                :: update_flag_ens(LIS_rc%ngrid(n))
  logical                :: update_flag_new(LIS_rc%ngrid(n))
  real                   :: tmp1(LIS_rc%nensem(n))
  integer                :: pcount
  logical                :: bounds_violation
  real                   :: MinEnsSM1 ,MaxEnsSM1
  real                   :: tmpval
  integer                :: nIter
  real                   :: smc_tmp 

  call ESMF_StateGet(LSM_State,"LAI",laiField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     if(NoahMP50_struc(n)%noahmp50(t)%param%sla.ne.0) then 
        NoahMP50_struc(n)%noahmp50(t)%lai = lai(t)
        lfmass = lai(t)*1000.0/(NoahMP50_struc(n)%noahmp50(t)%param%sla)
        NoahMP50_struc(n)%noahmp50(t)%lfmass = lfmass
     endif
  enddo

end subroutine noahmp50_setvegvars

