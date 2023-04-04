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
! !ROUTINE: noahmp401_updatetws
!  \label{noahmp401_updatetws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li; created for Noah-MP4.0.1 
!
! !INTERFACE:
subroutine noahmp401_updatetws(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use noahmp401_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!  
!  This routine assigns the soil moisture prognostic variables to noah's
!  model space. 
! 
!EOP

  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  type(ESMF_Field)       :: sm1IncrField
  type(ESMF_Field)       :: sm2IncrField
  type(ESMF_Field)       :: sm3IncrField
  type(ESMF_Field)       :: sm4IncrField
  type(ESMF_Field)       :: sweField, sweIncrField
  type(ESMF_Field)       :: snodField, snodIncrField

  !Wanshu
  type(ESMF_Field)     :: gwField
  type(ESMF_Field)     :: gwIncrField
  real, pointer        :: gws(:)
  real, pointer        :: gwsIncr(:)  

  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: soilmIncr1(:)
  real, pointer          :: soilmIncr2(:)
  real, pointer          :: soilmIncr3(:)
  real, pointer          :: soilmIncr4(:)
  real, pointer          :: swe(:), sweincr(:)
  real, pointer          :: snod(:), snodincr(:)
  integer                :: t,i,m,gid
  integer                :: status
  real                   :: swetmp, snodtmp,sndens
  logical                :: update_flag(LIS_rc%ngrid(n))
  real                   :: perc_violation(LIS_rc%ngrid(n))

  real                   :: snodmean(LIS_rc%ngrid(n))
  integer                :: nsnodmean(LIS_rc%ngrid(n))
  

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in noahmp401_updatetws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 2 failed in noahmp401_updatetws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 3 failed in noahmp401_updatetws")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 4 failed in noahmp401_updatetws")

  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp401_updatetws")
 
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)



  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 1",sm1IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 1 failed in noahmp401_updatetws")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 2",sm2IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 2 failed in noahmp401_updatetws")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 3",sm3IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 3 failed in noahmp401_updatetws")
  call ESMF_StateGet(LSM_Incr_State,"Soil Moisture Layer 4",sm4IncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet: Soil Moisture Layer 4 failed in noahmp401_updatetws")

  call ESMF_StateGet(LSM_Incr_State, "Groundwater Storage",gwIncrField,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp401_updatetws")
  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)


  !-------------------
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp401_updatetws")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in noahmp401_updatetws")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in noahmp401_updatetws")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in noahmp401_updatetws")
  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Groundwater Storage failed in noahmp401_updatetws")
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)



  call ESMF_FieldGet(sm1IncrField,localDE=0,farrayPtr=soilmIncr1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noahmp401_updatetws")
  call ESMF_FieldGet(sm2IncrField,localDE=0,farrayPtr=soilmIncr2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in noahmp401_updatetws")
  call ESMF_FieldGet(sm3IncrField,localDE=0,farrayPtr=soilmIncr3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in noahmp401_updatetws")
  call ESMF_FieldGet(sm4IncrField,localDE=0,farrayPtr=soilmIncr4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in noahmp401_updatetws")
  call ESMF_FieldGet(gwIncrField, localDE=0, farrayPtr=gwsIncr,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Groundwater Storage failed in noahmp401_updatetws")
  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)



!  write(*,fmt='(I4.4, 1x, I2.2, 1x, I2.2, 1x, I2.2, 1x, I2.2,1x,2E14.6)') &
!       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr,LIS_rc%mn,&
!       sum(swe(991:1000))/10.0, sum(sweincr(991:1000))/10.0
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     soilm1(t) = soilm1(t) + soilmIncr1(t)
     soilm2(t) = soilm2(t) + soilmIncr2(t)
     soilm3(t) = soilm3(t) + soilmIncr3(t)
     soilm4(t) = soilm4(t) + soilmIncr4(t)
     gws(t)    = gws(t)    + gwsIncr(t)
     swe(t) = swe(t)  + sweIncr(t)
  enddo

#if 0   
  
  update_flag    = .true.
  perc_violation = 0.0
  snodmean       = 0.0
  nsnodmean      = 0

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)     
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)

     swetmp = swe(t) + sweincr(t)

     if((swetmp.lt.0)) then
        update_flag(gid) = .false.
        perc_violation(gid) = perc_violation(gid) +1
     endif

  enddo

  do gid=1,LIS_rc%ngrid(n)
     perc_violation(gid) = perc_violation(gid) / real(LIS_rc%nensem(n))
  enddo

! For ensembles that are unphysical, compute the ensemble average after excluding them. This
! is done only if the majority of the ensemble members are good (>80%)



! If the update is unphysical, simply set to the average of
! the good ensemble members. If all else fails, do not update.

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     gid = LIS_domain(n)%gindex(&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%col,&
          LIS_surface(n,LIS_rc%lsm_index)%tile(t)%row)



     swetmp  = swe(t) + sweincr(t)

!Use the model's snow density from the previous timestep
     sndens = 0.0
     if(noahmp401_struc(n)%noahmp401(t)%snowh.gt.0) then
       sndens = noahmp401_struc(n)%noahmp401(t)%sneqv/noahmp401_struc(n)%noahmp401(t)%snowh
     endif

     if(update_flag(gid)) then
        snod(t) = snodtmp
        swe(t) = swetmp
     elseif(perc_violation(gid).lt.0.2) then
       if(snodtmp.lt.0.0) then  ! average of the good ensemble members
          snod(t) = snodmean(gid)
          swe(t) =  snodmean(gid)*sndens
       else
          snod(t) = snodtmp
          swe(t) = swetmp
       endif
     else            ! do not update
       snod(t) = noahmp401_struc(n)%noahmp401(t)%snowh
       swe(t)  = noahmp401_struc(n)%noahmp401(t)%sneqv
     end if
     
  enddo
#endif
  
!  write(*,fmt='(I4.4, 1x, I2.2, 1x, I2.2, 1x, I2.2, 1x, I2.2,1x,2E14.6)') &
!       LIS_rc%yr, LIS_rc%mo, LIS_rc%da, LIS_rc%hr,LIS_rc%mn,&
!       sum(swe(991:1000))/10.0, sum(sweincr(991:1000))/10.0


end subroutine noahmp401_updatetws
