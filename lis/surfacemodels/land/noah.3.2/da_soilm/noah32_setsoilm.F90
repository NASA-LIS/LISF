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
! !ROUTINE: noah32_setsoilm
!  \label{noah32_setsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah32_setsoilm(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify, LIS_logunit
  use noah32_lsmMod

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

  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  integer                :: t,i,m
  integer                :: status
  real                   :: delta1,delta2,delta3,delta4
  real                   :: sm_threshold_lo
  real                   :: sm_threshold_hi
  integer                :: cutoff_exceeded

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 1 failed in noah32_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 2 failed in noah32_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 3 failed in noah32_setsoilm")
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateSet: Soil Moisture Layer 4 failed in noah32_setsoilm")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 1 failed in noah32_setsoilm")
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 2 failed in noah32_setsoilm")
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 3 failed in noah32_setsoilm")
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet: Soil Moisture Layer 4 failed in noah32_setsoilm")

  do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)/LIS_rc%nensem(n)
     cutoff_exceeded=0
     do m=1,LIS_rc%nensem(n)
        t=(i-1)*LIS_rc%nensem(n)+m

        sm_threshold_lo = noah32_struc(n)%noah(t)%smcdry
        sm_threshold_hi = noah32_struc(n)%noah(t)%smcmax  !- 0.05

        delta1 = soilm1(t)-noah32_struc(n)%noah(t)%smc(1)
        delta2 = soilm2(t)-noah32_struc(n)%noah(t)%smc(2)
        delta3 = soilm3(t)-noah32_struc(n)%noah(t)%smc(3)
        delta4 = soilm4(t)-noah32_struc(n)%noah(t)%smc(4)

        ! test to see if cutoffs are exceeded anywhere
        if(noah32_struc(n)%noah(t)%sh2o(1)+delta1.lt.sm_threshold_lo .or.&
             noah32_struc(n)%noah(t)%sh2o(1)+delta1.gt.&
             sm_threshold_hi) then 
           cutoff_exceeded=1
        endif

        ! test to see if cutoffs are exceeded anywhere
        if(noah32_struc(n)%noah(t)%sh2o(2)+delta2.lt.sm_threshold_lo .or.&
             noah32_struc(n)%noah(t)%sh2o(2)+delta2.gt.&
             sm_threshold_hi) then 
           cutoff_exceeded=1
        endif

        ! test to see if cutoffs are exceeded anywhere
        if(noah32_struc(n)%noah(t)%sh2o(3)+delta3.lt.sm_threshold_lo .or.&
             noah32_struc(n)%noah(t)%sh2o(3)+delta3.gt.&
             sm_threshold_hi) then 
           cutoff_exceeded=1
        endif

        ! test to see if cutoffs are exceeded anywhere
        if(noah32_struc(n)%noah(t)%sh2o(4)+delta4.lt.sm_threshold_lo .or.&
             noah32_struc(n)%noah(t)%sh2o(4)+delta4.gt.&
             sm_threshold_hi) then 
           cutoff_exceeded=1
        endif
     enddo
     if(cutoff_exceeded.eq.0) then
        do m=1,LIS_rc%nensem(n)
           t=(i-1)*LIS_rc%nensem(n)+m
           sm_threshold_lo = noah32_struc(n)%noah(t)%smcdry
           sm_threshold_hi = noah32_struc(n)%noah(t)%smcmax  !- 0.05
          
           delta1 = soilm1(t)-noah32_struc(n)%noah(t)%smc(1)
           delta2 = soilm2(t)-noah32_struc(n)%noah(t)%smc(2)
           delta3 = soilm3(t)-noah32_struc(n)%noah(t)%smc(3)
           delta4 = soilm4(t)-noah32_struc(n)%noah(t)%smc(4)

           noah32_struc(n)%noah(t)%sh2o(1) = noah32_struc(n)%noah(t)%sh2o(1)+ delta1
           noah32_struc(n)%noah(t)%sh2o(2) = noah32_struc(n)%noah(t)%sh2o(2)+ delta2
           noah32_struc(n)%noah(t)%sh2o(3) = noah32_struc(n)%noah(t)%sh2o(3)+ delta3
           noah32_struc(n)%noah(t)%sh2o(4) = noah32_struc(n)%noah(t)%sh2o(4)+ delta4
           noah32_struc(n)%noah(t)%smc(1) = soilm1(t)
           noah32_struc(n)%noah(t)%smc(2) = soilm2(t)
           noah32_struc(n)%noah(t)%smc(3) = soilm3(t)
           noah32_struc(n)%noah(t)%smc(4) = soilm4(t)
        enddo
     else
        !do not apply increments
     endif
  enddo

end subroutine noah32_setsoilm

