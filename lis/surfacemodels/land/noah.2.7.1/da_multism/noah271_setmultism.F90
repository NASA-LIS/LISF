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
! !ROUTINE: noah271_setmultism
!  \label{noah271_setmultism}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah271_setmultism(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noah271_lsmMod

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
  real, allocatable          :: soilm1(:)
  real, allocatable          :: soilm2(:)
  real, allocatable          :: soilm3(:)
  real, allocatable          :: soilm4(:)
  integer                :: t
  integer                :: status
  real                   :: delta1,delta2,delta3,delta4

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     delta1 = soilm1(t)-noah271_struc(n)%noah(t)%smc(1)
     delta2 = soilm2(t)-noah271_struc(n)%noah(t)%smc(2)
     delta3 = soilm3(t)-noah271_struc(n)%noah(t)%smc(3)
     delta4 = soilm4(t)-noah271_struc(n)%noah(t)%smc(4)

     if(noah271_struc(n)%noah(t)%sh2o(1)+delta1.gt.0.10 .and.&
          noah271_struc(n)%noah(t)%sh2o(1)+delta1.lt.0.55) then 
        noah271_struc(n)%noah(t)%sh2o(1) = noah271_struc(n)%noah(t)%sh2o(1)+&
             delta1
        noah271_struc(n)%noah(t)%smc(1) = soilm1(t)
     endif
     if(noah271_struc(n)%noah(t)%sh2o(2)+delta2.gt.0.10 .and.&
          noah271_struc(n)%noah(t)%sh2o(2)+delta2.lt.0.55) then 
        noah271_struc(n)%noah(t)%sh2o(2) = noah271_struc(n)%noah(t)%sh2o(2)+&
             soilm2(t)-noah271_struc(n)%noah(t)%smc(2)
        noah271_struc(n)%noah(t)%smc(2) = soilm2(t)
     endif
     if(noah271_struc(n)%noah(t)%sh2o(3)+delta3.gt.0.10 .and.&
          noah271_struc(n)%noah(t)%sh2o(3)+delta3.lt.0.55) then 
        noah271_struc(n)%noah(t)%sh2o(3) = noah271_struc(n)%noah(t)%sh2o(3)+&
             soilm3(t)-noah271_struc(n)%noah(t)%smc(3)
        noah271_struc(n)%noah(t)%smc(3) = soilm3(t)
     endif
     if(noah271_struc(n)%noah(t)%sh2o(4)+delta4.gt.0.10 .and.&
          noah271_struc(n)%noah(t)%sh2o(4)+delta4.lt.0.55) then 
        noah271_struc(n)%noah(t)%sh2o(4) = noah271_struc(n)%noah(t)%sh2o(4)+&
             soilm4(t)-noah271_struc(n)%noah(t)%smc(4)
        noah271_struc(n)%noah(t)%smc(4) = soilm4(t)
     endif
!     if(noah271_struc(n)%noah(t)%sh2o(1).le.0) then 
!        print*, 'problem1 ',t,noah271_struc(n)%noah(t)%sh2o(1), &
!             noah271_struc(n)%noah(t)%smc(1)
!        stop
!     endif
!     if(noah271_struc(n)%noah(t)%sh2o(2).le.0) then 
!        print*, 'problem2 ',t,noah271_struc(n)%noah(t)%sh2o(2), &
!             noah271_struc(n)%noah(t)%smc(2)
!        stop
!     endif
!     if(noah271_struc(n)%noah(t)%sh2o(3).le.0) then 
!        print*, 'problem3 ',t,noah271_struc(n)%noah(t)%sh2o(3), &
!             noah271_struc(n)%noah(t)%smc(3)
!        stop
!     endif
!     if(noah271_struc(n)%noah(t)%sh2o(4).le.0) then 
!        print*, 'problem4 ',t,noah271_struc(n)%noah(t)%sh2o(4), &
!             noah271_struc(n)%noah(t)%smc(4)
!        stop
!     endif
  enddo
end subroutine noah271_setmultism

