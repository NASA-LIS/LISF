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
! !ROUTINE: GLS_qcdec_raint
!  \label{GLS_qcdec_raint}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine GLS_qcdec_raint(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify


  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  

! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: rft1Field,rft2Field,rft3Field,stField,siField
  type(ESMF_Field)       :: feasField
  integer, pointer       :: mod_flag(:)
  real, pointer          :: rft1(:),rft2(:),rft3(:),st(:),si(:)
  real                   :: rft1_min, rft1_max
  real                   :: rft2_min, rft2_max
  real                   :: rft3_min, rft3_max
  real                   :: st_min, st_max
  real                   :: si_min, si_max
  integer                :: t
  integer                :: status
  
  n = 1

  call ESMF_StateGet(DEC_State,"RF_THRESHOLD1",rft1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(DEC_State,"RF_THRESHOLD2",rft2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(DEC_State,"RF_THRESHOLD3",rft3Field,rc=status)
  call LIS_verify(status)
!  call ESMF_StateGet(DEC_State,"SUSCEP_INDEX",siField,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(DEC_State,"SUSCEP_INDEX_THRESHOLD",stField,rc=status)
!  call LIS_verify(status)

  call ESMF_FieldGet(rft1Field,localDE=0,farrayPtr=rft1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(rft1Field,'MinRange',rft1_min,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(rft1Field,'MaxRange',rft1_max,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(rft2Field,localDE=0,farrayPtr=rft2,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(rft2Field,'MinRange',rft2_min,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(rft2Field,'MaxRange',rft2_max,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(rft3Field,localDE=0,farrayPtr=rft3,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(rft3Field,'MinRange',rft3_min,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(rft3Field,'MaxRange',rft3_max,rc=status)
  call LIS_verify(status)

!  call ESMF_FieldGet(stField,localDE=0,farrayPtr=st,rc=status)
!  call LIS_verify(status)
!  call ESMF_AttributeGet(stField,'MinRange',st_min,rc=status)
!  call LIS_verify(status)
!  call ESMF_AttributeGet(stField,'MaxRange',st_max,rc=status)
!  call LIS_verify(status)
#if 0 
  call ESMF_FieldGet(siField,localDE=0,farrayPtr=si,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(siField,'MinRange',si_min,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(siField,'MaxRange',si_max,rc=status)
  call LIS_verify(status)

#endif
  call ESMF_StateGet(Feas_State, "Feasibility Flag", feasField, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(feasField,localDE=0,farrayPtr=mod_flag,rc=status)
  call LIS_verify(status)

#if 0 
  n = 1
  mod_flag = 0 

  do t=1,LIS_rc%ntiles(n)
     if(rft1(t).lt.rft1_min) then 
        rft1(t) = rft1_min
        mod_flag(t) = 1
     endif
     if(rft1(t).gt.rft1_max) then 
        rft1(t) = rft1_max
        mod_flag(t) = 1
     endif

     if(rft2(t).lt.rft2_min) then 
        rft2(t) = rft2_min
        mod_flag(t) = 1
     endif
     if(rft2(t).gt.rft2_max) then 
        rft2(t) = rft2_max
        mod_flag(t) = 1
     endif

     if(rft3(t).lt.rft3_min) then 
        rft3(t) = rft3_min
        mod_flag(t) = 1
     endif
     if(rft3(t).gt.rft3_max) then 
        rft3(t) = rft3_max
        mod_flag(t) = 1
     endif

     if(st(t).lt.st_min) then 
        st(t) = st_min
        mod_flag(t) = 1
        print*, 'st min exceeded ',st(t), st_min
        stop
     endif
     if(st(t).gt.st_max) then 
        st(t) = st_max
        mod_flag(t) = 1
        print*, 'st max exceeded ',st(t), st_max
        stop
     endif

     if(si(t).lt.si_min) then 
        si(t) = si_min
        mod_flag(t) = 1
        print*, 'si min exceeded ',si(t), si_min
        stop
     endif
     if(si(t).gt.si_max) then 
        si(t) = si_max
        mod_flag(t) = 1
        print*, 'si max exceeded ',si(t), si_max
        stop
     endif
!constraints
     if(rft2(t).lt.rft1(t).or.&
          rft3(t).lt.rft1(t).or.&
          rft3(t).lt.rft2(t).or.&
          si(t).lt.st(t)) then 
        mod_flag(t) = 1
     endif

  enddo
#endif
  do t=1, LIS_rc%ntiles(n)
     !constraints
     if(rft2(t).lt.rft1(t).or.&
          rft3(t).lt.rft1(t).or.&
          rft3(t).lt.rft2(t)) then 
        mod_flag(t) = 1
     endif
     
  enddo

end subroutine GLS_qcdec_raint

