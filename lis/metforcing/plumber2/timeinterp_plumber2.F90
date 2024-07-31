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
!
! !ROUTINE: timeinterp_plumber2
! \label{timeinterp_plumber2}
! 
!
! !REVISION HISTORY:
! 16 Sep 2021: Mark Beauharnois, based on Bondville timeinterp
!
! !INTERFACE:
subroutine timeinterp_plumber2(n,findex)
  ! !USES:
  use ESMF
  use LIS_logMod, only           : LIS_logunit,LIS_verify,LIS_endrun
  use LIS_coreMod, only          : LIS_rc,LIS_domain
  use LIS_metforcingMod, only    : LIS_FORC_Base_State, LIS_forc
  use plumber2_forcingMod, only  : plumber2_struc
  use LIS_FORC_AttributesMod

  implicit none
  ! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

  ! !DESCRIPTION:
  ! Temporally interpolates the PLUMBER2 forcing data to
  ! the model timestep.  All variables except precipitation
  ! are linearly interpolated. 
  !
  ! Index order in 'metdata'
    !  1. prec rain (total precip)
    !  2. psurf
    !  3. qair
    !  4. tair
    !  5. swdown
    !  6. lwdown
    !  7. wind u
    !  8. wind v
    !  9. lai    currently not used
  ! 
  !EOP

  real :: wt1,wt2
  integer :: t
  integer :: index1,tid
  integer :: status
  type(ESMF_Field)   :: lwdField, pcpField, psurfField, qairField
  type(ESMF_Field)   :: swdField, tairField, uField, vField

  real,pointer       :: tair(:),qair(:),uwind(:),vwind(:)
  real,pointer       :: swdown(:),lwdown(:),psurf(:),pcp(:)
! real,pointer       :: lai(:)
  integer, save      :: count

  !      write(LIS_logunit,*) 'starting timeinterp_plumber2'

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),   &
       tairField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),   &
       qairField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1), &
       swdField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1), &
       lwdField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable LWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_E%varname(1), &
       uField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Wind_N%varname(1), &
       vField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Wind_N in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Psurf%varname(1),  &
       psurfField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Psurf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Rainf%varname(1),  &
       pcpField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Rainf in the forcing variables list')

! call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LAI%varname(1), &
!      laiField,rc=status)
! call LIS_verify(status,                                          &
!      'Error: Enable LAI in the forcing variables list')


  call ESMF_FieldGet(tairField,localDE=0,farrayPtr=tair,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(qairField,localDE=0,farrayPtr=qair,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swdown,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwdown,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

! call ESMF_FieldGet(laiField,localDE=0,farrayPtr=lai,rc=status)
! call LIS_verify(status)

  !write(LIS_logunit,*) 'Ptime1: ',plumber2_struc(n)%plumber2time1
  !write(LIS_logunit,*) 'Ptime2: ',plumber2_struc(n)%plumber2time2
  !write(LIS_logunit,*) 'realtime: ',LIS_rc%time
  wt1 = (plumber2_struc(n)%plumber2time2-LIS_rc%time) /        &
        (plumber2_struc(n)%plumber2time2-                      &
         plumber2_struc(n)%plumber2time1)
  wt2 = 1.0 - wt1
  !      write(LIS_logunit,*) 'timeinter_plumber2 wts: ',wt1,wt2

  if (LIS_rc%startcode.eq."restart") then
     count = count + 1
     if (wt1.gt.1.5) then
        wt1 = 1.0
        wt2 = 0.0
     endif
     if (count.eq.2) then
        wt1 = 0.5
        wt2 = 0.5
     endif
  endif

    !  1. prec rain (total precip)
    !  2. psurf
    !  3. qair
    !  4. tair
    !  5. swdown
    !  6. lwdown
    !  7. wind u
    !  8. wind v
    !  9. lai    currently not used

  do t = 1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index

     !! -- pcp (1)
     ! Precip rate is not interpolated, but carried forward
     if (plumber2_struc(n)%metdata1(1,t).ne.plumber2_struc(n)%undef) then
        pcp(t) = plumber2_struc(n)%metdata1(1,index1)
     endif

     !! -- psurf (2)
     if ((plumber2_struc(n)%metdata1(2,t).ne.plumber2_struc(n)%undef).and.        &
          (plumber2_struc(n)%metdata2(2,t).ne.plumber2_struc(n)%undef)) then
        psurf(t) = wt1 * plumber2_struc(n)%metdata1(2,index1) +                      &
             wt2 * plumber2_struc(n)%metdata2(2,index1)
     endif

     !! -- qair (3)
     if ((plumber2_struc(n)%metdata1(3,t).ne.plumber2_struc(n)%undef).and.        &
          (plumber2_struc(n)%metdata2(3,t).ne.plumber2_struc(n)%undef)) then
        qair(t) = wt1 * plumber2_struc(n)%metdata1(3,index1) +                      &
             wt2 * plumber2_struc(n)%metdata2(3,index1)
     endif

     !! -- tair (4)
     if ((plumber2_struc(n)%metdata1(4,t).ne.plumber2_struc(n)%undef).and.        &
          (plumber2_struc(n)%metdata2(4,t).ne.plumber2_struc(n)%undef)) then
        tair(t) = wt1 * plumber2_struc(n)%metdata1(4,index1) +                    &
             wt2 * plumber2_struc(n)%metdata2(4,index1)
     endif

     !! -- swdown (5)
     if ((plumber2_struc(n)%metdata1(5,t).ne.plumber2_struc(n)%undef).and.        &
          (plumber2_struc(n)%metdata2(5,t).ne.plumber2_struc(n)%undef)) then
        swdown(t) = wt1 * plumber2_struc(n)%metdata1(5,index1) +                    &
             wt2 * plumber2_struc(n)%metdata2(5,index1)
     endif
     
     !! -- lwdown (6)
     if ((plumber2_struc(n)%metdata1(6,t).ne.plumber2_struc(n)%undef).and.        &
          (plumber2_struc(n)%metdata2(6,t).ne.plumber2_struc(n)%undef)) then
        lwdown(t) = wt1 * plumber2_struc(n)%metdata1(6,index1) +                     &
             wt2 * plumber2_struc(n)%metdata2(6,index1)
     endif

     !! -- wind u (7)
     if ((plumber2_struc(n)%metdata1(7,t).ne.plumber2_struc(n)%undef).and.        &
          (plumber2_struc(n)%metdata2(7,t).ne.plumber2_struc(n)%undef)) then
        uwind(t) = wt1 * plumber2_struc(n)%metdata1(7,index1) +                     &
             wt2 * plumber2_struc(n)%metdata2(7,index1)
     endif

     !! -- wind v (8)
     if ((plumber2_struc(n)%metdata1(8,t).ne.plumber2_struc(n)%undef).and.        &
          (plumber2_struc(n)%metdata2(8,t).ne.plumber2_struc(n)%undef)) then
        vwind(t) = wt1 * plumber2_struc(n)%metdata1(8,index1) +                     &
             wt2 * plumber2_struc(n)%metdata2(8,index1)
     endif

     !! -- lai (9)
!    if ((plumber2_struc(n)%metdata1(9,t).ne.plumber2_struc(n)%undef).and.        &
!         (plumber2_struc(n)%metdata2(9,t).ne.plumber2_struc(n)%undef)) then
!       lai(t) = wt1 * plumber2_struc(n)%metdata1(9,index1) +                    &
!            wt2 * plumber2_struc(n)%metdata2(9,index1)
!       !            write(LIS_logunit,*) plumber2_struc(n)%metdata1(9,index1),plumber2_struc(n)%metdata2(9,index1)
!    endif

  enddo

end subroutine timeinterp_plumber2
