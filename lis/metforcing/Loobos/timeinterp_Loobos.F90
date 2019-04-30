!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: timeinterp_Loobos
! \label{timeinterp_Loobos}
! 
!
! !REVISION HISTORY:
! 05 Oct 2010: David Mocko, Updated for Loobos test case
!
! !INTERFACE:
subroutine timeinterp_Loobos(n,findex)
  ! !USES:
  use ESMF
  use LIS_logMod, only           : LIS_logunit,LIS_verify,LIS_endrun
  use LIS_coreMod, only          : LIS_rc,LIS_domain
  use LIS_metforcingMod, only   : LIS_FORC_Base_State, LIS_forc
  use Loobos_forcingMod, only : Loobos_struc
  use LIS_FORC_AttributesMod

  implicit none
  ! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

  ! !DESCRIPTION:
  ! Temporally interpolates the Loobos forcing data to
  ! the model timestep.  All variables except precipitation
  ! is linearly interpolated. 
  ! 
  !EOP

  real :: wt1,wt2
  integer :: t
  integer :: index1,tid
  integer :: status
  type(ESMF_Field)   :: tairField,qairField,uField,vField,swdownField
  type(ESMF_Field)   :: lwdownField,psurfField,rainField,snowField
  real,pointer       :: tair(:),qair(:),uwind(:),vwind(:)
  real,pointer       :: swdown(:),lwdown(:),psurf(:),rain(:),snow(:)
  real, parameter    :: eps = 0.622
  real               :: svp,qs,E

  !      write(LIS_logunit,*) 'starting timeinterp_Loobos'

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Tair%varname(1),   &
       tairField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Qair%varname(1),   &
       qairField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Swdown%varname(1), &
       swdownField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1), &
       lwdownField,rc=status)
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
       rainField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Rainf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_Snowf%varname(1), &
       snowField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Snowf in the forcing variables list')

  call ESMF_FieldGet(tairField,localDE=0,farrayPtr=tair,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(qairField,localDE=0,farrayPtr=qair,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swdownField,localDE=0,farrayPtr=swdown,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdownField,localDE=0,farrayPtr=lwdown,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(rainField,localDE=0,farrayPtr=rain,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(snowField,localDE=0,farrayPtr=snow,rc=status)
  call LIS_verify(status)

  !      write(LIS_logunit,*) 'Btime1: ',Loobos_struc(n)%Loobostime1
  !      write(LIS_logunit,*) 'Btime2: ',Loobos_struc(n)%Loobostime2
  !      write(LIS_logunit,*) 'realtime: ',LIS_rc%time
  wt1 = (Loobos_struc(n)%Loobostime2-LIS_rc%time) /        &
       (Loobos_struc(n)%Loobostime2-                      &
       Loobos_struc(n)%Loobostime1)
  wt2 = 1.0 - wt1
  !      write(LIS_logunit,*) wt1,wt2

  do t = 1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
        psurf(t) = wt1 * Loobos_struc(n)%metdata1(1,index1) +                     &
             wt2 * Loobos_struc(n)%metdata2(1,index1)
        tair(t) = wt1 * Loobos_struc(n)%metdata1(2,index1) +                      &
             wt2 * Loobos_struc(n)%metdata2(2,index1)
        qair(t) = wt1 * Loobos_struc(n)%metdata1(3,index1) +                      &
             wt2 * Loobos_struc(n)%metdata2(3,index1)
        uwind(t) = wt1 * Loobos_struc(n)%metdata1(4,index1) +                     &
             wt2 * Loobos_struc(n)%metdata2(4,index1)
        vwind(t) = 0.0
        swdown(t) = wt1 * Loobos_struc(n)%metdata1(6,index1) +                    &
             wt2 * Loobos_struc(n)%metdata2(6,index1)
        lwdown(t) = wt1 * Loobos_struc(n)%metdata1(7,index1) +                    &
             wt2 * Loobos_struc(n)%metdata2(7,index1)
        rain(t) = Loobos_struc(n)%metdata1(8,index1)
        snow(t) = Loobos_struc(n)%metdata1(9,index1)
  enddo
end subroutine timeinterp_Loobos


