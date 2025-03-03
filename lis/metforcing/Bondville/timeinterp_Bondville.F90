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
! !ROUTINE: timeinterp_Bondville
! \label{timeinterp_Bondville}
! 
!
! !REVISION HISTORY:
! 05 Oct 2010: David Mocko, Updated for Bondville test case
!
! !INTERFACE:
subroutine timeinterp_Bondville(n,findex)
  ! !USES:
  use ESMF
  use LIS_logMod, only           : LIS_logunit,LIS_verify,LIS_endrun
  use LIS_coreMod, only          : LIS_rc,LIS_domain
  use LIS_metforcingMod, only   : LIS_FORC_Base_State, LIS_forc
  use Bondville_forcingMod, only : Bondville_struc
  use LIS_FORC_AttributesMod

  implicit none
  ! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex

  ! !DESCRIPTION:
  ! Temporally interpolates the Bondville forcing data to
  ! the model timestep.  All variables except precipitation
  ! is linearly interpolated. 
  ! 
  !EOP

  real :: wt1,wt2
  integer :: t
  integer :: index1,tid
  integer :: status
  type(ESMF_Field)   :: tairField,qairField,uField,vField,swdownField
  type(ESMF_Field)   :: lwdownField,psurfField,pcpField,cpcpField
  real,pointer       :: tair(:),qair(:),uwind(:),vwind(:)
  real,pointer       :: swdown(:),lwdown(:),psurf(:),pcp(:),cpcp(:)
  real, parameter    :: eps = 0.622
  real               :: svp,qs,E,EsFuncT
  integer, save      :: count

  !      write(LIS_logunit,*) 'starting timeinterp_Bondville'

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
       pcpField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable Rainf in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_CRainf%varname(1), &
       cpcpField,rc=status)
  call LIS_verify(status,                                          &
       'Error: Enable CRainf in the forcing variables list')

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

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  call LIS_verify(status)

!  write(LIS_logunit,*) 'Btime1: ',Bondville_struc(n)%Bondvilletime1
!  write(LIS_logunit,*) 'Btime2: ',Bondville_struc(n)%Bondvilletime2
!  write(LIS_logunit,*) 'realtime: ',LIS_rc%time
  wt1 = (Bondville_struc(n)%Bondvilletime2-LIS_rc%time) /        &
        (Bondville_struc(n)%Bondvilletime2-                      &
         Bondville_struc(n)%Bondvilletime1)
  wt2 = 1.0 - wt1
  !      write(LIS_logunit,*) wt1,wt2

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

  do t = 1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     if ((Bondville_struc(n)%metdata1(1,t).ne.Bondville_struc(n)%undef).and.        &
          (Bondville_struc(n)%metdata2(1,t).ne.Bondville_struc(n)%undef)) then
        psurf(t) = wt1 * Bondville_struc(n)%metdata1(1,index1) +                     &
             wt2 * Bondville_struc(n)%metdata2(1,index1)
        !            write(LIS_logunit,*) Bondville_struc(n)%metdata1(1,index1),Bondville_struc(n)%metdata2(1,index1)
     endif
     if ((Bondville_struc(n)%metdata1(2,t).ne.Bondville_struc(n)%undef).and.        &
          (Bondville_struc(n)%metdata2(2,t).ne.Bondville_struc(n)%undef)) then
        tair(t) = wt1 * Bondville_struc(n)%metdata1(2,index1) +                      &
             wt2 * Bondville_struc(n)%metdata2(2,index1)
        !            write(LIS_logunit,*) Bondville_struc(n)%metdata1(2,index1),Bondville_struc(n)%metdata2(2,index1)
     endif
     if ((Bondville_struc(n)%metdata1(3,t).ne.Bondville_struc(n)%undef).and.        &
          (Bondville_struc(n)%metdata2(3,t).ne.Bondville_struc(n)%undef)) then
        qair(t) = wt1 * Bondville_struc(n)%metdata1(3,index1) +                      &
             wt2 * Bondville_struc(n)%metdata2(3,index1)
        if (Bondville_struc(n)%MP.eq.0) then
        ! Conversion of psurf and specific humidity are from the
        ! module_ascii_io.F code included in the Noah3.1 release package.
           psurf(t) = psurf(t) * 100.0
        ! Convert RH [ % ] to Specific Humidity [ kg kg{-1} ]
        ! This computation from NCEP's Noah v2.7.1 driver.
           qair(t) = qair(t) * 1.E-2
           svp = EsFuncT(tair(t))
!           print *,'SFCTMP: ',tair(t)
!           print *,'SVP: ',SVP
!           print *,'EPS: ',EPS
!           print *,'SFCPRS: ',psurf(t)
           qs = eps * svp / (psurf(t) - (1.-eps) * svp)
!           print *,'QS: ',QS
!           print *,'RHF: ',qair(t)
           E = (psurf(t)*svp*qair(t))/(psurf(t) - svp*(1. - qair(t)))
!           print *,'E: ',E
           qair(t) = (eps*E)/(psurf(t)-(1.0-eps)*E)
!           print *,'SPECHUMD: ',qair(t)
           if (qair(t).lt.0.1E-5) qair(t) = 0.1E-5
           if (qair(t).ge.qs)     qair(t) = qs*0.99
        else
        ! Conversion of psurf and specific humidity are from the
        ! HRLDAS_forcing/run/examples/single_point/create_point_data.f90
        ! code included in the HRLDAS Noah-MP.4.0.1 release package.
           psurf(t) = psurf(t) * 100.0
           svp = 611.2*exp(17.67*(tair(t)-273.15)/(tair(t)-29.65)) ! [Pa]
           e   = qair(t)/100.0 * svp                               ! [Pa]
           qair(t) = (0.622*e)/(psurf(t)-(1.0-0.622)*e)            ! now it is specific humidity
        endif
        !            write(LIS_logunit,*) Bondville_struc(n)%metdata1(3,index1),Bondville_struc(n)%metdata2(3,index1)
     endif
     if ((Bondville_struc(n)%metdata1(4,t).ne.Bondville_struc(n)%undef).and.        &
          (Bondville_struc(n)%metdata2(4,t).ne.Bondville_struc(n)%undef)) then
        uwind(t) = wt1 * Bondville_struc(n)%metdata1(4,index1) +                     &
             wt2 * Bondville_struc(n)%metdata2(4,index1)
        !            write(LIS_logunit,*) Bondville_struc(n)%metdata1(4,index1),Bondville_struc(n)%metdata2(4,index1)
     endif
     if ((Bondville_struc(n)%metdata1(5,t).ne.Bondville_struc(n)%undef).and.        &
          (Bondville_struc(n)%metdata2(5,t).ne.Bondville_struc(n)%undef)) then
        vwind(t) = 0.0
        !            write(LIS_logunit,*) Bondville_struc(n)%metdata1(5,index1),Bondville_struc(n)%metdata2(5,index1)
     endif
     if ((Bondville_struc(n)%metdata1(6,t).ne.Bondville_struc(n)%undef).and.        &
          (Bondville_struc(n)%metdata2(6,t).ne.Bondville_struc(n)%undef)) then
        swdown(t) = wt1 * Bondville_struc(n)%metdata1(6,index1) +                    &
             wt2 * Bondville_struc(n)%metdata2(6,index1)
        !            write(LIS_logunit,*) Bondville_struc(n)%metdata1(6,index1),Bondville_struc(n)%metdata2(6,index1)
     endif
     if ((Bondville_struc(n)%metdata1(7,t).ne.Bondville_struc(n)%undef).and.        &
          (Bondville_struc(n)%metdata2(7,t).ne.Bondville_struc(n)%undef)) then
        lwdown(t) = wt1 * Bondville_struc(n)%metdata1(7,index1) +                    &
             wt2 * Bondville_struc(n)%metdata2(7,index1)
        !            write(LIS_logunit,*) Bondville_struc(n)%metdata1(7,index1),Bondville_struc(n)%metdata2(7,index1)
     endif
     ! Precip rate is not interpolated, but carried forward
     if (Bondville_struc(n)%metdata1(8,t).ne.Bondville_struc(n)%undef) then
        pcp(t) = Bondville_struc(n)%metdata1(8,index1)
     endif
  enddo
end subroutine timeinterp_Bondville

!------------------------------------------------------------------------------------------------

real function EsFuncT(T) result(E)

  implicit none

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !C  PURPOSE:  TO CALCULATE VALUES OF SAT. VAPOR PRESSURE E [ Pa ]
  !C            FORMULAS AND CONSTANTS FROM ROGERS AND YAU, 1989.
  !C
  !C                         ADDED BY PABLO J. GRUNMANN, 7/9/97.
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  real, intent(in) :: T     ! Temperature [ K ]

  real, parameter  :: TO    = 273.15
  real, parameter  :: CPV   = 1870.0  ! Specific heat of water vapor  [ J kg{-1} K{-1} ]
  real, parameter  :: RV    = 461.5   ! Water vapor gas constant      [ J kg{-1} K{-1} ]
  real, parameter  :: CW    = 4187.0  ! Specific heat of liquid water [ J kg{-1} K{-1} ]
  real, parameter  :: ESO   = 611.2   ! Sat. vapor pres. at T = T0    [ Pa ]
  real, parameter  :: LVH2O = 2.501E6 ! Latent Heat of Vaporization   [ J kg{-1} ]

  real :: LW

  !     CLAUSIUS-CLAPEYRON: DES/DT = L*ES/(RV*T^2)

  LW = LVH2O - (CW - CPV) * (T - TO)
  E = ESO * EXP(LW * (1.0/TO - 1.0/T) /RV)

end function esfunct

