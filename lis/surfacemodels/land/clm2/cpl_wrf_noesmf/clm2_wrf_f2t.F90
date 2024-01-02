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
! !ROUTINE: clm2_wrf_f2t
! \label{clm2_wrf_f2t}
!
! !REVISION HISTORY:
!  15 Oct 2002: Sujay Kumar; Initial Code
!
! !INTERFACE:
subroutine clm2_wrf_f2t(n)
! !USES:      
  use ESMF
  use LIS_precisionMod
  use LIS_coreMod ,      only : LIS_rc
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_FORC_AttributesMod
  use LIS_logMod,         only : LIS_verify
  use clm2_lsmMod        
  use clm2_varcon  , only : rair, cpair, po2, pco2, tcrit, tfrz    

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the Clm
!  model tiles in a coupled mode to WRF
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

!  lisimpdataname(1) = 'Incident Shortwave Radiation'
!  lisimpdataname(2) = 'Incident Direct Surface Shortwave Radiation'
!  lisimpdataname(3) = 'Incident Diffuse Surface Shortwave Radiation'
!  lisimpdataname(4) = 'Incident Longwave Radiation'
!  lisimpdataname(5) = 'Near Surface Specific Humidity'
!  lisimpdataname(6) = 'Surface Pressure'
!  lisimpdataname(7) = 'Near Surface Air Temperature'
!  lisimpdataname(8) = 'Rainfall Rate'
!  lisimpdataname(9) = 'Snowfall Rate'
!  lisimpdataname(10) = 'Northward Wind'
!  lisimpdataname(11) = 'Eastward Wind'
!  lisimpdataname(12) = 'Height of Forcing Variables'
!  lisimpdataname(13) = 'Surface Exchange Coefficient for Heat'
!  lisimpdataname(14) = 'Surface Exchange Coefficient for Momentum'
!  lisimpdataname(15) = '2m Surface Exchange Coefficient for Heat'
!  lisimpdataname(16) = '2m Surface Exchange Coefficient for Moisture'
!  lisimpdataname(17) = 'Surface Emissivity'
!  lisimpdataname(18) = 'Saturated Mixing Ratio'
!  lisimpdataname(19) = 'Cosine of Zenith Angle'
!  lisimpdataname(20) = 'Surface Albedo'

  real(r8), allocatable :: solar(:)
  real(r8), allocatable :: prcp(:)
  real(r8) :: forc_vp
  integer            :: t,status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
  type(ESMF_Field)   :: chField, chs2Field, cqs2Field,q2satField,coszField, zField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: ch(:),chs2(:),cqs2(:),q2sat(:),cosz(:),zval(:)


  allocate(solar(LIS_rc%ntiles(n)))
  allocate(prcp(LIS_rc%ntiles(n)))

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: Tair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: Qair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: uwind')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: vwind')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: Psurf')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: Rainf')

  !  call ESMF_StateGet(LIS_FORC_State(n),"Convective Rainfall Rate",cpcpField,&
  !       rc=status)
  !  call LIS_verify(status, 'Error in ESMF_StateGet: Rainf_cp')

  !  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),snowfField,&
  !       rc=status)
  !  call LIS_verify(status, 'Error in ESMF_StateGet: Snowf_f')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Ch%varname(1)),chField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: ch')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Chs2%varname(1)),chs2Field,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: chs2')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Cqs2%varname(1)),cqs2Field,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: cqs2')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Q2sat%varname(1)),q2satField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: q2sat')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_cosz%varname(1)),coszField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: cosz')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Forc_Hgt%varname(1)),zField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: zval')

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status)

  !  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  !  call LIS_verify(status)
  !  call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
  !  call LIS_verify(status)

  call ESMF_FieldGet(chField,localDE=0,farrayPtr=ch,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(chs2Field,localDE=0,farrayPtr=chs2,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(cqs2Field,localDE=0,farrayPtr=cqs2,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2satField,localDE=0,farrayPtr=q2sat,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(coszField,localDE=0,farrayPtr=cosz,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(zField,localDE=0,farrayPtr=zval,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     clm2_struc(n)%clm(t)%forc_t        = tmp(t)
     clm2_struc(n)%clm(t)%forc_q        = q2(t)
     solar(t)             = swd(t)
     clm2_struc(n)%clm(t)%forc_solad(1) = solar(t)*35.0/100.0
     clm2_struc(n)%clm(t)%forc_solad(2) = solar(t)*35.0/100.0
     clm2_struc(n)%clm(t)%forc_solai(1) = solar(t)*15.0/100.0
     clm2_struc(n)%clm(t)%forc_solai(2) = solar(t)*15.0/100.0
     clm2_struc(n)%clm(t)%forc_lwrad    = lwd(t)
     clm2_struc(n)%clm(t)%forc_u        = uwind(t)
     clm2_struc(n)%clm(t)%forc_v        = vwind(t)
     clm2_struc(n)%clm(t)%forc_pbot     = psurf(t)
     prcp(t)              = pcp(t)

     clm2_struc(n)%clm(t)%forc_hgt      = 0.5*zval(t)
     clm2_struc(n)%clm(t)%forc_ch       = ch(t)
     clm2_struc(n)%clm(t)%forc_chs2       = chs2(t)
     clm2_struc(n)%clm(t)%forc_cqs2       = cqs2(t)
     clm2_struc(n)%clm(t)%forc_q2sat    = q2sat(t)
     clm2_struc(n)%clm(t)%forc_cosz     = cosz(t)

     !       clm2_struc(n)%clm(t)%forc_hgt_u    = clm2_struc(n)%clm(t)%forc_hgt+&
     !           clm2_struc(n)%clm(t)%displa+clm2_struc(n)%clm(t)%z0m
     !       clm2_struc(n)%clm(t)%forc_hgt_t    = clm2_struc(n)%clm(t)%forc_hgt+&
     !           clm2_struc(n)%clm(t)%displa+clm2_struc(n)%clm(t)%z0m
     !       clm2_struc(n)%clm(t)%forc_hgt_q    = clm2_struc(n)%clm(t)%forc_hgt+&
     !           clm2_struc(n)%clm(t)%displa+clm2_struc(n)%clm(t)%z0m

     if (prcp(t) > 0.) then
        if (clm2_struc(n)%clm(t)%forc_t > (tfrz + tcrit)) then
           clm2_struc(n)%clm(t)%itypprc   = 1
           clm2_struc(n)%clm(t)%forc_rain = prcp(t)
           clm2_struc(n)%clm(t)%forc_snow = 0.
        else
           clm2_struc(n)%clm(t)%itypprc   = 2
           clm2_struc(n)%clm(t)%forc_rain = 0.
           clm2_struc(n)%clm(t)%forc_snow = prcp(t)
        endif
     else
        clm2_struc(n)%clm(t)%itypprc   = 0
        clm2_struc(n)%clm(t)%forc_rain = 0.
        clm2_struc(n)%clm(t)%forc_snow = 0
     endif

     ! Derive new fields (potential temperature, vapor pressure, 
     ! air density, CO2, and O2) and copy solar radiations

     !=== Cuurently potential temperature set to 2 m temperature since 
     !=== we only get surface pressure in our forcing and elevation differences 
     !=== are accounted for in for.f

     !===LDAS modification: Slight change to be consistent with our forcing dataset
     !          clm2_struc(n)%clm(t)%forc_th  = clm2_struc(n)%clm(t)%forc_t * (clm2_struc(n)%clm(t)%forc_psrf/clm2_struc(n)%clm(t)%forc_pbot)**(rair/cpair)
     clm2_struc(n)%clm(t)%forc_th  = clm2_struc(n)%clm(t)%forc_t * &
          (clm2_struc(n)%clm(t)%forc_pbot/clm2_struc(n)%clm(t)%forc_pbot)**(rair/cpair)
     forc_vp  = clm2_struc(n)%clm(t)%forc_q*clm2_struc(n)%clm(t)%forc_pbot / &
          (0.622+0.378*clm2_struc(n)%clm(t)%forc_q)
     clm2_struc(n)%clm(t)%forc_rho = &
          (clm2_struc(n)%clm(t)%forc_pbot-0.378*forc_vp) / &
          (rair*clm2_struc(n)%clm(t)%forc_t)
  enddo

  deallocate(solar)
  deallocate(prcp)
end subroutine clm2_wrf_f2t
