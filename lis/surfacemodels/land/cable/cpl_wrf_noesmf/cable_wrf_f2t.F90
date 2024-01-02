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
! !ROUTINE: cable_wrf_f2t
! \label{cable_wrf_f2t}
!
! !REVISION HISTORY:
!  15 Oct 2002: Sujay Kumar; Initial Code
!  17 Oct 2011: Claire Carouge (ccc), CABLE LSM improvements
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
!
! !INTERFACE:
subroutine cable_wrf_f2t(n)
! !USES:      
  use ESMF
  use LIS_precisionMod
  use LIS_coreMod ,      only : LIS_rc, LIS_surface
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_FORC_AttributesMod
  use LIS_logMod,         only : LIS_verify
  use cable_lsmMod        
  use cable_physical_constants, only : tfrz   

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the CABLE
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

  real(r8), pointer  :: prcp(:)
  integer            :: tid
  integer            :: t,status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
  type(ESMF_Field)   :: chField,chs2Field,cqs2Field,q2satField,coszField,zField,co2Field
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: ch(:),chs2(:),cqs2(:),q2sat(:),cosz(:),zval(:),co2(:)

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

  !ccc
  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_co2%varname(1)),co2Field,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: co2')

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

  call ESMF_FieldGet(co2Field,localDE=0,farrayPtr=co2,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
     cable_struc(n)%cable(t)%tair          = tmp(tid)
     cable_struc(n)%cable(t)%qair         = q2(tid)
     cable_struc(n)%cable(t)%swdown  = swd(tid)
     cable_struc(n)%cable(t)%lwdown   = lwd(tid)
     cable_struc(n)%cable(t)%uwind        = uwind(tid)
     cable_struc(n)%cable(t)%vwind        = vwind(tid)
     cable_struc(n)%cable(t)%psurf     = psurf(tid)/100.
     prcp(t)              = pcp(tid)

     cable_struc(n)%cable(t)%za      = 0.5*zval(tid)
     cable_struc(n)%cable(t)%forc_ch       = ch(tid)
     cable_struc(n)%cable(t)%forc_chs2       = chs2(tid)
     cable_struc(n)%cable(t)%forc_cqs2       = cqs2(tid)
     cable_struc(n)%cable(t)%q2sat    = q2sat(tid)
     cable_struc(n)%cable(t)%coszen     = cosz(tid)

     if (prcp(tid) > 0.) then
        ! Threshold taken from CLM2. 
        if (cable_struc(n)%cable(t)%tair > (tfrz + 2.5)) then
           cable_struc(n)%cable(t)%rainf = prcp(tid)
           cable_struc(n)%cable(t)%snowf = 0.
        else
           cable_struc(n)%cable(t)%rainf = 0.
           cable_struc(n)%cable(t)%snowf = prcp(tid)
        endif
     else
        cable_struc(n)%cable(t)%rainf = 0.
        cable_struc(n)%cable(t)%snowf = 0
     endif

     ! CO2
     if (co2(1) <= 0.) then
        cable_struc(n)%cable(t)%co2 = cable_struc(n)%fixedco2/1000000.0
     else
        cable_struc(n)%cable(t)%co2 = co2(1)
     end if
  enddo

  deallocate(prcp)

end subroutine cable_wrf_f2t
