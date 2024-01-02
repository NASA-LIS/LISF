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
! !ROUTINE: noah271_wrf_esmff2t
! \label{noah271_wrf_esmff2t}
!
! !REVISION HISTORY:
!  15 Oct 2002: Sujay Kumar; Initial Code
!
! !INTERFACE:
subroutine noah271_wrf_esmff2t(n)
! !USES:      
  use ESMF
  use LIS_coreMod ,      only : LIS_rc
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_FORC_AttributesMod
  use LIS_logMod,         only : LIS_verify
  use noah271_lsmMod,        only : noah271_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the Noah
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
!  lisimpdataname(15) = 'Surface Emissivity'
!  lisimpdataname(16) = 'Saturated Mixing Ratio'
!  lisimpdataname(17) = 'Cosine of Zenith Angle'
!  lisimpdataname(18) = 'Surface Albedo'
!  lisimpdataname(19) = 'Sea Ice Mask'

  integer            :: t,status
  type(ESMF_Field)   :: tmpField,qairField,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
  type(ESMF_Field)   :: chField, q2satField,emissField, zField, xiceField
  type(ESMF_Field)   :: qsfcField, chs2Field, cqs2Field,q2Field
  type(ESMF_Field)   :: t2Field, th2Field, tmnField
  real,allocatable       :: tmp(:),qair(:),uwind(:),vwind(:),snowf(:)
  real,allocatable       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,allocatable       :: ch(:),q2sat(:),emiss(:),zval(:), xice(:)
  real,allocatable       :: qsfc(:),chs2(:),cqs2(:),t2(:),th2(:),q2(:),tmn(:)

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: Tair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),qairField,&
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

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),snowfField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: Snowf_f')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Ch%varname(1)),chField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: ch')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Q2sat%varname(1)),q2satField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: q2sat')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Emiss%varname(1)),emissField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: emiss')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Forc_Hgt%varname(1)),zField,&
       rc=status)
  call LIS_verify(status, 'Error in ESMF_StateGet: zval')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Xice%varname(1)),&
       xiceField,rc=status)
  call LIS_verify(status,'noah271_f2t: error getting xice')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_QSFC%varname(1)),&
       qsfcField,rc=status)
  call LIS_verify(status,'noah271_f2t: error getting qsfc')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_chs2%varname(1)),&
       chs2Field,rc=status)
  call LIS_verify(status,'noah271_f2t: error getting chs2')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_cqs2%varname(1)),&
       cqs2Field,rc=status)
  call LIS_verify(status,'noah271_f2t: error getting cqs2')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_t2%varname(1)),&
       t2Field,rc=status)
  call LIS_verify(status,'noah271_f2t: error getting t2')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_q2%varname(1)),&
       q2Field,rc=status)
  call LIS_verify(status,'noah271_f2t: error getting q2')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_th2%varname(1)),&
       th2Field,rc=status)
  call LIS_verify(status,'noah271_f2t: error getting th2')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_tmn%varname(1)),&
       tmnField,rc=status)
  call LIS_verify(status,'noah271_f2t: error getting TMN')



  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(qairField,localDE=0,farrayPtr=qair,rc=status)
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

  call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(chField,localDE=0,farrayPtr=ch,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2satField,localDE=0,farrayPtr=q2sat,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(emissField,localDE=0,farrayPtr=emiss,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(zField,localDE=0,farrayPtr=zval,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(xiceField,localDE=0, farrayPtr=xice,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving xice')

  call ESMF_FieldGet(qsfcField,localDE=0, farrayPtr=qsfc,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving qsfc')

  call ESMF_FieldGet(chs2Field,localDE=0, farrayPtr=chs2,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving chs2')

  call ESMF_FieldGet(cqs2Field,localDE=0, farrayPtr=cqs2,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving cqs2')

  call ESMF_FieldGet(t2Field,localDE=0, farrayPtr=t2,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving t2')

  call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving q2')

  call ESMF_FieldGet(th2Field,localDE=0, farrayPtr=th2,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving th2')

  call ESMF_FieldGet(tmnField,localDE=0, farrayPtr=tmn,rc=status)
  call LIS_verify(status,'noah271_f2t: error retrieving tmn')

  do t=1,LIS_rc%ntiles(n)
     noah271_struc(n)%noah(t)%tair=tmp(t)
     noah271_struc(n)%noah(t)%qair=qair(t)
     noah271_struc(n)%noah(t)%swdown=swd(t)
     noah271_struc(n)%noah(t)%lwdown=lwd(t)
     noah271_struc(n)%noah(t)%uwind=uwind(t)
     noah271_struc(n)%noah(t)%vwind=vwind(t)
     noah271_struc(n)%noah(t)%psurf=psurf(t)
     noah271_struc(n)%noah(t)%xice=xice(t)

     noah271_struc(n)%noah(t)%qsfc=qsfc(t)
     noah271_struc(n)%noah(t)%chs2=chs2(t)
     noah271_struc(n)%noah(t)%cqs2=cqs2(t)

     noah271_struc(n)%noah(t)%t2=t2(t)
     noah271_struc(n)%noah(t)%th2=th2(t)
     noah271_struc(n)%noah(t)%tempbot=tmn(t)
     noah271_struc(n)%noah(t)%q2=q2(t)

     if(pcp(t).ne.LIS_rc%udef) then 
        noah271_struc(n)%noah(t)%rainf=pcp(t)
     else
        noah271_struc(n)%noah(t)%rainf=0.0
     endif
     if(snowf(t).ne.LIS_rc%udef) then 
        noah271_struc(n)%noah(t)%snowf=snowf(t)
     else
        noah271_struc(n)%noah(t)%snowf=0.0
     endif
     if(ch(t).ne.LIS_rc%udef) then 
        noah271_struc(n)%noah(t)%ch=ch(t)
     endif
     if(q2sat(t).ne.LIS_rc%udef) then 
        noah271_struc(n)%noah(t)%q2sat=q2sat(t)
     endif
     if(emiss(t).ne.LIS_rc%udef) then 
        noah271_struc(n)%noah(t)%emiss=emiss(t)
     endif
     if(zval(t).ne.LIS_rc%udef) then 
        noah271_struc(n)%noah(t)%z=zval(t)
     endif
     
  enddo

end subroutine noah271_wrf_esmff2t
