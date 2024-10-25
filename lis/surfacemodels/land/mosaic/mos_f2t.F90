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
! !ROUTINE: mos_f2t
! \label{mos_f2t}
!
! !REVISION HISTORY:
!  28 Jan 2002: Jon Gottschalck; Added option for different number of forcing variables
! 25 Sep 2007: Sujay Kumar; Upgraded for LIS5.0
!
! !INTERFACE:
subroutine mos_f2t(n)
! !USES:      
  use ESMF
  use LIS_coreMod ,       only : LIS_rc, LIS_surface
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod,  only : LIS_FORC_State
  use LIS_logMod,         only : LIS_verify
  use mos_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the mosaic
!  model tiles. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

  integer            :: t,tid, status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField,fhgtField,acondField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: fheight(:),acond(:)

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),&
                     tmpField, rc=status)
  call LIS_verify(status,'mos_f2t: error getting Tair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),&
                     q2Field, rc=status)
  call LIS_verify(status,'mos_f2t: error getting Qair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),&
                     swdField, rc=status)
  call LIS_verify(status,'mos_f2t: error getting SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),&
                     lwdField, rc=status)
  call LIS_verify(status,'mos_f2t: error getting LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),&
                     uField, rc=status)
  call LIS_verify(status,'mos_f2t: error getting Wind_E')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),&
                     vField, rc=status)
  call LIS_verify(status,'mos_f2t: error getting Wind_N')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),&
                     psurfField, rc=status)
  call LIS_verify(status, 'mos_f2t: error getting PSurf')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),&
                     pcpField, rc=status)
  call LIS_verify(status,'mos_f2t: error getting Rainf')

  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),&
                        cpcpField, rc=status)
     call LIS_verify(status,'mos_f2t: error getting CRainf')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),&
                        snowfField, rc=status)
     call LIS_verify(status,'mos_f2t: error getting Snowf')
  endif

  if ( mos_struc(n)%forcing_z > 0 ) then
     call ESMF_StateGet(LIS_FORC_State(n),"Height of Atmospheric Forcing",&
                             fhgtField, rc=status)
     call LIS_verify(status)
  endif

  if ( mos_struc(n)%forcing_ch > 0 ) then
     call ESMF_StateGet(LIS_FORC_State(n),&
                        "Surface Exchange Coefficient for Heat",acondField,&
                        rc=status)
     call LIS_verify(status)
  endif

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

  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LIS_verify(status, 'mos_f2t: error retrieving cpcp')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
     call LIS_verify(status,'mos_f2t: error retrieving snowf')
  endif

  if ( mos_struc(n)%forcing_z > 0 ) then
     call ESMF_FieldGet(fhgtField,localDE=0,farrayPtr=fheight,rc=status)
     call LIS_verify(status)
  endif
  
  if ( mos_struc(n)%forcing_ch > 0 ) then
     call ESMF_FieldGet(acondField,localDE=0,farrayPtr=acond,rc=status)
     call LIS_verify(status)
  endif

  mos_struc(n)%forc_count = mos_struc(n)%forc_count +1
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id

     mos_struc(n)%mos(t)%tair=mos_struc(n)%mos(t)%tair + tmp(tid)
     mos_struc(n)%mos(t)%qair=mos_struc(n)%mos(t)%qair + q2(tid)
     mos_struc(n)%mos(t)%swdown=mos_struc(n)%mos(t)%swdown + swd(tid)
     mos_struc(n)%mos(t)%lwdown=mos_struc(n)%mos(t)%lwdown + lwd(tid)
     mos_struc(n)%mos(t)%psurf=mos_struc(n)%mos(t)%psurf + psurf(tid)

     if ( mos_struc(n)%forcing_ch > 0 ) then
        if(acond(tid).le.0) then
           mos_struc(n)%mos(t)%acond=mos_struc(n)%mos(t)%acond + 0.00001     
        else
           mos_struc(n)%mos(t)%acond=mos_struc(n)%mos(t)%acond + acond(tid)
        endif
     else
        mos_struc(n)%mos(t)%acond=LIS_rc%udef
     endif

     if ( mos_struc(n)%forcing_z > 0 ) then
        mos_struc(n)%mos(t)%obsz = mos_struc(n)%mos(t)%obsz + fheight(tid)
     else  
        !set default mosaic value
        mos_struc(n)%mos(t)%obsz = mos_struc(n)%mos(t)%obsz + 50.0  
     endif

     if(uwind(tid).ne.LIS_rc%udef) then 
        mos_struc(n)%mos(t)%uwind=mos_struc(n)%mos(t)%uwind + uwind(tid)
     else
        mos_struc(n)%mos(t)%uwind=mos_struc(n)%mos(t)%uwind + 0.0
     endif
     if(vwind(tid).ne.LIS_rc%udef) then 
        mos_struc(n)%mos(t)%vwind=mos_struc(n)%mos(t)%vwind + vwind(tid)
     else
        mos_struc(n)%mos(t)%vwind=mos_struc(n)%mos(t)%vwind + 0.0
     endif

     if(pcp(tid).ne.LIS_rc%udef) then 
        mos_struc(n)%mos(t)%rainf=mos_struc(n)%mos(t)%rainf + pcp(tid)
     else
        mos_struc(n)%mos(t)%rainf=mos_struc(n)%mos(t)%rainf + 0.0
     endif
     if ( LIS_FORC_CRainf%selectOpt.eq.1) then
        if(cpcp(tid).ne.LIS_rc%udef) then 
           mos_struc(n)%mos(t)%rainf_c=mos_struc(n)%mos(t)%rainf_c+ cpcp(tid)
        else
           mos_struc(n)%mos(t)%rainf_c=mos_struc(n)%mos(t)%rainf_c + 0.0
        endif
     else
        mos_struc(n)%mos(t)%rainf_c=mos_struc(n)%mos(t)%rainf_c + 0.0
     endif
     ! If there is snowf add it to precipitation.  Mosaic does not use
     ! separate rainf and snowf.  It determines what to do with precipitation.
     if(LIS_FORC_Snowf%selectOpt.eq.1) then 
        if(snowf(tid).ne.LIS_rc%udef) then 
           mos_struc(n)%mos(t)%rainf=mos_struc(n)%mos(t)%rainf+snowf(tid)
        endif
     endif
     mos_struc(n)%mos(t)%snowf=mos_struc(n)%mos(t)%snowf+ 0.0
  enddo

end subroutine mos_f2t
