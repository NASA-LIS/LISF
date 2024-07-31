!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: template_f2t
! \label{template_f2t}
!
! !REVISION HISTORY: 
! 21 Jul 2004: Sujay Kumar   Initial Specification
! 23 Oct 2007: Kristi Arsenault, Implemented code for LISv5.0
! 
! !INTERFACE:
subroutine template_f2t(n)

! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_surface
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_FORC_AttributesMod 
  use LIS_logMod,         only : LIS_verify
  use template_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!
!  Forcing-only option (template) for calling the forcing transfer routines.
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

  integer            :: t,tid,status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,snowfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)


  if(LIS_FORC_Tair%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_Qair%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_SWdown%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
    call LIS_verify(status)
  endif

  if(LIS_FORC_LWdown%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
    call LIS_verify(status)
  endif

  if(LIS_FORC_Wind_E%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_Wind_N%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_Psurf%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_Rainf%selectOpt.eq.1) then
    call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)
  endif
  
  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LIS_verify(status)
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),snowfField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
     call LIS_verify(status)
  endif

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tid=LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id

     if(LIS_FORC_Tair%selectOpt.eq.1) then
       template_struc(n)%template(t)%tair=tmp(tid)
     endif
     if(LIS_FORC_Qair%selectOpt.eq.1) then
       template_struc(n)%template(t)%qair=q2(tid)
     endif
     if(LIS_FORC_SWdown%selectOpt.eq.1) then
       template_struc(n)%template(t)%swdown=swd(tid)
     endif
     if(LIS_FORC_LWdown%selectOpt.eq.1) then
       template_struc(n)%template(t)%lwdown=lwd(tid)
     endif
     if(LIS_FORC_Wind_E%selectOpt.eq.1) then
       template_struc(n)%template(t)%uwind=uwind(tid)
     endif
     if(LIS_FORC_Wind_N%selectOpt.eq.1) then
       template_struc(n)%template(t)%vwind=vwind(tid)
     endif
     if(LIS_FORC_Psurf%selectOpt.eq.1) then
       template_struc(n)%template(t)%psurf=psurf(tid)
     endif
     if(pcp(tid).ne.LIS_rc%udef) then
        template_struc(n)%template(t)%rainf=pcp(tid)
     else
        template_struc(n)%template(t)%rainf=0.0
     endif
     if(LIS_FORC_CRainf%selectOpt.eq.1) then 
        if(cpcp(tid).ne.LIS_rc%udef) then 
           template_struc(n)%template(t)%rainf_c=cpcp(tid)
        else
           template_struc(n)%template(t)%rainf_c=0.0
        endif
     endif
     if(LIS_FORC_Snowf%selectOpt.eq.1) then 
        if(snowf(tid).ne.LIS_rc%udef) then
           template_struc(n)%template(t)%snowf=snowf(tid)
        else
           template_struc(n)%template(t)%snowf=0.0
        endif
     endif
  enddo


end subroutine template_f2t
