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
! !ROUTINE: jules52_f2t
! \label{jules52_f2t}
!
! !REVISION HISTORY: 
! 21 Jul 2004: Sujay Kumar   Initial Specification
! 23 Oct 2007: Kristi Arsenault, Implemented code for LISv5.0
! 16 May 2016; Shugong Wang; initial implementation for JULES 4.3
! 01 Feb 2018; Shugong Wang; updated for JULES 5.0 
! 28 Nov 2018; Shugong Wang; updated for JULES 5.2 
! !INTERFACE:
subroutine jules52_f2t(n)

! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_surface, LIS_domain 
  use LIS_metforcingMod,  only : LIS_FORC_State
  use LIS_FORC_AttributesMod 
  use LIS_logMod,         only : LIS_verify
  use jules52_lsmMod
  use debug_latlon 
  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
!
! !DESCRIPTION:
!
!  Forcing-only option (jules52) for calling the forcing transfer routines.
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP

  integer            :: t,tid,status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,rainField,crainField,snowfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),rain(:),crain(:)
  
  real                 :: lat, lon
  integer              :: row, col, pft 
  INTEGER    COMM, RANK, IERROR

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status)
  
  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),rainField,&
       rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(rainField,localDE=0,farrayPtr=rain,rc=status)
  call LIS_verify(status)
  
  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),crainField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(crainField,localDE=0,farrayPtr=crain,rc=status)
     call LIS_verify(status)
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),snowfField,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
     call LIS_verify(status)
  endif
  
  !!! set the forcing counter 
  jules52_struc(n)%forc_count = jules52_struc(n)%forc_count + 1 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
      col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
      lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
      lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

      lat_d = lat
      lon_d = lon 
     
     tid=LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id

     jules52_struc(n)%jules52(t)%tair   = jules52_struc(n)%jules52(t)%tair   + tmp(tid)
     jules52_struc(n)%jules52(t)%qair   = jules52_struc(n)%jules52(t)%qair   + q2(tid)
     jules52_struc(n)%jules52(t)%swdown = jules52_struc(n)%jules52(t)%swdown + swd(tid)
     jules52_struc(n)%jules52(t)%lwdown = jules52_struc(n)%jules52(t)%lwdown + lwd(tid)
     jules52_struc(n)%jules52(t)%wind_e = jules52_struc(n)%jules52(t)%wind_e + uwind(tid)
     jules52_struc(n)%jules52(t)%wind_n = jules52_struc(n)%jules52(t)%wind_n + vwind(tid)
     jules52_struc(n)%jules52(t)%psurf  = jules52_struc(n)%jules52(t)%psurf  + psurf(tid)
     jules52_struc(n)%jules52(t)%rainf  = jules52_struc(n)%jules52(t)%rainf  + rain(tid)
     
     if(LIS_FORC_CRainf%selectOpt .eq. 1) then
       if(crain(tid) .ne. LIS_rc%udef) then 
          jules52_struc(n)%jules52(t)%rainf_c = jules52_struc(n)%jules52(t)%rainf_c + crain(tid)
       else
          jules52_struc(n)%jules52(t)%rainf_c = jules52_struc(n)%jules52(t)%rainf_c + 0.0
       endif
     else
       jules52_struc(n)%jules52(t)%rainf_c = 0.0
     endif

     if(LIS_FORC_Snowf%selectOpt .eq. 1) then
        if(snowf(tid) .ne. LIS_rc%udef) then 
          jules52_struc(n)%jules52(t)%snowf = jules52_struc(n)%jules52(t)%snowf + snowf(tid)
        else
          jules52_struc(n)%jules52(t)%snowf = jules52_struc(n)%jules52(t)%snowf + 0.0
        endif
     else
       jules52_struc(n)%jules52(t)%snowf = 0.0
     endif
  enddo

end subroutine jules52_f2t
