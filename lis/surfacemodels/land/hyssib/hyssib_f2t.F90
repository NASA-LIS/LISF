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
! !ROUTINE: hyssib_f2t
! \label{hyssib_f2t}

! !REVISION HISTORY:
!  21 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
!  25 Aug 2007: Chuck Alonge, Updates for LIS 5.0 Compliance
!  27 Oct 2010: David Mocko, changes for HY-SSiB in LIS6.1
!
! !INTERFACE:
subroutine hyssib_f2t(n)
! !USES:
  use ESMF
  use LIS_coreMod       , only : LIS_rc, LIS_surface
  use LIS_FORC_AttributesMod 
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_logMod,         only : LIS_verify
  use hyssib_lsmMod

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
!  This routine transfers the LIS provided forcing onto the Hyssib
!  model tiles. 
! 
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!EOP

  integer            :: t,tid,status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField,fhgtField,snowfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:),fheight(:)

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),&
                     tmpField, rc=status)
  call LIS_verify(status,'hyssib_f2t: error getting Tair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),&
                     q2Field, rc=status)
  call LIS_verify(status,'hyssib_f2t: error getting Qair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),&
                     swdField, rc=status)
  call LIS_verify(status,'hyssib_f2t: error getting SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),&
                     lwdField, rc=status)
  call LIS_verify(status,'hyssib_f2t: error getting LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),&
                     uField, rc=status)
  call LIS_verify(status,'hyssib_f2t: error getting Wind_E')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),&
                     vField, rc=status)
  call LIS_verify(status,'hyssib_f2t: error getting Wind_N')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),&
                     psurfField, rc=status)
  call LIS_verify(status, 'hyssib_f2t: error getting PSurf')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),&
                     pcpField, rc=status)
  call LIS_verify(status,'hyssib_f2t: error getting Rainf')

  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),&
                        cpcpField, rc=status)
     call LIS_verify(status,'hyssib_f2t: error getting CRainf')
  endif

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Forc_Hgt%varname(1)),&
                        fhgtField, rc=status)
     call LIS_verify(status,'hyssib_f2t: error getting Forc_Hgt')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Snowf%varname(1)),&
                        snowfField, rc=status)
     call LIS_verify(status,'hyssib_f2t: error getting Snowf')
  endif

  call ESMF_FieldGet(tmpField, localDE=0, farrayPtr=tmp, rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(q2Field, localDE=0, farrayPtr=q2, rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(swdField, localDE=0, farrayPtr=swd, rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField, localDE=0, farrayPtr=lwd, rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(uField, localDE=0, farrayPtr=uwind, rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(vField, localDE=0, farrayPtr=vwind, rc=status)
  call LIS_verify(status)
        
  call ESMF_FieldGet(psurfField, localDE=0, farrayPtr=psurf, rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(pcpField, localDE=0, farrayPtr=pcp, rc=status)
  call LIS_verify(status)

  if(LIS_FORC_CRainf%selectOpt.eq.1) then 
     call ESMF_FieldGet(cpcpField, localDE=0, farrayPtr=cpcp, rc=status)
     call LIS_verify(status, 'hyssib_f2t: error retrieving cpcp')
  endif

  if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then 
     call ESMF_FieldGet(fhgtField, localDE=0, farrayPtr=fheight, rc=status)
     call LIS_verify(status, 'hyssib_f2t: error retrieving forc_hgt')
  endif

  if(LIS_FORC_Snowf%selectOpt.eq.1) then 
     call ESMF_FieldGet(snowfField, localDE=0, farrayPtr=snowf, rc=status)
     call LIS_verify(status,'hyssib_f2t: error retrieving snowf')
  endif

  hyssib_struc(n)%forc_count = 0 
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
     hyssib_struc(n)%hyssib(t)%tair=hyssib_struc(n)%hyssib(t)%tair + tmp(tid)
     hyssib_struc(n)%hyssib(t)%qair=hyssib_struc(n)%hyssib(t)%qair + q2(tid)
     hyssib_struc(n)%hyssib(t)%swdown=hyssib_struc(n)%hyssib(t)%swdown + swd(tid)
     hyssib_struc(n)%hyssib(t)%lwdown=hyssib_struc(n)%hyssib(t)%lwdown + lwd(tid)
     hyssib_struc(n)%hyssib(t)%uwind=hyssib_struc(n)%hyssib(t)%uwind + uwind(tid)
     hyssib_struc(n)%hyssib(t)%vwind=hyssib_struc(n)%hyssib(t)%vwind + vwind(tid)
     hyssib_struc(n)%hyssib(t)%psurf=hyssib_struc(n)%hyssib(t)%psurf + psurf(tid)
     if(LIS_FORC_Forc_Hgt%selectOpt.eq.1) then
        hyssib_struc(n)%hyssib(t)%fheight=fheight(tid)
     endif

     if(pcp(tid).ne.LIS_rc%udef) then 
        hyssib_struc(n)%hyssib(t)%rainf_in=&
             hyssib_struc(n)%hyssib(t)%rainf_in + pcp(tid)
     else
        hyssib_struc(n)%hyssib(t)%rainf_in= &
              hyssib_struc(n)%hyssib(t)%rainf_in + 0.0
     endif
     if ( LIS_FORC_CRainf%selectOpt.eq.1) then
        if(cpcp(tid).ne.LIS_rc%udef) then 
           hyssib_struc(n)%hyssib(t)%rainf_cp= & 
                hyssib_struc(n)%hyssib(t)%rainf_cp + cpcp(tid)
        else
           hyssib_struc(n)%hyssib(t)%rainf_cp= & 
                 hyssib_struc(n)%hyssib(t)%rainf_cp + 0.0
        endif
     else
        hyssib_struc(n)%hyssib(t)%rainf_cp= & 
              hyssib_struc(n)%hyssib(t)%rainf_cp + 0.0
     endif
     if ( LIS_FORC_Snowf%selectOpt.eq.1) then
        if(snowf(tid).ne.LIS_rc%udef) then 
           hyssib_struc(n)%hyssib(t)%snowf_in= & 
                 hyssib_struc(n)%hyssib(t)%snowf_in + snowf(tid)
        else
           hyssib_struc(n)%hyssib(t)%snowf_in= & 
                 hyssib_struc(n)%hyssib(t)%snowf_in + 0.0
        endif
     else
        hyssib_struc(n)%hyssib(t)%snowf_in= & 
              hyssib_struc(n)%hyssib(t)%snowf_in + 0.0
     endif
  enddo
end subroutine hyssib_f2t

