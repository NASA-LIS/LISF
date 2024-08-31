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
! !ROUTINE: cable_f2t
! \label{cable_f2t}
!
! !REVISION HISTORY:
!  21 Jul 2004: Sujay Kumar, Initial Specification
!  23 Oct 2007: Kristi Arsenault, Implemented code for LISv5.0
!  26 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!
! !INTERFACE:
subroutine cable_f2t(n)
! !USES:
  use ESMF
  use LIS_coreMod,        only : LIS_rc, LIS_surface
  use LIS_metforcingMod, only : LIS_FORC_State
  use LIS_FORC_AttributesMod
  use LIS_logMod,         only : LIS_verify
  use cable_lsmMod

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
!
! !DESCRIPTION:
!  This subroutine transfers the LIS provided forcing onto the
!  CABLE model tiles.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!  \end{description}
!EOP
  integer            :: t,tid,status
  type(ESMF_Field)   :: tmpField,q2Field,swdField,lwdField
  type(ESMF_Field)   :: uField,vField,psurfField
  type(ESMF_Field)   :: pcpField,cpcpField,snowfField
  real, pointer      :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real, pointer      :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  call ESMF_StateGet(LIS_FORC_State(n),                            &
       trim(LIS_FORC_Tair%varname(1)),tmpField,rc=status)
  call LIS_verify(status,'cable_f2t: error getting Tair')

  call ESMF_StateGet(LIS_FORC_State(n),                            &
       trim(LIS_FORC_Qair%varname(1)),q2Field,rc=status)
  call LIS_verify(status,'cable_f2t: error getting Qair')

  call ESMF_StateGet(LIS_FORC_State(n),                            &
       trim(LIS_FORC_SWdown%varname(1)),swdField,rc=status)
  call LIS_verify(status,'cable_f2t: error getting SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),                            &
       trim(LIS_FORC_LWdown%varname(1)),lwdField,rc=status)
  call LIS_verify(status,'cable_f2t: error getting LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),                            &
       trim(LIS_FORC_Wind_E%varname(1)),uField,rc=status)
  call LIS_verify(status,'cable_f2t: error getting Wind_E')

  call ESMF_StateGet(LIS_FORC_State(n),                            &
       trim(LIS_FORC_Wind_N%varname(1)),vField,rc=status)
  call LIS_verify(status,'cable_f2t: error getting Wind_N')

  call ESMF_StateGet(LIS_FORC_State(n),                            &
       trim(LIS_FORC_Psurf%varname(1)),psurfField,rc=status)
  call LIS_verify(status,'cable_f2t: error getting PSurf')

  call ESMF_StateGet(LIS_FORC_State(n),                            &
       trim(LIS_FORC_Rainf%varname(1)),pcpField,rc=status)
  call LIS_verify(status,'cable_f2t: error getting Rainf')

  if (LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),                         &
          trim(LIS_FORC_Snowf%varname(1)),snowfField,rc=status)
     call LIS_verify(status,'cable_f2t: error getting Snowf')
  endif

  if (LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_StateGet(LIS_FORC_State(n),                         &
          trim(LIS_FORC_CRainf%varname(1)),cpcpField,rc=status)
     call LIS_verify(status,'cable_f2t: error getting CRainf')
  endif

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status,'cable_f2t: error retrieving Tair')

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status,'cable_f2t: error retrieving Qair')

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status,'cable_f2t: error retrieving SWdown')

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status,'cable_f2t: error retrieving LWdown')

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status,'cable_f2t: error retrieving Wind_E')

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status,'cable_f2t: error retrieving Wind_N')

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status,'cable_f2t: error retrieving PSurf')

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status,'cable_f2t: error retrieving Rainf')

  if (LIS_FORC_Snowf%selectOpt.eq.1) then
     call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
     call LIS_verify(status,'cable_f2t: error retrieving Snowf')
  endif

  if (LIS_FORC_CRainf%selectOpt.eq.1) then
     call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
     call LIS_verify(status,'cable_f2t: error retrieving CRainf')
  endif

  do t = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     tid = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%tile_id
     cable_struc(n)%cable(t)%tair = tmp(tid)
     cable_struc(n)%cable(t)%qair = q2(tid)
     cable_struc(n)%cable(t)%swdown = swd(tid)
     cable_struc(n)%cable(t)%lwdown = lwd(tid)
     cable_struc(n)%cable(t)%uwind = uwind(tid)
     cable_struc(n)%cable(t)%vwind = vwind(tid)
     ! Convert PSurf to mb
     cable_struc(n)%cable(t)%psurf = psurf(tid) / 100.0
     if (pcp(tid).ne.LIS_rc%udef) then
        cable_struc(n)%cable(t)%rainf = pcp(tid)
     else
        cable_struc(n)%cable(t)%rainf = 0.0
     endif
     if (LIS_FORC_Snowf%selectOpt.eq.1) then
        if (snowf(tid).ne.LIS_rc%udef) then
           cable_struc(n)%cable(t)%snowf = snowf(tid)
        else
           cable_struc(n)%cable(t)%snowf = 0.0
        endif
     else
        cable_struc(n)%cable(t)%snowf = 0.0
     endif
     if (LIS_FORC_CRainf%selectOpt.eq.1) then
        if (cpcp(tid).ne.LIS_rc%udef) then
           cable_struc(n)%cable(t)%rainf_c = cpcp(tid)
        else
           cable_struc(n)%cable(t)%rainf_c = 0.0
        endif
     else
        cable_struc(n)%cable(t)%rainf_c = 0.0
     endif
     ! Add an "if ..._CO2%selectOpt" line here if CO2 forcing available
     ! Convert from ppmv in lis.config file
     cable_struc(n)%cable(t)%co2 = cable_struc(n)%fixedco2/1000000.0
  enddo

end subroutine cable_f2t
