!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine lis_gce_run()
  
  use ESMF
  use LIS_coreMod,         only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,      only : LIS_timemgr_set
  use LIS_FORC_AttributesMod  
  use LIS_metforcingMod,  only : LIS_FORC_state, LIS_perturb_forcing
  use LIS_surfaceModelMod, only : LIS_surfaceModel_f2t, LIS_surfaceModel_run,&
       LIS_surfaceModel_output, LIS_surfaceModel_writerestart, &
       LIS_surfaceModel_perturb_states, LIS_surfaceModel_setexport
  use LIS_paramsMod,        only : LIS_setDynparams
  use LIS_DAobservationsMod, only : LIS_readDAobservations, &
       LIS_perturb_DAobservations
  use LIS_dataAssimMod,    only : LIS_dataassim_run, LIS_dataassim_output
  use LIS_logMod,          only : LIS_logunit, LIS_verify
  use lisgceGridCompMod,   only : lisgce_import, lisgce_export

  implicit none
  
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField
  integer             :: status 
  integer             :: n , c,r, t
! only one nest
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  n = 1
  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tmpField,&
       rc=status)
  call LIS_verify(status, 'ESMF_StateGet: Tair: lis_gce_run')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status, 'ESMF_StateGet: Qair: lis_gce_run')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status, 'ESMF_StateGet: SWdown: lis_gce_run')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status, 'ESMF_StateGet: LWdown: lis_gce_run')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uField,&
       rc=status)
  call LIS_verify(status, 'ESMF_StateGet: Wind_E: lis_gce_run')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vField,&
       rc=status)
  call LIS_verify(status, 'ESMF_StateGet: Wind_N: lis_gce_run')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status, 'ESMF_StateGet: Psurf: lis_gce_run')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
       rc=status)
  call LIS_verify(status, 'ESMF_StateGet: Rainf: lis_gce_run')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CRainf%varname(1)),cpcpField,&
       rc=status)
  call LIS_verify(status, 'ESMF_StateGet: CRainf: lis_gce_run')

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LIS_verify(status, 'FieldGet: Tair: lis_gce_run')

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LIS_verify(status, 'FieldGet: Qair: lis_gce_run')
  
  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status, 'FieldGet: swd: lis_gce_run')

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status, 'FieldGet: lwd: lis_gce_run')

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LIS_verify(status, 'FieldGet: u: lis_gce_run')

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LIS_verify(status, 'FieldGet: v: lis_gce_run')
        
  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LIS_verify(status, 'FieldGet: psurf: lis_gce_run')

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LIS_verify(status, 'FieldGet: pcp: lis_gce_run')

  call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
  call LIS_verify(status, 'FieldGet: cpcp: lis_gce_run')

  
  do t=1,LIS_rc%ntiles(n)

     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row

     tmp(t)     = lisgce_import%data_tair(c,r)
     q2(t)      = lisgce_import%data_qair(c,r)
     swd(t)     = lisgce_import%data_swd(c,r)
     lwd(t)     = lisgce_import%data_lwd(c,r)
     uwind(t)   = lisgce_import%data_u(c,r)
     vwind(t)   = lisgce_import%data_v(c,r)
     psurf(t)   = lisgce_import%data_psurf(c,r)
     pcp(t)     = lisgce_import%data_prcp(c,r)
     cpcp(t)    = -9999.0
  enddo
  call LIS_timemgr_set(LIS_rc,lisgce_import%yr,&
       lisgce_import%mo,lisgce_import%da,&
       lisgce_import%hr,lisgce_import%mn,lisgce_import%ss,0,0.0) 

!map the import state to forc_state
  call LIS_setDynparams(n)
  call LIS_perturb_forcing(n)
  call LIS_surfaceModel_f2t(n)  
  call LIS_surfaceModel_run(n)
  call LIS_surfaceModel_output(n)  
  call LIS_surfaceModel_writerestart(n)
  flush(LIS_logunit)

  call LIS_surfaceModel_setexport(n)
  
  print*, 'Done LIS_RUN'
end subroutine lis_gce_run
