!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
subroutine lisgeosrun_coupled(n)
! 
!BOP
! !ROUTINE: lisgeosrun_coupled
! 
! !REVISION HISTORY: 
! 9 Nov 2026: Ronnie Abolafia-Rosenzweig, Initial Specification
! !USES: 
  use ESMF
  use LISGEOSGridCompMod
  use LIS_coreMod
  use LIS_metforcingMod
  use LIS_timeMgrMod
  use LIS_surfaceModelMod
  use LIS_paramsMod
  use LIS_metforcingMod
  use LIS_perturbMod
  use LIS_DAobservationsMod
  use LIS_dataAssimMod
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use LIS_irrigationMod,     only : LIS_irrigation_run ! EMK
!EOP
 
  implicit none 
! !ARGUMENTS:   
  integer, intent(in)   :: n 
!
! !DESCRIPTION: 
!  This routine defines the set of steps required from LIS during a coupled
!  LIS-GEOS run. The routine translates the GEOS import fields to LIS forcing
!  variables, and then calls the LIS LSM physics routines. Finally the 
!  export state from LIS to GEOS is generated. 
!
! The order of operations in GEOS is:
!
! 1. write output and restart data
! 2. run physics to update state
! 3. advance time
!
! Say GEOS has processed times $t_1$ and $t_2$, and that GEOS is now
! at step 1 with time $t_3$.  The above sequence of steps results in:
!
! 1. write output and restart data valid for time $t_3$
! 2. run physics to update state to time $t_4$
! 3. advance time to $t_4$
!
! To be consistent with GEOS, LIS accepts GEOS's time as is, and LIS runs
! with the same order of operations.
!
! Note that this is different from how offline LIS runs.  For an offline
! LIS run, LIS first advances time, then it runs physics to update state,
! and then it writes output.  So for an offline run, at the top of LIS'
! control loop, at time $t_3$, LIS first advances time to $t_4$.  Then
! it runs physics to update state to time $t_4$.  Then it writes output
! valid for time $t_4$.
! 
!EOP  
  integer            :: t,c,r
  integer            :: status
  type(ESMF_Field)   :: tairField, q2Field, swdField, lwdField, uwindField
  type(ESMF_Field)   :: vwindField, psurfField,rainfField,emissField,chField
  type(ESMF_Field)   :: chs2Field, cqs2Field, forchgtField, coszField, &
                        q2satField, xiceField, co2Field
  type(ESMF_Field)   :: tmnField ! EMK 
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:),snowf(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)
  real,pointer       :: ch(:),chs2(:), cqs2(:), q2sat(:),emiss(:),zval(:),&
                        cosz(:), xice(:), co2(:)
  real,pointer       :: tmn(:) ! EMK

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Tair%varname(1)),tairField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Tair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Qair%varname(1)),q2Field,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Qair')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_SWdown%varname(1)),swdField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting SWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_LWdown%varname(1)),lwdField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting LWdown')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_E%varname(1)),uwindField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Wind_E')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Wind_N%varname(1)),vwindField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Wind_N')
  
  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Psurf%varname(1)),psurfField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Psurf')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),rainfField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Rainf')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Emiss%varname(1)),emissField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Emiss')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Ch%varname(1)),chField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Ch')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Chs2%varname(1)),chs2Field,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Chs2')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Cqs2%varname(1)),cqs2Field,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Cqs2')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Q2sat%varname(1)),q2satField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting q2sat')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Forc_Hgt%varname(1)),forchgtField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Forc Hgt')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Cosz%varname(1)),coszField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Cosz')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_XICE%varname(1)),xiceField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting Xice')

  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_TMN%varname(1)),tmnField,&
       rc=status)
  call LIS_verify(status,'lisgeosrun_coupled: error getting TMN')

!  call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_CO2%varname(1)),co2Field,&
!       rc=status)
!  call LIS_verify(status,'lisgeosrun_coupled: error getting CO2')

  call ESMF_FieldGet(tairField, localDE=0, farrayPtr= tmp,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving tmp')

  call ESMF_FieldGet(q2Field, localDE=0, farrayPtr= q2,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving q2')

  call ESMF_FieldGet(swdField, localDE=0, farrayPtr= swd,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving swd')

  call ESMF_FieldGet(lwdField, localDE=0, farrayPtr= lwd,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving lwd')

  call ESMF_FieldGet(uwindField, localDE=0, farrayPtr= uwind,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving uwind')

  call ESMF_FieldGet(vwindField, localDE=0, farrayPtr= vwind,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving vwind')

  call ESMF_FieldGet(psurfField, localDE=0, farrayPtr=psurf,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving psurf')

  call ESMF_FieldGet(rainfField, localDE=0, farrayPtr=pcp,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving pcp')

  call ESMF_FieldGet(emissField, localDE=0, farrayPtr=emiss,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving emiss')

  call ESMF_FieldGet(chField, localDE=0, farrayPtr=ch,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving ch')

 call ESMF_FieldGet(chs2Field, localDE=0, farrayPtr=chs2,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving chs2')

 call ESMF_FieldGet(cqs2Field, localDE=0, farrayPtr=cqs2,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving cqs2')

  call ESMF_FieldGet(q2satField, localDE=0, farrayPtr=q2sat,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving q2sat')

  call ESMF_FieldGet(forchgtField, localDE=0, farrayPtr=zval,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving zval')

  call ESMF_FieldGet(coszField, localDE=0, farrayPtr=cosz,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving cosz')

  call ESMF_FieldGet(xiceField, localDE=0, farrayPtr=xice,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving xice')

  call ESMF_FieldGet(tmnField, localDE=0, farrayPtr=tmn,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving tmn')

!  call ESMF_FieldGet(co2Field, localDE=0, farrayPtr=co2,rc=status)
!  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving co2')

call ESMF_FieldGet(tmnField, localDE=0, farrayPtr=tmn,rc=status)
  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving tmn')

!  call ESMF_FieldGet(co2Field, localDE=0, farrayPtr=co2,rc=status)
!  call LIS_verify(status, 'lisgeosrun_coupled: error retrieving co2')

  ! ====================================================================
  ! DEFENSIVE POINTER VALIDATION
  ! ====================================================================
  ! Ensure all pointers are properly associated before dereferencing
  if (.not. associated(tmp)) then
     print *, 'ERROR lisgeosrun_coupled: tmp pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(q2)) then
     print *, 'ERROR lisgeosrun_coupled: q2 pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(swd)) then
     print *, 'ERROR lisgeosrun_coupled: swd pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(lwd)) then
     print *, 'ERROR lisgeosrun_coupled: lwd pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(uwind)) then
     print *, 'ERROR lisgeosrun_coupled: uwind pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(vwind)) then
     print *, 'ERROR lisgeosrun_coupled: vwind pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(psurf)) then
     print *, 'ERROR lisgeosrun_coupled: psurf pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(pcp)) then
     print *, 'ERROR lisgeosrun_coupled: pcp pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(emiss)) then
     print *, 'ERROR lisgeosrun_coupled: emiss pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(ch)) then
     print *, 'ERROR lisgeosrun_coupled: ch pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(chs2)) then
     print *, 'ERROR lisgeosrun_coupled: chs2 pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(cqs2)) then
     print *, 'ERROR lisgeosrun_coupled: cqs2 pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(q2sat)) then
     print *, 'ERROR lisgeosrun_coupled: q2sat pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(zval)) then
     print *, 'ERROR lisgeosrun_coupled: zval pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(cosz)) then
     print *, 'ERROR lisgeosrun_coupled: cosz pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(xice)) then
     print *, 'ERROR lisgeosrun_coupled: xice pointer not associated'
     call LIS_endrun()
  endif
  if (.not. associated(tmn)) then
     print *, 'ERROR lisgeosrun_coupled: tmn pointer not associated'
     call LIS_endrun()
  endif

  ! Verify array sizes match expected tile count
  if (size(tmp) /= LIS_rc%ntiles(n)) then
     write(LIS_logunit,*) 'ERROR: size(tmp)=', size(tmp), &
          ' but LIS_rc%ntiles(n)=', LIS_rc%ntiles(n)
     print *, 'ERROR lisgeosrun_coupled: tmp array size mismatch'
     call LIS_endrun()
  endif

  do t=1,LIS_rc%ntiles(n)
     c = LIS_domain(n)%tile(t)%col
     r = LIS_domain(n)%tile(t)%row

     tmp(t)     = LISGEOS_import(n)%data_tair(c,r)
     q2(t)      = LISGEOS_import(n)%data_qair(c,r)
     swd(t)     = LISGEOS_import(n)%data_swd(c,r)
     lwd(t)     = LISGEOS_import(n)%data_lwd(c,r)
     uwind(t)   = LISGEOS_import(n)%data_u(c,r)
     vwind(t)   = LISGEOS_import(n)%data_v(c,r)
     psurf(t)   = LISGEOS_import(n)%data_psurf(c,r)
     pcp(t)     = LISGEOS_import(n)%data_prcp(c,r) / LIS_rc%nts(n) ! pcp from 
                  ! GEOS is total in mm.  convert to rate (mm/s)

     emiss(t)   = LISGEOS_import(n)%data_emiss(c,r)
     ch(t)      = LISGEOS_import(n)%data_ch(c,r)
     chs2(t)      = LISGEOS_import(n)%data_chs2(c,r)
     cqs2(t)      = LISGEOS_import(n)%data_cqs2(c,r)
     q2sat(t)   = LISGEOS_import(n)%data_q2sat(c,r)
     zval(t)    = LISGEOS_import(n)%data_z(c,r)
     cosz(t)    = LISGEOS_import(n)%data_cosz(c,r)
     xice(t)    = LISGEOS_import(n)%data_xice(c,r)

     tmn(t)     = LISGEOS_import(n)%data_tmn(c,r) ! EMK

!     co2(t)     = LISGEOS_import(n)%data_co2(1)
!     print*, 'LISF ',t, tmp(t), q2(t), swd(t), lwd(t), uwind(t), vwind(t), &
!          psurf(t), pcp(t), emiss(t), ch(t), q2sat(t), zval(t), cosz(t)
  enddo

  call LIS_timemgr_set(LIS_rc,LISGEOS_import(n)%yr,&
       LISGEOS_import(n)%mo,LISGEOS_import(n)%da,&
       LISGEOS_import(n)%hr,LISGEOS_import(n)%mn,&
       LISGEOS_import(n)%ss,0,0.0)
!  print*, 'GEOSTIME ',LISGEOS_import%yr,LISGEOS_import%mo,LISGEOS_import%da,&
!       LISGEOS_import%hr,LISGEOS_import%mn,LISGEOS_import%ss
!  print*, 'LISTIME ',LIS_rc%yr,LIS_rc%mo,LIS_rc%da,&
!       LIS_rc%hr,LIS_rc%mn,LIS_rc%ss

  if ( .not. LISGEOS_import(n)%startflag ) then 
     call LIS_surfaceModel_output(n)
     call LIS_surfaceModel_writerestart(n)
  endif

  if(LISGEOS_import(n)%startflag) then 
     LISGEOS_import(n)%startflag = .false.
  endif

  call LIS_setDynparams(n) 
  call LIS_perturb_forcing(n)
  call LIS_irrigation_run(n) ! EMK
  call LIS_surfaceModel_f2t(n)     
  call LIS_surfaceModel_run(n)
  call LIS_surfaceModel_perturb_states(n)
  call LIS_readDAobservations(n)
  call LIS_perturb_DAobservations(n)
  call LIS_dataassim_run(n)
  call LIS_dataassim_output(n)

  call LIS_surfaceModel_setexport(n)

  flush(LIS_logunit)

end subroutine lisgeosrun_coupled
