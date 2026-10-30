!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.8
!
! Copyright (c) 2026 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP52_set_pedecvars
!  \label{NoahMP52_set_pedecvars}
!
! !REVISION HISTORY:
!  May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!  Oct 2026: David Mocko: Modified for Noah-MP-5.2
!
! !INTERFACE:
subroutine NoahMP52_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc
  use LIS_logMod,       only : LIS_verify
  use NoahMP52_lsmMod, only : NoahMP52_struc
  use NoahMP52_peMod,  only : NoahMP52_pe_struc

  implicit none
! !ARGUMENTS:
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!
!  This routine assigns the decision space to NoahMP model variables.
!
!EOP
  integer                :: n
  real, pointer          :: vdata(:)
  character*100          :: vname
  integer, pointer       :: mod_flag_NoahMP52(:)
  integer                :: i
  integer                :: status

  external :: NoahMP52_getvardata
  external :: NoahMP52_checkBounds
  external :: NoahMP52_checkConstraints
  external :: NoahMP52_setVars
  external :: NoahMP52_setModFlag

  n = 1

  allocate(mod_flag_NoahMP52(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_NoahMP52 = 0

  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,NoahMP52_pe_struc(n)%nparams
     if(NoahMP52_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(NoahMP52_pe_struc(n)%param_name(i))
        call NoahMP52_getvardata(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call NoahMP52_checkBounds(n,DEC_State,vname, vdata, mod_flag_NoahMP52)
     endif
  enddo
  deallocate(vdata)

  !update modflags based on constraints
  call NoahMP52_checkConstraints(n,DEC_State, mod_flag_NoahMP52)

  !set variables given modflag; if flag set will leave values alone
  call NoahMP52_setVars(n,DEC_State,mod_flag_NoahMP52)

  !send mod flag to ESMF state (feasibility flag)
  call NoahMP52_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMP52)
end subroutine NoahMP52_set_pedecvars

!BOP
!
! !ROUTINE: NoahMP52_getvardata
! \label{NoahMP52_getvardata}
!
! !INTERFACE:
subroutine NoahMP52_getvardata(n,DEC_State,vname, vdata, statusStateGet)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_verify

  implicit none
! !ARGUMENTS:
  integer                :: n
  type(ESMF_State)       :: DEC_State
  character*100          :: vname
  real                   :: vdata(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!
!  This routine assigns the decision space to NoahMP model variables.
!
!EOP
  real, pointer          :: vardata(:)
  type(ESMF_Field)       :: varField
  integer                :: statusStateGet, statusFieldGet

  call ESMF_StateGet(DEC_State,vname,varField,rc=statusStateGet)

  if(statusStateGet.eq.0) then
     call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
          rc=statusFieldGet)
     call LIS_verify(statusFieldGet)
     vdata=vardata
  endif

end subroutine NoahMP52_getvardata

subroutine NoahMP52_checkBounds(n,DEC_State,vname, vardata, mod_flag_NoahMP52)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_verify

  implicit none
! !ARGUMENTS:
  integer               :: n
  type(ESMF_State)      :: DEC_State
  character*100         :: vname
  real                  :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  integer               :: mod_flag_NoahMP52(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!
!  This routine assigns the decision space to NoahMP model variables.
!
!EOP
  type(ESMF_Field)       :: varField
  real                   :: vardata_min, vardata_max
  integer                :: status
  integer                :: t

  call ESMF_StateGet(DEC_State,vname,varField,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(varField,'MinRange',vardata_min,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(varField,'MaxRange',vardata_max,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata(t).lt.vardata_min) then
        mod_flag_NoahMP52(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then
        mod_flag_NoahMP52(t) = 1
     endif
  enddo
end subroutine NoahMP52_checkBounds

subroutine NoahMP52_checkConstraints(n,DEC_State,mod_flag_NoahMP52)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc
  use NoahMP52_lsmMod, only : NoahMP52_struc

  implicit none
! !ARGUMENTS:
  integer               :: n
  type(ESMF_State)      :: DEC_State
  integer               :: mod_flag_NoahMP52(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!
!  This routine assigns the decision space to NoahMP model variables.
!
!EOP

  character*100          :: vname
  integer                :: t
  integer                :: status1, status2
  real, allocatable      :: vardata1(:)
  real, allocatable      :: vardata2(:)
  real, allocatable      :: vardata3(:)

  external :: NoahMP52_getvardata

  allocate(vardata1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(vardata2(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(vardata3(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  vname='SMCMAX'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata1, status1)

  if(status1.ne.0) vardata1=NoahMP52_struc(n)%noahmp52(:)%param%smcmax(1)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.NoahMP52_struc(n)%noahmp52(t)%param%smcdry(1)) then
        mod_flag_NoahMP52(t) = 1
     endif
  enddo

  !SMCREF > SMCWLT
  vname='SMCREF'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCWLT'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP52_struc(n)%noahmp52(:)%param%smcref(1)
  if(status2.ne.0) vardata2=NoahMP52_struc(n)%noahmp52(:)%param%smcwlt(1)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMP52(t) = 1
     endif
  enddo

  !SMCMAX > SMCREF
  vname='SMCMAX'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCREF'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP52_struc(n)%noahmp52(:)%param%smcmax(1)
  if(status2.ne.0) vardata2=NoahMP52_struc(n)%noahmp52(:)%param%smcref(1)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMP52(t) = 1
     endif
  enddo

  !HVT > HVB
  vname='HVT'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='HVB'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP52_struc(n)%noahmp52(:)%param%HVT
  if(status2.ne.0) vardata2=NoahMP52_struc(n)%noahmp52(:)%param%HVB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).lt.vardata2(t)) then
        mod_flag_NoahMP52(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='HVT'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='Z0MVT'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP52_struc(n)%noahmp52(:)%param%HVT
  if(status2.ne.0) vardata2=NoahMP52_struc(n)%noahmp52(:)%param%Z0MVT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMP52(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='MNSNALB'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='MXSNALB'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP52_struc(n)%noahmp52(:)%param%MNSNALB
  if(status2.ne.0) vardata2=NoahMP52_struc(n)%noahmp52(:)%param%MXSNALB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).ge.vardata2(t)) then
        mod_flag_NoahMP52(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='T_ULIMIT'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='T_MLIMIT'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata2, status2)
  vname='T_LLIMIT'
  call NoahMP52_getvardata(n,DEC_State,vname,vardata3, status2)
  if(status1.ne.0) vardata1=NoahMP52_struc(n)%noahmp52(:)%param%T_ULIMIT
  if(status2.ne.0) vardata2=NoahMP52_struc(n)%noahmp52(:)%param%T_MLIMIT
  if(status2.ne.0) vardata3=NoahMP52_struc(n)%noahmp52(:)%param%T_LLIMIT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if((vardata3(t).gt.vardata2(t)).or.&
          (vardata2(t).gt.vardata1(t)).or.&
          (vardata3(t).gt.vardata1(t))) then
        mod_flag_NoahMP52(t) = 1
     endif
  enddo

  deallocate(vardata1)
  deallocate(vardata2)
  deallocate(vardata3)

end subroutine NoahMP52_checkConstraints

subroutine NoahMP52_setVars(n,DEC_State,mod_flag_NoahMP52)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc
  use LIS_logMod,       only : LIS_verify
  use NoahMP52_lsmMod, only : NoahMP52_struc
  use NoahMP52_peMod,  only : NoahMP52_pe_struc

  implicit none
! !ARGUMENTS:
  integer               :: n
  integer               :: mod_flag_NoahMP52(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)      :: DEC_State
!
! !DESCRIPTION:
!
!  This routine assigns the decision space to NoahMP model variables.
!  Only does so if the proposed parameter set is feasible (meets bounds and
!  constraints)
!
!EOP
  real                   :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  character*100          :: vname
  integer                :: i,t, status

  external :: NoahMP52_getvardata

  do i=1,NoahMP52_pe_struc(n)%nparams
     if(NoahMP52_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(NoahMP52_pe_struc(n)%param_name(i))
        call NoahMP52_getvardata(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if(mod_flag_NoahMP52(t).eq.0) then
              if(vname.eq."TOPT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%topt   = vardata(t)
              if(vname.eq."RGL") &
                   NoahMP52_struc(n)%noahmp52(t)%param%rgl    = vardata(t)
              if(vname.eq."RSMAX") &
                   NoahMP52_struc(n)%noahmp52(t)%param%rsmax  = vardata(t)
              if(vname.eq."RSMIN") &
                   NoahMP52_struc(n)%noahmp52(t)%param%rsmin  = vardata(t)
              if(vname.eq."HS") &
                   NoahMP52_struc(n)%noahmp52(t)%param%hs     = vardata(t)
              if(vname.eq."NROOT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%nroot  = vardata(t)
              if(vname.eq."CSOIL") &
                   NoahMP52_struc(n)%noahmp52(t)%param%csoil  = vardata(t)
              if(vname.eq."BEXP") &
                   NoahMP52_struc(n)%noahmp52(t)%param%bexp   = vardata(t)
              if(vname.eq."DKSAT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%dksat  = vardata(t)
              if(vname.eq."DWSAT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%dwsat  = vardata(t)
              if(vname.eq."PSISAT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%psisat = vardata(t)
              if(vname.eq."QUARTZ") &
                   NoahMP52_struc(n)%noahmp52(t)%param%quartz = vardata(t)
              if(vname.eq."SMCMAX") &
                   NoahMP52_struc(n)%noahmp52(t)%param%smcmax = vardata(t)
              if(vname.eq."SMCREF") &
                   NoahMP52_struc(n)%noahmp52(t)%param%smcref = vardata(t)
              if(vname.eq."SMCWLT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%smcwlt = vardata(t)
              if(vname.eq."CZIL") &
                   NoahMP52_struc(n)%noahmp52(t)%param%czil   = vardata(t)
              if(vname.eq."SLOPE") &
                   NoahMP52_struc(n)%noahmp52(t)%param%slope  = vardata(t)
              if(vname.eq."CH2OP") &
                   NoahMP52_struc(n)%noahmp52(t)%param%CH2OP  = vardata(t)
              if(vname.eq."DLEAF") &
                   NoahMP52_struc(n)%noahmp52(t)%param%DLEAF  = vardata(t)
              if(vname.eq."Z0MVT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%Z0MVT  = vardata(t)
              if(vname.eq."HVT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%HVT    = vardata(t)
              if(vname.eq."HVB") &
                   NoahMP52_struc(n)%noahmp52(t)%param%HVB    = vardata(t)
              if(vname.eq."RC") &
                   NoahMP52_struc(n)%noahmp52(t)%param%RC     = vardata(t)
              if(vname.eq."MFSNO") &
                   NoahMP52_struc(n)%noahmp52(t)%param%MFSNO  = vardata(t)
              if(vname.eq."ALBSAT1") &
                   NoahMP52_struc(n)%noahmp52(t)%param%ALBSAT(1) = vardata(t)
              if(vname.eq."ALBSAT2") &
                   NoahMP52_struc(n)%noahmp52(t)%param%ALBSAT(2) = vardata(t)
              if(vname.eq."ALBDRY1") &
                   NoahMP52_struc(n)%noahmp52(t)%param%ALBDRY(1) = vardata(t)
              if(vname.eq."ALBDRY2") &
                   NoahMP52_struc(n)%noahmp52(t)%param%ALBDRY(2) = vardata(t)
              if(vname.eq."ALBICE1") &
                   NoahMP52_struc(n)%noahmp52(t)%param%ALBICE(1) = vardata(t)
              if(vname.eq."ALBICE2") &
                   NoahMP52_struc(n)%noahmp52(t)%param%ALBICE(2) = vardata(t)
              if(vname.eq."OMEGAS1") &
                   NoahMP52_struc(n)%noahmp52(t)%param%OMEGAS(1) = vardata(t)
              if(vname.eq."OMEGAS2") &
                   NoahMP52_struc(n)%noahmp52(t)%param%OMEGAS(2) = vardata(t)
              if(vname.eq."BETADS") &
                   NoahMP52_struc(n)%noahmp52(t)%param%BETADS = vardata(t)
              if(vname.eq."BETAIS") &
                   NoahMP52_struc(n)%noahmp52(t)%param%BETAIS = vardata(t)
              if(vname.eq."EG1") &
                   NoahMP52_struc(n)%noahmp52(t)%param%EG(1)  = vardata(t)
              if(vname.eq."EG2") &
                   NoahMP52_struc(n)%noahmp52(t)%param%EG(2)  = vardata(t)
              if(vname.eq."Z0SNO") &
                   NoahMP52_struc(n)%noahmp52(t)%param%Z0SNO  = vardata(t)
              if(vname.eq."SSI") &
                   NoahMP52_struc(n)%noahmp52(t)%param%SSI    = vardata(t)
              if(vname.eq."SWEMX") &
                   NoahMP52_struc(n)%noahmp52(t)%param%SWEMX  = vardata(t)
              if(vname.eq."RSURF_SNOW") &
                   NoahMP52_struc(n)%noahmp52(t)%param%RSURF_SNOW = vardata(t)
              if(vname.eq."MNSNALB") &
                   NoahMP52_struc(n)%noahmp52(t)%param%MNSNALB    = vardata(t)
              if(vname.eq."MXSNALB") &
                   NoahMP52_struc(n)%noahmp52(t)%param%MXSNALB    = vardata(t)
              if(vname.eq."SNDECAYEXP") &
                   NoahMP52_struc(n)%noahmp52(t)%param%SNDECAYEXP = vardata(t)
              if(vname.eq."T_ULIMIT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%T_ULIMIT   = vardata(t)
              if(vname.eq."T_MLIMIT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%T_MLIMIT   = vardata(t)
              if(vname.eq."T_LLIMIT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%T_LLIMIT   = vardata(t)
              if(vname.eq."SNOWF_SCALEF") &
                   NoahMP52_struc(n)%noahmp52(t)%param%snowf_scalef = &
                   vardata(t)
              if(vname.eq."RHOL1") &
                   NoahMP52_struc(n)%noahmp52(t)%param%RHOL(1) = vardata(t)
              if(vname.eq."RHOL2") &
                   NoahMP52_struc(n)%noahmp52(t)%param%RHOL(2) = vardata(t)
              if(vname.eq."RHOS1") &
                   NoahMP52_struc(n)%noahmp52(t)%param%RHOS(1) = vardata(t)
              if(vname.eq."RHOS2") &
                   NoahMP52_struc(n)%noahmp52(t)%param%RHOS(2) = vardata(t)
              if(vname.eq."TAUL1") &
                   NoahMP52_struc(n)%noahmp52(t)%param%TAUL(1) = vardata(t)
              if(vname.eq."TAUL2") &
                   NoahMP52_struc(n)%noahmp52(t)%param%TAUL(2) = vardata(t)
              if(vname.eq."TAUS1") &
                   NoahMP52_struc(n)%noahmp52(t)%param%TAUS(1) = vardata(t)
              if(vname.eq."TAUS2") &
                   NoahMP52_struc(n)%noahmp52(t)%param%TAUS(2) = vardata(t)
              if(vname.eq."XL") &
                   NoahMP52_struc(n)%noahmp52(t)%param%XL      = vardata(t)
              if(vname.eq."CWPVT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%CWPVT   = vardata(t)
              if(vname.eq."C3PSN") &
                   NoahMP52_struc(n)%noahmp52(t)%param%C3PSN   = vardata(t)
              if(vname.eq."KC25") &
                   NoahMP52_struc(n)%noahmp52(t)%param%KC25    = vardata(t)
              if(vname.eq."AKC") &
                   NoahMP52_struc(n)%noahmp52(t)%param%AKC     = vardata(t)
              if(vname.eq."KO25") &
                   NoahMP52_struc(n)%noahmp52(t)%param%KO25    = vardata(t)
              if(vname.eq."AKO") &
                   NoahMP52_struc(n)%noahmp52(t)%param%AKO     = vardata(t)
              if(vname.eq."AVCMX") &
                   NoahMP52_struc(n)%noahmp52(t)%param%AVCMX   = vardata(t)
              if(vname.eq."AQE") &
                   NoahMP52_struc(n)%noahmp52(t)%param%AQE     = vardata(t)
              if(vname.eq."LTOVRC") &
                   NoahMP52_struc(n)%noahmp52(t)%param%LTOVRC  = vardata(t)
              if(vname.eq."DILEFC") &
                   NoahMP52_struc(n)%noahmp52(t)%param%DILEFC  = vardata(t)
              if(vname.eq."DILEFW") &
                   NoahMP52_struc(n)%noahmp52(t)%param%DILEFW  = vardata(t)
              if(vname.eq."RMF25") &
                   NoahMP52_struc(n)%noahmp52(t)%param%RMF25   = vardata(t)
              if(vname.eq."SLA") &
                   NoahMP52_struc(n)%noahmp52(t)%param%SLA     = vardata(t)
              if(vname.eq."FRAGR") &
                   NoahMP52_struc(n)%noahmp52(t)%param%FRAGR   = vardata(t)
              if(vname.eq."TMIN") &
                   NoahMP52_struc(n)%noahmp52(t)%param%TMIN    = vardata(t)
              if(vname.eq."VCMX25") &
                   NoahMP52_struc(n)%noahmp52(t)%param%VCMX25  = vardata(t)
              if(vname.eq."TDLEF") &
                   NoahMP52_struc(n)%noahmp52(t)%param%TDLEF   = vardata(t)
              if(vname.eq."BP") &
                   NoahMP52_struc(n)%noahmp52(t)%param%BP      = vardata(t)
              if(vname.eq."MP") &
                   NoahMP52_struc(n)%noahmp52(t)%param%MP      = vardata(t)
              if(vname.eq."QE25") &
                   NoahMP52_struc(n)%noahmp52(t)%param%QE25    = vardata(t)
              if(vname.eq."RMS25") &
                   NoahMP52_struc(n)%noahmp52(t)%param%RMS25   = vardata(t)
              if(vname.eq."RMR25") &
                   NoahMP52_struc(n)%noahmp52(t)%param%RMR25   = vardata(t)
              if(vname.eq."ARM") &
                   NoahMP52_struc(n)%noahmp52(t)%param%ARM     = vardata(t)
              if(vname.eq."FOLNMX") &
                   NoahMP52_struc(n)%noahmp52(t)%param%FOLNMX  = vardata(t)
              if(vname.eq."WDPOOL") &
                   NoahMP52_struc(n)%noahmp52(t)%param%WDPOOL  = vardata(t)
              if(vname.eq."WRRAT") &
                   NoahMP52_struc(n)%noahmp52(t)%param%WRRAT   = vardata(t)
              if(vname.eq."MRP") &
                   NoahMP52_struc(n)%noahmp52(t)%param%MRP     = vardata(t)
           endif
        enddo
     endif
  enddo
end subroutine NoahMP52_setVars

subroutine NoahMP52_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMP52)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_verify

  implicit none
! !ARGUMENTS:
  integer               :: n
  integer               :: mod_flag_NoahMP52(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)      :: DEC_State
  type(ESMF_State)      :: Feas_State
!
! !DESCRIPTION:
!
!  This routine sets the feasibility flag
!
!EOP
  type(ESMF_Field)       :: feasField
  integer                :: t
  integer                :: status
  integer, pointer       :: modflag(:)

  call ESMF_StateGet(Feas_State, "Feasibility Flag", feasField, rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(feasField,localDE=0,farrayPtr=modflag,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(mod_flag_NoahMP52(t).eq.1) then
        modflag(t)=1
     endif
  enddo

end subroutine NoahMP52_setModFlag
