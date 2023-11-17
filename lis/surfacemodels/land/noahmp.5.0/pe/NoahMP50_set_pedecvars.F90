!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP50_set_pedecvars
!  \label{NoahMP50_set_pedecvars}
!
! !REVISION HISTORY:
! 02 Feb 2018: Soni Yatheendradas; Initial Specification
! May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!


! !INTERFACE:
subroutine NoahMP50_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_surface
  use LIS_logMod,       only : LIS_logunit,LIS_verify
  use NoahMP50_lsmMod, only : NoahMP50_struc
  use NoahMP50_peMod,  only : NoahMP50_pe_struc

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
  integer, pointer       :: mod_flag_NoahMP50(:)
  integer                :: i,t
  integer                :: status

  n = 1

  allocate(mod_flag_NoahMP50(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_NoahMP50 = 0
    
  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,NoahMP50_pe_struc(n)%nparams
     if(NoahMP50_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(NoahMP50_pe_struc(n)%param_name(i))
        call NoahMP50_getvardata(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call NoahMP50_checkBounds(n,DEC_State,vname, vdata, mod_flag_NoahMP50)
     endif
  enddo
  deallocate(vdata) 

  !update modflags based on constraints
  call NoahMP50_checkConstraints(n,DEC_State, mod_flag_NoahMP50)

  !set variables given modflag; if flag set will leave values alone
  call NoahMP50_setVars(n,DEC_State,mod_flag_NoahMP50)

  !send mod flag to ESMF state (feasibility flag)
  call NoahMP50_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMP50)
end subroutine NoahMP50_set_pedecvars

!BOP
! 
! !ROUTINE: randArray
! \label{randArray}
!
! !INTERFACE: 
subroutine NoahMP50_getvardata(n,DEC_State,vname, vdata, statusStateGet)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_logunit,LIS_verify

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
  integer                :: statusStateGet, statusFieldGet,i
  
  call ESMF_StateGet(DEC_State,vname,varField,rc=statusStateGet)
!  call LIS_verify(status)
  
  if(statusStateGet.eq.0) then
     call ESMF_FieldGet(varField,localDE=0,farrayPtr=vardata,&
          rc=statusFieldGet)
     call LIS_verify(statusFieldGet)
     vdata=vardata
  endif
  
end subroutine NoahMP50_getvardata

subroutine NoahMP50_checkBounds(n,DEC_State,vname, vardata, mod_flag_NoahMP50)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  character*100          :: vname
  real                   :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  integer                :: mod_flag_NoahMP50(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
        mod_flag_NoahMP50(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then 
        mod_flag_NoahMP50(t) = 1
     endif
  enddo
end subroutine NoahMP50_checkBounds

subroutine NoahMP50_checkConstraints(n,DEC_State,mod_flag_NoahMP50)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_surface
  use NoahMP50_lsmMod, only : NoahMP50_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  integer                :: mod_flag_NoahMP50(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP model variables. 
! 
!EOP
  type(ESMF_Field)       :: varField
  real                   :: vardata_min, vardata_max
  character*100          :: vname
  integer                :: t
  integer                :: status1, status2
  real, allocatable      :: vardata1(:)
  real, allocatable      :: vardata2(:)
  real, allocatable      :: vardata3(:)

  allocate(vardata1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(vardata2(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(vardata3(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  vname='SMCMAX'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata1, status1)

  if(status1.ne.0) vardata1=NoahMP50_struc(n)%noahmp50(:)%param%smcmax(1)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.NoahMP50_struc(n)%noahmp50(t)%param%smcdry(1)) then
        mod_flag_NoahMP50(t) = 1
     endif
  enddo

  !SMCREF > SMCWLT
  vname='SMCREF'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCWLT'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP50_struc(n)%noahmp50(:)%param%smcref(1)
  if(status2.ne.0) vardata2=NoahMP50_struc(n)%noahmp50(:)%param%smcwlt(1)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMP50(t) = 1
     endif
  enddo

  !SMCMAX > SMCREF
  vname='SMCMAX'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCREF'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP50_struc(n)%noahmp50(:)%param%smcmax(1)
  if(status2.ne.0) vardata2=NoahMP50_struc(n)%noahmp50(:)%param%smcref(1)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMP50(t) = 1
     endif
  enddo

  !HVT > HVB
  vname='HVT'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='HVB'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP50_struc(n)%noahmp50(:)%param%HVT
  if(status2.ne.0) vardata2=NoahMP50_struc(n)%noahmp50(:)%param%HVB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).lt.vardata2(t)) then ! SY: Note .lt. instead of .le., following some entries with HVT=HVB in MPTABLE_UMD.TBL
        mod_flag_NoahMP50(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='HVT'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='Z0MVT'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP50_struc(n)%noahmp50(:)%param%HVT
  if(status2.ne.0) vardata2=NoahMP50_struc(n)%noahmp50(:)%param%Z0MVT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then 
        mod_flag_NoahMP50(t) = 1
     endif
  enddo


  !HVT > Z0MVT
  vname='MNSNALB'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='MXSNALB'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP50_struc(n)%noahmp50(:)%param%MNSNALB
  if(status2.ne.0) vardata2=NoahMP50_struc(n)%noahmp50(:)%param%MXSNALB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).ge.vardata2(t)) then 
        mod_flag_NoahMP50(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='T_ULIMIT'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='T_MLIMIT'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata2, status2)
  vname='T_LLIMIT'
  call NoahMP50_getvardata(n,DEC_State,vname,vardata3, status2)
  if(status1.ne.0) vardata1=NoahMP50_struc(n)%noahmp50(:)%param%T_ULIMIT
  if(status2.ne.0) vardata2=NoahMP50_struc(n)%noahmp50(:)%param%T_MLIMIT
  if(status2.ne.0) vardata3=NoahMP50_struc(n)%noahmp50(:)%param%T_LLIMIT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if((vardata3(t).gt.vardata2(t)).or.&
          (vardata2(t).gt.vardata1(t)).or.&
          (vardata3(t).gt.vardata1(t))) then 
        mod_flag_NoahMP50(t) = 1
     endif
  enddo

  deallocate(vardata1)
  deallocate(vardata2)
  deallocate(vardata3)

end subroutine NoahMP50_checkConstraints

subroutine NoahMP50_setVars(n,DEC_State,mod_flag_NoahMP50)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_surface
  use LIS_logMod,       only : LIS_logunit,LIS_verify
  use NoahMP50_lsmMod, only : NoahMP50_struc
  use NoahMP50_peMod,  only : NoahMP50_pe_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer                :: mod_flag_NoahMP50(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)       :: DEC_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP model variables. 
!  Only does so if the proposed parameter set is feasible (meets bounds and constraints)
! 
!EOP
  real                   :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  character*100          :: vname
  integer                :: i,t, status

  do i=1,NoahMP50_pe_struc(n)%nparams
     if(NoahMP50_pe_struc(n)%param_select(i).eq.1) then 
        vname=trim(NoahMP50_pe_struc(n)%param_name(i))
        call NoahMP50_getvardata(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if(mod_flag_NoahMP50(t).eq.0) then 
              if(vname.eq."TOPT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%topt   = vardata(t) 
              if(vname.eq."RGL") & 
                   NoahMP50_struc(n)%noahmp50(t)%param%rgl    = vardata(t) 
              if(vname.eq."RSMAX") &
                   NoahMP50_struc(n)%noahmp50(t)%param%rsmax  = vardata(t) 
              if(vname.eq."RSMIN") & 
                   NoahMP50_struc(n)%noahmp50(t)%param%rsmin  = vardata(t) 
              if(vname.eq."HS") & 
                   NoahMP50_struc(n)%noahmp50(t)%param%hs     = vardata(t) 
              if(vname.eq."NROOT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%nroot  = vardata(t) 
              if(vname.eq."CSOIL") &
                   NoahMP50_struc(n)%noahmp50(t)%param%csoil  = vardata(t) 
              if(vname.eq."BEXP") & 
                   NoahMP50_struc(n)%noahmp50(t)%param%bexp   = vardata(t) 
              if(vname.eq."DKSAT") & 
                   NoahMP50_struc(n)%noahmp50(t)%param%dksat  = vardata(t) 
              if(vname.eq."DWSAT") & 
                   NoahMP50_struc(n)%noahmp50(t)%param%dwsat  = vardata(t) 
              if(vname.eq."PSISAT") & 
                   NoahMP50_struc(n)%noahmp50(t)%param%psisat = vardata(t) 
              if(vname.eq."QUARTZ") & 
                   NoahMP50_struc(n)%noahmp50(t)%param%quartz = vardata(t) 
              if(vname.eq."SMCMAX") &
                   NoahMP50_struc(n)%noahmp50(t)%param%smcmax = vardata(t)
              if(vname.eq."SMCREF") &
                   NoahMP50_struc(n)%noahmp50(t)%param%smcref = vardata(t) 
              if(vname.eq."SMCWLT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%smcwlt = vardata(t) 
              if(vname.eq."CZIL") &
                   NoahMP50_struc(n)%noahmp50(t)%param%czil   = vardata(t) 
              if(vname.eq."SLOPE") &
                   NoahMP50_struc(n)%noahmp50(t)%param%slope  = vardata(t) 
              if(vname.eq."CH2OP") &
                   NoahMP50_struc(n)%noahmp50(t)%param%CH2OP  = vardata(t)
              if(vname.eq."DLEAF") &
                   NoahMP50_struc(n)%noahmp50(t)%param%DLEAF  = vardata(t)
              if(vname.eq."Z0MVT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%Z0MVT  = vardata(t)
              if(vname.eq."HVT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%HVT    = vardata(t)
              if(vname.eq."HVB") &
                   NoahMP50_struc(n)%noahmp50(t)%param%HVB    = vardata(t)
              if(vname.eq."RC") &
                   NoahMP50_struc(n)%noahmp50(t)%param%RC     = vardata(t)
              if(vname.eq."MFSNO") &
                   NoahMP50_struc(n)%noahmp50(t)%param%MFSNO  = vardata(t)
              if(vname.eq."ALBSAT1") &
                   NoahMP50_struc(n)%noahmp50(t)%param%ALBSAT(1) = vardata(t)
              if(vname.eq."ALBSAT2") &
                   NoahMP50_struc(n)%noahmp50(t)%param%ALBSAT(2) = vardata(t)
              if(vname.eq."ALBDRY1") &
                   NoahMP50_struc(n)%noahmp50(t)%param%ALBDRY(1) = vardata(t)
              if(vname.eq."ALBDRY2") &
                   NoahMP50_struc(n)%noahmp50(t)%param%ALBDRY(2) = vardata(t)
              if(vname.eq."ALBICE1") &
                   NoahMP50_struc(n)%noahmp50(t)%param%ALBICE(1) = vardata(t)
              if(vname.eq."ALBICE2") &
                   NoahMP50_struc(n)%noahmp50(t)%param%ALBICE(2) = vardata(t)
              if(vname.eq."OMEGAS1") &
                   NoahMP50_struc(n)%noahmp50(t)%param%OMEGAS(1) = vardata(t)
              if(vname.eq."OMEGAS2") &
                   NoahMP50_struc(n)%noahmp50(t)%param%OMEGAS(2) = vardata(t)
              if(vname.eq."BETADS") &
                   NoahMP50_struc(n)%noahmp50(t)%param%BETADS = vardata(t)
              if(vname.eq."BETAIS") &
                   NoahMP50_struc(n)%noahmp50(t)%param%BETAIS = vardata(t)
              if(vname.eq."EG1") &
                   NoahMP50_struc(n)%noahmp50(t)%param%EG(1)  = vardata(t)
              if(vname.eq."EG2") &
                   NoahMP50_struc(n)%noahmp50(t)%param%EG(2)  = vardata(t)
              if(vname.eq."Z0SNO") &
                   NoahMP50_struc(n)%noahmp50(t)%param%Z0SNO  = vardata(t)
              if(vname.eq."SSI") &
                   NoahMP50_struc(n)%noahmp50(t)%param%SSI    = vardata(t)
              if(vname.eq."SWEMX") &
                   NoahMP50_struc(n)%noahmp50(t)%param%SWEMX  = vardata(t)
              if(vname.eq."RSURF_SNOW") &
                   NoahMP50_struc(n)%noahmp50(t)%param%RSURF_SNOW = vardata(t)
              if(vname.eq."MNSNALB") &
                   NoahMP50_struc(n)%noahmp50(t)%param%MNSNALB    = vardata(t)
              if(vname.eq."MXSNALB") &
                   NoahMP50_struc(n)%noahmp50(t)%param%MXSNALB    = vardata(t)
              if(vname.eq."SNDECAYEXP") &
                   NoahMP50_struc(n)%noahmp50(t)%param%SNDECAYEXP = vardata(t)
              if(vname.eq."T_ULIMIT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%T_ULIMIT   = vardata(t)
              if(vname.eq."T_MLIMIT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%T_MLIMIT   = vardata(t)
              if(vname.eq."T_LLIMIT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%T_LLIMIT   = vardata(t)
              if(vname.eq."SNOWF_SCALEF") &
                   NoahMP50_struc(n)%noahmp50(t)%param%snowf_scalef = vardata(t)              
              if(vname.eq."RHOL1") &
                   NoahMP50_struc(n)%noahmp50(t)%param%RHOL(1) = vardata(t)
              if(vname.eq."RHOL2") &
                   NoahMP50_struc(n)%noahmp50(t)%param%RHOL(2) = vardata(t)
              if(vname.eq."RHOS1") &
                   NoahMP50_struc(n)%noahmp50(t)%param%RHOS(1) = vardata(t)
              if(vname.eq."RHOS2") &
                   NoahMP50_struc(n)%noahmp50(t)%param%RHOS(2) = vardata(t)
              if(vname.eq."TAUL1") &
                   NoahMP50_struc(n)%noahmp50(t)%param%TAUL(1) = vardata(t)
              if(vname.eq."TAUL2") &
                   NoahMP50_struc(n)%noahmp50(t)%param%TAUL(2) = vardata(t)
              if(vname.eq."TAUS1") &
                   NoahMP50_struc(n)%noahmp50(t)%param%TAUS(1) = vardata(t)
              if(vname.eq."TAUS2") &
                   NoahMP50_struc(n)%noahmp50(t)%param%TAUS(2) = vardata(t)
              if(vname.eq."XL") &
                   NoahMP50_struc(n)%noahmp50(t)%param%XL      = vardata(t)
              if(vname.eq."CWPVT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%CWPVT   = vardata(t)
              if(vname.eq."C3PSN") &
                   NoahMP50_struc(n)%noahmp50(t)%param%C3PSN   = vardata(t)
              if(vname.eq."KC25") &
                   NoahMP50_struc(n)%noahmp50(t)%param%KC25    = vardata(t)
              if(vname.eq."AKC") &
                   NoahMP50_struc(n)%noahmp50(t)%param%AKC     = vardata(t)
              if(vname.eq."KO25") &
                   NoahMP50_struc(n)%noahmp50(t)%param%KO25    = vardata(t)
              if(vname.eq."AKO") &
                   NoahMP50_struc(n)%noahmp50(t)%param%AKO     = vardata(t)
              if(vname.eq."AVCMX") &
                   NoahMP50_struc(n)%noahmp50(t)%param%AVCMX   = vardata(t)
              if(vname.eq."AQE") &
                   NoahMP50_struc(n)%noahmp50(t)%param%AQE     = vardata(t)
              if(vname.eq."LTOVRC") &
                   NoahMP50_struc(n)%noahmp50(t)%param%LTOVRC  = vardata(t)
              if(vname.eq."DILEFC") &
                   NoahMP50_struc(n)%noahmp50(t)%param%DILEFC  = vardata(t)
              if(vname.eq."DILEFW") &
                   NoahMP50_struc(n)%noahmp50(t)%param%DILEFW  = vardata(t)
              if(vname.eq."RMF25") &
                   NoahMP50_struc(n)%noahmp50(t)%param%RMF25   = vardata(t)
              if(vname.eq."SLA") &
                   NoahMP50_struc(n)%noahmp50(t)%param%SLA     = vardata(t)
              if(vname.eq."FRAGR") &
                   NoahMP50_struc(n)%noahmp50(t)%param%FRAGR   = vardata(t)
              if(vname.eq."TMIN") &
                   NoahMP50_struc(n)%noahmp50(t)%param%TMIN    = vardata(t)
              if(vname.eq."VCMX25") &
                   NoahMP50_struc(n)%noahmp50(t)%param%VCMX25  = vardata(t)
              if(vname.eq."TDLEF") &
                   NoahMP50_struc(n)%noahmp50(t)%param%TDLEF   = vardata(t)
              if(vname.eq."BP") &
                   NoahMP50_struc(n)%noahmp50(t)%param%BP      = vardata(t)
              if(vname.eq."MP") &
                   NoahMP50_struc(n)%noahmp50(t)%param%MP      = vardata(t)
              if(vname.eq."QE25") &
                   NoahMP50_struc(n)%noahmp50(t)%param%QE25    = vardata(t)
              if(vname.eq."RMS25") &
                   NoahMP50_struc(n)%noahmp50(t)%param%RMS25   = vardata(t)
              if(vname.eq."RMR25") &
                   NoahMP50_struc(n)%noahmp50(t)%param%RMR25   = vardata(t)
              if(vname.eq."ARM") &
                   NoahMP50_struc(n)%noahmp50(t)%param%ARM     = vardata(t)
              if(vname.eq."FOLNMX") &
                   NoahMP50_struc(n)%noahmp50(t)%param%FOLNMX  = vardata(t)
              if(vname.eq."WDPOOL") &
                   NoahMP50_struc(n)%noahmp50(t)%param%WDPOOL  = vardata(t)
              if(vname.eq."WRRAT") &
                   NoahMP50_struc(n)%noahmp50(t)%param%WRRAT   = vardata(t)
              if(vname.eq."MRP") &
                   NoahMP50_struc(n)%noahmp50(t)%param%MRP     = vardata(t)
           endif
        enddo
     endif
  enddo
end subroutine NoahMP50_setVars

subroutine NoahMP50_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMP50)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer                :: mod_flag_NoahMP50(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
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
     if(mod_flag_NoahMP50(t).eq.1) then 
        modflag(t)=1
     endif
  enddo

end subroutine NoahMP50_setModFlag
