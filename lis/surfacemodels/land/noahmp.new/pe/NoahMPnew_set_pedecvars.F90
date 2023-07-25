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
! !ROUTINE: NoahMPnew_set_pedecvars
!  \label{NoahMPnew_set_pedecvars}
!
! !REVISION HISTORY:
! 02 Feb 2018: Soni Yatheendradas; Initial Specification
! May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!


! !INTERFACE:
subroutine NoahMPnew_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_surface
  use LIS_logMod,       only : LIS_logunit,LIS_verify
  use NoahMPnew_lsmMod, only : NoahMPnew_struc
  use NoahMPnew_peMod,  only : NoahMPnew_pe_struc

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
  integer, pointer       :: mod_flag_NoahMPnew(:)
  integer                :: i,t
  integer                :: status

  n = 1

  allocate(mod_flag_NoahMPnew(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_NoahMPnew = 0
    
  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,NoahMPnew_pe_struc(n)%nparams
     if(NoahMPnew_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(NoahMPnew_pe_struc(n)%param_name(i))
        call NoahMPnew_getvardata(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call NoahMPnew_checkBounds(n,DEC_State,vname, vdata, mod_flag_NoahMPnew)
     endif
  enddo
  deallocate(vdata) 

  !update modflags based on constraints
  call NoahMPnew_checkConstraints(n,DEC_State, mod_flag_NoahMPnew)

  !set variables given modflag; if flag set will leave values alone
  call NoahMPnew_setVars(n,DEC_State,mod_flag_NoahMPnew)

  !send mod flag to ESMF state (feasibility flag)
  call NoahMPnew_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMPnew)
end subroutine NoahMPnew_set_pedecvars

!BOP
! 
! !ROUTINE: randArray
! \label{randArray}
!
! !INTERFACE: 
subroutine NoahMPnew_getvardata(n,DEC_State,vname, vdata, statusStateGet)
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
  
end subroutine NoahMPnew_getvardata

subroutine NoahMPnew_checkBounds(n,DEC_State,vname, vardata, mod_flag_NoahMPnew)
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
  integer                :: mod_flag_NoahMPnew(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
        mod_flag_NoahMPnew(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then 
        mod_flag_NoahMPnew(t) = 1
     endif
  enddo
end subroutine NoahMPnew_checkBounds

subroutine NoahMPnew_checkConstraints(n,DEC_State,mod_flag_NoahMPnew)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_surface
  use NoahMPnew_lsmMod, only : NoahMPnew_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  integer                :: mod_flag_NoahMPnew(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata1, status1)

  if(status1.ne.0) vardata1=NoahMPnew_struc(n)%noahmpnew(:)%param%smcmax(1)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.NoahMPnew_struc(n)%noahmpnew(t)%param%smcdry(1)) then
        mod_flag_NoahMPnew(t) = 1
     endif
  enddo

  !SMCREF > SMCWLT
  vname='SMCREF'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCWLT'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMPnew_struc(n)%noahmpnew(:)%param%smcref(1)
  if(status2.ne.0) vardata2=NoahMPnew_struc(n)%noahmpnew(:)%param%smcwlt(1)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMPnew(t) = 1
     endif
  enddo

  !SMCMAX > SMCREF
  vname='SMCMAX'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCREF'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMPnew_struc(n)%noahmpnew(:)%param%smcmax(1)
  if(status2.ne.0) vardata2=NoahMPnew_struc(n)%noahmpnew(:)%param%smcref(1)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMPnew(t) = 1
     endif
  enddo

  !HVT > HVB
  vname='HVT'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='HVB'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMPnew_struc(n)%noahmpnew(:)%param%HVT
  if(status2.ne.0) vardata2=NoahMPnew_struc(n)%noahmpnew(:)%param%HVB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).lt.vardata2(t)) then ! SY: Note .lt. instead of .le., following some entries with HVT=HVB in MPTABLE_UMD.TBL
        mod_flag_NoahMPnew(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='HVT'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='Z0MVT'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMPnew_struc(n)%noahmpnew(:)%param%HVT
  if(status2.ne.0) vardata2=NoahMPnew_struc(n)%noahmpnew(:)%param%Z0MVT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then 
        mod_flag_NoahMPnew(t) = 1
     endif
  enddo


  !HVT > Z0MVT
  vname='MNSNALB'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='MXSNALB'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMPnew_struc(n)%noahmpnew(:)%param%MNSNALB
  if(status2.ne.0) vardata2=NoahMPnew_struc(n)%noahmpnew(:)%param%MXSNALB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).ge.vardata2(t)) then 
        mod_flag_NoahMPnew(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='T_ULIMIT'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='T_MLIMIT'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata2, status2)
  vname='T_LLIMIT'
  call NoahMPnew_getvardata(n,DEC_State,vname,vardata3, status2)
  if(status1.ne.0) vardata1=NoahMPnew_struc(n)%noahmpnew(:)%param%T_ULIMIT
  if(status2.ne.0) vardata2=NoahMPnew_struc(n)%noahmpnew(:)%param%T_MLIMIT
  if(status2.ne.0) vardata3=NoahMPnew_struc(n)%noahmpnew(:)%param%T_LLIMIT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if((vardata3(t).gt.vardata2(t)).or.&
          (vardata2(t).gt.vardata1(t)).or.&
          (vardata3(t).gt.vardata1(t))) then 
        mod_flag_NoahMPnew(t) = 1
     endif
  enddo

  deallocate(vardata1)
  deallocate(vardata2)
  deallocate(vardata3)

end subroutine NoahMPnew_checkConstraints

subroutine NoahMPnew_setVars(n,DEC_State,mod_flag_NoahMPnew)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_surface
  use LIS_logMod,       only : LIS_logunit,LIS_verify
  use NoahMPnew_lsmMod, only : NoahMPnew_struc
  use NoahMPnew_peMod,  only : NoahMPnew_pe_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer                :: mod_flag_NoahMPnew(LIS_rc%npatch(n,LIS_rc%lsm_index))
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

  do i=1,NoahMPnew_pe_struc(n)%nparams
     if(NoahMPnew_pe_struc(n)%param_select(i).eq.1) then 
        vname=trim(NoahMPnew_pe_struc(n)%param_name(i))
        call NoahMPnew_getvardata(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if(mod_flag_NoahMPnew(t).eq.0) then 
              if(vname.eq."TOPT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%topt   = vardata(t) 
              if(vname.eq."RGL") & 
                   NoahMPnew_struc(n)%noahmpnew(t)%param%rgl    = vardata(t) 
              if(vname.eq."RSMAX") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%rsmax  = vardata(t) 
              if(vname.eq."RSMIN") & 
                   NoahMPnew_struc(n)%noahmpnew(t)%param%rsmin  = vardata(t) 
              if(vname.eq."HS") & 
                   NoahMPnew_struc(n)%noahmpnew(t)%param%hs     = vardata(t) 
              if(vname.eq."NROOT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%nroot  = vardata(t) 
              if(vname.eq."CSOIL") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%csoil  = vardata(t) 
              if(vname.eq."BEXP") & 
                   NoahMPnew_struc(n)%noahmpnew(t)%param%bexp   = vardata(t) 
              if(vname.eq."DKSAT") & 
                   NoahMPnew_struc(n)%noahmpnew(t)%param%dksat  = vardata(t) 
              if(vname.eq."DWSAT") & 
                   NoahMPnew_struc(n)%noahmpnew(t)%param%dwsat  = vardata(t) 
              if(vname.eq."PSISAT") & 
                   NoahMPnew_struc(n)%noahmpnew(t)%param%psisat = vardata(t) 
              if(vname.eq."QUARTZ") & 
                   NoahMPnew_struc(n)%noahmpnew(t)%param%quartz = vardata(t) 
              if(vname.eq."SMCMAX") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%smcmax = vardata(t)
              if(vname.eq."SMCREF") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%smcref = vardata(t) 
              if(vname.eq."SMCWLT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%smcwlt = vardata(t) 
              if(vname.eq."CZIL") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%czil   = vardata(t) 
              if(vname.eq."SLOPE") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%slope  = vardata(t) 
              if(vname.eq."CH2OP") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%CH2OP  = vardata(t)
              if(vname.eq."DLEAF") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%DLEAF  = vardata(t)
              if(vname.eq."Z0MVT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%Z0MVT  = vardata(t)
              if(vname.eq."HVT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%HVT    = vardata(t)
              if(vname.eq."HVB") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%HVB    = vardata(t)
              if(vname.eq."RC") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%RC     = vardata(t)
              if(vname.eq."MFSNO") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%MFSNO  = vardata(t)
              if(vname.eq."ALBSAT1") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%ALBSAT(1) = vardata(t)
              if(vname.eq."ALBSAT2") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%ALBSAT(2) = vardata(t)
              if(vname.eq."ALBDRY1") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%ALBDRY(1) = vardata(t)
              if(vname.eq."ALBDRY2") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%ALBDRY(2) = vardata(t)
              if(vname.eq."ALBICE1") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%ALBICE(1) = vardata(t)
              if(vname.eq."ALBICE2") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%ALBICE(2) = vardata(t)
              if(vname.eq."OMEGAS1") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%OMEGAS(1) = vardata(t)
              if(vname.eq."OMEGAS2") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%OMEGAS(2) = vardata(t)
              if(vname.eq."BETADS") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%BETADS = vardata(t)
              if(vname.eq."BETAIS") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%BETAIS = vardata(t)
              if(vname.eq."EG1") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%EG(1)  = vardata(t)
              if(vname.eq."EG2") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%EG(2)  = vardata(t)
              if(vname.eq."Z0SNO") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%Z0SNO  = vardata(t)
              if(vname.eq."SSI") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%SSI    = vardata(t)
              if(vname.eq."SWEMX") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%SWEMX  = vardata(t)
              if(vname.eq."RSURF_SNOW") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%RSURF_SNOW = vardata(t)
              if(vname.eq."MNSNALB") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%MNSNALB    = vardata(t)
              if(vname.eq."MXSNALB") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%MXSNALB    = vardata(t)
              if(vname.eq."SNDECAYEXP") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%SNDECAYEXP = vardata(t)
              if(vname.eq."T_ULIMIT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%T_ULIMIT   = vardata(t)
              if(vname.eq."T_MLIMIT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%T_MLIMIT   = vardata(t)
              if(vname.eq."T_LLIMIT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%T_LLIMIT   = vardata(t)
              if(vname.eq."SNOWF_SCALEF") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%snowf_scalef = vardata(t)              
              if(vname.eq."RHOL1") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%RHOL(1) = vardata(t)
              if(vname.eq."RHOL2") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%RHOL(2) = vardata(t)
              if(vname.eq."RHOS1") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%RHOS(1) = vardata(t)
              if(vname.eq."RHOS2") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%RHOS(2) = vardata(t)
              if(vname.eq."TAUL1") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%TAUL(1) = vardata(t)
              if(vname.eq."TAUL2") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%TAUL(2) = vardata(t)
              if(vname.eq."TAUS1") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%TAUS(1) = vardata(t)
              if(vname.eq."TAUS2") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%TAUS(2) = vardata(t)
              if(vname.eq."XL") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%XL      = vardata(t)
              if(vname.eq."CWPVT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%CWPVT   = vardata(t)
              if(vname.eq."C3PSN") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%C3PSN   = vardata(t)
              if(vname.eq."KC25") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%KC25    = vardata(t)
              if(vname.eq."AKC") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%AKC     = vardata(t)
              if(vname.eq."KO25") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%KO25    = vardata(t)
              if(vname.eq."AKO") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%AKO     = vardata(t)
              if(vname.eq."AVCMX") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%AVCMX   = vardata(t)
              if(vname.eq."AQE") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%AQE     = vardata(t)
              if(vname.eq."LTOVRC") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%LTOVRC  = vardata(t)
              if(vname.eq."DILEFC") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%DILEFC  = vardata(t)
              if(vname.eq."DILEFW") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%DILEFW  = vardata(t)
              if(vname.eq."RMF25") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%RMF25   = vardata(t)
              if(vname.eq."SLA") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%SLA     = vardata(t)
              if(vname.eq."FRAGR") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%FRAGR   = vardata(t)
              if(vname.eq."TMIN") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%TMIN    = vardata(t)
              if(vname.eq."VCMX25") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%VCMX25  = vardata(t)
              if(vname.eq."TDLEF") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%TDLEF   = vardata(t)
              if(vname.eq."BP") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%BP      = vardata(t)
              if(vname.eq."MP") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%MP      = vardata(t)
              if(vname.eq."QE25") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%QE25    = vardata(t)
              if(vname.eq."RMS25") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%RMS25   = vardata(t)
              if(vname.eq."RMR25") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%RMR25   = vardata(t)
              if(vname.eq."ARM") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%ARM     = vardata(t)
              if(vname.eq."FOLNMX") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%FOLNMX  = vardata(t)
              if(vname.eq."WDPOOL") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%WDPOOL  = vardata(t)
              if(vname.eq."WRRAT") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%WRRAT   = vardata(t)
              if(vname.eq."MRP") &
                   NoahMPnew_struc(n)%noahmpnew(t)%param%MRP     = vardata(t)
           endif
        enddo
     endif
  enddo
end subroutine NoahMPnew_setVars

subroutine NoahMPnew_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMPnew)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer                :: mod_flag_NoahMPnew(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
     if(mod_flag_NoahMPnew(t).eq.1) then 
        modflag(t)=1
     endif
  enddo

end subroutine NoahMPnew_setModFlag
