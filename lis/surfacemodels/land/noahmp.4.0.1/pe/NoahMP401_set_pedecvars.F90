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
! !ROUTINE: NoahMP401_set_pedecvars
!  \label{NoahMP401_set_pedecvars}
!
! !REVISION HISTORY:
! 02 Feb 2018: Soni Yatheendradas; Initial Specification
!


! !INTERFACE:
subroutine NoahMP401_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use NoahMP401_lsmMod, only : NoahMP401_struc
  use NoahMP401_peMod,  only : NoahMP401_pe_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
! 
!EOP
  integer                :: n
  real, pointer          :: vdata(:)
  character*100          :: vname
  integer, pointer       :: mod_flag_NoahMP401(:)
  integer                :: i,t
  integer                :: status

  n = 1

  allocate(mod_flag_NoahMP401(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_NoahMP401 = 0
    
  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,NoahMP401_pe_struc(n)%nparams
     if(NoahMP401_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(NoahMP401_pe_struc(n)%param_name(i))
        call NoahMP401_getvardata(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call NoahMP401_checkBounds(n,DEC_State,vname, vdata, mod_flag_NoahMP401)
     endif
  enddo
  deallocate(vdata) 

  !update modflags based on constraints
  call NoahMP401_checkConstraints(n,DEC_State, mod_flag_NoahMP401)

  !set variables given modflag; if flag set will leave values alone
  call NoahMP401_setVars(n,DEC_State,mod_flag_NoahMP401)

  !send mod flag to ESMF state (feasibility flag)
  call NoahMP401_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMP401)
end subroutine NoahMP401_set_pedecvars

!BOP
! 
! !ROUTINE: randArray
! \label{randArray}
!
! !INTERFACE: 
subroutine NoahMP401_getvardata(n,DEC_State,vname, vdata, statusStateGet)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  character*100          :: vname
  real          :: vdata(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
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
  
end subroutine NoahMP401_getvardata

subroutine NoahMP401_checkBounds(n,DEC_State,vname, vardata, mod_flag_NoahMP401)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  character*100          :: vname
  real          :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  integer       :: mod_flag_NoahMP401(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
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
        mod_flag_NoahMP401(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then 
        mod_flag_NoahMP401(t) = 1
     endif
  enddo
end subroutine NoahMP401_checkBounds

subroutine NoahMP401_checkConstraints(n,DEC_State,mod_flag_NoahMP401)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use NoahMP401_lsmMod, only : NoahMP401_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  integer       :: mod_flag_NoahMP401(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
! 
!EOP
  type(ESMF_Field)       :: varField
  real                   :: vardata_min, vardata_max
  character*100          :: vname
  integer                :: t
  integer       :: status1, status2
  real, allocatable :: vardata1(:)
  real, allocatable :: vardata2(:)
  real, allocatable :: vardata3(:)

  allocate(vardata1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(vardata2(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(vardata3(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  vname='SMCMAX'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata1, status1)

  if(status1.ne.0) vardata1=NoahMP401_struc(n)%noahmp401(:)%param%smcmax(1)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.NoahMP401_struc(n)%noahmp401(t)%param%smcdry(1)) then
        mod_flag_NoahMP401(t) = 1
     endif
  enddo

  !SMCREF > SMCWLT
  vname='SMCREF'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCWLT'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP401_struc(n)%noahmp401(:)%param%smcref(1)
  if(status2.ne.0) vardata2=NoahMP401_struc(n)%noahmp401(:)%param%smcwlt(1)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMP401(t) = 1
     endif
  enddo

  !SMCMAX > SMCREF
  vname='SMCMAX'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCREF'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP401_struc(n)%noahmp401(:)%param%smcmax(1)
  if(status2.ne.0) vardata2=NoahMP401_struc(n)%noahmp401(:)%param%smcref(1)
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMP401(t) = 1
     endif
  enddo

  !HVT > HVB
  vname='HVT'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='HVB'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP401_struc(n)%noahmp401(:)%param%HVT
  if(status2.ne.0) vardata2=NoahMP401_struc(n)%noahmp401(:)%param%HVB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).lt.vardata2(t)) then ! SY: Note .lt. instead of .le., following some entries with HVT=HVB in MPTABLE_UMD.TBL
        mod_flag_NoahMP401(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='HVT'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='Z0MVT'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP401_struc(n)%noahmp401(:)%param%HVT
  if(status2.ne.0) vardata2=NoahMP401_struc(n)%noahmp401(:)%param%Z0MVT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then 
        mod_flag_NoahMP401(t) = 1
     endif
  enddo


  !HVT > Z0MVT
  vname='MNSNALB'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='MXSNALB'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP401_struc(n)%noahmp401(:)%param%MNSNALB
  if(status2.ne.0) vardata2=NoahMP401_struc(n)%noahmp401(:)%param%MXSNALB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).ge.vardata2(t)) then 
        mod_flag_NoahMP401(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='T_ULIMIT'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='T_MLIMIT'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata2, status2)
  vname='T_LLIMIT'
  call NoahMP401_getvardata(n,DEC_State,vname,vardata3, status2)
  if(status1.ne.0) vardata1=NoahMP401_struc(n)%noahmp401(:)%param%T_ULIMIT
  if(status2.ne.0) vardata2=NoahMP401_struc(n)%noahmp401(:)%param%T_MLIMIT
  if(status2.ne.0) vardata3=NoahMP401_struc(n)%noahmp401(:)%param%T_LLIMIT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if((vardata3(t).gt.vardata2(t)).or.&
          (vardata2(t).gt.vardata1(t)).or.&
          (vardata3(t).gt.vardata1(t))) then 
        mod_flag_NoahMP401(t) = 1
     endif
  enddo

  deallocate(vardata1)
  deallocate(vardata2)
  deallocate(vardata3)

end subroutine NoahMP401_checkConstraints

subroutine NoahMP401_setVars(n,DEC_State,mod_flag_NoahMP401)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use NoahMP401_lsmMod, only : NoahMP401_struc
  use NoahMP401_peMod,  only : NoahMP401_pe_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_NoahMP401(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)       :: DEC_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to NoahMP3.6 model variables. 
!  Only does so if the proposed parameter set is feasible (meets bounds and constraints)
! 
!EOP
  real          :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  character*100          :: vname
  integer                :: i,t, status

  do i=1,NoahMP401_pe_struc(n)%nparams
     if(NoahMP401_pe_struc(n)%param_select(i).eq.1) then 
        vname=trim(NoahMP401_pe_struc(n)%param_name(i))
        call NoahMP401_getvardata(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if(mod_flag_NoahMP401(t).eq.0) then 
              if(vname.eq."TOPT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%topt = vardata(t) 
              if(vname.eq."RGL") & 
                   NoahMP401_struc(n)%noahmp401(t)%param%rgl = vardata(t) 
              if(vname.eq."RSMAX") &
                   NoahMP401_struc(n)%noahmp401(t)%param%rsmax = vardata(t) 
              if(vname.eq."RSMIN") & 
                   NoahMP401_struc(n)%noahmp401(t)%param%rsmin = vardata(t) 
              if(vname.eq."HS") & 
                   NoahMP401_struc(n)%noahmp401(t)%param%hs = vardata(t) 
              if(vname.eq."NROOT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%nroot = vardata(t) 
              if(vname.eq."CSOIL") &
                   NoahMP401_struc(n)%noahmp401(t)%param%csoil = vardata(t) 
              if(vname.eq."BEXP") & 
                   NoahMP401_struc(n)%noahmp401(t)%param%bexp = vardata(t) 
              if(vname.eq."DKSAT") & 
                   NoahMP401_struc(n)%noahmp401(t)%param%dksat = vardata(t) 
              if(vname.eq."DWSAT") & 
                   NoahMP401_struc(n)%noahmp401(t)%param%dwsat = vardata(t) 
              if(vname.eq."PSISAT") & 
                   NoahMP401_struc(n)%noahmp401(t)%param%psisat = vardata(t) 
              if(vname.eq."QUARTZ") & 
                   NoahMP401_struc(n)%noahmp401(t)%param%quartz = vardata(t) 
              if(vname.eq."SMCMAX") &
                   NoahMP401_struc(n)%noahmp401(t)%param%smcmax = vardata(t)
              if(vname.eq."SMCREF") &
                   NoahMP401_struc(n)%noahmp401(t)%param%smcref = vardata(t) 
              if(vname.eq."SMCWLT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%smcwlt = vardata(t) 
              if(vname.eq."CZIL") &
                   NoahMP401_struc(n)%noahmp401(t)%param%czil = vardata(t) 
              if(vname.eq."SLOPE") &
                   NoahMP401_struc(n)%noahmp401(t)%param%slope = vardata(t) 
              if(vname.eq."CH2OP") &
                   NoahMP401_struc(n)%noahmp401(t)%param%CH2OP = vardata(t)
              if(vname.eq."DLEAF") &
                   NoahMP401_struc(n)%noahmp401(t)%param%DLEAF = vardata(t)
              if(vname.eq."Z0MVT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%Z0MVT = vardata(t)
              if(vname.eq."HVT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%HVT = vardata(t)
              if(vname.eq."HVB") &
                   NoahMP401_struc(n)%noahmp401(t)%param%HVB = vardata(t)
              if(vname.eq."RC") &
                   NoahMP401_struc(n)%noahmp401(t)%param%RC = vardata(t)
              if(vname.eq."MFSNO") &
                   NoahMP401_struc(n)%noahmp401(t)%param%MFSNO = vardata(t)
              if(vname.eq."ALBSAT1") &
                   NoahMP401_struc(n)%noahmp401(t)%param%ALBSAT(1) = vardata(t)
              if(vname.eq."ALBSAT2") &
                   NoahMP401_struc(n)%noahmp401(t)%param%ALBSAT(2) = vardata(t)
              if(vname.eq."ALBDRY1") &
                   NoahMP401_struc(n)%noahmp401(t)%param%ALBDRY(1) = vardata(t)
              if(vname.eq."ALBDRY2") &
                   NoahMP401_struc(n)%noahmp401(t)%param%ALBDRY(2) = vardata(t)
              if(vname.eq."ALBICE1") &
                   NoahMP401_struc(n)%noahmp401(t)%param%ALBICE(1) = vardata(t)
              if(vname.eq."ALBICE2") &
                   NoahMP401_struc(n)%noahmp401(t)%param%ALBICE(2) = vardata(t)
              if(vname.eq."OMEGAS1") &
                   NoahMP401_struc(n)%noahmp401(t)%param%OMEGAS(1) = vardata(t)
              if(vname.eq."OMEGAS2") &
                   NoahMP401_struc(n)%noahmp401(t)%param%OMEGAS(2) = vardata(t)
              if(vname.eq."BETADS") &
                   NoahMP401_struc(n)%noahmp401(t)%param%BETADS = vardata(t)
              if(vname.eq."BETAIS") &
                   NoahMP401_struc(n)%noahmp401(t)%param%BETAIS = vardata(t)
              if(vname.eq."EG1") &
                   NoahMP401_struc(n)%noahmp401(t)%param%EG(1) = vardata(t)
              if(vname.eq."EG2") &
                   NoahMP401_struc(n)%noahmp401(t)%param%EG(2) = vardata(t)
              if(vname.eq."Z0SNO") &
                   NoahMP401_struc(n)%noahmp401(t)%param%Z0SNO = vardata(t)
              if(vname.eq."SSI") &
                   NoahMP401_struc(n)%noahmp401(t)%param%SSI = vardata(t)
              if(vname.eq."SWEMX") &
                   NoahMP401_struc(n)%noahmp401(t)%param%SWEMX = vardata(t)
              if(vname.eq."RSURF_SNOW") &
                   NoahMP401_struc(n)%noahmp401(t)%param%RSURF_SNOW = vardata(t)
              if(vname.eq."MNSNALB") &
                   NoahMP401_struc(n)%noahmp401(t)%param%MNSNALB = vardata(t)
              if(vname.eq."MXSNALB") &
                   NoahMP401_struc(n)%noahmp401(t)%param%MXSNALB= vardata(t)
              if(vname.eq."SNDECAYEXP") &
                   NoahMP401_struc(n)%noahmp401(t)%param%SNDECAYEXP = vardata(t)
              if(vname.eq."T_ULIMIT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%T_ULIMIT = vardata(t)
              if(vname.eq."T_MLIMIT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%T_MLIMIT = vardata(t)
              if(vname.eq."T_LLIMIT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%T_LLIMIT = vardata(t)
              if(vname.eq."SNOWF_SCALEF") &
                   NoahMP401_struc(n)%noahmp401(t)%param%snowf_scalef = vardata(t)              
              if(vname.eq."RHOL1") &
                   NoahMP401_struc(n)%noahmp401(t)%param%RHOL(1) = vardata(t)
              if(vname.eq."RHOL2") &
                   NoahMP401_struc(n)%noahmp401(t)%param%RHOL(2) = vardata(t)
              if(vname.eq."RHOS1") &
                   NoahMP401_struc(n)%noahmp401(t)%param%RHOS(1) = vardata(t)
              if(vname.eq."RHOS2") &
                   NoahMP401_struc(n)%noahmp401(t)%param%RHOS(2) = vardata(t)
              if(vname.eq."TAUL1") &
                   NoahMP401_struc(n)%noahmp401(t)%param%TAUL(1) = vardata(t)
              if(vname.eq."TAUL2") &
                   NoahMP401_struc(n)%noahmp401(t)%param%TAUL(2) = vardata(t)
              if(vname.eq."TAUS1") &
                   NoahMP401_struc(n)%noahmp401(t)%param%TAUS(1) = vardata(t)
              if(vname.eq."TAUS2") &
                   NoahMP401_struc(n)%noahmp401(t)%param%TAUS(2) = vardata(t)
              if(vname.eq."XL") &
                   NoahMP401_struc(n)%noahmp401(t)%param%XL = vardata(t)
              if(vname.eq."CWPVT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%CWPVT = vardata(t)
              if(vname.eq."C3PSN") &
                   NoahMP401_struc(n)%noahmp401(t)%param%C3PSN = vardata(t)
              if(vname.eq."KC25") &
                   NoahMP401_struc(n)%noahmp401(t)%param%KC25 = vardata(t)
              if(vname.eq."AKC") &
                   NoahMP401_struc(n)%noahmp401(t)%param%AKC = vardata(t)
              if(vname.eq."KO25") &
                   NoahMP401_struc(n)%noahmp401(t)%param%KO25 = vardata(t)
              if(vname.eq."AKO") &
                   NoahMP401_struc(n)%noahmp401(t)%param%AKO = vardata(t)
              if(vname.eq."AVCMX") &
                   NoahMP401_struc(n)%noahmp401(t)%param%AVCMX = vardata(t)
              if(vname.eq."AQE") &
                   NoahMP401_struc(n)%noahmp401(t)%param%AQE = vardata(t)
              if(vname.eq."LTOVRC") &
                   NoahMP401_struc(n)%noahmp401(t)%param%LTOVRC = vardata(t)
              if(vname.eq."DILEFC") &
                   NoahMP401_struc(n)%noahmp401(t)%param%DILEFC = vardata(t)
              if(vname.eq."DILEFW") &
                   NoahMP401_struc(n)%noahmp401(t)%param%DILEFW = vardata(t)
              if(vname.eq."RMF25") &
                   NoahMP401_struc(n)%noahmp401(t)%param%RMF25 = vardata(t)
              if(vname.eq."SLA") &
                   NoahMP401_struc(n)%noahmp401(t)%param%SLA = vardata(t)
              if(vname.eq."FRAGR") &
                   NoahMP401_struc(n)%noahmp401(t)%param%FRAGR = vardata(t)
              if(vname.eq."TMIN") &
                   NoahMP401_struc(n)%noahmp401(t)%param%TMIN = vardata(t)
              if(vname.eq."VCMX25") &
                   NoahMP401_struc(n)%noahmp401(t)%param%VCMX25 = vardata(t)
              if(vname.eq."TDLEF") &
                   NoahMP401_struc(n)%noahmp401(t)%param%TDLEF = vardata(t)
              if(vname.eq."BP") &
                   NoahMP401_struc(n)%noahmp401(t)%param%BP = vardata(t)
              if(vname.eq."MP") &
                   NoahMP401_struc(n)%noahmp401(t)%param%MP = vardata(t)
              if(vname.eq."QE25") &
                   NoahMP401_struc(n)%noahmp401(t)%param%QE25 = vardata(t)
              if(vname.eq."RMS25") &
                   NoahMP401_struc(n)%noahmp401(t)%param%RMS25 = vardata(t)
              if(vname.eq."RMR25") &
                   NoahMP401_struc(n)%noahmp401(t)%param%RMR25 = vardata(t)
              if(vname.eq."ARM") &
                   NoahMP401_struc(n)%noahmp401(t)%param%ARM = vardata(t)
              if(vname.eq."FOLNMX") &
                   NoahMP401_struc(n)%noahmp401(t)%param%FOLNMX = vardata(t)
              if(vname.eq."WDPOOL") &
                   NoahMP401_struc(n)%noahmp401(t)%param%WDPOOL = vardata(t)
              if(vname.eq."WRRAT") &
                   NoahMP401_struc(n)%noahmp401(t)%param%WRRAT = vardata(t)
              if(vname.eq."MRP") &
                   NoahMP401_struc(n)%noahmp401(t)%param%MRP = vardata(t)
           endif
        enddo
     endif
  enddo
end subroutine NoahMP401_setVars

subroutine NoahMP401_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMP401)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_NoahMP401(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
     if(mod_flag_NoahMP401(t).eq.1) then 
        modflag(t)=1
     endif
  enddo

end subroutine NoahMP401_setModFlag
