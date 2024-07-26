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
! !ROUTINE: NoahMP36_set_pedecvars
!  \label{NoahMP36_set_pedecvars}
!
! !REVISION HISTORY:
! 02 Feb 2018: Soni Yatheendradas; Initial Specification
!


! !INTERFACE:
subroutine NoahMP36_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use NoahMP36_lsmMod, only : NoahMP36_struc
  use NoahMP36_peMod,  only : NoahMP36_pe_struc

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
  integer, pointer       :: mod_flag_NoahMP36(:)
  integer                :: i,t
  integer                :: status

  n = 1

  allocate(mod_flag_NoahMP36(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_NoahMP36 = 0
    
  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,NoahMP36_pe_struc(n)%nparams
     if(NoahMP36_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(NoahMP36_pe_struc(n)%param_name(i))
        call NoahMP36_getvardata(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call NoahMP36_checkBounds(n,DEC_State,vname, vdata, mod_flag_NoahMP36)
     endif
  enddo
  deallocate(vdata) 

  !update modflags based on constraints
  call NoahMP36_checkConstraints(n,DEC_State, mod_flag_NoahMP36)

  !set variables given modflag; if flag set will leave values alone
  call NoahMP36_setVars(n,DEC_State,mod_flag_NoahMP36)

  !send mod flag to ESMF state (feasibility flag)
  call NoahMP36_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMP36)
end subroutine NoahMP36_set_pedecvars

!BOP
! 
! !ROUTINE: randArray
! \label{randArray}
!
! !INTERFACE: 
subroutine NoahMP36_getvardata(n,DEC_State,vname, vdata, statusStateGet)
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
  
end subroutine NoahMP36_getvardata

subroutine NoahMP36_checkBounds(n,DEC_State,vname, vardata, mod_flag_NoahMP36)
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
  integer       :: mod_flag_NoahMP36(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
        mod_flag_NoahMP36(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then 
        mod_flag_NoahMP36(t) = 1
     endif
  enddo
end subroutine NoahMP36_checkBounds

subroutine NoahMP36_checkConstraints(n,DEC_State,mod_flag_NoahMP36)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use NoahMP36_lsmMod, only : NoahMP36_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  integer       :: mod_flag_NoahMP36(LIS_rc%npatch(n,LIS_rc%lsm_index))
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

  allocate(vardata1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(vardata2(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  !SMCMAX > SMCDRY
  vname='SMCMAX'
  call NoahMP36_getvardata(n,DEC_State,vname,vardata1, status1)
!  vname='SMCDRY'
!  call NoahMP36_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP36_struc(n)%noahmp36(:)%smcmax
!  if(status2.ne.0) vardata2=NoahMP36_struc(n)%noahmp36(:)%smcdry
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     if(vardata1(t).le.vardata2(t)) then
     if(vardata1(t).le.NoahMP36_struc(n)%noahmp36(t)%smcdry) then
        mod_flag_NoahMP36(t) = 1
     endif
  enddo

  !SMCREF > SMCWLT
  vname='SMCREF'
  call NoahMP36_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCWLT'
  call NoahMP36_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP36_struc(n)%noahmp36(:)%smcref
  if(status2.ne.0) vardata2=NoahMP36_struc(n)%noahmp36(:)%smcwlt
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMP36(t) = 1
     endif
  enddo

  !SMCMAX > SMCREF
  vname='SMCMAX'
  call NoahMP36_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCREF'
  call NoahMP36_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP36_struc(n)%noahmp36(:)%smcmax
  if(status2.ne.0) vardata2=NoahMP36_struc(n)%noahmp36(:)%smcref
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_NoahMP36(t) = 1
     endif
  enddo

  !HVT > HVB
  vname='HVT'
  call NoahMP36_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='HVB'
  call NoahMP36_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP36_struc(n)%noahmp36(:)%HVT
  if(status2.ne.0) vardata2=NoahMP36_struc(n)%noahmp36(:)%HVB
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).lt.vardata2(t)) then ! SY: Note .lt. instead of .le., following some entries with HVT=HVB in MPTABLE_UMD.TBL
        mod_flag_NoahMP36(t) = 1
     endif
  enddo

  !HVT > Z0MVT
  vname='HVT'
  call NoahMP36_getvardata(n,DEC_State,vname,vardata1, status1)
  vname='Z0MVT'
  call NoahMP36_getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=NoahMP36_struc(n)%noahmp36(:)%HVT
  if(status2.ne.0) vardata2=NoahMP36_struc(n)%noahmp36(:)%Z0MVT
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then 
        mod_flag_NoahMP36(t) = 1
     endif
  enddo

  deallocate(vardata1)
  deallocate(vardata2)

end subroutine NoahMP36_checkConstraints

subroutine NoahMP36_setVars(n,DEC_State,mod_flag_NoahMP36)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use NoahMP36_lsmMod, only : NoahMP36_struc
  use NoahMP36_peMod,  only : NoahMP36_pe_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_NoahMP36(LIS_rc%npatch(n,LIS_rc%lsm_index))
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

  do i=1,NoahMP36_pe_struc(n)%nparams
     if(NoahMP36_pe_struc(n)%param_select(i).eq.1) then 
        vname=trim(NoahMP36_pe_struc(n)%param_name(i))
        call NoahMP36_getvardata(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if(mod_flag_NoahMP36(t).eq.0) then 
              if(vname.eq."TOPT") &
                   NoahMP36_struc(n)%noahmp36(t)%topt = vardata(t) 
              if(vname.eq."RGL") & 
                   NoahMP36_struc(n)%noahmp36(t)%rgl = vardata(t) 
              if(vname.eq."RSMAX") &
                   NoahMP36_struc(n)%noahmp36(t)%rsmax = vardata(t) 
              if(vname.eq."RSMIN") & 
                   NoahMP36_struc(n)%noahmp36(t)%rsmin = vardata(t) 
              if(vname.eq."HS") & 
                   NoahMP36_struc(n)%noahmp36(t)%hs = vardata(t) 
              if(vname.eq."NROOT") &
                   NoahMP36_struc(n)%noahmp36(t)%nroot = vardata(t) 
              if(vname.eq."CSOIL") &
                   NoahMP36_struc(n)%noahmp36(t)%csoil = vardata(t) 
              if(vname.eq."BEXP") & 
                   NoahMP36_struc(n)%noahmp36(t)%bexp = vardata(t) 
              if(vname.eq."DKSAT") & 
                   NoahMP36_struc(n)%noahmp36(t)%dksat = vardata(t) 
              if(vname.eq."DWSAT") & 
                   NoahMP36_struc(n)%noahmp36(t)%dwsat = vardata(t) 
              if(vname.eq."PSISAT") & 
                   NoahMP36_struc(n)%noahmp36(t)%psisat = vardata(t) 
              if(vname.eq."QUARTZ") & 
                   NoahMP36_struc(n)%noahmp36(t)%quartz = vardata(t) 
              if(vname.eq."SMCMAX") &
                   NoahMP36_struc(n)%noahmp36(t)%smcmax = vardata(t)
              if(vname.eq."SMCREF") &
                   NoahMP36_struc(n)%noahmp36(t)%smcref = vardata(t) 
              if(vname.eq."SMCWLT") &
                   NoahMP36_struc(n)%noahmp36(t)%smcwlt = vardata(t) 
              if(vname.eq."CZIL") &
                   NoahMP36_struc(n)%noahmp36(t)%czil = vardata(t) 
              if(vname.eq."FRZK") &
                   NoahMP36_struc(n)%noahmp36(t)%frzk = vardata(t) 
              if(vname.eq."REFDK") &
                   NoahMP36_struc(n)%noahmp36(t)%refdk = vardata(t) 
              if(vname.eq."REFKDT") &
                   NoahMP36_struc(n)%noahmp36(t)%refkdt = vardata(t) 
              if(vname.eq."SLOPE") &
                   NoahMP36_struc(n)%noahmp36(t)%slope = vardata(t) 
              if(vname.eq."CH2OP") &
                   NoahMP36_struc(n)%noahmp36(t)%CH2OP = vardata(t)
              if(vname.eq."DLEAF") &
                   NoahMP36_struc(n)%noahmp36(t)%DLEAF = vardata(t)
              if(vname.eq."Z0MVT") &
                   NoahMP36_struc(n)%noahmp36(t)%Z0MVT = vardata(t)
              if(vname.eq."HVT") &
                   NoahMP36_struc(n)%noahmp36(t)%HVT = vardata(t)
              if(vname.eq."HVB") &
                   NoahMP36_struc(n)%noahmp36(t)%HVB = vardata(t)
              if(vname.eq."RC") &
                   NoahMP36_struc(n)%noahmp36(t)%RC = vardata(t)
              if(vname.eq."RHOL1") &
                   NoahMP36_struc(n)%noahmp36(t)%RHOL1 = vardata(t)
              if(vname.eq."RHOL2") &
                   NoahMP36_struc(n)%noahmp36(t)%RHOL2 = vardata(t)
              if(vname.eq."RHOS1") &
                   NoahMP36_struc(n)%noahmp36(t)%RHOS1 = vardata(t)
              if(vname.eq."RHOS2") &
                   NoahMP36_struc(n)%noahmp36(t)%RHOS2 = vardata(t)
              if(vname.eq."TAUL1") &
                   NoahMP36_struc(n)%noahmp36(t)%TAUL1 = vardata(t)
              if(vname.eq."TAUL2") &
                   NoahMP36_struc(n)%noahmp36(t)%TAUL2 = vardata(t)
              if(vname.eq."TAUS1") &
                   NoahMP36_struc(n)%noahmp36(t)%TAUS1 = vardata(t)
              if(vname.eq."TAUS2") &
                   NoahMP36_struc(n)%noahmp36(t)%TAUS2 = vardata(t)
              if(vname.eq."XL") &
                   NoahMP36_struc(n)%noahmp36(t)%XL = vardata(t)
              if(vname.eq."CWPVT") &
                   NoahMP36_struc(n)%noahmp36(t)%CWPVT = vardata(t)
              if(vname.eq."C3PSN") &
                   NoahMP36_struc(n)%noahmp36(t)%C3PSN = vardata(t)
              if(vname.eq."KC25") &
                   NoahMP36_struc(n)%noahmp36(t)%KC25 = vardata(t)
              if(vname.eq."AKC") &
                   NoahMP36_struc(n)%noahmp36(t)%AKC = vardata(t)
              if(vname.eq."KO25") &
                   NoahMP36_struc(n)%noahmp36(t)%KO25 = vardata(t)
              if(vname.eq."AKO") &
                   NoahMP36_struc(n)%noahmp36(t)%AKO = vardata(t)
              if(vname.eq."AVCMX") &
                   NoahMP36_struc(n)%noahmp36(t)%AVCMX = vardata(t)
              if(vname.eq."AQE") &
                   NoahMP36_struc(n)%noahmp36(t)%AQE = vardata(t)
              if(vname.eq."LTOVRC") &
                   NoahMP36_struc(n)%noahmp36(t)%LTOVRC = vardata(t)
              if(vname.eq."DILEFC") &
                   NoahMP36_struc(n)%noahmp36(t)%DILEFC = vardata(t)
              if(vname.eq."DILEFW") &
                   NoahMP36_struc(n)%noahmp36(t)%DILEFW = vardata(t)
              if(vname.eq."RMF25") &
                   NoahMP36_struc(n)%noahmp36(t)%RMF25 = vardata(t)
              if(vname.eq."SLA") &
                   NoahMP36_struc(n)%noahmp36(t)%SLA = vardata(t)
              if(vname.eq."FRAGR") &
                   NoahMP36_struc(n)%noahmp36(t)%FRAGR = vardata(t)
              if(vname.eq."TMIN") &
                   NoahMP36_struc(n)%noahmp36(t)%TMIN = vardata(t)
              if(vname.eq."VCMX25") &
                   NoahMP36_struc(n)%noahmp36(t)%VCMX25 = vardata(t)
              if(vname.eq."TDLEF") &
                   NoahMP36_struc(n)%noahmp36(t)%TDLEF = vardata(t)
              if(vname.eq."BP") &
                   NoahMP36_struc(n)%noahmp36(t)%BP = vardata(t)
              if(vname.eq."MP") &
                   NoahMP36_struc(n)%noahmp36(t)%MP = vardata(t)
              if(vname.eq."QE25") &
                   NoahMP36_struc(n)%noahmp36(t)%QE25 = vardata(t)
              if(vname.eq."RMS25") &
                   NoahMP36_struc(n)%noahmp36(t)%RMS25 = vardata(t)
              if(vname.eq."RMR25") &
                   NoahMP36_struc(n)%noahmp36(t)%RMR25 = vardata(t)
              if(vname.eq."ARM") &
                   NoahMP36_struc(n)%noahmp36(t)%ARM = vardata(t)
              if(vname.eq."FOLNMX") &
                   NoahMP36_struc(n)%noahmp36(t)%FOLNMX = vardata(t)
              if(vname.eq."WDPOOL") &
                   NoahMP36_struc(n)%noahmp36(t)%WDPOOL = vardata(t)
              if(vname.eq."WRRAT") &
                   NoahMP36_struc(n)%noahmp36(t)%WRRAT = vardata(t)
              if(vname.eq."MRP") &
                   NoahMP36_struc(n)%noahmp36(t)%MRP = vardata(t)
           endif
        enddo
     endif
  enddo
end subroutine NoahMP36_setVars

subroutine NoahMP36_setModFlag(n,DEC_State,Feas_State,mod_flag_NoahMP36)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_NoahMP36(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
!  write(LIS_logunit,*) 'NoahMP36_setModFlag 1 modflag:', modflag

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(mod_flag_NoahMP36(t).eq.1) then 
        modflag(t)=1
     endif
  enddo
!  write(LIS_logunit,*) 'NoahMP36_setModFlag 2 modflag:', modflag
end subroutine NoahMP36_setModFlag
