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
! !ROUTINE: noah33_set_pedecvars
!  \label{noah33_set_pedecvars}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!


! !INTERFACE:
subroutine noah33_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use noah33_lsmMod, only : noah33_struc
  use noah33_peMod,  only : noah33_pe_struc

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Noah's model variables. 
! 
!EOP
  integer                :: n
  real, pointer          :: vdata(:)
  character*100          :: vname
  integer, pointer       :: mod_flag_noah33(:)
  integer                :: i,t
  integer                :: status

  n = 1

  allocate(mod_flag_noah33(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_noah33 = 0
    
  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,noah33_pe_struc(n)%nparams
     if(noah33_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(noah33_pe_struc(n)%param_name(i))
        call getvardata(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call checkBounds(n,DEC_State,vname, vdata, mod_flag_noah33)
     endif
  enddo
  deallocate(vdata) 

  !update modflags based on constraints
  call checkConstraints(n,DEC_State, mod_flag_noah33)

  !set variables given modflag; if flag set will leave values alone
  call setVars(n,DEC_State,mod_flag_noah33)

  !send mod flag to ESMF state (feasibility flag)
  call setModFlag(n,DEC_State,Feas_State,mod_flag_noah33)
end subroutine noah33_set_pedecvars

!BOP
! 
! !ROUTINE: randArray
! \label{randArray}
!
! !INTERFACE: 
subroutine getvardata(n,DEC_State,vname, vdata, statusStateGet)
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
!  This routine assigns the decision space to Noah's model variables. 
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
  
end subroutine getvardata

subroutine checkBounds(n,DEC_State,vname, vardata, mod_flag_noah33)
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
  integer       :: mod_flag_noah33(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Noah's model variables. 
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
        mod_flag_noah33(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then 
        mod_flag_noah33(t) = 1
     endif
  enddo
end subroutine checkBounds

subroutine checkConstraints(n,DEC_State,mod_flag_noah33)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use noah33_lsmMod, only : noah33_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  integer       :: mod_flag_noah33(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Noah's model variables. 
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
  call getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCDRY'
  call getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=noah33_struc(n)%noah(:)%smcmax
  if(status2.ne.0) vardata2=noah33_struc(n)%noah(:)%smcdry
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_noah33(t) = 1
     endif
  enddo

  !SMCREF > SMCWLT
  vname='SMCREF'
  call getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCWLT'
  call getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=noah33_struc(n)%noah(:)%smcref
  if(status2.ne.0) vardata2=noah33_struc(n)%noah(:)%smcwlt
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_noah33(t) = 1
     endif
  enddo

  !SMCMAX > SMCREF
  vname='SMCMAX'
  call getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMCREF'
  call getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=noah33_struc(n)%noah(:)%smcmax
  if(status2.ne.0) vardata2=noah33_struc(n)%noah(:)%smcref
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        mod_flag_noah33(t) = 1
     endif
  enddo

  !SMCMAX > SMC1
  vname='SMCMAX'
  call getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMC1'
  call getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=noah33_struc(n)%noah(:)%smcmax
  if(status2.ne.0) then
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        vardata2(t)=noah33_struc(n)%noah(t)%smc(1)
     enddo
  endif
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        if(status2.ne.0) then !smc1 not optimized, just fix smc1 as spinup relied on
           noah33_struc(n)%noah(t)%smc(1)=vardata1(t)
        else  !both optimized, need to declare infeasible           
           mod_flag_noah33(t) = 1
        endif
     endif
  enddo

  !SMCMAX > SMC2
  vname='SMCMAX'
  call getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMC2'
  call getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=noah33_struc(n)%noah(:)%smcmax
  if(status2.ne.0) then
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        vardata2(t)=noah33_struc(n)%noah(t)%smc(2)
     enddo
  endif
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        if(status2.ne.0) then !smc2 not optimized, just fix smc2 as spinup relied on
           noah33_struc(n)%noah(t)%smc(2)=vardata1(t)
        else  !both optimized, need to declare infeasible           
           mod_flag_noah33(t) = 1
        endif
     endif
  enddo

  !SMCMAX > SMC3
  vname='SMCMAX'
  call getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMC3'
  call getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=noah33_struc(n)%noah(:)%smcmax
  if(status2.ne.0) then
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        vardata2(t)=noah33_struc(n)%noah(t)%smc(3)
     enddo
  endif
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        if(status2.ne.0) then !smc3 not optimized, just fix smc3 as spinup relied on
           noah33_struc(n)%noah(t)%smc(3)=vardata1(t)
        else  !both optimized, need to declare infeasible           
           mod_flag_noah33(t) = 1
        endif
     endif
  enddo

  !SMCMAX > SMC4
  vname='SMCMAX'
  call getvardata(n,DEC_State,vname,vardata1, status1)
  vname='SMC4'
  call getvardata(n,DEC_State,vname,vardata2, status2)
  if(status1.ne.0) vardata1=noah33_struc(n)%noah(:)%smcmax
  if(status2.ne.0) then
     do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        vardata2(t)=noah33_struc(n)%noah(t)%smc(4)
     enddo
  endif
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(vardata1(t).le.vardata2(t)) then
        if(status2.ne.0) then !smc4 not optimized, just fix smc4 as spinup relied on
           noah33_struc(n)%noah(t)%smc(4)=vardata1(t)
        else  !both optimized, need to declare infeasible           
           mod_flag_noah33(t) = 1
        endif
     endif
  enddo
  deallocate(vardata1)
  deallocate(vardata2)

end subroutine checkConstraints

subroutine setVars(n,DEC_State,mod_flag_noah33)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
  use noah33_lsmMod, only : noah33_struc
  use noah33_peMod,  only : noah33_pe_struc

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_noah33(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)       :: DEC_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to Noah's model variables. 
!  Only does so if the proposed parameter set is feasible (meets bounds and constraints)
! 
!EOP
  real          :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  character*100          :: vname
  integer                :: i,t, status

  do i=1,noah33_pe_struc(n)%nparams
     if(noah33_pe_struc(n)%param_select(i).eq.1) then 
        vname=trim(noah33_pe_struc(n)%param_name(i))
        call getvardata(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if(mod_flag_noah33(t).eq.0) then 
              if(vname.eq."SMCMAX") &
                   noah33_struc(n)%noah(t)%smcmax = vardata(t)
              if(vname.eq."PSISAT") & 
                   noah33_struc(n)%noah(t)%psisat = vardata(t) 
              if(vname.eq."DKSAT") & 
                   noah33_struc(n)%noah(t)%dksat = vardata(t) 
              if(vname.eq."DWSAT") & 
                   noah33_struc(n)%noah(t)%dwsat = vardata(t) 
              if(vname.eq."BEXP") & 
                   noah33_struc(n)%noah(t)%bexp = vardata(t) 
              if(vname.eq."QUARTZ") & 
                   noah33_struc(n)%noah(t)%quartz = vardata(t) 
              if(vname.eq."RSMIN") & 
                   noah33_struc(n)%noah(t)%rsmin = vardata(t) 
              if(vname.eq."RGL") & 
                   noah33_struc(n)%noah(t)%rgl = vardata(t) 
              if(vname.eq."HS") & 
                   noah33_struc(n)%noah(t)%hs = vardata(t) 
              if(vname.eq."Z0") &
                   noah33_struc(n)%noah(t)%z0 = vardata(t) 
              if(vname.eq."LAI") &
                   noah33_struc(n)%noah(t)%lai = vardata(t) 
              if(vname.eq."CFACTR") &
                   noah33_struc(n)%noah(t)%cfactr = vardata(t) 
              if(vname.eq."CMCMAX") &
                   noah33_struc(n)%noah(t)%cmcmax = vardata(t) 
              if(vname.eq."SBETA") &
                   noah33_struc(n)%noah(t)%sbeta = vardata(t) 
              if(vname.eq."RSMAX") &
                   noah33_struc(n)%noah(t)%rsmax = vardata(t) 
              if(vname.eq."TOPT") &
                   noah33_struc(n)%noah(t)%topt = vardata(t) 
              if(vname.eq."REFDK") &
                   noah33_struc(n)%noah(t)%refdk = vardata(t) 
              if(vname.eq."FXEXP") &
                   noah33_struc(n)%noah(t)%fxexp = vardata(t) 
              if(vname.eq."REFKDT") &
                   noah33_struc(n)%noah(t)%refkdt = vardata(t) 
              if(vname.eq."CZIL") &
                   noah33_struc(n)%noah(t)%czil = vardata(t) 
              if(vname.eq."CSOIL") &
                   noah33_struc(n)%noah(t)%csoil = vardata(t) 
              if(vname.eq."FRZK") &
                   noah33_struc(n)%noah(t)%frzk = vardata(t) 
              if(vname.eq."SNUP") &
                   noah33_struc(n)%noah(t)%snup = vardata(t) 
              if(vname.eq."SMCREF") &
                   noah33_struc(n)%noah(t)%smcref = vardata(t) 
              if(vname.eq."SMCDRY") &
                   noah33_struc(n)%noah(t)%smcdry = vardata(t) 
              if(vname.eq."SMCWLT") &
                   noah33_struc(n)%noah(t)%smcwlt = vardata(t) 
              if(vname.eq."F1") &
                   noah33_struc(n)%noah(t)%f1 = vardata(t) 
              if(vname.eq."SLOPE") &
                   noah33_struc(n)%noah(t)%slope = vardata(t) 
              if(vname.eq."EMISS") &
                   noah33_struc(n)%noah(t)%emiss = vardata(t)
              if(vname.eq."SIGMA_FLX") &
                   noah33_struc(n)%noah(t)%sigma_flx = vardata(t) 
              if(vname.eq."SIGMA_SM") &
                   noah33_struc(n)%noah(t)%sigma_sm = vardata(t) 
              if(vname.eq."SMC1") &
                   noah33_struc(n)%noah(t)%smc(1) = vardata(t) 
              if(vname.eq."SMC2") &
                   noah33_struc(n)%noah(t)%smc(2) = vardata(t) 
              if(vname.eq."SMC3") &
                   noah33_struc(n)%noah(t)%smc(3) = vardata(t) 
              if(vname.eq."SMC4") &
                   noah33_struc(n)%noah(t)%smc(4) = vardata(t) 
              if(vname.eq."STC1") &
                   noah33_struc(n)%noah(t)%stc(1) = vardata(t) 
              if(vname.eq."STC2") &
                   noah33_struc(n)%noah(t)%stc(2) = vardata(t) 
              if(vname.eq."STC3") &
                   noah33_struc(n)%noah(t)%stc(3) = vardata(t) 
              if(vname.eq."STC4") &
                   noah33_struc(n)%noah(t)%stc(4) = vardata(t) 
           endif
        enddo
     endif
  enddo
end subroutine setVars

subroutine setModFlag(n,DEC_State,Feas_State,mod_flag_noah33)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_noah33(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
     if(mod_flag_noah33(t).eq.1) then 
        modflag(t)=1
     endif
  enddo
end subroutine setModFlag
