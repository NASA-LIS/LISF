!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: cmem3_set_pedecvars
!  \label{cmem3_set_pedecvars}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!


! !INTERFACE:
subroutine cmem3_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
#if (defined RTMS) 
  use CMEM3_Mod, only : cmem3_struc
  use cmem3_peMod,  only : cmem3_pe_struc
#endif

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to CMEM3's model variables. 
! 
!EOP
  integer                :: n
  real, allocatable          :: vdata(:)
  character*100          :: vname
  integer, allocatable       :: mod_flag_cmem3(:)
  integer                :: i,t
  integer                :: status

#if (defined RTMS) 
  n = 1

  allocate(mod_flag_cmem3(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_cmem3 = 0
    
  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,cmem3_pe_struc(n)%nparams
     if(cmem3_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(cmem3_pe_struc(n)%param_name(i))
        call getvardata_cmem3(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call checkBounds_cmem3(n,DEC_State,vname, vdata, mod_flag_cmem3)
     endif
  enddo
  deallocate(vdata) 

  !update modflags based on constraints
  call checkConstraints_cmem3(n,DEC_State, mod_flag_cmem3)

  !set variables given modflag; if flag set will leave values alone
  call setVars_cmem3(n,DEC_State,mod_flag_cmem3)

  !send mod flag to ESMF state (feasibility flag)
  call setModFlag_cmem3(n,DEC_State,Feas_State,mod_flag_cmem3)
#endif
end subroutine cmem3_set_pedecvars

!BOP
! 
! !ROUTINE: randArray
! \label{randArray}
!
! !INTERFACE: 
subroutine getvardata_cmem3(n,DEC_State,vname, vdata, statusStateGet)
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
!  This routine assigns the decision space to CMEM3's model variables. 
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
  
end subroutine getvardata_cmem3

subroutine checkBounds_cmem3(n,DEC_State,vname, vardata, mod_flag_cmem3)
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
  integer       :: mod_flag_cmem3(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to CMEM3's model variables. 
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
        mod_flag_cmem3(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then 
        mod_flag_cmem3(t) = 1
     endif
  enddo
end subroutine checkBounds_cmem3

subroutine checkConstraints_cmem3(n,DEC_State,mod_flag_cmem3)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
#if (defined RTMS) 
  use CMEM3_Mod, only : cmem3_struc
#endif
  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  integer       :: mod_flag_cmem3(LIS_rc%npatch(n,LIS_rc%lsm_index))
#if (defined RTMS) 
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to CMEM3's model variables. 
! 
!EOP
  type(ESMF_Field)       :: varField
  real                   :: vardata_min, vardata_max
  character*100          :: vname
  integer                :: t
  real          :: vardata1(LIS_rc%npatch(n,LIS_rc%lsm_index))
  real          :: vardata2(LIS_rc%npatch(n,LIS_rc%lsm_index))
  integer       :: status1, status2

!  !SMCMAX > SMCDRY
!  vname='SMCMAX'
!  call getvardata_cmem3(n,DEC_State,vname,vardata1, status1)
!  vname='SMCDRY'
!  call getvardata_cmem3(n,DEC_State,vname,vardata2, status2)
!  if(status1.ne.0) vardata1=cmem3_struc(n)%noah(:)%smcmax
!  if(status2.ne.0) vardata2=cmem3_struc(n)%noah(:)%smcdry
!  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     if(vardata1(t).le.vardata2(t)) then
!        mod_flag_cmem3(t) = 1
!     endif
!  enddo
#endif
end subroutine checkConstraints_cmem3

subroutine setVars_cmem3(n,DEC_State,mod_flag_cmem3)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
#if (defined RTMS) 
  use CMEM3_Mod, only : cmem3_struc
  use cmem3_peMod,  only : cmem3_pe_struc
#endif
  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_cmem3(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)       :: DEC_State
#if (defined RTMS) 
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to CMEM3's model variables. 
!  Only does so if the proposed parameter set is feasible (meets bounds and constraints)
! 
!EOP
  real          :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  character*100          :: vname
  integer                :: i,t, status

  do i=1,cmem3_pe_struc(n)%nparams
     if(cmem3_pe_struc(n)%param_select(i).eq.1) then 
        vname=trim(cmem3_pe_struc(n)%param_name(i))
        call getvardata_cmem3(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if(mod_flag_cmem3(t).eq.0) then 
              if(vname.eq."SRMAX") &
                   cmem3_struc(n)%srmax(t) = vardata(t)
              if(vname.eq."SRMAX2SRMIN") &
                   cmem3_struc(n)%srmax2srmin(t) = vardata(t)
              if(vname.eq."D_LEAF") &
                   cmem3_struc(n)%d_leaf(t) = vardata(t)
              if(vname.eq."BGF_FIXED") &
                   cmem3_struc(n)%bgf_fixed(t) = vardata(t)
              if(vname.eq."K_LAI2VGF") &
                   cmem3_struc(n)%k_lai2vgf(t) = vardata(t)
              if(vname.eq."VWC2LAI") &
                   cmem3_struc(n)%vwc2lai(t) = vardata(t)
              if(vname.eq."M_D") &
                   cmem3_struc(n)%m_d(t) = vardata(t)
           endif
        enddo
     endif
  enddo
#endif
end subroutine setVars_cmem3

subroutine setModFlag_cmem3(n,DEC_State,Feas_State,mod_flag_cmem3)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_cmem3(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
     if(mod_flag_cmem3(t).eq.1) then 
        modflag(t)=1
     endif
  enddo
end subroutine setModFlag_cmem3
