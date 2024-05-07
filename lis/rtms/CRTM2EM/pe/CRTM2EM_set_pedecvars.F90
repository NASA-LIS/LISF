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
! !ROUTINE: crtm2em_set_pedecvars
!  \label{crtm2em_set_pedecvars}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!


! !INTERFACE:
subroutine crtm2em_set_pedecvars(DEC_State, Feas_State)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
#if (defined RTMS) 
  use CRTM2_EMMod, only : crtm_struc
  use crtm2em_peMod,  only : crtm2em_pe_struc
#endif

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: DEC_State
  type(ESMF_State)       :: Feas_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to CRTM2EM's model variables. 
! 
!EOP
  integer                :: n
  real, allocatable          :: vdata(:)
  character*100          :: vname
  integer, allocatable       :: mod_flag_crtm2em(:)
  integer                :: i,t
  integer                :: status

#if (defined RTMS) 
  n = 1

  allocate(mod_flag_crtm2em(LIS_rc%npatch(n,LIS_rc%lsm_index)))

  mod_flag_crtm2em = 0
    
  !set modflag based on bounds
  allocate(vdata(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  do i=1,crtm2em_pe_struc(n)%nparams
     if(crtm2em_pe_struc(n)%param_select(i).eq.1) then
        vname=trim(crtm2em_pe_struc(n)%param_name(i))
        call getvardata_crtm2em(n,DEC_State,vname, vdata, status)
        call LIS_verify(status)
        call checkBounds_crtm2em(n,DEC_State,vname, vdata, mod_flag_crtm2em)
     endif
  enddo
  deallocate(vdata) 

  !update modflags based on constraints
  call checkConstraints_crtm2em(n,DEC_State, mod_flag_crtm2em)

  !set variables given modflag; if flag set will leave values alone
  call setVars_crtm2em(n,DEC_State,mod_flag_crtm2em)

  !send mod flag to ESMF state (feasibility flag)
  call setModFlag_crtm2em(n,DEC_State,Feas_State,mod_flag_crtm2em)
#endif
end subroutine crtm2em_set_pedecvars

!BOP
! 
! !ROUTINE: randArray
! \label{randArray}
!
! !INTERFACE: 
subroutine getvardata_crtm2em(n,DEC_State,vname, vdata, statusStateGet)
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
!  This routine assigns the decision space to CRTM2EM's model variables. 
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
  
end subroutine getvardata_crtm2em

subroutine checkBounds_crtm2em(n,DEC_State,vname, vardata, mod_flag_crtm2em)
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
  integer       :: mod_flag_crtm2em(LIS_rc%npatch(n,LIS_rc%lsm_index))
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to CRTM2EM's model variables. 
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
        mod_flag_crtm2em(t) = 1
     endif
     if(vardata(t).gt.vardata_max) then 
        mod_flag_crtm2em(t) = 1
     endif
  enddo
end subroutine checkBounds_crtm2em

subroutine checkConstraints_crtm2em(n,DEC_State,mod_flag_crtm2em)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
#if (defined RTMS) 
  use CRTM2_EMMod, only : crtm_struc
#endif
  implicit none
! !ARGUMENTS: 
  integer                :: n
  type(ESMF_State)       :: DEC_State
  integer       :: mod_flag_crtm2em(LIS_rc%npatch(n,LIS_rc%lsm_index))
#if (defined RTMS) 
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to CRTM2EM's model variables. 
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
!  call getvardata_crtm2em(n,DEC_State,vname,vardata1, status1)
!  vname='SMCDRY'
!  call getvardata_crtm2em(n,DEC_State,vname,vardata2, status2)
!  if(status1.ne.0) vardata1=crtm_struc(n)%noah(:)%smcmax
!  if(status2.ne.0) vardata2=crtm_struc(n)%noah(:)%smcdry
!  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     if(vardata1(t).le.vardata2(t)) then
!        mod_flag_crtm2em(t) = 1
!     endif
!  enddo
#endif
end subroutine checkConstraints_crtm2em

subroutine setVars_crtm2em(n,DEC_State,mod_flag_crtm2em)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify
#if (defined RTMS) 
  use CRTM2_EMMod, only : crtm_struc
  use crtm2em_peMod,  only : crtm2em_pe_struc
#endif
  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_crtm2em(LIS_rc%npatch(n,LIS_rc%lsm_index))
  type(ESMF_State)       :: DEC_State
#if (defined RTMS) 
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to CRTM2EM's model variables. 
!  Only does so if the proposed parameter set is feasible (meets bounds and constraints)
! 
!EOP
  real          :: vardata(LIS_rc%npatch(n,LIS_rc%lsm_index))
  character*100          :: vname
  integer                :: i,t, status

  do i=1,crtm2em_pe_struc(n)%nparams
     if(crtm2em_pe_struc(n)%param_select(i).eq.1) then 
        vname=trim(crtm2em_pe_struc(n)%param_name(i))
        call getvardata_crtm2em(n,DEC_State,vname,vardata, status)
        call LIS_verify(status)
        do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           if(mod_flag_crtm2em(t).eq.0) then 
              if(vname.eq."SIGMA") &
                   crtm_struc(n)%SFC(t)%sigma = vardata(t)
              if(vname.eq."LEAF_THICK") &
                   crtm_struc(n)%SFC(t)%leaf_thick = vardata(t)
              if(vname.eq."BGF_FIXED") &
                   crtm_struc(n)%bgf_fixed(t) = vardata(t)
              if(vname.eq."K_LAI2VGF") &
                   crtm_struc(n)%k_lai2vgf(t) = vardata(t)
              if(vname.eq."WATER_CONTENT_PER_LAI") &
                   crtm_struc(n)%SFC(t)%water_content_per_lai = vardata(t)
              if(vname.eq."SSALB_FACTOR") &
                   crtm_struc(n)%SFC(t)%ssalb_factor = vardata(t)
           endif
        enddo
     endif
  enddo
#endif
end subroutine setVars_crtm2em

subroutine setModFlag_crtm2em(n,DEC_State,Feas_State,mod_flag_crtm2em)
! !USES:
  use ESMF
  use LIS_coreMod,   only : LIS_rc, LIS_surface
  use LIS_logMod,    only : LIS_logunit,LIS_verify

  implicit none
! !ARGUMENTS: 
  integer                :: n
  integer       :: mod_flag_crtm2em(LIS_rc%npatch(n,LIS_rc%lsm_index))
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
     if(mod_flag_crtm2em(t).eq.1) then 
        modflag(t)=1
     endif
  enddo
end subroutine setModFlag_crtm2em
