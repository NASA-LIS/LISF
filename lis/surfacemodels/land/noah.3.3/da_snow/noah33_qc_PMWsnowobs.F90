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
! !ROUTINE: noah33_qc_PMWsnowobs
! \label{noah33_qc_PMWsnowobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 13Mar2014: Yuqiong Liu, modified for qc PMW SWE and snow depth obs
! !INTERFACE:
subroutine noah33_qc_PMWsnowobs(n,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use LIS_constantsMod, only : LIS_CONST_TKFRZ
  use noah33_lsmMod
  use PMW_snow_Mod, only : PMW_snow_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  type(ESMF_State)         :: OBS_State
!
! !DESCRIPTION:
!
!  This subroutine performs any model-based QC of the observation 
!  prior to data assimilation. Here the snow observations
!  are flagged when LSM indicates that (1) rain is falling (2)
!  ground is fully or partially covered with snow. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[OBS\_State] ESMF state container for observations \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)         :: obs_snow_field

  real, pointer            :: snowobs(:)
  integer                  :: t
  integer                  :: gid
  integer                  :: status

  
  call ESMF_StateGet(OBS_State,"Observation01",obs_snow_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in noah33_qc_PMWsnowobs")
  call ESMF_FieldGet(obs_snow_field,localDE=0,farrayPtr=snowobs,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet failed in noah33_qc_PMWsnowobs")
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     gid  = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index

     if(snowobs(gid).ne.LIS_rc%udef) then 
        !GVF mask
        if(PMW_snow_struc(n)%gvf_mask.eq.1 .and. noah33_struc(n)%noah(t)%shdfac.gt.0.7) then 
           snowobs(gid) = LIS_rc%udef
        !forest mask        
        elseif(PMW_snow_struc(n)%lctype_mask.eq.1 .and. noah33_struc(n)%noah(t)%vegt.le.4) then
           snowobs(gid) = LIS_rc%udef
        !LSM temperature mask
        elseif(PMW_snow_struc(n)%temp_mask.eq.1 .and. noah33_struc(n)%noah(t)%t1.ge.278.15) then 
           snowobs(gid) = LIS_rc%udef
        !LSM temperature mask
        elseif(PMW_snow_struc(n)%temp_mask.eq.1 .and. noah33_struc(n)%noah(t)%stc(1).ge.278.15) then 
           snowobs(gid) = LIS_rc%udef
        endif
     endif
  enddo

end subroutine noah33_qc_PMWsnowobs

