!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noah39_qc_ldtsiobs
! \label{noah39_qc_ldtsiobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 09Apr2019: Eric Kemp: Updated to Noah 3.9 and LDT-SI
!
! !INTERFACE:
subroutine noah39_qc_ldtsiobs(n,OBS_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only : LIS_verify
  use noah39_lsmMod

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

  real, pointer            :: ldtsiobs(:)
  integer                  :: status

  
  call ESMF_StateGet(OBS_State,"Observation01",obs_snow_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in noah39_qc_ldtsiobs")
  call ESMF_FieldGet(obs_snow_field,localDE=0,farrayPtr=ldtsiobs,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet failed in noah39_qc_ldtsiobs")
  
  !do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

!     gid  = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index

!     if(noah39_struc(n)%noah(t)%shdfac.gt.0.7) then 
!        ldtsiobs(gid) = LIS_rc%udef        
!     endif
  !enddo

end subroutine noah39_qc_ldtsiobs

