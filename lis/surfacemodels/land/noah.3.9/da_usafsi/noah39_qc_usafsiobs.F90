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
! !ROUTINE: noah39_qc_usafsiobs
! \label{noah39_qc_usafsiobs}
!
! !REVISION HISTORY:
! 25Feb2008: Sujay Kumar: Initial Specification
! 09Apr2019: Eric Kemp: Updated to Noah 3.9 and LDT-SI
! 13Dec2019: Eric Kemp: Replaced LDTSI with USAFSI
!
! !INTERFACE:
subroutine noah39_qc_usafsiobs(n,OBS_State)
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

  real, pointer            :: usafsiobs(:)
  integer                  :: status

  
  call ESMF_StateGet(OBS_State,"Observation01",obs_snow_field,&
       rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet failed in noah39_qc_usafsiobs")
  call ESMF_FieldGet(obs_snow_field,localDE=0,farrayPtr=usafsiobs,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet failed in noah39_qc_usafsiobs")
  
  !do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

!     gid  = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index

!     if(noah39_struc(n)%noah(t)%shdfac.gt.0.7) then 
!        usafsiobs(gid) = LIS_rc%udef        
!     endif
  !enddo

end subroutine noah39_qc_usafsiobs

