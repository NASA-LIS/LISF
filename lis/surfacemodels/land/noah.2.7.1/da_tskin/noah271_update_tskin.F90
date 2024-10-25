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
! !ROUTINE: noah271_update_tskin
! \label{noah271_update_tskin}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!
! !INTERFACE:
subroutine noah271_update_tskin(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use noah271_lsmMod
  use LIS_logMod, only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)   :: soiltField, soiltIncrField
  integer            :: status, t
  real, allocatable      :: soilt(:), soiltIncr(:)
  
!  call ESMF_StateGet(LSM_State, "Skin Temperature", tskinField, rc=status)
!  call LIS_verify(status)
  
!  call ESMF_FieldGet(tskinField,localDE=0,farrayPtr= tskin, rc=status)
!  call LIS_verify(status)

  call ESMF_StateGet(LSM_State, "Soil Temperature", soiltField, rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(soiltField,localDE=0,farrayPtr= soilt, rc=status)
  call LIS_verify(status)

!  call ESMF_StateGet(LSM_Incr_State, "Skin Temperature", tskinIncrField, rc=status)
!  call LIS_verify(status)
  
!  call ESMF_FieldGet(tskinIncrField,localDE=0,farrayPtr= tskinIncr, rc=status)
!  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State, "Soil Temperature", soiltIncrField, rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(soiltIncrField,localDE=0,farrayPtr= soiltIncr, rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
!     tskin(t) = tskin(t) + tskinIncr(t)
     soilt(t) = soilt(t) + soiltIncr(t)
!     gid = LIS_domain(n)%tile(t)%index
!     if(gid.eq.1) write(LIS_logunit,*)'upd ',t,tskin(t), tskinIncr(t)
  enddo

end subroutine noah271_update_tskin

