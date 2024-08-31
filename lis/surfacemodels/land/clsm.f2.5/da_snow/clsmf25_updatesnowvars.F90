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
! !ROUTINE: clsmf25_updatesnowvars
! \label{clsmf25_updatesnowvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!
! !INTERFACE:
subroutine clsmf25_updatesnowvars(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use clsmf25_lsmMod
  use LIS_logMod,   only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Returns the snow related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \item[LSM\_Incr\_State] ESMF State container for LSM state increments \newline
!  \end{description}
!
!EOP

  type(ESMF_Field)       :: swe1Field,swe2Field,swe3Field
  type(ESMF_Field)       :: swe1IncrField,swe2IncrField,swe3IncrField

  integer                :: t
  integer                :: status
  real, pointer          :: swe1(:),swe2(:),swe3(:)
  real, pointer          :: swe1incr(:),swe2incr(:),swe3incr(:)
 
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 1",swe1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 2",swe2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 3",swe3Field,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"Water Equivalent Snow 1",swe1incrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Water Equivalent Snow 2",swe2incrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Water Equivalent Snow 3",swe3incrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swe1Field,localDE=0,farrayPtr=swe1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(swe2Field,localDE=0,farrayPtr=swe2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(swe3Field,localDE=0,farrayPtr=swe3,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swe1incrField,localDE=0,farrayPtr=swe1incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(swe2incrField,localDE=0,farrayPtr=swe2incr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(swe3incrField,localDE=0,farrayPtr=swe3incr,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if(swe1(t)+ swe1Incr(t).ge.0) then 
        swe1(t) = swe1(t) + swe1Incr(t)
     endif
     if(swe2(t) + swe2Incr(t).ge.0) then 
        swe2(t) = swe2(t) + swe2Incr(t)
     endif
     if(swe3(t) + swe3Incr(t).ge.0) then 
        swe3(t) = swe3(t) + swe3Incr(t)
     endif
     if(swe1(t).lt.0) then 
        print*, 'problem in upd ',t, swe1(t),swe2(t), swe3(t)
     endif
  enddo
  
end subroutine clsmf25_updatesnowvars

