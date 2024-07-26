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
! !ROUTINE: clsmf25_getsnowvars
! \label{clsmf25_getsnowvars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!
! !INTERFACE:
subroutine clsmf25_getsnowvars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use clsmf25_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
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
!  \end{description}
!
!EOP
  type(ESMF_Field)       :: swe1Field,swe2Field,swe3Field

  integer                :: t
  integer                :: status
  real, pointer          :: swe1(:),swe2(:),swe3(:)
 
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 1",swe1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 2",swe2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 3",swe3Field,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swe1Field,localDE=0,farrayPtr=swe1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(swe2Field,localDE=0,farrayPtr=swe2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(swe3Field,localDE=0,farrayPtr=swe3,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     swe1(t) = clsmf25_struc(n)%cat_progn(t)%wesn(1)
     swe2(t) = clsmf25_struc(n)%cat_progn(t)%wesn(2)
     swe3(t) = clsmf25_struc(n)%cat_progn(t)%wesn(3)
  enddo
  
end subroutine clsmf25_getsnowvars

