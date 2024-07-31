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
! !ROUTINE: noah32_updatesnodep
! \label{noah32_updatesnodep}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 21 Jul 2011: James Geiger; Modified for Noah 3.2
! 01 May 2014: Yuqiong Liu; modified for better QC
!
! !INTERFACE:
subroutine noah32_updatesnodep(n, LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use noah32_lsmMod
  use LIS_logMod,   only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Updates the related state prognostic variable objects for
!  SNODEP data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \item[LSM\_Incr\_State] ESMF State container for LSM state increments \newline
!  \end{description}
!
!EOP
  type(ESMF_Field)       :: sweField, sweIncrField
  type(ESMF_Field)       :: snodField, snodIncrField

  integer                :: t
  integer                :: status
  real, allocatable          :: swe(:), sweincr(:)
  real, allocatable          :: snod(:), snodincr(:)
 
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State,"Snowdepth",snodIncrField,rc=status)
  call LIS_verify(status)

 
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodIncrField,localDE=0,farrayPtr=snodincr,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      swe(t) = swe(t) + sweincr(t)
      snod(t) = snod(t) + snodincr(t)
      if (swe(t).lt.0 .or. snod(t).lt.0) then
         swe(t)=0.0
         snod(t)=0.0
     endif
     if (swe(t).gt.snod(t)) then
        swe(t)=snod(t)
     endif

  enddo

!  stop
end subroutine noah32_updatesnodep

