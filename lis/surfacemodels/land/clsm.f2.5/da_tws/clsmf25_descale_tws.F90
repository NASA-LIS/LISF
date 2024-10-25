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
!
! !ROUTINE: clsmf25_descale_tws
! \label{clsmf25_descale_tws}
!
! !REVISION HISTORY:
! 12 Feb 2006; Sujay Kumar, Initial Code
!
! !INTERFACE:    
subroutine clsmf25_descale_tws(n,LSM_State, LSM_Incr_State)
! !USES:
  use ESMF
  use LIS_coreMod, only  : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use clsmf25_lsmMod

  implicit none
! !ARGUMENTS:   
  integer, intent(IN)       :: n 
  type(ESMF_State)          :: LSM_State
  type(ESMF_State)          :: LSM_Incr_State
!
! !DESCRIPTION:
!  Routine to retrieve soil moisture prognostic variables in catchment.
!
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP  
  type(ESMF_Field)          :: catdefField,surfField, rzoneField
  real, pointer             :: catdef(:),surf(:),rzone(:)
  type(ESMF_Field)          :: catdefIncField,surfIncField,rzoneIncField
  real, pointer             :: catdefInc(:),surfInc(:),rzoneInc(:)
  integer                   :: t
  integer                   :: status

  call ESMF_StateGet(LSM_State, "Catchment Deficit",catdefField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State, "Root Zone Excess",catdefField,rc=status)
  call LIS_verify(status)  
  call ESMF_StateGet(LSM_State, "Surface Excess",surfField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(catdefField, localDE=0, farrayPtr=catdef,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(catdefField, localDE=0, farrayPtr=rzone,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(surfField,localDE=0, farrayPtr=surf,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_Incr_State, "Catchment Deficit",catdefIncField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_Incr_State, "Root Zone Excess",rzoneIncField,rc=status)
  call LIS_verify(status)  
  call ESMF_StateGet(LSM_Incr_State, "Surface Excess",surfIncField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(catdefIncField, localDE=0, farrayPtr=catdefInc,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(rzoneIncField, localDE=0, farrayPtr=rzoneInc,rc=status)
  call LIS_verify(status)  
  call ESMF_FieldGet(surfIncField,localDE=0, farrayPtr=surfInc,rc=status)
  call LIS_verify(status)


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
    catdef(t) = catdef(t)*500.0  
    surf(t)  =  surf(t)*5.0 
    
    catdefInc(t) = catdefInc(t)*500.0
    surfInc(t) = surfInc(t)*5.0
  enddo

end subroutine clsmf25_descale_tws
