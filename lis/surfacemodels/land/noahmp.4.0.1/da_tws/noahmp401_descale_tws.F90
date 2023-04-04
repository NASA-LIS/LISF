!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_descale_tws
! \label{noahmp401_descale_tws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li; created for Noah-MP4.0.1
!
! !INTERFACE:
subroutine noahmp401_descale_tws(n, LSM_State, LSM_Incr_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noahmp401_lsmMod
  use module_sf_noahmplsm_401
  use LIS_constantsMod,  only : LIS_CONST_RHOFW ! Natt

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
  type(ESMF_State)       :: LSM_Incr_State
!
! !DESCRIPTION:
!
!  Descales twsoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP

!Wanshu

  type(ESMF_Field)     :: gwsField
  real, pointer        :: gws(:)
  type(ESMF_Field)     :: gwsIncrField
  real, pointer        :: gwsIncr(:)
  integer              :: t
  integer              :: status
  
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: sweIncrField
  real, pointer          :: swe(:)
  real, pointer          :: sweincr(:)
  integer                :: SOILTYP           ! soil type index [-]
  real                   :: MAX_THRESHOLD , MIN_THRESHOLD


  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for SWE in noahmp401_descale_tws')

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for SWE in noahmp401_descale_tws')
  call ESMF_StateGet(LSM_Incr_State,"SWE",sweIncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweIncrField,localDE=0,farrayPtr=sweincr,rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwsField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for GWS in noahmp401_descale_tws')

  call ESMF_FieldGet(gwsField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for GWS in noahmp401_descale_tws')
  call ESMF_StateGet(LSM_Incr_State,"Groundwater Storage",gwsIncrField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(gwsIncrField,localDE=0,farrayPtr=gwsincr,rc=status)
  call LIS_verify(status)

  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     swe(t)    = swe(t)!*100.0
     sweincr(t)    = sweincr(t)!*100.0

     gws(t)    = gws(t)!*1000.0
     gwsincr(t)    = gwsincr(t)!*10000.0

  enddo

end subroutine noahmp401_descale_tws
