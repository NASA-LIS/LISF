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
! !ROUTINE: noahmpnew_scale_tws
! \label{noahmpnew_scale_tws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li; created for Noah-MP4.0.1
! May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!
! !INTERFACE:
subroutine noahmpnew_scale_tws(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noahmpnew_lsmMod
  use LIS_constantsMod,  only : LIS_CONST_RHOFW ! Natt

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Scales twsoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP

!Wanshu
 
  type(ESMF_Field)     :: gwField
  real, pointer        :: gws(:)
  integer              :: t
  integer              :: status
  
! Natt
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)


#if 0
  ! Natt
  ! Scale TWS states to mm (Note, GWS is already in mm)
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm1 in noahmpnew_scale_tws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm2 in noahmpnew_scale_tws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm3 in noahmpnew_scale_tws')
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for sm4 in noahmpnew_scale_tws')
  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for SWE in noahmpnew_scale_tws')
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for Snowdepth in noahmpnew_scale_tws')
  
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm1 in noahmpnew_scale_tws')
  call ESMF_FieldGet(sm2Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm2 in noahmpnew_scale_tws')
  call ESMF_FieldGet(sm3Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm3 in noahmpnew_scale_tws')
  call ESMF_FieldGet(sm4Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for sm4 in noahmpnew_scale_tws')
  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for SWE in noahmpnew_scale_tws')
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for Snowdepth in noahmpnew_scale_tws')
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     soilm1(t) = soilm1(t) * NoahMPnew_struc(n)%sldpth(1)*1000.0 ! m3m-3 -> mm
     soilm2(t) = soilm2(t) * NoahMPnew_struc(n)%sldpth(2)*1000.0 ! m3m-3 -> mm
     soilm3(t) = soilm3(t) * NoahMPnew_struc(n)%sldpth(3)*1000.0 ! m3m-3 -> mm
     soilm4(t) = soilm4(t) * NoahMPnew_struc(n)%sldpth(4)*1000.0 ! m3m-3 -> mm
     swe(t)    = swe(t) * 1000.0  ! m -> mm
     snod(t)   = snod(t) * 1000.0 ! m -> mm
  enddo
#endif  
  
end subroutine noahmpnew_scale_tws
