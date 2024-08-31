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
! !ROUTINE: jules52_qcsoilm
! \label{jules52_qcsoilm}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 20 Dec 2018: Mahdi Navari; Modified for JULES 5.2
!
! !INTERFACE:
subroutine jules52_qcsoilm(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use jules52_lsmMod
  !use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use jules_soil_mod, only:  dzsoil ! EMK
  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: sm1Field
!  type(ESMF_Field)       :: sm2Field
!  type(ESMF_Field)       :: sm3Field
!  type(ESMF_Field)       :: sm4Field
  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
!  real, pointer          :: soilm2(:)
!  real, pointer          :: soilm3(:)
!  real, pointer          :: soilm4(:)
  real                   :: smmax1!,smmax2,smmax3,smmax4
  real                   :: smmin1!,smmin2,smmin3,smmin4
 
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 1 failed in jules52_qcsoilm")
 
  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 1 failed in jules52_qcsoilm")

  call ESMF_AttributeGet(sm1Field,"Max Value",smmax1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Max Value failed in jules52_qcsoilm")
  call ESMF_AttributeGet(sm1Field,"Min Value",smmin1,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Min Value failed in jules52_qcsoilm")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

! MN: the unit of the smmax1, smmin1 are [m3/m3] they are comming from the  jules_sm_attribs.txt 
!     if(soilm1(t) * 1/LIS_sfmodel_struc(n)%lyrthk(1)*1/1000.gt.smmax1) & 
!	soilm1(t) = smmax1 / (1/LIS_sfmodel_struc(n)%lyrthk(1)*1/1000) ! [m3w/m3s] / [[1/m1s]*[1/kg/m3w]]=[kg/m2s]
!     if(soilm1(t) * 1/LIS_sfmodel_struc(n)%lyrthk(1)*1/1000.lt.smmin1) &
!	soilm1(t) = smmin1/ (1/LIS_sfmodel_struc(n)%lyrthk(1)*1/1000)

!EMK...Need to use MKS units.  Since LIS_sfmodel_struc%lyrthk is in cm, we
!use dzoil instead.
     if(soilm1(t) * 1/dzsoil(1)*1/1000.gt.smmax1) & 
	soilm1(t) = smmax1 / (1/dzsoil(1)*1/1000) ! [m3w/m3s] / [[1/m1s]*[1/kg/m3w]]=[kg/m2s]
     if(soilm1(t) * 1/dzsoil(1)*1/1000.lt.smmin1) &
	soilm1(t) = smmin1/ (1/dzsoil(1)*1/1000)
  enddo

end subroutine jules52_qcsoilm

