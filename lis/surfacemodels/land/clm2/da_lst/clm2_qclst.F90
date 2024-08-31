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
! !ROUTINE: clm2_qclst
! \label{clm2_qclst}
!
! !REVISION HISTORY:
! 1 Apr 2007: Sujay Kumar: Initial Specification
!
! !INTERFACE:
subroutine clm2_qclst(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_verify
  use clm2_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Qc's the LST related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[LSM\_State] ESMF State container for LSM state variables
!  \end{description}
!EOP
  type(ESMF_Field)       :: lst1Field
  type(ESMF_Field)       :: tgField
  type(ESMF_Field)       :: ts1Field 
  type(ESMF_Field)       :: trField 

  integer                :: t
  integer                :: status
  
  real, pointer          :: lst1(:)
  real, pointer          :: tground(:)
  real, pointer          :: tsoil1(:)
  real, pointer          :: trad(:)

  real                   :: lstmax1
  real                   :: lstmin1  
  real                   :: tgmax1
  real                   :: tgmin1 
  real                   :: ts1max1
  real                   :: ts1min1
  real                   :: trmax1
  real                   :: trmin1
  
  
  call ESMF_StateGet(LSM_State,"Vegetation Temperature",lst1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Ground Surface Temperature",tgField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Temperature 1",ts1Field,rc=status)
  call LIS_verify(status)  
  call ESMF_StateGet(LSM_State,"Radiometric Temperature",trField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lst1Field,localDE=0,farrayPtr=lst1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(tgField,localDE=0,farrayPtr=tground,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(ts1Field,localDE=0,farrayPtr=tsoil1,rc=status)
  call LIS_verify(status)  
  call ESMF_FieldGet(trField,localDE=0,farrayPtr=trad,rc=status)
  call LIS_verify(status)

 
  call ESMF_AttributeGet(lst1Field,"Max Value",lstmax1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(lst1Field,"Min Value",lstmin1,rc=status)
  call LIS_verify(status)

  call ESMF_AttributeGet(tgField,"Max Value",tgmax1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(tgField,"Min Value",tgmin1,rc=status)
  call LIS_verify(status)
  
  call ESMF_AttributeGet(ts1Field,"Max Value",ts1max1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(ts1Field,"Min Value",ts1min1,rc=status)
  call LIS_verify(status)
  
  call ESMF_AttributeGet(trField,"Max Value",trmax1,rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(trField,"Min Value",trmin1,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     if(lst1(t).gt.lstmax1) lst1(t) = lstmax1
     if(lst1(t).lt.lstmin1) lst1(t) = lstmin1
     
     if(tground(t).gt.tgmax1) tground(t) = tgmax1
     if(tground(t).lt.tgmin1) tground(t) = tgmin1
     
     if(tsoil1(t).gt.ts1max1) tsoil1(t) = ts1max1
     if(tsoil1(t).lt.ts1min1) tsoil1(t) = ts1min1
     
     if(trad(t).gt.trmax1) trad(t) = trmax1
     if(trad(t).lt.trmin1) trad(t) = trmin1
     
  enddo

end subroutine clm2_qclst
