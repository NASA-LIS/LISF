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
! !ROUTINE: clm2_setlst
! \label{clm2_setlst}
!
! !REVISION HISTORY:
! 1 Apr 2007: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine clm2_setlst(n, LSM_State)

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
!  sets the land surface temperature related state prognostic variables for
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
 
  call ESMF_StateGet(LSM_State,"Vegetation Temperature",lst1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Ground Surface Temperature",&
       tgField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Soil Temperature 1",&
       ts1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Radiometric Temperature",&
       trField,rc=status)
  call LIS_verify(status)
 

  call ESMF_FieldGet(lst1Field,localDE=0,farrayPtr=lst1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(tgField,localDE=0,farrayPtr=tground,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(ts1Field,localDE=0,farrayPtr=tsoil1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(trField,localDE=0,farrayPtr=trad,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%ntiles(n)
     clm2_struc(n)%clm(t)%t_veg = lst1(t)
     clm2_struc(n)%clm(t)%t_grnd = tground(t)
!     clm2_struc(n)%clm(t)%t_soisno(1) = tsoil1(t)
     clm2_struc(n)%clm(t)%t_rad = trad(t)
  enddo

end subroutine clm2_setlst

