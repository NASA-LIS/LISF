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
! !ROUTINE: clsmf25_gettws
! \label{clsmf25_gettws}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 28Sep2011: Ben Zaitchik: Applied to GRACE
!
! !INTERFACE:
subroutine clsmf25_gettws(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,  only  : LIS_verify
  use clsmf25_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the tws related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  type(ESMF_Field)       :: catdefField
  type(ESMF_Field)       :: srfexcField
  type(ESMF_FIeld)       :: rzexcField
  type(ESMF_Field)       :: wesn1Field
  type(ESMF_Field)       :: wesn2Field
  type(ESMF_Field)       :: wesn3Field
  integer                :: t
  integer                :: status
  real, pointer          :: catdef(:)
  real, pointer          :: srfexc(:)
  real, pointer          :: rzexc(:)
  real, pointer          :: wesn1(:)
  real, pointer          :: wesn2(:)
  real, pointer          :: wesn3(:)        
 
  call ESMF_StateGet(LSM_State,"Catchment Deficit",catdefField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Surface Excess",srfexcField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Root Zone Excess",rzexcField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 1",wesn1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 2",wesn2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Water Equivalent Snow 3",wesn3Field,rc=status)
  call LIS_verify(status)       
    
  call ESMF_FieldGet(catdefField,localDE=0,farrayPtr=catdef,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(srfexcField,localDE=0,farrayPtr=srfexc,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(rzexcField,localDE=0,farrayPtr=rzexc,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn1Field,localDE=0,farrayPtr=wesn1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn2Field,localDE=0,farrayPtr=wesn2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(wesn3Field,localDE=0,farrayPtr=wesn3,rc=status)
  call LIS_verify(status)      

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
    catdef(t) = clsmf25_struc(n)%cat_progn(t)%catdef
    srfexc(t) = clsmf25_struc(n)%cat_progn(t)%srfexc
    rzexc(t) = clsmf25_struc(n)%cat_progn(t)%rzexc
    wesn1(t) = clsmf25_struc(n)%cat_progn(t)%wesn(1)
    wesn2(t) = clsmf25_struc(n)%cat_progn(t)%wesn(2)
    wesn3(t) = clsmf25_struc(n)%cat_progn(t)%wesn(3)

!    if(LIS_localPet.eq.24.and.t.eq.7311) then 
!       print*, 'get ',t, wesn1(t), wesn2(t), wesn3(t)
!     endif
  enddo

end subroutine clsmf25_gettws

