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
! !ROUTINE: noahmp36_getsws_hymap2
!  \label{noahmp36_getsws_hymap2}
!
! !REVISION HISTORY:
! 12 Sep 2019: Augusto Getirana; implementation of two-way coupling
!
! !INTERFACE:
subroutine noahmp36_getsws_hymap2(n)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_routingMod, only : LIS_runoff_state
  use LIS_logMod
  use LIS_historyMod
  use noahmp36_lsmMod, only : noahmp36_struc

  implicit none
! !ARGUMENTS: 
  integer,  intent(in)   :: n 
!
! !DESCRIPTION:
!   This routine defines the surface water storage variables in NoahMP
!   to be updated based on feedback from HYMAP2
!  
!EOP
  type(ESMF_Field)       :: rivsto_field
  type(ESMF_Field)       :: fldsto_field
  type(ESMF_Field)       :: fldfrc_field
  real, pointer          :: rivstotmp(:)
  real, pointer          :: fldstotmp(:)
  real, pointer          :: fldfrctmp(:)
  integer                :: t
  integer                :: c,r
  integer                :: status
  integer                :: enable2waycpl
  
  call ESMF_AttributeGet(LIS_runoff_state(n),"2 way coupling",&
       enable2waycpl, rc=status)
  call LIS_verify(status)

  if(enable2waycpl==1) then 
     ! River Storage
     call ESMF_StateGet(LIS_runoff_state(n),"River Storage",rivsto_field,rc=status)
     call LIS_verify(status,'ESMF_StateGet failed for River Storage')
     
     call ESMF_FieldGet(rivsto_field,localDE=0,farrayPtr=rivstotmp,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for River Storage')
     where(rivstotmp/=LIS_rc%udef) &
          NOAHMP36_struc(n)%noahmp36(:)%rivsto=rivstotmp/NOAHMP36_struc(n)%dt

     ! Flood Storage
     call ESMF_StateGet(LIS_runoff_state(n),"Flood Storage",fldsto_field,rc=status)
     call LIS_verify(status,'ESMF_StateGet failed for Flood Storage')
     
     call ESMF_FieldGet(fldsto_field,localDE=0,farrayPtr=fldstotmp,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for Flood Storage')
     where(fldstotmp/=LIS_rc%udef)&
          NOAHMP36_struc(n)%noahmp36(:)%fldsto=fldstotmp/NOAHMP36_struc(n)%dt

     ! Flooded Fraction Flag
     call ESMF_StateGet(LIS_runoff_state(n),"Flooded Fraction",fldfrc_field,rc=status)
     call LIS_verify(status,'ESMF_StateGet failed for Flooded Fraction')
     
     call ESMF_FieldGet(fldfrc_field,localDE=0,farrayPtr=fldfrctmp,rc=status)
     call LIS_verify(status,'ESMF_FieldGet failed for Flooded Fraction')
     NOAHMP36_struc(n)%noahmp36(:)%fldfrc=fldfrctmp
  endif  

end subroutine noahmp36_getsws_hymap2
