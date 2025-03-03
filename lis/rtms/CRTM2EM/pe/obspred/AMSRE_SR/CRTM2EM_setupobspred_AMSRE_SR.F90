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
! !ROUTINE: CRTM2EM_setupobspred_AMSRE_SR
!  \label{CRTM2EM_setupobspred_AMSRE_SR}
!
! !REVISION HISTORY:
! 16 Jul 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine CRTM2EM_setupobspred_AMSRE_SR(OBSPred)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_vecTile
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: OBSPred
!
! !DESCRIPTION:
!  
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: emField
  integer                :: status
  character*2, parameter :: fnames(6)=(/'07', '11', '19', '24', '37', '89'/)
  character*1, parameter :: pnames(2)=(/'V','H'/)
  integer                   ::  f !freq counter
  integer                   ::  p !polarization counter
  integer                   :: numfreqs=6        
  integer                   :: numchannels=12    
  integer                   :: numpolarizations=2

  n = 1


  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

  do f=1, numfreqs
     do p=1, numpolarizations 
        emField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecTile(n), &
             name="Emissivity" // fnames(f) // pnames(p), rc=status)
        call LIS_verify(status)
        
        call ESMF_StateAdd(OBSPred,(/emField/),rc=status)
        call LIS_verify(status)
     enddo
  enddo
end subroutine CRTM2EM_setupobspred_AMSRE_SR

