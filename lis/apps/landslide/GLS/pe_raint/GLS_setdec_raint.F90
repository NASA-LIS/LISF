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
! !ROUTINE: GLS_setdec_raint
!  \label{GLS_setdec_raint}
!
! !REVISION HISTORY:
! 06 Feb 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine GLS_setdec_raint(DEC_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_soilsMod,  only : LIS_soils
  use LIS_logMod,       only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_verify, LIS_endrun
  use GLSMod,           only : GLS_data

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: DEC_State
!
! !DESCRIPTION:
!  
!  This routine assigns the decision space to GLS  model variables. 
! 
!EOP
  integer                :: n
  type(ESMF_Field)       :: rft1Field,rft2Field,rft3Field,stField,siField
  real, pointer          :: rft1(:),rft2(:),rft3(:),st(:),si(:)
  integer                :: t,g,m
  integer                :: status


  n = 1
  call ESMF_StateGet(DEC_State,"RF_THRESHOLD1",rft1Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(DEC_State,"RF_THRESHOLD2",rft2Field,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(DEC_State,"RF_THRESHOLD3",rft3Field,rc=status)
  call LIS_verify(status)
!  call ESMF_StateGet(DEC_State,"SUSCEP_INDEX",siField,rc=status)
!  call LIS_verify(status)
!  call ESMF_StateGet(DEC_State,"SUSCEP_INDEX_THRESHOLD",stField,rc=status)
!  call LIS_verify(status)

  call ESMF_FieldGet(rft1Field,localDE=0,farrayPtr=rft1,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(rft2Field,localDE=0,farrayPtr=rft2,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(rft3Field,localDE=0,farrayPtr=rft3,rc=status)
  call LIS_verify(status)
!  call ESMF_FieldGet(siField,localDE=0,farrayPtr=si,rc=status)
!  call LIS_verify(status)
!  call ESMF_FieldGet(stField,localDE=0,farrayPtr=st,rc=status)
!  call LIS_verify(status)

  do g=1,LIS_rc%ntiles(n)/LIS_rc%nensem(n)
     do m=1,LIS_rc%nensem(n)
        t = (g-1)*LIS_rc%nensem(n)+m
        GLS_data(m)%rf_threshold1 = rft1(t)
        GLS_data(m)%rf_threshold2 = rft2(t)
        GLS_data(m)%rf_threshold3 = rft3(t)
!        GLS_data(m)%suscep_index_threshold = st(t)
!        print*, 'param ',t, rft1(t), rft2(t), rft3(t)
!        GLS_data%suscep_index(t)  = si(t)
     enddo
  enddo

end subroutine GLS_setdec_raint



