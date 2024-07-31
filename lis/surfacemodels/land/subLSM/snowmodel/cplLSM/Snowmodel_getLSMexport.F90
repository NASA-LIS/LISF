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
! !ROUTINE: snowmodel_getLSMexport
! \label{snowmodel_getLSMexport}
!
! !REVISION HISTORY:
! 19 Sep 2020: Sujay Kumar; Initial Specification
!  9 Dec 2020: Mahdi Navari; edited to take into account the Crocus slope correction 
!  2 Aug 2021: Kristi Arsenault; Edited to support SnowModel implementation
!
! !INTERFACE:
subroutine snowmodel_getLSMexport(n, SubLSM2LSM_State)
! !USES:

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use snowmodel_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: SubLSM2LSM_State
! 
! !DESCRIPTION:
! 
! 
!EOP
  type(ESMF_Field)   :: snwdField, sweField
  real, pointer      :: swe(:), snwd(:)
  integer            :: t
  integer            :: status
  real               :: tmp_SLOPE
 
  call ESMF_StateGet(SubLSM2LSM_State,"Total SWE",sweField,rc=status)
  call LIS_verify(status,'Snowmodel_getLSMexport: SM SWE state error get')
  call ESMF_StateGet(SubLSM2LSM_State,"Total snowdepth",snwdField,rc=status)
  call LIS_verify(status,'Snowmodel_getLSMexport: SM snowd state error get')

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,'Snowmodel_getLSMexport: SM swe field error get')
  call ESMF_FieldGet(snwdField,localDE=0,farrayPtr=snwd,rc=status)
  call LIS_verify(status,'Snowmodel_getLSMexport: SM snwd field error get')

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
   ! SnowModel states
     swe(t)  = snowmodel_struc(n)%sm(t)%swe_depth   ! SWE:  meters
     snwd(t) = snowmodel_struc(n)%sm(t)%snow_depth  ! Snow depth:  meters
     ! If values are very small -- set to 0.
     if( swe(t) < 0.0001 ) then 
       swe(t) = 0.
     endif
     if( snwd(t) < 0.001 ) then 
       snwd(t) = 0.
     endif 
  enddo

end subroutine snowmodel_getLSMexport


