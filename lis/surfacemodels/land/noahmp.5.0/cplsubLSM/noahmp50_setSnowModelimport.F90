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
! !ROUTINE: noahmp50_setSnowModelimport
! \label{noahmp50_setSnowModelimport}
!
! !REVISION HISTORY:
! 19 Sep 2020: Sujay Kumar; Initial Specification
! 12 Aug 2021: Kristi Arsenault; Added SnowModel
! May 2023: Cenlin He; modified for refactored NoahMP v5 and later
!
! !INTERFACE:
subroutine noahmp50_setSnowModelimport(n, SubLSM2LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use NoahMP50_lsmMod

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
  real               :: dsneqv,dsnowh
  integer            :: t
  integer            :: status

  call ESMF_StateGet(SubLSM2LSM_State,"Total SWE",sweField,rc=status)
  call LIS_verify(status,"noahmp50_setSnowModelimport: error in swe state get from SnowModel")
  call ESMF_StateGet(SubLSM2LSM_State,"Total snowdepth",snwdField,rc=status)
  call LIS_verify(status,"noahmp50_setSnowModelimport: error in snwd state get from SnowModel")

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status,"noahmp50_setSnowModelimport: error in swe data values")
  call ESMF_FieldGet(snwdField,localDE=0,farrayPtr=snwd,rc=status)
  call LIS_verify(status,"noahmp50_setSnowModelimport: error in snwd data values")

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     ! SnowModel has its SWE in meters -- conversion here to mm:
     dsneqv = (swe(t)*1000.) - NoahMP50_struc(n)%noahmp50(t)%sneqv   !in mm
     dsnowh = snwd(t) - NoahMP50_struc(n)%noahmp50(t)%snowh  !in m

!     if( dsneqv > 0. .or. dsnowh > 0. ) then
!        write(501,*) t, swe(t), noahmp50_struc(n)%noahmp50(t)%sneqv, &
!                        snwd(t), noahmp50_struc(n)%noahmp50(t)%snowh 
!     endif

     ! Update NoahMP's Snow states:
     call noahmp50_snow_update(n, t, dsneqv, dsnowh)

  enddo

end subroutine noahmp50_setSnowModelimport


