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
! !ROUTINE: noah39_snow_update
! \label{noah39_snow_update}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
!  02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
!   7 Jun 2021: Mahdi Navari; Modified for Noah39
!
! !INTERFACE:
subroutine noah39_snow_update(n, t, dsneqv, dsnowh)
! !USES:
  use LIS_coreMod, only : LIS_rc
  use noah39_lsmMod
  use LIS_logMod, only : LIS_logunit, LIS_verify

  implicit none
! !ARGUMENTS: 
  integer, intent(in)  :: n
  integer, intent(in)  :: t
  real                 :: dsneqv !mm
  real                 :: dsnowh !m
! 
! !DESCRIPTION:
! 
!  This routine assigns the snow progognostic variables to noah's
!  model space. 
! 
!EOP

     noah39_struc(n)%noah(t)%sneqv = dsneqv ! swe(t)
     noah39_struc(n)%noah(t)%snowh = dsnowh ! snod(t)
end subroutine noah39_snow_update

