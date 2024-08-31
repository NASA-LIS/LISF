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
! !ROUTINE: jules5x_setusafsivars
! \label{jules5x_setusafsivars}
!
! !REVISION HISTORY:
! 27Feb2005: Sujay Kumar; Initial Specification
! 25Jun2006: Sujay Kumar: Updated for the ESMF design
! 02 Mar 2010: Sujay Kumar; Modified for Noah 3.1
! 21 Jul 2011: James Geiger; Modified for Noah 3.2
! 05 Nov 2018: Yeosang Yoon; Modified for Jules 5.0 and SNODEP data
! 08 Jul 2019: Yeosang Yoon; Modified for Jules.5.0 and LDT-SI data
! 13 Dec 2019: Eric Kemp; Replaced LDTSI with USAFSI.
! 17 Feb 2020: Yeosang Yoon; Modified for Jules 5.x
!
! !INTERFACE:
subroutine jules5x_setusafsivars(n, LSM_State)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_domain,LIS_surface
  use LIS_snowMod, only : LIS_snow_struc
  use LIS_logMod,  only : LIS_logunit, LIS_verify, LIS_endrun
  use jules5x_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
! 
! !DESCRIPTION:
! 
!  This routine assigns the snow progognostic variables to noah's
!  model space. 
! 
!EOP
  type(ESMF_State)       :: LSM_State
  type(ESMF_Field)       :: sweField
  type(ESMF_Field)       :: snodField
  real, pointer          :: swe(:)
  real, pointer          :: snod(:)
  real                   :: dsneqv,dsnowh
  integer                :: t, pft
  integer                :: status

  call ESMF_StateGet(LSM_State,"SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(LSM_State,"Snowdepth",snodField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snodField,localDE=0,farrayPtr=snod,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

     pft = jules5x_struc(n)%jules5x(t)%pft

     dsneqv = swe(t) - jules5x_struc(n)%jules5x(t)%snow_mass_ij     ! unit: mm
     dsnowh = snod(t) - jules5x_struc(n)%jules5x(t)%snowdepth(pft)  ! unit: m

     ! update snow prognostic variables using jules's snow physics 
     call jules5x_usafsi_update(n, t, dsneqv, dsnowh)

  enddo
  
end subroutine jules5x_setusafsivars

