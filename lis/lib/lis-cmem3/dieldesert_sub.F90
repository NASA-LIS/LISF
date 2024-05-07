!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!!!! Insert banner here !!!!
!---------------------------------------------------------------------------
! Calculate the dielectric constant of desert sand 
!  Matzler, 1998, IEEE TGARS, 36, 317-319. 
!---------------------------------------------------------------------------

!SUBROUTINE diel_desert 

! Purpose :
!   Calculate the dielectric constant of desert sand 

! Reference:
!   Matzler, 1998, IEEE TGARS, 36, 317-319.

! Input Variables:
!   fghz	 frequency (GHz) 

! Output variables: 
!   eps		 dielectric constant of desert sand 

!---------------------------------------------------------------------------

SUBROUTINE diel_desert (fghz, eps)

IMPLICIT NONE

real :: fghz
COMPLEX :: eps

REAL, parameter :: eInf = 2.53    
REAL, parameter :: eS   = 2.79 
REAL, parameter :: f0 = 0.27    ! Ghz 
REAL, parameter :: app = 0.002    ! a'' in Eq. (1) 
COMPLEX, parameter  ::  ci = (0, 1.0) 
!---------------------------------------------------------------------------

eps = eInf + (eS - eInf ) / ( 1 - ci * fghz / f0 ) + ci * app         ! Eq. (1) 

!write(*, *) "====", real(eps), aimag(eps)

return 
END SUBROUTINE diel_desert 

