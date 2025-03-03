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

!SUBROUTINE CMEM_VEG

! Purpose :
! -------
!  Vegetation model

! Author :
! ------
!   23-Aug-2006 Thomas Holmes   *ECMWF*

! Input Variables:
!   fghz         frequency (GHz)
!   theta        zenith angle (degree)
!   vtype        vegetation cover type, real (UMD class numbers, 0-13)
!   t_veg        vegetation temperature (K)
!   h_veg        vegetation height (m) 
!   lai          leaf area index 
!   wc_veg       vegetation water content per unit area (kg/m2)
!   rho_veg      vegetation density (kg/m3)
!   m_d          dry mass fraction of vegetation
!   d_leaf       leaf thickness (mm)
!
! Output Variables:
!   tau_veg      microwave emission from vegetation (K) (index 1.h-pol, 2.v-pol)
!   tb_veg       microwave emission from vegetation (K) (index 1.h-pol, 2.v-pol)
!
! Internal variables 
!   omega	 single scattering albedo of vegetation (index 1.h-pol, 2.v-pol)
!   omega_eff	 effective single scattering albedo of vegetation (1-h, 2-v) 
!   a_geo        vegetation structure coefficient (index 1.h-pol, 2.v-pol)
!
! CMEM default values
! rho_veg = 950.
! m_d = 0.3
! d_leaf = 0.2
!---------------------------------------------------------------------------

SUBROUTINE CMEM_VEG(fghz, theta, vtype, t_veg, h_veg, lai, &
                    wc_veg,                &    ! veg params
                    rho_veg, m_d, d_leaf,  &
	            tau_veg, tb_veg) 

IMPLICIT NONE

real, parameter :: tfreeze = 273.15
real, parameter :: sal_vw = 6.0 
real :: fghz, theta
integer :: vtype
real :: t_veg, h_veg, lai, wc_veg, rho_veg, m_d, d_leaf
real :: tau_veg(2), tb_veg(2)

real :: omega(2), omega_eff(2), a_geo(2)

complex :: eps_vw

INTEGER :: i
INTEGER :: medium, isal
!---------------------------------------------------------------------------

if (vtype .GE. 12 ) then   ! 12 and 13 are bare land
  tau_veg = 0.0
  tb_veg = 0.0 
  return 
end if 

call get_a_geo(theta, vtype, a_geo)

    ! 2. Calculate vegetation opacity (tau_veg)
    !    -------------------------------------
    medium = 0
    isal = 0
    CALL DIEL_WAT (fghz, (t_veg-tfreeze), sal_vw, 0.0, 0.0, 0.0, medium, isal, eps_vw) 
    CALL VEGWEGM(fghz, theta, &            ! sensor params
                   eps_vw, wc_veg, a_geo, &    ! veg params
                   rho_veg, m_d, d_leaf,  &
                   tau_veg, omega)                  ! output

    call get_omega_eff(fghz, vtype, omega, omega_eff)

    ! 3. Calculate microwave emission from vegetation (tb_veg)
    !    -----------------------------------------------------
    DO i = 1, 2
      tb_veg(i) = (1. - omega_eff(i)) * (1. - exp(-tau_veg(i))) * t_veg
    ENDDO

END SUBROUTINE CMEM_VEG

! Compute single scattering albedo of vegetation (index 1.h-pol, 2.v-pol)
! 2/19/2011, Y.Tian: added weak frequency-dependency

! Input variable:
! fghz          frequency (GHz)
! vtype		UMD land cover type (0-13)

! Output variable: 
! w_eff 	single scattering albedo of vegetation (index 1.h-pol, 2.v-pol) 

subroutine get_omega_eff(fghz, vtype, omega, omega_eff) 
 implicit none
 integer vtype 
 real :: fghz, omega(2), omega_eff(2), Zw_eff(2, 0:13) 
 
! tentative numbers 
Zw_eff(:, 0)=(0.0,0.0)          ! 0     WATER
Zw_eff(:, 1)=(0.07,0.07)        ! 1     EVERGREEN NEEDLELEAF FOREST
Zw_eff(:, 2)=(0.08,0.08)        ! 2     EVERGREEN BROADLEAF FOREST
Zw_eff(:, 3)=(0.095,0.095)      ! 3     DECIDUOUS NEEDLELEAF FOREST 
Zw_eff(:, 4)=(0.05,0.05)        ! 4     DECIDUOUS BROADLEAF FOREST 
Zw_eff(:, 5)=(0.05,0.05)        ! 5     MIXED FOREST 
Zw_eff(:, 6)=(0.05,0.05)          ! 6     WOODLAND 
Zw_eff(:, 7)=(0.05,0.05)        ! 7     WOODED GRASSLAND 
Zw_eff(:, 8)=(0.05,0.05)        ! 8     CLOSED SHRUBLAND 
Zw_eff(:, 9)=(0.05,0.05)        ! 9     OPEN SHRUBLAND 
Zw_eff(:, 10)=(0.05,0.05)       ! 10     GRASSLAND 
Zw_eff(:, 11)=(0.05,0.05)       ! 11     CROPLAND 
Zw_eff(:, 12)=(0.0,0.0)         ! 12     BARE GROUND 
Zw_eff(:, 13)=(0.0,0.0)         ! 13     URBAN AND BUILT-UP 

  select case (vtype)
    case (0:13)
      omega_eff(:) = Zw_eff(:, vtype) 
    case default
      omega_eff(:) = (0.0, 0.0) 
  end select

return 
end subroutine get_omega_eff

! Compute vegetation structure coefficient 

! Input variable:
! theta         zenith angle (degree)
! vtype         UMD land cover type (0-13)

! Output variable:
! a_geo          vegetation structure coefficient (index 1.h-pol, 2.v-pol)
!---------------------------------------------------------------------------
!  About vegetation structure coefficient a_geo:
!
!    a_geo(pol) = u(pol) / 3
!
!    u takes vegetation geometry into account
!
!    Theoretical examples:
!       horizontal leaves:  a_geo(h) = 1     a_geo(v) = (cos(theta))**2
!       isotropic leaves:   a_geo(h) = 2/3   a_geo(v) = 2/3
!       vertical stalks:    a_geo(h) = 0     a_geo(v) = (sin(theta))**2
!       isotropic stalks:   a_geo(h) = 1/3   a_geo(v) = 1/3
!---------------------------------------------------------------------------

subroutine get_a_geo(theta, vtype, a_geo)
 implicit none
 real :: theta, c2, s2
 integer :: vtype
 real :: pi,  a_geo(2), Za_geo(2, 0:13)

 pi = acos(-1.0)
 c2 = cos(theta*pi/180.0)**2
 s2 = sin(theta*pi/180.0)**2

! tentative numbers 
Za_geo(:, 0)=(0.0,0.0)          ! 0     WATER
Za_geo(:, 1)=(0.66, 0.66)       ! 1     EVERGREEN NEEDLELEAF FOREST
Za_geo(1, 2)=1.00               ! 2     EVERGREEN BROADLEAF FOREST
Za_geo(2, 2)=c2                
Za_geo(:, 3)=(0.66, 0.66)       ! 3     DECIDUOUS NEEDLELEAF FOREST
Za_geo(1, 4)=1.00               ! 4     DECIDUOUS BROADLEAF FOREST
Za_geo(2, 4)=c2          
Za_geo(:, 5)=(0.66, 0.66)       ! 5     MIXED FOREST
Za_geo(:, 6)=(0.66, 0.66)       ! 6     WOODLAND
Za_geo(1, 7)=1.00               ! 7     WOODED GRASSLAND
Za_geo(2, 7)=c2  
Za_geo(:, 8)=(0.33, 0.33)       ! 8     CLOSED SHRUBLAND
Za_geo(:, 9)=(0.33, 0.33)       ! 9     OPEN SHRUBLAND
Za_geo(1,10)=1.00               ! 10     GRASSLAND
Za_geo(2,10)=c2 
Za_geo(1,11)=1.00               ! 11     CROPLAND
Za_geo(2,11)=c2  
Za_geo(:,12)=(0.0,0.0)         ! 12     BARE GROUND
Za_geo(:,13)=(0.0,0.0)         ! 13     URBAN AND BUILT-UP

  select case (vtype)
    case (0:13)
      a_geo(:) = Za_geo(:, vtype)
    case default
      a_geo(:) = (0.0, 0.0)
  end select

return
end subroutine get_a_geo

