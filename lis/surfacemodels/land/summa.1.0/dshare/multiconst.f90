! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE multiconst
 USE nrtype
 ! define physical constants
 REAL(DP), PARAMETER           :: ave_slp        =  101325.0_dp ! mean sea level pressure              (Pa)
 REAL(DP), PARAMETER           :: vkc            =  0.4_dp      ! von Karman constant                  (-)
 REAL(DP), PARAMETER           :: satvpfrz       =  610.8_dp    ! sat vapour pressure at 273.16K       (Pa)
 REAL(DP), PARAMETER           :: w_ratio        =  0.622_dp    ! molecular ratio water to dry air     (-)
 REAL(DP), PARAMETER           :: R_da           =  287.053_dp  ! gas constant for dry air             (Pa K-1 m3 kg-1; J kg-1 K-1)
 REAL(DP), PARAMETER           :: R_wv           = 461.285_dp   ! gas constant for water vapor         (Pa K-1 m3 kg-1; J kg-1 K-1)
 !REAL(DP), PARAMETER           :: Rgas           =  8.314462_dp ! universal gas constant               (J mol-1 K-1)
 ! use same digits as Noah-MP -- chasing the difference
 REAL(DP), PARAMETER           :: Rgas           = 8.314_dp     ! universal gas constant               (J mol-1 K-1)
 REAL(DP), PARAMETER           :: gravity        = 9.80616_dp   ! acceleration of gravity              (m s-2)
 REAL(DP), PARAMETER           :: Cp_air         = 1005._dp     ! specific heat of air                 (J kg-1 K-1)
 REAL(DP), PARAMETER           :: Cp_ice         = 2114._dp     ! specific heat of ice                 (J kg-1 K-1)
 REAL(DP), PARAMETER           :: Cp_soil        = 850._dp      ! specific heat of soil                (J kg-1 K-1)
 REAL(DP), PARAMETER           :: Cp_water       = 4181._dp     ! specific heat of liquid water        (J kg-1 K-1)
 REAL(DP), PARAMETER           :: Tfreeze        = 273.16_dp    ! temperature at freezing              (K)
 REAL(DP), PARAMETER           :: TriplPt        = 273.16_dp    ! triple point of water                (K)
 REAL(DP), PARAMETER           :: LH_fus         = 333700.0_dp  ! latent heat of fusion                (J kg-1)
 REAL(DP), PARAMETER           :: LH_vap         = 2501000.0_dp ! latent heat of vaporization          (J kg-1)
 REAL(DP), PARAMETER           :: LH_sub         = 2834700.0_dp ! latent heat of sublimation           (J kg-1)
 REAL(DP), PARAMETER           :: sb             = 5.6705d-8    ! Stefan Boltzman constant             (W m-2 K-4)
 REAL(DP), PARAMETER           :: em_sno         = 0.99_dp      ! emissivity of snow                   (-)
 REAL(DP), PARAMETER           :: lambda_air     = 0.026_dp     ! thermal conductivity of air          (W m-1 K-1)
 REAL(DP), PARAMETER           :: lambda_ice     = 2.50_dp      ! thermal conductivity of ice          (W m-1 K-1)
 REAL(DP), PARAMETER           :: lambda_water   = 0.60_dp      ! thermal conductivity of liquid water (W m-1 K-1)
 REAL(DP), PARAMETER           :: iden_air       = 1.293_dp     ! intrinsic density of air             (kg m-3)
 REAL(DP), PARAMETER           :: iden_ice       = 917.0_dp     ! intrinsic density of ice             (kg m-3)
 REAL(DP), PARAMETER           :: iden_water     = 1000.0_dp    ! intrinsic density of liquid water    (kg m-3)
 REAL(DP), PARAMETER           :: secprday       = 86400._dp    ! number of seconds in a day
 REAL(DP), PARAMETER           :: secprhour      = 3600._dp     ! number of seconds in an hour
 REAL(DP), PARAMETER           :: secprmin       = 60._dp       ! number of seconds in a minute
 REAL(DP), PARAMETER           :: minprhour      = 60._dp       ! number of minutes in an hour

! define missing values
real(dp),parameter,public      :: quadMissing    = -9999._qp    ! missing quadruple precision number
real(dp),parameter,public      :: realMissing    = -9999._dp    ! missing double precision number
integer(i4b),parameter,public  :: integerMissing = -9999        ! missing integer 

END MODULE multiconst
