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

module snow_utils_module
USE nrtype
implicit none
private
public::fracliquid
public::templiquid
public::dFracLiq_dTk
public::tcond_snow
contains


 ! ***********************************************************************************************************
 ! public function fracliquid: compute fraction of liquid water
 ! ***********************************************************************************************************
 function fracliquid(Tk,fc_param)
 USE multiconst,only:Tfreeze
 implicit none
 real(dp),intent(in) :: Tk         ! temperature (K)
 real(dp),intent(in) :: fc_param   ! freezing curve parameter (K-1)
 real(dp)            :: fracliquid ! fraction of liquid water (-)
 ! compute fraction of liquid water (-)
 fracliquid = 1._dp / ( 1._dp + (fc_param*( Tfreeze - min(Tk,Tfreeze) ))**2._dp )
 end function fracliquid


 ! ***********************************************************************************************************
 ! public function templiquid: invert the fraction of liquid water function
 ! ***********************************************************************************************************
 function templiquid(fracliquid,fc_param)
 USE multiconst,only:Tfreeze
 implicit none
 real(dp),intent(in) :: fracliquid ! fraction of liquid water (-)
 real(dp),intent(in) :: fc_param   ! freezing curve parameter (K-1)
 real(dp)            :: templiquid ! temperature (K)
 ! compute temperature based on the fraction of liquid water (K)
 templiquid = Tfreeze - ((1._dp/fracliquid - 1._dp)/fc_param**2._dp)**(0.5_dp)
 end function templiquid


 ! ***********************************************************************************************************
 ! public function dFracLiq_dTk: differentiate the freezing curve
 ! ***********************************************************************************************************
 function dFracLiq_dTk(Tk,fc_param)
 USE multiconst,only:Tfreeze
 implicit none
 ! dummies
 real(dp),intent(in) :: Tk           ! temperature (K)
 real(dp),intent(in) :: fc_param     ! freezing curve parameter (K-1)
 real(dp)            :: dFracLiq_dTk ! differentiate the freezing curve (K-1)
 ! locals
 real(dp)            :: Tdep         ! temperature depression (K)
 real(dp)            :: Tdim         ! dimensionless temperature (-)
 ! compute local variables (just to make things more efficient)
 Tdep = Tfreeze - min(Tk,Tfreeze)
 Tdim = fc_param*Tdep
 ! differentiate the freezing curve w.r.t temperature
 dFracLiq_dTk = (fc_param*2._dp*Tdim) / ( ( 1._dp + Tdim**2._dp)**2._dp )
 end function dFracLiq_dTk


 ! ***********************************************************************************************************
 ! public subroutine tcond_snow: compute thermal conductivity of snow
 ! ***********************************************************************************************************
 subroutine tcond_snow(BulkDenIce,thermlcond,err,message)
 USE multiconst,only:lambda_air,lambda_ice  ! thermal conductivity of air and ice
 USE globalData,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 USE mDecisions_module,only:Yen1965,Mellor1977,Jordan1991 ! named variables defining thermal conductivity options
 implicit none
 real(dp),intent(in)      :: BulkDenIce     ! bulk density of ice (kg m-3)
 real(dp),intent(out)     :: thermlcond     ! thermal conductivity of snow (W m-1 K-1)
 integer(i4b),intent(out) :: err            ! error code
 character(*),intent(out) :: message        ! error message
 ! initialize error control
 err=0; message="tcond_snow/"
 ! compute thermal conductivity of snow
 select case(model_decisions(iLookDECISIONS%thCondSnow)%iDecision)
  case(Yen1965);      thermlcond = 3.217d-6 * BulkDenIce**2._dp               ! Yen (1965)
  case(Mellor1977);   thermlcond = 2.576d-6 * BulkDenIce**2._dp + 7.4d-2      ! Mellor (1977)
  case(Jordan1991);   thermlcond = lambda_air + (7.75d-5*BulkDenIce + 1.105d-6*(BulkDenIce**2._dp)) &
                                     * (lambda_ice-lambda_air)                ! Jordan (1991)
  case default
   err=10; message=trim(message)//"unknownOption"; return
 end select
 end subroutine tcond_snow


end module snow_utils_module
