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

module convE2Temp_module
USE nrtype
implicit none
private
public::E2T_lookup
public::E2T_nosoil
public::temp2ethpy
! define the look-up table used to compute temperature based on enthalpy
integer(i4b),parameter            :: nlook=10001       ! number of elements in the lookup table
real(dp),dimension(nlook),public  :: E_lookup          ! enthalpy values (J kg-1)
real(dp),dimension(nlook),public  :: T_lookup          ! temperature values (K)
contains


 ! ************************************************************************************************************************
 ! public subroutine E2T_lookup: define a look-up table to compute specific enthalpy based on temperature, assuming no soil
 ! ************************************************************************************************************************
 subroutine E2T_lookup(mpar_data,err,message)
 USE nr_utility_module,only:arth                       ! use to build vectors with regular increments
 USE spline_int_module,only:spline,splint              ! use for cubic spline interpolation
 USE multiconst,only:Tfreeze                           ! freezing point (K)
 USE var_lookup,only:iLookPARAM                        ! named variables to define structure element
 USE data_types,only:var_dlength                       ! data vector with variable length dimension (dp): x%var(:)%dat(:)
 implicit none
 ! declare dummy variables
 type(var_dlength),intent(in)  :: mpar_data            ! model parameters
 integer(i4b),intent(out)      :: err                  ! error code
 character(*),intent(out)      :: message              ! error message
 ! declare local variables
 character(len=128)            :: cmessage             ! error message in downwind routine
 real(dp),parameter            :: T_start=260.0_dp     ! start temperature value where all liquid water is assumed frozen (K)
 real(dp)                      :: T_incr,E_incr        ! temperature/enthalpy increments
 real(dp),dimension(nlook)     :: Tk                   ! initial temperature vector
 real(dp),dimension(nlook)     :: Ey                   ! initial enthalpy vector
 real(dp),parameter            :: waterWght=1._dp      ! weight applied to total water (kg m-3) --- cancels out
 real(dp),dimension(nlook)     :: T2deriv              ! 2nd derivatives of the interpolating function at tabulated points
 integer(i4b)                  :: ilook                ! loop through lookup table
 ! initialize error control
 err=0; message="E2T_lookup/"
 ! associate
 associate( snowfrz_scale => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1) )
 ! define initial temperature vector
 T_incr = (Tfreeze - T_start) / real(nlook-1, kind(dp))  ! temperature increment
 Tk     = arth(T_start,T_incr,nlook)
 ! ***** compute specific enthalpy (NOTE: J m-3 --> J kg-1) *****
 do ilook=1,nlook
  Ey(ilook) = temp2ethpy(Tk(ilook),waterWght,snowfrz_scale)/waterWght  ! (J m-3 --> J kg-1)
 end do
 ! define the final enthalpy vector
 E_incr   = (-Ey(1)) / real(nlook-1, kind(dp))  ! enthalpy increment
 E_lookup = arth(Ey(1),E_incr,nlook)
 ! use cubic spline interpolation to obtain temperature values at the desired values of enthalpy
 call spline(Ey,Tk,1.e30_dp,1.e30_dp,T2deriv,err,cmessage)  ! get the second derivatives
 if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
 do ilook=1,nlook
  call splint(Ey,Tk,T2deriv,E_lookup(ilook),T_lookup(ilook),err,cmessage)
  if(err/=0) then; message=trim(message)//trim(cmessage); return; end if
  !write(*,'(i6,1x,2(f20.4,1x))') ilook, E_lookup(ilook), T_lookup(ilook)
 end do
 end associate
 end subroutine E2T_lookup


 ! ************************************************************************************************************************
 ! public subroutine E2T_nosoil: compute temperature based on specific enthalpy -- appropriate when no dry mass, as in snow
 ! ************************************************************************************************************************
 subroutine E2T_nosoil(Ey,BulkDenWater,fc_param,Tk,err,message)
 ! compute temperature based on enthalpy -- appropriate when no dry mass, as in snow
 USE multiconst, only: Cp_ice ! specific heat of ice (J kg-1 K-1)
 implicit none
 ! declare dummy variables
 real(dp),intent(in)      :: Ey            ! total enthalpy (J m-3)
 real(dp),intent(in)      :: BulkDenWater  ! bulk density of water (kg m-3)
 real(dp),intent(in)      :: fc_param      ! freezing curve parameter (K-1)
 real(dp),intent(out)     :: Tk            ! initial temperature guess / final temperature value (K)
 integer(i4b),intent(out) :: err           ! error code
 character(*),intent(out) :: message       ! error message
 ! declare local variables
 real(dp),parameter       :: dx=1.d-8      ! finite difference increment (J kg-1)
 real(dp),parameter       :: atol=1.d-12   ! convergence criteria (J kg-1)
 real(dp)                 :: E_spec        ! specific enthalpy (J kg-1)
 real(dp)                 :: E_incr        ! enthalpy increment
 integer(i4b)             :: niter=15      ! maximum number of iterations
 integer(i4b)             :: iter          ! iteration index
 integer(i4b)             :: i0            ! position in lookup table
 real(dp)                 :: Tg0,Tg1       ! trial temperatures (K)
 real(dp)                 :: Ht0,Ht1       ! specific enthalpy, based on the trial temperatures (J kg-1)
 real(dp)                 :: f0,f1         ! function evaluations (difference between enthalpy guesses)
 real(dp)                 :: dh            ! enthalpy derivative
 real(dp)                 :: dT            ! temperature increment
 ! initialize error control
 err=0; message="E2T_nosoil/"
 ! convert input of total enthalpy (J m-3) to total specific enthalpy (J kg-1)
 E_spec = Ey/BulkDenWater ! (NOTE: no soil)
 !write(*,'(a,1x,10(e20.10,1x))') 'E_spec, E_lookup(1)', E_spec, E_lookup(1)

 ! ***** get initial guess and derivative assuming all water is frozen
 if(E_spec<E_lookup(1))then ! process cases below the limit of the look-up tab;e
  ! get temperature guess
  Tg0 = (E_spec - E_lookup(1))/Cp_ice + T_lookup(1)
  Tg1 = Tg0+dx
  ! compute enthalpy
  Ht0 = temp2ethpy(Tg0,1._dp,fc_param)
  Ht1 = temp2ethpy(Tg1,1._dp,fc_param)
  ! compute function evaluations
  f0  = Ht0 - E_spec
  f1  = Ht1 - E_spec

 ! ***** get initial guess and derivative from the look-up table
 else
  ! get enthalpy increment
  E_incr = E_lookup(2) - E_lookup(1)
  ! get position in lookup table
  i0 = ceiling( (E_spec - E_lookup(1)) / E_incr, kind(i4b) )
  ! check found the appropriate value in the look-up table
  if(E_spec < E_lookup(i0) .or. E_spec > E_lookup(i0+1) .or. &
     i0 < 1 .or. i0+1 > nlook)then
   err=10; message=trim(message)//'problem finding appropriate value in lookup table'; return
  end if
  ! get temperature guess
  Tg0 = T_lookup(i0)
  Tg1 = T_lookup(i0+1)
  ! compute function evaluations
  f0  = E_lookup(i0) - E_spec
  f1  = E_lookup(i0+1) - E_spec
 end if

 ! compute initial derivative
 dh  = (f1 - f0) / (Tg1 - Tg0)
 ! compute initial change in T
 dT  = -f0/dh
 !write(*,'(a,1x,f12.5,1x,10(e20.10,1x))') 'Tg1, f0, f1, dh, dT = ', Tg1, f0, f1, dh, dT
 ! exit if already close enough
 if(abs(dT)<atol)then
  Tk = Tg0+dT
  return
 end if

 ! **** iterate a little
 do iter=1,niter
  ! comute new value of Tg
  Tg1 = Tg0+dT
  ! get new function evaluation
  Ht1 = temp2ethpy(Tg1,1._dp,fc_param)
  f1  = Ht1 - E_spec
  ! compute derivative if dT
  dh  = (f1 - f0)/dT
  ! compute change in T
  dT  = -f1/dh
  ! print progress
  !write(*,'(a,1x,i4,1x,f12.5,1x,10(e20.10,1x))') 'iter, Tg1, Ht1, f1, dh, dT = ', iter, Tg1, Ht1, f1, dh, dT
  ! exit if converged
  if(abs(dT)<atol)then
   Tk = Tg1+dT
   return
  end if
  ! get ready for next iteration -- save old function evaluation and temperature
  f0  = f1
  Tg0 = Tg1
  ! and check for convergence
  if(iter==niter)then; err=20; message=trim(message)//"failedToConverge"; return; end if
 end do  ! (iteration loop)
 end subroutine E2T_nosoil


 ! ************************************************************************************************************************
 ! public function temp2ethpy: compute total enthalpy based on temperature and mass (J m-3)
 ! ************************************************************************************************************************
 function temp2ethpy(Tk,BulkDenWater,fc_param)
 ! used to compute enthalpy based on temperature and total mass in layer (snow or soil)
 ! NOTE: enthalpy is a relative value, defined as zero at Tfreeze where all water is liquid
 USE multiconst, only: Tfreeze, &                   ! freezing point of water (K)
                       Cp_soil,Cp_water,Cp_ice,&    ! specific heat of soil, water and ice (J kg-1 K-1)
                       LH_fus                       ! latent heat of fusion (J kg-1)
 implicit none
 ! declare dummy variables
 real(dp),intent(in)  :: Tk            ! layer temperature (K)
 real(dp),intent(in)  :: BulkDenWater  ! bulk density of water (kg m-3)
 real(dp),intent(in)  :: fc_param      ! freezing curve parameter (K-1)
 real(dp)             :: temp2ethpy    ! return value of the function, total specific enthalpy (J m-3)
 ! declare local variables
 real(dp)             :: frac_liq      ! fraction of liquid water
 real(dp)             :: enthTempWater ! temperature component of specific enthalpy for total water (liquid and ice) (J kg-1)
 real(dp)             :: enthMass      ! mass component of specific enthalpy (J kg-1)
 ! NOTE: this function assumes the freezing curve for snow ... it needs modification to use vanGenuchten functions for soil
 ! compute the fraction of liquid water in the given layer
 frac_liq     = 1._dp / ( 1._dp + ( fc_param*( Tfreeze - min(Tk,Tfreeze) ) )**2._dp )
 ! compute the temperature component of enthalpy for the soil constituent (J kg-1)
 !enthTempSoil = Cp_soil*(Tk - Tfreeze)
 ! compute the temperature component of enthalpy for total water (J kg-1)
 ! NOTE: negative enthalpy means require energy to bring to Tfreeze
 if(Tk< Tfreeze) enthTempWater =   Cp_ice*(Tk - Tfreeze) - (Cp_water - Cp_ice)*(atan(fc_param*(Tfreeze - Tk))/fc_param)
 if(Tk>=Tfreeze) enthTempWater = Cp_water*(Tk - Tfreeze)
 ! compute the mass component of enthalpy -- energy required to melt ice (J kg-1)
 ! NOTE: negative enthalpy means require energy to bring to Tfreeze
 enthMass     = -LH_fus*(1._dp - frac_liq)
 ! finally, compute the total enthalpy (J m-3)
 ! NOTE: this is the case for snow (no soil).. function needs modification to use vanGenuchten functions for soil
 temp2ethpy   = BulkDenWater*(enthTempWater + enthMass) !+ BulkDenSoil*enthTempSoil
 end function temp2ethpy


end module ConvE2Temp_module
