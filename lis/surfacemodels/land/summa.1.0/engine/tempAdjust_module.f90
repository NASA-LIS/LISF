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

module tempAdjust_module
! data types
USE nrtype
! physical constants
USE multiconst,only:Tfreeze         ! freezing point of pure water (K)
USE multiconst,only:LH_fus          ! latent heat of fusion (J kg-1)
USE multiconst,only:Cp_ice          ! specific heat of ice (J kg-1 K-1)
USE multiconst,only:Cp_water        ! specific heat of liquid water (J kg-1 K-1)
USE multiconst,only:iden_water      ! intrinsic density of water (kg m-3)
implicit none
private
public::tempAdjust

contains


 ! ************************************************************************************************
 ! public subroutine tempAdjust: compute change in snow stored on the vegetation canopy
 ! ************************************************************************************************
 subroutine tempAdjust(&
                       ! input: derived parameters
                       canopyDepth,                 & ! intent(in): canopy depth (m)
                       ! input/output: data structures
                       mpar_data,                   & ! intent(in):    model parameters
                       prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,                   & ! intent(out):   model diagnostic variables for a local HRU
                       ! output: error control
                       err,message)                   ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_types,only:&
                     var_d,              & ! data vector (dp)
                     var_dlength           ! data vector with variable length dimension (dp)
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookPARAM,iLookPROG,iLookDIAG  ! named variables for structure elements
 ! utility routines
 USE snow_utils_module,only:fracliquid     ! compute fraction of liquid water
 USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
 implicit none
 ! ------------------------------------------------------------------------------------------------
 ! input: derived parameters
 real(dp),intent(in)             :: canopyDepth         ! depth of the vegetation canopy (m)
 ! input/output: data structures
 type(var_dlength),intent(in)    :: mpar_data           ! model parameters
 type(var_dlength),intent(inout) :: prog_data           ! model prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data           ! model diagnostic variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! ------------------------------------------------------------------------------------------------
 ! local variables for canopy thermodynamics
 integer(i4b)                  :: iTry                       ! trial index
 integer(i4b)                  :: iter                       ! iteration index
 integer(i4b),parameter        :: maxiter=100                ! maximum number of iterations
 real(dp)                      :: fLiq                       ! fraction of liquid water (-)
 real(dp)                      :: tempMin,tempMax            ! solution constraints for temperature (K)
 real(dp)                      :: nrgMeltFreeze              ! energy required to melt-freeze the water to the current canopy temperature (J m-3)
 real(dp)                      :: scalarCanopyWat            ! total canopy water (kg m-2)
 real(dp)                      :: scalarCanopyIceOld         ! canopy ice content after melt-freeze to the initial temperature (kg m-2)
 real(dp),parameter            :: resNrgToler=0.1_dp         ! tolerance for the energy residual (J m-3)
 real(dp)                      :: f1,f2,x1,x2,fTry,xTry,fDer,xInc ! iteration variables
 logical(lgt)                  :: fBis                       ! .true. if bisection
 ! -------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='tempAdjust/'
 ! ------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&

 ! model parameters for canopy thermodynamics (input)
 snowfrz_scale             => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),              & ! intent(in): [dp] scaling factor for snow freezing curve (K)
 specificHeatVeg           => mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1),            & ! intent(in): [dp] specific heat of vegetation mass (J kg-1 K-1)
 maxMassVegetation         => mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1),          & ! intent(in): [dp] maximum mass of vegetation (full foliage) (kg m-2)

 ! state variables (input/output)
 scalarCanopyLiq           => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),             & ! intent(inout): [dp] mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyIce           => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1),             & ! intent(inout): [dp] mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyTemp          => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1),            & ! intent(inout): [dp] temperature of the vegetation canopy (K)

 ! diagnostic variables (output)
 scalarBulkVolHeatCapVeg   => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1)      & ! intent(out): [dp] volumetric heat capacity of the vegetation (J m-3 K-1)

 )  ! associate variables in the data structures
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! ** preliminaries

 ! compute the total canopy water (state variable: will not change)
 scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce
 !write(*,'(a,1x,3(f20.10,1x))') 'scalarCanopyWat, scalarCanopyLiq, scalarCanopyIce = ', scalarCanopyWat, scalarCanopyLiq, scalarCanopyIce

 ! compute the fraction of liquid water associated with the canopy temperature
 fLiq = fracliquid(scalarCanopyTemp,snowfrz_scale)

 ! compute the new volumetric ice content
 ! NOTE: new value; iterations will adjust this value for consistency with temperature
 scalarCanopyIceOld = (1._dp - fLiq)*scalarCanopyWat

 ! compute volumetric heat capacity of vegetation (J m-3 K-1)
 scalarBulkVolHeatCapVeg = specificHeatVeg*maxMassVegetation/canopyDepth + & ! vegetation component
                           Cp_water*scalarCanopyLiq/canopyDepth          + & ! liquid water component
                           Cp_ice*scalarCanopyIce/canopyDepth                ! ice component

 ! compute the energy required to melt-freeze the water to the current canopy temperature (J m-3)
 nrgMeltFreeze = LH_fus*(scalarCanopyIceOld - scalarCanopyIce)/canopyDepth

 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! ** get ready for iterating

 ! compute initial function and derivative
 x1   = scalarCanopyTemp
 f1   = nrgMeltFreeze
 fDer = resNrgDer(x1,scalarBulkVolHeatCapVeg,snowfrz_scale)

 ! compute new function based on newton step from the first function
 x2 = x1 + f1 / fDer
 f2 = resNrgFunc(x2,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 !print*, 'x1, x2 = ', x1, x2
 !print*, 'f1, f2 = ', f1, f2

 ! ensure that we bracket the root
 if(f1*f2 > 0._dp)then
  xInc = f1 / fDer
  x2   = 1._dp
  do iter=1,maxiter
   ! successively expand limit in order to bracket the root
   x2 = x1 + sign(x2,xInc)*2._dp
   f2 = resNrgFunc(x2,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
   if(f1*f2 < 0._dp)exit
   ! check that we bracketed the root
   ! (should get here in just a couple of expansions)
   if(iter==maxiter)then
    message=trim(message)//'unable to bracket the root'
    err=20; return
   end if
  end do ! trying to bracket the root
 end if  ! first check that we bracketed the root
 !print*, 'x1, x2 = ', x1, x2
 !print*, 'f1, f2 = ', f1, f2

 ! define initial constraints
 if(x1 < x2)then
  tempMin = x1
  tempMax = x2
 else
  tempMin = x2
  tempMax = x1
 end if
 !print*, 'tempMin, tempMax = ', tempMin, tempMax

 ! get starting trial
 xInc = huge(1._dp)
 xTry = 0.5_dp*(x1 + x2)
 fTry = resNrgFunc(xTry,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 fDer = resNrgDer(xTry,scalarBulkVolHeatCapVeg,snowfrz_scale)
 !print*, 'xTry = ', xTry
 !print*, 'fTry = ', fTry

 ! check the functions at the limits (should be of opposing sign)
 !f1 = resNrgFunc(tempMax,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 !f2 = resNrgFunc(tempMin,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
 !print*, 'f1, f2 = ', f1, f2

 ! -----------------------------------------------------------------------------------------------------------------------------------------------------
 ! iterate
 do iter=1,maxiter

  ! bisect if out of range
  if(xTry <= tempMin .or. xTry >= tempMax)then
   xTry = 0.5_dp*(tempMin + tempMax)  ! new value
   fBis = .true.

  ! value in range; use the newton step
  else
   xInc = fTry/fDer
   xTry = xTry + xInc
   fBis = .false.

  end if  ! (switch between bi-section and newton)

  ! compute new function and derivative
  fTry = resNrgFunc(xTry,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
  fDer = resNrgDer(xTry,scalarBulkVolHeatCapVeg,snowfrz_scale)
  !print*, 'tempMin, tempMax = ', tempMin, tempMax

  ! update limits
  if(fTry < 0._dp)then
   tempMax = min(xTry,tempMax)
  else
   tempMin = max(tempMin,xTry)
  end if

  ! check the functions at the limits (should be of opposing sign)
  !f1 = resNrgFunc(tempMax,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
  !f2 = resNrgFunc(tempMin,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
  !print*, 'f1, f2 = ', f1, f2

  ! print progress
  !write(*,'(a,1x,i4,1x,l1,1x,e20.10,1x,4(f20.10,1x))') 'iter, fBis, fTry, xTry, xInc, tempMin, tempMax = ', iter, fBis, fTry, xTry, xInc, tempMin, tempMax

  ! check convergence
  if(abs(fTry) < resNrgToler) exit

  ! check non-convergence
  if(iter==maxiter)then
   ! (print out a 1-d x-section)
   do iTry=1,maxiter
    xTry = 1.0_dp*real(iTry,kind(1._dp))/real(maxiter,kind(1._dp)) + 272.5_dp
    fTry = resNrgFunc(xTry,scalarCanopyTemp,scalarBulkVolHeatCapVeg,snowfrz_scale)
    write(*,'(a,1x,i4,1x,e20.10,1x,4(f20.10,1x))') 'iTry, fTry, xTry = ', iTry, fTry, xTry
   end do
   ! (return with error)
   message=trim(message)//'unable to converge'
   err=20; return
  end if

 end do  ! iterating
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------

 ! update state variables
 scalarCanopyTemp = xTry
 scalarCanopyIce  = (1._dp - fracliquid(xTry,snowfrz_scale))*scalarCanopyWat
 scalarCanopyLiq  = scalarCanopyWat - scalarCanopyIce

 ! end association to variables in the data structure
 end associate

 contains


  ! ************************************************************************************************
  ! internal function resNrgFunc: calculate the residual in energy (J m-3)
  ! ************************************************************************************************
  function resNrgFunc(xTemp,xTemp0,bulkVolHeatCapVeg,snowfrz_scale)
  !
  implicit none
  real(dp),intent(in) :: xTemp              ! temperature (K)
  real(dp),intent(in) :: xTemp0             ! initial temperature (K)
  real(dp),intent(in) :: bulkVolHeatCapVeg  ! volumetric heat capacity of veg (J m-3 K-1)
  real(dp),intent(in) :: snowfrz_scale      ! scaling factor in freezing curve (K-1)
  real(dp)            :: xIce               ! canopy ice content (kg m-2)
  real(dp)            :: resNrgFunc         ! residual in energy (J m-3)
  xIce       = (1._dp - fracliquid(xTemp,snowfrz_scale))*scalarCanopyWat
  resNrgFunc = -bulkVolHeatCapVeg*(xTemp - xTemp0) + LH_fus*(xIce - scalarCanopyIceOld)/canopyDepth + nrgMeltFreeze
  return
  end function resNrgFunc


  ! ************************************************************************************************
  ! internal function resNrgDer: calculate the derivative (J m-3 K-1)
  ! ************************************************************************************************
  function resNrgDer(xTemp,bulkVolHeatCapVeg,snowfrz_scale)
  implicit none
  real(dp),intent(in) :: xTemp              ! temperature (K)
  real(dp),intent(in) :: bulkVolHeatCapVeg  ! volumetric heat capacity of veg (J m-3 K-1)
  real(dp),intent(in) :: snowfrz_scale      ! scaling factor in freezing curve (K-1)
  real(dp)            :: dW_dT              ! derivative in canopy ice content w.r.t. temperature (kg m-2 K-1)
  real(dp)            :: resNrgDer          ! derivative (J m-3 K-1)
  dW_dT     = -scalarCanopyWat*dFracLiq_dTk(xTemp,snowfrz_scale)
  resNrgDer = bulkVolHeatCapVeg - dW_dT*LH_fus/canopyDepth
  return
  end function resNrgDer

 end subroutine tempAdjust


end module tempAdjust_module
