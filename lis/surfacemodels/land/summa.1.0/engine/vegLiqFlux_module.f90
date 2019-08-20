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

module vegLiqFlux_module
USE nrtype
! look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:         &
                      unDefined,    & ! original model (no flexibility in canopy interception): 100% of rainfall is intercepted by the vegetation canopy
                      sparseCanopy, & ! fraction of rainfall that never hits the canopy (throughfall); drainage above threshold
                      storageFunc     ! throughfall a function of canopy storage; 100% throughfall when canopy is at capacity
implicit none
private
public::vegLiqFlux
contains


 ! ************************************************************************************************
 ! public subroutine vegLiqFlux: compute water balance for the vegetation canopy
 ! ************************************************************************************************
 subroutine vegLiqFlux(&
                       ! input
                       computeVegFlux,               & ! intent(in): flag to denote if computing energy flux over vegetation
                       scalarCanopyLiqTrial,         & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       scalarRainfall,               & ! intent(in): rainfall rate (kg m-2 s-1)
                       ! input-output: data structures
                       mpar_data,                    & ! intent(in): model parameters
                       diag_data,                    & ! intent(in): local HRU model diagnostic variables
                       ! output
                       scalarThroughfallRain,        & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                       scalarCanopyLiqDrainage,      & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                       scalarThroughfallRainDeriv,   & ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
                       scalarCanopyLiqDrainageDeriv, & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                       err,message)                    ! intent(out): error control
 ! model decisions
 USE globalData,only:model_decisions                              ! model decision structure
 USE var_lookup,only:iLookDECISIONS                               ! named variables for elements of the decision structure
 ! named variables 
 USE var_lookup,only:iLookPARAM,iLookDIAG ! named variables for structure elements
 ! data types
 USE data_types,only:var_d           ! x%var(:)       (dp)
 USE data_types,only:var_dlength     ! x%var(:)%dat   (dp)
 implicit none
 ! input
 logical(lgt),intent(in)         :: computeVegFlux               ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(dp),intent(in)             :: scalarCanopyLiqTrial         ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)             :: scalarRainfall               ! rainfall (kg m-2 s-1)
 ! input-output: data structures
 type(var_dlength),intent(in)    :: mpar_data                    ! model parameters
 type(var_dlength),intent(inout) :: diag_data                    ! model diagnostic variables for the local basin
 ! output
 real(dp),intent(out)            :: scalarThroughfallRain        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)            :: scalarCanopyLiqDrainage      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp),intent(out)            :: scalarThroughfallRainDeriv   ! derivative in throughfall w.r.t. canopy liquid water (s-1)
 real(dp),intent(out)            :: scalarCanopyLiqDrainageDeriv ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
 integer(i4b),intent(out)        :: err                          ! error code
 character(*),intent(out)        :: message                      ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! make association of local variables with information in the data structures
 associate(&
  ixCanopyInterception       => model_decisions(iLookDECISIONS%cIntercept)%iDecision, & ! intent(in): index defining choice of parameterization for canopy interception
  scalarCanopyLiqMax         => diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1),   & ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
  scalarThroughfallScaleRain => mpar_data%var(iLookPARAM%throughfallScaleRain)%dat(1),& ! intent(in): fraction of rain that hits the ground without touching the canopy (-)
  scalarCanopyDrainageCoeff  => mpar_data%var(iLookPARAM%canopyDrainageCoeff)%dat(1)  & ! intent(in): canopy drainage coefficient (s-1)
 ) ! associating local variables with information in the data structures 
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="vegLiqFlux/"

 ! set throughfall to inputs if vegetation is completely buried with snow
 if(.not.computeVegFlux)then
  scalarThroughfallRain        = scalarRainfall
  scalarCanopyLiqDrainage      = 0._dp
  scalarThroughfallRainDeriv   = 0._dp
  scalarCanopyLiqDrainageDeriv = 0._dp
  return
 end if

 ! compute throughfall
 select case(ixCanopyInterception)

  ! original model (no flexibility in canopy interception): 100% of rainfall is intercepted by the vegetation canopy
  ! NOTE: this could be done with scalarThroughfallScaleRain=0, though requires setting scalarThroughfallScaleRain in all test cases
  case(unDefined)
   scalarThroughfallRain      = 0._dp
   scalarThroughfallRainDeriv = 0._dp

  ! fraction of rainfall hits the ground without ever touching the canopy
  case(sparseCanopy)
   scalarThroughfallRain      = scalarThroughfallScaleRain*scalarRainfall
   scalarThroughfallRainDeriv = 0._dp

  ! throughfall a function of canopy storage
  case(storageFunc)

   ! throughfall during wetting-up phase
   if(scalarCanopyLiqTrial < scalarCanopyLiqMax)then
    scalarThroughfallRain      = scalarRainfall*(scalarCanopyLiqTrial/scalarCanopyLiqMax)
    scalarThroughfallRainDeriv = scalarRainfall/scalarCanopyLiqMax

   ! all rain falls through the canopy when the canopy is at capacity
   else
    scalarThroughfallRain      = scalarRainfall
    scalarThroughfallRainDeriv = 0._dp
   end if

  case default; err=20; message=trim(message)//'unable to identify option for canopy interception'; return

 end select ! (option for canopy interception)

 ! compute canopy drainage
 if(scalarCanopyLiqTrial > scalarCanopyLiqMax)then
  scalarCanopyLiqDrainage       = scalarCanopyDrainageCoeff*(scalarCanopyLiqTrial - scalarCanopyLiqMax)
  scalarCanopyLiqDrainageDeriv  = scalarCanopyDrainageCoeff
 else
  scalarCanopyLiqDrainage       = 0._dp
  scalarCanopyLiqDrainageDeriv  = 0._dp
 end if

 !write(*,'(a,1x,f25.15)') 'scalarCanopyLiqDrainage = ', scalarCanopyLiqDrainage

 ! end association of local variables with information in the data structures
 end associate

 end subroutine vegLiqFlux


end module vegLiqFlux_module
