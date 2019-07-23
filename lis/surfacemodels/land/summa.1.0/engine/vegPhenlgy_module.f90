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

module vegPhenlgy_module
! data types
USE nrtype
! look-up values for the boundary conditions
USE mDecisions_module,only:      &
 prescribedHead,                  &         ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 prescribedTemp,                  &         ! prescribed temperature
 zeroFlux                                   ! zero flux
! look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:      &
 noah_mp,                        &         ! full Noah-MP implementation (including albedo)
 CLM_2stream,                    &         ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,                    &         ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,                     &         ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                                  ! Beer's Law (as implemented in VIC)
implicit none
private
public::vegPhenlgy
! algorithmic parameters
real(dp),parameter     :: valueMissing=-9999._dp  ! missing value, used when diagnostic or state variables are undefined
real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers
contains


 ! ************************************************************************************************
 ! public subroutine vegPhenlgy: compute vegetation phenology
 ! ************************************************************************************************
 subroutine vegPhenlgy(&
                       ! input/output: data structures
                       model_decisions,             & ! intent(in):    model decisions
                       type_data,                   & ! intent(in):    type of vegetation and soil
                       attr_data,                   & ! intent(in):    spatial attributes
                       mpar_data,                   & ! intent(in):    model parameters
                       prog_data,                   & ! intent(in):    prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): diagnostic variables for a local HRU
                       ! output
                       computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                       canopyDepth,                 & ! intent(out): canopy depth (m)
                       exposedVAI,                  & ! intent(out): exposed vegetation area index (LAI + SAI)
                       err,message)                   ! intent(out): error control
 ! -------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_types,only:&
                     var_i,            & ! data vector (i4b)
                     var_d,            & ! data vector (dp)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookTYPE,iLookATTR,iLookPARAM,iLookDIAG,iLookPROG  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                                      ! named variables for elements of the decision structure
 ! modules
 USE NOAHMP_ROUTINES,only:phenology         ! determine vegetation phenology
 ! common variables
 USE globalData,only:urbanVegCategory       ! vegetation category for urban areas
 USE globalData,only:fracJulday             ! fractional julian days since the start of year
 USE globalData,only:yearLength             ! number of days in the current year
 implicit none
 ! -------------------------------------------------------------------------------------------------
 ! input/output
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 type(var_i),intent(in)          :: type_data           ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data           ! spatial attributes
 type(var_dlength),intent(in)    :: mpar_data           ! model parameters
 type(var_dlength),intent(in)    :: prog_data           ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data           ! diagnostic variables for a local HRU
 ! output
 logical(lgt),intent(out)        :: computeVegFlux      ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(dp),intent(out)            :: canopyDepth         ! canopy depth (m)
 real(dp),intent(out)            :: exposedVAI          ! exposed vegetation area index (LAI + SAI)
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! -------------------------------------------------------------------------------------------------
 ! local
 real(dp)                 :: notUsed_heightCanopyTop    ! height of the top of the canopy layer (m)
 real(dp)                 :: heightAboveSnow            ! height top of canopy is above the snow surface (m)
 ! initialize error control
 err=0; message="vegPhenlgy/"
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&

 ! input: model decisions
 ix_bcUpprTdyn                   => model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision,      & ! intent(in): [i4b] choice of upper boundary condition for thermodynamics
 ix_bcUpprSoiH                   => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,      & ! intent(in): [i4b] index of method used for the upper boundary condition for soil hydrology

 ! local attributes
 vegTypeIndex                    => type_data%var(iLookTYPE%vegTypeIndex),                     & ! intent(in): [i4b] vegetation type index
 latitude                        => attr_data%var(iLookATTR%latitude),                         & ! intent(in): [dp] latitude

 ! model state variables
 scalarSnowDepth                 => prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),           & ! intent(in):    [dp] snow depth on the ground surface (m)
 scalarCanopyTemp                => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1),          & ! intent(in):    [dp] temperature of the vegetation canopy at the start of the sub-step (K)

 ! diagnostic variables and parameters (input)
 heightCanopyTop                 => mpar_data%var(iLookPARAM%heightCanopyTop)%dat(1),          & ! intent(in):    [dp] height of the top of the canopy layer (m)
 heightCanopyBottom              => mpar_data%var(iLookPARAM%heightCanopyBottom)%dat(1),       & ! intent(in):    [dp] height of the bottom of the canopy layer (m)
 scalarRootZoneTemp              => diag_data%var(iLookDIAG%scalarRootZoneTemp)%dat(1),        & ! intent(in):    [dp] root zone temperature (K)

 ! diagnostic variables and parameters (input/output)
 scalarLAI                       => diag_data%var(iLookDIAG%scalarLAI)%dat(1),                 & ! intent(inout): [dp] one-sided leaf area index (m2 m-2)
 scalarSAI                       => diag_data%var(iLookDIAG%scalarSAI)%dat(1),                 & ! intent(inout): [dp] one-sided stem area index (m2 m-2)

 ! diagnostic variables and parameters (output)
 scalarExposedLAI                => diag_data%var(iLookDIAG%scalarExposedLAI)%dat(1),          & ! intent(out): [dp] exposed leaf area index after burial by snow (m2 m-2)
 scalarExposedSAI                => diag_data%var(iLookDIAG%scalarExposedSAI)%dat(1),          & ! intent(out): [dp] exposed stem area index after burial by snow (m2 m-2)
 scalarGrowingSeasonIndex        => diag_data%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1)   & ! intent(out): [dp] growing season index (0=off, 1=on)

 ) ! associate variables in data structure
 ! ----------------------------------------------------------------------------------------------------------------------------------

 ! check if we have isolated the snow-soil domain (used in test cases)
 if(ix_bcUpprTdyn == prescribedTemp .or. ix_bcUpprTdyn == zeroFlux .or. ix_bcUpprSoiH == prescribedHead)then

  ! isolated snow-soil domain: do not compute fluxes over vegetation
  computeVegFlux = .false.

  ! set vegetation phenology variables to missing
  scalarLAI                = valueMissing    ! one-sided leaf area index (m2 m-2)
  scalarSAI                = valueMissing    ! one-sided stem area index (m2 m-2)
  scalarExposedLAI         = valueMissing    ! exposed leaf area index after burial by snow (m2 m-2)
  scalarExposedSAI         = valueMissing    ! exposed stem area index after burial by snow (m2 m-2)
  scalarGrowingSeasonIndex = valueMissing    ! growing season index (0=off, 1=on)
  exposedVAI               = valueMissing    ! exposed vegetation area index (m2 m-2)
  canopyDepth              = valueMissing    ! canopy depth (m)
  heightAboveSnow          = valueMissing    ! height top of canopy is above the snow surface (m)

 ! compute vegetation phenology (checks for complete burial of vegetation)
 else

  ! determine vegetation phenology
  ! NOTE: recomputing phenology every sub-step accounts for changes in exposed vegetation associated with changes in snow depth
  call phenology(&
                 ! input
                 vegTypeIndex,                & ! intent(in): vegetation type index
                 urbanVegCategory,            & ! intent(in): vegetation category for urban areas
                 scalarSnowDepth,             & ! intent(in): snow depth on the ground surface (m)
                 scalarCanopyTemp,            & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                 latitude,                    & ! intent(in): latitude
                 yearLength,                  & ! intent(in): number of days in the current year
                 fracJulday,                  & ! intent(in): fractional julian days since the start of year
                 scalarLAI,                   & ! intent(inout): one-sided leaf area index (m2 m-2)
                 scalarSAI,                   & ! intent(inout): one-sided stem area index (m2 m-2)
                 scalarRootZoneTemp,          & ! intent(in): root zone temperature (K)
                 ! output
                 notUsed_heightCanopyTop,     & ! intent(out): height of the top of the canopy layer (m)
                 scalarExposedLAI,            & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                 scalarExposedSAI,            & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                 scalarGrowingSeasonIndex     ) ! intent(out): growing season index (0=off, 1=on)

  ! determine additional phenological variables
  exposedVAI      = scalarExposedLAI + scalarExposedSAI   ! exposed vegetation area index (m2 m-2)
  canopyDepth     = heightCanopyTop - heightCanopyBottom  ! canopy depth (m)
  heightAboveSnow = heightCanopyTop - scalarSnowDepth     ! height top of canopy is above the snow surface (m)

  ! determine if need to include vegetation in the energy flux routines
  computeVegFlux  = (exposedVAI > 0.05_dp .and. heightAboveSnow > 0.05_dp)

 end if  ! (check if the snow-soil column is isolated)

 ! end association to variables in the data structure
 end associate

 end subroutine vegPhenlgy


end module vegPhenlgy_module
