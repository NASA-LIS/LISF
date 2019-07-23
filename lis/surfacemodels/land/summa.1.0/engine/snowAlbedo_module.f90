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

module snowAlbedo_module

! data types
USE nrtype                        ! numerical recipes data types

! physical constants
USE multiconst,only:Tfreeze       ! freezing point of pure water (K)

! look-up values for the choice of snow albedo options
USE mDecisions_module,only:  &
 constantDecay,              &    ! constant decay in snow albedo (e.g., VIC, CLASS)
 variableDecay                    ! variable decay in snow albedo (e.g., BATS approach, with destructive metamorphism + soot content)

! look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:  &
 noah_mp,                    &    ! full Noah-MP implementation (including albedo)
 CLM_2stream,                &    ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,                &    ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,                 &    ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                         ! Beer's Law (as implemented in VIC)

! -------------------------------------------------------------------------------------------------
implicit none
private
public::snowAlbedo
! dimensions
integer(i4b),parameter        :: nBands=2      ! number of spectral bands for shortwave radiation
contains


 ! *******************************************************************************************************
 ! public subroutine snowAlbedo: muster program to compute energy fluxes at vegetation and ground surfaces
 ! *******************************************************************************************************
 subroutine snowAlbedo(&
                       ! input: model control
                       dt,                                    & ! intent(in): model time step (s)
                       snowPresence,                          & ! intent(in): logical flag to denote if snow is present
                       ! input/output: data structures
                       model_decisions,                       & ! intent(in):    model decisions
                       mpar_data,                             & ! intent(in):    model parameters
                       flux_data,                             & ! intent(in):    model flux variables
                       diag_data,                             & ! intent(inout): model diagnostic variables for a local HRU
                       prog_data,                             & ! intent(inout): model prognostic variables for a local HRU
                       ! output: error control
                       err,message)                             ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_types,only:&
                     var_i,            & ! data vector (i4b)
                     var_d,            & ! data vector (dp)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookPARAM,iLookFLUX,iLookDIAG,iLookPROG   ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                             ! named variables for elements of the decision structure
 ! provide access to desired modules
 USE snow_utils_module,only:fracliquid                          ! compute fraction of liquid water at a given temperature
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)             :: dt                          ! model time step
 logical(lgt),intent(in)         :: snowPresence                ! logical flag to denote if snow is present
 ! input/output: data structures
 type(model_options),intent(in)  :: model_decisions(:)          ! model decisions
 type(var_dlength),intent(in)    :: mpar_data                   ! model parameters
 type(var_dlength),intent(in)    :: flux_data                   ! model flux variables
 type(var_dlength),intent(inout) :: diag_data                   ! model diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: prog_data                   ! model prognostic variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                         ! error code
 character(*),intent(out)        :: message                     ! error message
 ! local variables
 integer(i4b),parameter          :: ixVisible=1                  ! named variable to define index in array of visible part of the spectrum
 integer(i4b),parameter          :: ixNearIR=2                   ! named variable to define index in array of near IR part of the spectrum
 real(dp),parameter              :: valueMissing=-9999._dp       ! missing value -- will cause problems if snow albedo is ever used for the non-snow case
 real(dp),parameter              :: slushExp=10._dp              ! "slush" exponent, to increase decay when snow is near Tfreeze
 real(dp),parameter              :: fractionLiqThresh=0.001_dp   ! threshold for the fraction of liquid water to switch to spring albedo minimum
 real(dp)                        :: fractionLiq                  ! fraction of liquid water (-)
 real(dp)                        :: age1,age2,age3               ! aging factors (-)
 real(dp)                        :: decayFactor                  ! albedo decay factor (-)
 real(dp)                        :: refreshFactor                ! albedo refreshment factor, representing albedo increase due to snowfall (-)
 real(dp)                        :: albedoMin                    ! minimum albedo -- depends if in winter or spring conditions (-)
 real(dp)                        :: fZen                         ! factor to modify albedo at low zenith angles (-)
 real(dp),parameter              :: bPar=2._dp                   ! empirical parameter in fZen
 ! initialize error control
 err=0; message='snowAlbedo/'
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&
 ! input: model decisions
 ixCanopySrad              => model_decisions(iLookDECISIONS%canopySrad)%iDecision,   & ! intent(in): index of method used for canopy sw radiation
 ixAlbedoMethod            => model_decisions(iLookDECISIONS%alb_method)%iDecision,   & ! intent(in): index of method used for snow albedo
 ! input: model parameters
 Frad_vis                  => mpar_data%var(iLookPARAM%Frad_vis)%dat(1),              & ! intent(in): fraction of radiation in visible part of spectrum (-)
 Frad_direct               => mpar_data%var(iLookPARAM%Frad_direct)%dat(1),           & ! intent(in): fraction direct solar radiation (-)
 albedoMax                 => mpar_data%var(iLookPARAM%albedoMax)%dat(1),             & ! intent(in): maximum snow albedo for a single spectral band (-)
 albedoMinWinter           => mpar_data%var(iLookPARAM%albedoMinWinter)%dat(1),       & ! intent(in): minimum snow albedo during winter for a single spectral band (-)
 albedoMinSpring           => mpar_data%var(iLookPARAM%albedoMinSpring)%dat(1),       & ! intent(in): minimum snow albedo during spring for a single spectral band (-)
 albedoMaxVisible          => mpar_data%var(iLookPARAM%albedoMaxVisible)%dat(1),      & ! intent(in): maximum snow albedo in the visible part of the spectrum (-)
 albedoMinVisible          => mpar_data%var(iLookPARAM%albedoMinVisible)%dat(1),      & ! intent(in): minimum snow albedo in the visible part of the spectrum (-)
 albedoMaxNearIR           => mpar_data%var(iLookPARAM%albedoMaxNearIR)%dat(1),       & ! intent(in): maximum snow albedo in the near infra-red part of the spectrum (-)
 albedoMinNearIR           => mpar_data%var(iLookPARAM%albedoMinNearIR)%dat(1),       & ! intent(in): minimum snow albedo in the near infra-red part of the spectrum (-)
 albedoDecayRate           => mpar_data%var(iLookPARAM%albedoDecayRate)%dat(1),       & ! intent(in): albedo decay rate (s)
 tempScalGrowth            => mpar_data%var(iLookPARAM%tempScalGrowth)%dat(1),        & ! intent(in): temperature scaling factor for grain growth (K-1)
 albedoSootLoad            => mpar_data%var(iLookPARAM%albedoSootLoad)%dat(1),        & ! intent(in): soot load factor (-)
 albedoRefresh             => mpar_data%var(iLookPARAM%albedoRefresh)%dat(1),         & ! intent(in): critical mass necessary for albedo refreshment (kg m-2)
 snowfrz_scale             => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1),         & ! intent(in): scaling parameter for the freezing curve for snow (K-1)
 ! input: model variables
 surfaceTemp               => prog_data%var(iLookPROG%mLayerTemp)%dat(1),             & ! intent(in): surface temperature
 snowfallRate              => flux_data%var(iLookFLUX%scalarSnowfall)%dat(1),         & ! intent(in): snowfall rate (kg m-2 s-1)
 cosZenith                 => diag_data%var(iLookDIAG%scalarCosZenith)%dat(1),        & ! intent(in): cosine of the zenith angle (-)
 ! input/output: model variables
 scalarSnowAlbedo          => prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1),       & ! intent(inout): snow albedo for the entire spectral band (-)
 spectralSnowAlbedoDirect  => diag_data%var(iLookDIAG%spectralSnowAlbedoDirect)%dat,  & ! intent(inout): direct snow albedo in each spectral band (-)
 spectralSnowAlbedoDiffuse => prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat  & ! intent(inout): diffuse snow albedo in each spectral band (-)
 ) ! end associate statement
 ! --------------------------------------------------------------------------------------------------------------------------------------

 ! return early if computing radiation in noah-MP
 if(ixCanopySrad==noah_mp) return

 ! return early if no snow
 if(.not. snowPresence)then
  scalarSnowAlbedo             = valueMissing
  spectralSnowAlbedoDirect(:)  = valueMissing
  spectralSnowAlbedoDiffuse(:) = valueMissing
  return
 end if

 ! compute fractional increase in albedo associated with snowfall
 refreshFactor = dt*snowfallRate/albedoRefresh

 ! identify option for snow albedo
 select case(ixAlbedoMethod)


  ! *** constant decay rate
  case(constantDecay)
   ! compute decay rate
   decayFactor = dt/albedoDecayRate
   ! compute minimum albedo
   fractionLiq = fracliquid(surfaceTemp,snowfrz_scale) ! fraction of liquid water
   if(scalarSnowAlbedo < albedoMinWinter .or. fractionLiq > fractionLiqThresh)then
    albedoMin = albedoMinSpring
   else
    albedoMin = albedoMinWinter
   end if
   ! compute average albedo
   call computeAlbedo(scalarSnowAlbedo,refreshFactor,decayFactor,albedoMax,albedoMin)
   ! assume albedo is the same in visible and near infra-red bands, and for direct and diffuse radiation
   spectralSnowAlbedoDiffuse(ixVisible) = scalarSnowAlbedo
   spectralSnowAlbedoDiffuse(ixNearIR)  = scalarSnowAlbedo
   spectralSnowAlbedoDirect(ixVisible)  = scalarSnowAlbedo
   spectralSnowAlbedoDirect(ixNearIR)   = scalarSnowAlbedo


  ! *** variable decay rate
  case(variableDecay)
   ! compute decay factor
   age1 = exp(-tempScalGrowth*(Tfreeze - surfaceTemp ))  ! temperature dependence
   age2 = age1**slushExp                                 ! increase with liquid water
   age3 = albedoSootLoad                                 ! soot loading
   decayFactor = dt*(age1 + age2 + age3)/albedoDecayRate
   ! compute diffuse albedo for the different spectral bands
   call computeAlbedo(spectralSnowAlbedoDiffuse(ixVisible),refreshFactor,decayFactor,albedoMaxVisible,albedoMinVisible)
   call computeAlbedo(spectralSnowAlbedoDiffuse(ixNearIR), refreshFactor,decayFactor,albedoMaxNearIR, albedoMinNearIR)
   ! compute factor to modify direct albedo at low zenith angles
   if(cosZenith < 0.5_dp)then
    fZen = (1._dp/bPar)*( ((1._dp + bPar)/(1._dp + 2._dp*bPar*cosZenith)) - 1._dp)
   else
    fZen = 0._dp
   end if
   ! compute direct albedo
   spectralSnowAlbedoDirect(ixVisible) = spectralSnowAlbedoDiffuse(ixVisible) + 0.4_dp*fZen*(1._dp - spectralSnowAlbedoDiffuse(ixVisible))
   spectralSnowAlbedoDirect(ixNearIR)  = spectralSnowAlbedoDiffuse(ixNearIR)  + 0.4_dp*fZen*(1._dp - spectralSnowAlbedoDiffuse(ixNearIR))

   ! compute average albedo
   scalarSnowAlbedo = (        Frad_direct)*(Frad_vis*spectralSnowAlbedoDirect(ixVisible) + (1._dp - Frad_vis)*spectralSnowAlbedoDirect(ixNearIR) ) + &
                      (1._dp - Frad_direct)*(Frad_vis*spectralSnowAlbedoDirect(ixVisible) + (1._dp - Frad_vis)*spectralSnowAlbedoDirect(ixNearIR) )

  ! check that we identified the albedo option
  case default; err=20; message=trim(message)//'unable to identify option for snow albedo'; return

 end select  ! identify option for snow albedo

 ! check
 if(scalarSnowAlbedo < 0._dp)then; err=20; message=trim(message)//'unable to identify option for snow albedo'; return; end if

 ! end association to data structures
 end associate

 end subroutine snowAlbedo


 ! *******************************************************************************************************
 ! private subroutine computeAlbedo: compute change in albedo -- implicit solution
 ! *******************************************************************************************************
 subroutine computeAlbedo(snowAlbedo,refreshFactor,decayFactor,albedoMax,albedoMin)
 implicit none
 ! dummy variables
 real(dp),intent(inout)   :: snowAlbedo    ! snow albedo (-)
 real(dp),intent(in)      :: refreshFactor ! albedo refreshment factor (-)
 real(dp),intent(in)      :: decayFactor   ! albedo decay factor (-)
 real(dp),intent(in)      :: albedoMax     ! maximum albedo (-)
 real(dp),intent(in)      :: albedoMin     ! minimum albedo (-)
 ! local variables
 real(dp)                 :: albedoChange ! change in albedo over the time step (-)
 ! compute change in albedo
 albedoChange = refreshFactor*(albedoMax - snowAlbedo) - (decayFactor*(snowAlbedo - albedoMin)) / (1._dp + decayFactor)
 snowAlbedo   = snowAlbedo + albedoChange
 if(snowAlbedo > albedoMax) snowAlbedo = albedoMax
 end subroutine computeAlbedo


end module snowAlbedo_module
