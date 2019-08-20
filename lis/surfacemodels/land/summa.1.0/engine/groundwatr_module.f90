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

module groundwatr_module
! numerical recipes data types
USE nrtype
! model constants
USE multiconst,only:iden_water ! density of water (kg m-3)
! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization
! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_d,     & ! data vector (dp)
                    var_dlength  ! data vector with variable length dimension (dp)
! provide access to named variables defining elements in the data structures
USE var_lookup,only:iLookATTR,iLookPARAM,iLookPROG,iLookFLUX
 ! utility modules
implicit none
! constant parameters
real(dp),parameter     :: valueMissing=-9999._dp    ! missing value parameter
real(dp),parameter     :: verySmall=epsilon(1.0_dp) ! a very small number (used to avoid divide by zero)
real(dp),parameter     :: dx=1.e-8_dp               ! finite difference increment
private
public::groundwatr
contains


 ! ************************************************************************************************
 ! public subroutine groundwatr: compute the groundwater sink term in Richards' equation
 ! ************************************************************************************************
 !
 ! Method
 ! ------
 !
 ! Here we assume that water avaialble for shallow groundwater flow includes is all water above
 ! "field capacity" below the depth zCrit, where zCrit is defined as the lowest point in the soil
 ! profile where the volumetric liquid water content is less than field capacity.
 !
 ! We further assume that transmssivity (m2 s-1) for each layer is defined asuming that the water
 ! available for saturated flow is located at the bottom of the soil profile. Specifically:
 !  trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL
 !  trSoil(iLayer)  = trTotal(iLayer) - trTotal(iLayer+1)
 ! where zActive(iLayer) is the effective water table thickness for all layers up to and including
 ! the current layer (working from the bottom to the top).
 !
 ! The outflow from each layer is then (m3 s-1)
 !  mLayerOutflow(iLayer) = trSoil(iLayer)*tan_slope*contourLength
 ! where contourLength is the width of a hillslope (m) parallel to a stream
 !
 ! ************************************************************************************************
 subroutine groundwatr(&

                       ! input: model control
                       nSnow,                                  & ! intent(in): number of snow layers
                       nSoil,                                  & ! intent(in): number of soil layers
                       nLayers,                                & ! intent(in): total number of layers
                       getSatDepth,                            & ! intent(in): logical flag to compute index of the lowest saturated layer

                       ! input: state and diagnostic variables
                       mLayerdTheta_dPsi,                      & ! intent(in): derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
                       mLayerMatricHeadLiq,                    & ! intent(in): liquid water matric potential (m)
                       mLayerVolFracLiq,                       & ! intent(in): volumetric fraction of liquid water (-)
                       mLayerVolFracIce,                       & ! intent(in): volumetric fraction of ice (-)

                       ! input/output: data structures
                       attr_data,                              & ! intent(in):    spatial attributes
                       mpar_data,                              & ! intent(in):    model parameters
                       prog_data,                              & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,                              & ! intent(in):    model diagnostic variables for a local HRU
                       flux_data,                              & ! intent(inout): model fluxes for a local HRU

                       ! output: baseflow
                       ixSaturation,                           & ! intent(inout) index of lowest saturated layer (NOTE: only computed on the first iteration)
                       mLayerBaseflow,                         & ! intent(out): baseflow from each soil layer (m s-1)
                       dBaseflow_dMatric,                      & ! intent(out): derivative in baseflow w.r.t. matric head (s-1)

                       ! output: error control
                       err,message)                              ! intent(out): error control
 ! ---------------------------------------------------------------------------------------
 ! provide access to the derived types to define the data structures
 USE data_types,only:&
                     var_d,            & ! data vector (dp)
                     var_dlength         ! data vector with variable length dimension (dp)
 ! provide access to named variables defining elements in the data structures
 USE var_lookup,only:iLookATTR           ! named variables for structure elements
 USE var_lookup,only:iLookPROG           ! named variables for structure elements
 USE var_lookup,only:iLookDIAG           ! named variables for structure elements
 USE var_lookup,only:iLookFLUX           ! named variables for structure elements
 USE var_lookup,only:iLookPARAM          ! named variables for structure elements
 ! utility modules
 USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water as a function of matric head
 USE soil_utils_module,only:hydCond_psi         ! compute hydraulic conductivity as a function of matric head
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: model control
 integer(i4b),intent(in)          :: nSnow                        ! number of snow layers
 integer(i4b),intent(in)          :: nSoil                        ! number of soil layers
 integer(i4b),intent(in)          :: nLayers                      ! total number of layers
 logical(lgt),intent(in)          :: getSatDepth                  ! logical flag to compute index of the lowest saturated layer
 ! input: state and diagnostic variables
 real(dp),intent(in)              :: mLayerdTheta_dPsi(:)         ! derivative in the soil water characteristic w.r.t. matric head in each layer (m-1)
 real(dp),intent(in)              :: mLayerMatricHeadLiq(:)       ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)              :: mLayerVolFracLiq(:)          ! volumetric fraction of liquid water (-)
 real(dp),intent(in)              :: mLayerVolFracIce(:)          ! volumetric fraction of ice (-)
 ! input/output: data structures
 type(var_d),intent(in)           :: attr_data                    ! spatial attributes
 type(var_dlength),intent(in)     :: mpar_data                    ! model parameters
 type(var_dlength),intent(in)     :: prog_data                    ! prognostic variables for a local HRU
 type(var_dlength),intent(in)     :: diag_data                    ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout)  :: flux_data                    ! model fluxes for a local HRU
 ! output: baseflow
 integer(i4b),intent(inout)       :: ixSaturation                 ! index of lowest saturated layer (NOTE: only computed on the first iteration)
 real(dp),intent(out)             :: mLayerBaseflow(:)            ! baseflow from each soil layer (m s-1)
 real(dp),intent(out)             :: dBaseflow_dMatric(:,:)       ! derivative in baseflow w.r.t. matric head (s-1)
 ! output: error control
 integer(i4b),intent(out)         :: err                          ! error code
 character(*),intent(out)         :: message                      ! error message
 ! ---------------------------------------------------------------------------------------
 ! * local variables
 ! ---------------------------------------------------------------------------------------
 ! general local variables
 integer(i4b)                    :: iLayer                        ! index of soil layer
 real(dp),dimension(nSoil,nSoil) :: dBaseflow_dVolLiq             ! derivative in the baseflow flux w.r.t. volumetric liquid water content (m s-1)
 ! local variables to compute the numerical Jacobian
 logical(lgt),parameter          :: doNumericalJacobian=.false.   ! flag to compute the numerical Jacobian
 real(dp),dimension(nSoil)       :: mLayerMatricHeadPerturbed     ! perturbed matric head (m)
 real(dp),dimension(nSoil)       :: mLayerVolFracLiqPerturbed     ! perturbed volumetric fraction of liquid water (-)
 real(dp),dimension(nSoil)       :: mLayerBaseflowPerturbed       ! perturbed baseflow (m s-1)
 real(dp),dimension(nSoil,nSoil) :: nJac                          ! numerical Jacobian (s-1)
 ! ***************************************************************************************
 ! ***************************************************************************************
 ! initialize error control
 err=0; message='groundwatr/'
 ! ---------------------------------------------------------------------------------------
 ! ---------------------------------------------------------------------------------------
 ! associate variables in data structures
 associate(&

 ! input: baseflow parameters
 fieldCapacity           => mpar_data%var(iLookPARAM%fieldCapacity)%dat(1),         & ! intent(in): [dp] field capacity (-)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                & ! intent(in): [dp] soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat,                & ! intent(in): [dp] residual volumetric water content (-)

 ! input: van Genuchten soil parametrers
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat,                & ! intent(in): [dp] van Genutchen "alpha" parameter (m-1)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat,                    & ! intent(in): [dp] van Genutchen "n" parameter (-)
 vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat,               & ! intent(in): [dp] van Genutchen "m" parameter (-)

 ! output: diagnostic variables
 scalarExfiltration      => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1),     & ! intent(out):[dp]    exfiltration from the soil profile (m s-1)
 mLayerColumnOutflow     => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat        & ! intent(out):[dp(:)] column outflow from each soil layer (m3 s-1)

 )  ! end association to variables in data structures

 ! ************************************************************************************************
 ! (1) compute the "active" portion of the soil profile
 ! ************************************************************************************************

 ! get index of the lowest saturated layer
 if(getSatDepth)then  ! NOTE: only compute for the first flux call
  ixSaturation = nSoil+1  ! unsaturated profile when ixSaturation>nSoil
  do iLayer=nSoil,1,-1  ! start at the lowest soil layer and work upwards to the top layer
   if(mLayerVolFracLiq(iLayer) > fieldCapacity)then; ixSaturation = iLayer  ! index of saturated layer -- keeps getting over-written as move upwards
   else; exit; end if                                                        ! (only consider saturated layer at the bottom of the soil profile)
  end do  ! (looping through soil layers)
 end if

 ! check for an early return (no layers are "active")
 if(ixSaturation > nSoil)then
  scalarExfiltration     = 0._dp   ! exfiltration from the soil profile (m s-1)
  mLayerColumnOutflow(:) = 0._dp   ! column outflow from each soil layer (m3 s-1)
  mLayerBaseflow(:)      = 0._dp   ! baseflow from each soil layer (m s-1)
  dBaseflow_dMatric(:,:) = 0._dp   ! derivative in baseflow w.r.t. matric head (s-1)
  return
 end if  ! if some layers are saturated

 ! ************************************************************************************************
 ! (2) compute the baseflow flux and its derivative w.r.t volumetric liquid water content
 ! ************************************************************************************************

 ! use private subroutine to compute baseflow (for multiple calls for numerical Jacobian)
 call computeBaseflow(&
                      ! input: control and state variables
                      nSnow,                   & ! intent(in): number of snow layers
                      nSoil,                   & ! intent(in): number of soil layers
                      nLayers,                 & ! intent(in): total number of layers
                      .true.,                  & ! intent(in): .true. if analytical derivatives are desired
                      ixSaturation,            & ! intent(in): index of upper-most "saturated" layer
                      mLayerVolFracLiq,        & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                      mLayerVolFracIce,        & ! intent(in): volumetric fraction of ice in each soil layer (-)
                      ! input/output: data structures
                      attr_data,               & ! intent(in):    spatial attributes
                      mpar_data,               & ! intent(in):    model parameters
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      flux_data,               & ! intent(inout): model fluxes for a local HRU
                      ! output: fluxes and derivatives
                      mLayerBaseflow,          & ! intent(out): baseflow flux in each soil layer (m s-1)
                      dBaseflow_dVolLiq)         ! intent(out): derivative in baseflow w.r.t. volumetric liquid water content (s-1)

 ! use the chain rule to compute the baseflow derivative w.r.t. matric head (s-1)
 do iLayer=1,nSoil
  dBaseflow_dMatric(1:iLayer,iLayer) = dBaseflow_dVolLiq(1:iLayer,iLayer)*mLayerdTheta_dPsi(iLayer)
  if(iLayer<nSoil) dBaseflow_dMatric(iLayer+1:nSoil,iLayer) = 0._dp
 end do

 ! ************************************************************************************************
 ! (4) compute the numerical Jacobian
 ! ************************************************************************************************
 if(doNumericalJacobian)then

  ! first, print the analytical Jacobian
  print*, '** analytical Jacobian:'
  write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=1,nSoil)
  do iLayer=1,nSoil; write(*,'(i4,1x,100(e12.5,1x))') iLayer, dBaseflow_dMatric(1:nSoil,iLayer); end do

  ! get a copy of the state vector to perturb
  mLayerMatricHeadPerturbed(:) = mLayerMatricHeadLiq(:)
  mLayerVolFracLiqPerturbed(:) = mLayerVolFracLiq(:)

  ! loop through the state variables
  do iLayer=1,nSoil

   ! perturb state vector
   mLayerMatricHeadPerturbed(iLayer) = mLayerMatricHeadPerturbed(iLayer) + dx

   ! compute the columetruc liquid water content
   mLayerVolFracLiqPerturbed(iLayer) = volFracLiq(mLayerMatricHeadPerturbed(iLayer),vGn_alpha(iLayer),theta_res(iLayer),theta_sat(iLayer),vGn_n(iLayer),vGn_m(iLayer))

   ! compute baseflow flux
   ! NOTE: This is an optional second call to computeBaseflow that is invoked when computing numerical derivatives.
   !       Since the purpose here is to compute the numerical derivatives, we do not need to compute analytical derivatives also.
   !       Hence, analytical derivatives are not desired
   call computeBaseflow(&
                        ! input: control and state variables
                        nSnow,                     & ! intent(in): number of snow layers
                        nSoil,                     & ! intent(in): number of soil layers
                        nLayers,                   & ! intent(in): total number of layers
                        .false.,                   & ! intent(in): .true. if analytical derivatives are desired
                        ixSaturation,              & ! intent(in): index of upper-most "saturated" layer
                        mLayerVolFracLiqPerturbed, & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                        mLayerVolFracIce,          & ! intent(in): volumetric fraction of ice in each soil layer (-)
                        ! input/output: data structures
                        attr_data,                 & ! intent(in):    spatial attributes
                        mpar_data,                 & ! intent(in):    model parameters
                        prog_data,                 & ! intent(in):    model prognostic variables for a local HRU
                        flux_data,                 & ! intent(inout): model fluxes for a local HRU
                        ! output: fluxes and derivatives
                        mLayerBaseflowPerturbed,   & ! intent(out): baseflow flux in each soil layer (m s-1)
                        dBaseflow_dVolLiq)           ! intent(out): ** NOT USED ** derivative in baseflow w.r.t. volumetric liquid water content (s-1)

   ! compute the numerical Jacobian
   nJac(:,iLayer) = (mLayerBaseflowPerturbed(:) - mLayerBaseflow(:))/dx    ! compute the Jacobian

   ! set the state back to the input value
   mLayerMatricHeadPerturbed(iLayer) = mLayerMatricHeadLiq(iLayer)
   mLayerVolFracLiqPerturbed(iLayer) = mLayerVolFracLiq(iLayer)

  end do  ! looping through state variables

  ! print the numerical Jacobian
  print*, '** numerical Jacobian:'
  write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=1,nSoil)
  do iLayer=1,nSoil; write(*,'(i4,1x,100(e12.5,1x))') iLayer, nJac(1:nSoil,iLayer); end do
  !pause 'testing Jacobian'

 end if  ! if desire to compute the Jacobian

 ! end association to variables in data structures
 end associate

 end subroutine groundwatr


 ! ***********************************************************************************************************************
 ! * private subroutine computeBaseflow: private subroutine so can be used to test the numerical jacobian
 ! ***********************************************************************************************************************
 subroutine computeBaseflow(&
                            ! input: control and state variables
                            nSnow,                         & ! intent(in): number of snow layers
                            nSoil,                         & ! intent(in): number of soil layers
                            nLayers,                       & ! intent(in): total number of layers
                            derivDesired,                  & ! intent(in): .true. if derivatives are desired
                            ixSaturation,                  & ! intent(in): index of upper-most "saturated" layer
                            mLayerVolFracLiq,              & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                            mLayerVolFracIce,              & ! intent(in): volumetric fraction of ice in each soil layer (-)
                            ! input/output: data structures
                            attr_data,                     & ! intent(in):    spatial attributes
                            mpar_data,                     & ! intent(in):    model parameters
                            prog_data,                     & ! intent(in):    model prognostic variables for a local HRU
                            flux_data,                     & ! intent(inout): model fluxes for a local HRU
                            ! output: fluxes and derivatives
                            mLayerBaseflow,                & ! intent(out): baseflow flux in each soil layer (m s-1)
                            dBaseflow_dVolLiq)               ! intent(out): derivative in baseflow w.r.t. volumetric liquid water content (s-1)
 implicit none
 ! ---------------------------------------------------------------------------------------
 ! * dummy variables
 ! ---------------------------------------------------------------------------------------
 ! input: control and state variables
 integer(i4b),intent(in)          :: nSnow                   ! number of snow layers
 integer(i4b),intent(in)          :: nSoil                   ! number of soil layers
 integer(i4b),intent(in)          :: nLayers                 ! total number of layers
 logical(lgt),intent(in)          :: derivDesired            ! .true. if derivatives are desired
 integer(i4b),intent(in)          :: ixSaturation            ! index of upper-most "saturated" layer
 real(dp),intent(in)              :: mLayerVolFracLiq(:)     ! volumetric fraction of liquid water (-)
 real(dp),intent(in)              :: mLayerVolFracIce(:)     ! volumetric fraction of ice (-)
 ! input/output: data structures
 type(var_d),intent(in)           :: attr_data               ! spatial attributes
 type(var_dlength),intent(in)     :: mpar_data               ! model parameters
 type(var_dlength),intent(in)     :: prog_data               ! prognostic variables for a local HRU
 type(var_dlength),intent(inout)  :: flux_data               ! model fluxes for a local HRU
 ! output: baseflow
 real(dp),intent(out)             :: mLayerBaseflow(:)       ! baseflow from each soil layer (m s-1)
 real(dp),intent(out)             :: dBaseflow_dVolLiq(:,:)  ! derivative in baseflow w.r.t. matric head (s-1)
 ! ---------------------------------------------------------------------------------------
 ! * local variables
 ! ---------------------------------------------------------------------------------------
 ! general local variables
 integer(i4b)                    :: iLayer,jLayer            ! index of model layer
 ! local variables for the exfiltration
 real(dp)                        :: totalColumnInflow        ! total column inflow (m s-1)
 real(dp)                        :: totalColumnOutflow       ! total column outflow (m s-1)
 real(dp)                        :: availStorage             ! available storage (m)
 real(dp),parameter              :: xMinEval=0.002_dp        ! minimum value to evaluate the exfiltration function (m)
 real(dp),parameter              :: xCenter=0.001_dp         ! center of the exfiltration function (m)
 real(dp),parameter              :: xWidth=0.0001_dp         ! width of the exfiltration function (m)
 real(dp)                        :: expF,logF                ! logistic smoothing function (-)
 ! local variables for the lateral flux among soil columns
 real(dp)                        :: activePorosity           ! "active" porosity associated with storage above a threshold (-)
 real(dp)                        :: drainableWater           ! drainable water in eaxch layer (m)
 real(dp)                        :: tran0                    ! maximum transmissivity (m2 s-1)
 real(dp),dimension(nSoil)       :: zActive                  ! water table thickness associated with storage below and including the given layer (m)
 real(dp),dimension(nSoil)       :: trTotal                  ! total transmissivity associated with total water table depth zActive (m2 s-1)
 real(dp),dimension(nSoil)       :: trSoil                   ! transmissivity of water in a given layer (m2 s-1)
 ! local variables for the derivatives
 real(dp)                        :: qbTotal                  ! total baseflow (m s-1)
 real(dp)                        :: length2area              ! ratio of hillslope width to hillslope area (m m-2)
 real(dp),dimension(nSoil)       :: depth2capacity           ! ratio of layer depth to total subsurface storage capacity (-)
 real(dp),dimension(nSoil)       :: dXdS                     ! change in dimensionless flux w.r.t. change in dimensionless storage (-)
 real(dp),dimension(nSoil)       :: dLogFunc_dLiq            ! derivative in the logistic function w.r.t. volumetric liquid water content (-)
 real(dp),dimension(nSoil)       :: dExfiltrate_dVolLiq      ! derivative in exfiltration w.r.t. volumetric liquid water content (-)
 ! local variables for testing (debugging)
 logical(lgt),parameter          :: printFlag=.false.        ! flag for printing (debugging)
 real(dp)                        :: xDepth,xTran,xFlow       ! temporary variables (depth, transmissivity, flow)
 ! ---------------------------------------------------------------------------------------
 ! * association to data in structures
 ! ---------------------------------------------------------------------------------------
 associate(&

 ! input: coordinate variables
 soilDepth               => prog_data%var(iLookPROG%iLayerHeight)%dat(nLayers),       & ! intent(in): [dp]    total soil depth (m)
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat(nSnow+1:nLayers),& ! intent(in): [dp(:)] depth of each soil layer (m)

 ! input: diagnostic variables
 surfaceHydCond          => flux_data%var(iLookFLUX%mLayerSatHydCondMP)%dat(1),       & ! intent(in): [dp]    saturated hydraulic conductivity at the surface (m s-1)
 mLayerColumnInflow      => flux_data%var(iLookFLUX%mLayerColumnInflow)%dat,          & ! intent(in): [dp(:)] inflow into each soil layer (m3/s)

 ! input: local attributes
 HRUarea                 => attr_data%var(iLookATTR%HRUarea),                         & ! intent(in): [dp] HRU area (m2)
 tan_slope               => attr_data%var(iLookATTR%tan_slope),                       & ! intent(in): [dp] tan water table slope, taken as tan local ground surface slope (-)
 contourLength           => attr_data%var(iLookATTR%contourLength),                   & ! intent(in): [dp] length of contour at downslope edge of HRU (m)

 ! input: baseflow parameters
 zScale_TOPMODEL         => mpar_data%var(iLookPARAM%zScale_TOPMODEL)%dat(1),         & ! intent(in): [dp]    TOPMODEL exponent (-)
 kAnisotropic            => mpar_data%var(iLookPARAM%kAnisotropic)%dat(1),            & ! intent(in): [dp]    anisotropy factor for lateral hydraulic conductivity (-
 fieldCapacity           => mpar_data%var(iLookPARAM%fieldCapacity)%dat(1),           & ! intent(in): [dp]    field capacity (-)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat,                  & ! intent(in): [dp(:)] soil porosity (-)

 ! output: diagnostic variables
 scalarExfiltration      => flux_data%var(iLookFLUX%scalarExfiltration)%dat(1),       & ! intent(out):[dp]    exfiltration from the soil profile (m s-1)
 mLayerColumnOutflow     => flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat          & ! intent(out):[dp(:)] column outflow from each soil layer (m3 s-1)

 )  ! end association to variables in data structures
 ! ***********************************************************************************************************************
 ! ***********************************************************************************************************************
 ! start routine here

 ! ***********************************************************************************************************************
 ! (1) compute the baseflow flux in each soil layer
 ! ***********************************************************************************************************************

 ! compute the maximum transmissivity
 ! NOTE: this can be done as a pre-processing step
 tran0 = kAnisotropic*surfaceHydCond*soilDepth/zScale_TOPMODEL   ! maximum transmissivity (m2 s-1)

 ! compute the water table thickness (m) and transmissivity in each layer (m2 s-1)
 do iLayer=nSoil,ixSaturation,-1  ! loop through "active" soil layers, from lowest to highest
  ! define drainable water in each layer (m)
  activePorosity = theta_sat(iLayer) - fieldCapacity ! "active" porosity (-)
  drainableWater = mLayerDepth(iLayer)*(max(0._dp,mLayerVolFracLiq(iLayer) - fieldCapacity))/activePorosity
  ! compute layer transmissivity
  if(iLayer==nSoil)then
   zActive(iLayer) = drainableWater                                       ! water table thickness associated with storage in a given layer (m)
   trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL   ! total transmissivity for total depth zActive (m2 s-1)
   trSoil(iLayer)  = trTotal(iLayer)                                      ! transmissivity of water in a given layer (m2 s-1)
  else
   zActive(iLayer) = zActive(iLayer+1) + drainableWater
   trTotal(iLayer) = tran0*(zActive(iLayer)/soilDepth)**zScale_TOPMODEL
   trSoil(iLayer)  = trTotal(iLayer) - trTotal(iLayer+1)
  end if
  !write(*,'(a,1x,i4,1x,10(f20.15,1x))') 'iLayer, mLayerMatricHeadLiq(iLayer), mLayerVolFracLiq(iLayer), zActive(iLayer), trTotal(iLayer), trSoil(iLayer) = ', &
  !                                       iLayer, mLayerMatricHeadLiq(iLayer), mLayerVolFracLiq(iLayer), zActive(iLayer), trTotal(iLayer), trSoil(iLayer)
 end do  ! looping through soil layers

 ! set un-used portions of the vectors to zero
 if(ixSaturation>1)then
  zActive(1:ixSaturation-1) = 0._dp
  trTotal(1:ixSaturation-1) = 0._dp
  trSoil(1:ixSaturation-1)  = 0._dp
 end if

 ! compute the outflow from each layer (m3 s-1)
 mLayerColumnOutflow(1:nSoil) = trSoil(1:nSoil)*tan_slope*contourLength

 ! compute total column inflow and total column outflow (m s-1)
 totalColumnInflow  = sum(mLayerColumnInflow(1:nSoil))/HRUarea
 totalColumnOutflow = sum(mLayerColumnOutflow(1:nSoil))/HRUarea

 ! compute the available storage (m)
 availStorage = sum(mLayerDepth(1:nSoil)*(theta_sat - (mLayerVolFracLiq(1:nSoil)+mLayerVolFracIce(1:nSoil))) )

 ! compute the smoothing function (-)
 if(availStorage < xMinEval)then
  ! (compute the logistic function)
  expF = exp((availStorage - xCenter)/xWidth)
  logF = 1._dp / (1._dp + expF)
  ! (compute the derivative in the logistic function w.r.t. volumetric liquid water content in each soil layer)
  dLogFunc_dLiq(1:nSoil) = mLayerDepth(1:nSoil)*(expF/xWidth)/(1._dp + expF)**2._dp
 else
  logF             = 0._dp
  dLogFunc_dLiq(:) = 0._dp
 end if

 ! compute the exfiltartion (m s-1)
 if(totalColumnInflow > totalColumnOutflow .and. logF > tiny(1._dp))then
  scalarExfiltration = logF*(totalColumnInflow - totalColumnOutflow)  ! m s-1
  !write(*,'(a,1x,10(f30.20,1x))') 'scalarExfiltration = ', scalarExfiltration
 else
  scalarExfiltration = 0._dp
 end if

 ! check
 !write(*,'(a,1x,10(f30.20,1x))') 'zActive(1), soilDepth, availStorage, logF, scalarExfiltration = ', &
 !                                 zActive(1), soilDepth, availStorage, logF, scalarExfiltration
 !if(scalarExfiltration > tiny(1.0_dp))  pause 'exfiltrating'

 ! compute the baseflow in each layer (m s-1)
 mLayerBaseflow(1:nSoil) = (mLayerColumnOutflow(1:nSoil) - mLayerColumnInflow(1:nSoil))/HRUarea

 ! compute the total baseflow
 qbTotal = sum(mLayerBaseflow)

 ! add exfiltration to the baseflow flux at the top layer
 mLayerBaseflow(1)      = mLayerBaseflow(1) + scalarExfiltration
 mLayerColumnOutflow(1) = mLayerColumnOutflow(1) + scalarExfiltration*HRUarea

 ! test
 if(printFlag)then
  xDepth = sum(mLayerDepth(ixSaturation:nSoil)*(mLayerVolFracLiq(ixSaturation:nSoil) - fieldCapacity))/sum(theta_sat(ixSaturation:nSoil) - fieldCapacity)  ! "effective" water table thickness (m)
  xTran  = tran0*(xDepth/soilDepth)**zScale_TOPMODEL  ! transmissivity for the entire aquifer (m2 s-1)
  xFlow  = xTran*tan_slope*contourLength/HRUarea   ! total column outflow (m s-1)
  print*, 'ixSaturation = ', ixSaturation
  write(*,'(a,1x,5(f30.20,1x))') 'tran0, soilDepth                 = ', tran0, soilDepth
  write(*,'(a,1x,5(f30.20,1x))') 'surfaceHydCond, zScale_TOPMODEL  = ', surfaceHydCond, zScale_TOPMODEL
  write(*,'(a,1x,5(f30.20,1x))') 'xDepth, zActive(ixSaturation)    = ', xDepth, zActive(ixSaturation)
  write(*,'(a,1x,5(f30.20,1x))') 'xTran, trTotal(ixSaturation)     = ', xTran, trTotal(ixSaturation)
  write(*,'(a,1x,5(f30.20,1x))') 'xFlow, totalColumnOutflow        = ', xFlow, sum(mLayerColumnOutflow(:))/HRUarea
  !pause 'check groundwater'
 end if

 ! ***********************************************************************************************************************
 ! (2) compute the derivative in the baseflow flux w.r.t. volumetric liquid water content (m s-1)
 ! ***********************************************************************************************************************

 ! initialize the derivative matrix
 dBaseflow_dVolLiq(:,:) = 0._dp

 ! check if derivatives are actually required
 if(.not.derivDesired) return

 ! compute ratio of hillslope width to hillslope area (m m-2)
 length2area = tan_slope*contourLength/HRUarea

 ! compute the ratio of layer depth to maximum water holding capacity (-)
 depth2capacity(1:nSoil) = mLayerDepth(1:nSoil)/sum( (theta_sat(1:nSoil) - fieldCapacity)*mLayerDepth(1:nSoil) )

 ! compute the change in dimensionless flux w.r.t. change in dimensionless storage (-)
 dXdS(1:nSoil) = zScale_TOPMODEL*(zActive(1:nSoil)/SoilDepth)**(zScale_TOPMODEL - 1._dp)

 ! loop through soil layers
 do iLayer=1,nSoil
  ! compute diagonal terms (s-1)
  dBaseflow_dVolLiq(iLayer,iLayer) = tran0*dXdS(iLayer)*depth2capacity(iLayer)*length2area
  ! compute off-diagonal terms
  do jLayer=iLayer+1,nSoil  ! (only dependent on layers below)
   dBaseflow_dVolLiq(iLayer,jLayer) = tran0*(dXdS(iLayer) - dXdS(iLayer+1))*depth2capacity(jLayer)*length2area
  end do  ! looping through soil layers
 end do  ! looping through soil layers

 ! compute the derivative in the exfiltration flux w.r.t. volumetric liquid water content (m s-1)
 if(qbTotal < 0._dp)then
  do iLayer=1,nSoil
   dExfiltrate_dVolLiq(iLayer) = dBaseflow_dVolLiq(iLayer,iLayer)*logF + dLogFunc_dLiq(iLayer)*qbTotal
  end do  ! looping through soil layers
  dBaseflow_dVolLiq(1,1:nSoil) = dBaseflow_dVolLiq(1,1:nSoil) - dExfiltrate_dVolLiq(1:nSoil)
 end if

 ! end association to data in structures
 end associate

 end subroutine computeBaseflow


end module groundwatr_module
