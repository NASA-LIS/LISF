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

module computJacob_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! domain types
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands
USE globalData,only: ixDiag         ! index for the diagonal band
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (dp)

implicit none
! define constants
real(dp),parameter     :: verySmall=tiny(1.0_dp)     ! a very small number
integer(i4b),parameter :: ixBandOffset=kl+ku+1       ! offset in the band Jacobian matrix

private
public::computJacob
contains

 ! **********************************************************************************************************
 ! public subroutine computJacob: compute the Jacobian matrix
 ! **********************************************************************************************************
 subroutine computJacob(&
                        ! input: model control
                        dt,                         & ! intent(in):    length of the time step (seconds)
                        nSnow,                      & ! intent(in):    number of snow layers
                        nSoil,                      & ! intent(in):    number of soil layers
                        nLayers,                    & ! intent(in):    total number of layers
                        computeVegFlux,             & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                        computeBaseflow,            & ! intent(in):    flag to indicate if we need to compute baseflow
                        ixMatrix,                   & ! intent(in):    form of the Jacobian matrix
                        ! input: data structures
                        indx_data,                  & ! intent(in):    index data
                        prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,                  & ! intent(in):    model diagnostic variables for a local HRU
                        deriv_data,                 & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                        dBaseflow_dMatric,          & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                        ! input-output: Jacobian and its diagonal
                        dMat,                       & ! intent(inout): diagonal of the Jacobian matrix
                        aJac,                       & ! intent(out):   Jacobian matrix
                        ! output: error control
                        err,message)                  ! intent(out):   error code and error message
 ! named variables for structure elements
 USE var_lookup,only:iLookPROG                        ! named variables for structure elements
 USE var_lookup,only:iLookDIAG                        ! named variables for structure elements
 USE var_lookup,only:iLookINDEX                       ! named variables for structure elements
 USE var_lookup,only:iLookDERIV                       ! named variables for structure elements
 ! -----------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 real(dp),intent(in)               :: dt              ! length of the time step (seconds)
 integer(i4b),intent(in)           :: nSnow           ! number of snow layers
 integer(i4b),intent(in)           :: nSoil           ! number of soil layers
 integer(i4b),intent(in)           :: nLayers         ! total number of layers in the snow+soil domain
 logical(lgt),intent(in)           :: computeVegFlux  ! flag to indicate if computing fluxes over vegetation
 logical(lgt),intent(in)           :: computeBaseflow ! flag to indicate if computing baseflow
 integer(i4b),intent(in)           :: ixMatrix        ! form of the Jacobian matrix
 ! input: data structures
 type(var_ilength),intent(in)      :: indx_data       ! indices defining model states and layers
 type(var_dlength),intent(in)      :: prog_data       ! prognostic variables for a local HRU
 type(var_dlength),intent(in)      :: diag_data       ! diagnostic variables for a local HRU
 type(var_dlength),intent(in)      :: deriv_data      ! derivatives in model fluxes w.r.t. relevant state variables
 real(dp),intent(in)               :: dBaseflow_dMatric(:,:) ! derivative in baseflow w.r.t. matric head (s-1)
 ! input-output: Jacobian and its diagonal
 real(dp),intent(inout)            :: dMat(:)         ! diagonal of the Jacobian matrix
 real(dp),intent(out)              :: aJac(:,:)       ! Jacobian matrix
 ! output variables
 integer(i4b),intent(out)          :: err             ! error code
 character(*),intent(out)          :: message         ! error message
 ! --------------------------------------------------------------
 ! * local variables
 ! --------------------------------------------------------------
 ! indices of model state variables
 integer(i4b)                      :: jState          ! index of state within the state subset
 integer(i4b)                      :: qState          ! index of cross-derivative state variable for baseflow
 integer(i4b)                      :: nrgState        ! energy state variable
 integer(i4b)                      :: watState        ! hydrology state variable
 integer(i4b)                      :: nState          ! number of state variables
 ! indices of model layers
 integer(i4b)                      :: iLayer          ! index of model layer
 integer(i4b)                      :: jLayer          ! index of model layer within the full state vector (hydrology)
 integer(i4b)                      :: pLayer          ! indices of soil layers (used for the baseflow derivatives)
 ! conversion factors
 real(dp)                          :: convLiq2tot     ! factor to convert liquid water derivative to total water derivative
 ! --------------------------------------------------------------
 ! associate variables from data structures
 associate(&
 ! indices of model state variables
 ixCasNrg                     => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)                       ,& ! intent(in): [i4b(:)] index of canopy air space energy state variable
 ixVegNrg                     => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)                       ,& ! intent(in): [i4b(:)] index of canopy energy state variable
 ixVegHyd                     => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)                       ,& ! intent(in): [i4b(:)] index of canopy hydrology state variable (mass)
 ixTopNrg                     => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)                       ,& ! intent(in): [i4b(:)] index of upper-most energy state in the snow+soil subdomain
 ixTopHyd                     => indx_data%var(iLookINDEX%ixTopHyd)%dat(1)                       ,& ! intent(in): [i4b(:)] index of upper-most hydrology state in the snow+soil subdomain
 ! vectors of indices for specfic state types within specific sub-domains IN THE FULL STATE VECTOR
 ixNrgLayer                   => indx_data%var(iLookINDEX%ixNrgLayer)%dat                        ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer                   => indx_data%var(iLookINDEX%ixHydLayer)%dat                        ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ! vector of energy indices for the snow and soil domains
 ! NOTE: states not in the subset are equal to integerMissing
 ixSnowSoilNrg                => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowOnlyNrg                => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow domain
 ixSoilOnlyNrg                => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
 ! vector of hydrology indices for the snow and soil domains
 ! NOTE: states not in the subset are equal to integerMissing
 ixSnowSoilHyd                => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 ixSnowOnlyHyd                => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow domain
 ixSoilOnlyHyd                => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
 ! number of state variables of a specific type
 nSnowSoilNrg                 => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
 nSnowOnlyNrg                 => indx_data%var(iLookINDEX%nSnowOnlyNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
 nSoilOnlyNrg                 => indx_data%var(iLookINDEX%nSoilOnlyNrg )%dat(1)                  ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
 nSnowSoilHyd                 => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
 nSnowOnlyHyd                 => indx_data%var(iLookINDEX%nSnowOnlyHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
 nSoilOnlyHyd                 => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)                  ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
 ! type and index of model control volume
 ixHydType                    => indx_data%var(iLookINDEX%ixHydType)%dat                         ,& ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
 ixDomainType                 => indx_data%var(iLookINDEX%ixDomainType)%dat                      ,& ! intent(in): [i4b(:)] indices defining the type of the domain (iname_veg, iname_snow, iname_soil)
 ixControlVolume              => indx_data%var(iLookINDEX%ixControlVolume)%dat                   ,& ! intent(in): [i4b(:)] index of the control volume for specific model domains
 ! mapping between states and model layers
 ixMapSubset2Full             => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat                  ,& ! intent(in): [i4b(:)] list of indices in the full state vector that are in the state subset
 ixMapFull2Subset             => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat                  ,& ! intent(in): [i4b(:)] list of indices in the state subset in each element of the full state vector
 ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
 dCanairNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanairTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. canopy air temperature
 dCanairNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dCanopyTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. canopy temperature
 dCanairNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanairNetFlux_dGroundTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. ground temperature
 dCanopyNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanairTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. canopy air temperature
 dCanopyNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanopyTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. canopy temperature
 dCanopyNetFlux_dGroundTemp   => deriv_data%var(iLookDERIV%dCanopyNetFlux_dGroundTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. ground temperature
 dCanopyNetFlux_dCanLiq       => deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanLiq      )%dat(1)  ,& ! intent(in): [dp]     derivative in net canopy fluxes w.r.t. canopy liquid water content
 dGroundNetFlux_dCanairTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanairTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net ground flux w.r.t. canopy air temperature
 dGroundNetFlux_dCanopyTemp   => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanopyTemp  )%dat(1)  ,& ! intent(in): [dp]     derivative in net ground flux w.r.t. canopy temperature
 dGroundNetFlux_dCanLiq       => deriv_data%var(iLookDERIV%dGroundNetFlux_dCanLiq      )%dat(1)  ,& ! intent(in): [dp]     derivative in net ground fluxes w.r.t. canopy liquid water content
 ! derivatives in evaporative fluxes w.r.t. relevant state variables
 dCanopyEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanair )%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy air temperature
 dCanopyEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanopy )%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy temperature
 dCanopyEvaporation_dTGround  => deriv_data%var(iLookDERIV%dCanopyEvaporation_dTGround )%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. ground temperature
 dCanopyEvaporation_dCanLiq   => deriv_data%var(iLookDERIV%dCanopyEvaporation_dCanLiq  )%dat(1)  ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy liquid water content
 dGroundEvaporation_dTCanair  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanair )%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy air temperature
 dGroundEvaporation_dTCanopy  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanopy )%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy temperature
 dGroundEvaporation_dTGround  => deriv_data%var(iLookDERIV%dGroundEvaporation_dTGround )%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. ground temperature
 dGroundEvaporation_dCanLiq   => deriv_data%var(iLookDERIV%dGroundEvaporation_dCanLiq  )%dat(1)  ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy liquid water content
 ! derivatives in canopy water w.r.t canopy temperature
 dCanLiq_dTcanopy             => deriv_data%var(iLookDERIV%dCanLiq_dTcanopy            )%dat(1)  ,& ! intent(in): [dp]     derivative of canopy liquid storage w.r.t. temperature
 dTheta_dTkCanopy             => deriv_data%var(iLookDERIV%dTheta_dTkCanopy            )%dat(1)  ,& ! intent(in): [dp]     derivative of volumetric liquid water content w.r.t. temperature
 ! derivatives in canopy liquid fluxes w.r.t. canopy water
 scalarCanopyLiqDeriv         => deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv        )%dat(1)  ,& ! intent(in): [dp]     derivative in (throughfall + drainage) w.r.t. canopy liquid water
 ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
 dNrgFlux_dTempAbove          => deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove         )%dat     ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer above
 dNrgFlux_dTempBelow          => deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow         )%dat     ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer below
 ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
 iLayerLiqFluxSnowDeriv       => deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv      )%dat     ,& ! intent(in): [dp(:)]  derivative in vertical liquid water flux at layer interfaces
 ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
 dVolTot_dPsi0                => deriv_data%var(iLookDERIV%dVolTot_dPsi0               )%dat     ,& ! intent(in): [dp(:)]  derivative in total water content w.r.t. total water matric potential
 dq_dHydStateAbove            => deriv_data%var(iLookDERIV%dq_dHydStateAbove           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer above
 dq_dHydStateBelow            => deriv_data%var(iLookDERIV%dq_dHydStateBelow           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer below
 dCompress_dPsi               => deriv_data%var(iLookDERIV%dCompress_dPsi              )%dat     ,& ! intent(in): [dp(:)]  derivative in compressibility w.r.t matric head
 ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
 dq_dNrgStateAbove            => deriv_data%var(iLookDERIV%dq_dNrgStateAbove           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer above
 dq_dNrgStateBelow            => deriv_data%var(iLookDERIV%dq_dNrgStateBelow           )%dat     ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer below
 mLayerdTheta_dTk             => deriv_data%var(iLookDERIV%mLayerdTheta_dTk            )%dat     ,& ! intent(in): [dp(:)]  derivative of volumetric liquid water content w.r.t. temperature
 ! diagnostic variables
 scalarFracLiqVeg             => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)                ,& ! intent(in): [dp]     fraction of liquid water on vegetation (-)
 scalarBulkVolHeatCapVeg      => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1)         ,& ! intent(in): [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
 mLayerFracLiqSnow            => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat                  ,& ! intent(in): [dp(:)]  fraction of liquid water in each snow layer (-)
 mLayerVolHtCapBulk           => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat                 ,& ! intent(in): [dp(:)]  bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)          
 scalarSoilControl            => diag_data%var(iLookDIAG%scalarSoilControl)%dat(1)               ,& ! intent(in): [dp]     soil control on infiltration, zero or one
 ! canopy and layer depth
 canopyDepth                  => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)               ,& ! intent(in): [dp   ]  canopy depth (m)
 mLayerDepth                  => prog_data%var(iLookPROG%mLayerDepth)%dat                         & ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
 ) ! making association with data in structures
 ! --------------------------------------------------------------
 ! initialize error control
 err=0; message='computJacob/'

 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************
 ! * PART 0: PRELIMINARIES (INITIALIZE JACOBIAN AND COMPUTE TIME-VARIABLE DIAGONAL TERMS)
 ! *********************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************

 ! get the number of state variables
 nState = size(dMat)

 ! initialize the Jacobian
 ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver
 aJac(:,:) = 0._dp  ! analytical Jacobian matrix

 ! compute terms in the Jacobian for vegetation (excluding fluxes)
 ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
 if(ixVegNrg/=integerMissing) dMat(ixVegNrg) = scalarBulkVolHeatCapVeg + LH_fus*iden_water*dTheta_dTkCanopy       ! volumetric heat capacity of the vegetation (J m-3 K-1)

 ! compute additional terms for the Jacobian for the snow-soil domain (excluding fluxes)
 ! NOTE: energy for snow+soil is computed *within* the iteration loop as it includes phase change
 do iLayer=1,nLayers
  if(ixSnowSoilNrg(iLayer)/=integerMissing) dMat(ixSnowSoilNrg(iLayer)) = mLayerVolHtCapBulk(iLayer) + LH_fus*iden_water*mLayerdTheta_dTk(iLayer)
 end do

 ! compute additional terms for the Jacobian for the soil domain (excluding fluxes)
 do iLayer=1,nSoil
  if(ixSoilOnlyHyd(iLayer)/=integerMissing) dMat(ixSoilOnlyHyd(iLayer)) = dVolTot_dPsi0(iLayer) + dCompress_dPsi(iLayer)
 end do

 ! define the form of the matrix
 select case(ixMatrix)

  ! *********************************************************************************************************************************************************
  ! *********************************************************************************************************************************************************
  ! * PART 1: BAND MATRIX
  ! *********************************************************************************************************************************************************
  ! *********************************************************************************************************************************************************
  case(ixBandMatrix)

   ! check
   if(size(aJac,1)/=nBands .or. size(aJac,2)/=size(dMat))then
    message=trim(message)//'unexpected shape of the Jacobian matrix: expect aJac(nBands,nState)'
    err=20; return
   end if

   ! -----
   ! * energy and liquid fluxes over vegetation...
   ! ---------------------------------------------
   if(computeVegFlux)then  ! (derivatives only defined when vegetation protrudes over the surface)
   
    ! * diagonal elements for the vegetation canopy (-)
    if(ixCasNrg/=integerMissing) aJac(ixDiag,ixCasNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + dMat(ixCasNrg)
    if(ixVegNrg/=integerMissing) aJac(ixDiag,ixVegNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + dMat(ixVegNrg)
    if(ixVegHyd/=integerMissing) aJac(ixDiag,ixVegHyd) = -scalarFracLiqVeg*(dCanopyEvaporation_dCanLiq - scalarCanopyLiqDeriv)*dt + 1._dp     ! ixVegHyd: CORRECT

    ! * cross-derivative terms w.r.t. canopy water
    if(ixVegHyd/=integerMissing)then
     ! cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
     if(ixCasNrg/=integerMissing) aJac(ixOffDiag(ixVegHyd,ixCasNrg),ixCasNrg) = -dCanopyEvaporation_dTCanair*dt                                                        ! ixCasNrg: CORRECT
     if(ixVegNrg/=integerMissing) aJac(ixOffDiag(ixVegHyd,ixVegNrg),ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy             ! ixVegNrg: CORRECT
     if(ixTopNrg/=integerMissing) aJac(ixOffDiag(ixVegHyd,ixTopNrg),ixTopNrg) = -dCanopyEvaporation_dTGround*dt                                                        ! ixTopNrg: CORRECT
     ! cross-derivative terms w.r.t. canopy water (kg-1 m2)
     if(ixTopHyd/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarFracLiqVeg*scalarCanopyLiqDeriv)/iden_water
     ! cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
     ! NOTE: dIce/dLiq = (1 - scalarFracLiqVeg); dIce*LH_fus/canopyDepth = J m-3; dLiq = kg m-2
     if(ixVegNrg/=integerMissing) aJac(ixOffDiag(ixVegNrg,ixVegHyd),ixVegHyd) = (dt/canopyDepth)   *(-dCanopyNetFlux_dCanLiq) - (1._dp - scalarFracLiqVeg)*LH_fus/canopyDepth   ! dF/dLiq
     if(ixTopNrg/=integerMissing) aJac(ixOffDiag(ixTopNrg,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanLiq)
    endif     

    ! cross-derivative terms between surface hydrology and the temperature of the vegetation canopy (K-1)
    if(ixVegNrg/=integerMissing)then
     if(ixTopHyd/=integerMissing) aJac(ixOffDiag(ixTopHyd,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water
    endif

    ! cross-derivative terms w.r.t. the temperature of the canopy air space (J m-3 K-1)
    if(ixCasNrg/=integerMissing)then
     if(ixVegNrg/=integerMissing) aJac(ixOffDiag(ixCasNrg,ixVegNrg),ixVegNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)
     if(ixTopNrg/=integerMissing) aJac(ixOffDiag(ixCasNrg,ixTopNrg),ixTopNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dGroundTemp)
    endif

    ! cross-derivative terms w.r.t. the temperature of the vegetation canopy (J m-3 K-1)
    if(ixVegNrg/=integerMissing)then
     if(ixCasNrg/=integerMissing) aJac(ixOffDiag(ixVegNrg,ixCasNrg),ixCasNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)
     if(ixTopNrg/=integerMissing) aJac(ixOffDiag(ixVegNrg,ixTopNrg),ixTopNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)
    endif

    ! cross-derivative terms w.r.t. the temperature of the surface (J m-3 K-1)
    if(ixTopNrg/=integerMissing)then
     if(ixCasNrg/=integerMissing) aJac(ixOffDiag(ixTopNrg,ixCasNrg),ixCasNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)
     if(ixVegNrg/=integerMissing) aJac(ixOffDiag(ixTopNrg,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)
    endif

   endif  ! if there is a need to compute energy fluxes within vegetation

   ! -----
   ! * energy fluxes for the snow+soil domain...
   ! -------------------------------------------
   if(nSnowSoilNrg>0)then
    do iLayer=1,nLayers  ! loop through all layers in the snow+soil domain

     ! check if the state is in the subset
     if(ixSnowSoilNrg(iLayer)==integerMissing) cycle

     ! - define index within the state subset and the full state vector
     jState = ixSnowSoilNrg(iLayer)        ! index within the state subset

     ! - diagonal elements
     aJac(ixDiag,jState)   = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(jState)

     ! - lower-diagonal elements
     if(iLayer > 1)then
      if(ixSnowSoilNrg(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSnowSoilNrg(iLayer-1),jState),jState) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
     endif

     ! - upper diagonal elements
     if(iLayer < nLayers)then
      if(ixSnowSoilNrg(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSnowSoilNrg(iLayer+1),jState),jState) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
     endif

    end do  ! (looping through energy states in the snow+soil domain)
   endif   ! (if the subset includes energy state variables in the snow+soil domain)

   ! -----
   ! * liquid water fluxes for the snow domain...
   ! --------------------------------------------
   if(nSnowOnlyHyd>0)then
    do iLayer=1,nSnow  ! loop through layers in the snow domain

     ! - check that the snow layer is desired
     if(ixSnowOnlyHyd(iLayer)==integerMissing) cycle

     ! - define state indices for the current layer
     watState = ixSnowOnlyHyd(iLayer)   ! hydrology state index within the state subset

     ! compute factor to convert liquid water derivative to total water derivative
     select case( ixHydType(iLayer) )
      case(iname_watLayer); convLiq2tot = mLayerFracLiqSnow(iLayer)
      case default;         convLiq2tot = 1._dp
     end select

     ! - diagonal elements
     aJac(ixDiag,watState) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot + dMat(watState)

     ! - lower-diagonal elements
     if(iLayer > 1)then
      if(ixSnowOnlyHyd(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSnowOnlyHyd(iLayer-1),watState),watState) = 0._dp  ! sub-diagonal: no dependence on other layers
     endif

     ! - upper diagonal elements
     if(iLayer < nSnow)then
      if(ixSnowOnlyHyd(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSnowOnlyHyd(iLayer+1),watState),watState) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot       ! dVol(below)/dLiq(above) -- (-)
     endif

     ! - compute cross-derivative terms for energy
     ! NOTE: increase in volumetric liquid water content balanced by a decrease in volumetric ice content
     if(nSnowOnlyNrg>0)then

      ! (define the energy state)
      nrgState = ixSnowOnlyNrg(iLayer)       ! index within the full state vector
      if(nrgstate/=integerMissing)then       ! (energy state for the current layer is within the state subset)

       ! (cross-derivative terms for the current layer)
       aJac(ixOffDiag(nrgState,watState),watState) = -(1._dp - mLayerFracLiqSnow(iLayer))*LH_fus*iden_water     ! (dF/dLiq)
       aJac(ixOffDiag(watState,nrgState),nrgState) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! (dVol/dT)

       ! (cross-derivative terms for the layer below)
       if(iLayer < nSnow)then
        aJac(ixOffDiag(ixSnowOnlyHyd(iLayer+1),nrgState),nrgState) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)        ! dVol(below)/dT(above) -- K-1
       endif ! (if there is a water state in the layer below the current layer in the given state subset)

      endif ! (if the energy state for the current layer is within the state subset)
     endif ! (if state variables exist for energy in snow+soil layers)

    end do  ! (looping through liquid water states in the snow domain)
   endif   ! (if the subset includes hydrology state variables in the snow domain)

   ! -----
   ! * liquid water fluxes for the soil domain...
   ! --------------------------------------------
   if(nSoilOnlyHyd>0)then
    do iLayer=1,nSoil

     ! - check that the soil layer is desired
     if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle

     ! - define state indices
     watState = ixSoilOnlyHyd(iLayer)         ! hydrology state index within the state subset

     ! - define indices of the soil layers
     jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

     ! - compute the diagonal elements
     ! all terms *excluding* baseflow
     aJac(ixDiag,watState) = (dt/mLayerDepth(jLayer))*(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(watState)

     ! - compute the lower-diagonal elements
     if(iLayer > 1)then
      if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(iLayer-1),watState),watState) = (dt/mLayerDepth(jLayer-1))*( dq_dHydStateBelow(iLayer-1))
     endif

     ! - compute the upper-diagonal elements
     if(iLayer<nSoil)then
      if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(iLayer+1),watState),watState) = (dt/mLayerDepth(jLayer+1))*(-dq_dHydStateAbove(iLayer))
     endif

    end do  ! (looping through hydrology states in the soil domain)
   endif   ! (if the subset includes hydrology state variables in the soil domain)

   ! -----
   ! * derivative in liquid water fluxes w.r.t. temperature for the soil domain...
   ! -----------------------------------------------------------------------------
   if(nSoilOnlyHyd>0 .and. nSoilOnlyNrg>0)then
    do iLayer=1,nSoilOnlyHyd

     ! - check that the soil layer is desired
     if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle

     ! - define index of hydrology state variable within the state subset
     watState = ixSoilOnlyHyd(iLayer)

     ! - define indices of the soil layers
     jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

     ! - define the energy state variable
     nrgState = ixNrgLayer(jLayer)       ! index within the full state vector

     ! only compute derivatives if the energy state for the current layer is within the state subset
     if(nrgstate/=integerMissing)then

      ! - compute the Jacobian for the layer itself
      aJac(ixOffDiag(watState,nrgState),nrgState) = (dt/mLayerDepth(jLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance

      ! - include derivatives w.r.t. ground evaporation
      if(nSnow==0 .and. iLayer==1)then  ! upper-most soil layer
       if(computeVegFlux)then
        aJac(ixOffDiag(watState,ixVegHyd),ixVegHyd) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dCanLiq/iden_water)  ! dVol/dLiq (kg m-2)-1
        aJac(ixOffDiag(watState,ixCasNrg),ixCasNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanair/iden_water) ! dVol/dT (K-1)
        aJac(ixOffDiag(watState,ixVegNrg),ixVegNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanopy/iden_water) ! dVol/dT (K-1)
       endif
       aJac(ixOffDiag(watState,ixTopNrg),ixTopNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTGround/iden_water) + aJac(ixOffDiag(watState,ixTopNrg),ixTopNrg) ! dVol/dT (K-1)
      endif

      ! melt-freeze: compute derivative in energy with respect to mass
      if(mLayerdTheta_dTk(jLayer) > verySmall)then  ! ice is present
       aJac(ixOffDiag(nrgState,watState),watState) = -dVolTot_dPsi0(iLayer)*LH_fus*iden_water    ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
      else
       aJac(ixOffDiag(nrgState,watState),watState) = 0._dp
      endif

      ! - compute lower diagonal elements
      if(iLayer>1)then
       if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(iLayer-1),nrgState),nrgState) = (dt/mLayerDepth(jLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
      endif

      ! compute upper-diagonal elements
      if(iLayer<nSoil)then
       if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixOffDiag(ixSoilOnlyHyd(iLayer+1),nrgState),nrgState) = (dt/mLayerDepth(jLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1
      endif

     endif   ! (if the energy state for the current layer is within the state subset)

    end do  ! (looping through soil layers)
   endif   ! (if there are state variables for both water and energy in the soil domain)
   
   if(globalPrintFlag)then
    print*, '** banded analytical Jacobian:'
    write(*,'(a4,1x,100(i17,1x))') 'xCol', (iLayer, iLayer=min(iJac1,nState),min(iJac2,nState))
    do iLayer=kl+1,nBands
     write(*,'(i4,1x,100(e17.10,1x))') iLayer, (aJac(iLayer,jLayer),jLayer=min(iJac1,nState),min(iJac2,nState))
    end do
   end if
   !print*, 'PAUSE: banded analytical Jacobian'; read(*,*)

  ! *********************************************************************************************************************************************************
  ! *********************************************************************************************************************************************************
  ! * PART 2: FULL MATRIX
  ! *********************************************************************************************************************************************************
  ! *********************************************************************************************************************************************************
  case(ixFullMatrix)

   ! check
   if(size(aJac,1)/=size(dMat) .or. size(aJac,2)/=size(dMat))then
    message=trim(message)//'unexpected shape of the Jacobian matrix: expect aJac(nState,nState)'
    err=20; return
   end if

   ! -----
   ! * energy and liquid fluxes over vegetation...
   ! ---------------------------------------------
   if(computeVegFlux)then  ! (derivatives only defined when vegetation protrudes over the surface)
  
    ! * liquid water fluxes for vegetation canopy (-)
    if(ixVegHyd/=integerMissing) aJac(ixVegHyd,ixVegHyd) = -scalarFracLiqVeg*(dCanopyEvaporation_dCanLiq - scalarCanopyLiqDeriv)*dt + 1._dp

    ! * cross-derivative terms for canopy water
    if(ixVegHyd/=integerMissing)then
     ! cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
     if(ixCasNrg/=integerMissing) aJac(ixVegHyd,ixCasNrg) = -dCanopyEvaporation_dTCanair*dt
     if(ixVegNrg/=integerMissing) aJac(ixVegHyd,ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy
     if(ixTopNrg/=integerMissing) aJac(ixVegHyd,ixTopNrg) = -dCanopyEvaporation_dTGround*dt
     ! cross-derivative terms w.r.t. canopy water (kg-1 m2)
     if(ixTopHyd/=integerMissing) aJac(ixTopHyd,ixVegHyd) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarFracLiqVeg*scalarCanopyLiqDeriv)/iden_water
     ! cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
     ! NOTE: dIce/dLiq = (1 - scalarFracLiqVeg); dIce*LH_fus/canopyDepth = J m-3; dLiq = kg m-2
     if(ixVegNrg/=integerMissing) aJac(ixVegNrg,ixVegHyd) = (dt/canopyDepth)   *(-dCanopyNetFlux_dCanLiq) - (1._dp - scalarFracLiqVeg)*LH_fus/canopyDepth   ! dF/dLiq
     if(ixTopNrg/=integerMissing) aJac(ixTopNrg,ixVegHyd) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanLiq)
    endif  

    ! cross-derivative terms w.r.t. canopy temperature (K-1)
    if(ixVegNrg/=integerMissing)then
     if(ixTopHyd/=integerMissing) aJac(ixTopHyd,ixVegNrg) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water
    endif  

    ! energy fluxes with the canopy air space (J m-3 K-1)
    if(ixCasNrg/=integerMissing)then
                                  aJac(ixCasNrg,ixCasNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + dMat(ixCasNrg)
     if(ixVegNrg/=integerMissing) aJac(ixCasNrg,ixVegNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)
     if(ixTopNrg/=integerMissing) aJac(ixCasNrg,ixTopNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dGroundTemp)
    endif 
 
    ! energy fluxes with the vegetation canopy (J m-3 K-1)
    if(ixVegNrg/=integerMissing)then
     if(ixCasNrg/=integerMissing) aJac(ixVegNrg,ixCasNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)
                                  aJac(ixVegNrg,ixVegNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + dMat(ixVegNrg)
     if(ixTopNrg/=integerMissing) aJac(ixVegNrg,ixTopNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)
    endif  

    ! energy fluxes with the surface (J m-3 K-1)
    if(ixTopNrg/=integerMissing)then
     if(ixCasNrg/=integerMissing) aJac(ixTopNrg,ixCasNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)
     if(ixVegNrg/=integerMissing) aJac(ixTopNrg,ixVegNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)
    endif  

   endif  ! if there is a need to compute energy fluxes within vegetation
  
   ! -----
   ! * energy fluxes for the snow+soil domain...
   ! -------------------------------------------
   if(nSnowSoilNrg>0)then
    do iLayer=1,nLayers  ! loop through all layers in the snow+soil domain 

     ! check if the state is in the subset
     if(ixSnowSoilNrg(iLayer)==integerMissing) cycle

     ! - define index within the state subset and the full state vector
     jState = ixSnowSoilNrg(iLayer)        ! index within the state subset

     ! - diagonal elements
     aJac(jState,jState)   = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(jState)

     ! - lower-diagonal elements
     if(iLayer > 1)then
      if(ixSnowSoilNrg(iLayer-1)/=integerMissing) aJac(ixSnowSoilNrg(iLayer-1),jState) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
     endif

     ! - upper diagonal elements
     if(iLayer < nLayers)then
      if(ixSnowSoilNrg(iLayer+1)/=integerMissing) aJac(ixSnowSoilNrg(iLayer+1),jState) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
     endif

    end do  ! (looping through energy states in the snow+soil domain)
   endif   ! (if the subset includes energy state variables in the snow+soil domain)
  
   ! -----
   ! * liquid water fluxes for the snow domain...
   ! --------------------------------------------
   if(nSnowOnlyHyd>0)then
    do iLayer=1,nSnow  ! loop through layers in the snow domain

     ! - check that the snow layer is desired
     if(ixSnowOnlyHyd(iLayer)==integerMissing) cycle

     ! - define state indices for the current layer
     watState = ixSnowOnlyHyd(iLayer)   ! hydrology state index within the state subset

     ! compute factor to convert liquid water derivative to total water derivative
     select case( ixHydType(iLayer) )
      case(iname_watLayer); convLiq2tot = mLayerFracLiqSnow(iLayer)
      case default;         convLiq2tot = 1._dp
     end select

     ! - diagonal elements
     aJac(watState,watState) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot + dMat(watState)

     ! - lower-diagonal elements
     if(iLayer > 1)then
      if(ixSnowOnlyHyd(iLayer-1)/=integerMissing) aJac(ixSnowOnlyHyd(iLayer-1),watState) = 0._dp  ! sub-diagonal: no dependence on other layers
     endif

     ! - upper diagonal elements
     if(iLayer < nSnow)then
      if(ixSnowOnlyHyd(iLayer+1)/=integerMissing) aJac(ixSnowOnlyHyd(iLayer+1),watState) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot       ! dVol(below)/dLiq(above) -- (-)
     endif

     ! - compute cross-derivative terms for energy
     ! NOTE: increase in volumetric liquid water content balanced by a decrease in volumetric ice content
     if(nSnowOnlyNrg>0)then

      ! (define the energy state)
      nrgState = ixSnowOnlyNrg(iLayer)       ! index within the full state vector
      if(nrgstate/=integerMissing)then       ! (energy state for the current layer is within the state subset)

       ! (cross-derivative terms for the current layer)
       aJac(nrgState,watState) = -(1._dp - mLayerFracLiqSnow(iLayer))*LH_fus*iden_water     ! (dF/dLiq)
       aJac(watState,nrgState) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! (dVol/dT)

       ! (cross-derivative terms for the layer below)
       if(iLayer < nSnow)then
        aJac(ixSnowOnlyHyd(iLayer+1),nrgState) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)        ! dVol(below)/dT(above) -- K-1
       endif ! (if there is a water state in the layer below the current layer in the given state subset)

      endif ! (if the energy state for the current layer is within the state subset)
     endif ! (if state variables exist for energy in snow+soil layers)

    end do  ! (looping through liquid water states in the snow domain)
   endif   ! (if the subset includes hydrology state variables in the snow domain)
   
   ! -----
   ! * liquid water fluxes for the soil domain...
   ! --------------------------------------------
   if(nSoilOnlyHyd>0)then
    do iLayer=1,nSoil

     ! - check that the soil layer is desired
     if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle

     ! - define state indices
     watState = ixSoilOnlyHyd(iLayer)         ! hydrology state index within the state subset

     ! - define indices of the soil layers
     jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

     ! - compute the diagonal elements
     ! all terms *excluding* baseflow
     aJac(watState,watState) = (dt/mLayerDepth(jLayer))*(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(watState)

     ! - compute the lower-diagonal elements
     if(iLayer > 1)then
      if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer-1),watState) = (dt/mLayerDepth(jLayer-1))*( dq_dHydStateBelow(iLayer-1))
     endif

     ! - compute the upper-diagonal elements
     if(iLayer<nSoil)then
      if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer+1),watState) = (dt/mLayerDepth(jLayer+1))*(-dq_dHydStateAbove(iLayer))
     endif
   
     ! - include terms for baseflow
     if(computeBaseflow)then
      do pLayer=1,nSoil
       qState = ixSoilOnlyHyd(pLayer)  ! hydrology state index within the state subset
       aJac(watState,qState) = aJac(watState,qState) + (dt/mLayerDepth(jLayer))*dBaseflow_dMatric(iLayer,pLayer)
      end do
     endif
   
    end do  ! (looping through hydrology states in the soil domain)
   endif   ! (if the subset includes hydrology state variables in the soil domain)

   ! -----
   ! * derivative in liquid water fluxes w.r.t. temperature for the soil domain...
   ! -----------------------------------------------------------------------------
   if(nSoilOnlyHyd>0 .and. nSoilOnlyNrg>0)then
    do iLayer=1,nSoilOnlyHyd

     ! - check that the soil layer is desired
     if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle

     ! - define index of hydrology state variable within the state subset
     watState = ixSoilOnlyHyd(iLayer)

     ! - define indices of the soil layers
     jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

     ! - define the energy state variable
     nrgState = ixNrgLayer(jLayer)       ! index within the full state vector
     
     ! only compute derivatives if the energy state for the current layer is within the state subset
     if(nrgstate/=integerMissing)then

      ! - compute the Jacobian for the layer itself
      aJac(watState,nrgState) = (dt/mLayerDepth(jLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance
    
      ! - include derivatives w.r.t. ground evaporation
      if(nSnow==0 .and. iLayer==1)then  ! upper-most soil layer
       if(computeVegFlux)then
        aJac(watState,ixVegHyd) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dCanLiq/iden_water)  ! dVol/dLiq (kg m-2)-1
        aJac(watState,ixCasNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanair/iden_water) ! dVol/dT (K-1)
        aJac(watState,ixVegNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanopy/iden_water) ! dVol/dT (K-1)
       endif
       aJac(watState,ixTopNrg) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTGround/iden_water) + aJac(watState,ixTopNrg) ! dVol/dT (K-1)
      endif
    
      ! melt-freeze: compute derivative in energy with respect to mass
      if(mLayerdTheta_dTk(jLayer) > verySmall)then  ! ice is present
       aJac(nrgState,watState) = -dVolTot_dPsi0(iLayer)*LH_fus*iden_water    ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
      else
       aJac(nrgState,watState) = 0._dp
      endif
    
      ! - compute lower diagonal elements
      if(iLayer>1)then
       if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer-1),nrgState) = (dt/mLayerDepth(jLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
      endif

      ! compute upper-diagonal elements
      if(iLayer<nSoil)then
       if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac(ixSoilOnlyHyd(iLayer+1),nrgState) = (dt/mLayerDepth(jLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1
      endif

     endif   ! (if the energy state for the current layer is within the state subset) 
 
    end do  ! (looping through soil layers)
   endif   ! (if there are state variables for both water and energy in the soil domain)
  
   ! print the Jacobian
   if(globalPrintFlag)then
    print*, '** analytical Jacobian:'
    write(*,'(a4,1x,100(i12,1x))') 'xCol', (iLayer, iLayer=min(iJac1,nState),min(iJac2,nState))
    do iLayer=iJac1,iJac2; write(*,'(i4,1x,100(e12.5,1x))') iLayer, aJac(min(iJac1,nState):min(iJac2,nState),iLayer); end do
   end if

  ! ***
  ! check
  case default; err=20; message=trim(message)//'unable to identify option for the type of matrix'; return

 end select  ! type of matrix

 ! end association to variables in the data structures
 end associate

 end subroutine computJacob

 
 ! **********************************************************************************************************
 ! private function: get the off-diagonal index in the band-diagonal matrix
 ! **********************************************************************************************************
 function ixOffDiag(jState,iState)
 implicit none
 integer(i4b),intent(in)  :: jState    ! off-diagonal state
 integer(i4b),intent(in)  :: iState    ! diagonal state
 integer(i4b)             :: ixOffDiag ! off-diagonal index in gthe band-diagonal matrix
 ixOffDiag = ixBandOffset + jState - iState
 end function ixOffDiag

end module computJacob_module
