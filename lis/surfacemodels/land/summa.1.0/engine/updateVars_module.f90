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

module updateVars_module

! data types
USE nrtype

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_cas       ! named variables for canopy air space
USE globalData,only:iname_veg       ! named variables for vegetation canopy
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! metadata for information in the data structures
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! constants
USE multiconst,only:&
                    gravity,      & ! acceleration of gravity              (m s-2)
                    Tfreeze,      & ! temperature at freezing              (K)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (dp)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (dp)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookDIAG             ! named variables for structure elements
USE var_lookup,only:iLookPROG             ! named variables for structure elements
USE var_lookup,only:iLookDERIV            ! named variables for structure elements
USE var_lookup,only:iLookPARAM            ! named variables for structure elements
USE var_lookup,only:iLookINDEX            ! named variables for structure elements

! provide access to routines to update states
USE updatState_module,only:updateSnow     ! update snow states
USE updatState_module,only:updateSoil     ! update soil states

! provide access to functions for the constitutive functions and derivatives
USE snow_utils_module,only:fracliquid     ! compute the fraction of liquid water (snow) 
USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
USE soil_utils_module,only:dTheta_dTk     ! differentiate the freezing curve w.r.t. temperature (soil)
USE soil_utils_module,only:dTheta_dPsi    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:dPsi_dTheta    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric water content
USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water
USE soil_utils_module,only:crit_soilT     ! compute critical temperature below which ice exists
USE soil_utils_module,only:liquidHead     ! compute the liquid water matric potential 

! IEEE checks
USE, intrinsic :: ieee_arithmetic            ! check values (NaN, etc.)

implicit none
private
public::updateVars

contains

 ! **********************************************************************************************************
 ! public subroutine updateVars: compute diagnostic variables
 ! **********************************************************************************************************
 subroutine updateVars(&
                       ! input
                       do_adjustTemp,                             & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                       explicitEuler,                             & ! intent(in):    flag to denote computing the explicit Euler solution
                       mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,                                 & ! intent(inout): model diagnostic variables for a local HRU
                       deriv_data,                                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       ! output: variables for the vegetation canopy
                       scalarCanopyTempTrial,                     & ! intent(inout): trial value of canopy temperature (K)
                       scalarCanopyWatTrial,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                       scalarCanopyLiqTrial,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                       scalarCanopyIceTrial,                      & ! intent(inout): trial value of canopy ice content (kg m-2)
                       ! output: variables for the snow-soil domain
                       mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                       mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-) 
                       mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-) 
                       mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-) 
                       mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                       mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                       ! output: error control 
                       err,message)                                 ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 implicit none 
 ! input
 logical(lgt)     ,intent(in)    :: do_adjustTemp                   ! flag to adjust temperature to account for the energy used in melt+freeze
 logical(lgt)     ,intent(in)    :: explicitEuler                   ! flag to denote computing the explicit Euler solution
 type(var_dlength),intent(in)    :: mpar_data                       ! model parameters for a local HRU
 type(var_ilength),intent(in)    :: indx_data                       ! indices defining model states and layers                 
 type(var_dlength),intent(in)    :: prog_data                       ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data                       ! diagnostic variables for a local HRU
 type(var_dlength),intent(inout) :: deriv_data                      ! derivatives in model fluxes w.r.t. relevant state variables
 ! output: variables for the vegetation canopy
 real(dp),intent(inout)          :: scalarCanopyTempTrial           ! trial value of canopy temperature (K)
 real(dp),intent(inout)          :: scalarCanopyWatTrial            ! trial value of canopy total water (kg m-2)
 real(dp),intent(inout)          :: scalarCanopyLiqTrial            ! trial value of canopy liquid water (kg m-2)
 real(dp),intent(inout)          :: scalarCanopyIceTrial            ! trial value of canopy ice content (kg m-2)
 ! output: variables for the snow-soil domain
 real(dp),intent(inout)          :: mLayerTempTrial(:)              ! trial vector of layer temperature (K)
 real(dp),intent(inout)          :: mLayerVolFracWatTrial(:)        ! trial vector of volumetric total water content (-)
 real(dp),intent(inout)          :: mLayerVolFracLiqTrial(:)        ! trial vector of volumetric liquid water content (-)
 real(dp),intent(inout)          :: mLayerVolFracIceTrial(:)        ! trial vector of volumetric ice water content (-)
 real(dp),intent(inout)          :: mLayerMatricHeadTrial(:)        ! trial vector of total water matric potential (m)
 real(dp),intent(inout)          :: mLayerMatricHeadLiqTrial(:)     ! trial vector of liquid water matric potential (m)
 ! output: error control 
 integer(i4b),intent(out)        :: err                             ! error code
 character(*),intent(out)        :: message                         ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! general local variables
 integer(i4b)                    :: iState                          ! index of model state variable
 integer(i4b)                    :: iLayer                          ! index of layer within the snow+soil domain
 integer(i4b)                    :: ixFullVector                    ! index within full state vector
 integer(i4b)                    :: ixDomainType                    ! name of a given model domain
 integer(i4b)                    :: ixControlIndex                  ! index within a given model domain
 integer(i4b)                    :: ixOther,ixOtherLocal            ! index of the coupled state variable within the (full, local) vector
 logical(lgt)                    :: isCoupled                       ! .true. if a given variable shared another state variable in the same control volume
 logical(lgt)                    :: isNrgState                      ! .true. if a given variable is an energy state 
 logical(lgt),allocatable        :: computedCoupling(:)             ! .true. if computed the coupling for a given state variable
 real(dp)                        :: scalarVolFracLiq                ! volumetric fraction of liquid water (-)
 real(dp)                        :: scalarVolFracIce                ! volumetric fraction of ice (-)
 real(dp)                        :: Tcrit                           ! critical soil temperature below which ice exists (K)
 real(dp)                        :: xTemp                           ! temporary temperature (K)
 real(dp)                        :: effSat                          ! effective saturation (-)
 real(dp)                        :: avPore                          ! available pore space (-)
 character(len=256)              :: cMessage                        ! error message of downwind routine
 logical(lgt),parameter          :: printFlag=.false.               ! flag to turn on printing
 ! iterative solution for temperature
 real(dp)                        :: meltNrg                         ! energy for melt+freeze (J m-3)
 real(dp)                        :: residual                        ! residual in the energy equation (J m-3)
 real(dp)                        :: derivative                      ! derivative in the energy equation (J m-3 K-1)
 real(dp)                        :: tempInc                         ! iteration increment (K)
 integer(i4b)                    :: iter                            ! iteration index
 integer(i4b)                    :: niter                           ! number of iterations
 integer(i4b),parameter          :: maxiter=100                     ! maximum number of iterations
 real(dp),parameter              :: nrgConvTol=1.e-4_dp             ! convergence tolerance for energy (J m-3)
 real(dp),parameter              :: tempConvTol=1.e-6_dp            ! convergence tolerance for temperature (K)
 real(dp)                        :: critDiff                        ! temperature difference from critical (K)
 real(dp)                        :: tempMin                         ! minimum bracket for temperature (K)
 real(dp)                        :: tempMax                         ! maximum bracket for temperature (K)
 logical(lgt)                    :: bFlag                           ! flag to denote that iteration increment was constrained using bi-section
 real(dp),parameter              :: epsT=1.e-7_dp                   ! small interval above/below critical temperature (K)
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of model layers, and layer type
 nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
 nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
 nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
 mLayerDepth             => prog_data%var(iLookPROG%mLayerDepth)%dat               ,& ! intent(in):  [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
 ! indices defining model states and layers
 ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy energy state variable
 ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
 ! indices in the full vector for specific domains
 ixNrgCanair             => indx_data%var(iLookINDEX%ixNrgCanair)%dat              ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in canopy air space domain
 ixNrgCanopy             => indx_data%var(iLookINDEX%ixNrgCanopy)%dat              ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
 ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat              ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
 ixNrgLayer              => indx_data%var(iLookINDEX%ixNrgLayer)%dat               ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer              => indx_data%var(iLookINDEX%ixHydLayer)%dat               ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ! mapping between the full state vector and the state subset
 ixMapFull2Subset        => indx_data%var(iLookINDEX%ixMapFull2Subset)%dat         ,& ! intent(in):  [i4b(:)] list of indices in the state subset for each state in the full state vector
 ixMapSubset2Full        => indx_data%var(iLookINDEX%ixMapSubset2Full)%dat         ,& ! intent(in):  [i4b(:)] [state subset] list of indices of the full state vector in the state subset
 ! type of domain, type of state variable, and index of control volume within domain
 ixDomainType_subset     => indx_data%var(iLookINDEX%ixDomainType_subset)%dat      ,& ! intent(in):  [i4b(:)] [state subset] id of domain for desired model state variables
 ixControlVolume         => indx_data%var(iLookINDEX%ixControlVolume)%dat          ,& ! intent(in):  [i4b(:)] index of the control volume for different domains (veg, snow, soil)
 ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat              ,& ! intent(in):  [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
 ! snow parameters
 snowfrz_scale           => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)         ,& ! intent(in):  [dp] scaling parameter for the snow freezing curve (K-1)
 ! depth-varying model parameters
 vGn_m                   => diag_data%var(iLookDIAG%scalarVGn_m)%dat               ,& ! intent(in):  [dp(:)] van Genutchen "m" parameter (-)
 vGn_n                   => mpar_data%var(iLookPARAM%vGn_n)%dat                    ,& ! intent(in):  [dp(:)] van Genutchen "n" parameter (-)
 vGn_alpha               => mpar_data%var(iLookPARAM%vGn_alpha)%dat                ,& ! intent(in):  [dp(:)] van Genutchen "alpha" parameter (m-1)
 theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat                ,& ! intent(in):  [dp(:)] soil porosity (-)
 theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat                ,& ! intent(in):  [dp(:)] soil residual volumetric water content (-)
 ! model diagnostic variables (heat capacity)
 canopyDepth             => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)      ,& ! intent(in):  [dp   ] canopy depth (m)
 scalarBulkVolHeatCapVeg => diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1),& ! intent(in):  [dp   ] volumetric heat capacity of the vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat        ,& ! intent(in):  [dp(:)] volumetric heat capacity in each layer (J m-3 K-1)
 ! model diagnostic variables (fraction of liquid water)
 scalarFracLiqVeg        => diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)       ,& ! intent(out): [dp]    fraction of liquid water on vegetation (-)
 mLayerFracLiqSnow       => diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat         ,& ! intent(out): [dp(:)] fraction of liquid water in each snow layer (-)
 ! model states for the vegetation canopy
 scalarCanairTemp        => prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)       ,& ! intent(in):  [dp] temperature of the canopy air space (K)
 scalarCanopyTemp        => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)       ,& ! intent(in):  [dp] temperature of the vegetation canopy (K)
 scalarCanopyWat         => prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)        ,& ! intent(in):  [dp] mass of total water on the vegetation canopy (kg m-2)
 ! model state variable vectors for the snow-soil layers
 mLayerTemp              => prog_data%var(iLookPROG%mLayerTemp)%dat                ,& ! intent(in):  [dp(:)] temperature of each snow/soil layer (K)
 mLayerVolFracWat        => prog_data%var(iLookPROG%mLayerVolFracWat)%dat          ,& ! intent(in):  [dp(:)] volumetric fraction of total water (-)
 mLayerMatricHead        => prog_data%var(iLookPROG%mLayerMatricHead)%dat          ,& ! intent(in):  [dp(:)] total water matric potential (m)
 mLayerMatricHeadLiq     => diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat       ,& ! intent(in):  [dp(:)] liquid water matric potential (m)
 ! model diagnostic variables from a previous solution
 scalarCanopyLiq         => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)        ,& ! intent(in):  [dp(:)] mass of liquid water on the vegetation canopy (kg m-2)
 scalarCanopyIce         => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)        ,& ! intent(in):  [dp(:)] mass of ice on the vegetation canopy (kg m-2)
 mLayerVolFracLiq        => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat          ,& ! intent(in):  [dp(:)] volumetric fraction of liquid water (-)
 mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat          ,& ! intent(in):  [dp(:)] volumetric fraction of ice (-)
 ! derivatives
 dVolTot_dPsi0           => deriv_data%var(iLookDERIV%dVolTot_dPsi0   )%dat        ,& ! intent(out): [dp(:)] derivative in total water content w.r.t. total water matric potential
 dPsiLiq_dPsi0           => deriv_data%var(iLookDERIV%dPsiLiq_dPsi0   )%dat        ,& ! intent(out): [dp(:)] derivative in liquid water matric pot w.r.t. the total water matric pot (-)
 dPsiLiq_dTemp           => deriv_data%var(iLookDERIV%dPsiLiq_dTemp   )%dat        ,& ! intent(out): [dp(:)] derivative in the liquid water matric potential w.r.t. temperature
 mLayerdTheta_dTk        => deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat        ,& ! intent(out): [dp(:)] derivative of volumetric liquid water content w.r.t. temperature
 dTheta_dTkCanopy        => deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)      & ! intent(out): [dp]    derivative of volumetric liquid water content w.r.t. temperature
 ) ! association with variables in the data structures

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! initialize error control
 err=0; message='updateVars/'

 ! allocate space and assign values to the flag vector
 allocate(computedCoupling(size(ixMapSubset2Full)),stat=err)        ! .true. if computed the coupling for a given state variable
 if(err/=0)then; message=trim(message)//'problem allocating computedCoupling'; return; endif
 computedCoupling(:)=.false.

 ! loop through model state variables
 do iState=1,size(ixMapSubset2Full)

  ! check the need for the computations
  if(computedCoupling(iState)) cycle

  ! -----
  ! - compute indices...
  ! --------------------

  ! get domain type, and index of the control volume within the domain
  ixFullVector   = ixMapSubset2Full(iState)       ! index within full state vector
  ixDomainType   = ixDomainType_subset(iState)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
  ixControlIndex = ixControlVolume(ixFullVector)  ! index within a given domain

  ! get the layer index
  select case(ixDomainType)
   case(iname_cas);  cycle ! canopy air space: do nothing
   case(iname_veg);  iLayer = 0 
   case(iname_snow); iLayer = ixControlIndex
   case(iname_soil); iLayer = ixControlIndex + nSnow
   case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil'; return
  end select
   
  ! get the index of the other (energy or mass) state variable within the full state vector
  select case(ixDomainType)
   case(iname_veg)             ; ixOther = merge(ixHydCanopy(1),    ixNrgCanopy(1),    ixStateType(ixFullVector)==iname_nrgCanopy)
   case(iname_snow, iname_soil); ixOther = merge(ixHydLayer(iLayer),ixNrgLayer(iLayer),ixStateType(ixFullVector)==iname_nrgLayer)
   case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
  end select

  ! get the index in the local state vector
  ixOtherLocal = ixMapFull2Subset(ixOther)  ! ixOtherLocal could equal integerMissing
  if(ixOtherLocal/=integerMissing) computedCoupling(ixOtherLocal)=.true.

  ! check if we have a coupled solution
  isCoupled    = (ixOtherLocal/=integerMissing)

  ! check if we are an energy state
  isNrgState   = (ixStateType(ixFullVector)==iname_nrgCanopy .or. ixStateType(ixFullVector)==iname_nrgLayer)

  if(printFlag)then
   print*, 'iState         = ', iState, size(ixMapSubset2Full) 
   print*, 'ixFullVector   = ', ixFullVector
   print*, 'ixDomainType   = ', ixDomainType
   print*, 'ixControlIndex = ', ixControlIndex
   print*, 'ixOther        = ', ixOther
   print*, 'ixOtherLocal   = ', ixOtherLocal
   print*, 'do_adjustTemp  = ', do_adjustTemp
   print*, 'isCoupled      = ', isCoupled
   print*, 'isNrgState     = ', isNrgState
  endif

  ! =======================================================================================================================================
  ! =======================================================================================================================================
  ! =======================================================================================================================================
  ! =======================================================================================================================================
  ! =======================================================================================================================================
  ! =======================================================================================================================================

  ! update hydrology state variables for the uncoupled solution
  if(.not.isNrgState .and. .not.isCoupled)then

   ! update the total water from volumetric liquid water
   if(ixStateType(ixFullVector)==iname_liqCanopy .or. ixStateType(ixFullVector)==iname_liqLayer)then
    select case(ixDomainType)
     case(iname_veg);    scalarCanopyWatTrial          = scalarCanopyLiqTrial          + scalarCanopyIceTrial
     case(iname_snow);   mLayerVolFracWatTrial(iLayer) = mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer)*iden_ice/iden_water
     case(iname_soil);   mLayerVolFracWatTrial(iLayer) = mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer) ! no volume expansion
     case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, or iname_soil'; return
    end select
   endif

   ! update the total water and the total water matric potential
   if(ixDomainType==iname_soil)then
    select case( ixStateType(ixFullVector) )
     ! --> update the total water from the liquid water matric potential
     case(iname_lmpLayer)
      effSat = volFracLiq(mLayerMatricHeadLiqTrial(ixControlIndex),vGn_alpha(ixControlIndex),0._dp,1._dp,vGn_n(ixControlIndex),vGn_m(ixControlIndex))  ! effective saturation
      avPore = theta_sat(ixControlIndex) - mLayerVolFracIceTrial(iLayer) - theta_res(ixControlIndex)  ! available pore space
      mLayerVolFracLiqTrial(iLayer) = effSat*avPore + theta_res(ixControlIndex)
      mLayerVolFracWatTrial(iLayer) = mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer) ! no volume expansion
      mLayerMatricHeadTrial(ixControlIndex) = matricHead(mLayerVolFracWatTrial(iLayer),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex)) 
      !write(*,'(a,1x,i4,1x,3(f20.10,1x))') 'mLayerVolFracLiqTrial(iLayer) 1 = ', iLayer, mLayerVolFracLiqTrial(iLayer), mLayerVolFracIceTrial(iLayer), mLayerVolFracWatTrial(iLayer)
     ! --> update the total water from the total water matric potential
     case(iname_matLayer)
      mLayerVolFracWatTrial(iLayer) = volFracLiq(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
     ! --> update the total water matric potential (assume already have mLayerVolFracWatTrial given block above)
     case(iname_liqLayer, iname_watLayer)
      mLayerMatricHeadTrial(ixControlIndex) = matricHead(mLayerVolFracWatTrial(iLayer),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex)) 
     case default; err=20; message=trim(message)//'expect iname_lmpLayer, iname_matLayer, iname_liqLayer, or iname_watLayer'; return
    end select
   endif  ! if in the soil domain

  endif  ! if hydrology state variable or uncoupled solution

  ! compute the critical soil temperature below which ice exists
  select case(ixDomainType)
   case(iname_veg, iname_snow); Tcrit = Tfreeze
   case(iname_soil);            Tcrit = crit_soilT( mLayerMatricHeadTrial(ixControlIndex) )
   case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
  end select
  
  ! initialize temperature 
  select case(ixDomainType)
   case(iname_veg);              xTemp = scalarCanopyTempTrial
   case(iname_snow, iname_soil); xTemp = mLayerTempTrial(iLayer)
   case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
  end select

  ! define brackets for the root
  ! NOTE: start with an enormous range; updated quickly in the iterations
  tempMin = xTemp - 10._dp
  tempMax = xTemp + 10._dp

  ! get iterations (set to maximum iterations if adjusting the temperature)
  niter = merge(maxiter, 1, do_adjustTemp)

  ! iterate
  iterations: do iter=1,niter

   ! restrict temperature
   if(xTemp <= tempMin .or. xTemp >= tempMax)then
    xTemp = 0.5_dp*(tempMin + tempMax)  ! new value
    bFlag = .true.
   else
    bFlag = .false.
   endif

   ! -----
   ! - compute derivatives...
   ! ------------------------
  
   ! compute the derivative in total water content w.r.t. total water matric potential (m-1)
   ! NOTE 1: valid for frozen and unfrozen conditions
   ! NOTE 2: for case "iname_lmpLayer", dVolTot_dPsi0 = dVolLiq_dPsi
   if(ixDomainType==iname_soil)then
    select case( ixStateType(ixFullVector) )
     case(iname_lmpLayer); dVolTot_dPsi0(ixControlIndex) = dTheta_dPsi(mLayerMatricHeadLiqTrial(ixControlIndex),vGn_alpha(ixControlIndex),0._dp,1._dp,vGn_n(ixControlIndex),vGn_m(ixControlIndex))*avPore
     case default;         dVolTot_dPsi0(ixControlIndex) = dTheta_dPsi(mLayerMatricHeadTrial(ixControlIndex),vGn_alpha(ixControlIndex),theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
    end select
   endif

   ! compute the derivative in liquid water content w.r.t. temperature
   ! --> partially frozen: dependence of liquid water on temperature
   if(xTemp<Tcrit)then
    select case(ixDomainType)
     case(iname_veg);  dTheta_dTkCanopy         = dFracLiq_dTk(xTemp,snowfrz_scale)*scalarCanopyWat/(iden_water*canopyDepth)
     case(iname_snow); mLayerdTheta_dTk(iLayer) = dFracLiq_dTk(xTemp,snowfrz_scale)*mLayerVolFracWatTrial(iLayer)
     case(iname_soil); mLayerdTheta_dTk(iLayer) = dTheta_dTk(xTemp,theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_alpha(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))
     case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
    end select  ! domain type

   ! --> unfrozen: no dependence of liquid water on temperature
   else
    select case(ixDomainType)
     case(iname_veg);              dTheta_dTkCanopy         = 0._dp 
     case(iname_snow, iname_soil); mLayerdTheta_dTk(iLayer) = 0._dp
     case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
    end select  ! domain type
   endif
  
   ! -----
   ! - update volumetric fraction of liquid water and ice...
   !    => case of hydrology state uncoupled with energy (and when not adjusting the temperature)...
   ! -----------------------------------------------------------------------------------------------
  
   ! case of hydrology state uncoupled with energy (and when not adjusting the temperature)
   if(.not.do_adjustTemp .and. .not.isNrgState .and. .not.isCoupled)then
  
    ! compute the fraction of snow
    select case(ixDomainType)
     case(iname_veg);  scalarFracLiqVeg          = fracliquid(xTemp,snowfrz_scale) 
     case(iname_snow); mLayerFracLiqSnow(iLayer) = fracliquid(xTemp,snowfrz_scale)
     case(iname_soil)  ! do nothing 
     case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
    end select  ! domain type
 
   ! -----
   ! - update volumetric fraction of liquid water and ice...
   !    => case of energy state or coupled solution (or adjusting the temperature)... 
   ! --------------------------------------------------------------------------------
  
   ! case of energy state OR coupled solution (or adjusting the temperature)
   ! NOTE: do not go in here if we have the explicit Euler solution for energy state variables (isNrgState or isCoupled)
   !         --> the explicit Euler updates for ice are done separately, based on the energy partitioning into melt and temperature at the start of the step
   elseif(do_adjustTemp .or. ( (isNrgState .or. isCoupled) .and. .not.explicitEuler ) )then

    ! identify domain type
    select case(ixDomainType)
   
     ! *** vegetation canopy
     case(iname_veg)
   
      ! compute volumetric fraction of liquid water and ice
      call updateSnow(xTemp,                                        & ! intent(in)   : temperature (K)
                      scalarCanopyWatTrial/(iden_water*canopyDepth),& ! intent(in)   : volumetric fraction of total water (-)
                      snowfrz_scale,                                & ! intent(in)   : scaling parameter for the snow freezing curve (K-1)
                      scalarVolFracLiq,                             & ! intent(out)  : trial volumetric fraction of liquid water (-)
                      scalarVolFracIce,                             & ! intent(out)  : trial volumetric fraction if ice (-)
                      scalarFracLiqVeg,                             & ! intent(out)  : fraction of liquid water (-)
                      err,cmessage)                                   ! intent(out)  : error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   
      ! compute mass of water on the canopy
      ! NOTE: possibilities for speed-up here
      scalarCanopyLiqTrial =          scalarFracLiqVeg *scalarCanopyWatTrial 
      scalarCanopyIceTrial = (1._dp - scalarFracLiqVeg)*scalarCanopyWatTrial
      
     ! *** snow layers
     case(iname_snow)

      ! compute volumetric fraction of liquid water and ice
      call updateSnow(xTemp,                                        & ! intent(in)   : temperature (K)
                      mLayerVolFracWatTrial(iLayer),                & ! intent(in)   : mass state variable = trial volumetric fraction of water (-)
                      snowfrz_scale,                                & ! intent(in)   : scaling parameter for the snow freezing curve (K-1)
                      mLayerVolFracLiqTrial(iLayer),                & ! intent(out)  : trial volumetric fraction of liquid water (-)
                      mLayerVolFracIceTrial(iLayer),                & ! intent(out)  : trial volumetric fraction if ice (-)
                      mLayerFracLiqSnow(iLayer),                    & ! intent(out)  : fraction of liquid water (-)
                      err,cmessage)                                   ! intent(out)  : error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   
     ! *** soil layers
     case(iname_soil)
   
      ! compute volumetric fraction of liquid water and ice
      call updateSoil(xTemp,                                                                                                                        & ! intent(in)   : temperature (K)
                      mLayerMatricHeadTrial(ixControlIndex),                                                                                        & ! intent(in)   : total water matric potential (m)
                      vGn_alpha(ixControlIndex),vGn_n(ixControlIndex),theta_sat(ixControlIndex),theta_res(ixControlIndex),vGn_m(ixControlIndex),    & ! intent(in)   : soil parameters
                      mLayerVolFracWatTrial(iLayer),                                                                                                & ! intent(in)   : mass state variable = trial volumetric fraction of water (-)
                      mLayerVolFracLiqTrial(iLayer),                                                                                                & ! intent(out)  : trial volumetric fraction of liquid water (-)
                      mLayerVolFracIceTrial(iLayer),                                                                                                & ! intent(out)  : trial volumetric fraction if ice (-)
                      err,cmessage)                                                                                                                   ! intent(out)  : error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! check
     case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return
   
    end select  ! domain type
  
   ! explicit Euler solution where energy state variables exist
   else

    ! do nothing (input = output) -- but check that we got here correctly
    if(explicitEuler .and. (isNrgState .or. isCoupled) )then
     scalarVolFracLiq = realMissing
     scalarVolFracIce = realMissing
    else
     message=trim(message)//'unexpected else branch: expect explicit Euler solution where energy state variables exist'
     err=20; return
    endif

   endif  ! if energy state or solution is coupled 

   ! -----
   ! - update temperatures...
   ! ------------------------

   ! check the need to adjust temperature
   if(do_adjustTemp)then

    ! get the melt energy
    meltNrg = merge(LH_fus*iden_ice, LH_fus*iden_water, ixDomainType==iname_snow)

    ! compute the residual and the derivative
    select case(ixDomainType)

     ! * vegetation
     case(iname_veg)
      call xTempSolve(&
                      ! constant over iterations
                      meltNrg         = meltNrg                                 ,&  ! intent(in)    : energy for melt+freeze (J m-3)
                      heatCap         = scalarBulkVolHeatCapVeg                 ,&  ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                      tempInit        = scalarCanopyTemp                        ,&  ! intent(in)    : initial temperature (K)
                      volFracIceInit  = scalarCanopyIce/(iden_water*canopyDepth),&  ! intent(in)    : initial volumetric fraction of ice (-)
                      ! trial values
                      xTemp           = xTemp                                   ,&  ! intent(inout) : trial value of temperature
                      dLiq_dT         = dTheta_dTkCanopy                        ,&  ! intent(in)    : derivative in liquid water content w.r.t. temperature (K-1)
                      volFracIceTrial = scalarVolFracIce                        ,&  ! intent(in)    : trial value for volumetric fraction of ice
                      ! residual and derivative
                      residual        = residual                                ,&  ! intent(out)   : residual (J m-3)
                      derivative      = derivative                               )  ! intent(out)   : derivative (J m-3 K-1)

     ! * snow and soil
     case(iname_snow, iname_soil)
      call xTempSolve(&
                      ! constant over iterations
                      meltNrg         = meltNrg                                 ,&  ! intent(in)    : energy for melt+freeze (J m-3)
                      heatCap         = mLayerVolHtCapBulk(iLayer)              ,&  ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                      tempInit        = mLayerTemp(iLayer)                      ,&  ! intent(in)    : initial temperature (K)
                      volFracIceInit  = mLayerVolFracIce(iLayer)                ,&  ! intent(in)    : initial volumetric fraction of ice (-)
                      ! trial values
                      xTemp           = xTemp                                   ,&  ! intent(inout) : trial value of temperature
                      dLiq_dT         = mLayerdTheta_dTk(iLayer)                ,&  ! intent(in)    : derivative in liquid water content w.r.t. temperature (K-1)
                      volFracIceTrial = mLayerVolFracIceTrial(iLayer)           ,&  ! intent(in)    : trial value for volumetric fraction of ice
                      ! residual and derivative
                      residual        = residual                                ,&  ! intent(out)   : residual (J m-3)
                      derivative      = derivative                               )  ! intent(out)   : derivative (J m-3 K-1)

     ! * check
     case default; err=20; message=trim(message)//'expect case to be iname_veg, iname_snow, iname_soil'; return

    end select  ! domain type

    ! check validity of residual
    if( ieee_is_nan(residual) )then
     message=trim(message)//'residual is not valid'
     err=20; return
    endif

    ! update bracket
    if(residual < 0._dp)then
     tempMax = min(xTemp,tempMax)
    else
     tempMin = max(tempMin,xTemp)
    end if

    ! compute iteration increment
    tempInc    = residual/derivative  ! K

    ! check
    if(globalPrintFlag)& 
    write(*,'(i4,1x,e20.10,1x,5(f20.10,1x),L1)') iter, residual, xTemp-Tcrit, tempInc, Tcrit, tempMin, tempMax, bFlag
 
    ! check convergence
    if(abs(residual) < nrgConvTol .or. abs(tempInc) < tempConvTol) exit iterations

    ! add constraints for snow temperature
    if(ixDomainType==iname_veg .or. ixDomainType==iname_snow)then
     if(tempInc > Tcrit - xTemp) tempInc=(Tcrit - xTemp)*0.5_dp  ! simple bi-section method
    endif  ! if the domain is vegetation or snow  

    ! deal with the discontinuity between partially frozen and unfrozen soil
    if(ixDomainType==iname_soil)then
     ! difference from the temperature below which ice exists
     critDiff = Tcrit - xTemp
     ! --> initially frozen (T < Tcrit)
     if(critDiff > 0._dp)then
      if(tempInc > critDiff) tempInc = critDiff + epsT  ! set iteration increment to slightly above critical temperature
     ! --> initially unfrozen (T > Tcrit)
     else
      if(tempInc < critDiff) tempInc = critDiff - epsT  ! set iteration increment to slightly below critical temperature
     endif
    endif  ! if the domain is soil
   
    ! update the temperature trial
    xTemp = xTemp + tempInc

    ! check failed convergence
    if(iter==maxiter)then
     message=trim(message)//'failed to converge'
     err=-20; return ! negative error code = try to recover
    endif
  
   endif   ! if adjusting the temperature

  end do iterations ! iterating

  ! save temperature 
  select case(ixDomainType)
   case(iname_veg);              scalarCanopyTempTrial   = xTemp
   case(iname_snow, iname_soil); mLayerTempTrial(iLayer) = xTemp
  end select
  
  ! =======================================================================================================================================
  ! =======================================================================================================================================

  ! -----
  ! - compute the liquid water matric potential (and necessay derivatives)...
  ! -------------------------------------------------------------------------

  ! only for soil
  if(ixDomainType==iname_soil)then

   ! check liquid water
   if(mLayerVolFracLiqTrial(iLayer) > theta_sat(ixControlIndex) )then
    message=trim(message)//'liquid water greater than porosity'
    err=20; return
   endif

   ! case of hydrology state uncoupled with energy
   if(.not.isNrgState .and. .not.isCoupled)then

    ! derivatives relating liquid water matric potential to total water matric potential and temperature
    dPsiLiq_dPsi0(ixControlIndex) = 1._dp  ! exact correspondence (psiLiq=psi0)
    dPsiLiq_dTemp(ixControlIndex) = 0._dp  ! no relationship between liquid water matric potential and temperature
  
   ! case of energy state or coupled solution 
   else
   
    ! compute the liquid matric potential (and the derivatives w.r.t. total matric potential and temperature)
    call liquidHead(&
                    ! input
                    mLayerMatricHeadTrial(ixControlIndex)                                                                                     ,& ! intent(in) : total water matric potential (m)
                    mLayerVolFracLiqTrial(iLayer)                                                                                             ,& ! intent(in) : volumetric fraction of liquid water (-)
                    mLayerVolFracIceTrial(iLayer)                                                                                             ,& ! intent(in) : volumetric fraction of ice (-)
                    vGn_alpha(ixControlIndex),vGn_n(ixControlIndex),theta_sat(ixControlIndex),theta_res(ixControlIndex),vGn_m(ixControlIndex), & ! intent(in) : soil parameters
                    dVolTot_dPsi0(ixControlIndex)                                                                                             ,& ! intent(in) : derivative in the soil water characteristic (m-1)
                    mLayerdTheta_dTk(iLayer)                                                                                                  ,& ! intent(in) : derivative in volumetric total water w.r.t. temperature (K-1)
                    ! output
                    mLayerMatricHeadLiqTrial(ixControlIndex)                                                                                  ,& ! intent(out): liquid water matric potential (m)
                    dPsiLiq_dPsi0(ixControlIndex)                                                                                             ,& ! intent(out): derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
                    dPsiLiq_dTemp(ixControlIndex)                                                                                             ,& ! intent(out): derivative in the liquid water matric potential w.r.t. temperature (m K-1)
                    err,cmessage)                                                                                                                ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   endif  ! switch between hydrology and energy state

  endif  ! if domain is soil

 end do ! looping through state variables

 ! deallocate space
 deallocate(computedCoupling,stat=err)        ! .true. if computed the coupling for a given state variable
 if(err/=0)then; message=trim(message)//'problem deallocating computedCoupling'; return; endif

 ! end association to the variables in the data structures
 end associate

 end subroutine updateVars


 ! **********************************************************************************************************
 ! private subroutine xTempSolve: compute residual and derivative for temperature
 ! **********************************************************************************************************
 subroutine xTempSolve(&
                       ! input: constant over iterations
                       meltNrg          ,&  ! intent(in)    : energy for melt+freeze (J m-3)
                       heatCap          ,&  ! intent(in)    : volumetric heat capacity (J m-3 K-1)
                       tempInit         ,&  ! intent(in)    : initial temperature (K)
                       volFracIceInit   ,&  ! intent(in)    : initial volumetric fraction of ice (-)
                       ! input-output: trial values
                       xTemp            ,&  ! intent(inout) : trial value of temperature
                       dLiq_dT          ,&  ! intent(in)    : derivative in liquid water content w.r.t. temperature (K-1)
                       volFracIceTrial  ,&  ! intent(in)    : trial value for volumetric fraction of ice
                       ! output: residual and derivative
                       residual         ,&  ! intent(out)   : residual (J m-3)
                       derivative        )  ! intent(out)   : derivative (J m-3 K-1)
 implicit none
 ! input: constant over iterations
 real(dp),intent(in)             :: meltNrg                         ! energy for melt+freeze (J m-3)
 real(dp),intent(in)             :: heatCap                         ! volumetric heat capacity (J m-3 K-1)
 real(dp),intent(in)             :: tempInit                        ! initial temperature (K)
 real(dp),intent(in)             :: volFracIceInit                  ! initial volumetric fraction of ice (-)
 ! input-output: trial values
 real(dp),intent(inout)          :: xTemp                           ! trial value for temperature
 real(dp),intent(in)             :: dLiq_dT                         ! derivative in liquid water content w.r.t. temperature (K-1)
 real(dp),intent(in)             :: volFracIceTrial                 ! trial value for the volumetric fraction of ice (-)
 ! output: residual and derivative
 real(dp),intent(out)            :: residual         ! residual (J m-3)
 real(dp),intent(out)            :: derivative       ! derivative (J m-3 K-1)
 ! subroutine starts here
 residual   = -heatCap*(xTemp - tempInit) + meltNrg*(volFracIceTrial - volFracIceInit)  ! J m-3
 derivative = heatCap + LH_fus*iden_water*dLiq_dT  ! J m-3 K-1
 end subroutine xTempSolve 

end module updateVars_module
