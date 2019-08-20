! cable_types.f90
!
! This file declares the main derived type variables for 
! CABLE and contains the interface that allocates them.
!
! Development by (alphabetically) Gab Abramowitz, Harvey Davies, 
! Martin Dix, Vanessa Haverd, Eva Kowalczyk, Ray Leuning,
! Mike Raupach, Ying-Ping Wang.
!
! bugs to bernard.pak@csiro.au
!
! The modules included are:
!   cable_types,
!     which contains subroutines:
!      alloc_*_type, and
!      dealloc_*_type

MODULE cable_types
  ! Contains all variables which are not subroutine-internal
  USE cable_dimensions, ONLY: r_1,r_2,i_d,ms,ncp,ncs,nrb,mf,ncstringlength
  IMPLICIT NONE
  PUBLIC
  PRIVATE r_1,r_2,i_d,ms,ncp,ncs,nrb,mf
  ! Energy and water balance variables:
  TYPE balances_type 
     REAL(r_2), DIMENSION(:), POINTER :: canopy_drybal ! energy balance for dry canopy
     REAL(r_2), DIMENSION(:), POINTER :: canopy_wetbal ! energy balance for wet canopy
     REAL(r_2), DIMENSION(:), POINTER :: evap_tot ! cumulative evapotranspiration (mm/dels)
     REAL(r_2), DIMENSION(:), POINTER :: osnowd0  ! snow depth, first time step
     REAL(r_2), DIMENSION(:), POINTER :: precip_tot ! cumulative precipitation (mm/dels)
     REAL(r_2), DIMENSION(:), POINTER :: rnoff_tot  ! cumulative runoff (mm/dels)
     REAL(r_2), DIMENSION(:), POINTER :: wbtot0 ! total soil water (mm), first time step
	 REAL(r_2), DIMENSION(:), POINTER :: Radbal ! radiation balance
	 REAL(r_2), DIMENSION(:), POINTER :: RadbalSum ! cumulative radiation balance
	 REAL(r_2), DIMENSION(:), POINTER :: Ebal ! energy balance
	 REAL(r_2), DIMENSION(:), POINTER :: EbalSum ! cumulative energy balance
	 REAL(r_2), DIMENSION(:), POINTER :: EbalSoil ! soil energy balance
	 REAL(r_2), DIMENSION(:), POINTER :: EbalSnow ! snow energy balance
	 REAL(r_2), DIMENSION(:), POINTER :: EbalVeg ! veg energy balance
	 REAL(r_2), DIMENSION(:), POINTER :: Wbal ! water balance (mm/dels)
	 REAL(r_2), DIMENSION(:), POINTER :: WbalSum ! water balance (mm/dels)
	 REAL(r_2), DIMENSION(:), POINTER :: WbalSoil ! soil water balance (mm/dels)
	 REAL(r_2), DIMENSION(:), POINTER :: WbalVeg ! water balance (wet canopy) (mm/dels)
	 REAL(r_2), DIMENSION(:), POINTER :: WbalSnow ! soil water balance (mm/dels)
  END TYPE balances_type
  ! Soil parameters:
  TYPE soil_parameter_type
     ! Default / SoilSnow parameters:
     REAL(r_1), DIMENSION(:), POINTER :: albsoil ! soil reflectance
     REAL(r_1), DIMENSION(:), POINTER :: bch  ! parameter b in Campbell equation
     REAL(r_1), DIMENSION(:), POINTER :: c3   ! c3 drainage coeff (fraction) (EK nov 2007)
     REAL(r_1), DIMENSION(:), POINTER :: clay ! fraction of soil which is clay
     REAL(r_1), DIMENSION(:), POINTER :: cnsd ! thermal conductivity of dry soil [W/m/K]
     REAL(r_1), DIMENSION(:), POINTER :: css  ! soil specific heat capacity [kJ/kg/K]
     REAL(r_1), DIMENSION(:), POINTER :: hsbh  ! difsat * etasat (=hyds*abs(sucs)*bch)
     REAL(r_1), DIMENSION(:), POINTER :: hyds  ! hydraulic conductivity @ saturation [m/s], Ksat
     INTEGER(i_d), DIMENSION(:), POINTER :: i2bp3  ! parameter one in K vis suction (=nint(bch)+2)
     INTEGER(i_d), DIMENSION(:), POINTER :: ibp2   ! parameter two in K vis suction (function of pbch)
     INTEGER(i_d), DIMENSION(:), POINTER :: isoilm ! integer soil type
     REAL(r_1), DIMENSION(:), POINTER :: rhosoil ! soil density [kg/m3]
     REAL(r_1), DIMENSION(:), POINTER :: rs20  ! soil respiration at 20 C [mol m-2 s-1]
     REAL(r_1), DIMENSION(:), POINTER :: sand  ! fraction of soil which is sand
     REAL(r_1), DIMENSION(:), POINTER :: sfc   ! vol H2O @ field capacity
     REAL(r_1), DIMENSION(:), POINTER :: silt  ! fraction of soil which is silt
     REAL(r_1), DIMENSION(:), POINTER :: ssat  ! vol H2O @ saturation
     REAL(r_1), DIMENSION(:), POINTER :: sucs  ! suction at saturation (m)
     REAL(r_1), DIMENSION(:), POINTER :: swilt ! vol H2O @ wilting
     REAL(r_1), DIMENSION(ms) :: zse   ! thickness of each soil layer (1=top) in m
     REAL(r_1), DIMENSION(ms+1) :: zshh ! distance between consecutive layer midpoints (m)
     ! Additional SLI parameters: 'B' suffix implies 2nd soil horizon parameters
     INTEGER(i_d), DIMENSION(:), POINTER :: nhorizons ! number of soil horizons
     REAL(r_1), DIMENSION(:), POINTER :: bchB  ! parameter b in Campbell equation  
     REAL(r_1), DIMENSION(:), POINTER :: clayB ! fraction of soil which is clay    
     REAL(r_1), DIMENSION(:), POINTER :: cnsdB ! thermal conductivity of dry soil [W/m/K]  
     REAL(r_1), DIMENSION(:), POINTER :: cssB  ! soil specific heat capacity [kJ/kg/K]      
     REAL(r_1), DIMENSION(:), POINTER :: hydsB  ! hydraulic conductivity @ saturation [m/s], Ksat    
     REAL(r_1), DIMENSION(:), POINTER :: rhosoilB ! soil density [kg/m3]   
     REAL(r_1), DIMENSION(:), POINTER :: sandB  ! fraction of soil which is sand  
     REAL(r_1), DIMENSION(:), POINTER :: siltB  ! fraction of soil which is silt   
     REAL(r_1), DIMENSION(:), POINTER :: ssatB  ! vol H2O @ saturation          
     REAL(r_1), DIMENSION(:), POINTER :: sucsB  ! suction at saturation (m)     
     REAL(r_1), DIMENSION(:), POINTER :: swiltB ! vol H2O @ wilting             
     REAL(r_1), DIMENSION(:), POINTER :: depthA ! depth of A horizon             
     REAL(r_1), DIMENSION(:), POINTER :: depthB ! depth of B horizon           
     REAL(r_1), DIMENSION(:), POINTER :: sfcB    ! field capacity of B horizon
     REAL(r_2), DIMENSION(:), POINTER  :: clitt      ! litter (tC/ha)
     REAL(r_1), DIMENSION(:,:), POINTER :: swilt_vec ! vol H2O @ wilting
  END TYPE soil_parameter_type
  ! Soil and snow variables:
  TYPE soil_snow_type 
     ! Default / SoilSnow variables:
     REAL(r_2), DIMENSION(:,:), POINTER :: albsoilsn ! soil + snow reflectance
     REAL(r_2), DIMENSION(:), POINTER :: cls     ! factor for latent heat
     REAL(r_2), DIMENSION(:), POINTER :: dfn_dtg ! d(canopy%fns)/d(ssoil%tgg)
     REAL(r_2), DIMENSION(:), POINTER :: dfh_dtg ! d(canopy%fhs)/d(ssoil%tgg)
     REAL(r_2), DIMENSION(:), POINTER :: dfe_ddq ! d(canopy%fes)/d(dq)
     REAL(r_2), DIMENSION(:), POINTER :: ddq_dtg ! d(dq)/d(ssoil%tgg)
     REAL(r_2), DIMENSION(:), POINTER :: evapsn  ! snow evaporation
     REAL(r_2), DIMENSION(:), POINTER :: fwtop   ! water flux to the soil (EK nov 2007)
     REAL(r_2), DIMENSION(:,:), POINTER :: gammzz ! heat capacity for each soil layer
     INTEGER(i_d), DIMENSION(:), POINTER :: isflag ! 0 => no snow 1 => snow
     REAL(r_2), DIMENSION(:), POINTER :: osnowd  ! snow depth from previous time step
     REAL(r_2), DIMENSION(:), POINTER :: potev   ! potential evapotranspiration
     REAL(r_2), DIMENSION(:), POINTER :: pwb_min
     REAL(r_2), DIMENSION(:), POINTER :: runoff  ! total runoff (mm/dels)
     REAL(r_2), DIMENSION(:), POINTER :: rnof1   ! surface runoff (mm/dels)
     REAL(r_2), DIMENSION(:), POINTER :: rnof2   ! deep drainage (mm/dels)
     REAL(r_2), DIMENSION(:), POINTER :: rtsoil  ! turbulent resistance for soil
     REAL(r_2), DIMENSION(:,:), POINTER :: sconds ! EK nov 2007
     REAL(r_2), DIMENSION(:,:), POINTER :: sdepth ! snow depth
     REAL(r_2), DIMENSION(:,:), POINTER  :: smass ! snow mass
     REAL(r_2), DIMENSION(:), POINTER :: snage   ! snow age
     REAL(r_2), DIMENSION(:), POINTER :: snowd   ! snow depth (liquid water)
     REAL(r_2), DIMENSION(:), POINTER :: smelt   ! snow melt (EK nov 2007)
     REAL(r_2), DIMENSION(:,:), POINTER  :: ssdn ! snow densities
     REAL(r_2), DIMENSION(:), POINTER :: ssdnn   ! average snow density
     REAL(r_2), DIMENSION(:,:), POINTER :: tgg   ! soil temperature in K
     REAL(r_2), DIMENSION(:,:), POINTER  :: tggsn ! snow temperature in K
     REAL(r_1), DIMENSION(:), POINTER :: tss     ! surface temperature (weighted soil, snow)
     REAL(r_1), DIMENSION(:), POINTER :: tssold  ! previous time step surface temperature (weighted soil, snow)
     REAL(r_2), DIMENSION(:,:), POINTER :: wb    ! volumetric soil moisture (solid+liq)
     REAL(r_2), DIMENSION(:,:), POINTER :: wbfice ! fraction of ssat that is ice
     REAL(r_2), DIMENSION(:,:), POINTER :: wbice  ! volumentric soil ice
     REAL(r_2), DIMENSION(:,:), POINTER :: wblf  ! fraction of ssat that is liquid
     REAL(r_2), DIMENSION(:), POINTER :: wbtot   ! total soil water (mm)
     REAL(r_1), DIMENSION(:), POINTER :: wetfac ! surface wetness fact. at current time step
     REAL(r_1), DIMENSION(:), POINTER :: owetfac ! surface wetness fact. at previous time step
     ! Additional SLI variables:
     REAL(r_2), DIMENSION(:,:), POINTER  :: S      ! moisture content relative to sat value    (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:), POINTER :: SL        ! litter moisture content relative to sat value (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:), POINTER :: TL        ! litter temperature in K     (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:), POINTER :: h0        ! pond height in m            (edit vh 23/01/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: rex    ! root extraction from each layer (mm/dels)
     REAL(r_2), DIMENSION(:,:), POINTER :: wflux    ! water flux at layer boundaries (mm s-1)
     REAL(r_2), DIMENSION(:,:), POINTER :: hflux    ! heat flux at layer boundaries (J m-2 s-1)
     REAL(r_2), DIMENSION(:), POINTER :: delwcol    ! change in water column (mm / dels)
     REAL(r_2), DIMENSION(:), POINTER :: zdelta        ! water table depth           (edit vh 23/06/08)
     REAL(r_2), DIMENSION(:,:), POINTER :: kth        ! thermal conductivity           (edit vh 29/07/08)
     REAL(r_2), DIMENSION(:), POINTER :: Tsurface       !  tepmerature at surface (soil, pond or litter) (edit vh 22/10/08)
     REAL(r_2), DIMENSION(:), POINTER :: rh0       !  relative humidity at top of soil column (edit vh 22/10/08)
     REAL(r_2), DIMENSION(:), POINTER :: rhsurface       ! relative at surface (soil, pond or litter) (edit vh 22/10/08)
     REAL(r_2), DIMENSION(:), POINTER :: lE              ! soil latent heat flux 
     REAL(r_2), DIMENSION(:), POINTER  :: evap			! soil evaporation (mm / dels)
     REAL(r_2), DIMENSION(:,:), POINTER  :: ciso      ! concentration of minor isotopologue in soil water (kg m-3 water)   (edit vh 8/01/09)
     REAL(r_2), DIMENSION(:), POINTER  :: cisoL      ! concentration of minor isotopologue in litter water (kg m-3 water)   (edit vh 8/01/09)
     REAL(r_2), DIMENSION(:), POINTER :: delwcolA    ! change in water column (A horizon) (mm / dels)
     REAL(r_2), DIMENSION(:), POINTER :: delwcolB    ! change in water column (B horizon) (mm / dels)
     REAL(r_2), DIMENSION(:), POINTER :: rexA    ! root extraction (A horizon) (mm / dels)
     REAL(r_2), DIMENSION(:), POINTER :: rexB    ! root extraction (B horizon) (mm / dels)
     REAL(r_2), DIMENSION(:), POINTER  :: SA      ! degree of saturation  (A horizon)
     REAL(r_2), DIMENSION(:), POINTER  :: SB      ! degree of saturation  (B horizon)
     REAL(r_2), DIMENSION(:), POINTER  :: leachAB      ! downward flux between A and B horizons (mm/dels)
     REAL(r_2), DIMENSION(:), POINTER  :: rlitt      ! resistance to heat/moisture transfer through litter (m-1 s)
     REAL(r_2), DIMENSION(:,:), POINTER  :: thetai  ! volumetric ice content (MC)
  END TYPE soil_snow_type
  ! Vegetation parameters:
  TYPE veg_parameter_type
     ! Original CABLE / CBM parameters:
     REAL(r_1), DIMENSION(:), POINTER :: canst1 ! max intercepted water by canopy (mm/LAI)
     LOGICAL,   DIMENSION(:), POINTER :: deciduous ! flag used for phenology fix
     REAL(r_1), DIMENSION(:), POINTER :: dleaf  ! chararacteristc legnth of leaf (m)
     REAL(r_1), DIMENSION(:), POINTER :: ejmax  ! max pot. electron transport rate top leaf(mol/m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: extkn  ! extinction coef for vertical
     REAL(r_1), DIMENSION(:), POINTER :: frac4  ! fraction of c4 plants
     REAL(r_1), DIMENSION(:,:), POINTER :: froot  ! fraction of root in each soil layer
     REAL(r_1), DIMENSION(:), POINTER :: hc     ! roughness height of canopy
     INTEGER(i_d),DIMENSION(:), POINTER :: iveg ! vegetation type
     INTEGER(i_d),DIMENSION(:), POINTER :: meth ! method for calculation of canopy fluxes and temp.
     REAL(r_1), DIMENSION(:), POINTER :: rp20   ! plant respiration coefficient at 20 C
     REAL(r_1), DIMENSION(:), POINTER :: rpcoef ! temperature coef nonleaf plant respiration (1/C)
     REAL(r_1), DIMENSION(:), POINTER :: shelrb ! sheltering factor (dimensionless) 
     REAL(r_1), DIMENSION(:), POINTER :: tmaxvj ! max temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(:), POINTER :: tminvj ! min temperature of the start of photosynthesis
     REAL(r_1), DIMENSION(:), POINTER :: vbeta  ! stomatal sensitivity to soil water
     REAL(r_1), DIMENSION(:), POINTER :: vegcf  ! biome-specific soil respiration rate
     REAL(r_1), DIMENSION(:), POINTER :: vcmax  ! maximum RuBP carboxylation rate top leaf (mol/m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: vlai   ! leaf area index
     REAL(r_1), DIMENSION(:), POINTER :: wai   ! wood area index (stem+branches+twigs)
     REAL(r_1), DIMENSION(:), POINTER :: xfang  ! leaf angle PARAMETER
     ! Additional SLI veg parameters:
     REAL(r_1), DIMENSION(:), POINTER :: rootbeta  ! parameter for estimating vertical root mass distribution (froot)
     REAL(r_1), DIMENSION(:), POINTER :: gamma ! parameter in root efficiency function (Lai and Katul 2000)
     REAL(r_1), DIMENSION(:), POINTER :: ZR ! maximum rooting depth (cm)
     REAL(r_1), DIMENSION(:), POINTER :: F10 ! fraction of roots in top 10 cm
  END TYPE veg_parameter_type
  ! Canopy/vegetation variables:
  TYPE canopy_type
     REAL(r_2), DIMENSION(:), POINTER :: cansto ! canopy water storage (mm)
     REAL(r_1), DIMENSION(:), POINTER :: cduv  ! drag coefficient for momentum
     REAL(r_2), DIMENSION(:), POINTER :: delwc ! change in canopy water store (mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: dewmm ! dewfall (mm)
     REAL(r_2), DIMENSION(:), POINTER :: dgdtg ! derivative of gflux wrt soil temp
     REAL(r_1), DIMENSION(:), POINTER :: fe    ! total latent heat (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fh    ! total sensible heat (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fpn   ! plant photosynthesis (g C s-1)
     REAL(r_1), DIMENSION(:), POINTER :: frp   ! plant respiration (g C m-2 s-1)
     REAL(r_1), DIMENSION(:), POINTER :: frpw  ! plant respiration (g C m-2 s-1)???
     REAL(r_1), DIMENSION(:), POINTER :: frpr  ! plant respiration (g C m-2 s-1)???
     REAL(r_1), DIMENSION(:), POINTER :: frs   ! soil respiration (g C m-2 s-1)
     REAL(r_1), DIMENSION(:), POINTER :: fnee  ! net carbon flux (g C m-2 s-1)
     REAL(r_1), DIMENSION(:), POINTER :: frday ! daytime leaf resp
     REAL(r_1), DIMENSION(:), POINTER :: fnv   ! net rad. avail. to canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fev   ! latent hf from canopy (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fevc  ! dry canopy transpiration (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fevw  ! lat heat fl wet canopy (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: potev_c ! canopy potential evapotranspitation (YP & Mao, jun08)
     REAL(r_1), DIMENSION(:), POINTER :: fhv   ! sens heatfl from canopy (W/m2)
     REAL(r_2), DIMENSION(:), POINTER :: fhvw  ! sens heatfl from wet canopy (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fns   ! net rad avail to soil (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fes   ! latent heatfl from soil (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fhs   ! sensible heat flux from soil
     REAL(r_1), DIMENSION(:), POINTER :: fwet   ! fraction of canopy wet
     REAL(r_1), DIMENSION(:), POINTER :: ga    ! ground heat flux (W/m2) ???
     REAL(r_1), DIMENSION(:), POINTER :: ghflux  ! ground heat flux (W/m2) ???
     REAL(r_2), DIMENSION(:), POINTER :: precis ! throughfall to soil, after snow (mm)
     REAL(r_1), DIMENSION(:), POINTER :: qscrn ! specific humudity at screen height (g/g)
     REAL(r_1), DIMENSION(:), POINTER :: rnet  ! net radiation absorbed by surface (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: segg  ! latent heatfl from soil mm (EK nov 2007)
     REAL(r_1), DIMENSION(:), POINTER :: sghflux ! ground heat flux (W/m2) ???
     REAL(r_1), DIMENSION(:), POINTER :: spill ! can.storage excess after dewfall (mm)
     REAL(r_1), DIMENSION(:), POINTER :: through ! canopy throughfall (mm)
     REAL(r_1), DIMENSION(:), POINTER :: tscrn ! air temperature at screen height (oC)
     REAL(r_1), DIMENSION(:), POINTER :: tv    ! vegetation temp (K)
     REAL(r_1), DIMENSION(:), POINTER :: us    ! friction velocity
     REAL(r_1), DIMENSION(:), POINTER :: uscrn ! wind speed at screen height (m/s)
     REAL(r_1), DIMENSION(:), POINTER :: vlaiw  ! lai adjusted for snow depth for calculation of resistances
     REAL(r_2), DIMENSION(:), POINTER :: wcint ! canopy rainfall interception (mm)
     ! Additional variables:
     REAL(r_1), DIMENSION(:,:), POINTER :: gw    ! dry canopy conductance (ms-1) edit vh 6/7/09
     REAL(r_1), DIMENSION(:,:), POINTER :: gswx   ! stomatal conductance (ms-1) edit vh 6/7/09
     REAL(r_1), DIMENSION(:,:,:), POINTER :: ancj    ! limiting photosynthetic rates (Rubisco,RuBP,sink) vh 6/7/09
     REAL(r_1), DIMENSION(:,:), POINTER :: tlfy    ! sunlit and shaded leaf temperatures
     REAL(r_1), DIMENSION(:,:), POINTER :: ecy    ! sunlit and shaded leaf transpiration (dry canopy)
     REAL(r_1), DIMENSION(:,:), POINTER :: ecx    ! sunlit and shaded leaf latent heat flux
     REAL(r_1), DIMENSION(:,:,:), POINTER :: ci    ! intra-cellular CO2 vh 6/7/09
     REAL(r_1), DIMENSION(:), POINTER :: fwsoil    !
  END TYPE canopy_type
  ! Radiation variables:
  TYPE radiation_type
     REAL(r_1), DIMENSION(:,:), POINTER :: albedo ! canopy+soil albedo
     REAL(r_1), DIMENSION(:), POINTER     :: extkb  ! beam radiation extinction coeff
     REAL(r_1), DIMENSION(:), POINTER     :: extkd2 ! diffuse 2D radiation extinction coeff
     REAL(r_1), DIMENSION(:), POINTER ::  extkd ! diffuse radiation extinction coeff (-)
     REAL(r_1), DIMENSION(:), POINTER :: flws   ! soil long-wave radiation
     REAL(r_1), DIMENSION(:,:), POINTER  :: fvlai  ! leaf area index of big leaf
     REAL(r_2), DIMENSION(:,:), POINTER  :: gradis ! radiative conductance
     REAL(r_1), DIMENSION(:), POINTER :: latitude  ! latitude
     REAL(r_1), DIMENSION(:), POINTER :: lwabv ! non-isothermal long wave emission by vegetation
     REAL(r_1), DIMENSION(:,:,:), POINTER :: qcan ! absorbed radiation for canopy (W/m^2)
     REAL(r_1), DIMENSION(:), POINTER     :: qssabs ! absorbed short-wave radiation for soil
     REAL(r_1), DIMENSION(:,:), POINTER ::  rhocdf ! canopy diffuse reflectance (-)
     REAL(r_1), DIMENSION(:,:), POINTER  :: rniso  !  sum(rad%qcan, 3) total abs by canopy (W/m2)
     REAL(r_1), DIMENSION(:,:), POINTER  :: scalex ! scaling PARAMETER for big leaf
     REAL(r_1), DIMENSION(:), POINTER     :: transd ! fraction SW diffuse transmitted through canopy
     REAL(r_1), DIMENSION(:), POINTER ::  trad  !  radiative temperature (soil and veg)
     ! new types, ypw 11/july/2008
     REAL(r_1),DIMENSION(:,:), POINTER  :: reffdf  !effective conopy diffuse reflectance
     REAL(r_1),DIMENSION(:,:), POINTER  :: reffbm  !effective conopy beam reflectance
     REAL(r_1), DIMENSION(:,:), POINTER :: extkbm  !modified k beam(6.20)(for leaf scattering)
     REAL(r_1), DIMENSION(:,:), POINTER :: extkdm  !modified k diffuse(6.20)(for leaf scattering)
     REAL(r_1), DIMENSION(:), POINTER   :: fbeam   !beam fraction
     REAL(r_1), DIMENSION(:,:),POINTER  :: cexpkbm ! canopy beam transmittance
     REAL(r_1), DIMENSION(:,:),POINTER  :: cexpkdm ! canopy diffuse transmittance
     ! ********* ypw 11/july/2008
  END TYPE radiation_type
  ! Roughness variables:
  TYPE roughness_type
     ! "coexp": coefficient in exponential in-canopy wind profile
     ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
     ! canopy and roughness-sublayer U(z) at z=h
     REAL(r_1), DIMENSION(:), POINTER   :: coexp ! Extinction coefficient for wind profile in canopy
     REAL(r_1), DIMENSION(:), POINTER   :: disp  ! zero-plane displacement
     REAL(r_1), DIMENSION(:), POINTER   :: hruff ! canopy height above snow level
     REAL(r_1), DIMENSION(:), POINTER   :: rt0us ! eq. 3.54, SCAM manual (CSIRO tech report 132)
     REAL(r_1), DIMENSION(:), POINTER   :: rt1usa ! resistance from disp to hruf
     REAL(r_1), DIMENSION(:), POINTER   :: rt1usb ! resistance from hruf to zruffs (or zref if zref<zruffs)
     REAL(r_1), DIMENSION(:), POINTER   :: rt1 ! 1/aerodynamic conductance
     REAL(r_1), DIMENSION(:), POINTER   :: term2, term3, term5, term6 ! for aerodynamic resistance calc.
     ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
     REAL(r_1), DIMENSION(:), POINTER   :: usuh ! Friction velocity/windspeed at canopy height
     REAL(r_1), DIMENSION(:), POINTER   :: za   ! level of lowest atmospheric model layer
     REAL(r_1), DIMENSION(:), POINTER  :: za_uv ! level of lowest uv grid
     REAL(r_1), DIMENSION(:), POINTER  :: za_tq ! level of lowest tq grid
     REAL(r_1), DIMENSION(:), POINTER   :: z0m  ! roughness length
     REAL(r_1), DIMENSION(:), POINTER   :: zref ! Reference height for met forcing
     REAL(r_1), DIMENSION(:), POINTER  :: zref_uv ! Ref height wrt uv grid
     REAL(r_1), DIMENSION(:), POINTER  :: zref_tq ! Ref height wrt tq grid
     REAL(r_1), DIMENSION(:), POINTER   :: zruffs ! SCALAR Roughness sublayer depth (ground=origin)
     REAL(r_1), DIMENSION(:), POINTER   :: z0soilsn ! roughness length of bare soil surface
     REAL(r_1), DIMENSION(:), POINTER   :: z0soil ! roughness length of bare soil surface
  END TYPE roughness_type
  ! Air variables:
  TYPE air_type
     REAL(r_1), DIMENSION(:), POINTER :: rho  ! dry air density (kg m-3)
     REAL(r_1), DIMENSION(:), POINTER :: volm ! molar volume (m3 mol-1)
     REAL(r_1), DIMENSION(:), POINTER :: rlam ! latent heat for water (j/kg)
     REAL(r_1), DIMENSION(:), POINTER :: qsat ! saturation specific humidity
     REAL(r_1), DIMENSION(:), POINTER :: epsi ! d(qsat)/dT ((kg/kg)/K)
     REAL(r_1), DIMENSION(:), POINTER :: visc ! air kinematic viscosity (m2/s)
     REAL(r_1), DIMENSION(:), POINTER :: psyc ! psychrometric constant
     REAL(r_1), DIMENSION(:), POINTER :: dsatdk ! d(es)/dT (mb/K)
     REAL(r_1), DIMENSION(:), POINTER :: cmolar ! conv. from m/s to mol/m2/s
  END TYPE air_type
  ! Flux variables from a homogenous patch (subset of a grid cell)
  TYPE patch_flux_type
     REAL(r_1), DIMENSION(:), POINTER :: fe ! latent heat (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fh ! sensible heat (W/m2)
  END TYPE patch_flux_type
  ! Meterological data:
  TYPE met_type
     REAL(r_1), DIMENSION(:), POINTER :: ca   ! CO2 concentration (mol/mol)
     INTEGER(i_d), DIMENSION(:), POINTER :: year ! local time year AD 
     INTEGER(i_d), DIMENSION(:), POINTER :: moy  ! local time month of year 
     REAL(r_1), DIMENSION(:), POINTER :: doy  ! local time day of year = days since 0 hr 1st Jan 
     REAL(r_1), DIMENSION(:), POINTER :: hod  ! local hour of day
     REAL(r_1), DIMENSION(:), POINTER :: fsd  ! downward short-wave radiation (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: fld  ! downward long-wave radiation (W/m2)
     REAL(r_1), DIMENSION(:), POINTER :: precip  ! rainfall (liquid+solid)(mm/dels)
     REAL(r_1), DIMENSION(:), POINTER :: precip_s ! solid preipitation only (mm/dels) (EK nov 2007)
     REAL(r_1), DIMENSION(:), POINTER :: tc    ! surface air temperature (oC)
     REAL(r_1), DIMENSION(:), POINTER :: tk    ! surface air temperature (oK)
     REAL(r_1), DIMENSION(:), POINTER :: tvair ! within canopy air temperature (oK)
     REAL(r_1), DIMENSION(:), POINTER :: tvrad ! radiative vegetation temperature (K)
     REAL(r_1), DIMENSION(:), POINTER :: pmb   ! surface air pressure (mbar)
     REAL(r_1), DIMENSION(:), POINTER :: ua    ! surface wind speed (m/s)
     REAL(r_1), DIMENSION(:), POINTER :: qv    ! surface specific humidity (g/g)
     REAL(r_1), DIMENSION(:), POINTER :: qvair ! within canopy specific humidity (g/g)
     REAL(r_1), DIMENSION(:), POINTER :: da    ! water vap pressure deficit at ref height (Pa)
     REAL(r_1), DIMENSION(:), POINTER :: dva   ! in canopy water vap pressure deficit (Pa)
     REAL(r_1), DIMENSION(:), POINTER :: coszen  ! cos(zenith angle of sun)
     REAL(r_1), DIMENSION(:), POINTER :: tk_old   ! surface air temperature of time step before (oK)
  END TYPE met_type
  ! Cumulative flux variables:
  TYPE sum_flux_type
     REAL(r_1), DIMENSION(:), POINTER :: sumpn  ! sum of canopy photosynthesis (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER :: sumrp  ! sum of plant respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER :: sumrpw ! sum of plant respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER :: sumrpr ! sum of plant respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER :: sumrs  ! sum of soil respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER :: sumrd  ! sum of daytime respiration (g C m-2)
     REAL(r_1), DIMENSION(:), POINTER :: dsumpn ! daily sumpn
     REAL(r_1), DIMENSION(:), POINTER :: dsumrp ! daily sumrp
     REAL(r_1), DIMENSION(:), POINTER :: dsumrs ! daily sumrs
     REAL(r_1), DIMENSION(:), POINTER :: dsumrd ! daily sumrd
     REAL(r_1), DIMENSION(:), POINTER :: sumxrp ! sum plant resp. modifier
     REAL(r_1), DIMENSION(:), POINTER :: sumxrs ! sum soil resp. modifier
  END TYPE sum_flux_type
  TYPE bgc_pool_type
     ! Default BGC vars:
     REAL(r_1), DIMENSION(:,:), POINTER :: cplant ! plant carbon (g C/m2))
     REAL(r_1), DIMENSION(:,:), POINTER :: csoil  ! soil carbon (g C/m2)
     REAL(r_1), DIMENSION(ncp) :: ratecp ! plant carbon rate constant (1/year)
     REAL(r_1), DIMENSION(ncs) :: ratecs ! soil carbon rate constant (1/year)
     ! Additional SLI bgc variables
     REAL(r_1), DIMENSION(:), POINTER :: clitt  ! litter carbon (g C/m2)
  END TYPE bgc_pool_type
  TYPE model_structure_type ! All model structure variables below accept 'default'
     CHARACTER(LEN=ncstringlength),DIMENSION(:),POINTER :: canopy  ! optional canopy structures - see cable.nml
     CHARACTER(LEN=ncstringlength),DIMENSION(:),POINTER :: soil    ! optional canopy structures - see cable.nml
     CHARACTER(LEN=ncstringlength),DIMENSION(:),POINTER :: sli_isotope ! optional canopy structures - see cable.nml
     CHARACTER(LEN=ncstringlength),DIMENSION(:),POINTER :: photosynthesis ! optional canopy structures - see cable.nml
     CHARACTER(LEN=ncstringlength),DIMENSION(:),POINTER :: sli_litter ! optional soil litter structures - see cable.nml
     CHARACTER(LEN=ncstringlength),DIMENSION(:),POINTER :: sli_coupled ! optional soil energy/water structures - see cable.nml
  END TYPE model_structure_type

  ! Functions for allocating these types
  ! All overloaded so code only needs to call alloc_cbm_var
  ! Alloc routines could all initialise to NaN or zero for debugging?
  PUBLIC :: alloc_cbm_var
  PRIVATE :: alloc_bgc_pool_type, dealloc_bgc_pool_type
  INTERFACE alloc_cbm_var
     MODULE PROCEDURE alloc_balances_type,            &
          alloc_soil_parameter_type,      &
          alloc_soil_snow_type,           &
          alloc_veg_parameter_type,       &
          alloc_canopy_type,              &
          alloc_radiation_type,           &
          alloc_roughness_type,           &
          alloc_air_type,                 &
          alloc_met_type,                 &
          alloc_sum_flux_type,            &
          alloc_bgc_pool_type,            &
          alloc_model_structure_type 
  END INTERFACE
  INTERFACE dealloc_cbm_var
     MODULE PROCEDURE dealloc_balances_type,            &
          dealloc_soil_parameter_type,      &
          dealloc_soil_snow_type,           &
          dealloc_veg_parameter_type,       &
          dealloc_canopy_type,              &
          dealloc_radiation_type,           &
          dealloc_roughness_type,           &
          dealloc_air_type,                 &
          dealloc_met_type,                 &
          dealloc_sum_flux_type,            &
          dealloc_bgc_pool_type,            &
          dealloc_model_structure_type
  END INTERFACE
CONTAINS

  SUBROUTINE alloc_balances_type(var,landunits,model_structure)
    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % canopy_drybal(landunits) )
    ALLOCATE ( var % evap_tot(landunits) )
    ALLOCATE ( var % osnowd0(landunits) )
    ALLOCATE ( var % precip_tot(landunits) )
    ALLOCATE ( var % rnoff_tot(landunits) )
    ALLOCATE ( var % wbtot0(landunits) )
    ALLOCATE ( var % canopy_wetbal(landunits) )
	ALLOCATE ( var % Radbal(landunits) )
	ALLOCATE ( var % RadbalSum(landunits) )
	ALLOCATE ( var % Ebal(landunits) )
	ALLOCATE ( var % EbalSum(landunits) )
	ALLOCATE ( var % EbalSoil(landunits) )
	ALLOCATE ( var % EbalSnow(landunits) )
	ALLOCATE ( var % EbalVeg(landunits) )
	ALLOCATE ( var % Wbal(landunits) )
	ALLOCATE ( var % WbalSum(landunits) )
	ALLOCATE ( var % WbalSoil(landunits) )
	ALLOCATE ( var % WbalVeg(landunits) )
	ALLOCATE ( var % WbalSnow(landunits) )
  END SUBROUTINE alloc_balances_type

  SUBROUTINE alloc_soil_parameter_type(var,landunits,model_structure)
    TYPE(soil_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % albsoil(landunits) )
    ALLOCATE ( var % bch(landunits) )
    ALLOCATE ( var % c3(landunits) )
    ALLOCATE ( var % clay(landunits) )
    ALLOCATE ( var % cnsd(landunits) )
    ALLOCATE ( var % css(landunits) )
    ALLOCATE ( var % hsbh(landunits) )
    ALLOCATE ( var % hyds(landunits) )
    ALLOCATE ( var % i2bp3(landunits) )
    ALLOCATE ( var % ibp2(landunits) )
    ALLOCATE ( var % isoilm(landunits) )
    ALLOCATE ( var % rhosoil(landunits) )
    ALLOCATE ( var % rs20(landunits) )
    ALLOCATE ( var % sand(landunits) )
    ALLOCATE ( var % sfc(landunits) )
    ALLOCATE ( var % silt(landunits) )
    ALLOCATE ( var % ssat(landunits) )
    ALLOCATE ( var % sucs(landunits) )
    ALLOCATE ( var % swilt(landunits) )
    ! Allocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       ALLOCATE ( var % nhorizons(landunits) )
       ALLOCATE ( var % bchB(landunits) )
       ALLOCATE ( var % clayB(landunits) )
       ALLOCATE ( var % cnsdB(landunits) )
       ALLOCATE ( var % cssB(landunits) )
       ALLOCATE ( var % hydsB(landunits) )
       ALLOCATE ( var % rhosoilB(landunits) )
       ALLOCATE ( var % sandB(landunits) )
       ALLOCATE ( var % siltB(landunits) )
       ALLOCATE ( var % ssatB(landunits) )
       ALLOCATE ( var % sucsB(landunits) )
       ALLOCATE ( var % swiltB(landunits) )
       ALLOCATE ( var % depthA(landunits) )
       ALLOCATE ( var % depthB(landunits) )
       ALLOCATE ( var % sfcB(landunits) )
       ALLOCATE ( var % clitt(landunits) )
       ALLOCATE ( var % swilt_vec(landunits,ms) )
    END IF
    ! Allocate variables for canopy_vh soil model:
    IF(ANY(model_structure%canopy=='canopy_vh')) THEN
       IF(.NOT.(ASSOCIATED(var % swilt_vec))) ALLOCATE ( var % swilt_vec(landunits,ms) )
    END IF

  END SUBROUTINE alloc_soil_parameter_type

  SUBROUTINE alloc_soil_snow_type(var,landunits,model_structure)
    TYPE(soil_snow_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % albsoilsn(landunits,nrb) )
    ALLOCATE ( var % cls(landunits) )
    ALLOCATE ( var % dfn_dtg(landunits) )
    ALLOCATE ( var % dfh_dtg(landunits) )
    ALLOCATE ( var % dfe_ddq(landunits) )
    ALLOCATE ( var % ddq_dtg(landunits) )
    ALLOCATE ( var % evapsn(landunits) )
    ALLOCATE ( var % fwtop(landunits) )
    ALLOCATE ( var % gammzz(landunits,ms) )
    ALLOCATE ( var % isflag(landunits) )
    ALLOCATE ( var % osnowd(landunits) )
    ALLOCATE ( var % potev(landunits) )
    ALLOCATE ( var % pwb_min(landunits) )
    ALLOCATE ( var % runoff(landunits) )
    ALLOCATE ( var % rnof1(landunits) )
    ALLOCATE ( var % rnof2(landunits) )
    ALLOCATE ( var % rtsoil(landunits) )
    ALLOCATE ( var % sconds(landunits,3) )
    ALLOCATE ( var % sdepth(landunits,3) )
    ALLOCATE ( var % smass(landunits,3) )
    ALLOCATE ( var % snage(landunits) )
    ALLOCATE ( var % snowd(landunits) )
    ALLOCATE ( var % smelt(landunits) )
    ALLOCATE ( var % ssdn(landunits,3) )
    ALLOCATE ( var % ssdnn(landunits) )
    ALLOCATE ( var % tgg(landunits,ms) )
    ALLOCATE ( var % tggsn(landunits,3) )
    ALLOCATE ( var % tss(landunits) )
    ALLOCATE ( var % tssold(landunits) )
    ALLOCATE ( var % wb(landunits,ms) )
    ALLOCATE ( var % wbfice(landunits,ms) )
    ALLOCATE ( var % wbice(landunits,ms) )
    ALLOCATE ( var % wblf(landunits,ms) )
    ALLOCATE ( var % wbtot(landunits) )
    ALLOCATE ( var % wetfac(landunits) )
    ALLOCATE ( var % owetfac(landunits) )
    ALLOCATE ( var % evap(landunits) )
    ALLOCATE ( var % delwcol(landunits) )
    ! Allocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       ALLOCATE ( var % S(landunits,ms) )	
       ALLOCATE ( var % SL(landunits) )	
       ALLOCATE ( var % TL(landunits) )	
       ALLOCATE ( var % h0(landunits) )	
       ALLOCATE ( var % rex(landunits,ms) )	
       ALLOCATE ( var % wflux(landunits,0:ms) )	
       ALLOCATE ( var % hflux(landunits,0:ms) )	       	
       ALLOCATE ( var % zdelta(landunits) )
       ALLOCATE ( var % kth(landunits,ms) )	
       ALLOCATE ( var % Tsurface(landunits) )	
       ALLOCATE ( var % rh0(landunits) )	
       ALLOCATE ( var % rhsurface(landunits) )	
       ALLOCATE ( var % lE(landunits) )	 	    
       ALLOCATE ( var % ciso(landunits,ms+1) )	
       ALLOCATE ( var % cisoL(landunits) )	   
       ALLOCATE ( var % SA(landunits) )
       ALLOCATE ( var % SB(landunits) )	
       ALLOCATE ( var % delwcolA(landunits) )	
       ALLOCATE ( var % delwcolB(landunits) )	
       ALLOCATE ( var % rexA(landunits) )	
       ALLOCATE ( var % rexB(landunits) )	
       ALLOCATE ( var % leachAB(landunits) ) 
       ALLOCATE ( var % rlitt(landunits) ) 
       ALLOCATE ( var % thetai(landunits,ms) ) 
    END IF
  END SUBROUTINE alloc_soil_snow_type

  SUBROUTINE alloc_veg_parameter_type(var,landunits,model_structure)
    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % canst1(landunits) )
    ALLOCATE ( var % deciduous(landunits) ) 
    ALLOCATE ( var % dleaf(landunits) )
    ALLOCATE ( var % ejmax(landunits) )
    ALLOCATE ( var % extkn(landunits) )  
    ALLOCATE ( var % frac4(landunits) )
    ALLOCATE ( var % froot(landunits,ms) )
    ALLOCATE ( var % hc(landunits) )
    ALLOCATE ( var % iveg(landunits) )
    ALLOCATE ( var % meth(landunits) )
    ALLOCATE ( var % rp20(landunits) )
    ALLOCATE ( var % rpcoef(landunits) )
    ALLOCATE ( var % shelrb(landunits) )
    ALLOCATE ( var % tmaxvj(landunits) )
    ALLOCATE ( var % tminvj(landunits) )
    ALLOCATE ( var % vbeta(landunits) )
    ALLOCATE ( var % vcmax(landunits) )
    ALLOCATE ( var % vegcf(landunits) )  
    ALLOCATE ( var % vlai(landunits) )
    ALLOCATE ( var % wai(landunits) )  
    ALLOCATE ( var % xfang(landunits) )
    ! Allocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       ALLOCATE ( var % rootbeta(landunits) )
       ALLOCATE ( var % gamma(landunits) )  
       ALLOCATE ( var % F10(landunits) )  
       ALLOCATE ( var % ZR(landunits) )  
    END IF
    ! Allocate variables for canopy_vh model:
    IF(ANY(model_structure%canopy=='canopy_vh')) THEN
      IF(.NOT.ASSOCIATED(var % gamma)) ALLOCATE ( var % gamma(landunits) ) 
    END IF
  END SUBROUTINE alloc_veg_parameter_type

  SUBROUTINE alloc_canopy_type(var,landunits,model_structure)
    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % cansto(landunits) )
    ALLOCATE ( var % delwc(landunits) )
    ALLOCATE ( var % dewmm(landunits) )
    ALLOCATE ( var % fe(landunits) )
    ALLOCATE ( var % fh(landunits) )
    ALLOCATE ( var % fpn(landunits) )
    ALLOCATE ( var % frp(landunits) )
    ALLOCATE ( var % frpw(landunits) )
    ALLOCATE ( var % frpr(landunits) )
    ALLOCATE ( var % frs(landunits) )
    ALLOCATE ( var % fnee(landunits) )
    ALLOCATE ( var % frday(landunits) )
    ALLOCATE ( var % fnv(landunits) )
    ALLOCATE ( var % fev(landunits) )
    ALLOCATE ( var % fevc(landunits) )
    ALLOCATE ( var % fevw(landunits) )
    ALLOCATE ( var % fhv(landunits) )
    ALLOCATE ( var % fhvw(landunits) )
    ALLOCATE ( var % fns(landunits) )
    ALLOCATE ( var % fes(landunits) )
    ALLOCATE ( var % fhs(landunits) )
    ALLOCATE ( var % vlaiw(landunits) )
    ALLOCATE ( var % fwet(landunits) )
    ALLOCATE ( var % tv(landunits) )
    ALLOCATE ( var % ga(landunits) )
    ALLOCATE ( var % ghflux(landunits) )
    ALLOCATE ( var % segg(landunits) )
    ALLOCATE ( var % sghflux(landunits) )
    ALLOCATE ( var % dgdtg(landunits) )
    ALLOCATE ( var % through(landunits) )
    ALLOCATE ( var % precis(landunits) )
    ALLOCATE ( var % rnet(landunits) )
    ALLOCATE ( var % spill(landunits) )
    ALLOCATE ( var % wcint(landunits) )
    ALLOCATE ( var % us(landunits) )
    ALLOCATE ( var % tscrn(landunits) )
    ALLOCATE ( var % qscrn(landunits) )
    ALLOCATE ( var % uscrn(landunits) )
    ALLOCATE ( var % cduv(landunits) )
    ALLOCATE ( var % potev_c(landunits) )
    ALLOCATE ( var % gw(landunits,mf) )
    ALLOCATE ( var % ancj(landunits,mf,3) )
    ALLOCATE ( var % ci(landunits,mf,3) )
    ALLOCATE ( var % gswx(landunits,mf) )
    ALLOCATE ( var % ecy(landunits,mf) )
    ALLOCATE ( var % ecx(landunits,mf) )
    ALLOCATE ( var % tlfy(landunits,mf) )
    ALLOCATE ( var % fwsoil(landunits) )
  END SUBROUTINE alloc_canopy_type

  SUBROUTINE alloc_radiation_type(var,landunits,model_structure)
    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % albedo(landunits,nrb) )
    ALLOCATE ( var % extkb(landunits) )
    ALLOCATE ( var % extkd2(landunits) )
    ALLOCATE ( var % extkd(landunits) )
    ALLOCATE ( var % flws(landunits) )
    ALLOCATE ( var % fvlai(landunits,mf) )
    ALLOCATE ( var % gradis(landunits,mf) )
    ALLOCATE ( var % latitude(landunits) )
    ALLOCATE ( var % lwabv(landunits) )
    ALLOCATE ( var % qcan(landunits,mf,nrb) )
    ALLOCATE ( var % qssabs(landunits) )
    ALLOCATE ( var % rhocdf(landunits,nrb) )
    ALLOCATE ( var % rniso(landunits,mf) )
    ALLOCATE ( var % scalex(landunits,mf) )
    ALLOCATE ( var % transd(landunits) )
    ALLOCATE ( var % trad(landunits) )
    ! new types, ypw 11/july/2008
    ALLOCATE ( var % reffdf(landunits,nrb) )
    ALLOCATE ( var % reffbm(landunits,nrb) )
    ALLOCATE ( var % extkbm(landunits,nrb) )
    ALLOCATE ( var % extkdm(landunits,nrb) )
    ALLOCATE ( var % fbeam(landunits) )
    ALLOCATE ( var % cexpkbm(landunits,nrb) )
    ALLOCATE ( var % cexpkdm(landunits,nrb) )
    ! ********** ypw 11/july/2008
  END SUBROUTINE alloc_radiation_type

  SUBROUTINE alloc_roughness_type(var,landunits,model_structure)
    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % coexp(landunits) )
    ALLOCATE ( var % disp(landunits) )
    ALLOCATE ( var % hruff(landunits) )
    ALLOCATE ( var % rt0us(landunits) )
    ALLOCATE ( var % rt1usa(landunits) )
    ALLOCATE ( var % rt1usb(landunits) )
    ALLOCATE ( var % rt1(landunits) )
    ALLOCATE ( var % term2(landunits) )
    ALLOCATE ( var % term3(landunits) )
    ALLOCATE ( var % term5(landunits) )
    ALLOCATE ( var % term6(landunits) )
    ALLOCATE ( var % usuh(landunits) )
    ALLOCATE ( var % za(landunits) )
    ALLOCATE ( var % z0m(landunits) )
    ALLOCATE ( var % zref(landunits) )
    ALLOCATE ( var % zruffs(landunits) )
    ALLOCATE ( var % z0soilsn(landunits) )
    ALLOCATE ( var % z0soil(landunits) )
    ALLOCATE ( var % zref_uv(landunits) )
    ALLOCATE ( var % zref_tq(landunits) )
    ALLOCATE ( var % za_uv(landunits) )
    ALLOCATE ( var % za_tq(landunits) )
  END SUBROUTINE alloc_roughness_type

  SUBROUTINE alloc_air_type(var,landunits,model_structure)
    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % rho(landunits) )
    ALLOCATE ( var % volm(landunits) )
    ALLOCATE ( var % rlam(landunits) )
    ALLOCATE ( var % qsat(landunits) )
    ALLOCATE ( var % epsi(landunits) )
    ALLOCATE ( var % visc(landunits) )
    ALLOCATE ( var % psyc(landunits) )
    ALLOCATE ( var % dsatdk(landunits) )
    ALLOCATE ( var % cmolar(landunits) )
  END SUBROUTINE alloc_air_type

  SUBROUTINE alloc_met_type(var,landunits,model_structure)
    TYPE(met_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % ca(landunits) )
    ALLOCATE ( var % year(landunits) )
    ALLOCATE ( var % moy(landunits) )
    ALLOCATE ( var % doy(landunits) )
    ALLOCATE ( var % hod(landunits) )
    ALLOCATE ( var % fsd(landunits) )
    ALLOCATE ( var % fld(landunits) )
    ALLOCATE ( var % precip(landunits) )
    ALLOCATE ( var % precip_s(landunits) )
    ALLOCATE ( var % tc(landunits) )
    ALLOCATE ( var % tk(landunits) )
    ALLOCATE ( var % tvair(landunits) )
    ALLOCATE ( var % tvrad(landunits) )
    ALLOCATE ( var % pmb(landunits) )
    ALLOCATE ( var % ua(landunits) )
    ALLOCATE ( var % qv(landunits) )
    ALLOCATE ( var % qvair(landunits) )
    ALLOCATE ( var % da(landunits) )
    ALLOCATE ( var % dva(landunits) )
    ALLOCATE ( var % coszen(landunits) )
    ! Deallocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       ALLOCATE ( var % tk_old(landunits) )
    END IF
  END SUBROUTINE alloc_met_type

  SUBROUTINE alloc_sum_flux_type(var,landunits,model_structure)
    TYPE(sum_flux_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % sumpn(landunits) )
    ALLOCATE ( var % sumrp(landunits) )
    ALLOCATE ( var % sumrpw(landunits) )
    ALLOCATE ( var % sumrpr(landunits) )
    ALLOCATE ( var % sumrs(landunits) )
    ALLOCATE ( var % sumrd(landunits) )
    ALLOCATE ( var % dsumpn(landunits) )
    ALLOCATE ( var % dsumrp(landunits) )
    ALLOCATE ( var % dsumrs(landunits) )
    ALLOCATE ( var % dsumrd(landunits) )
    ALLOCATE ( var % sumxrp(landunits) )
    ALLOCATE ( var % sumxrs(landunits) )
  END SUBROUTINE alloc_sum_flux_type

  SUBROUTINE alloc_bgc_pool_type(var,landunits,model_structure)
    TYPE(bgc_pool_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    ALLOCATE ( var % cplant(landunits,ncp) )
    ALLOCATE ( var % csoil(landunits,ncs) )
    ! Allocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       ALLOCATE ( var % clitt(landunits) )
    END IF
  END SUBROUTINE alloc_bgc_pool_type

  SUBROUTINE alloc_model_structure_type(var,landunits)
    TYPE(model_structure_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var%canopy(landunits) )
    ALLOCATE ( var%soil(landunits) )
    ALLOCATE ( var%sli_isotope(landunits) )
    ALLOCATE ( var%photosynthesis(landunits) )
    ALLOCATE ( var%sli_litter(landunits) )
    ALLOCATE ( var%sli_coupled(landunits) )
  END SUBROUTINE alloc_model_structure_type

  ! Begin deallocation routines:
  SUBROUTINE dealloc_balances_type(var,landunits,model_structure)
    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % canopy_drybal )
    NULLIFY ( var % evap_tot )
    NULLIFY ( var % osnowd0 )
    NULLIFY ( var % precip_tot )
    NULLIFY ( var % rnoff_tot )
    NULLIFY ( var % wbtot0 )
    NULLIFY ( var % canopy_wetbal )
	NULLIFY ( var % Radbal)
	NULLIFY ( var % RadbalSum)
	NULLIFY ( var % Ebal)
	NULLIFY ( var % EbalSum)
	NULLIFY ( var % EbalSoil)
	NULLIFY ( var % EbalSnow)
	NULLIFY ( var % EbalVeg)
	NULLIFY ( var % Wbal)
	NULLIFY ( var % WbalSum)
	NULLIFY ( var % WbalSoil)
	NULLIFY ( var % WbalVeg)
	NULLIFY ( var % WbalSnow)
  END SUBROUTINE dealloc_balances_type
  
  SUBROUTINE dealloc_soil_parameter_type(var,landunits,model_structure)
    TYPE(soil_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % albsoil )
    NULLIFY ( var % bch )
    NULLIFY ( var % c3 )
    NULLIFY ( var % clay )
    NULLIFY ( var % cnsd )
    NULLIFY ( var % css )
    NULLIFY ( var % hsbh )
    NULLIFY ( var % hyds )
    NULLIFY ( var % i2bp3 )
    NULLIFY ( var % ibp2 )
    NULLIFY ( var % isoilm )
    NULLIFY ( var % rhosoil )
    NULLIFY ( var % rs20 )
    NULLIFY ( var % sand )
    NULLIFY ( var % sfc )
    NULLIFY ( var % silt )
    NULLIFY ( var % ssat )
    NULLIFY ( var % sucs )
    NULLIFY ( var % swilt )
    ! Deallocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       NULLIFY ( var % nhorizons)
       NULLIFY ( var % bchB )
       NULLIFY ( var % clayB )
       NULLIFY ( var % cnsdB)
       NULLIFY ( var % cssB)
       NULLIFY ( var % hydsB)
       NULLIFY ( var % rhosoilB)
       NULLIFY ( var % sandB)
       NULLIFY ( var % siltB)
       NULLIFY ( var % ssatB)
       NULLIFY ( var % sucsB )
       NULLIFY ( var % swiltB )
       NULLIFY ( var % depthA)
       NULLIFY ( var % depthB)
       NULLIFY ( var % sfcB)
       NULLIFY ( var % clitt )
       NULLIFY ( var % swilt_vec )
    END IF
    IF(ANY(model_structure%canopy=='canopy_vh')) THEN
       IF(ASSOCIATED(var % swilt_vec)) NULLIFY ( var % swilt_vec )
    END IF
  END SUBROUTINE dealloc_soil_parameter_type

  SUBROUTINE dealloc_soil_snow_type(var,landunits,model_structure)
    TYPE(soil_snow_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % albsoilsn )
    NULLIFY ( var % cls )
    NULLIFY ( var % dfn_dtg )
    NULLIFY ( var % dfh_dtg )
    NULLIFY ( var % dfe_ddq )
    NULLIFY ( var % ddq_dtg )
    NULLIFY ( var % evapsn )
    NULLIFY ( var % fwtop )
    NULLIFY ( var % gammzz )
    NULLIFY ( var % isflag )
    NULLIFY ( var % osnowd )
    NULLIFY ( var % potev )
    NULLIFY ( var % pwb_min )
    NULLIFY ( var % runoff )
    NULLIFY ( var % rnof1 )
    NULLIFY ( var % rnof2 )
    NULLIFY ( var % rtsoil )
    NULLIFY ( var % sconds )
    NULLIFY ( var % sdepth )
    NULLIFY ( var % smass )
    NULLIFY ( var % snage )
    NULLIFY ( var % snowd )
    NULLIFY ( var % smelt )
    NULLIFY ( var % ssdn )
    NULLIFY ( var % ssdnn )
    NULLIFY ( var % tgg )
    NULLIFY ( var % tggsn )
    NULLIFY ( var % tss )
    NULLIFY ( var % tssold )
    NULLIFY ( var % wb )
    NULLIFY ( var % wbfice )
    NULLIFY ( var % wbice )
    NULLIFY ( var % wblf )
    NULLIFY ( var % wbtot )
    NULLIFY ( var % evap)
    NULLIFY ( var % delwcol )
    ! Deallocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       NULLIFY ( var % S )	
       NULLIFY ( var % SL )	
       NULLIFY ( var % TL )	
       NULLIFY ( var % h0)	
       NULLIFY ( var % rex )	
       NULLIFY ( var % wflux )
       NULLIFY ( var % hflux )		
       NULLIFY ( var % zdelta )	
       NULLIFY ( var % kth )	
       NULLIFY ( var % Tsurface )
       NULLIFY ( var % rh0 )	
       NULLIFY ( var % rhsurface )	
       NULLIFY ( var % lE )		
       NULLIFY ( var % ciso )	
       NULLIFY ( var % cisoL )	
       NULLIFY ( var % SA )	
       NULLIFY ( var % SB )	
       NULLIFY ( var % delwcolA )	
       NULLIFY ( var % delwcolB )	
       NULLIFY ( var % rexA )
       NULLIFY ( var % rexB )	
       NULLIFY ( var % leachAB ) 
       NULLIFY ( var % rlitt )
       NULLIFY ( var % thetai )
    END IF
  END SUBROUTINE dealloc_soil_snow_type

  SUBROUTINE dealloc_veg_parameter_type(var,landunits,model_structure)
    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % iveg )
    NULLIFY ( var % meth )
    NULLIFY ( var % vlai )
    NULLIFY ( var % froot )
    NULLIFY ( var % canst1 )
    NULLIFY ( var % ejmax )
    NULLIFY ( var % frac4 )
    NULLIFY ( var % wai )  ! new addition in Oct 2007 (YP)
    NULLIFY ( var % vegcf )  ! new addition in Oct 2007 (YP)
    NULLIFY ( var % tminvj )
    NULLIFY ( var % tmaxvj )
    NULLIFY ( var % vbeta )
    NULLIFY ( var % hc )
    NULLIFY ( var % shelrb )
    NULLIFY ( var % vcmax )
    NULLIFY ( var % xfang )
    NULLIFY ( var % dleaf )
    NULLIFY ( var % rp20 )
    NULLIFY ( var % rpcoef )
    NULLIFY ( var % extkn )  ! new addition in Oct 2007 (YP)
    NULLIFY ( var % deciduous )  ! rml addition 22/10/07
    ! Deallocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       NULLIFY ( var % rootbeta )
       NULLIFY ( var % gamma ) ! vh 20/07/09
       NULLIFY ( var % F10 )
       NULLIFY ( var % ZR )
    END IF
    ! Deallocate variables for canopy_vh model:
    IF(ANY(model_structure%canopy=='canopy_vh')) THEN
      IF(ASSOCIATED(var % gamma)) NULLIFY ( var % gamma )  
    END IF

  END SUBROUTINE dealloc_veg_parameter_type

  SUBROUTINE dealloc_canopy_type(var,landunits,model_structure)
    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % cansto )
    NULLIFY ( var % delwc )
    NULLIFY ( var % dewmm )
    NULLIFY ( var % fe )
    NULLIFY ( var % fh )
    NULLIFY ( var % fpn )
    NULLIFY ( var % frp )
    NULLIFY ( var % frpw )
    NULLIFY ( var % frpr )
    NULLIFY ( var % frs )
    NULLIFY ( var % fnee )
    NULLIFY ( var % frday )
    NULLIFY ( var % fnv )
    NULLIFY ( var % fev )
    NULLIFY ( var % fevc )
    NULLIFY ( var % fevw )
    NULLIFY ( var % fhv )
    NULLIFY ( var % fhvw )
    NULLIFY ( var % fns )
    NULLIFY ( var % fes )
    NULLIFY ( var % fhs )
    NULLIFY ( var % vlaiw )
    NULLIFY ( var % fwet )
    NULLIFY ( var % tv )
    NULLIFY ( var % ga )
    NULLIFY ( var % ghflux )
    NULLIFY ( var % segg )
    NULLIFY ( var % sghflux )
    NULLIFY ( var % dgdtg )
    NULLIFY ( var % through )
    NULLIFY ( var % precis )
    NULLIFY ( var % rnet )
    NULLIFY ( var % spill )
    NULLIFY ( var % wcint )
    NULLIFY ( var % us )
    NULLIFY ( var % tscrn )
    NULLIFY ( var % qscrn )
    NULLIFY ( var % uscrn )
    NULLIFY ( var % cduv )
    NULLIFY ( var % potev_c )
    NULLIFY ( var % gw )
    NULLIFY ( var % ancj )
    NULLIFY ( var % gswx )
    NULLIFY ( var % tlfy )
    NULLIFY ( var % ecy )
    NULLIFY ( var % ecx )
    NULLIFY ( var % ci )
    NULLIFY ( var % fwsoil )
  END SUBROUTINE dealloc_canopy_type

  SUBROUTINE dealloc_radiation_type(var,landunits,model_structure)
    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % albedo )
    NULLIFY ( var % extkb )
    NULLIFY ( var % extkd2 )
    NULLIFY ( var % extkd )
    NULLIFY ( var % flws )
    NULLIFY ( var % fvlai )
    NULLIFY ( var % gradis )
    NULLIFY ( var % latitude )
    NULLIFY ( var % lwabv )
    NULLIFY ( var % qcan )
    NULLIFY ( var % qssabs )
    NULLIFY ( var % rhocdf )
    NULLIFY ( var % rniso )
    NULLIFY ( var % scalex )
    NULLIFY ( var % transd )
    NULLIFY ( var % trad )
    ! new types, ypw 11/july/2008
    NULLIFY ( var % reffdf )
    NULLIFY ( var % reffbm )
    NULLIFY ( var % extkbm )
    NULLIFY ( var % extkdm )
    NULLIFY ( var % fbeam )
    NULLIFY ( var % cexpkbm )
    NULLIFY ( var % cexpkdm )
    ! ********** ypw 11/july/2008
  END SUBROUTINE dealloc_radiation_type

  SUBROUTINE dealloc_roughness_type(var,landunits,model_structure)
    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % coexp )
    NULLIFY ( var % disp )
    NULLIFY ( var % hruff )
    NULLIFY ( var % rt0us )
    NULLIFY ( var % rt1usa )
    NULLIFY ( var % rt1usb )
    NULLIFY ( var % rt1 )
    NULLIFY ( var % term2 )
    NULLIFY ( var % term3 )
    NULLIFY ( var % term5 )
    NULLIFY ( var % term6 )
    NULLIFY ( var % usuh )
    NULLIFY ( var % za )
    NULLIFY ( var % z0m )
    NULLIFY ( var % zref )
    NULLIFY ( var % zruffs )
    NULLIFY ( var % z0soilsn )
    NULLIFY ( var % z0soil )
  END SUBROUTINE dealloc_roughness_type

  SUBROUTINE dealloc_air_type(var,landunits,model_structure)
    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % rho )
    NULLIFY ( var % volm )
    NULLIFY ( var % rlam )
    NULLIFY ( var % qsat )
    NULLIFY ( var % epsi )
    NULLIFY ( var % visc )
    NULLIFY ( var % psyc )
    NULLIFY ( var % dsatdk )
    NULLIFY ( var % cmolar )
  END SUBROUTINE dealloc_air_type

  SUBROUTINE dealloc_met_type(var,landunits,model_structure)
    TYPE(met_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % ca )
    NULLIFY ( var % year )
    NULLIFY ( var % moy )
    NULLIFY ( var % doy )
    NULLIFY ( var % hod )
    NULLIFY ( var % fsd )
    NULLIFY ( var % fld )
    NULLIFY ( var % precip )
    NULLIFY ( var % precip_s )
    NULLIFY ( var % tc )
    NULLIFY ( var % tk )
    NULLIFY ( var % tvair )
    NULLIFY ( var % tvrad )
    NULLIFY ( var % pmb )
    NULLIFY ( var % ua )
    NULLIFY ( var % qv )
    NULLIFY ( var % qvair )
    NULLIFY ( var % da )
    NULLIFY ( var % dva )
    NULLIFY ( var % coszen )
    ! Deallocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       NULLIFY ( var % tk_old )	
    END IF
  END SUBROUTINE dealloc_met_type

  SUBROUTINE dealloc_sum_flux_type(var,landunits,model_structure)
    TYPE(sum_flux_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % sumpn )
    NULLIFY ( var % sumrp )
    NULLIFY ( var % sumrpw )
    NULLIFY ( var % sumrpr )
    NULLIFY ( var % sumrs )
    NULLIFY ( var % sumrd )
    NULLIFY ( var % dsumpn )
    NULLIFY ( var % dsumrp )
    NULLIFY ( var % dsumrs )
    NULLIFY ( var % dsumrd )
    NULLIFY ( var % sumxrp )
    NULLIFY ( var % sumxrs )
  END SUBROUTINE dealloc_sum_flux_type

  SUBROUTINE dealloc_bgc_pool_type(var,landunits,model_structure)
    TYPE(bgc_pool_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    TYPE(model_structure_type), INTENT(IN) :: model_structure
    NULLIFY ( var % cplant )
    NULLIFY ( var % csoil )
    ! Deallocate variables for SLI soil model:
    IF(ANY(model_structure%soil=='sli')) THEN
       NULLIFY ( var % clitt )	
    END IF
  END SUBROUTINE dealloc_bgc_pool_type

  SUBROUTINE dealloc_model_structure_type(var,landunits)
    TYPE(model_structure_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    NULLIFY ( var%canopy )
    NULLIFY ( var%soil )
    NULLIFY ( var%sli_isotope )
    NULLIFY ( var%photosynthesis )
    NULLIFY ( var%sli_litter )
    NULLIFY ( var%sli_coupled )
  END SUBROUTINE dealloc_model_structure_type

END MODULE cable_types
