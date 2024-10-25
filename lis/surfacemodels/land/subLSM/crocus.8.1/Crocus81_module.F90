!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module Crocus81_module
!BOP
!
! !MODULE: Crocus81_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the Crocus81 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[n]
!     nest id. unit: -
!   \item[nsnow]
!     number of snow layer. unit: 
!   \item[nimpur]
!     number of impurtites. unit: 
!   \item[year]
!     year of the currrent time step. unit: -
!   \item[month]
!     month of the current time step. unit: -
!   \item[day]
!     day of the current time step. unit: -
!   \item[hour]
!     hour of the current time step. unit: -
!   \item[minute]
!     minute of the current time step. unit: -
!   \item[SNOWRES\_opt]
!     SNOWRES\_opt  = ISBA-SNOW3L turbulant exchange option
!   'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!   'RIL' = Limit Richarson number under very stable
! conditions (currently testing)
!   'M98'  = Martin et Lejeune 1998 : older computation for turbulent fluxes coefficents in Crocus. unit: 
!   \item[OMEB\_BOOL]
!     !True = coupled to MEB. surface fluxes are IMPOSED
! as an upper boundary condition to the explicit snow schemes. 
!If = False, then energy budget and fluxes are computed herein.. unit: -
!   \item[GLACIER\_BOOL]
!     !True = Over permanent snow and ice, initialise WGI=WSAT, Hsnow$>$=10m and allow 0.8$<$SNOALB$<$0.85
! False = No specific treatment. unit: -
!   \item[HIMPLICIT\_WIND\_opt]
!     wind implicitation option  'OLD' = direct , 'NEW' = Taylor serie, order 1. unit: -
!   \item[SNOWSWE]
!     Snow layer(s) liquid Water Equivalent (SWE:kg m-2). unit: kg/m2
!   \item[SNOWRHO]
!     Snow layer(s) averaged density (kg/m3). unit: kg/m3
!   \item[SNOWHEAT]
!     Snow layer(s) heat content (J/m2). unit: J/m2
!   \item[SNOWALB]
!     snow surface albedo. unit: -
!   \item[SNOWGRAN1]
!     Snow layers grain feature 1. unit: -
!   \item[SNOWGRAN2]
!     Snow layer grain feature 2. unit: -
!   \item[SNOWHIST]
!     Snow layer grain historical parameter (only for non dendritic snow) (-) in {0-5}. unit: -
!   \item[SNOWAGE]
!     Age since snowfall (day). unit: day
!   \item[PTSTEP]
!     time step of the integration. unit: s
!   \item[PPS]
!     pressure at atmospheric model surface (Pa). unit: Pa
!   \item[SRSNOW]
!     snow rate (SWE) [kg/(m2 s)]. unit: kg/(m2 s)
!   \item[RRSNOW]
!     rain rate [kg/(m2 s)]. unit: kg/(m2 s)
!   \item[TA]
!     atmospheric temperature at level za (K). unit: K
!   \item[TG]
!     Surface soil temperature (effective temperature the of layer lying below snow) (K)  (for snowcro.F90 we only use the surface layer ZP\_TG(:,1))  (#nsoil depends on 2-L, 3-L DIF). unit: K
!   \item[SW\_RAD]
!     incoming solar radiation (W/m2). unit: W/m2
!   \item[QA]
!     atmospheric specific humidity at level za. unit: -
!   \item[Wind\_E]
!     Eastward Wind. unit: m/s
!   \item[Wind\_N]
!     Northward Wind. unit: m/s
!   \item[LW\_RAD]
!     atmospheric infrared radiation (W/m2). unit: W/m2
!   \item[UREF]
!     reference height of the wind. unit: m
!   \item[SLOPE]
!     angle between the normal to the surface and the vertical  (MN:  replaced PDIRCOSZW with slope and computed the cosine in the driver). unit: -
!   \item[ZREF]
!     Reference height of the first atmospheric level (m). unit: m
!   \item[Z0NAT]
!     grid box average roughness length (m) (roughness length for momentum). unit: m
!   \item[Z0EFF]
!     roughness length for momentum (modd\_diagn.F90 effective roughness length for heat(!?)). unit: m
!   \item[Z0HNAT]
!     grid box average roughness length for heat. unit: m
!   \item[ALB]
!     soil/vegetation albedo (monthly values). unit: -
!   \item[SOILCOND]
!     soil thermal conductivity (W m-1 K-1). unit: W /(m K)
!   \item[D\_G]
!     !Assumed first soil layer thickness (m)
!Used to calculate ground/snow heat flux   (D\_G(:,1)). unit: m
!   \item[SNOWLIQ]
!     Snow layer(s) liquid water content (m). unit: m
!   \item[SNOWTEMP]
!     Snow layer(s) temperature (K). unit: K
!   \item[SNOWDZ]
!     Snow layer(s) thickness (m). unit: m
!   \item[THRUFAL]
!     Rate that liquid water leaves snow pack: paritioned into soil infiltration/runoff  by ISBA [kg/(m2 s)]. unit: kg/(m2 s)
!   \item[GRNDFLUX]
!     Soil/snow interface heat flux (W/m2). unit: W/m2
!   \item[SNDRIFT]
!     Blowing snow sublimation (kg/m2/s) NOTE: Snow compaction and metamorphism due to drift, Mass is unchanged  (Assistance #1592). unit: kg/m2/s
!   \item[RI\_n]
!     Richardson number (-)  NOTE: RI has not been initialized in CALL\_MODEL (If not OMED initalized to undefined in the snow3L\_isba.F90). unit: -
!   \item[EMISNOW]
!     Snow surface emissivity (initialize to XEMISSN
!XEMISSN  comes from MODD\_SNOW\_PAR). unit: -
!   \item[CDSNOW]
!     Drag coefficient for momentum over snow (-). unit: -
!   \item[USTARSNOW]
!     Friction velocity over snow (m/s);. unit: m/s
!   \item[CHSNOW]
!     Drag coefficient for heat over snow  (-). unit: -
!   \item[SNOWHMASS]
!     Heat content change due to mass changes in snowpack: for budget calculations only.. unit: J/m2
!   \item[QS]
!     Surface humidity (kg/kg). unit: kg/kg
!   \item[PERMSNOWFRAC]
!     Fraction of permanet snow/ice. unit: -
!   \item[LAT]
!     Latitude in decimal degree  (latitude (degrees +North)). unit: degrees
!   \item[LON]
!     Longitude in decimal year    (longitude (degrees +East)). unit: degrees
!   \item[SNOWDRIFT\_opt]
!     Mechanical transformation of snow grain and compaction + effect of wind on falling snow properties
!	'NONE': No snowdrift scheme
!	'DFLT': falling snow falls as purely dendritic
! 	'GA01': Gallee et al 2001
!	'VI13': Vionnet et al 2013. unit: -
!   \item[SNOWDRIFT\_SUBLIM\_BOOL]
!     Logicals for snowdrift sublimation. unit: -
!   \item[SNOW\_ABS\_ZENITH\_BOOL]
!     Activate parametrization of solar absorption for polar regions (If True modify solar absorption as a function of solar zenithal angle). unit: -
!   \item[SNOWMETAMO\_opt]
!     Metamorphism scheme: B92 (historical version, Brun et al 92), C13, T07, F06 (see Carmagnola et al 2014). unit: -
!   \item[SNOWRAD\_opt]
!     Radiative transfer scheme. HSNOWRAD=B92 Brun et al 1992.  HSNOWRAD=T17 (Tuzet et al. 2017) (Libois et al. 2013) TARTES with impurities content scheme. unit: -
!   \item[ATMORAD\_BOOL]
!     Activate atmotartes scheme  (default=.FALSE. # This option is not stable yet, but it is supposed to compute the direct/diffuse ratio directly from atmospheric informations (AOD, Ozone column, Water column..)). unit: -
!   \item[IMPWET]
!     [init\_surf\_atmn.F90   --$>$ wet deposit coefficient for each impurity typeÊ (g/m_/s)  ,   snowcro.F90 --$>$ Dry and wet deposit coefficient from Forcing File(g/m_/s)   (# 1553 You can either feed the model with prescribed and constant deposition fluxes or introduce a wet and dry deposition field directly in the forcing file. ). unit: g/m_/s
!   \item[IMPDRY]
!     [init\_surf\_atmn.F90   --$>$ wet deposit coefficient for each impurity typeÊ (g)  ,   snowcro.F90 --$>$ Dry and wet deposit coefficient from Forcing File(g/m_/s). unit: g/m_/s
!   \item[SNOWFALL\_opt]
!     New options for multiphysics version (Cluzet et al 2016). Falling snow scheme: V12 (Vionnet et al. 2012) , A76 (Anderson 1976), S02 (Lehning and al. 2002), P75 (Pahaut 1975). unit: -
!   \item[SNOWCOND\_opt]
!     Thermal conductivity scheme: Y81 (Yen 1981), I02 (Boone et al. 2002) C11 (Calonne et al. 2011). unit: -
!   \item[SNOWHOLD\_opt]
!     liquid water content scheme: B92 (Brun et al. 1992) O04 (Oleson et al., 2004) S02 (SNOWPACK, Lehning et al, 2002) B02 (ISBA\_ES, Boone et al. 2002). unit: -
!   \item[SNOWCOMP\_opt]
!     B92 snow compaction basis version and B93 for slightly different parameters   (NOTE: IN THE CODE S14, T11, ). unit: -
!   \item[SNOWZREF\_opt]
!     reference height is constant or variable from the snow surface: CST (constant from snow surface, i.e. Col de Porte) or VAR (variable from snow surface = snow depth has to be removed from reference height). unit: -
!   \item[SNOWMAK\_dz]
!     Snowmaking thickness (m). unit: m
!   \item[SNOWCOMPACT\_BOOL]
!     Snowmaking and Grooming options. unit: -
!   \item[SNOWMAK\_BOOL]
!     Snowmaking and Grooming options. unit: -
!   \item[SNOWTILLER\_BOOL]
!     Snowmaking and Grooming options. unit: -
!   \item[SELF\_PROD\_BOOL]
!     Snowmaking and Grooming options. unit: -
!   \item[SNOWMAK\_PROP\_BOOL]
!     Snowmaking and Grooming options. unit: -
!   \item[PRODSNOWMAK\_BOOL]
!     Snowmaking and Grooming options. unit: -
!   \item[SLOPE\_DIR]
!     !Typical slope aspect in the grid  (deg from N clockwise). unit: degrees
!   \item[SAND]
!     Soil sand fraction (-). unit: -
!   \item[SILT]
!     Soil silt fraction (-). unit: -
!   \item[CLAY]
!     Soil clay fraction (-). unit: -
!   \item[POROSITY]
!     Soil porosity (m3 m-3). unit: m3/m3
!   \item[SD\_1D]
!     Total snow depth, temporally added. unit: m
!   \item[SWE\_1D]
!     Total SWE, temporally added. unit: kg/m2
!   \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  10/18/19: Mahdi Navari, Shugong Wang Initial implementation for LIS 7 and Crocus81
!
!EOP
    implicit none
    private
    type, public :: crocus81dec
        !-------------------------------------------------------------------------
        ! forcing
        !-------------------------------------------------------------------------
        real               :: PPS
        real               :: SRSNOW
        real               :: RRSNOW
        real               :: TA
        real               :: SW_RAD
        real               :: QA
        real               :: Wind_E
        real               :: Wind_N
        real               :: LW_RAD
        !-------------------------------------------------------------------------
        ! spatial parameter
        !-------------------------------------------------------------------------
        !logical            :: GLACIER_BOOL ! will be computed in the driver using PERMSNOWFRAC
        real               :: SLOPE
        real               :: SOILCOND ! will be computed in the drvier using soil parameters
        real               :: PERMSNOWFRAC
        real               :: SLOPE_DIR
        real               :: SAND
        real               :: SILT
        real               :: CLAY
        real               :: POROSITY
        real               :: XWGI
        real               :: XWG
        !-------------------------------------------------------------------------
        ! multilevel spatial parameter
        !-------------------------------------------------------------------------
        real               :: TG ! for now assume there is no energy exchange in the soil 
                                     ! snow interface and read the variable from lis.config
        real, pointer      :: ALB(:) 
        !-------------------------------------------------------------------------
        ! state
        !-------------------------------------------------------------------------
        real, pointer      :: SNOWSWE(:)
        real, pointer      :: SNOWRHO(:)
        real, pointer      :: SNOWHEAT(:)
        real               :: SNOWALB
        real, pointer      :: SNOWGRAN1(:)
        real, pointer      :: SNOWGRAN2(:)
        real, pointer      :: SNOWHIST(:)
        real, pointer      :: SNOWAGE(:)
        real, pointer      :: SNOWLIQ(:)
        real, pointer      :: SNOWTEMP(:)
        real, pointer      :: SNOWDZ(:)
        real               :: GRNDFLUX
        real               :: SNDRIFT
        real               :: RI_n
        real               :: CDSNOW
        real               :: USTARSNOW
        real               :: CHSNOW
        real               :: SNOWMAK_dz
        !-------------------------------------------------------------------------
        ! output
        !-------------------------------------------------------------------------
        real               :: THRUFAL
        real               :: EMISNOW
        real               :: SNOWHMASS
        real               :: QS
        real               :: SD_1D
        real               :: SWE_1D
    end type crocus81dec
end module Crocus81_module
