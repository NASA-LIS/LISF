!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module NoahMP36_module
!BOP
!
! !MODULE: NoahMP36_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the NoahMP36 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[landuse\_tbl\_name]
!     Noah model landuse parameter table. unit: -
!   \item[soil\_tbl\_name]
!     Noah model soil parameter table. unit: -
!   \item[gen\_tbl\_name]
!     Noah model general parameter table. unit: -
!   \item[noahmp\_tbl\_name]
!     NoahMP parameter table. unit: -
!   \item[landuse\_scheme\_name]
!     landuse classification scheme. unit: -
!   \item[soil\_scheme\_name]
!     soil classification scheme. unit: -
!   \item[dveg\_opt]
!     vegetation model. unit: -
!   \item[crs\_opt]
!     canopy stomatal resistance. unit: -
!   \item[btr\_opt]
!     soil moisture factor for stomatal resistance. unit: -
!   \item[run\_opt]
!     runoff and groundwater. unit: -
!   \item[sfc\_opt]
!     surface layer drag coefficients (CH \& CM). unit: -
!   \item[frz\_opt]
!     supercooled liquid water. unit: -
!   \item[inf\_opt]
!     frozen soil permeability. unit: -
!   \item[rad\_opt]
!     radiation transfer. unit: -
!   \item[alb\_opt]
!     snow surface albedo. unit: -
!   \item[snf\_opt]
!     rainfall \& snowfall. unit: -
!   \item[tbot\_opt]
!     lower boundary of soil temperature. unit: -
!   \item[stc\_opt]
!     snow/soil temperature time scheme. unit: -
!   \item[nslcats]
!     the number of total soil types in parameter table. unit: -
!   \item[nlucats]
!     the number of total land cover types in parameter table. unit: -
!   \item[nslpcats]
!     the number of total slope category for Noah baseflow. unit: -
!   \item[latitude]
!     latitude in decimal degree. unit: -
!   \item[longitude]
!     longitude in decimal degree. unit: -
!   \item[year]
!     year of the current time step. unit: -
!   \item[month]
!     month of the current time step. unit: -
!   \item[day]
!     day of the current time step. unit: -
!   \item[hour]
!     hour of the current time step. unit: -
!   \item[minute]
!     minute of the current time step. unit: -
!   \item[dt]
!     time step in seconds. unit: s
!   \item[nsoil]
!     number of soil layers. unit: -
!   \item[sldpth]
!     thickness of soil layers. unit: -
!   \item[nsnow]
!     maximum number of snow layers. unit: -
!   \item[shdfac\_monthly]
!     monthly values for green vegetation fraction. unit: -
!   \item[vegetype]
!     land cover type index. unit: -
!   \item[soiltype]
!     soil type index. unit: -
!   \item[slopetype]
!     slope type for Noah baseflow. unit: -
!   \item[urban\_vegetype]
!     urban land cover type index. unit: -
!   \item[ice\_flag]
!     ice flag: 0 = no ice, 1 = ice. unit: -
!   \item[st\_flag]
!     surface type 1=soil, 2=lake. unit: -
!   \item[sc\_idx]
!     soil color type. unit: -
!   \item[iz0tlnd]
!     option of Chen adjustment of Czil. unit: -
!   \item[smceq]
!     equilibrium soil water content. unit: m^3 m-3
!   \item[tair]
!     air temperature. unit: K
!   \item[psurf]
!     air pressure. unit: Pa
!   \item[wind\_e]
!     eastward wind speed. unit: m s-1
!   \item[wind\_n]
!     northward wind speed. unit: m s-1
!   \item[qair]
!     near Surface Specific Humidity. unit: kg kg-1
!   \item[swdown]
!     downward solar radiation. unit: w/m2
!   \item[lwdown]
!     downward longwave radiation. unit: w/m2
!   \item[prcp]
!     total precipitation Rate. unit: kg m-2 s-1
!   \item[tbot]
!     deep-layer soil temperature. unit: K
!   \item[pblh]
!     planetary boundary layer height. unit: m
!   \item[zlvl]
!     reference height of temperature and humidity. unit: m
!   \item[albold]
!     snow albedo at last time step. unit: -
!   \item[sneqvo]
!     snow mass at the last time step. unit: mm
!   \item[sstc]
!     snow/soil temperature. unit: K
!   \item[sh2o]
!     volumetric liquid soil moisture. unit: m^3 m-3
!   \item[smc]
!     volumetric soil moisture, ice + liquid. unit: m^3 m-3
!   \item[tah]
!     canopy air temperature. unit: K
!   \item[eah]
!     canopy air vapor pressure. unit: Pa
!   \item[fwet]
!     wetted or snowed fraction of canopy. unit: -
!   \item[canliq]
!     intercepted liquid water. unit: mm
!   \item[canice]
!     intercepted ice mass. unit: mm
!   \item[tv]
!     vegetation temperature. unit: K
!   \item[tg]
!     ground temperature (skin temperature). unit: K
!   \item[qsnow]
!     snowfall on the ground. unit: mm s-1
!   \item[isnow]
!     actual number of snow layers. unit: -
!   \item[zss]
!     snow/soil layer-bottom depth from snow surface. unit: m
!   \item[snowh]
!     snow height. unit: m
!   \item[sneqv]
!     snow water equivalent. unit: mm
!   \item[snowice]
!     snow-layer ice. unit: mm
!   \item[snowliq]
!     snow-layer liquid water. unit: mm
!   \item[zwt]
!     depth to water table. unit: m
!   \item[wa]
!     water storage in aquifer. unit: mm
!   \item[wt]
!     water in aquifer and saturated soil. unit: mm
!   \item[wslake]
!     lake water storage. unit: mm
!   \item[lfmass]
!     leaf mass. unit: g/m2
!   \item[rtmass]
!     mass of fine roots. unit: g/m2
!   \item[stmass]
!     stem mass. unit: g/m2
!   \item[wood]
!     mass of wood including woody roots. unit: g/m2
!   \item[stblcp]
!     stable carbon in deep soil. unit: g/m2
!   \item[fastcp]
!     short-lived carbon in shallow soil. unit: g/m2
!   \item[lai]
!     leaf area index. unit: -
!   \item[sai]
!     stem area index. unit: -
!   \item[cm]
!     momentum drag coefficient. unit: s/m
!   \item[ch]
!     sensible heat exchange coefficient. unit: s/m
!   \item[tauss]
!     snow aging term. unit: -
!   \item[smcwtd]
!     soil water content between bottom of the soil and water table. unit: m^3 m-3
!   \item[deeprech]
!     recharge to or from the water table when deep. unit: m
!   \item[rech]
!     recharge to or from the water table when shallow. unit: m
!   \item[fsa]
!     total absorbed solar radiation. unit: W m-2
!   \item[fsr]
!     total reflected solar radiation. unit: W m-2
!   \item[fira]
!     total net longwave radiation to atmosphere. unit: W m-2
!   \item[fsh]
!     total sensible heat to atmosphere. unit: W m-2
!   \item[ssoil]
!     ground heat flux to soil. unit: W m-2
!   \item[fcev]
!     canopy evaporative heat to atmosphere. unit: W m-2
!   \item[fgev]
!     ground evaporative heat to atmosphere. unit: W m-2
!   \item[fctr]
!     transpiration heat to atmosphere. unit: W m-2
!   \item[ecan]
!     evaporation rate of canopy water. unit: kg m-2 s-1
!   \item[etran]
!     transpiration rate. unit: kg m-2 s-1
!   \item[edir]
!     direct evaporation rate from surface. unit: kg m-2 s-1
!   \item[trad]
!     surface radiative temperature. unit: K
!   \item[tgb]
!     ground temperature. unit: K
!   \item[tgv]
!     ground surface temperature. unit: K
!   \item[t2mv]
!     2-m air temperature over vegetated part. unit: K
!   \item[t2mb]
!     2-m height air temperature. unit: K
!   \item[q2v]
!     2-m specific humidity over vegetation. unit: kg kg-1
!   \item[q2b]
!     2-m air specific humidity. unit: kg kg-1
!   \item[runsrf]
!     surface runoff. unit: kg m-2 s-1
!   \item[runsub]
!     baseflow (saturation excess). unit: kg m-2 s-1
!   \item[apar]
!     photosynthesis active energy by canopy. unit: W m-2
!   \item[psn]
!     total photosynthesis of CO2. unit: umol m-2 s-1
!   \item[sav]
!     solar radiation absorbed by vegetation. unit: W m-2
!   \item[sag]
!     solar radiation absorbed by ground. unit: W m-2
!   \item[fsno]
!     snow-cover fraction on the ground. unit: -
!   \item[nee]
!     net ecosystem exchange of CO2. unit: g/m2s
!   \item[gpp]
!     net instantaneous assimilation of carbon. unit: g/m2s
!   \item[npp]
!     net primary productivity of carbon. unit: g/m2s
!   \item[fveg]
!     green vegetation fraction. unit: -
!   \item[albedo]
!     surface albedo. unit: -
!   \item[qsnbot]
!     melting water out of snow bottom. unit: kg m-2 s-1
!   \item[ponding]
!     surface ponding. unit: mm
!   \item[ponding1]
!     surface ponding1. unit: mm
!   \item[ponding2]
!     surface ponding2. unit: mm
!   \item[rssun]
!     sunlit stomatal resistance. unit: s/m
!   \item[rssha]
!     shaded stomatal resistance. unit: s/m
!   \item[bgap]
!     between canopy gap fraction for beam. unit: -
!   \item[wgap]
!     within canopy gap fraction for beam. unit: -
!   \item[chv]
!     sensible heat exchange coefficient over vegetated fraction. unit: s/m
!   \item[chb]
!     sensible heat exchange coefficient over bare-ground fraction. unit: s/m
!   \item[emissi]
!     surface emissivity. unit: -
!   \item[shg]
!     ground sensible heat. unit: W m-2
!   \item[shc]
!     canopy sensible heat. unit: W m-2
!   \item[shb]
!     bare ground sensible heat. unit: W m-2
!   \item[evg]
!     ground evaporation heat. unit: W m-2
!   \item[evb]
!     bare ground evaporation heat. unit: W m-2
!   \item[ghv]
!     ground heat flux. unit: W m-2
!   \item[ghb]
!     bare ground heat flux. unit: W m-2
!   \item[irg]
!     ground net long wave radiation. unit: W m-2
!   \item[irc]
!     canopy net long wave radiation. unit: W m-2
!   \item[irb]
!     bare ground net long wave radiation. unit: W m-2
!   \item[tr]
!     transpiration heat. unit: W m-2
!   \item[evc]
!     canopy evaporation heat. unit: W m-2
!   \item[chleaf]
!     leaf exchange coefficient. unit: -
!   \item[chuc]
!     under canopy exchange coefficient. unit: -
!   \item[chv2]
!     sensible heat exchange coefficient over vegetated fraction. unit: -
!   \item[chb2]
!     sensible heat exchange coefficient over bare-ground. unit: -
!   \item[fpice]
!     snow fraction in precipitation. unit: -
!   \item[sfcheadrt]
!     extra output for WRF-HYDRO. unit: m
!   \item[sigma_sm]
!     standard deviation of soil moisture. unit: %
!   \item[topt]
!     optimum transpiration air temperature.
!   \item[rgl]
!     parameter used in radiation stress function.
!   \item[rsmax]
!     maximum stomatal resistance.
!   \item[rsmin]
!     minimum Canopy Resistance. unit: s/m
!   \item[hs]
!     parameter used in vapor pressure deficit function.
!   \item[csoil]
!     vol. soil heat capacity. unit: j/m3/K
!   \item[bexp]
!     B parameter.
!   \item[dksat]
!     saturated soil hydraulic conductivity.
!   \item[dwsat]
!     saturated soil hydraulic diffusivity.
!   \item[psisat]
!     saturated soil matric potential.
!   \item[quartz]
!     soil quartz content.
!   \item[smcmax]
!     porosity, saturated value of soil moisture (volumetric).
!   \item[smcref]
!     reference soil moisture (field capacity).
!   \item[smcwlt]
!     wilting point soil moisture (volumetric).
!   \item[czil]
!     Calculate roughness length of heat. 
!   \item[frzk]
!     frozen ground parameter. 
!   \item[refdk]
!     parameters in the surface runoff parameteriz. 
!   \item[refkdt]
!     parameters in the surface runoff parameteriz. 
!   \item[slope]
!     slope index (0 - 1).
!   \item[CH2OP]
!     maximum intercepted h2o per unit lai+sai. unit: mm
!   \item[DLEAF]
!     characteristic leaf dimension. unit: m
!   \item[Z0MVT]
!     momentum roughness length. unit: m
!   \item[HVT]
!     top of canopy. unit: m
!   \item[HVB]
!     bottom of canopy. unit: m
!   \item[RC]
!     tree crown radius. unit: m
!   \item[RHOL1]
!     leaf reflectance (1=vis)
!   \item[RHOL2]
!     leaf reflectance (2=nir)
!   \item[RHOS1]
!     stem reflectance (1=vis)
!   \item[RHOS2]
!     stem reflectance (2=nir)
!   \item[TAUL1]
!     leaf transmittance (1=vis)
!   \item[TAUL2]
!     leaf transmittance (2=nir)
!   \item[TAUS1]
!     stem transmittance (1=vis)
!   \item[TAUS2]
!     stem transmittance (2=nir)
!   \item[XL]
!     leaf/stem orientation index
!   \item[CWPVT]
!     empirical canopy wind parameter
!   \item[C3PSN]
!     photosynthetic pathway (0. = c4, 1. = c3)
!   \item[KC25]
!     co2 michaelis-menten constant at 25c. unit: pa
!   \item[AKC]
!     q10 for kc25
!   \item[KO25]
!     o2 michaelis-menten constant at 25c. unit: pa
!   \item[AKO]
!     q10 for ko25
!   \item[AVCMX]
!     q10 for vcmx25
!   \item[AQE]
!     q10 for qe25
!   \item[LTOVRC]
!     leaf turnover. unit: 1/s
!   \item[DILEFC]
!     coeficient for leaf stress death. unit: 1/s
!   \item[DILEFW]
!     coeficient for leaf stress death. unit: 1/s
!   \item[RMF25]
!     leaf maintenance respiration at 25c. unit: umol co2/m**2/s
!   \item[SLA]
!     single-side leaf area per Kg. unit: m2/kg
!   \item[FRAGR]
!     fraction of growth respiration. original was 0.3
!   \item[TMIN]
!     minimum temperature for photosynthesis. unit: k
!   \item[VCMX25]
!     maximum rate of carboxylation at 25c. unit: umol co2/m**2/s
!   \item[TDLEF]
!     characteristic T for leaf freezing. unit: K
!   \item[BP]
!     minimum leaf conductance. unit: umol/m**2/s
!   \item[MP]
!     slope of conductance-to-photosynthesis relationship
!   \item[QE25]
!     quantum efficiency at 25c. unit: umol co2 / umol photon
!   \item[RMS25]
!     stem maintenance respiration at 25c. unit: umol co2/kg bio/s
!   \item[RMR25]
!     root maintenance respiration at 25c. unit: umol co2/kg bio/s
!   \item[ARM]
!     q10 for maintenance respiration
!   \item[FOLNMX]
!     foliage nitrogen concentration when f(n)=1. unit: %
!   \item[WDPOOL]
!     wood pool (switch 1 or 0) depending on woody or not. unit: -
!   \item[WRRAT]
!     wood to non-wood ratio
!   \item[MRP]
!     microbial respiration parameter. unit: umol co2 /kg c/ s
!   \item[wtrflx]
!    total water flux. unit: kg m-2 s-1
!   \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  9/4/14: Shugong Wang Initial implementation for LIS 7 and NoahMP36
!  2/1/18: Soni Yatheendradas: Added calibratable parameters for OPT
!
!EOP
  implicit none
  private
  type, public :: noahmp36dec
     !-------------------------------------------------------------------------
     ! forcing
     !-------------------------------------------------------------------------
     real               :: tair
     real               :: sfctmp
     real               :: psurf
     real               :: wind_e
     real               :: wind_n
     real               :: qair
     real               :: swdown
     real               :: lwdown
     real               :: prcp
     !-------------------------------------------------------------------------
     ! spatial parameter
     !-------------------------------------------------------------------------
     integer            :: vegetype
     integer            :: soiltype
     integer            :: slopetype
     real               :: tbot
     real               :: pblh
     !-------------------------------------------------------------------------
     ! multilevel spatial parameter
     !-------------------------------------------------------------------------
     real, pointer      :: shdfac_monthly(:)
     real, pointer      :: smceq(:)
     !-------------------------------------------------------------------------
     ! state
     !-------------------------------------------------------------------------
     real               :: albold
     real               :: sneqvo
     real, pointer      :: sstc(:)
     real, pointer      :: sh2o(:)
     real, pointer      :: smc(:)
     real               :: tah
     real               :: eah
     real               :: fwet
     real               :: canliq
     real               :: canice
     real               :: tv
     real               :: tg
     real               :: qsnow
     integer            :: isnow
     real, pointer      :: zss(:)
     real               :: snowh
     real               :: sneqv
     real, pointer      :: snowice(:)
     real, pointer      :: snowliq(:)
     real               :: zwt
     real               :: wa
     real               :: wt
     real               :: wslake
     real               :: lfmass
     real               :: rtmass
     real               :: stmass
     real               :: wood
     real               :: stblcp
     real               :: fastcp
     real               :: lai
     real               :: sai
     real               :: cm
     real               :: ch
     real               :: tauss
     real               :: smcwtd
     real               :: deeprech
     real               :: rech
     real               :: zlvl 
     real               :: albd(2)
     real               :: albi(2)
     logical            :: alb_upd_flag
     !ag (12Sep2019)
     ! 2-way coupling parameters
     real               :: rivsto
     real               :: fldsto
     real               :: fldfrc
     !-------------------------------------------------------------------------
     ! output
     !-------------------------------------------------------------------------
     real               :: fsa
     real               :: fsr
     real               :: fira
     real               :: fsh
     real               :: ssoil
     real               :: fcev
     real               :: fgev
     real               :: fctr
     real               :: ecan
     real               :: etran
     real               :: edir
     real               :: trad
     real               :: tgb
     real               :: tgv
     real               :: t2mv
     real               :: t2mb
     real               :: q2v
     real               :: q2b
     real               :: runsrf
     real               :: runsub
     real               :: apar
     real               :: psn
     real               :: sav
     real               :: sag
     real               :: fsno
     real               :: nee
     real               :: gpp
     real               :: npp
     real               :: fveg
     real               :: albedo
     real               :: qsnbot
     real               :: ponding
     real               :: ponding1
     real               :: ponding2
     real               :: rssun
     real               :: rssha
     real               :: bgap
     real               :: wgap
     real               :: chv
     real               :: chb
     real               :: emissi
     real               :: shg
     real               :: shc
     real               :: shb
     real               :: evg
     real               :: evb
     real               :: ghv
     real               :: ghb
     real               :: irg
     real               :: irc
     real               :: irb
     real               :: tr
     real               :: evc
     real               :: chleaf
     real               :: chuc
     real               :: chv2
     real               :: chb2
     real               :: fpice
     real               :: sfcheadrt
     !Added by Chandana Gangodagamage
     real :: sfchead1rt
     real :: infxs1rt
     real :: soldrain1rt

     !-------------------------------------------------------------------------
     ! read in from ARS sm data files ! SY
     !-------------------------------------------------------------------------
     !real               :: smc_std
     !-------------------------------------------------------------------------
     ! calibratable parameters for OPTUE ! SY
     !-------------------------------------------------------------------------
     ! SY: Begin from REDPRM
     !SY: begin vegetation parameters
     real               :: topt
     real               :: rgl
     real               :: rsmax
     real               :: rsmin
     real               :: hs
     real               :: nroot
     !SY: end vegetation parameters
     !SY: begin soil parameters
     real               :: csoil
     real               :: bexp
     real               :: dksat
     real               :: dwsat
     !real               :: f1 ! SY: Not used by NoahMP3.6 from REDPRM 
     real               :: psisat
     real               :: quartz
     !real               :: smcdry ! SY: Not used by NoahMP3.6 from REDPRM
     real               :: smcmax
     real               :: smcref
     real               :: smcwlt
     !SY: end soil parameters
     !SY: begin universal parameters (not dependent on SOILTYP, VEGTYP)
     real               :: czil
     real               :: frzk
     real               :: refdk
     real               :: refkdt
     real               :: slope
     !SY: end universal parameters (not dependent on SOILTYP, VEGTYP)
     ! SY: End from REDPRM
     ! SY: Begin from read_mp_veg_parameters
     real               :: CH2OP
     real               :: DLEAF
     real               :: Z0MVT
     real               :: HVT
     real               :: HVB
     real               :: RC
     real               :: RHOL1
     real               :: RHOL2
     real               :: RHOS1
     real               :: RHOS2
     real               :: TAUL1
     real               :: TAUL2
     real               :: TAUS1
     real               :: TAUS2
     real               :: XL
     real               :: CWPVT
     real               :: C3PSN
     real               :: KC25
     real               :: AKC
     real               :: KO25
     real               :: AKO
     real               :: AVCMX
     real               :: AQE
     real               :: LTOVRC
     real               :: DILEFC
     real               :: DILEFW
     real               :: RMF25
     real               :: SLA
     real               :: FRAGR
     real               :: TMIN
     real               :: VCMX25
     real               :: TDLEF
     real               :: BP
     real               :: MP
     real               :: QE25
     real               :: RMS25
     real               :: RMR25
     real               :: ARM
     real               :: FOLNMX
     real               :: WDPOOL
     real               :: WRRAT
     real               :: MRP
     ! SY: End from read_mp_veg_parameters
     !-------------------------------------------------------------------------
     ! used for constraints on calibratable parameters for OPTUE ! SY
     !-------------------------------------------------------------------------
     real               :: smcdry ! SY: Not used by NoahMP3.6 from REDPRM, but read in from table
#ifdef PARFLOW
     real, pointer      :: wtrflx(:)
#endif

  end type noahmp36dec
end module NoahMP36_module
