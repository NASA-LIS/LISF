!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module noahmpglacier3911_module
!BOP
!
! !MODULE: noahmpglacier3911_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the noahmpglacier3911 1-d variables.
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
!   \end{description}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
!EOP
  implicit none
  private
  type, public :: noahmpgldec
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
     ! state
     !-------------------------------------------------------------------------
     real               :: zlvl 
     real               :: qsnow
     real               :: sneqvo
     integer            :: isnow
     real, pointer      :: snowice(:)
     real, pointer      :: snowliq(:)
     real, pointer      :: zss(:)
     real, pointer      :: sstc(:)
     real, pointer      :: smc(:)
     real, pointer      :: sh2o(:)
     real               :: tg
     real               :: sneqv
     real               :: snowh
     real               :: albold
     real               :: cm
     real               :: ch
     real               :: tauss
     !-------------------------------------------------------------------------
     ! spatial parameter
     !-------------------------------------------------------------------------
     real               :: tbot
     !-------------------------------------------------------------------------
     ! output
     !-------------------------------------------------------------------------
     real               :: ponding
     real               :: ponding1
     real               :: ponding2

     real               :: fsa
     real               :: fsr
     real               :: fira
     real               :: fsh
     real               :: ssoil
     real               :: fgev
     real               :: edir
     real               :: trad
     real               :: runsrf
     real               :: runsub
     real               :: sag
     real               :: albedo
     real               :: qsnbot
     real               :: emissi
     real               :: chb2
     real               :: fpice
     real               :: sfcheadrt
#if 0 
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





#endif
  end type noahmpgldec
end module noahmpglacier3911_module
