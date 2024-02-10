!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module Ac71_module
!BOP
!
! !MODULE: Ac71_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the Ac71 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[landuse\_tbl\_name]
!     Noah model landuse parameter table. unit: -
!   \item[soil\_tbl\_name]
!     Noah model soil parameter table. unit: -
!   \item[gen\_tbl\_name]
!     Noah model general parameter table. unit: -
!   \item[ac\_tbl\_name]
!     Ac parameter table. unit: -
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
!   \end{description}
!
! !REVISION HISTORY:
!  This module is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the module is defined by Sujay Kumar. 
!  9/4/14: Shugong Wang Initial implementation for LIS 7 and Ac71
!  2/1/18: Soni Yatheendradas: Added calibratable parameters for OPT
!
!EOP
  
  use ac_global, only: rep_RootZoneWC,&
                       rep_Content,&
                       rep_EffectiveRain,&
                       rep_sum,&
                       rep_RootZoneSalt,&
                       rep_sim,&
                       rep_DayEventInt,&
                       rep_Shapes,&
                       rep_soil,&
                       rep_Assimilates,&
                       rep_Onset,&
                       rep_EndSeason,&
                       rep_EffectStress,&
                       rep_IrriECw,&
                       rep_clim,&
                       rep_CropFileSet,&
                       rep_Cuttings,&
                       rep_Manag,&
                       rep_param,&
                       rep_IniSWC,&
                       rep_storage,&
                       rep_DayEventDbl,&
                       rep_Crop,&
                       rep_PerennialPeriod,&
                       CompartmentIndividual,&
                       SoilLayerIndividual,&
                       GetNrCompartments,&
                       GetSoil_NrSoilLayers

  use ac_run, only:    repIrriInfoRecord,&
                       rep_GwTable,&
                       rep_plotPar,&
                       rep_StressTot,&
                       repCutInfoRecord,&
                       rep_Transfer

  use ac_kinds, only: dp,&
                      int8,&
                      int32,&
                      intEnum,&
                      sp

  implicit none
  private
  type, public :: ac71dec
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
     real               :: PREC_ac
     real               :: TMIN_ac
     real               :: TMAX_ac
     real               :: ETo_ac
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
     real, pointer      :: smc(:)
     !-------------------------------------------------------------------------

     !!! MB: AC71
     !!!
     integer            :: daynri
     real               :: RootZoneWC_Actual
     real               :: RootZoneWC_FC
     real               :: RootZoneWC_WP
     real               :: RootZoneWC_SAT
     real               :: RootZoneWC_Leaf
     real               :: RootZoneWC_Thresh
     real               :: RootZoneWC_Sen
     real               :: RootZoneWC_ZtopAct
     real               :: RootZoneWC_ZtopFC
     real               :: RootZoneWC_ZtopWP
     real               :: RootZoneWC_ZtopThresh
     logical            :: HarvestNow
     type(rep_RootZoneWC) :: RootZoneWC
     type(rep_Content) :: TotalSaltContent
     type(rep_Content) :: TotalWaterContent
     type(rep_EffectiveRain) :: effectiverain
     type(CompartmentIndividual), dimension(12) :: Compartment
     type(SoilLayerIndividual), dimension(5) :: soillayer
     type(rep_sum) :: SumWaBal
     type(repIrriInfoRecord) :: IrriInfoRecord1, IrriInfoRecord2
     type(rep_DayEventInt), dimension(5) :: IrriBeforeSeason
     type(rep_DayEventInt), dimension(5) :: IrriAfterSeason
     type(rep_IrriECw) :: IrriECw
     type(rep_Manag) :: Management
     type(rep_PerennialPeriod) :: perennialperiod
     type(rep_param) :: simulparam
     type(rep_Cuttings) :: Cuttings
     type(rep_Onset) :: onset
     type(rep_EndSeason) :: endseason
     type(rep_Crop) :: crop
     type(rep_soil) :: Soil
     type(rep_clim)  :: TemperatureRecord, ClimRecord, RainRecord, EToRecord
     integer :: IrriInterval
     type(rep_RootZoneSalt) :: RootZoneSalt
     type(rep_sim) :: Simulation

    integer(intEnum) :: GenerateTimeMode
    integer(intEnum) :: GenerateDepthMode
    integer(intEnum) :: IrriMode
    integer(intEnum) :: IrriMethod
    integer(int32) :: DaySubmerged
    integer(int32) :: MaxPlotNew
    integer(int32) :: NrCompartments
    integer(int32) :: IrriFirstDayNr
    integer(int32) :: ZiAqua ! Depth of Groundwater table below
                             ! soil surface in centimeter 
    integer(int8) :: IniPercTAW ! Default Value for Percentage TAW for Initial
                                ! Soil Water Content Menu
    integer(int8) :: MaxPlotTr
    integer(int8) :: OutputAggregate

    integer(int32) :: NrRuns
    integer(int8) :: irun
    integer(int32) :: InitializeRun
    integer(intEnum) :: TheProjectType

    logical :: EvapoEntireSoilSurface ! True of soil wetted by RAIN (false = IRRIGATION and fw < 1)
    logical :: PreDay, OutDaily
    logical :: Out1Wabal
    logical :: Out2Crop
    logical :: Out3Prof
    logical :: Out4Salt
    logical :: Out5CompWC
    logical :: Out6CompEC
    logical :: Out7Clim
    logical :: Part1Mult,Part2Eval

    real(dp) :: CCiActual
    real(dp) :: CCiprev
    real(dp) :: CCiTopEarlySen
    real(dp) :: CRsalt ! gram/m2
    real(dp) :: CRwater ! mm/day
    real(dp) :: ECdrain ! EC drain water dS/m
    real(dp) :: ECiAqua ! EC of the groundwater table in dS/m
    real(dp) :: ECstorage !EC surface storage dS/m
    real(dp) :: Eact ! mm/day
    real(dp) :: Epot ! mm/day
    !real(dp) :: ETo_ac ! mm/day
    real(dp) :: Drain  ! mm/day
    real(dp) :: Infiltrated ! mm/day
    real(dp) :: Irrigation ! mm/day
    !real(dp) :: PREC_ac  ! mm/day
    real(dp) :: RootingDepth
    real(dp) :: Runoff  ! mm/day
    real(dp) :: SaltInfiltr ! salt infiltrated in soil profile Mg/ha
    real(dp) :: Surf0 ! surface water [mm] begin day
    real(dp) :: SurfaceStorage !mm/day
    real(dp) :: Tact ! mm/day
    real(dp) :: Tpot ! mm/day
    real(dp) :: TactWeedInfested !mm/day
    !real(dp) :: Tmax_ac ! degC
    !real(dp) :: Tmin_ac ! degC

    ! variables from run.f90
    type(rep_GwTable) :: GwTable
    type(rep_DayEventDbl), dimension(31) :: EToDataSet
    type(rep_DayEventDbl), dimension(31) :: RainDataSet
    type(rep_plotPar) :: PlotVarCrop
    type(rep_StressTot) :: StressTot
    type(repCutInfoRecord) :: CutInfoRecord1, CutInfoRecord2
    type(rep_Transfer) :: Transfer
    type(rep_DayEventDbl), dimension(31) :: TminDataSet, TmaxDataSet
    type(rep_sum) :: PreviousSum

    integer(int32) :: Tadj, GDDTadj
    integer(int32) :: DayLastCut,NrCut,SumInterval
    integer(int8)  :: PreviousStressLevel, StressSFadjNEW

    real(dp) :: Bin
    real(dp) :: Bout
    real(dp) :: GDDayi
    real(dp) :: CO2i
    real(dp) :: FracBiomassPotSF
    real(dp) :: SumETo,SumGDD, Ziprev,SumGDDPrev
    real(dp) :: CCxWitheredTpotNoS
    real(dp) :: Coeffb0,Coeffb1,Coeffb2
    real(dp) :: Coeffb0Salt,Coeffb1Salt,Coeffb2Salt
    real(dp) :: StressLeaf,StressSenescence !! stress for leaf expansion and senescence
    real(dp) :: DayFraction,GDDayFraction
    real(dp) :: CGCref,GDDCGCref 
    real(dp) :: TimeSenescence !! calendar days or GDDays
    real(dp) :: SumKcTop, SumKcTopStress, SumKci
    real(dp) :: CCoTotal, CCxTotal, CDCTotal, GDDCDCTotal, CCxCropWeedsNoSFstress
    real(dp) :: WeedRCi, CCiActualWeedInfested, fWeedNoS, Zeval
    real(dp) :: BprevSum, YprevSum, SumGDDcuts, HItimesBEF
    real(dp) :: ScorAT1, ScorAT2, HItimesAT1, HItimesAT2, HItimesAT
    real(dp) :: alfaHI, alfaHIAdj
    real(dp) :: WPi
    !! DelayedGermination
    integer(int32) :: NextSimFromDayNr !! the Simulation.FromDayNr for next run if delayed germination and KeepSWC

    !! Evaluation
    integer(int32) :: DayNr1Eval,DayNrEval
    integer(int8)  :: LineNrEval

    !! specific for StandAlone
    real(dp) :: PreviousSumETo, PreviousSumGDD, PreviousBmob,PreviousBsto
    integer(int8)  :: StageCode
    integer(int32) :: PreviousDayNr
    logical :: NoYear

    character(len=:), allocatable :: fEval_filename

    logical :: WaterTableInProfile, StartMode, NoMoreCrop
    logical :: GlobalIrriECw ! for versions before 3.2 where EC of 
                             ! irrigation water was not yet recorded


character(len=:), allocatable :: RainFile
character(len=:), allocatable :: RainFileFull
character(len=:), allocatable :: RainDescription
character(len=:), allocatable :: EToFile
character(len=:), allocatable :: EToFileFull
character(len=:), allocatable :: EToDescription
character(len=:), allocatable :: CalendarFile
character(len=:), allocatable :: CalendarFileFull
character(len=:), allocatable :: CalendarDescription
character(len=:), allocatable :: CO2File
character(len=:), allocatable :: CO2FileFull
character(len=:), allocatable :: CO2Description
character(len=:), allocatable :: IrriFile
character(len=:), allocatable :: IrriFileFull
character(len=:), allocatable :: CropFile
character(len=:), allocatable :: CropFileFull
character(len=:), allocatable :: CropDescription
character(len=:), allocatable :: PathNameProg
character(len=:), allocatable :: PathNameOutp
character(len=:), allocatable :: PathNameSimul
character(len=:), allocatable :: Climate_Filename
character(len=:), allocatable :: ETo_Filename
character(len=:), allocatable :: Rain_Filename
character(len=:), allocatable :: CO2_Filename
character(len=:), allocatable :: Crop_Filename
character(len=:), allocatable :: Management_Filename
character(len=:), allocatable :: Irrigation_Filename
character(len=:), allocatable :: ProfFile
character(len=:), allocatable :: ProfFilefull
character(len=:), allocatable :: ProfDescription
character(len=:), allocatable :: ManFile
character(len=:), allocatable :: ManFilefull
character(len=:), allocatable :: ObservationsFile
character(len=:), allocatable :: ObservationsFilefull
character(len=:), allocatable :: ObservationsDescription
character(len=:), allocatable :: OffSeasonFile
character(len=:), allocatable :: OffSeasonFilefull
character(len=:), allocatable :: OutputName
character(len=:), allocatable :: GroundWaterFile
character(len=:), allocatable :: GroundWaterFilefull
character(len=:), allocatable :: ClimateFile
character(len=:), allocatable :: ClimateFileFull
character(len=:), allocatable :: ClimateDescription
character(len=:), allocatable :: IrriDescription
character(len=:), allocatable :: ClimFile
character(len=:), allocatable :: SWCiniFile
character(len=:), allocatable :: SWCiniFileFull
character(len=:), allocatable :: SWCiniDescription
character(len=:), allocatable :: ProjectDescription
character(len=:), allocatable :: ProjectFile
character(len=:), allocatable :: ProjectFileFull
character(len=:), allocatable :: MultipleProjectDescription
character(len=:), allocatable :: MultipleProjectFile
character(len=:), allocatable :: TemperatureFile
character(len=:), allocatable :: TemperatureFileFull
character(len=:), allocatable :: TemperatureDescription
character(len=:), allocatable :: MultipleProjectFileFull
character(len=:), allocatable :: FullFileNameProgramParameters
character(len=:), allocatable :: ManDescription
character(len=:), allocatable :: ClimDescription
character(len=:), allocatable :: OffSeasonDescription
character(len=:), allocatable :: GroundwaterDescription

  end type ac71dec
end module Ac71_module
