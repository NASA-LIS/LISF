!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module AC72_module
!BOP
!
! !MODULE: AC72_module
!
! !DESCRIPTION:
!  The code in this file provides a description of the
!  data structure containing the AC72 1-d variables.
!  The variables specified in the data structure include:
!
!  \begin{description}      
!   \item[landuse\_tbl\_name]
!     Noah model landuse parameter table. unit: -
!   !!!!!LB provide documentation and units
!   \end{description}
!
! !REVISION HISTORY:
!  06 MAR 2024; Louise Busschaert, initial implementation
!
!EOP
  
  use ac_global, only: &
          CompartmentIndividual,&
          SoilLayerIndividual, &
          rep_Content,&
          rep_Crop,&
          rep_CropFileSet,&
          rep_Cuttings,&
          rep_DayEventDbl,&
          rep_DayEventInt,&
          rep_EffectiveRain,&
          rep_EndSeason,&
          rep_IrriECw,&
          rep_Manag,&
          rep_Onset,&
          rep_PerennialPeriod,&
          rep_RootZoneSalt,&
          rep_RootZoneWC,&
          rep_clim,&
          rep_param,&
          rep_sim,&
          rep_soil,&
          rep_sum

  use ac_run, only: &
          repCutInfoRecord,&
          repIrriInfoRecord,&
          rep_GwTable,&
          rep_StressTot,&
          rep_Transfer,&
          rep_plotPar

  use ac_kinds, only: dp,&
                int8,&
                int32,&
                intEnum,&
                sp
  implicit none
  private
  type, public :: ac72dec
    !-------------------------------------------------------------------------
    ! forcing
    !-------------------------------------------------------------------------
    real               :: tair
    real               :: tmax
    real               :: tmin
    real               :: tdew
    real               :: psurf
    real               :: wndspd
    real               :: swdown
    real               :: prcp
    real               :: eto
    !-------------------------------------------------------------------------
    ! spatial parameter
    !-------------------------------------------------------------------------
    character(len=100) :: cropt
    integer            :: soiltype
    !-------------------------------------------------------------------------
    ! multilevel spatial parameter
    !-------------------------------------------------------------------------
    real, pointer      :: tmincli_monthly(:)
    real, pointer      :: tmaxcli_monthly(:)
    real, pointer      :: smceq(:)
    !-------------------------------------------------------------------------
    ! state
    !-------------------------------------------------------------------------
    real, pointer      :: smc(:)
    !-------------------------------------------------------------------------
    ! AC specific
    !-------------------------------------------------------------------------
     integer            :: daynri
     !real, pointer      :: ac70_soilwc(:)
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
     type(rep_RootZoneWC) :: RootZoneWC
     !type(CompartmentIndividual), dimension(GetNrCompartments()) :: Compartment
     !type(SoilLayerIndividual), dimension(GetSoil_NrSoilLayers()) :: soillayer
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
     type(rep_CropFileSet) :: CropFileSet
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
    integer :: irun
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
    logical :: Out8Irri
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
    real(dp) :: CCxWitheredTpot,CCxWitheredTpotNoS
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

    logical :: WaterTableInProfile, StartMode, NoMoreCrop, CGCadjustmentAfterCutting
    logical :: GlobalIrriECw ! for versions before 3.2 where EC of 
                             ! irrigation water was not yet recorded
    logical  :: HarvestNow
    real(sp), pointer :: Tmax_record(:)
    real(sp), pointer :: Tmin_record(:)  
    logical :: germ_reached, harv_reached
    logical :: maxR_reached, sene_reached, flowr_reached
    integer :: irri_lnr


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
     ! OUTPUT
  end type ac72dec
end module AC72_module
