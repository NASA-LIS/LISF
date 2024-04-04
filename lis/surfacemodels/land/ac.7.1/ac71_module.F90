!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
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
  type, public :: ac71dec
    !-------------------------------------------------------------------------
    ! forcing
    !-------------------------------------------------------------------------
    real               :: tair
    real               :: tmax
    real               :: tmin
    real               :: tdew
    real               :: sfctmp
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
    real, pointer      :: smceq(:)
    !-------------------------------------------------------------------------
    ! state
    !-------------------------------------------------------------------------
    real, pointer      :: smc(:)
    !-------------------------------------------------------------------------
    ! AC specific
    !-------------------------------------------------------------------------
    character(len=:), allocatable :: CalendarDescription
    character(len=:), allocatable :: CalendarFile
    character(len=:), allocatable :: CalendarFileFull
    character(len=:), allocatable :: ClimateDescription
    character(len=:), allocatable :: ClimateFile
    character(len=:), allocatable :: ClimateFileFull
    character(len=:), allocatable :: Climate_Filename
    character(len=:), allocatable :: ClimDescription
    character(len=:), allocatable :: ClimFile
    character(len=:), allocatable :: CO2Description
    character(len=:), allocatable :: CO2File
    character(len=:), allocatable :: CO2FileFull
    character(len=:), allocatable :: CO2_Filename
    character(len=:), allocatable :: CropDescription
    character(len=:), allocatable :: CropFile
    character(len=:), allocatable :: CropFileFull
    character(len=:), allocatable :: Crop_Filename
    character(len=:), allocatable :: EToDescription
    character(len=:), allocatable :: EToFile
    character(len=:), allocatable :: EToFileFull
    character(len=:), allocatable :: ETo_Filename
    character(len=:), allocatable :: fEval_filename
    character(len=:), allocatable :: FullFileNameProgramParameters
    character(len=:), allocatable :: GroundwaterDescription
    character(len=:), allocatable :: GroundWaterFile
    character(len=:), allocatable :: GroundWaterFilefull
    character(len=:), allocatable :: IrriDescription
    character(len=:), allocatable :: IrriFile
    character(len=:), allocatable :: IrriFileFull
    character(len=:), allocatable :: Irrigation_Filename
    character(len=:), allocatable :: Management_Filename
    character(len=:), allocatable :: ManDescription
    character(len=:), allocatable :: ManFile
    character(len=:), allocatable :: ManFilefull
    character(len=:), allocatable :: MultipleProjectDescription
    character(len=:), allocatable :: MultipleProjectFile
    character(len=:), allocatable :: MultipleProjectFileFull
    character(len=:), allocatable :: ObservationsDescription
    character(len=:), allocatable :: ObservationsFile
    character(len=:), allocatable :: ObservationsFilefull
    character(len=:), allocatable :: OffSeasonDescription
    character(len=:), allocatable :: OffSeasonFile
    character(len=:), allocatable :: OffSeasonFilefull
    character(len=:), allocatable :: OutputName
    character(len=:), allocatable :: PathNameOutp
    character(len=:), allocatable :: PathNameProg
    character(len=:), allocatable :: PathNameSimul
    character(len=:), allocatable :: ProfDescription
    character(len=:), allocatable :: ProfFile
    character(len=:), allocatable :: ProfFilefull
    character(len=:), allocatable :: ProjectDescription
    character(len=:), allocatable :: ProjectFile
    character(len=:), allocatable :: ProjectFileFull
    character(len=:), allocatable :: RainDescription
    character(len=:), allocatable :: RainFile
    character(len=:), allocatable :: RainFileFull
    character(len=:), allocatable :: Rain_Filename
    character(len=:), allocatable :: SWCiniDescription
    character(len=:), allocatable :: SWCiniFile
    character(len=:), allocatable :: SWCiniFileFull
    character(len=:), allocatable :: TemperatureDescription
    character(len=:), allocatable :: TemperatureFile
    character(len=:), allocatable :: TemperatureFileFull
    integer        :: daynri
    integer        :: IrriInterval
    integer(int32) :: DayLastCut
    integer(int32) :: DayNr1Eval
    integer(int32) :: DayNrEval
    integer(int32) :: DaySubmerged
    integer(int32) :: GDDTadj
    integer(int32) :: InitializeRun
    integer(int32) :: IrriFirstDayNr
    integer(int32) :: MaxPlotNew
    integer(int32) :: NextSimFromDayNr
    integer(int32) :: NrCompartments
    integer(int32) :: NrCut
    integer(int32) :: NrRuns
    integer(int32) :: PreviousDayNr
    integer(int32) :: SumInterval
    integer(int32) :: Tadj
    integer(int32) :: ZiAqua
    integer        :: LineNrEval
    integer        :: PreviousStressLevel
    integer        :: StageCode
    integer        :: StressSFadjNEW
    integer       :: IniPercTAW
    integer       :: irun
    integer       :: MaxPlotTr
    integer       :: OutputAggregate
    integer       :: NoMoreCrop
    integer(intEnum) :: GenerateDepthMode
    integer(intEnum) :: GenerateTimeMode
    integer(intEnum) :: IrriMethod
    integer(intEnum) :: IrriMode
    integer(intEnum) :: TheProjectType
    logical  :: HarvestNow
    logical :: EvapoEntireSoilSurface
    logical :: GlobalIrriECw
    !logical :: NoMoreCrop
    logical :: NoYear
    logical :: Out1Wabal
    logical :: Out2Crop
    logical :: Out3Prof
    logical :: Out4Salt
    logical :: Out5CompWC
    logical :: Out6CompEC
    logical :: Out7Clim
    logical :: OutDaily
    logical :: Part1Mult
    logical :: Part2Eval
    logical :: PreDay
    logical :: StartMode
    logical :: WaterTableInProfile
    real    :: RootZoneWC_Actual
    real    :: RootZoneWC_FC
    real    :: RootZoneWC_Leaf
    real    :: RootZoneWC_SAT
    real    :: RootZoneWC_Sen
    real    :: RootZoneWC_Thresh
    real    :: RootZoneWC_WP
    real    :: RootZoneWC_ZtopAct
    real    :: RootZoneWC_ZtopFC
    real    :: RootZoneWC_ZtopThresh
    real    :: RootZoneWC_ZtopWP
    real :: alfaHI
    real :: alfaHIAdj
    real :: Bin
    real :: Bout
    real :: BprevSum
    real :: CCiActual
    real :: CCiActualWeedInfested
    real :: CCiprev
    real :: CCiTopEarlySen
    real :: CCoTotal
    real :: CCxCropWeedsNoSFstress
    real :: CCxTotal
    real :: CCxWitheredTpotNoS
    real :: CDCTotal
    real :: CGCref,GDDCGCref 
    real :: CO2i
    real :: Coeffb0
    real :: Coeffb0Salt
    real :: Coeffb1
    real :: Coeffb1Salt
    real :: Coeffb2
    real :: Coeffb2Salt
    real :: CRsalt
    real :: CRwater
    real :: DayFraction
    real :: Drain  
    real :: Eact
    real :: ECdrain 
    real :: ECiAqua
    real :: ECstorage
    real :: Epot 
    real :: FracBiomassPotSF
    real :: fWeedNoS
    real :: GDDayFraction
    real :: GDDayi
    real :: GDDCDCTotal
    real :: HItimesAT
    real :: HItimesAT1
    real :: HItimesAT2
    real :: HItimesBEF
    real :: Infiltrated 
    real :: Irrigation 
    real :: PreviousBmob
    real :: PreviousBsto
    real :: PreviousSumETo
    real :: PreviousSumGDD
    real :: RootingDepth
    real :: Runoff  
    real :: SaltInfiltr
    real :: ScorAT1
    real :: ScorAT2
    real :: StressLeaf
    real :: StressSenescence
    real :: SumETo
    real :: SumGDD
    real :: SumGDDcuts
    real :: SumGDDPrev
    real :: SumKci
    real :: SumKcTop
    real :: SumKcTopStress
    real :: Surf0
    real :: SurfaceStorage
    real :: Tact 
    real :: TactWeedInfested
    real :: TimeSenescence
    real :: Tpot 
    real :: WeedRCi
    real :: WPi
    real :: YprevSum
    real :: Zeval
    real :: Ziprev
    type(CompartmentIndividual), dimension(12) :: Compartment
    type(repCutInfoRecord) :: CutInfoRecord1
    type(repCutInfoRecord) :: CutInfoRecord2
    type(repIrriInfoRecord) :: IrriInfoRecord1
    type(repIrriInfoRecord) :: IrriInfoRecord2
    type(rep_clim) :: ClimRecord
    type(rep_clim) :: EToRecord
    type(rep_clim) :: RainRecord
    type(rep_clim) :: TemperatureRecord
    type(rep_Content) :: TotalSaltContent
    type(rep_Content) :: TotalWaterContent
    type(rep_Crop) :: Crop
    type(rep_Cuttings) :: Cuttings
    type(rep_DayEventDbl), dimension(31) :: EToDataSet
    type(rep_DayEventDbl), dimension(31) :: RainDataSet
    type(rep_DayEventDbl), dimension(31) :: TminDataSet
    type(rep_DayEventDbl), dimension(31) ::TmaxDataSet
    type(rep_DayEventInt), dimension(5) :: IrriAfterSeason
    type(rep_DayEventInt), dimension(5) :: IrriBeforeSeason
    type(rep_EffectiveRain) :: effectiverain
    type(rep_EndSeason) :: endseason
    type(rep_GwTable) :: GwTable
    type(rep_IrriECw) :: IrriECw
    type(rep_Manag) :: Management
    type(rep_Onset) :: onset
    type(rep_param) :: simulparam
    type(rep_PerennialPeriod) :: perennialperiod
    type(rep_plotPar) :: PlotVarCrop
    type(rep_RootZoneSalt) :: RootZoneSalt
    type(rep_RootZoneWC) :: RootZoneWC
    type(rep_sim) :: Simulation
    type(rep_soil) :: Soil
    type(rep_StressTot) :: StressTot
    type(rep_sum) :: PreviousSum
    type(rep_sum) :: SumWaBal
    type(rep_Transfer) :: Transfer
    type(SoilLayerIndividual), dimension(5) :: soillayer

  end type ac71dec
end module Ac71_module
