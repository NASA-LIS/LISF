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
    integer(int32) :: DayLastCut,NrCut,SumInterval
    integer(int32) :: DayNr1Eval,DayNrEval
    integer(int32) :: DaySubmerged
    integer(int32) :: InitializeRun
    integer(int32) :: IrriFirstDayNr
    integer(int32) :: MaxPlotNew
    integer(int32) :: NextSimFromDayNr
    integer(int32) :: NrCompartments
    integer(int32) :: NrRuns
    integer(int32) :: PreviousDayNr
    integer(int32) :: Tadj, GDDTadj
    integer(int32) :: ZiAqua
    integer(int8)  :: LineNrEval
    integer(int8)  :: PreviousStressLevel, StressSFadjNEW
    integer(int8)  :: StageCode
    integer(int8) :: IniPercTAW
    integer(int8) :: irun
    integer(int8) :: MaxPlotTr
    integer(int8) :: OutputAggregate
    integer(intEnum) :: GenerateDepthMode
    integer(intEnum) :: GenerateTimeMode
    integer(intEnum) :: IrriMethod
    integer(intEnum) :: IrriMode
    integer(intEnum) :: TheProjectType
    logical  :: HarvestNow
    logical :: EvapoEntireSoilSurface
    logical :: GlobalIrriECw
    logical :: NoYear
    logical :: Out1Wabal
    logical :: Out2Crop
    logical :: Out3Prof
    logical :: Out4Salt
    logical :: Out5CompWC
    logical :: Out6CompEC
    logical :: Out7Clim
    logical :: Part1Mult,Part2Eval
    logical :: PreDay, OutDaily
    logical :: WaterTableInProfile, StartMode, NoMoreCrop
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
    real(dp) :: alfaHI, alfaHIAdj
    real(dp) :: Bin
    real(dp) :: Bout
    real(dp) :: BprevSum, YprevSum, SumGDDcuts, HItimesBEF
    real(dp) :: CCiActual
    real(dp) :: CCiprev
    real(dp) :: CCiTopEarlySen
    real(dp) :: CCoTotal, CCxTotal, CDCTotal, GDDCDCTotal, CCxCropWeedsNoSFstress
    real(dp) :: CCxWitheredTpotNoS
    real(dp) :: CGCref,GDDCGCref 
    real(dp) :: CO2i
    real(dp) :: Coeffb0,Coeffb1,Coeffb2
    real(dp) :: Coeffb0Salt,Coeffb1Salt,Coeffb2Salt
    real(dp) :: CRsalt
    real(dp) :: CRwater
    real(dp) :: DayFraction,GDDayFraction
    real(dp) :: Drain  
    real(dp) :: Eact
    real(dp) :: ECdrain 
    real(dp) :: ECiAqua
    real(dp) :: ECstorage
    real(dp) :: Epot 
    real(dp) :: FracBiomassPotSF
    real(dp) :: GDDayi
    real(dp) :: Infiltrated 
    real(dp) :: Irrigation 
    real(dp) :: PreviousSumETo, PreviousSumGDD, PreviousBmob,PreviousBsto
    real(dp) :: RootingDepth
    real(dp) :: Runoff  
    real(dp) :: SaltInfiltr
    real(dp) :: ScorAT1, ScorAT2, HItimesAT1, HItimesAT2, HItimesAT
    real(dp) :: StressLeaf,StressSenescence
    real(dp) :: SumETo,SumGDD, Ziprev,SumGDDPrev
    real(dp) :: SumKcTop, SumKcTopStress, SumKci
    real(dp) :: Surf0
    real(dp) :: SurfaceStorage
    real(dp) :: Tact 
    real(dp) :: TactWeedInfested
    real(dp) :: TimeSenescence
    real(dp) :: Tpot 
    real(dp) :: WeedRCi, CCiActualWeedInfested, fWeedNoS, Zeval
    real(dp) :: WPi
    type(CompartmentIndividual), dimension(12) :: Compartment
    type(repCutInfoRecord) :: CutInfoRecord1, CutInfoRecord2
    type(repIrriInfoRecord) :: IrriInfoRecord1, IrriInfoRecord2
    type(rep_clim)  :: TemperatureRecord, ClimRecord, RainRecord, EToRecord
    type(rep_Content) :: TotalSaltContent
    type(rep_Content) :: TotalWaterContent
    type(rep_Crop) :: crop
    type(rep_Cuttings) :: Cuttings
    type(rep_DayEventDbl), dimension(31) :: EToDataSet
    type(rep_DayEventDbl), dimension(31) :: RainDataSet
    type(rep_DayEventDbl), dimension(31) :: TminDataSet, TmaxDataSet
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
