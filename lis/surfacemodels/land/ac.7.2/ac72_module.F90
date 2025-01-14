!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
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
! !REVISION HISTORY:
!  04 NOV 2024; Louise Busschaert, initial implementation
!
!EOP
  
  use ac_global, only: &
        CompartmentIndividual,&
        rep_clim,&
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
        rep_param,&
        rep_PerennialPeriod,&
        rep_RootZoneSalt,&
        rep_RootZoneWC,&
        rep_sim,&
        rep_soil,&
        rep_sum, &
        SoilLayerIndividual

  use ac_run, only: &
          repCutInfoRecord,&
          repIrriInfoRecord,&
          rep_GwTable,&
          rep_StressTot,&
          rep_Transfer,&
          rep_plotPar

  use ac_kinds, only: int8,&
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
        ! AC specific (all AC global variables)
        !-------------------------------------------------------------------------
        real               :: cycle_complete
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
        integer(int32) :: ZiAqua
        integer(int8) :: IniPercTAW
        integer(int8) :: MaxPlotTr
        integer(int8) :: OutputAggregate

        integer(int32) :: NrRuns
        integer :: irun
        integer(int32) :: InitializeRun
        integer(int32) :: read_Trecord
        integer(intEnum) :: TheProjectType

        logical :: EvapoEntireSoilSurface
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

        real :: CCiActual
        real :: CCiprev
        real :: CCiTopEarlySen
        real :: CRsalt
        real :: CRwater
        real :: ECdrain
        real :: ECiAqua
        real :: ECstorage
        real :: Eact
        real :: Epot
        real :: Drain
        real :: Infiltrated
        real :: Irrigation
        real :: RootingDepth
        real :: Runoff
        real :: SaltInfiltr
        real :: Surf0
        real :: SurfaceStorage
        real :: Tact
        real :: Tpot
        real :: TactWeedInfested

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
        integer(int32)  :: PreviousStressLevel, StressSFadjNEW

        real :: Bin
        real :: Bout
        real :: GDDayi
        real :: CO2i
        real :: FracBiomassPotSF
        real :: SumETo,SumGDD, Ziprev,SumGDDPrev
        real :: CCxWitheredTpot,CCxWitheredTpotNoS
        real :: Coeffb0,Coeffb1,Coeffb2
        real :: Coeffb0Salt,Coeffb1Salt,Coeffb2Salt
        real :: StressLeaf,StressSenescence
        real :: DayFraction,GDDayFraction
        real :: CGCref,GDDCGCref 
        real :: TimeSenescence
        real :: SumKcTop, SumKcTopStress, SumKci
        real :: CCoTotal, CCxTotal, CDCTotal, GDDCDCTotal, CCxCropWeedsNoSFstress
        real :: WeedRCi, CCiActualWeedInfested, fWeedNoS, Zeval
        real :: BprevSum, YprevSum, SumGDDcuts, HItimesBEF
        real :: ScorAT1, ScorAT2, HItimesAT1, HItimesAT2, HItimesAT
        real :: alfaHI, alfaHIAdj
        real :: WPi
        integer(int32) :: NextSimFromDayNr

        !! Evaluation
        integer(int32) :: DayNr1Eval,DayNrEval
        integer(int8)  :: LineNrEval

        !! specific for StandAlone
        real :: PreviousSumETo, PreviousSumGDD, PreviousBmob,PreviousBsto
        integer(int8)  :: StageCode
        integer(int32) :: PreviousDayNr
        logical :: NoYear

        character(len=:), allocatable :: fEval_filename

        logical :: WaterTableInProfile, StartMode, CGCadjustmentAfterCutting
        logical :: GlobalIrriECw
        logical  :: HarvestNow
        real(sp), pointer :: Tmax_record(:)
        real(sp), pointer :: Tmin_record(:)  
        integer :: irri_lnr

        logical :: NoMoreCrop

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
  end type ac72dec
end module AC72_module
