module ac_run

use ac_climprocessing, only:    GetDecadeEToDataset, &
                                GetDecadeRainDataSet, &
                                GetMonthlyEToDataset, &
                                GetMonthlyRainDataset
use ac_global, only:    AdjustSizeCompartments, &
                        CompartmentIndividual, &
                        datatype_daily, &
                        datatype_decadely, &
                        datatype_monthly, &
                        DaysInMonth, &
                        DegreesDay, &
                        DetermineDate, &
                        DetermineDayNr, &
                        DetermineRootZoneWC, &
                        DetermineSaltContent, &
                        ECeComp, &
                        Equiv, &
                        FileExists, &
                        GetCCiActual, &
                        GetClimRecord_FromY, &
                        GetCompartment_i, &
                        GetCompartment_i, &
                        GetCompartment_Layer, &
                        GetCompartment_Theta, &
                        GetCompartment_Thickness, &
                        GetCrop_CCEffectEvapLate, &
                        GetCrop_CCo, &
                        GetPart2Eval, &
                        getcrop_ccsaltdistortion, &
                        GetCrop_CCx, &
                        GetCrop_CDC, &
                        GetCrop_CGC, &
                        GetCrop_CGC, &
                        GetCrop_Day1, &
                        GetCrop_DayN, &
                        GetCrop_DaysToCCini, &
                        GetCrop_DaysToFlowering, &
                        GetCrop_DaysToFullCanopy, &
                        GetCrop_DaysToGermination, &
                        GetCrop_DaysToHarvest, &
                        GetCrop_DaysToSenescence, &
                        GetCrop_DeterminancyLinked, &
                        GetCrop_dHIdt, &
                        GetCrop_DryMatter, &
                        GetCrop_GDDaysToCCini, &
                        GetCrop_GDDaysToFlowering, &
                        GetCrop_GDDaysToFullCanopy, &
                        GetCrop_GDDaysToGermination, &
                        GetCrop_GDDaysToHarvest, &
                        GetCrop_GDDaysToSenescence, &
                        GetCrop_GDDCDC, &
                        GetCrop_GDDCGC, &
                        GetCrop_GDDLengthFlowering, &
                        GetCrop_GDtranspLow, &
                        GetCrop_GDtranspLow, &
                        GetCrop_HI, &
                        GetCrop_KcDecline, &
                        GetCrop_KcTop, &
                        GetCrop_Length_i, &
                        GetCrop_LengthFlowering, &
                        GetCrop_ModeCycle, &
                        GetCrop_StressResponse, &
                        GetCrop_StressResponse_Calibrated, &
                        GetCrop_subkind, &
                        GetCrop_ModeCycle, &
                        GetCrop_RootMax, &
                        GetCrop_Tbase, &
                        GetCrop_Tupper, &
                        GetCrop_Length_i, &
                        GetCrop_WP, &
                        GetCrop_WPy, &
                        GetCRSalt, &
                        GetCRWater, &
                        GetDrain, &
                        GetEact, &
                        GetECdrain, &
                        GetECiAqua, &
                        GetEpot, &
                        GetETo, &
                        GetEToFile, &
                        GetEToFilefull, &
                        GetEToRecord_DataType, &
                        GetEToRecord_FromDayNr, &
                        GetManagement_FertilityStress, &
                        GetGroundWaterFile, &
                        GetGroundWaterFileFull, &
                        GetInfiltrated, &
                        GetIrriFile, &
                        GetIrriFilefull, &
                        GetIrriFirstDayNr, &
                        GetIrrigation, &
                        GetIrriMode, &
                        GetManagement_Cuttings_Criterion, &
                        GetManagement_Cuttings_Day1, &
                        GetManagement_Cuttings_FirstDayNr, &
                        GetManagement_Cuttings_Generate, &
                        GetManagement_Cuttings_HarvestEnd, &
                        GetManagement_Cuttings_NrDays, &
                        GetManagement_FertilityStress, &
                        GetNrCompartments, &
                        GetObservationsFile, &
                        GetObservationsFilefull, &
                        GetOutputAggregate, &
                        GetOutputName, &
                        getpart1mult, &
                        GetPathNameOutp, &
                        GetPathNameProg, &
                        GetPathNameSimul, &
                        GetSimulation_DelayedDays, &
                        GetRain, &
                        GetRainFile, &
                        GetRainFilefull, &
                        GetRainRecord_DataType, &
                        GetRainRecord_FromDayNr, &
                        GetRootingDepth, &
                        GetRootZoneSalt_ECe, &
                        GetRootZoneSalt_ECsw, &
                        GetRootZoneSalt_KsSalt, &
                        GetRootZoneWC_actual, &
                        GetRootZoneWC_FC, &
                        GetRootZoneWC_Leaf, &
                        GetRootZoneWC_SAT, &
                        GetRootZoneWC_Sen, &
                        GetRootZoneWC_Thresh, &
                        GetRootZoneWC_WP, &
                        GetRunoff, &
                        GetSaltInfiltr, &
                        GetSimulation_FromDayNr, &
                        GetSimulation_IrriECw, &
                        GetSimulation_SalinityConsidered, &
                        GetSimulation_Storage_Btotal, &
                        GetSimulation_SumGDD, &
                        GetSimulation_SWCtopsoilConsidered, &
                        GetSimulation_ToDayNr, &
                        GetSimulParam_GDDMethod, &
                        GetSimulParam_Tmax, &
                        GetSimulParam_Tmin, &
                        GetSoil_RootMax, &
                        GetSoilLayer_SAT, &
                        GetSumWaBal_Biomass, &
                        GetSumWaBal_BiomassPot, &
                        GetSumWaBal_BiomassUnlim, &
                        GetSumWaBal_SaltIn, &
                        GetSumWaBal_SaltOut, &
                        GetSumWaBal_CRsalt, &
                        GetSumWaBal_CRwater, &
                        GetSumWaBal_Drain, &
                        GetSumWaBal_ECropCycle, &
                        GetSumWaBal_Infiltrated, &
                        GetSumWaBal_Irrigation, &
                        GetSumWaBal_Rain, &
                        GetSumWaBal_Runoff, &
                        GetSumWaBal_Eact, &
                        GetSUmWaBal_Epot, &
                        GetSUmWaBal_Tact, &
                        GetSumWaBal_Tpot, &
                        GetSumWaBal_TrW, &
                        GetSumWaBal_YieldPart, &
                        GetSurfaceStorage, &
                        GetTact, &
                        GetTactWeedInfested, &
                        GetTemperatureFile, &
                        GetTemperatureFilefull, &
                        GetTemperatureRecord_DataType, &
                        GetTemperatureRecord_FromDayNr, &
                        GetTmax, &
                        GetTmin, &
                        GetTotalSaltContent_endDay, &
                        GetTotalWaterContent_EndDay, &
                        GetTpot, &
                        GetZiAqua, &
                        IrriMode_Generate, &
                        IrriMode_Manual, &
                        KsAny, &
                        KsTemperature, &
                        rep_DayEventDbl, &
                        LeapYear, &
                        GetOut1Wabal, &
                        GetOut2Crop, &
                        GetOut3Prof, &
                        GetOut4Salt, &
                        GetOut5CompWC, &
                        GetOut6CompEC, &
                        GetOut7Clim, &
                        rep_DayEventDbl, &
                        rep_sum, &
                        SetCompartment_i, &
                        SetCompartment_Theta, &
                        SetETo, &
                        SetRain, &
                        SetRootZoneSalt_ECe, &
                        SetRootZoneSalt_ECsw, &
                        SetRootZoneSalt_KsSalt, &
                        SetRootZoneWC_actual, &
                        SetRootZoneWC_FC, &
                        SetRootZoneWC_Leaf, &
                        SetRootZoneWC_SAT, &
                        SetRootZoneWC_Sen, &
                        SetRootZoneWC_Thresh, &
                        SetRootZoneWC_WP, &
                        SetSimulation_IrriECw, &
                        SetSimulation_SumGDD, &
                        SetSimulation_SWCtopSoilConsidered, &
                        GetSimulation_DelayedDays, &
                        GetSimulation_MultipleRun, &
                        GetSimulation_NrRuns, &
                        SetTmax, &
                        SetTmin, &
                        SplitStringInThreeParams, &
                        SplitStringInTwoParams, &
                        subkind_Grain, &
                        subkind_Tuber, &
                        undef_int, &
                        rep_EffectStress, &
                        EvapZmin, &
                        IrriMode_Inet, &
                        modeCycle_GDDays, &
                        GetOutDaily, &
                        plant_Seed, &
                        subkind_Forage, &
                        GetCompartment, &
                        GetSimulParam_ConstGwt, &
                        GetSWCIniFile, &
                        GetCrop_pdef, &
                        SetCrop_pActStom, &
                        SetCrop_pSenAct, &
                        GetCrop_pSenescence, &
                        SetCrop_pLeafAct, &
                        GetCrop_pLeafDefUL, &
                        CO2ForSimulationPeriod, &
                        GetCrop_ECemax, &
                        GetCrop_ECemin, &
                        GetSimulation_EffectStress, &
                        GetSimulation_EffectStress_RedCGC, &
                        GetSimulation_EffectStress_RedCCX, &
                        GetCrop_DaysToFullCanopySF, &
                        GetManagement_FertilityStress, &
                        GetCrop_DaysToFullCanopySF, &
                        SeasonalSumOfKcPot, &
                        SetSimulation_RCadj, &
                        GetManagement_WeedRC, &
                        MultiplierCCxSelfThinning, &
                        GetSimulation_YearSeason, &
                        GetManagement_WeedRC, &
                        GetManagement_WeedAdj, &
                        GetSimulation_RCadj, &
                        GetSimulation_SumGDDfromDay1, &
                        GetCrop_Planting, &
                        GetRootZoneSalt_KsSalt, &
                        GetRootingDepth, &
                        GetRootZoneSalt_KsSalt, &
                        GetRootZoneSalt_ECswFC, &
                        GetRootZoneSalt_ECsw, &
                        GetRootZoneSalt_ECe, &
                        GetManagement_BundHeight, &
                        GetManagement_Cuttings_Considered, &
                        GetCrop_DaysToMaxRooting, &
                        GetSoil_RootMax, &
                        GetCrop_RootMin, &
                        GetSimulation_Zrini, &
                        GetCCiPrev, &
                        GetCrop_GDDaysToMaxRooting, &
                        GetSoil_RootMax, &
                        GetCropFileFull, &
                        GetCrop_Assimilates_Mobilized, &
                        GetSimulation_Storage_Season, &
                        GetSimulation_Bini, &
                        GetSimulation_CCini, &
                        GetCrop_PlantingDens, &
                        GetCrop_GDDaysToFullCanopySF, &
                        GetSimulation_ResetIniSWC, &
                        GetManagement_WeedShape, &
                        GetCrop_YearCCx, &
                        GetCrop_CCxRoot, &
                        CCmultiplierWeed, &
                        CCmultiplierWeedAdjusted, &
                        GetSimulation_EffectStress_CDecline, &
                        GetCrop_SizeSeedling, &
                        GetSimulation_Storage_CropString, &
                        GetSimulation_Storage_Btotal, &
                        GetCrop_RootMax, &
                        GetCrop_RootShape, &
                        CCiNoWaterStressSF, &
                        MultiplierCCoSelfThinning, &
                        ActualRootingDepth, &
                        setcrop_gddaystofullcanopysf, &
                        setsimulation_sumgddfromday1, &
                        determinerootzonesaltcontent, &
                        setsimulation_storage_btotal, &
                        setcciprev, &
                        setsimulation_delayeddays, &
                        setcrop_ccoadjusted, &
                        setsumwabal_biomasstot, &
                        setsumwabal_biomasspot, &
                        checkforwatertableinprofile, &
                        setsimulation_dayanaero, &
                        setecstorage, &
                        setsumwabal_biomassunlim, &
                        setsimulation_hifinal, &
                        setcrop_ccxwithered, &
                        timetomaxcanopysf, &
                        setsimulation_fromdaynr, &
                        setsimulation_salinityconsidered, &
                        setsurfacestorage, &
                        setevapoentiresoilsurface, &
                        setsimulation_evapz, &
                        setmanagement_fertilitystress, &
                        setrootzonesalt_ecsw, &
                        setrootzonesalt_ece, &
                        setsimulation_evaplimiton, &
                        setrootingdepth, &
                        setcrop_ccxadjusted, &
                        setsimulation_effectstress_redcgc, &
                        setsimulation_evapwcsurf, &
                        setcrop_daystofullcanopysf, &
                        cropstressparameterssoilfertility, &
                        setsimulation_protectedseedling, &
                        setsumwabal_biomass, &
                        setcciactual, &
                        setsimulation_effectstress, &
                        setrootzonesalt_kssalt, &
                        setsimulation_storage_cropstring, &
                        setpreday, &
                        setsimulation_storage_season, &
                        setrootzonesalt_ecswfc, &
                        setsimulation_effectstress_redccx, &
                        setsimulation_sumetostress, &
                        setsimulation_germinate, &
                        TimeCuttings_DryB, &
                        TimeCuttings_DryY, &
                        TimeCuttings_FreshY, &
                        TimeCuttings_IntDay, &
                        TimeCuttings_IntGDD, &
                        typeproject_typenone, &
                        typeproject_typepro, &
                        typeproject_typeprm, &
                        undef_double, &
                        undef_int, &
                        GetManFile, &
                        GetManFileFull, &
                        GetCompartment_Theta, &
                        GetClimRecord_FromY, &
                        GetCCiActual, &
                        max_No_compartments, &
                        GetSimulation_MultipleRunWithKeepSWC, &
                        GetSimulation_MultipleRunConstZrx, &
                        GetSimulation_IniSWC_AtFC, &
                        GetMultipleProjectFileFull, &
                        GetProjectFileFull, &
                        GetSUmWaBal, &
                        GetPart1Mult, &
                        GetPart2Eval, &
                        GetObservationsFile, &
                        CalculateAdjustedFC, &
                        ResetSWCToFC, &
                        SetSumWaBal, &
                        GlobalZero, &
                        SetCompartment, &
                        calculateadjustedfc, &
                        setcompartment, &
                        Rep_DayEventInt, &
                        GetIrriBeforeSeason_i, &
                        GetIrriAfterSeason_i, &
                        GetGenerateTimeMode, &
                        GenerateTimeMode_WaterBetweenBunds, &
                        GenerateDepthMode_FixDepth, &
                        GenerateTimeMode_FixInt, &
                        GenerateTimeMode_AllRAW, &
                        GenerateTimeMode_AllDepl, &
                        GetGenerateDepthMode, &
                        GetSurfaceStorage, &
                        SetIrrigation, &
                        GetSoilLayer_WP, &
                        GetSoilLayer_FC, &
                        GetSimulParam_PercRAW, &
                        GetCompartment_Theta, &
                        GetCrop_Assimilates_Period, &
                        GetCrop_Assimilates_Stored, &
                        GetSimulation_SumEToStress, &
                        CanopyCoverNoStressSF, &
                        GetCCiActual, &
                        GetManagement_Cuttings_CCcut, &
                        GetPreDay, &
                        GetSimulation_SWCtopSoilConsidered, &
                        CCiniTotalFromTimeToCCini, &
                        modeCycle_CalendarDays, &
                        GetCrop_AnaeroPoint, &
                        KsTemperature, &
                        GetSumWaBal_BiomassTot, &
                        setsumwabal_irrigation, &
                        setsumwabal_yieldpart, &
                        seteciaqua, &
                        setsimulation_swctopsoilconsidered, &
                        settactweedinfested, &
                        setziaqua, &
                        determinerootzonewc
use ac_inforesults, only:       WriteAssessmentSimulation
use ac_kinds, only: dp, &
                    int8, &
                    int32, &
                    intenum
use ac_project_input, only: ProjectInput
use ac_rootunit, only:  AdjustedRootingDepth
use ac_simul,    only:  Budget_module, &
                        determinepotentialbiomass, &
                        determinebiomassandyield
use ac_tempprocessing, only:    CCxSaltStressRelationship, &
                                GetDecadeTemperatureDataSet, &
                                GetMonthlyTemperaturedataset, &
                                GrowingDegreeDays, &
                                LoadSimulationRunProject, &
                                StressBiomassRelationship, &
                                temperaturefilecoveringcropperiod
use ac_utils, only: assert, &
                    GetAquaCropDescriptionWithTimeStamp, &
                    roundc
use iso_fortran_env, only: iostat_end
implicit none


type rep_GwTable
    integer(int32) :: DNr1, DNr2
        !! Undocumented
    integer(int32) :: Z1, Z2
        !! cm
    real(dp) :: EC1, EC2
        !! dS/m
end type rep_GwTable


type rep_plotPar
    real(dp) :: PotVal, ActVal
        !! Undocumented
end type rep_plotPar


type repIrriInfoRecord
    logical :: NoMoreInfo
        !! Undocumented
    integer(int32) :: FromDay
        !! Undocumented
    integer(int32) :: ToDay
        !! Undocumented
    integer(int32) :: TimeInfo
        !! Undocumented
    integer(int32) :: DepthInfo
        !! Undocumented
end type repIrriInfoRecord


type rep_StressTot
    real(dp) :: Salt
        !! Undocumented
    real(dp) :: Temp
        !! Undocumented
    real(dp) :: Exp
        !! Undocumented
    real(dp) :: Sto
        !! Undocumented
    real(dp) :: Weed
        !! Undocumented
    integer(int32) :: NrD
        !! Undocumented
end type rep_StressTot


type repCutInfoRecord
    logical :: NoMoreInfo
        !! Undocumented
    integer(int32) :: FromDay
        !! Undocumented
    integer(int32) :: ToDay
        !! Undocumented
    integer(int32) :: IntervalInfo
        !! Undocumented
    real(dp) :: IntervalGDD
        !! Undocumented
    real(dp) :: MassInfo
        !! Undocumented
end type repCutInfoRecord


type rep_Transfer
    logical :: Store
        !! transfer of assimilates from above ground parts to root system is active
    logical :: Mobilize
        !! transfer of assimialtes from root system to above ground parts is active
    real(dp) :: ToMobilize
        !! Total mass of assimilates (ton/ha) to mobilize at start of the season
    real(dp) :: Bmobilized
        !! Cumulative sum of assimilates (ton/ha) mobilized form root system
end type rep_Transfer

character(len=:), allocatable :: TheProjectFile

integer :: fDaily  ! file handle
integer :: fDaily_iostat  ! IO status
integer :: fRun  ! file handle
integer :: fRun_iostat  ! IO status
integer :: fIrri  ! file handle
integer :: fIrri_iostat  ! IO status
integer :: fEToSIM ! file handle
integer :: fEToSIM_iostat ! IO status
integer :: fEval ! file handle
integer :: fEval_iostat ! IO status
integer :: fRainSIM ! file handle
integer :: fRainSIM_iostat ! IO status
integer :: fTempSIM ! file handle
integer :: fTempSIM_iostat ! IO status
integer :: fCuts ! file handle
integer :: fCuts_iostat ! IO status
integer :: fObs ! file handle
integer :: fObs_iostat ! IO status
integer :: fHarvest  ! file handle
integer :: fHarvest_iostat  ! IO status
character(len=:), allocatable :: fHarvest_filename  ! file name

type(rep_GwTable) :: GwTable
type(rep_DayEventDbl), dimension(31) :: EToDataSet
type(rep_DayEventDbl), dimension(31) :: RainDataSet
type(rep_plotPar) :: PlotVarCrop
type(repIrriInfoRecord) :: IrriInfoRecord1, IrriInfoRecord2
type(rep_StressTot) :: StressTot
type(repCutInfoRecord) :: CutInfoRecord1, CutInfoRecord2
type(rep_Transfer) :: Transfer
type(rep_DayEventDbl), dimension(31) :: TminDataSet, TmaxDataSet
type(rep_sum) :: PreviousSum

integer(int32) :: DayNri
integer(int32) :: IrriInterval
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



contains


subroutine open_file(fhandle, filename, mode, iostat)
    !! Opens a file in the given mode.
    integer, intent(out) :: fhandle
        !! file handle to be used for the open file
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')
    integer, intent(out) :: iostat
        !! IO status returned by open()

    logical :: file_exists

    inquire(file=filename, exist=file_exists)

    if (mode == 'r') then
        open(newunit=fhandle, file=trim(filename), status='old', &
             action='read', iostat=iostat)
    elseif (mode == 'a') then
        if (file_exists) then
            open(newunit=fhandle, file=trim(filename), status='old', &
                 position='append', action='write', iostat=iostat)
        else
            open(newunit=fhandle, file=trim(filename), status='new', &
                 action='write', iostat=iostat)
        end if
    elseif (mode == 'w') then
        open(newunit=fhandle, file=trim(filename), status='replace', &
             action='write', iostat=iostat)
    end if
end subroutine open_file


subroutine write_file(fhandle, line, advance, iostat)
    !! Writes one line to a file.
    integer, intent(in) :: fhandle
        !! file handle of an already-opened file
    character(len=*), intent(in) :: line
        !! line to write to the file
    logical, intent(in) :: advance
        !! whether or not to append a newline character
    integer, intent(out) :: iostat
        !! IO status returned by write()

    character(len=:), allocatable :: advance_str

    if (advance) then
        advance_str = 'yes'
    else
        advance_str = 'no'
    end if

    write(fhandle, '(a)', advance=advance_str, iostat=iostat) line
end subroutine write_file


function read_file(fhandle, iostat) result(line)
    !! Returns the next line read from the given file.
    integer, intent(in) :: fhandle
        !! file handle of an already-opened file
    integer, intent(out) :: iostat
        !! IO status returned by read()
    character(len=:), allocatable :: line
        !! string which will contain the content of the next line

    integer, parameter :: length = 1024  ! max. no of characters

    allocate(character(len=length) :: line)
    read(fhandle, '(a)', iostat=iostat) line
    line = trim(line)
end function read_file


!! Section for Getters and Setters for global variables

! TheProjectFile

function GetTheProjectFile() result(str)
    !! Getter for the "TheProjectFile" global variable.
    character(len=:), allocatable :: str

    str = TheProjectFile
end function GetTheProjectFile


subroutine SetTheProjectFile(str)
    !! Setter for the "TheProjectFile" global variable.
    character(len=*), intent(in) :: str

    TheProjectFile = str
end subroutine SetTheProjectFile

! fDaily


subroutine fDaily_open(filename, mode)
    !! Opens the given file, assigning it to the 'fDaily' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')

    call open_file(fDaily, filename, mode, fDaily_iostat)
end subroutine fDaily_open


subroutine fDaily_write(line, advance_in)
    !! Writes the given line to the fDaily file.
    character(len=*), intent(in) :: line
        !! line to write
    logical, intent(in), optional :: advance_in
        !! whether or not to append a newline character

    logical :: advance

    if (present(advance_in)) then
        advance = advance_in
    else
        advance = .true.
    end if
    call write_file(fDaily, line, advance, fDaily_iostat)
end subroutine fDaily_write


subroutine fDaily_close()
    close(fDaily)
end subroutine fDaily_close


! fRun

subroutine fRun_open(filename, mode)
    !! Opens the given file, assigning it to the 'fRun' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')

    call open_file(fRun, filename, mode, fRun_iostat)
end subroutine fRun_open


subroutine fRun_write(line, advance_in)
    !! Writes the given line to the fRun file.
    character(len=*), intent(in) :: line
        !! line to write
    logical, intent(in), optional :: advance_in
        !! whether or not to append a newline character

    logical :: advance

    if (present(advance_in)) then
        advance = advance_in
    else
        advance = .true.
    end if
    call write_file(fRun, line, advance, fRun_iostat)
end subroutine fRun_write


subroutine fRun_close()
    close(fRun)
end subroutine fRun_close


! fEval

subroutine fEval_open(filename, mode)
    !! Opens the given file, assigning it to the 'fEval' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')

    call open_file(fEval, filename, mode, fEval_iostat)
end subroutine fEval_open


subroutine fEval_write(line, advance_in)
    !! Writes the given line to the fEval file.
    character(len=*), intent(in) :: line
        !! line to write
    logical, intent(in), optional :: advance_in
        !! whether or not to append a newline character

    logical :: advance

    if (present(advance_in)) then
        advance = advance_in
    else
        advance = .true.
    end if
    call write_file(fEval, line, advance, fEval_iostat)
end subroutine fEval_write


subroutine fEval_close()
    close(fEval)
end subroutine fEval_close


subroutine fEval_erase()
    call unlink(GetfEval_filename())
end subroutine fEval_erase


function GetfEval_filename() result(filename)
    !! Getter for the fEval_filename

    character(len=:), allocatable :: filename

    filename = fEval_filename
end function GetfEval_filename


subroutine SetfEval_filename(filename)
    !! Setter for the fEval_filename
    character(len=*), intent(in) :: filename

    fEval_filename = filename
end subroutine SetfEval_filename


! fIrri

subroutine fIrri_open(filename, mode)
    !! Opens the given file, assigning it to the 'fIrri' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')

    call open_file(fIrri, filename, mode, fIrri_iostat)
end subroutine fIrri_open


function fIrri_read() result(line)
    !! Returns the next line read from the 'fIrri' file.
    character(len=:), allocatable :: line
        !! name of the file to assign the file handle to

    line = read_file(fIrri, fIrri_iostat)
end function fIrri_read


function fIrri_eof() result(eof)
    !! Returns whether the end of the 'fIrri' file has been reached.
    logical :: eof

    eof = fIrri_iostat == iostat_end
end function fIrri_eof


subroutine fIrri_close()
    close(fIrri)
end subroutine fIrri_close


! fEToSIM

subroutine fEToSIM_open(filename, mode)
    !! Opens the given file, assigning it to the 'fEToSIM' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')
    call open_file(fEToSIM, filename, mode, fEToSIM_iostat)
end subroutine fEToSIM_open


function fEToSIM_read() result(line)
    !! Returns the next line read from the 'fEToSIM' file.
    character(len=:), allocatable :: line
        !! name of the file to assign the file handle to

    line = read_file(fEToSIM, fEToSIM_iostat)
end function fEToSIM_read


subroutine fEToSIM_close()
    close(fEToSIM)
end subroutine fEToSIM_close


! fTempSIM

subroutine fTempSIM_open(filename, mode)
    !! Opens the given file, assigning it to the 'fTempSIM' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')
    call open_file(fTempSIM, filename, mode, fTempSIM_iostat)
end subroutine fTempSIM_open


function fTempSIM_read() result(line)
    !! Returns the next line read from the 'fTempSIM' file.
    character(len=:), allocatable :: line
        !! name of the file to assign the file handle to

    line = read_file(fTempSIM, fTempSIM_iostat)
end function fTempSIM_read


subroutine fTempSIM_close()
    close(fTempSIM)
end subroutine fTempSIM_close


! fRainSIM

subroutine fRainSIM_open(filename, mode)
    !! Opens the given file, assigning it to the 'fRainSIM' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')
    call open_file(fRainSIM, filename, mode, fRainSIM_iostat)
end subroutine fRainSIM_open


function fRainSIM_read() result(line)
    !! Returns the next line read from the 'fRainSIM' file.
    character(len=:), allocatable :: line
        !! name of the file to assign the file handle to

    line = read_file(fRainSIM, fRainSIM_iostat)
end function fRainSIM_read


subroutine fRainSIM_close()
    close(fRainSIM)
end subroutine fRainSIM_close


! fCuts

subroutine fCuts_open(filename, mode)
    !! Opens the given file, assigning it to the 'fCuts' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')
    call open_file(fCuts, filename, mode, fCuts_iostat)
end subroutine fCuts_open


function fCuts_read() result(line)
    !! Returns the next line read from the 'fCuts' file.
    character(len=:), allocatable :: line
        !! name of the file to assign the file handle to

    line = read_file(fCuts, fCuts_iostat)
end function fCuts_read


function fCuts_eof() result(eof)
    !! Returns whether the end of the 'fCuts' file has been reached.
    logical :: eof

    eof = fCuts_iostat == iostat_end
end function fCuts_eof


subroutine fCuts_close()
    close(fCuts)
end subroutine fCuts_close


! fObs

subroutine fObs_open(filename, mode)
    !! Opens the given file, assigning it to the 'fObs' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')

    call open_file(fObs, filename, mode, fObs_iostat)
end subroutine fObs_open


function fObs_read() result(line)
    !! Returns the next line read from the 'fObs' file.
    character(len=:), allocatable :: line
        !! name of the file to assign the file handle to

    line = read_file(fObs, fObs_iostat)
end function fObs_read


function fObs_eof() result(eof)
    !! Returns whether the end of the 'fObs' file has been reached.
    logical :: eof

    eof = fObs_iostat == iostat_end
end function fObs_eof


subroutine fObs_close()
    close(fObs)
end subroutine fObs_close


subroutine fObs_rewind()
    rewind(fObs)
end subroutine fObs_rewind


! fHarvest

function GetfHarvest_filename() result(str)
    !! Getter for the "fHarvest_filename" global variable.
    character(len=:), allocatable :: str

    str = fHarvest_filename
end function GetfHarvest_filename


subroutine SetfHarvest_filename(str)
    !! Setter for the "fHarvest_filename" global variable.
    character(len=*), intent(in) :: str

    fHarvest_filename = str
end subroutine SetfHarvest_filename


subroutine fHarvest_open(filename, mode)
    !! Opens the given file, assigning it to the 'fHarvest' file handle.
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')

    call open_file(fHarvest, filename, mode, fHarvest_iostat)
end subroutine fHarvest_open


subroutine fHarvest_write(line, advance_in)
    !! Writes the given line to the fHarvest file.
    character(len=*), intent(in) :: line
        !! line to write
    logical, intent(in), optional :: advance_in
        !! whether or not to append a newline character

    logical :: advance

    if (present(advance_in)) then
        advance = advance_in
    else
        advance = .true.
    end if
    call write_file(fHarvest, line, advance, fHarvest_iostat)

    if (advance) then
        ! For some reason (a compiler bug?) we need to explicitly flush
        ! the buffer here. Otherwise, with GNU Fortran v10.2.1 and DEBUG=0,
        ! the last line of test_defaultPROharvests.OUT in the Harvest test
        ! does not get written. This does not occur when either DEBUG=1 is
        ! applied or when using GNU Fortran v6.4.0.
        call flush(fHarvest)
    end if
end subroutine fHarvest_write


subroutine fHarvest_close()
    close(fHarvest)
end subroutine fHarvest_close


! Bin

real(dp) function GetBin()
    !! Getter for the "Bin" global variable.

    GetBin = Bin
end function GetBin


subroutine SetBin(Bin_in)
    !! Setter for the "Bin" global variable.
    real(dp), intent(in) :: Bin_in

    Bin = Bin_in
end subroutine SetBin


! Bout

real(dp) function GetBout()
    !! Getter for the "Bout" global variable.

    GetBout = Bout
end function GetBout


subroutine SetBout(Bout_in)
    !! Setter for the "Bout" global variable.
    real(dp), intent(in) :: Bout_in

    Bout = Bout_in
end subroutine SetBout


! GDDayi

real(dp) function GetGDDayi()
    !! Getter for the "GDDayi" global variable.

    GetGDDayi = GDDayi
end function GetGDDayi


subroutine SetGDDayi(GDDayi_in)
    !! Setter for the "GDDayi" global variable.
    real(dp), intent(in) :: GDDayi_in

    GDDayi = GDDayi_in
end subroutine SetGDDayi


! FracBiomass

real(dp) function GetFracBiomassPotSF()
    !! Getter for the "FracBiomassPotSF" global variable.

    GetFracBiomassPotSF = FracBiomassPotSF
end function GetFracBiomassPotSF


subroutine SetFracBiomassPotSF(FracBiomassPotSF_in)
    !! Setter for the "FracBiomassPotSF" global variable.
    real(dp), intent(in) :: FracBiomassPotSF_in

    FracBiomassPotSF = FracBiomassPotSF_in
end subroutine SetFracBiomassPotSF


! CO2i

real(dp) function GetCO2i()
    !! Getter for the "CO2i" global variable.

    GetCO2i = CO2i
end function GetCO2i


subroutine SetCO2i(CO2i_in)
    !! Setter for the "CO2i" global variable.
    real(dp), intent(in) :: CO2i_in

    CO2i = CO2i_in
end subroutine SetCO2i


! PreviousSum

type(rep_sum) function GetPreviousSum()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum = PreviousSum
end function GetPreviousSum


real(dp) function GetPreviousSum_Epot()
    !! Getter for the "PreviousSum" global variable.

     GetPreviousSum_Epot = PreviousSum%Epot
end function GetPreviousSum_Epot


real(dp) function GetPreviousSum_Tpot()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_Tpot = PreviousSum%Tpot
end function GetPreviousSum_Tpot


real(dp) function GetPreviousSum_Rain()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_Rain = PreviousSum%Rain
end function GetPreviousSum_Rain


real(dp) function GetPreviousSum_Irrigation()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_Irrigation = PreviousSum%Irrigation
end function GetPreviousSum_Irrigation


real(dp) function GetPreviousSum_Infiltrated()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_Infiltrated = PreviousSum%Infiltrated
end function GetPreviousSum_Infiltrated


real(dp) function GetPreviousSum_Runoff()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_Runoff = PreviousSum%Runoff
end function GetPreviousSum_Runoff


real(dp) function GetPreviousSum_Drain()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_Drain = PreviousSum%Drain
end function GetPreviousSum_Drain


real(dp) function GetPreviousSum_Eact()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_Eact = PreviousSum%Eact
end function GetPreviousSum_Eact


real(dp) function GetPreviousSum_Tact()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_Tact = PreviousSum%Tact
end function GetPreviousSum_Tact


real(dp) function GetPreviousSum_TrW()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_TrW = PreviousSum%TrW
end function GetPreviousSum_TrW


real(dp) function GetPreviousSum_ECropCycle()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_ECropCycle = PreviousSum%ECropCycle
end function GetPreviousSum_ECropCycle


real(dp) function GetPreviousSum_CRwater()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_CRwater = PreviousSum%CRwater
end function GetPreviousSum_CRwater


real(dp) function GetPreviousSum_Biomass()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_Biomass = PreviousSum%Biomass
end function GetPreviousSum_Biomass


real(dp) function GetPreviousSum_YieldPart()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_YieldPart = PreviousSum%YieldPart
end function GetPreviousSum_YieldPart


real(dp) function GetPreviousSum_BiomassPot()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_BiomassPot = PreviousSum%BiomassPot
end function GetPreviousSum_BiomassPot


real(dp) function GetPreviousSum_BiomassUnlim()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_BiomassUnlim = PreviousSum%BiomassUnlim
end function GetPreviousSum_BiomassUnlim


real(dp) function GetPreviousSum_BiomassTot()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_BiomassTot = PreviousSum%BiomassTot
end function GetPreviousSum_BiomassTot


real(dp) function GetPreviousSum_SaltIn()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_SaltIn = PreviousSum%SaltIn
end function GetPreviousSum_SaltIn


real(dp) function GetPreviousSum_SaltOut()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_SaltOut = PreviousSum%SaltOut
end function GetPreviousSum_SaltOut


real(dp) function GetPreviousSum_CRSalt()
    !! Getter for the "PreviousSum" global variable.

    GetPreviousSum_CRSalt = PreviousSum%CRSalt
end function GetPreviousSum_CRSalt


subroutine SetPreviousSum(PreviousSum_in)
    !! Setter for the "PreviousSum" global variable.
    type(rep_sum), intent(in) :: PreviousSum_in

    PreviousSum = PreviousSum_in
end subroutine SetPreviousSum


subroutine SetPreviousSum_Epot(Epot)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Epot

    PreviousSum%Epot = Epot
end subroutine SetPreviousSum_Epot


subroutine SetPreviousSum_Tpot(Tpot)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Tpot

    PreviousSum%Tpot = Tpot
end subroutine SetPreviousSum_Tpot


subroutine SetPreviousSum_Rain(Rain)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Rain

    PreviousSum%Rain = Rain
end subroutine SetPreviousSum_Rain


subroutine SetPreviousSum_Irrigation(Irrigation)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Irrigation

    PreviousSum%Irrigation = Irrigation
end subroutine SetPreviousSum_Irrigation


subroutine SetPreviousSum_Infiltrated(Infiltrated)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Infiltrated

    PreviousSum%Infiltrated = Infiltrated
end subroutine SetPreviousSum_Infiltrated


subroutine SetPreviousSum_Runoff(Runoff)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Runoff

    PreviousSum%Runoff = Runoff
end subroutine SetPreviousSum_Runoff


subroutine SetPreviousSum_Drain(Drain)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Drain

    PreviousSum%Drain = Drain
end subroutine SetPreviousSum_Drain


subroutine SetPreviousSum_Eact(Eact)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Eact

    PreviousSum%Eact = Eact
end subroutine SetPreviousSum_Eact


subroutine SetPreviousSum_Tact(Tact)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Tact

    PreviousSum%Tact = Tact
end subroutine SetPreviousSum_Tact


subroutine SetPreviousSum_TrW(TrW)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: TrW

    PreviousSum%TrW = TrW
end subroutine SetPreviousSum_TrW


subroutine SetPreviousSum_ECropCycle(ECropCycle)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: ECropCycle

    PreviousSum%ECropCycle = ECropCycle
end subroutine SetPreviousSum_ECropCycle


subroutine SetPreviousSum_CRwater(CRwater)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: CRwater

    PreviousSum%CRwater = CRwater
end subroutine SetPreviousSum_CRwater


subroutine SetPreviousSum_Biomass(Biomass)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: Biomass

    PreviousSum%Biomass = Biomass
end subroutine SetPreviousSum_Biomass


subroutine SetPreviousSum_YieldPart(YieldPart)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: YieldPart

    PreviousSum%YieldPart = YieldPart
end subroutine SetPreviousSum_YieldPart


subroutine SetPreviousSum_BiomassPot(BiomassPot)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: BiomassPot

    PreviousSum%BiomassPot = BiomassPot
end subroutine SetPreviousSum_BiomassPot


subroutine SetPreviousSum_BiomassUnlim(BiomassUnlim)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: BiomassUnlim

    PreviousSum%BiomassUnlim = BiomassUnlim
end subroutine SetPreviousSum_BiomassUnlim


subroutine SetPreviousSum_BiomassTot(BiomassTot)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: BiomassTot

    PreviousSum%BiomassTot = BiomassTot
end subroutine SetPreviousSum_BiomassTot


subroutine SetPreviousSum_SaltIn(SaltIn)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: SaltIn

    PreviousSum%SaltIn = SaltIn
end subroutine SetPreviousSum_SaltIn


subroutine SetPreviousSum_SaltOut(SaltOut)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: SaltOut

    PreviousSum%SaltOut = SaltOut
end subroutine SetPreviousSum_SaltOut


subroutine SetPreviousSum_CRSalt(CRSalt)
    !! Setter for the "PreviousSum" global variable.
    real(dp), intent(in) :: CRSalt

    PreviousSum%CRSalt = CRSalt
end subroutine SetPreviousSum_CRSalt

! GwTable

integer(int32) function GetGwTable_DNr1()
    !! Getter for the "GetGwTable" global variable.

    GetGwTable_DNr1 = GwTable%DNr1
end function GetGwTable_DNr1


integer(int32) function GetGwTable_DNr2()
    !! Getter for the "GetGwTable" global variable.

    GetGwTable_DNr2 = GwTable%DNr2
end function GetGwTable_DNr2


integer(int32) function GetGwTable_Z1()
    !! Getter for the "GetGwTable" global variable.

    GetGwTable_Z1 = GwTable%Z1
end function GetGwTable_Z1


integer(int32) function GetGwTable_Z2()
    !! Getter for the "GetGwTable" global variable.

    GetGwTable_Z2 = GwTable%Z2
end function GetGwTable_Z2


real(dp) function GetGwTable_EC1()
    !! Getter for the "GetGwTable" global variable.

    GetGwTable_EC1 = GwTable%EC1
end function GetGwTable_EC1


real(dp) function GetGwTable_EC2()
    !! Getter for the "GetGwTable" global variable.

    GetGwTable_EC2 = GwTable%EC2
end function GetGwTable_EC2


type(rep_GwTable) function GetGwTable()
    !! Getter for the "GetGwTable" global variable.

    GetGwTable = GwTable
end function GetGwTable


subroutine SetGwTable(GwT)
    !! Setter for the "GwTable" global variable.
    type(rep_GWTable), intent(in) :: GwT

    GwTable = GwT
end subroutine SetGwTable


subroutine SetGwTable_DNr1(DNr1)
    !! Setter for the "GwTable" global variable.
    integer(int32), intent(in) :: DNr1

    GwTable%DNr1 = DNr1
end subroutine SetGwTable_DNr1


subroutine SetGwTable_DNr2(DNr2)
    !! Setter for the "GwTable" global variable.
    integer(int32), intent(in) :: DNr2

    GwTable%DNr2 = DNr2
end subroutine SetGwTable_DNr2


subroutine SetGwTable_Z1(Z1)
    !! Setter for the "GwTable" global variable.
    integer(int32), intent(in) :: Z1

    GwTable%Z1 = Z1
end subroutine SetGwTable_Z1


subroutine SetGwTable_Z2(Z2)
    !! Setter for the "GwTable" global variable.
    integer(int32), intent(in) :: Z2

    GwTable%Z2 = Z2
end subroutine SetGwTable_Z2


subroutine SetGwTable_EC1(EC1)
    !! Setter for the "GwTable" global variable.
    real(dp), intent(in) :: EC1

    GwTable%EC1 = EC1
end subroutine SetGwTable_EC1


subroutine SetGwTable_EC2(EC2)
    !! Setter for the "GwTable" global variable.
    real(dp), intent(in) :: EC2

    GwTable%EC2 = EC2
end subroutine SetGwTable_EC2

! PlotVarCrop

type(rep_plotPar) function GetPlotVarCrop()
    !! Getter for the "PlotVarCrop" global variable.

    GetPlotVarCrop = PlotVarCrop
end function GetPlotVarCrop


subroutine SetPlotVarCrop(PlotVarCrop_in)
    !! Setter for the "PlotVarCrop" global variable.
    type(rep_plotPar), intent(in) :: PlotVarCrop_in

    PlotVarCrop = PlotVarCrop_in
end subroutine SetPlotVarCrop


subroutine SetPlotVarCrop_PotVal(PotVal)
    !! Setter for the "PlotVarCrop" global variable.
    real(dp), intent(in) :: PotVal

    PlotVarCrop%PotVal = PotVal
end subroutine SetPlotVarCrop_PotVal


subroutine SetPlotVarCrop_ActVal(ActVal)
    !! Setter for the "PlotVarCrop" global variable.
    real(dp), intent(in) :: ActVal

    PlotVarCrop%ActVal = ActVal
end subroutine SetPlotVarCrop_ActVal


real(dp) function GetPlotVarCrop_ActVal()
    !! Getter for the "PlotVarCrop_ActVal" global variable.

    GetPlotVarCrop_ActVal = PlotVarCrop%ActVal
end function GetPlotVarCrop_ActVal


real(dp) function GetPlotVarCrop_PotVal()
    !! Getter for the "PlotVarCrop_PotVal" global variable.

    GetPlotVarCrop_PotVal = PlotVarCrop%PotVal
end function GetPlotVarCrop_PotVal

! IrriInfoRecord1

function GetIrriInfoRecord1() result(IrriInfoRecord1_out)
    !! Getter for the "IrriInfoRecord1" global variable.
    type(repIrriInfoRecord) :: IrriInfoRecord1_out

    IrriInfoRecord1_out = IrriInfoRecord1
end function GetIrriInfoRecord1


subroutine SetIrriInfoRecord1(IrriInfoRecord1_in)
    !! Setter for the "TminDatSet" global variable.
    type(repIrriInfoRecord), intent(in) :: IrriInfoRecord1_in

    IrriInfoRecord1 = IrriInfoRecord1_in
end subroutine SetIrriInfoRecord1


logical function GetIrriInfoRecord1_NoMoreInfo()
    !! Getter for the "IrriInfoRecord1" global variable.

    GetIrriInfoRecord1_NoMoreInfo = IrriInfoRecord1%NoMoreInfo
end function GetIrriInfoRecord1_NoMoreInfo


integer(int32) function GetIrriInfoRecord1_FromDay()
    !! Getter for the "IrriInfoRecord1" global variable.

    GetIrriInfoRecord1_FromDay = IrriInfoRecord1%FromDay
end function GetIrriInfoRecord1_FromDay


integer(int32) function GetIrriInfoRecord1_ToDay()
    !! Getter for the "IrriInfoRecord1" global variable.

    GetIrriInfoRecord1_ToDay = IrriInfoRecord1%ToDay
end function GetIrriInfoRecord1_ToDay


integer(int32) function GetIrriInfoRecord1_TimeInfo()
    !! Getter for the "IrriInfoRecord1" global variable.

    GetIrriInfoRecord1_TimeInfo = IrriInfoRecord1%TimeInfo
end function GetIrriInfoRecord1_TimeInfo


integer(int32) function GetIrriInfoRecord1_DepthInfo()
    !! Getter for the "IrriInfoRecord1" global variable.

    GetIrriInfoRecord1_DepthInfo = IrriInfoRecord1%DepthInfo
end function GetIrriInfoRecord1_DepthInfo


subroutine SetIrriInfoRecord1_NoMoreInfo(NoMoreInfo)
    !! Setter for the "IrriInfoRecord1" global variable.
    logical, intent(in) :: NoMoreInfo

    IrriInfoRecord1%NoMoreInfo = NoMoreInfo
end subroutine SetIrriInfoRecord1_NoMoreInfo


subroutine SetIrriInfoRecord1_FromDay(FromDay)
    !! Setter for the "IrriInfoRecord1" global variable.
    integer(int32), intent(in) :: FromDay

    IrriInfoRecord1%FromDay = FromDay
end subroutine SetIrriInfoRecord1_FromDay


subroutine SetIrriInfoRecord1_ToDay(ToDay)
    !! Setter for the "IrriInfoRecord1" global variable.
    integer(int32), intent(in) :: ToDay

    IrriInfoRecord1%ToDay = ToDay
end subroutine SetIrriInfoRecord1_ToDay


subroutine SetIrriInfoRecord1_TimeInfo(TimeInfo)
    !! Setter for the "IrriInfoRecord1" global variable.
    integer(int32), intent(in) :: TimeInfo

    IrriInfoRecord1%TimeInfo = TimeInfo
end subroutine SetIrriInfoRecord1_TimeInfo


subroutine SetIrriInfoRecord1_DepthInfo(DepthInfo)
    !! Setter for the "IrriInfoRecord1" global variable.
    integer(int32), intent(in) :: DepthInfo

    IrriInfoRecord1%DepthInfo = DepthInfo
end subroutine SetIrriInfoRecord1_DepthInfo

! IrriInfoRecord2

function GetIrriInfoRecord2() result(IrriInfoRecord2_out)
    !! Getter for the "IrriInfoRecord2" global variable.
    type(repIrriInfoRecord) :: IrriInfoRecord2_out

    IrriInfoRecord2_out = IrriInfoRecord2
end function GetIrriInfoRecord2


subroutine SetIrriInfoRecord2(IrriInfoRecord2_in)
    !! Setter for the "TminDatSet" global variable.
    type(repIrriInfoRecord), intent(in) :: IrriInfoRecord2_in

    IrriInfoRecord2 = IrriInfoRecord2_in
end subroutine SetIrriInfoRecord2


logical function GetIrriInfoRecord2_NoMoreInfo()
    !! Getter for the "IrriInfoRecord2" global variable.

    GetIrriInfoRecord2_NoMoreInfo = IrriInfoRecord2%NoMoreInfo
end function GetIrriInfoRecord2_NoMoreInfo


integer(int32) function GetIrriInfoRecord2_FromDay()
    !! Getter for the "IrriInfoRecord2" global variable.

    GetIrriInfoRecord2_FromDay = IrriInfoRecord2%FromDay
end function GetIrriInfoRecord2_FromDay


integer(int32) function GetIrriInfoRecord2_ToDay()
    !! Getter for the "IrriInfoRecord2" global variable.

    GetIrriInfoRecord2_ToDay = IrriInfoRecord2%ToDay
end function GetIrriInfoRecord2_ToDay


integer(int32) function GetIrriInfoRecord2_TimeInfo()
    !! Getter for the "IrriInfoRecord2" global variable.

    GetIrriInfoRecord2_TimeInfo = IrriInfoRecord2%TimeInfo
end function GetIrriInfoRecord2_TimeInfo


integer(int32) function GetIrriInfoRecord2_DepthInfo()
    !! Getter for the "IrriInfoRecord2" global variable.

    GetIrriInfoRecord2_DepthInfo = IrriInfoRecord2%DepthInfo
end function GetIrriInfoRecord2_DepthInfo


subroutine SetIrriInfoRecord2_NoMoreInfo(NoMoreInfo)
    !! Setter for the "IrriInfoRecord2" global variable.
    logical, intent(in) :: NoMoreInfo

    IrriInfoRecord2%NoMoreInfo = NoMoreInfo
end subroutine SetIrriInfoRecord2_NoMoreInfo


subroutine SetIrriInfoRecord2_FromDay(FromDay)
    !! Setter for the "IrriInfoRecord2" global variable.
    integer(int32), intent(in) :: FromDay

    IrriInfoRecord2%FromDay = FromDay
end subroutine SetIrriInfoRecord2_FromDay


subroutine SetIrriInfoRecord2_ToDay(ToDay)
    !! Setter for the "IrriInfoRecord2" global variable.
    integer(int32), intent(in) :: ToDay

    IrriInfoRecord2%ToDay = ToDay
end subroutine SetIrriInfoRecord2_ToDay


subroutine SetIrriInfoRecord2_TimeInfo(TimeInfo)
    !! Setter for the "IrriInfoRecord2" global variable.
    integer(int32), intent(in) :: TimeInfo

    IrriInfoRecord2%TimeInfo = TimeInfo
end subroutine SetIrriInfoRecord2_TimeInfo


subroutine SetIrriInfoRecord2_DepthInfo(DepthInfo)
    !! Setter for the "IrriInfoRecord2" global variable.
    integer(int32), intent(in) :: DepthInfo

    IrriInfoRecord2%DepthInfo = DepthInfo
end subroutine SetIrriInfoRecord2_DepthInfo

! StressTot

type(rep_StressTot) function GetStressTot()
    !! Getter for the "StressTot" global variable.

    GetStressTot = StressTot
end function GetStressTot


real(dp) function GetStressTot_Salt()
    !! Getter for the "StressTot" global variable.

    GetStressTot_Salt = StressTot%Salt
end function GetStressTot_Salt


real(dp) function GetStressTot_Temp()
    !! Getter for the "StressTot" global variable.

    GetStressTot_Temp = StressTot%Temp
end function GetStressTot_Temp


real(dp) function GetStressTot_Exp()
    !! Getter for the "StressTot" global variable.

    GetStressTot_Exp = StressTot%Exp
end function GetStressTot_Exp


real(dp) function GetStressTot_Sto()
    !! Getter for the "StressTot" global variable.

    GetStressTot_Sto = StressTot%Sto
end function GetStressTot_Sto


real(dp) function GetStressTot_Weed()
    !! Getter for the "StressTot" global variable.

    GetStressTot_Weed = StressTot%Weed
end function GetStressTot_Weed


integer(int32) function GetStressTot_NrD()
    !! Getter for the "StressTot" global variable.

    GetStressTot_NrD = StressTot%NrD
end function GetStressTot_NrD


subroutine SetStressTot(StressTot_in)
    !! Setter for the "StressTot" global variable.
    type(rep_StressTot), intent(in) :: StressTot_in

    StressTot = StressTot_in
end subroutine SetStressTot


subroutine SetStressTot_Salt(Salt)
    !! Setter for the "StressTot" global variable.
    real(dp), intent(in) :: Salt

    StressTot%Salt = Salt
end subroutine SetStressTot_Salt


subroutine SetStressTot_Temp(Temp)
    !! Setter for the "StressTot" global variable.
    real(dp), intent(in) :: Temp

    StressTot%Temp = Temp
end subroutine SetStressTot_Temp


subroutine SetStressTot_Exp(Exp)
    !! Setter for the "StressTot" global variable.
    real(dp), intent(in) :: Exp

    StressTot%Exp = Exp
end subroutine SetStressTot_Exp


subroutine SetStressTot_Sto(Sto)
    !! Setter for the "StressTot" global variable.
    real(dp), intent(in) :: Sto

    StressTot%Sto = Sto
end subroutine SetStressTot_Sto


subroutine SetStressTot_Weed(Weed)
    !! Setter for the "StressTot" global variable.
    real(dp), intent(in) :: Weed

    StressTot%Weed = Weed
end subroutine SetStressTot_Weed


subroutine SetStressTot_NrD(NrD)
    !! Setter for the "StressTot" global variable.
    integer(int32), intent(in) :: NrD

    StressTot%NrD = NrD
end subroutine SetStressTot_NrD

! CutInfoRecord1

type(repCutInfoRecord) function GetCutInfoRecord1()
    !! Getter for the "CutInfoRecord1" global variable.

    GetCutInfoRecord1 = CutInfoRecord1
end function GetCutInfoRecord1


logical function GetCutInfoRecord1_NoMoreInfo()
    !! Getter for the "CutInfoRecord1" global variable.

    GetCutInfoRecord1_NoMoreInfo = CutInfoRecord1%NoMoreInfo
end function GetCutInfoRecord1_NoMoreInfo


integer(int32) function GetCutInfoRecord1_FromDay()
    !! Getter for the "CutInfoRecord1" global variable.

    GetCutInfoRecord1_FromDay = CutInfoRecord1%FromDay
end function GetCutInfoRecord1_FromDay


integer(int32) function GetCutInfoRecord1_ToDay()
    !! Getter for the "CutInfoRecord1" global variable.

    GetCutInfoRecord1_ToDay = CutInfoRecord1%ToDay
end function GetCutInfoRecord1_ToDay


integer(int32) function GetCutInfoRecord1_IntervalInfo()
    !! Getter for the "CutInfoRecord1" global variable.

    GetCutInfoRecord1_IntervalInfo = CutInfoRecord1%IntervalInfo
end function GetCutInfoRecord1_IntervalInfo


real(dp) function GetCutInfoRecord1_IntervalGDD()
    !! Getter for the "CutInfoRecord1" global variable.

    GetCutInfoRecord1_IntervalGDD = CutInfoRecord1%IntervalGDD
end function GetCutInfoRecord1_IntervalGDD


real(dp) function GetCutInfoRecord1_MassInfo()
    !! Getter for the "CutInfoRecord1" global variable.

    GetCutInfoRecord1_MassInfo = CutInfoRecord1%MassInfo
end function GetCutInfoRecord1_MassInfo


subroutine SetCutInfoRecord1(CutInfoRecord1_in)
    !! Setter for the "CutInfoRecord1" global variable.
    type(repCutInfoRecord), intent(in) :: CutInfoRecord1_in

    CutInfoRecord1 = CutInfoRecord1_in
end subroutine SetCutInfoRecord1


subroutine SetCutInfoRecord1_NoMoreInfo(NoMoreInfo)
    !! Setter for the "CutInfoRecord1" global variable.
    logical, intent(in) :: NoMoreInfo

    CutInfoRecord1%NoMoreInfo = NoMoreInfo
end subroutine SetCutInfoRecord1_NoMoreInfo


subroutine SetCutInfoRecord1_FromDay(FromDay)
    !! Setter for the "CutInfoRecord1" global variable.
    integer(int32), intent(in) :: FromDay

    CutInfoRecord1%FromDay = FromDay
end subroutine SetCutInfoRecord1_FromDay


subroutine SetCutInfoRecord1_ToDay(ToDay)
    !! Setter for the "CutInfoRecord1" global variable.
    integer(int32), intent(in) :: ToDay

    CutInfoRecord1%ToDay = ToDay
end subroutine SetCutInfoRecord1_ToDay


subroutine SetCutInfoRecord1_IntervalInfo(IntervalInfo)
    !! Setter for the "CutInfoRecord1" global variable.
    integer(int32), intent(in) :: IntervalInfo

    CutInfoRecord1%IntervalInfo = IntervalInfo
end subroutine SetCutInfoRecord1_IntervalInfo


subroutine SetCutInfoRecord1_IntervalGDD(IntervalGDD)
    !! Setter for the "CutInfoRecord1" global variable.
    real(dp), intent(in) :: IntervalGDD

    CutInfoRecord1%IntervalGDD = IntervalGDD
end subroutine SetCutInfoRecord1_IntervalGDD


subroutine SetCutInfoRecord1_MassInfo(MassInfo)
    !! Setter for the "CutInfoRecord1" global variable.
    real(dp), intent(in) :: MassInfo

    CutInfoRecord1%MassInfo = MassInfo
end subroutine SetCutInfoRecord1_MassInfo


! CutInfoRecord2

type(repCutInfoRecord) function GetCutInfoRecord2()
    !! Getter for the "CutInfoRecord2" global variable.

    GetCutInfoRecord2 = CutInfoRecord2
end function GetCutInfoRecord2


logical function GetCutInfoRecord2_NoMoreInfo()
    !! Getter for the "CutInfoRecord2" global variable.

    GetCutInfoRecord2_NoMoreInfo = CutInfoRecord2%NoMoreInfo
end function GetCutInfoRecord2_NoMoreInfo


integer(int32) function GetCutInfoRecord2_FromDay()
    !! Getter for the "CutInfoRecord2" global variable.

    GetCutInfoRecord2_FromDay = CutInfoRecord2%FromDay
end function GetCutInfoRecord2_FromDay


integer(int32) function GetCutInfoRecord2_ToDay()
    !! Getter for the "CutInfoRecord2" global variable.

    GetCutInfoRecord2_ToDay = CutInfoRecord2%ToDay
end function GetCutInfoRecord2_ToDay


integer(int32) function GetCutInfoRecord2_IntervalInfo()
    !! Getter for the "CutInfoRecord2" global variable.

    GetCutInfoRecord2_IntervalInfo = CutInfoRecord2%IntervalInfo
end function GetCutInfoRecord2_IntervalInfo


real(dp) function GetCutInfoRecord2_IntervalGDD()
    !! Getter for the "CutInfoRecord2" global variable.

    GetCutInfoRecord2_IntervalGDD = CutInfoRecord2%IntervalGDD
end function GetCutInfoRecord2_IntervalGDD


real(dp) function GetCutInfoRecord2_MassInfo()
    !! Getter for the "CutInfoRecord2" global variable.

    GetCutInfoRecord2_MassInfo = CutInfoRecord2%MassInfo
end function GetCutInfoRecord2_MassInfo


subroutine SetCutInfoRecord2(CutInfoRecord2_in)
    !! Setter for the "CutInfoRecord2" global variable.
    type(repCutInfoRecord), intent(in) :: CutInfoRecord2_in

    CutInfoRecord2 = CutInfoRecord2_in
end subroutine SetCutInfoRecord2


subroutine SetCutInfoRecord2_NoMoreInfo(NoMoreInfo)
    !! Setter for the "CutInfoRecord2" global variable.
    logical, intent(in) :: NoMoreInfo

    CutInfoRecord2%NoMoreInfo = NoMoreInfo
end subroutine SetCutInfoRecord2_NoMoreInfo


subroutine SetCutInfoRecord2_FromDay(FromDay)
    !! Setter for the "CutInfoRecord2" global variable.
    integer(int32), intent(in) :: FromDay

    CutInfoRecord2%FromDay = FromDay
end subroutine SetCutInfoRecord2_FromDay


subroutine SetCutInfoRecord2_ToDay(ToDay)
    !! Setter for the "CutInfoRecord2" global variable.
    integer(int32), intent(in) :: ToDay

    CutInfoRecord2%ToDay = ToDay
end subroutine SetCutInfoRecord2_ToDay


subroutine SetCutInfoRecord2_IntervalInfo(IntervalInfo)
    !! Setter for the "CutInfoRecord2" global variable.
    integer(int32), intent(in) :: IntervalInfo

    CutInfoRecord2%IntervalInfo = IntervalInfo
end subroutine SetCutInfoRecord2_IntervalInfo


subroutine SetCutInfoRecord2_IntervalGDD(IntervalGDD)
    !! Setter for the "CutInfoRecord2" global variable.
    real(dp), intent(in) :: IntervalGDD

    CutInfoRecord2%IntervalGDD = IntervalGDD
end subroutine SetCutInfoRecord2_IntervalGDD


subroutine SetCutInfoRecord2_MassInfo(MassInfo)
    !! Setter for the "CutInfoRecord2" global variable.
    real(dp), intent(in) :: MassInfo

    CutInfoRecord2%MassInfo = MassInfo
end subroutine SetCutInfoRecord2_MassInfo

! Transfer

type(rep_Transfer) function GetTransfer()
    !! Getter for the "Transfer" global variable.

    GetTransfer = Transfer
end function GetTransfer


logical function GetTransfer_Store()
    !! Getter for the "Transfer" global variable.

    GetTransfer_Store = Transfer%Store
end function GetTransfer_Store


logical function GetTransfer_Mobilize()
    !! Getter for the "Transfer" global variable.

    GetTransfer_Mobilize = Transfer%Mobilize
end function GetTransfer_Mobilize


real(dp) function GetTransfer_ToMobilize()
    !! Getter for the "Transfer" global variable.

    GetTransfer_ToMobilize = Transfer%ToMobilize
end function GetTransfer_ToMobilize


real(dp) function GetTransfer_Bmobilized()
    !! Getter for the "Transfer" global variable.

    GetTransfer_Bmobilized = Transfer%Bmobilized
end function GetTransfer_Bmobilized


subroutine SetTransfer(Transfer_in)
    !! Setter for the "Transfer" global variable.
    type(rep_Transfer), intent(in) :: Transfer_in

    Transfer = Transfer_in
end subroutine SetTransfer


subroutine SetTransfer_Store(Store)
    !! Setter for the "Transfer" global variable.
    logical, intent(in) :: Store

    Transfer%Store = Store
end subroutine SetTransfer_Store


subroutine SetTransfer_Mobilize(Mobilize)
    !! Setter for the "Transfer" global variable.
    logical, intent(in) :: Mobilize

    Transfer%Mobilize = Mobilize
end subroutine SetTransfer_Mobilize


subroutine SetTransfer_ToMobilize(ToMobilize)
    !! Setter for the "Transfer" global variable.
    real(dp), intent(in) :: ToMobilize

    Transfer%ToMobilize = ToMobilize
end subroutine SetTransfer_ToMobilize


subroutine SetTransfer_Bmobilized(Bmobilized)
    !! Setter for the "Transfer" global variable.
    real(dp), intent(in) :: Bmobilized

    Transfer%Bmobilized = Bmobilized
end subroutine SetTransfer_Bmobilized


!TminDatSet

function GetTminDataSet() result(TminDataSet_out)
    !! Getter for the "TminDataSet" global variable.
    type(rep_DayEventDbl), dimension(31) :: TminDataSet_out

    TminDataSet_out = TminDataSet
end function GetTminDataSet


function GetTminDataSet_i(i) result(TminDataSet_i)
    !! Getter for individual elements of the "TminDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl) :: TminDataSet_i

    TminDataSet_i = TminDataSet(i)
end function GetTminDataSet_i


integer(int32) function GetTminDataSet_DayNr(i)
    integer(int32), intent(in) :: i

    GetTminDataSet_DayNr = TminDataSet(i)%DayNr
end function GetTminDataSet_DayNr


real(dp) function GetTminDataSet_Param(i)
    integer(int32), intent(in) :: i

    GetTminDataSet_Param = TminDataSet(i)%Param
end function GetTminDataSet_Param


subroutine SetTminDataSet(TminDataSet_in)
    !! Setter for the "TminDatSet" global variable.
    type(rep_DayEventDbl), dimension(31), intent(in) :: TminDataSet_in

    TminDataSet = TminDataSet_in
end subroutine SetTminDataSet


subroutine SetTminDataSet_i(i, TminDataSet_i)
    !! Setter for individual element for the "TminDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl), intent(in) :: TminDataSet_i

    TminDataSet(i) = TminDataSet_i
end subroutine SetTminDataSet_i


subroutine SetTminDataSet_DayNr(i, DayNr_in)
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: DayNr_in

    TminDataSet(i)%DayNr = DayNr_in
end subroutine SetTminDataSet_DayNr


subroutine SetTminDataSet_Param(i, Param_in)
    integer(int32), intent(in) :: i
    real(dp), intent(in) :: Param_in

    TminDataSet(i)%Param = Param_in
end subroutine SetTminDataSet_Param

! TmaxDataSet

function GetTmaxDataSet() result(TmaxDataSet_out)
    !! Getter for the "TmaxDataSet" global variable.
    type(rep_DayEventDbl), dimension(31) :: TmaxDataSet_out

    TmaxDataSet_out = TmaxDataSet
end function GetTmaxDataSet


function GetTmaxDataSet_i(i) result(TmaxDataSet_i)
    !! Getter for individual elements of the "TmaxDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl) :: TmaxDataSet_i

    TmaxDataSet_i = TmaxDataSet(i)
end function GetTmaxDataSet_i


integer(int32) function GetTmaxDataSet_DayNr(i)
    integer(int32), intent(in) :: i

    GetTmaxDataSet_DayNr = TmaxDataSet(i)%DayNr
end function GetTmaxDataSet_DayNr


real(dp) function GetTmaxDataSet_Param(i)
    integer(int32), intent(in) :: i

    GetTmaxDataSet_Param = TmaxDataSet(i)%Param
end function GetTmaxDataSet_Param


subroutine SetTmaxDataSet(TmaxDataSet_in)
    !! Setter for the "TmaxDatSet" global variable.
    type(rep_DayEventDbl), dimension(31), intent(in) :: TmaxDataSet_in

    TmaxDataSet = TmaxDataSet_in
end subroutine SetTmaxDataSet


subroutine SetTmaxDataSet_i(i, TmaxDataSet_i)
    !! Setter for individual element for the "TmaxDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl), intent(in) :: TmaxDataSet_i

    TmaxDataSet(i) = TmaxDataSet_i
end subroutine SetTmaxDataSet_i


subroutine SetTmaxDataSet_DayNr(i, DayNr_in)
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: DayNr_in

    TmaxDataSet(i)%DayNr = DayNr_in
end subroutine SetTmaxDataSet_DayNr


subroutine SetTmaxDataSet_Param(i, Param_in)
    integer(int32), intent(in) :: i
    real(dp), intent(in) :: Param_in

    TmaxDataSet(i)%Param = Param_in
end subroutine SetTmaxDataSet_Param


integer(int32) function GetDayNri()

    GetDayNri = DayNri
end function GetDayNri


subroutine SetDayNri(DayNri_in)
    integer(int32), intent(in) :: DayNri_in

    DayNri = DayNri_in
end subroutine SetDayNri


! EToDataSet

function GetEToDataSet() result(EToDataSet_out)
    !! Getter for the "EToDataSet" global variable.
    type(rep_DayEventDbl), dimension(31) :: EToDataSet_out

    EToDataSet_out = EToDataSet
end function GetEToDataSet


function GetEToDataSet_i(i) result(EToDataSet_i)
    !! Getter for individual elements of the "EToDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl) :: EToDataSet_i

    EToDataSet_i = EToDataSet(i)
end function GetEToDataSet_i


integer(int32) function GetEToDataSet_DayNr(i)
    integer(int32), intent(in) :: i

    GetEToDataSet_DayNr = EToDataSet(i)%DayNr
end function GetEToDataSet_DayNr


real(dp) function GetEToDataSet_Param(i)
    integer(int32), intent(in) :: i

    GetEToDataSet_Param = EToDataSet(i)%Param
end function GetEToDataSet_Param


subroutine SetEToDataSet(EToDataSet_in)
    !! Setter for the "EToDatSet" global variable.
    type(rep_DayEventDbl), dimension(31), intent(in) :: EToDataSet_in

    EToDataSet = EToDataSet_in
end subroutine SetEToDataSet


subroutine SetEToDataSet_i(i, EToDataSet_i)
    !! Setter for individual element for the "EToDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl), intent(in) :: EToDataSet_i

    EToDataSet(i) = EToDataSet_i
end subroutine SetEToDataSet_i


subroutine SetEToDataSet_DayNr(i, DayNr_in)
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: DayNr_in

    EToDataSet(i)%DayNr = DayNr_in
end subroutine SetEToDataSet_DayNr


subroutine SetEToDataSet_Param(i, Param_in)
    integer(int32), intent(in) :: i
    real(dp), intent(in) :: Param_in

    EToDataSet(i)%Param = Param_in
end subroutine SetEToDataSet_Param

! RainDataSet

function GetRainDataSet() result(RainDataSet_out)
    !! Getter for the "RainDataSet" global variable.
    type(rep_DayEventDbl), dimension(31) :: RainDataSet_out

    RainDataSet_out = RainDataSet
end function GetRainDataSet


function GetRainDataSet_i(i) result(RainDataSet_i)
    !! Getter for individual elements of the "RainDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl) :: RainDataSet_i

    RainDataSet_i = RainDataSet(i)
end function GetRainDataSet_i


integer(int32) function GetRainDataSet_DayNr(i)
    integer(int32), intent(in) :: i

    GetRainDataSet_DayNr = RainDataSet(i)%DayNr
end function GetRainDataSet_DayNr


real(dp) function GetRainDataSet_Param(i)
    integer(int32), intent(in) :: i

    GetRainDataSet_Param = RainDataSet(i)%Param
end function GetRainDataSet_Param


subroutine SetRainDataSet(RainDataSet_in)
    !! Setter for the "RainDatSet" global variable.
    type(rep_DayEventDbl), dimension(31), intent(in) :: RainDataSet_in

    RainDataSet = RainDataSet_in
end subroutine SetRainDataSet


subroutine SetRainDataSet_i(i, RainDataSet_i)
    !! Setter for individual element for the "RainDataSet" global variable.
    integer(int32), intent(in) :: i
    type(rep_DayEventDbl), intent(in) :: RainDataSet_i

    RainDataSet(i) = RainDataSet_i
end subroutine SetRainDataSet_i


subroutine SetRainDataSet_DayNr(i, DayNr_in)
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: DayNr_in

    RainDataSet(i)%DayNr = DayNr_in
end subroutine SetRainDataSet_DayNr


subroutine SetRainDataSet_Param(i, Param_in)
    integer(int32), intent(in) :: i
    real(dp), intent(in) :: Param_in

    RainDataSet(i)%Param = Param_in
end subroutine SetRainDataSet_Param


logical function GetGlobalIrriECw()
    !! Getter for the GlobalIrriECw global variable

    GetGlobalIrriECw = GlobalIrriECw
end function GetGlobalIrriECw


subroutine SetGlobalIrriECw(GlobalIrriECw_in)
    !! Setter for the GlobalIrriECw global variable
    logical, intent(in) :: GlobalIrriECw_in

    GlobalIrriECw = GlobalIrriECw_in
end subroutine SetGlobalIrriECw


logical function GetWaterTableInProfile()
    !! Getter for the "WaterTableInProfile" global variable.

    GetWaterTableInProfile = WaterTableInProfile
end function GetWaterTableInProfile


subroutine SetWaterTableInProfile(WaterTableInProfile_in)
    !! Setter for the "WaterTableInProfile" global variable.
    logical, intent(in) :: WaterTableInProfile_in

    WaterTableInProfile = WaterTableInProfile_in
end subroutine SetWaterTableInProfile


logical function GetStartMode()
    !! Getter for the "StartMode" global variable.

    GetStartMode = StartMode
end function GetStartMode


subroutine SetStartMode(StartMode_in)
    !! Setter for the "StartMode" global variable.
    logical, intent(in) :: StartMode_in

    StartMode = StartMode_in
end subroutine SetStartMode


logical function GetNoMoreCrop()
    !! Getter for the "NoMoreCrop" global variable.

    GetNoMoreCrop = NoMoreCrop
end function GetNoMoreCrop


subroutine SetNoMoreCrop(NoMoreCrop_in)
    !! Setter for the "NoMoreCrop" global variable.
    logical, intent(in) :: NoMoreCrop_in

    NoMoreCrop = NoMoreCrop_in
end subroutine SetNoMoreCrop


integer(int32) function GetIrriInterval()
    !! Getter for the "IrriInterval" global variable.

    GetIrriInterval = IrriInterval
end function GetIrriInterval


subroutine SetIrriInterval(IrriInterval_in)
    !! Setter for the "IrriInterval" global variable.
    integer(int32), intent(in) :: IrriInterval_in

    IrriInterval = IrriInterval_in
end subroutine SetIrriInterval


integer(int32) function GetTadj()
    !! Getter for the "Tadj" global variable.

    GetTadj = Tadj
end function GetTadj


subroutine SetTadj(Tadj_in)
    !! Setter for the "Tadj" global variable.
    integer(int32), intent(in) :: Tadj_in

    Tadj = Tadj_in
end subroutine SetTadj


integer(int32) function GetGDDTadj()
    !! Getter for the "GDDTadj" global variable.

    GetGDDTadj = GDDTadj
end function GetGDDTadj


subroutine SetGDDTadj(GDDTadj_in)
    !! Setter for the "GDDTadj" global variable.
    integer(int32), intent(in) :: GDDTadj_in

    GDDTadj = GDDTadj_in
end subroutine SetGDDTadj


integer(int32) function GetDayLastCut()
    !! Getter for the "DayLastCut" global variable.

    GetDayLastCut = DayLastCut
end function GetDayLastCut


subroutine SetDayLastCut(DayLastCut_in)
    !! Setter for the "DayLastCut" global variable.
    integer(int32), intent(in) :: DayLastCut_in

    DayLastCut = DayLastCut_in
end subroutine SetDayLastCut


integer(int32) function GetNrCut()
    !! Getter for the "NrCut" global variable.

    GetNrCut = NrCut
end function GetNrCut


subroutine SetNrCut(NrCut_in)
    !! Setter for the "NrCut" global variable.
    integer(int32), intent(in) :: NrCut_in

    NrCut = NrCut_in
end subroutine SetNrCut


integer(int32) function GetSumInterval()
    !! Getter for the "SumInterval" global variable.

    GetSumInterval = SumInterval
end function GetSumInterval


subroutine SetSumInterval(SumInterval_in)
    !! Setter for the "SumInterval" global variable.
    integer(int32), intent(in) :: SumInterval_in

    SumInterval = SumInterval_in
end subroutine SetSumInterval


integer(int32) function GetPreviousStressLevel()
    !! Getter for the "PreviousStressLevel" global variable.

    GetPreviousStressLevel = int(PreviousStressLevel, kind=int32)
end function GetPreviousStressLevel


subroutine SetPreviousStressLevel(PreviousStressLevel_in)
    !! Setter for the "PreviousStressLevel" global variable.
    integer(int32), intent(in) :: PreviousStressLevel_in

    PreviousStressLevel = int(PreviousStressLevel_in, kind=int8)
end subroutine SetPreviousStressLevel


integer(int8) function GetStressSFadjNEW()
    !! Getter for the "StressSFadjNEW" global variable.

    GetStressSFadjNEW = int(StressSFadjNEW, kind=int32)
end function GetStressSFadjNEW


subroutine SetStressSFadjNEW(StressSFadjNEW_in)
    !! Setter for the "StressSFadjNEW" global variable.
    integer(int32), intent(in) :: StressSFadjNEW_in

    StressSFadjNEW = int(StressSFadjNEW_in, kind=int8)
end subroutine SetStressSFadjNEW


real(dp) function GetCCxWitheredTpotNoS()
    !! Getter for the "CCxWitheredTpotNoS" global variable.

    GetCCxWitheredTpotNoS = CCxWitheredTpotNoS
end function GetCCxWitheredTpotNoS


subroutine SetCCxWitheredTpotNoS(CCxWitheredTpotNoS_in)
    !! Setter for the "CCxWitheredTpotNoS" global variable.
    real(dp), intent(in) :: CCxWitheredTpotNoS_in

    CCxWitheredTpotNoS = CCxWitheredTpotNoS_in
end subroutine SetCCxWitheredTpotNoS


real(dp) function GetCoeffb0()
    !! Getter for the "Coeffb0" global variable.

    GetCoeffb0 = Coeffb0
end function GetCoeffb0


subroutine SetCoeffb0(Coeffb0_in)
    !! Setter for the "Coeffb0" global variable.
    real(dp), intent(in) :: Coeffb0_in

    Coeffb0 = Coeffb0_in
end subroutine SetCoeffb0


real(dp) function GetCoeffb1()
    !! Getter for the "Coeffb1" global variable.

    GetCoeffb1 = Coeffb1
end function GetCoeffb1


subroutine SetCoeffb1(Coeffb1_in)
    !! Setter for the "Coeffb1" global variable.
    real(dp), intent(in) :: Coeffb1_in

    Coeffb1 = Coeffb1_in
end subroutine SetCoeffb1


real(dp) function GetCoeffb2()
    !! Getter for the "Coeffb2" global variable.

    GetCoeffb2 = Coeffb2
end function GetCoeffb2


subroutine SetCoeffb2(Coeffb2_in)
    !! Setter for the "Coeffb2" global variable.
    real(dp), intent(in) :: Coeffb2_in

    Coeffb2 = Coeffb2_in
end subroutine SetCoeffb2


real(dp) function GetCoeffb0Salt()
    !! Getter for the "Coeffb0Salt" global variable.

    GetCoeffb0Salt = Coeffb0Salt
end function GetCoeffb0Salt


subroutine SetCoeffb0Salt(Coeffb0Salt_in)
    !! Setter for the "Coeffb0Salt" global variable.
    real(dp), intent(in) :: Coeffb0Salt_in

    Coeffb0Salt = Coeffb0Salt_in
end subroutine SetCoeffb0Salt


real(dp) function GetCoeffb1Salt()
    !! Getter for the "Coeffb1Salt" global variable.

    GetCoeffb1Salt = Coeffb1Salt
end function GetCoeffb1Salt


subroutine SetCoeffb1Salt(Coeffb1Salt_in)
    !! Setter for the "Coeffb1Salt" global variable.
    real(dp), intent(in) :: Coeffb1Salt_in

    Coeffb1Salt = Coeffb1Salt_in
end subroutine SetCoeffb1Salt


real(dp) function GetCoeffb2Salt()
    !! Getter for the "Coeffb2Salt" global variable.

    GetCoeffb2Salt = Coeffb2Salt
end function GetCoeffb2Salt


subroutine SetCoeffb2Salt(Coeffb2Salt_in)
    !! Setter for the "Coeffb2Salt" global variable.
    real(dp), intent(in) :: Coeffb2Salt_in

    Coeffb2Salt = Coeffb2Salt_in
end subroutine SetCoeffb2Salt


real(dp) function GetStressLeaf()
    !! Getter for the "StressLeaf" global variable.

    GetStressLeaf = StressLeaf
end function GetStressLeaf


subroutine SetStressLeaf(StressLeaf_in)
    !! Setter for the "StressLeaf" global variable.
    real(dp), intent(in) :: StressLeaf_in

    StressLeaf = StressLeaf_in
end subroutine SetStressLeaf


real(dp) function GetStressSenescence()
    !! Getter for the "StressSenescence" global variable.

    GetStressSenescence = StressSenescence
end function GetStressSenescence


subroutine SetStressSenescence(StressSenescence_in)
    !! Setter for the "StressSenescence" global variable.
    real(dp), intent(in) :: StressSenescence_in

    StressSenescence = StressSenescence_in
end subroutine SetStressSenescence


real(dp) function GetDayFraction()
    !! Getter for the "DayFraction" global variable.

    GetDayFraction = DayFraction
end function GetDayFraction


subroutine SetDayFraction(DayFraction_in)
    !! Setter for the "DayFraction" global variable.
    real(dp), intent(in) :: DayFraction_in

    DayFraction = DayFraction_in
end subroutine SetDayFraction


real(dp) function GetGDDayFraction()
    !! Getter for the "GDDayFraction" global variable.

    GetGDDayFraction = GDDayFraction
end function GetGDDayFraction


subroutine SetGDDayFraction(GDDayFraction_in)
    !! Setter for the "GDDayFraction" global variable.
    real(dp), intent(in) :: GDDayFraction_in

    GDDayFraction = GDDayFraction_in
end subroutine SetGDDayFraction


real(dp) function GetCGCref()
    !! Getter for the "CGCref" global variable.

    GetCGCref = CGCref
end function GetCGCref


subroutine SetCGCref(CGCref_in)
    !! Setter for the "CGCref" global variable.
    real(dp), intent(in) :: CGCref_in

    CGCref = CGCref_in
end subroutine SetCGCref


real(dp) function GetGDDCGCref()
    !! Getter for the "GDDCGCref" global variable.

    GetGDDCGCref = GDDCGCref
end function GetGDDCGCref


subroutine SetGDDCGCref(GDDCGCref_in)
    !! Setter for the "GDDCGCref" global variable.
    real(dp), intent(in) :: GDDCGCref_in

    GDDCGCref = GDDCGCref_in
end subroutine SetGDDCGCref


real(dp) function GetSumETo()
    !! Getter for the "SumETo" global variable.

    GetSumETo = SumETo
end function GetSumETo


subroutine SetSumETo(SumETo_in)
    !! Setter for the "SumETo" global variable.
    real(dp), intent(in) :: SumETo_in

    SumETo = SumETo_in
end subroutine SetSumETo


real(dp) function GetSumGDD()
    !! Getter for the "SumGDD" global variable.

    GetSumGDD = SumGDD
end function GetSumGDD


subroutine SetSumGDD(SumGDD_in)
    !! Setter for the "SumGDD" global variable.
    real(dp), intent(in) :: SumGDD_in

    SumGDD = SumGDD_in
end subroutine SetSumGDD


real(dp) function GetTimeSenescence()
    !! Getter for the "TimeSenescence" global variable.

    GetTimeSenescence = TimeSenescence
end function GetTimeSenescence


subroutine SetTimeSenescence(TimeSenescence_in)
    !! Setter for the "TimeSenescence" global variable.
    real(dp), intent(in) :: TimeSenescence_in

    TimeSenescence = TimeSenescence_in
end subroutine SetTimeSenescence


real(dp) function GetSumKcTop()
    !! Getter for the "SumKcTop" global variable.

    GetSumKcTop = SumKcTop
end function GetSumKcTop


subroutine SetSumKcTop(SumKcTop_in)
    !! Setter for the "SumKcTop" global variable.
    real(dp), intent(in) :: SumKcTop_in

    SumKcTop = SumKcTop_in
end subroutine SetSumKcTop


real(dp) function GetSumKcTopStress()
    !! Getter for the "SumKcTopStress" global variable.

    GetSumKcTopStress = SumKcTopStress
end function GetSumKcTopStress


subroutine SetSumKcTopStress(SumKcTopStress_in)
    !! Setter for the "SumKcTopStress" global variable.
    real(dp), intent(in) :: SumKcTopStress_in

    SumKcTopStress = SumKcTopStress_in
end subroutine SetSumKcTopStress


real(dp) function GetSumKci()
    !! Getter for the "SumKci" global variable.

    GetSumKci = SumKci
end function GetSumKci


subroutine SetSumKci(SumKci_in)
    !! Setter for the "SumKci" global variable.
    real(dp), intent(in) :: SumKci_in

    SumKci = SumKci_in
end subroutine SetSumKci


real(dp) function GetCCxCropWeedsNoSFstress()
    !! Getter for the "CCxCropWeedsNoSFstress" global variable.

    GetCCxCropWeedsNoSFstress = CCxCropWeedsNoSFstress
end function GetCCxCropWeedsNoSFstress


subroutine SetCCxCropWeedsNoSFstress(CCxCropWeedsNoSFstress_in)
    !! Setter for the "CCxCropWeedsNoSFstress" global variable.
    real(dp), intent(in) :: CCxCropWeedsNoSFstress_in

    CCxCropWeedsNoSFstress = CCxCropWeedsNoSFstress_in
end subroutine SetCCxCropWeedsNoSFstress


real(dp) function GetZiprev()
    !! Getter for the "Ziprev" global variable.

    GetZiprev = Ziprev
end function GetZiprev


subroutine SetZiprev(Ziprev_in)
    !! Setter for the "Ziprev" global variable.
    real(dp), intent(in) :: Ziprev_in

    Ziprev = Ziprev_in
end subroutine SetZiprev


real(dp) function GetSumGDDPrev()
    !! Getter for the "SumGDDPrev" global variable.

    GetSumGDDPrev = SumGDDPrev
end function GetSumGDDPrev


subroutine SetSumGDDPrev(SumGDDPrev_in)
    !! Setter for the "SumGDDPrev" global variable.
    real(dp), intent(in) :: SumGDDPrev_in

    SumGDDPrev = SumGDDPrev_in
end subroutine SetSumGDDPrev


real(dp) function GetCCoTotal()
    !! Getter for the "CCoTotal" global variable.

    GetCCoTotal = CCoTotal
end function GetCCoTotal


subroutine SetCCoTotal(CCoTotal_in)
    !! Setter for the "CCoTotal" global variable.
    real(dp), intent(in) :: CCoTotal_in

    CCoTotal = CCoTotal_in
end subroutine SetCCoTotal


real(dp) function GetCCxTotal()
    !! Getter for the "CCxTotal" global variable.

    GetCCxTotal = CCxTotal
end function GetCCxTotal


subroutine SetCCxTotal(CCxTotal_in)
    !! Setter for the "CCxTotal" global variable.
    real(dp), intent(in) :: CCxTotal_in

    CCxTotal = CCxTotal_in
end subroutine SetCCxTotal


real(dp) function GetCDCTotal()
    !! Getter for the "CDCTotal" global variable.

    GetCDCTotal = CDCTotal
end function GetCDCTotal


subroutine SetCDCTotal(CDCTotal_in)
    !! Setter for the "CDCTotal" global variable.
    real(dp), intent(in) :: CDCTotal_in

    CDCTotal = CDCTotal_in
end subroutine SetCDCTotal


real(dp) function GetGDDCDCTotal()
    !! Getter for the "GDDCDCTotal" global variable.

    GetGDDCDCTotal = GDDCDCTotal
end function GetGDDCDCTotal


subroutine SetGDDCDCTotal(GDDCDCTotal_in)
    !! Setter for the "GDDCDCTotal" global variable.
    real(dp), intent(in) :: GDDCDCTotal_in

    GDDCDCTotal = GDDCDCTotal_in
end subroutine SetGDDCDCTotal


real(dp) function GetWeedRCi()
    !! Getter for the "WeedRCi" global variable.

    GetWeedRCi = WeedRCi
end function GetWeedRCi


subroutine SetWeedRCi(WeedRCi_in)
    !! Setter for the "WeedRCi" global variable.
    real(dp), intent(in) :: WeedRCi_in

    WeedRCi = WeedRCi_in
end subroutine SetWeedRCi


real(dp) function GetCCiActualWeedInfested()
    !! Getter for the "CCiActualWeedInfested" global variable.

    GetCCiActualWeedInfested = CCiActualWeedInfested
end function GetCCiActualWeedInfested


subroutine SetCCiActualWeedInfested(CCiActualWeedInfested_in)
    !! Setter for the "CCiActualWeedInfested" global variable.
    real(dp), intent(in) :: CCiActualWeedInfested_in

    CCiActualWeedInfested = CCiActualWeedInfested_in
end subroutine SetCCiActualWeedInfested


real(dp) function GetfWeedNoS()
    !! Getter for the "fWeedNoS" global variable.

    GetfWeedNoS = fWeedNoS
end function GetfWeedNoS


subroutine SetfWeedNoS(fWeedNoS_in)
    !! Setter for the "fWeedNoS" global variable.
    real(dp), intent(in) :: fWeedNoS_in

    fWeedNoS = fWeedNoS_in
end subroutine SetfWeedNoS


real(dp) function GetBprevSum()
    !! Getter for the "BprevSum" global variable.

    GetBprevSum = BprevSum
end function GetBprevSum


subroutine SetBprevSum(BprevSum_in)
    !! Setter for the "BprevSum" global variable.
    real(dp), intent(in) :: BprevSum_in

    BprevSum = BprevSum_in
end subroutine SetBprevSum


real(dp) function GetYprevSum()
    !! Getter for the "YprevSum" global variable.

    GetYprevSum = YprevSum
end function GetYprevSum


subroutine SetYprevSum(YprevSum_in)
    !! Setter for the "YprevSum" global variable.
    real(dp), intent(in) :: YprevSum_in

    YprevSum = YprevSum_in
end subroutine SetYprevSum


real(dp) function GetSumGDDcuts()
    !! Getter for the "SumGDDcuts" global variable.

    GetSumGDDcuts = SumGDDcuts
end function GetSumGDDcuts


subroutine SetSumGDDcuts(SumGDDcuts_in)
    !! Setter for the "SumGDDcuts" global variable.
    real(dp), intent(in) :: SumGDDcuts_in

    SumGDDcuts = SumGDDcuts_in
end subroutine SetSumGDDcuts


real(dp) function GetHItimesBEF()
    !! Getter for the "HItimesBEF" global variable.

    GetHItimesBEF = HItimesBEF
end function GetHItimesBEF


subroutine SetHItimesBEF(HItimesBEF_in)
    !! Setter for the "HItimesBEF" global variable.
    real(dp), intent(in) :: HItimesBEF_in

    HItimesBEF = HItimesBEF_in
end subroutine SetHItimesBEF


real(dp) function GetScorAT1()
    !! Getter for the "ScorAT1" global variable.

    GetScorAT1 = ScorAT1
end function GetScorAT1


subroutine SetScorAT1(ScorAT1_in)
    !! Setter for the "ScorAT1" global variable.
    real(dp), intent(in) :: ScorAT1_in

    ScorAT1 = ScorAT1_in
end subroutine SetScorAT1


real(dp) function GetScorAT2()
    !! Getter for the "ScorAT2" global variable.

    GetScorAT2 = ScorAT2
end function GetScorAT2


subroutine SetScorAT2(ScorAT2_in)
    !! Setter for the "ScorAT2" global variable.
    real(dp), intent(in) :: ScorAT2_in

    ScorAT2 = ScorAT2_in
end subroutine SetScorAT2


real(dp) function GetHItimesAT1()
    !! Getter for the "HItimesAT1" global variable.

    GetHItimesAT1 = HItimesAT1
end function GetHItimesAT1


subroutine SetHItimesAT1(HItimesAT1_in)
    !! Setter for the "HItimesAT1" global variable.
    real(dp), intent(in) :: HItimesAT1_in

    HItimesAT1 = HItimesAT1_in
end subroutine SetHItimesAT1


real(dp) function GetHItimesAT2()
    !! Getter for the "HItimesAT2" global variable.

    GetHItimesAT2 = HItimesAT2
end function GetHItimesAT2


subroutine SetHItimesAT2(HItimesAT2_in)
    !! Setter for the "HItimesAT2" global variable.
    real(dp), intent(in) :: HItimesAT2_in

    HItimesAT2 = HItimesAT2_in
end subroutine SetHItimesAT2


real(dp) function GetHItimesAT()
    !! Getter for the "HItimesAT" global variable.

    GetHItimesAT = HItimesAT
end function GetHItimesAT


subroutine SetHItimesAT(HItimesAT_in)
    !! Setter for the "HItimesAT" global variable.
    real(dp), intent(in) :: HItimesAT_in

    HItimesAT = HItimesAT_in
end subroutine SetHItimesAT


real(dp) function GetalfaHI()
    !! Getter for the "alfaHI" global variable.

    GetalfaHI = alfaHI
end function GetalfaHI


subroutine SetalfaHI(alfaHI_in)
    !! Setter for the "alfaHI" global variable.
    real(dp), intent(in) :: alfaHI_in

    alfaHI = alfaHI_in
end subroutine SetalfaHI


real(dp) function GetalfaHIAdj()
    !! Getter for the "alfaHIAdj" global variable.

    GetalfaHIAdj = alfaHIAdj
end function GetalfaHIAdj


subroutine SetalfaHIAdj(alfaHIAdj_in)
    !! Setter for the "alfaHIAdj" global variable.
    real(dp), intent(in) :: alfaHIAdj_in

    alfaHIAdj = alfaHIAdj_in
end subroutine SetalfaHIAdj


real(dp) function GetPreviousSumETo()
    !! Getter for the "PreviousSumETo" global variable.

    GetPreviousSumETo = PreviousSumETo
end function GetPreviousSumETo


subroutine SetPreviousSumETo(PreviousSumETo_in)
    !! Setter for the "PreviousSumETo" global variable.
    real(dp), intent(in) :: PreviousSumETo_in

    PreviousSumETo = PreviousSumETo_in
end subroutine SetPreviousSumETo


real(dp) function GetPreviousSumGDD()
    !! Getter for the "PreviousSumGDD" global variable.

    GetPreviousSumGDD = PreviousSumGDD
end function GetPreviousSumGDD


subroutine SetPreviousSumGDD(PreviousSumGDD_in)
    !! Setter for the "PreviousSumGDD" global variable.
    real(dp), intent(in) :: PreviousSumGDD_in

    PreviousSumGDD = PreviousSumGDD_in
end subroutine SetPreviousSumGDD


real(dp) function GetPreviousBmob()
    !! Getter for the "PreviousBmob" global variable.

    GetPreviousBmob = PreviousBmob
end function GetPreviousBmob


subroutine SetPreviousBmob(PreviousBmob_in)
    !! Setter for the "PreviousBmob" global variable.
    real(dp), intent(in) :: PreviousBmob_in

    PreviousBmob = PreviousBmob_in
end subroutine SetPreviousBmob


real(dp) function GetPreviousBsto()
    !! Getter for the "PreviousBsto" global variable.

    GetPreviousBsto = PreviousBsto
end function GetPreviousBsto


subroutine SetPreviousBsto(PreviousBsto_in)
    !! Setter for the "PreviousBsto" global variable.
    real(dp), intent(in) :: PreviousBsto_in

    PreviousBsto = PreviousBsto_in
end subroutine SetPreviousBsto


integer(int32) function GetDayNr1Eval()
    !! Getter for the "DayNr1Eval" global variable.

    GetDayNr1Eval = DayNr1Eval
end function GetDayNr1Eval


subroutine SetDayNr1Eval(DayNr1Eval_in)
    !! Setter for the "DayNr1Eval" global variable.
    integer(int32), intent(in) :: DayNr1Eval_in

    DayNr1Eval = DayNr1Eval_in
end subroutine SetDayNr1Eval


integer(int32) function GetDayNrEval()
    !! Getter for the "DayNrEval" global variable.

    GetDayNrEval = DayNrEval
end function GetDayNrEval


subroutine SetDayNrEval(DayNrEval_in)
    !! Setter for the "DayNrEval" global variable.
    integer(int32), intent(in) :: DayNrEval_in

    DayNrEval = DayNrEval_in
end subroutine SetDayNrEval


integer(int32) function GetLineNrEval()
    !! Getter for the "LineNrEval" global variable.

    GetLineNrEval = int(LineNrEval, kind=int32)
end function GetLineNrEval


subroutine SetLineNrEval(LineNrEval_in)
    !! Setter for the "LineNrEval" global variable.
    integer(int32), intent(in) :: LineNrEval_in

    LineNrEval = int(LineNrEval_in, kind=int8)
end subroutine SetLineNrEval


real(dp) function GetZeval()
    !! Getter for the "Zeval" global variable.

    GetZeval = Zeval
end function GetZeval


subroutine SetZeval(Zeval_in)
    !! Setter for the "Zeval" global variable.
    real(dp), intent(in) :: Zeval_in

    Zeval = Zeval_in
end subroutine SetZeval


integer(int32) function GetNextSimFromDayNr()
    !! Getter for the "NextSimFromDayNr " global variable.

    GetNextSimFromDayNr = NextSimFromDayNr
end function GetNextSimFromDayNr


subroutine SetNextSimFromDayNr(NextSimFromDayNr_in)
    !! Setter for the "NextSimFromDayNr " global variable.
    integer(int32), intent(in) :: NextSimFromDayNr_in

    NextSimFromDayNr = NextSimFromDayNr_in
end subroutine SetNextSimFromDayNr


integer(int8) function GetStageCode()
    !! Getter for the "StageCode" global variable.

    GetStageCode = StageCode
end function GetStageCode


subroutine SetStageCode(StageCode_in)
    !! Setter for the "StageCode" global variable.
    integer(int8), intent(in) :: StageCode_in

    StageCode = StageCode_in
end subroutine SetStageCode


integer(int32) function GetPreviousDayNr()
    !! Getter for the "PreviousDayNr" global variable.

    GetPreviousDayNr = PreviousDayNr
end function GetPreviousDayNr


subroutine SetPreviousDayNr(PreviousDayNr_in)
    !! Setter for the "PreviousDayNr" global variable.
    integer(int32), intent(in) :: PreviousDayNr_in

    PreviousDayNr = PreviousDayNr_in
end subroutine SetPreviousDayNr


logical function GetNoYear()
    !! Getter for the NoYear global variable

    GetNoYear = NoYear
end function GetNoYear


subroutine SetNoYear(NoYear_in)
    !! Setter for the NoYear global variable
    logical, intent(in) :: NoYear_in

    NoYear = NoYear_in
end subroutine SetNoYear


!! END section global variables


subroutine AdjustForWatertable()

    real(dp) :: Ztot, Zi
    integer(int32) :: compi
    type(CompartmentIndividual) :: Compi_temp

    Ztot = 0.0_dp
    do compi = 1, GetNrCompartments()
        Ztot = Ztot + GetCompartment_Thickness(compi)
        Zi = Ztot - GetCompartment_Thickness(compi)/2.0_dp
        if (Zi >= (GetZiAqua()/100.0_dp)) then
            ! compartment at or below groundwater table
            call SetCompartment_Theta(compi, &
                GetSoilLayer_SAT(GetCompartment_Layer(compi))/100.0_dp)
            Compi_temp = GetCompartment_i(compi)
            call DetermineSaltContent(GetECiAqua(), Compi_temp)
            call SetCompartment_i(compi, Compi_temp)
        end if
    end do
end subroutine AdjustForWatertable


subroutine ResetPreviousSum(PreviousSum)
    type(rep_sum), intent(inout) :: PreviousSum

    PreviousSum%Epot = 0._dp
    PreviousSum%Tpot = 0._dp
    PreviousSum%Rain = 0._dp
    PreviousSum%Irrigation = 0._dp
    PreviousSum%Infiltrated = 0._dp
    PreviousSum%Runoff = 0._dp
    PreviousSum%Drain = 0._dp
    PreviousSum%Eact = 0._dp
    PreviousSum%Tact = 0._dp
    PreviousSum%TrW = 0._dp
    PreviousSum%ECropCycle = 0._dp
    PreviousSum%CRwater = 0._dp
    PreviousSum%Biomass = 0._dp
    PreviousSum%YieldPart = 0._dp
    PreviousSum%BiomassPot = 0._dp
    PreviousSum%BiomassUnlim = 0._dp
    PreviousSum%SaltIn = 0._dp
    PreviousSum%SaltOut = 0._dp
    PreviousSum%CRsalt = 0._dp
    call SetSumETo(0.0_dp)
    call SetSumGDD(0.0_dp)
    call SetPreviousSumETo(0.0_dp)
    call SetPreviousSumGDD(0.0_dp)
    call SetPreviousBmob(0.0_dp)
    call SetPreviousBsto(0.0_dp)
end subroutine ResetPreviousSum


subroutine CheckForPrint(TheProjectFile)
    character(len=*), intent(in) :: TheProjectFile

    integer(int32) :: DayN, MonthN, YearN, DayEndM
    real(dp) :: SaltIn, SaltOut, CRsalt, BiomassDay, BUnlimDay
    logical :: WriteNow

    call DetermineDate(GetDayNri(), DayN, MonthN, YearN)

    select case (GetOutputAggregate())
    case (1)
        ! 1: daily output
        BiomassDay = GetSumWaBal_Biomass() - GetPreviousSum_Biomass()
        BUnlimDay = GetSumWaBal_BiomassUnlim() - GetPreviousSum_BiomassUnlim()
        SaltIn = GetSumWaBal_SaltIn() - GetPreviousSum_SaltIn()
        SaltOut = GetSumWaBal_SaltOut() - GetPreviousSum_SaltOut()
        CRsalt = GetSumWaBal_CRsalt() - GetPreviousSum_CRsalt()
        call WriteTheResults(int(undef_int,kind=int8), DayN, MonthN, YearN, DayN, MonthN, &
                             YearN, GetRain(), GetETo(), GetGDDayi(), GetIrrigation(), &
                             GetInfiltrated(), GetRunoff(), GetDrain(), &
                             GetCRwater(), GetEact(), GetEpot(), GetTact(), &
                             GetTactWeedInfested(), GetTpot(), SaltIn, SaltOut, &
                             CRsalt, BiomassDay, BUnlimDay, GetBin(), GetBout(), &
                             TheProjectFile)
        call SetPreviousSum_Biomass(GetSumWaBal_Biomass())
        call SetPreviousSum_BiomassUnlim(GetSumWaBal_BiomassUnlim())
        call SetPreviousSum_SaltIn(GetSumWaBal_SaltIn())
        call SetPreviousSum_SaltOut(GetSumWaBal_SaltOut())
        call SetPreviousSum_CRsalt(GetSumWaBal_CRsalt())

    case (2,3)
        ! 2 or 3: 10-day or monthly output
        WriteNow = .false.
        DayEndM = DaysInMonth(MonthN)
        if (LeapYear(YearN) .and. (MonthN == 2)) then
            DayEndM = 29
        end if
        if (DayN == DayEndM) then
            WriteNow = .true.  ! 10-day and month
        end if
        if ((GetOutputAggregate() == 2) .and. ((DayN == 10) .or. (DayN == 20))) then
            WriteNow = .true. ! 10-day
        end if
        if (WriteNow) then
            call WriteIntermediatePeriod(TheProjectFile)
        end if
    end select
end subroutine CheckForPrint


subroutine GetGwtSet(DayNrIN, GwT)
    integer(int32), intent(in) :: DayNrIN
    type(rep_GwTable), intent(inout) :: GwT

    integer :: f0
    character(len=:), allocatable :: FileNameFull
    integer(int32) :: DayNr1Gwt, DNrini, rc
    integer(int32) :: i, dayi, monthi, yeari, Zini, yearACT
    real(dp) :: DayDouble, Zm, ECini
    character(len=255) :: StringREAD
    logical :: TheEnd

    ! FileNameFull
    if (GetGroundWaterFile() /= '(None)') then
        FileNameFull = GetGroundWaterFileFull()
    else
        FileNameFull = trim(GetPathNameProg())//'GroundWater.AqC'
    end if

    ! Get DayNr1Gwt
    open(newunit=f0, file=trim(FileNameFull), status='old', &
                 action='read', iostat=rc)
    read(f0, *, iostat=rc) ! description
    read(f0, *, iostat=rc) ! AquaCrop Version number
    read(f0, *, iostat=rc) ! Mode
    read(f0, *, iostat=rc) dayi
    read(f0, *, iostat=rc) monthi
    read(f0, *, iostat=rc) yeari
    call DetermineDayNr(dayi, monthi, yeari, DayNr1Gwt)

    ! Read first observation
    do i = 1, 3
        read(f0, *, iostat=rc)
    end do
    read(f0, '(a)', iostat=rc) StringREAD
    call SplitStringInThreeParams(StringREAD, DayDouble, Zm, GwT%EC2)
    GwT%DNr2 = DayNr1Gwt + roundc(DayDouble, mold=1_int32) - 1
    GwT%Z2 = roundc(Zm * 100, mold=1_int32)
    if (rc == iostat_end) then
        TheEnd = .true.
    else
        TheEnd = .false.
    end if

    ! Read next observations
    if (TheEnd) then
        ! only one observation
        GwT%DNr1 = GetSimulation_FromDayNr()
        GwT%Z1 = GwT%Z2
        GwT%EC1 = GwT%EC2
        GwT%DNr2 = GetSimulation_ToDayNr()
    else
        ! defined year
        if (DayNr1Gwt > 365) then
            if (DayNrIN < GwT%DNr2) then
                ! DayNrIN before 1st observation
                GwT%DNr1 = GetSimulation_FromDayNr()
                GwT%Z1 = GwT%Z2
                GwT%EC1 = GwT%EC2
            else
                ! DayNrIN after or at 1st observation
                loop1: do
                    GwT%DNr1 = GwT%DNr2
                    GwT%Z1 = GwT%Z2
                    GwT%EC1 = GwT%EC2
                    read(f0, '(a)', iostat=rc) StringREAD
                    call SplitStringInThreeParams(StringREAD, DayDouble, Zm, GwT%EC2)
                    GwT%DNr2 = DayNr1Gwt + roundc(DayDouble, mold=1_int32) - 1
                    GwT%Z2 = roundc(Zm * 100, mold=1_int32)
                    if (DayNrIN < GwT%DNr2) then
                        TheEnd = .true.
                    end if
                    if (TheEnd .or. (rc == iostat_end)) exit loop1
                end do loop1
                if (.not. TheEnd) then
                    ! DayNrIN after last observation
                    GwT%DNr1 = GwT%DNr2
                    GwT%Z1 = GwT%Z2
                    GwT%EC1 = GwT%EC2
                    GwT%DNr2 = GetSimulation_ToDayNr()
                end if
            end if
        end if ! defined year

        ! undefined year
        if (DayNr1Gwt <= 365) then
            call DetermineDate(DayNrIN, dayi, monthi, yearACT)
            if (yearACT /= 1901) then
                ! make 1st observation defined
                call DetermineDate(GwT%DNr2, dayi, monthi, yeari)
                call DetermineDayNr(dayi, monthi, yearACT, GwT%DNr2)
            end if
            if (DayNrIN < GwT%DNr2) then
                ! DayNrIN before 1st observation
                loop2: do
                    read(f0, '(a)', iostat=rc) StringREAD
                    call SplitStringInThreeParams(StringREAD, DayDouble, Zm, GwT%EC1)
                    GwT%DNr1 = DayNr1Gwt + roundc(DayDouble, mold=1_int32) - 1
                    call DetermineDate(GwT%DNr1, dayi, monthi, yeari)
                    call DetermineDayNr(dayi, monthi, yearACT, GwT%DNr1)
                    GwT%Z1 = roundc(Zm * 100, mold=1_int32)
                    if (rc == iostat_end) exit loop2
                end do loop2
                GwT%DNr1 = GwT%DNr1 - 365
            else
                ! save 1st observation
                DNrini = GwT%DNr2
                Zini = GwT%Z2
                ECini = GwT%EC2
                ! DayNrIN after or at 1st observation
                loop3: do
                    GwT%DNr1 = GwT%DNr2
                    GwT%Z1 = GwT%Z2
                    GwT%EC1 = GwT%EC2
                    read(f0, '(a)', iostat=rc) StringREAD
                    call SplitStringInThreeParams(StringREAD, DayDouble, Zm, GwT%EC2)
                    GwT%DNr2 = DayNr1Gwt + roundc(DayDouble, mold=1_int32) - 1
                    if (yearACT /= 1901) then
                        ! make observation defined
                        call DetermineDate(GwT%DNr2, dayi, monthi, yeari)
                        call DetermineDayNr(dayi, monthi, yearACT, GwT%DNr2)
                    end if
                    GwT%Z2 = roundc(Zm * 100, mold=1_int32)
                    if (DayNrIN < GwT%DNr2) then
                        TheEnd = .true.
                    end if
                    if (TheEnd .or. (rc == iostat_end)) exit loop3
                end do loop3
                if (.not. TheEnd) then
                    ! DayNrIN after last observation
                    GwT%DNr1 = GwT%DNr2
                    GwT%Z1 = GwT%Z2
                    GwT%EC1 = GwT%EC2
                    GwT%DNr2 = DNrini + 365
                    GwT%Z2 = Zini
                    GwT%EC2 = ECini
                end if
            end if
        end if ! undefined year
    end if ! more than 1 observation
    close(f0)
end subroutine GetGwtSet


subroutine GetNextHarvest()
    logical :: InfoLoaded
    integer(int32) :: DayNrXX
    integer(int32) :: FromDay_temp
    real(dp) :: IntervalInfo_temp, IntervalGDD_temp, MassInfo_temp
    character(len=:), allocatable :: TempString

    if (.not. GetManagement_Cuttings_Generate()) then
        TempString = fCuts_read()
        if (.not. fCuts_eof()) then
            read(TempString, *) FromDay_temp
            call SetCutInfoRecord1_FromDay(FromDay_temp)
            call SetCutInfoRecord1_NoMoreInfo(.false.)
            if (GetManagement_Cuttings_FirstDayNr() /= undef_int) then
                ! scroll to start growing cycle
                DayNrXX = GetManagement_Cuttings_FirstDayNr() + GetCutInfoRecord1_FromDay() -1
                 do while ((DayNrXX < GetCrop_Day1()) .and. (GetCutInfoRecord1_NoMoreInfo() .eqv. .false.))
                    TempString = fCuts_read()
                    if (.not. fCuts_eof()) then
                        read(TempString, *) FromDay_temp
                        call SetCutInfoRecord1_FromDay(FromDay_temp)
                        DayNrXX = GetManagement_Cuttings_FirstDayNr() + GetCutInfoRecord1_FromDay() -1
                    else
                        call SetCutInfoRecord1_NoMoreInfo(.true.)
                    end if
                end do
            end if
        else
            call SetCutInfoRecord1_NoMoreInfo(.true.)
        end if
    else
        if (GetNrCut() == 0) then
            if (GetManagement_Cuttings_Criterion() == TimeCuttings_IntDay) then
                TempString = fCuts_read()
                read(TempString, *) FromDay_temp, IntervalInfo_temp
                call SetCutInfoRecord1_FromDay(FromDay_temp)
                call SetCutInfoRecord1_IntervalInfo(roundc(IntervalInfo_temp, mold=1))
            elseif (GetManagement_Cuttings_Criterion() == TimeCuttings_IntGDD) then
                TempString = fCuts_read()
                read(TempString, *) FromDay_temp, IntervalGDD_temp
                call SetCutInfoRecord1_FromDay(FromDay_temp)
                call SetCutInfoRecord1_IntervalGDD(IntervalGDD_temp)
            elseif ((GetManagement_Cuttings_Criterion() == TimeCuttings_DryB) &
                    .or. (GetManagement_Cuttings_Criterion() == TimeCuttings_DryY) &
                    .or. (GetManagement_Cuttings_Criterion() == TimeCuttings_FreshY)) then
                TempString = fCuts_read()
                read(TempString, *) FromDay_temp, MassInfo_temp
                call SetCutInfoRecord1_FromDay(FromDay_temp)
                call SetCutInfoRecord1_MassInfo(MassInfo_temp)
            end if
            if (GetCutInfoRecord1_FromDay() < GetManagement_Cuttings_Day1()) then
                call SetCutInfoRecord1_FromDay(GetManagement_Cuttings_Day1())
            end if
            InfoLoaded = .false.
        end if
        loop2: do
            TempString = fCuts_read()
            if (.not. fCuts_eof()) then
                if (GetManagement_Cuttings_Criterion() == TimeCuttings_IntDay) then
                    read(TempString, *) FromDay_temp, IntervalInfo_temp
                    call SetCutInfoRecord2_FromDay(FromDay_temp)
                    call SetCutInfoRecord2_IntervalInfo(roundc(IntervalInfo_temp, mold=1))
                elseif (GetManagement_Cuttings_Criterion() == TimeCuttings_IntGDD) then
                    read(TempString, *) FromDay_temp, IntervalGDD_temp
                    call SetCutInfoRecord2_FromDay(FromDay_temp)
                    call SetCutInfoRecord2_IntervalGDD(IntervalGDD_temp)
                elseif ((GetManagement_Cuttings_Criterion() == TimeCuttings_DryB) &
                    .or. (GetManagement_Cuttings_Criterion() == TimeCuttings_DryY) &
                    .or. (GetManagement_Cuttings_Criterion() == TimeCuttings_FreshY)) then
                    read(TempString, *) FromDay_temp, MassInfo_temp
                    call SetCutInfoRecord2_FromDay(FromDay_temp)
                    call SetCutInfoRecord2_MassInfo(MassInfo_temp)
                end if
                if (GetCutInfoRecord2_FromDay() < GetManagement_Cuttings_Day1()) then
                    call SetCutInfoRecord2_FromDay(GetManagement_Cuttings_Day1())
                end if
                if (GetCutInfoRecord2_FromDay() <= GetCutInfoRecord1_FromDay()) then
                    ! CutInfoRecord2 becomes CutInfoRecord1
                    call SetCutInfoRecord1_FromDay(GetCutInfoRecord2_FromDay())
                    if (GetManagement_Cuttings_Criterion() == TimeCuttings_IntDay) then
                        call SetCutInfoRecord1_IntervalInfo(GetCutInfoRecord2_IntervalInfo())
                    elseif (GetManagement_Cuttings_Criterion() == TimeCuttings_IntGDD) then
                        call SetCutInfoRecord1_IntervalGDD(GetCutInfoRecord2_IntervalGDD())
                    elseif ((GetManagement_Cuttings_Criterion() == TimeCuttings_DryB) &
                            .or. (GetManagement_Cuttings_Criterion() == TimeCuttings_DryY) &
                            .or. (GetManagement_Cuttings_Criterion() == TimeCuttings_FreshY)) then
                        call SetCutInfoRecord1_MassInfo(GetCutInfoRecord2_MassInfo())
                    end if
                    call SetCutInfoRecord1_NoMoreInfo(.false.)
                else ! complete CutInfoRecord1
                    call SetCutInfoRecord1_ToDay(GetCutInfoRecord2_FromDay() - 1)
                    call SetCutInfoRecord1_NoMoreInfo(.false.)
                    if (GetManagement_Cuttings_NrDays() /= undef_int) then
                        if (GetCutInfoRecord1_ToDay() > (GetManagement_Cuttings_Day1() + GetManagement_Cuttings_NrDays() -1)) then
                            call SetCutInfoRecord1_ToDay(GetManagement_Cuttings_Day1() + GetManagement_Cuttings_NrDays() -1)
                            call SetCutInfoRecord1_NoMoreInfo(.true.)
                        end if
                    end if
                    InfoLoaded = .true.
                end if
            else ! Eof(fCuts)
                if (GetNrCut() > 0) then ! CutInfoRecord2 becomes CutInfoRecord1
                    call SetCutInfoRecord1_FromDay(GetCutInfoRecord2_FromDay())
                    if (GetManagement_Cuttings_Criterion() == TimeCuttings_IntDay) then
                        call SetCutInfoRecord1_IntervalInfo(GetCutInfoRecord2_IntervalInfo())
                    elseif (GetManagement_Cuttings_Criterion() == TimeCuttings_IntGDD) then
                        call SetCutInfoRecord1_IntervalGDD(GetCutInfoRecord2_IntervalGDD())
                    elseif ((GetManagement_Cuttings_Criterion() == TimeCuttings_DryB) &
                            .or. (GetManagement_Cuttings_Criterion() == TimeCuttings_DryY) &
                            .or. (GetManagement_Cuttings_Criterion() == TimeCuttings_FreshY)) then
                        call SetCutInfoRecord1_MassInfo(GetCutInfoRecord2_MassInfo())
                    end if
                end if
                call SetCutInfoRecord1_ToDay(GetCrop_DaysToHarvest())
                if (GetManagement_Cuttings_NrDays() /= undef_int) then
                    if (GetCutInfoRecord1_ToDay() > (GetManagement_Cuttings_Day1() + GetManagement_Cuttings_NrDays() -1)) then
                        call SetCutInfoRecord1_ToDay(GetManagement_Cuttings_Day1() + GetManagement_Cuttings_NrDays() -1)
                    end if
                end if
                call SetCutInfoRecord1_NoMoreInfo(.true.)
                InfoLoaded = .true.
            end if
            if (InfoLoaded .eqv. .true.) exit loop2
        end do loop2
    end if
end subroutine GetNextHarvest


subroutine GetSumGDDBeforeSimulation(SumGDDtillDay, SumGDDtillDayM1)
    real(dp), intent(inout) :: SumGDDtillDay
    real(dp), intent(inout) :: SumGDDtillDayM1

    character(len=:), allocatable :: totalname
    integer :: fTemp
    integer(int32) :: i
    character(len=255) :: StringREAD
    integer(int32) :: DayX
    real(dp) :: Tmin_temp, Tmax_temp
    type(rep_DayEventDbl), dimension(31) :: TmaxDataSet_temp, &
                                            TminDataSet_temp

    call SetSimulation_SumGDD(0._dp)
    if (GetTemperatureFile() /= '(None)') then
        totalname = GetTemperatureFilefull()

        if (FileExists(totalname)) then
            select case (GetTemperatureRecord_DataType())
            case (datatype_daily)
                open(newunit=fTemp, file=trim(totalname), status='old', &
                                                          action='read')
                read(fTemp, *) ! description
                read(fTemp, *) ! time step
                read(fTemp, *) ! day
                read(fTemp, *) ! month
                read(fTemp, *) ! year
                read(fTemp, *)
                read(fTemp, *)
                read(fTemp, *)
                ! days before first day of simulation (= DayNri)
                do i = GetTemperatureRecord_FromDayNr(), (DayNri - 1)
                    if (i < GetCrop_Day1()) then
                        read(fTemp, *)
                    else
                        read(fTemp, '(a)') StringREAD
                        Tmin_temp = GetTmin()
                        Tmax_temp = GetTmax()
                        call SplitStringInTwoParams(StringREAD, Tmin_temp, Tmax_temp)
                        call SetTmin(Tmin_temp)
                        call SetTmax(Tmax_temp)
                        call SetSimulation_SumGDD(GetSimulation_SumGDD() &
                                + DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), &
                                             GetTmin(), GetTmax(), &
                                             GetSimulParam_GDDMethod()))
                    end if
                end do
                close(fTemp)

            case (datatype_decadely)
                DayX = GetCrop_Day1()
                ! first day of cropping
                TminDataSet_temp = GetTminDataSet()
                TmaxDataSet_temp = GetTmaxDataSet()
                call GetDecadeTemperatureDataSet(DayX, TminDataSet_temp, &
                                                 TmaxDataSet_temp)
                call SetTminDataSet(TminDataSet_temp)
                call SetTmaxDataSet(TmaxDataSet_temp)
                i = 1
                do while (GetTminDataSet_DayNr(i) /= DayX)
                    i = i+1
                end do
                call SetTmin(GetTminDataSet_Param(i))
                call SetTmax(GetTmaxDataSet_Param(i))
                call SetSimulation_SumGDD(DegreesDay(GetCrop_Tbase(), &
                                GetCrop_Tupper(), GetTmin(), GetTmax(), &
                                GetSimulParam_GDDMethod()))
                ! next days
                do while (DayX < DayNri)
                    DayX = DayX + 1
                    if (DayX > GetTminDataSet_DayNr(31)) then
                        TminDataSet_temp = GetTminDataSet()
                        TmaxDataSet_temp = GetTmaxDataSet()
                        call GetDecadeTemperatureDataSet(DayX, &
                                TminDataSet_temp, TmaxDataSet_temp)
                        call SetTminDataSet(TminDataSet_temp)
                        call SetTmaxDataSet(TmaxDataSet_temp)
                        i = 0
                    end if
                    i = i+1
                    call SetTmin(GetTminDataSet_Param(i))
                    call SetTmax(GetTmaxDataSet_Param(i))
                    call SetSimulation_SumGDD(GetSimulation_SumGDD() &
                                + DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), &
                                             GetTmin(), GetTmax(), &
                                             GetSimulParam_GDDMethod()))
                end do
            case (datatype_monthly)
                DayX = GetCrop_Day1()
                ! first day of cropping
                TminDataSet_temp = GetTminDataSet()
                TmaxDataSet_temp = GetTmaxDataSet()
                call GetMonthlyTemperatureDataSet(DayX, TminDataSet_temp, &
                                                  TmaxDataSet_temp)
                call SetTminDataSet(TminDataSet_temp)
                call SetTmaxDataSet(TmaxDataSet_temp)
                i = 1
                do while (GetTminDataSet_DayNr(i) /= DayX)
                    i = i+1
                end do
                call SetTmin(GetTminDataSet_Param(i))
                call SetTmax(GetTmaxDataSet_Param(i))
                call SetSimulation_SumGDD(&
                        DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), &
                                   GetTmin(), GetTmax(), &
                                   GetSimulParam_GDDMethod()))
                ! next days
                do while (DayX < DayNri)
                    DayX = DayX + 1
                    if (DayX > GetTminDataSet_DayNr(31)) then
                        TminDataSet_temp = GetTminDataSet()
                        TmaxDataSet_temp = GetTmaxDataSet()
                        call GetMonthlyTemperatureDataSet(&
                                DayX, TminDataSet_temp, TmaxDataSet_temp)
                        call SetTminDataSet(TminDataSet_temp)
                        call SetTmaxDataSet(TmaxDataSet_temp)
                        i = 0
                    end if
                    i = i+1
                    call SetTmin(GetTminDataSet_Param(i))
                    call SetTmax(GetTmaxDataSet_Param(i))
                    call SetSimulation_SumGDD(GetSimulation_SumGDD() &
                            + DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), &
                                         GetTmin(), GetTmax(), &
                                         GetSimulParam_GDDMethod()))
                end do
            end select
        end if
    end if
    if (GetTemperatureFile() == '(None)') then
        call SetSimulation_SumGDD(DegreesDay(&
                                 GetCrop_Tbase(), GetCrop_Tupper(), &
                                 GetSimulParam_Tmin(), GetSimulParam_Tmax(), &
                                 GetSimulParam_GDDMethod()) &
                                 * (DayNri - GetCrop_Day1() + 1))
        if (GetSimulation_SumGDD() < 0._dp) then
            call SetSimulation_SumGDD(0._dp)
        end if
        SumGDDtillDay = GetSimulation_SumGDD()
        SumGDDtillDayM1 = DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), &
                                     GetSimulParam_Tmin(), GetSimulParam_Tmax(), &
                                     GetSimulParam_GDDMethod()) &
                          * (DayNri - GetCrop_Day1())
        if (SumGDDtillDayM1 < 0._dp) then
            SumGDDtillDayM1 = 0._dp
        end if
    else
        SumGDDtillDay = GetSimulation_SumGDD()
        SumGDDtillDayM1 = SumGDDtillDay &
                         - DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), &
                                      GetTmin(), GetTmax(), &
                                      GetSimulParam_GDDMethod())
    end if
end subroutine GetSumGDDBeforeSimulation


subroutine RelationshipsForFertilityAndSaltStress()

    real(dp) :: Coeffb0_temp
    real(dp) :: Coeffb1_temp
    real(dp) :: Coeffb2_temp
    real(dp) :: Coeffb0Salt_temp
    real(dp) :: Coeffb1Salt_temp
    real(dp) :: Coeffb2Salt_temp

    real(dp) :: X10, X20, X30, X40, X50, X60, X70, X80, X90
    integer(int8) :: BioTop, BioLow
    real(dp) :: StrTop, StrLow

    ! 1. Soil fertility
    call SetFracBiomassPotSF(1._dp)

    ! 1.a Soil fertility (Coeffb0,Coeffb1,Coeffb2 : Biomass-Soil Fertility stress)
    if (GetCrop_StressResponse_Calibrated()) then
        call StressBiomassRelationship(GetCrop_DaysToCCini(), GetCrop_GDDaysToCCini(), &
                                  GetCrop_DaysToGermination(), &
                                  GetCrop_DaysToFullCanopy(), &
                                  GetCrop_DaysToSenescence(), &
                                  GetCrop_DaysToHarvest(), &
                                  GetCrop_DaysToFlowering(), &
                                  GetCrop_LengthFlowering(), &
                                  GetCrop_GDDaysToGermination(), &
                                  GetCrop_GDDaysToFullCanopy(), &
                                  GetCrop_GDDaysToSenescence(), &
                                  GetCrop_GDDaysToHarvest(), &
                                  GetCrop_WPy(), GetCrop_HI(), &
                                  GetCrop_CCo(), GetCrop_CCx(), &
                                  GetCrop_CGC(), GetCrop_GDDCGC(), &
                                  GetCrop_CDC(), GetCrop_GDDCDC(), &
                                  GetCrop_KcTop(), GetCrop_KcDecline(), &
                                  real(GetCrop_CCEffectEvapLate(), kind= dp), &
                                  GetCrop_Tbase(), &
                                  GetCrop_Tupper(), GetSimulParam_Tmin(), &
                                  GetSimulParam_Tmax(), GetCrop_GDtranspLow(), &
                                  GetCrop_WP(), GetCrop_dHIdt(), GetCO2i(), &
                                  GetCrop_Day1(), GetCrop_DeterminancyLinked(), &
                                  GetCrop_StressResponse(),GetCrop_subkind(), &
                                  GetCrop_ModeCycle(), Coeffb0_temp, Coeffb1_temp, &
                                  Coeffb2_temp, X10, X20, X30, X40, X50, X60, X70)
        call SetCoeffb0(Coeffb0_temp)
        call SetCoeffb1(Coeffb1_temp)
        call SetCoeffb2(Coeffb2_temp)
    else
        call SetCoeffb0(real(undef_int, kind=dp))
        call SetCoeffb1(real(undef_int, kind=dp))
        call SetCoeffb2(real(undef_int, kind=dp))
    end if

    ! 1.b Soil fertility : FracBiomassPotSF
    if ((abs(GetManagement_FertilityStress()) > epsilon(0._dp)) .and. &
                                     GetCrop_StressResponse_Calibrated()) then
        BioLow = 100_int8
        StrLow = 0._dp
        loop: do
            BioTop = BioLow
            StrTop = StrLow
            BioLow = BioLow - 1_int8
            StrLow = GetCoeffb0() + GetCoeffb1()*BioLow + GetCoeffb2()*BioLow*BioLow
            if (((StrLow >= GetManagement_FertilityStress()) &
                         .or. (BioLow <= 0) .or. (StrLow >= 99.99_dp))) exit loop
        end do loop
        if (StrLow >= 99.99_dp) then
            StrLow = 100._dp
        end if
        if (abs(StrLow-StrTop) < 0.001_dp) then
            call SetFracBiomassPotSF(real(BioTop, kind=dp))
        else
            call SetFracBiomassPotSF(real(BioTop, kind=dp) - (GetManagement_FertilityStress() &
                                                    - StrTop)/(StrLow-StrTop))
        end if
    call SetFracBiomassPotSF(GetFracBiomassPotSF()/100._dp)
    end if

    ! 2. soil salinity (Coeffb0Salt,Coeffb1Salt,Coeffb2Salt : CCx/KsSto - Salt stress)
    if (GetSimulation_SalinityConsidered() .eqv. .true.) then
        call CCxSaltStressRelationship(GetCrop_DaysToCCini(), &
                                  GetCrop_GDDaysToCCini(), &
                                  GetCrop_DaysToGermination(), &
                                  GetCrop_DaysToFullCanopy(), &
                                  GetCrop_DaysToSenescence(), &
                                  GetCrop_DaysToHarvest(), &
                                  GetCrop_DaysToFlowering(), &
                                  GetCrop_LengthFlowering(), &
                                  GetCrop_GDDaysToFlowering(), &
                                  GetCrop_GDDLengthFlowering(), &
                                  GetCrop_GDDaysToGermination(), &
                                  GetCrop_GDDaysToFullCanopy(), &
                                  GetCrop_GDDaysToSenescence(), &
                                  GetCrop_GDDaysToHarvest(), &
                                  GetCrop_WPy(), GetCrop_HI(), &
                                  GetCrop_CCo(), GetCrop_CCx(), &
                                  GetCrop_CGC(), GetCrop_GDDCGC(), &
                                  GetCrop_CDC(), GetCrop_GDDCDC(), &
                                  GetCrop_KcTop(), GetCrop_KcDecline(), &
                                  real(GetCrop_CCEffectEvapLate(), kind=dp),  &
                                  GetCrop_Tbase(), GetCrop_Tupper(), &
                                  GetSimulParam_Tmin(), GetSimulParam_Tmax(), &
                                  GetCrop_GDtranspLow(), GetCrop_WP(), &
                                  GetCrop_dHIdt(), GetCO2i(), GetCrop_Day1(), &
                                  GetCrop_DeterminancyLinked(), &
                                  GetCrop_subkind(), GetCrop_ModeCycle(), &
                                  GetCrop_CCsaltDistortion(),Coeffb0Salt_temp, &
                                  Coeffb1Salt_temp, Coeffb2Salt_temp, X10, X20, X30, &
                                  X40, X50, X60, X70, X80, X90)
        call SetCoeffb0Salt(Coeffb0Salt_temp)
        call SetCoeffb1Salt(Coeffb1Salt_temp)
        call SetCoeffb2Salt(Coeffb2Salt_temp)
    else
        call SetCoeffb0Salt(real(undef_int, kind=dp))
        call SetCoeffb1Salt(real(undef_int, kind=dp))
        call SetCoeffb2Salt(real(undef_int, kind=dp))
    end if
end subroutine RelationshipsForFertilityAndSaltStress


! extra for output of daily results  -----------------------------

subroutine DetermineGrowthStage(Dayi, CCiPrev)
    integer(int32), intent(in) :: Dayi
    real(dp), intent(in) :: CCiPrev

    integer(int32) :: VirtualDay

    VirtualDay = Dayi - GetSimulation_DelayedDays() - GetCrop_Day1()
    if (VirtualDay < 0) then
        call SetStageCode(0_int8) ! before cropping period
    else
        if (VirtualDay < GetCrop_DaysToGermination()) then
            call SetStageCode(1_int8) ! sown --> emergence OR transplant recovering
        else
            call SetStageCode(2_int8) ! vegetative development
            if ((GetCrop_subkind() == subkind_Grain) .and. &
                (VirtualDay >= GetCrop_DaysToFlowering())) then
                if (VirtualDay < (GetCrop_DaysToFlowering() + &
                                  GetCrop_LengthFlowering())) then
                    call SetStageCode(3_int8) ! flowering
                else
                    call SetStageCode(4_int8) ! yield formation
                end if
            end if
            if ((GetCrop_subkind() == subkind_Tuber) .and. &
                (VirtualDay >= GetCrop_DaysToFlowering())) then
                call SetStageCode(4_int8) ! yield formation
            end if
            if ((VirtualDay > GetCrop_DaysToGermination()) .and.&
                (CCiPrev < epsilon(0._dp))) then
                call SetStageCode(int(undef_int, kind=int8))  ! no growth stage
            end if
            if (VirtualDay >= &
                (GetCrop_Length_i(1)+GetCrop_Length_i(2)+ &
                 GetCrop_Length_i(3)+GetCrop_Length_i(4))) then
                call SetStageCode(0_int8) ! after cropping period
            end if
        end if
    end if
end subroutine DetermineGrowthStage


subroutine WriteTitleDailyResults(TheProjectType, TheNrRun)
    integer(intenum), intent(in) :: TheProjectType
    integer(int8), intent(in) :: TheNrRun

    character(len=1025) :: Str1, Str2, tempstring
    real(dp) :: NodeD, Zprof
    integer(int32) :: Compi

    ! A. Run number
    call fDaily_write('')
    if (TheProjectType == typeproject_TypePRM) then
        write(Str1, '(i4)') TheNrRun
        call fDaily_write('   Run:'// trim(Str1))
    end if

    ! B. thickness of soil profile and root zone
    if ((GetOut1Wabal()) .or. (GetOut3Prof()) .or. (GetOut4Salt())) then
        Zprof = 0._dp
        do compi =1, GetNrCompartments()
            Zprof = Zprof + GetCompartment_Thickness(compi)
        end do
        write(Str1,'(f4.2)') Zprof
        if (roundc(GetSoil_RootMax()*1000._dp, mold=1) == &
                                roundc(GetCrop_RootMax()*1000._dp, mold=1)) then
            write(Str2, '(f4.2)') GetCrop_RootMax()
        else
            write(Str2,'(f4.2)') GetSoil_RootMax()
        end if
    end if

    ! C. 1st line title
    call fDaily_write(trim('   Day Month  Year   DAP Stage'), .false.)

    ! C1. Water balance
    if (GetOut1Wabal()) then
        if ((GetOut2Crop()) .or. (GetOut3Prof()) .or. (GetOut4Salt()) &
           .or. (GetOut5CompWC()) .or. (GetOut6CompEC()) .or. (GetOut7Clim())) then
            write(tempstring,'(4a)') '   WC(',trim(Str1),')   Rain     Irri   Surf'// &
                  '   Infilt   RO    Drain       CR    Zgwt', &
               '       Ex       E     E/Ex     Trx       Tr  Tr/Trx    ETx      ET  ET/ETx'
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring,'(4a)') '   WC(', trim(Str1), ')   Rain     Irri   Surf'// &
                          '   Infilt   RO    Drain       CR    Zgwt', &
                          '       Ex       E     E/Ex     Trx       Tr  Tr/Trx'// &
                                                      '    ETx      ET  ET/ETx'
            call fDaily_write(trim(tempstring))
        end if
    end if
    ! C2. Crop development and yield
    if (GetOut2Crop()) then
        if ((GetOut3Prof()) .or. (GetOut4Salt()) .or. (GetOut5CompWC()) &
                              .or. (GetOut6CompEC()) .or. (GetOut7Clim())) then
            write(tempstring, '(2a)') '      GD       Z     StExp  StSto  StSen'// &
                          ' StSalt StWeed   CC      CCw     StTr  Kc(Tr)'// &
                          '     Trx       Tr      TrW  Tr/Trx   WP', &
                          '    Biomass     HI    Y(dry)  Y(fresh)  Brelative'// &
                          '    WPet      Bin     Bout'
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(2a)')  '      GD       Z     StExp  StSto'// &
                          '  StSen StSalt StWeed   CC      CCw     StTr'// &
                          '  Kc(Tr)     Trx       Tr      TrW  Tr/Trx   WP', &
                          '    Biomass     HI    Y(dry)  Y(fresh)  Brelative'// &
                          '    WPet      Bin     Bout'
            call fDaily_write(trim(tempstring))
        end if
    end if
    ! C3. Profile/Root zone - Soil water content
    if (GetOut3Prof()) then
        if ((GetOut4Salt()) .or. (GetOut5CompWC()) .or. (GetOut6CompEC()) &
                                                    .or. (GetOut7Clim())) then
            write(tempstring, '(5a)') '  WC(', trim(Str1), ') Wr(', trim(Str2), ')     Z'// &
                  '       Wr    Wr(SAT)    Wr(FC)   Wr(exp)   Wr(sto)   Wr(sen)   Wr(PWP)'
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(5a)') '  WC(', trim(Str1), ') Wr(',trim(Str2), ')     Z'// &
                            '       Wr    Wr(SAT)    Wr(FC)   Wr(exp)   Wr(sto)'// &
                                '   Wr(sen)   Wr(PWP)'
            call fDaily_write(trim(tempstring))
        end if
    end if
    ! C4. Profile/Root zone - soil salinity
    if (GetOut4Salt()) then
        if ((GetOut5CompWC()) .or. (GetOut6CompEC()) .or. (GetOut7Clim())) then
            write(tempstring, '(3a)') '    SaltIn    SaltOut   SaltUp'// &
                            '   Salt(', trim(Str1), ')  SaltZ     Z       ECe'// &
                                '    ECsw   StSalt  Zgwt    ECgw'
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(3a)') '    SaltIn    SaltOut   SaltUp'// &
                            '   Salt(', trim(Str1), ')  SaltZ     Z       ECe'// &
                                '    ECsw   StSalt  Zgwt    ECgw'
            call fDaily_write(trim(tempstring))
        end if
    end if
    ! C5. Compartments - Soil water content  --!removed tempstring
    if (GetOut5CompWC()) then
        call fDaily_write(trim('       WC01'), .false.)
        do Compi = 2, (GetNrCompartments()-1)
            write(Str1, '(i2)') Compi
            call fDaily_write('       WC'// trim(Str1), .false.)
        end do
        write(Str1,'(i2)') GetNrCompartments()
        if ((GetOut6CompEC()) .or. (GetOut7Clim())) then
            call fDaily_write('       WC'// trim(Str1), .false.)
        else
            call fDaily_write('       WC'// trim(Str1))
        end if
    end if
    ! C6. Compartmens - Electrical conductivity of the saturated soil-paste extract
    if (GetOut6CompEC()) then
        call fDaily_write(trim('      ECe01'), .false.)
        do Compi = 2, (GetNrCompartments()-1)
            write(Str1, '(i2)') Compi
            call fDaily_write('      ECe'// trim(Str1), .false.)
        end do
        write(Str1, '(i2)') GetNrCompartments()
        if (GetOut7Clim()) then
            call fDaily_write('      ECe'// trim(Str1), .false.)
        else
            call fDaily_write('      ECe'// trim(Str1))
        end if
    end if
    ! C7. Climate input parameters
    if (GetOut7Clim()) then
        call fDaily_write(trim('     Rain       ETo      Tmin      Tavg      Tmax      CO2'))
    end if

    call fDaily_write('                              ', .false.)
    ! D1. Water balance
    if (GetOut1Wabal()) then
        if ((GetOut2Crop()) .or. (GetOut3Prof()) .or. (GetOut4Salt()) &
            .or. (GetOut5CompWC()) .or. (GetOut6CompEC()) .or. (GetOut7Clim())) then
            write(tempstring, '(2a)') '        mm      mm       mm     mm' // &
                  '     mm     mm       mm       mm      m ', &
                    '       mm       mm     %        mm       mm    %' // &
                            '        mm      mm       %'
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(2a)') '        mm      mm       mm     mm' // &
                  '     mm     mm       mm       mm      m ', &
                    '       mm       mm     %        mm       mm    %' // &
                        '        mm      mm       %'
            call fDaily_write(trim(tempstring))
        end if
    end if
    ! D2. Crop development and yield
    if (GetOut2Crop()) then
        if ((GetOut3Prof()) .or. (GetOut4Salt()) .or. (GetOut5CompWC()) .or. &
                                    (GetOut6CompEC()) .or. (GetOut7Clim())) then
            write(tempstring, '(2a)') '  degC-day     m       %      %      %'// &
                           '      %      %      %       %       %       -'// &
                           '        mm       mm       mm    %     g/m2', &
                           '    ton/ha      %    ton/ha   ton/ha'// &
                           '       %       kg/m3   ton/ha   ton/ha'
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(2a)') '  degC-day     m       %      %      %'// &
                          '      %      %      %       %       %       -'// &
                          '        mm       mm       mm    %     g/m2', &
                          '    ton/ha      %    ton/ha   ton/ha'// &
                          '       %       kg/m3   ton/ha   ton/ha'
            call fDaily_write(trim(tempstring))
        end if
    end if
    ! D3. Profile/Root zone - Soil water content
    if (GetOut3Prof()) then
        if ((GetOut4Salt()) .or. (GetOut5CompWC()) .or. (GetOut6CompEC()) &
                                                     .or. (GetOut7Clim())) then
            call fDaily_write(trim('      mm       mm       m       mm        mm'// &
                              '        mm        mm        mm        mm'// &
                              '         mm'), .false.)
        else
            call fDaily_write(trim('      mm       mm       m       mm        mm '// &
                              '       mm        mm        mm        mm        mm'))
        end if
    end if
    ! D4. Profile/Root zone - soil salinity
    if (GetOut4Salt()) then
        if ((GetOut5CompWC()) .or. (GetOut6CompEC()) .or. (GetOut7Clim())) then
            call fDaily_write(trim('    ton/ha    ton/ha    ton/ha    ton/ha'// &
                              '    ton/ha     m      dS/m    dS/m      %'// &
                              '     m      dS/m'), .false.)
        else
            call fDaily_write(trim('    ton/ha    ton/ha    ton/ha    ton/ha'// &
                              '    ton/ha     m      dS/m    dS/m      %'// &
                              '     m      dS/m'))
        end if
    end if
    ! D5. Compartments - Soil water content
    if (GetOut5CompWC()) then
        NodeD = GetCompartment_Thickness(1)/2._dp
        write(tempstring,'(f11.2)') NodeD
        call fDaily_write(trim(tempstring), .false.)
        do Compi = 2, (GetNrCompartments()-1)
            NodeD = NodeD + GetCompartment_Thickness(Compi-1)/2._dp &
                    + GetCompartment_Thickness(Compi)/2._dp
            write(tempstring,'(f11.2)') NodeD
            call fDaily_write(trim(tempstring), .false.)
        end do
        NodeD = NodeD + GetCompartment_Thickness(GetNrCompartments()-1)/2._dp &
                + GetCompartment_Thickness(GetNrCompartments())/2._dp
        if ((GetOut6CompEC()) .or. (GetOut7Clim())) then
            write(tempstring,'(f11.2)') NodeD
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring,'(f11.2)') NodeD
            call fDaily_write(trim(tempstring))
        end if
    end if
    ! D6. Compartmens - Electrical conductivity of the saturated soil-paste extract
    if (GetOut6CompEC()) then
        NodeD = GetCompartment_Thickness(1)/2._dp
        write(tempstring,'(f11.2)') NodeD
        call fDaily_write(trim(tempstring), .false.)
        do Compi = 2, (GetNrCompartments()-1)
            NodeD = NodeD + GetCompartment_Thickness(Compi-1)/2._dp &
                    + GetCompartment_Thickness(compi)/2._dp
            write(tempstring,'(f11.2)') NodeD
            call fDaily_write(trim(tempstring), .false.)
        end do
        NodeD = NodeD + GetCompartment_Thickness(GetNrCompartments()-1)/2._dp &
                + GetCompartment_Thickness(GetNrCompartments())/2._dp
        if (GetOut7Clim()) then
            write(tempstring,'(f11.2)') NodeD
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(f11.2)') NodeD
            call fDaily_write(trim(tempstring))
        end if
    end if
    ! D7. Climate input parameters
    if (GetOut7Clim()) then
        call fDaily_write(trim('       mm        mm     degC      degC      degC       ppm'))
    end if
end subroutine WriteTitleDailyResults


subroutine FinalizeRun2(NrRun, TheProjectType)
    integer(int8), intent(in) :: NrRun
    integer(intenum), intent(in) :: TheProjectType

    call CloseClimateFiles()
    call CloseIrrigationFile()
    call CloseManagementFile()

    if (GetPart2Eval() .and. (GetObservationsFile() /= '(None)')) then
        call CloseEvalDataPerformEvaluation(NrRun)
    end if


    contains


    subroutine CloseEvalDataPerformEvaluation(NrRun)
        integer(int8), intent(in) :: NrRun

        character(len=:), allocatable :: totalnameEvalStat
        character(len=1024) :: StrNr

        ! 1. Close Evaluation data file  and file with observations
        call fEval_close()
        if (GetLineNrEval() /= undef_int) then
            call fObs_close()
        end if

        ! 2. Specify File name Evaluation of simulation results - Statistics
        StrNr = ''
        if (GetSimulation_MultipleRun() .and. (GetSimulation_NrRuns() > 1)) then
            write(StrNr, '(i0)') NrRun
        end if

        select case (TheProjectType)
        case(typeproject_typepro)
            totalnameEvalStat = GetPathNameOutp() // GetOutputName() // 'PROevaluation.OUT'
        case(typeproject_typeprm)
            write(StrNr, '(i0)') NrRun
            totalnameEvalStat = GetPathNameOutp() // GetOutputName() // 'PRM' // trim(StrNr) // 'evaluation.OUT'
        end select

        ! 3. Create Evaluation statistics file
        call WriteAssessmentSimulation(trim(StrNr), totalnameEvalStat, &
                                       TheProjectType, &
                                       GetSimulation_FromDayNr(), &
                                       GetSimulation_ToDayNr())

        ! 4. Delete Evaluation data file
        call fEval_erase()
    end subroutine CloseEvalDataPerformEvaluation


    subroutine CloseClimateFiles()
        if (GetEToFile() /= '(None)') then
            call fEToSIM_close()
        end if
        if (GetRainFile() /= '(None)') then
            call fRainSIM_close()
        end if
        if (GetTemperatureFile() /= '(None)') then
            call fTempSIM_close()
        end if
    end subroutine CloseClimateFiles


    subroutine CloseIrrigationFile()
        if ((GetIrriMode() == IrriMode_Manual) .or. (GetIrriMode() == IrriMode_Generate)) then
            call fIrri_close()
        end if
    end subroutine CloseIrrigationFile


    subroutine CloseManagementFile()
        if (GetManagement_Cuttings_Considered()) then
            call fCuts_close()
        end if
    end subroutine CloseManagementFile
end subroutine FinalizeRun2


subroutine OpenIrrigationFile()

    character(len=:), allocatable :: totalname
    character(len=255) :: StringREAD
    integer(int32) :: i, DNr
    real(dp) :: Ir1, Ir2
    real(dp) :: VersionNr
    integer(int32) :: FromDay_temp, TimeInfo_temp, DepthInfo_temp
    real(dp) :: IrriECw_temp
    character(len=1025) :: TempString

    if ((GetIrriMode() == IrriMode_Manual) &
        .or. (GetIrriMode() == IrriMode_Generate)) then
        if (GetIrriFile() /= '(None)') then
            totalname = GetIrriFileFull()
        else
            totalname = GetPathNameProg() // 'IrriSchedule.AqC'
        end if
        call fIrri_open(totalname, 'r')
        TempString = fIrri_read() ! description
        TempString = fIrri_read() ! AquaCrop version
        read(TempString, *) VersionNr

        if (roundc(VersionNr*10, mold=1) < 32) then
            call SetGlobalIrriECw(.true.)
        else
            call SetGlobalIrriECw(.false.)
        end if
        do i = 1, 6
            TempString = fIrri_read()  ! irrigation info (already loaded)
        end do
        select case (GetIrriMode())
        case (IrriMode_Manual)
            if (GetIrriFirstDayNr() == undef_int) then
                DNr = GetDayNri() - GetCrop_Day1() + 1
            else
                DNr = GetDayNri() - GetIrriFirstDayNr() + 1
            end if
            loop: do
                StringREAD = fIrri_read()
                if (fIrri_eof()) then
                    call SetIrriInfoRecord1_NoMoreInfo(.true.)
                else
                    call SetIrriInfoRecord1_NoMoreInfo(.false.)
                    if (GetGlobalIrriECw()) then
                        call SplitStringInTwoParams(StringREAD, Ir1, Ir2)
                    else
                        IrriECw_temp = GetSimulation_IrriECw()
                        call SplitStringInThreeParams(StringREAD, Ir1, Ir2, &
                                                      IrriECw_temp)
                        call SetSimulation_IrriECw(IrriECw_temp)
                    end if
                    call SetIrriInfoRecord1_TimeInfo(roundc(Ir1, mold=1))
                    call SetIrriInfoRecord1_DepthInfo(roundc(Ir2, mold=1))
                end if
                if ((GetIrriInfoRecord1_NoMoreInfo()) &
                    .or. (GetIrriInfoRecord1_TimeInfo() >= DNr)) exit loop
            end do loop
        case(IrriMode_Generate)
            do i = 1, 2
                TempString = fIrri_read()
                ! time and depth criterion (already loaded)
            end do
            call SetIrriInfoRecord1_NoMoreInfo(.false.)
            if (roundc(VersionNr*10, mold=1) < 32) then
                TempString = fIrri_read()
                read(TempString, *) FromDay_temp, TimeInfo_temp, &
                                    DepthInfo_temp
                call SetIrriInfoRecord1_FromDay(FromDay_temp)
                call SetIrriInfoRecord1_TimeInfo(TimeInfo_temp)
                call SetIrriInfoRecord1_DepthInfo(DepthInfo_temp)
            else
                TempString = fIrri_read()
                read(TempString, *) FromDay_temp, TimeInfo_temp, &
                                    DepthInfo_temp, IrriECw_temp
                call SetIrriInfoRecord1_FromDay(FromDay_temp)
                call SetIrriInfoRecord1_TimeInfo(TimeInfo_temp)
                call SetIrriInfoRecord1_DepthInfo(DepthInfo_temp)
                call SetSimulation_IrriECw(IrriECw_temp)
            end if

            TempString = fIrri_read()
            if (fIrri_eof()) then
                call SetIrriInfoRecord1_ToDay(GetCrop_DayN() &
                                              - GetCrop_Day1() + 1)
            else
                call SetIrriInfoRecord2_NoMoreInfo(.false.)
                if (GetGlobalIrriECw()) then
                    read(TempString, *) FromDay_temp, TimeInfo_temp, &
                                        DepthInfo_temp
                    call SetIrriInfoRecord2_FromDay(FromDay_temp)
                    call SetIrriInfoRecord2_TimeInfo(TimeInfo_temp)
                    call SetIrriInfoRecord2_DepthInfo(DepthInfo_temp)
                else
                    read(TempString, *) FromDay_temp, TimeInfo_temp, &
                                        DepthInfo_temp, IrriEcw_temp
                    call SetIrriInfoRecord2_FromDay(FromDay_temp)
                    call SetIrriInfoRecord2_TimeInfo(TimeInfo_temp)
                    call SetIrriInfoRecord2_DepthInfo(DepthInfo_temp)
                    call SetSimulation_IrriECw(IrriECw_temp)
                end if
                call SetIrriInfoRecord1_ToDay(GetIrriInfoRecord2_FromDay() - 1)
            end if
        end select
    end if
end subroutine OpenIrrigationFile



subroutine WriteTitlePart1MultResults(TheProjectType, TheNrRun)
    integer(intEnum), intent(in) :: TheProjectType
    integer(int8), intent(in) :: TheNrRun

    integer(int32) :: Dayi, Monthi, Yeari
    real(dp) :: Nr
    character(len=1024) :: tempstring

    ! A. Run number
    call fHarvest_write('')
    if (TheProjectType == typeproject_TypePRM) then
        write(tempstring, '(i4)') TheNrRun
        call fHarvest_write('   Run:' // trim(tempstring))
    end if

    ! B. Title
    call fHarvest_write('    Nr   Day  Month Year   DAP Interval  Biomass    ' // &
        'Sum(B)   Dry-Yield  Sum(Y) Fresh-Yield  Sum(Y)')
    call fHarvest_write('                                 days     ton/ha    ' // &
        'ton/ha    ton/ha    ton/ha    ton/ha    ton/ha')

    ! C. start crop cycle
    call DetermineDate(GetCrop_Day1(), Dayi, Monthi, Yeari)
    call SetNoYear(Yeari == 1901)
    if (GetNoYear()) then
        if (Dayi == 0) then
            Dayi = 1
        end if
        Yeari = 9999
    end if
    Nr = 0._dp
    write(tempstring, '(i6, 3i6, f34.3, 2f20.3)') 0, Dayi, Monthi, Yeari, &
                                                  Nr, Nr, Nr
    call fHarvest_write(trim(tempstring))
end subroutine WriteTitlePart1MultResults


subroutine WriteTheResults(ANumber, Day1, Month1, Year1, DayN, MonthN, &
                           YearN, RPer, EToPer, GDDPer, IrriPer, InfiltPer, &
                           ROPer, DrainPer, CRwPer, EPer, ExPer, TrPer, TrWPer, &
                           TrxPer, SalInPer, SalOutPer, &
                           SalCRPer, BiomassPer, BUnlimPer, BmobPer, BstoPer, &
                           TheProjectFile)
    integer(int8), intent(in) :: ANumber
    integer(int32), intent(in) :: Day1
    integer(int32), intent(in) :: Month1
    integer(int32), intent(in) :: Year1
    integer(int32), intent(in) :: DayN
    integer(int32), intent(in) :: MonthN
    integer(int32), intent(in) :: YearN
    real(dp), intent(in) :: RPer
    real(dp), intent(in) :: EToPer
    real(dp), intent(in) :: GDDPer
    real(dp), intent(in) :: IrriPer
    real(dp), intent(in) :: InfiltPer
    real(dp), intent(in) :: ROPer
    real(dp), intent(in) :: DrainPer
    real(dp), intent(in) :: CRwPer
    real(dp), intent(in) :: EPer
    real(dp), intent(in) :: ExPer
    real(dp), intent(in) :: TrPer
    real(dp), intent(in) :: TrWPer
    real(dp), intent(in) :: TrxPer
    real(dp), intent(in) :: SalInPer
    real(dp), intent(in) :: SalOutPer
    real(dp), intent(in) :: SalCRPer
    real(dp), intent(in) :: BiomassPer
    real(dp), intent(in) :: BUnlimPer
    real(dp), intent(in) :: BmobPer
    real(dp), intent(in) :: BstoPer
    character(len=*), intent(in) :: TheProjectFile

    integer(int32) :: BrSF, RatioE, RatioT
    integer(int32) :: Year1_loc, YearN_loc
    real(dp) :: WPy, HI, tempreal
    character(len=1025) :: TempString

    ! copy intent(in) variables to locals
    Year1_loc = Year1
    YearN_loc = YearN

    ! start
    if (GetNoYear()) then
        Year1_loc = 9999
        YearN_loc = 9999
    end if
    if (ANumber == int(undef_int, int8)) then ! intermediate results
        select case (GetOutputAggregate())
        case(1)
            write(TempString, '(a, 3i9)') '      Day', Day1, Month1, Year1_loc
        case(2)
            write(TempString, '(a, 3i9)') '    10Day', Day1, Month1, Year1_loc
        case(3)
            write(TempString, '(a, 3i9)') '    Month', Day1, Month1, Year1_loc
        end select
        call fRun_write(trim(TempString), .false.)
    else
        write(TempString, '(i0)') ANumber
        TempString = 'Tot(' // trim(TempString) // ')'
        do while (len(trim(TempString)) < 9)
            TempString = ' ' // trim(TempString)
        end do
        call fRun_write(trim(TempString), .false.)
        write(TempString, '(3i9)') Day1, Month1, Year1_loc
        call fRun_write(trim(TempString), .false.)
    end if

    tempreal = roundc(GDDPer*10._dp, mold=1)
    ! Climatic conditions
    write(TempString, '(3f9.1, f9.2)') Rper, EToPer,&
          tempreal/10._dp, GetCO2i()
    call fRun_write(trim(TempString), .false.)
    ! Soil water parameters
    if (ExPer > 0._dp) then
        RatioE = roundc(100._dp*EPer/ExPer, mold=1)
    else
        RatioE = undef_int
    end if
    if (TrxPer > 0._dp) then
        RatioT = roundc(100._dp*TrPer/TrxPer, mold=1)
    else
        RatioT = undef_int
    end if

    write(TempString, '(6f9.1, i9, 2f9.1, i9)') IrriPer, InfiltPer, ROPer, &
                    DrainPer, CRwPer, EPer, RatioE, TrPer, TrWPer, RatioT
    call fRun_write(trim(TempString), .false.)

    ! Soil Salinity
    write(TempString, '(4f10.3)') SalInPer, SalOutPer, SalCRPer, &
                                  GetTotalSaltContent_EndDay()
    call fRun_write(trim(TempString), .false.)

    ! Seasonal stress
    write(TempString, '(7i9)') GetStressTot_NrD(), &
                               roundc(GetStressTot_Salt(), mold=1), &
                               GetManagement_FertilityStress(), &
                               roundc(GetStressTot_Weed(), mold=1), &
                               roundc(GetStressTot_Temp(), mold=1), &
                               roundc(GetStressTot_Exp(), mold=1), &
                               roundc(GetStressTot_Sto(), mold=1)
    call fRun_write(trim(TempString), .false.)

    ! Biomass production
    if ((BiomassPer > 0._dp) .and. (BUnlimPer > 0._dp)) then
        BrSF = roundc(100._dp*BiomassPer/BUnlimPer, mold=1)
        if (BrSF > 100) then
            BrSF = 100
        end if
    else
        BrSF = undef_int
    end if

    write(TempString, '(f10.3, i9)') BiomassPer, BrSF
    call fRun_write(trim(TempString), .false.)

    ! Crop yield
    ! Harvest Index
    if ((GetSumWaBal_Biomass() > epsilon(0._dp)) &
        .and. (GetSumWaBal_YieldPart() > epsilon(0._dp))) then
        HI = 100._dp*(GetSumWaBal_YieldPart())/(GetSumWaBal_Biomass())
    else
        if (GetSumWaBal_Biomass() > epsilon(0._dp)) then
            HI = 0._dp
        else
            HI = undef_double
        end if
    end if

    if (ANumber /= int(undef_int, int8)) then ! end of simulation run
        ! Water Use Efficiency yield
        if (((GetSumWaBal_Tact() > 0._dp) .or. (GetSumWaBal_ECropCycle() > 0._dp)) &
            .and. (GetSumWaBal_YieldPart() > 0._dp)) then
            WPy = (GetSumWaBal_YieldPart()*1000._dp) &
                    /((GetSumWaBal_Tact()+GetSumWaBal_ECropCycle())*10._dp)
        else
            WPy = 0.0_dp
        end if

        ! Fresh yield
        if ((GetCrop_DryMatter() == int(undef_int, int8)) &
            .or. (GetCrop_DryMatter() < epsilon(0._dp))) then
            write(TempString, '(f9.1, 2f9.3, f9.2)') HI, GetSumWaBal_YieldPart(), &
                                                    undef_double, WPy
        else
            write(TempString, '(f9.1, 2f9.3, f9.2)') HI, GetSumWaBal_YieldPart(), &
                (GetSumWaBal_YieldPart()/(GetCrop_DryMatter()/100._dp)), WPy
        end if
        call fRun_write(trim(TempString), .false.)

        ! Transfer of assimilates
        write(TempString, '(2f9.3)') GetTransfer_Bmobilized(), &
                                     GetSimulation_Storage_Btotal()
        call fRun_write(trim(TempString), .false.)
    else
        write(TempString, '(f9.1, 4i9, 2f9.3)') HI, undef_int, undef_int, undef_int, &
                                          undef_int, BmobPer, BstoPer
        call fRun_write(trim(TempString), .false.)
    end if

    ! End
    write(TempString, '(3i9)') DayN, MonthN, YearN_loc
    call fRun_write(trim(TempString), .false.)

    ! Project
    call fRun_write('  ' // TheProjectFile)
end subroutine WriteTheResults


subroutine InitializeSimulationRunPart1()
    !! Part1 (before reading the climate) of the initialization of a run
    !! Initializes parameters and states

    integer(int32) :: DNr1, DNr2
    real(dp) :: fWeed, fi
    integer(int8) :: Cweed
    integer(int32) :: Day1, Month1, Year1
    integer(int8) :: FertStress
    type(rep_GwTable) :: GwTable_temp
    integer(int8) :: RedCGC_temp, RedCCX_temp, RCadj_temp
    type(rep_EffectStress) :: EffectStress_temp
    logical :: bool_temp
    integer(int32) :: Crop_DaysToFullCanopySF_temp
    logical :: WaterTableInProfile_temp

    ! 1. Adjustments at start
    ! 1.1 Adjust soil water and salt content if water table IN soil profile
    WaterTableInProfile_temp = GetWaterTableInProfile()
    call CheckForWaterTableInProfile((GetZiAqua()/100._dp), &
               GetCompartment(), WaterTableInProfile_temp)
    call SetWaterTableInProfile(WaterTableInProfile_temp)
    if (GetWaterTableInProfile()) then
        call AdjustForWatertable
    end if
    if (.not. GetSimulParam_ConstGwt()) then
        GwTable_temp = GetGwTable()
        call GetGwtSet(GetSimulation_FromDayNr(), GwTable_temp)
        call SetGwTable(GwTable_temp)
    end if

    ! 1.2 Check if FromDayNr simulation needs to be adjusted
    ! from previous run if Keep initial SWC
    if ((GetSWCIniFile() == 'KeepSWC') .and. &
        (GetNextSimFromDayNr() /= undef_int)) then
        ! assign the adjusted DayNr defined in previous run
        if (GetNextSimFromDayNr() <= GetCrop_Day1()) then
            call SetSimulation_FromDayNr(GetNextSimFromDayNr())
        end if
    end if
    call SetNextSimFromDayNr(undef_int)

    ! 2. initial settings for Crop
    call SetCrop_pActStom(GetCrop_pdef())
    call SetCrop_pSenAct(GetCrop_pSenescence())
    call SetCrop_pLeafAct(GetCrop_pLeafDefUL())
    call SetEvapoEntireSoilSurface(.true.)
    call SetSimulation_EvapLimitON(.false.)
    call SetSimulation_EvapWCsurf(0._dp)
    call SetSimulation_EvapZ(EvapZmin/100._dp)
    call SetSimulation_SumEToStress(0._dp)
    call SetCCxWitheredTpotNoS(0._dp) ! for calculation Maximum Biomass
                                      ! unlimited soil fertility
    call SetSimulation_DayAnaero(0_int8) ! days of anaerobic conditions in
                                    ! global root zone
    ! germination
    if ((GetCrop_Planting() == plant_Seed) .and. &
        (GetSimulation_FromDayNr() <= GetCrop_Day1())) then
        call SetSimulation_Germinate(.false.)
    else
        call SetSimulation_Germinate(.true.)
        ! since already germinated no protection required
        call SetSimulation_ProtectedSeedling(.false.)
    end if
    ! delayed germination
    call SetSimulation_DelayedDays(0)

    ! 3. create temperature file covering crop cycle
    if (GetTemperatureFile() /= '(None)') then
        if (GetSimulation_ToDayNr() < GetCrop_DayN()) then
            call TemperatureFileCoveringCropPeriod(GetCrop_Day1(), &
                       GetSimulation_TodayNr())
        else
            call TemperatureFileCoveringCropPeriod(GetCrop_Day1(), &
                       GetCrop_DayN())
        end if
    end if

    ! 4. CO2 concentration during cropping period
    DNr1 = GetSimulation_FromDayNr()
    if (GetCrop_Day1() > GetSimulation_FromDayNr()) then
        DNr1 = GetCrop_Day1()
    end if
    DNr2 = GetSimulation_ToDayNr()
    if (GetCrop_DayN() < GetSimulation_ToDayNr()) then
        DNr2 = GetCrop_DayN()
    end if
    call SetCO2i(CO2ForSimulationPeriod(DNr1, DNr2))

    ! 5. seasonals stress coefficients
    bool_temp = ((GetCrop_ECemin() /= undef_int) .and. &
                 (GetCrop_ECemax() /= undef_int)) .and. &
                 (GetCrop_ECemin() < GetCrop_ECemax())
    call SetSimulation_SalinityConsidered(bool_temp)
    if (GetIrriMode() == IrriMode_Inet) then
        call SetSimulation_SalinityConsidered(.false.)
    end if
    call SetStressTot_NrD(undef_int)
    call SetStressTot_Salt(0._dp)
    call SetStressTot_Temp(0._dp)
    call SetStressTot_Exp(0._dp)
    call SetStressTot_Sto(0._dp)
    call SetStressTot_Weed(0._dp)

    ! 6. Soil fertility stress
    ! Coefficients for soil fertility - biomass relationship
    ! AND for Soil salinity - CCx/KsSto relationship
    call RelationshipsForFertilityAndSaltStress()

    ! No soil fertility stress
    if (GetManagement_FertilityStress() <= 0) then
        call SetManagement_FertilityStress(0_int8)
    end if

    ! Reset soil fertility parameters to selected value in management
    EffectStress_temp = GetSimulation_EffectStress()
    call CropStressParametersSoilFertility(GetCrop_StressResponse(), &
            GetManagement_FertilityStress(), EffectStress_temp)
    call SetSimulation_EffectStress(EffectStress_temp)
    FertStress = GetManagement_FertilityStress()
    RedCGC_temp = GetSimulation_EffectStress_RedCGC()
    RedCCX_temp = GetSimulation_EffectStress_RedCCX()
    Crop_DaysToFullCanopySF_temp = GetCrop_DaysToFullCanopySF()
    call TimeToMaxCanopySF(GetCrop_CCo(), GetCrop_CGC(), GetCrop_CCx(), &
           GetCrop_DaysToGermination(), GetCrop_DaysToFullCanopy(), &
           GetCrop_DaysToSenescence(), GetCrop_DaysToFlowering(), &
           GetCrop_LengthFlowering(), GetCrop_DeterminancyLinked(), &
           Crop_DaysToFullCanopySF_temp, RedCGC_temp, RedCCX_temp, FertStress)
    call SetCrop_DaysToFullCanopySF(Crop_DaysToFullCanopySF_temp)
    call SetManagement_FertilityStress(FertStress)
    call SetSimulation_EffectStress_RedCGC(RedCGC_temp)
    call SetSimulation_EffectStress_RedCCX(RedCCX_temp)
    call SetPreviousStressLevel(int(GetManagement_FertilityStress(),kind=int32))
    call SetStressSFadjNEW(int(GetManagement_FertilityStress(),kind=int32))
    ! soil fertility and GDDays
    if (GetCrop_ModeCycle() == modeCycle_GDDays) then
        if (GetManagement_FertilityStress() /= 0_int8) then
            call SetCrop_GDDaysToFullCanopySF(GrowingDegreeDays(&
                  GetCrop_DaysToFullCanopySF(), GetCrop_Day1(), &
                  GetCrop_Tbase(), GetCrop_Tupper(), GetSimulParam_Tmin(),&
                  GetSimulParam_Tmax()))
        else
            call SetCrop_GDDaysToFullCanopySF(GetCrop_GDDaysToFullCanopy())
        end if
    end if

    ! Maximum sum Kc (for reduction WP in season if soil fertility stress)
    call SetSumKcTop(SeasonalSumOfKcPot(GetCrop_DaysToCCini(), &
            GetCrop_GDDaysToCCini(), GetCrop_DaysToGermination(), &
            GetCrop_DaysToFullCanopy(), GetCrop_DaysToSenescence(), &
            GetCrop_DaysToHarvest(), GetCrop_GDDaysToGermination(), &
            GetCrop_GDDaysToFullCanopy(), GetCrop_GDDaysToSenescence(), &
            GetCrop_GDDaysToHarvest(), GetCrop_CCo(), GetCrop_CCx(), &
            GetCrop_CGC(), GetCrop_GDDCGC(), GetCrop_CDC(), GetCrop_GDDCDC(), &
            GetCrop_KcTop(), GetCrop_KcDecline(), real(GetCrop_CCEffectEvapLate(),kind=dp), &
            GetCrop_Tbase(), GetCrop_Tupper(), GetSimulParam_Tmin(), &
            GetSimulParam_Tmax(), GetCrop_GDtranspLow(), GetCO2i(), &
            GetCrop_ModeCycle()))
    call SetSumKcTopStress( GetSumKcTop() * GetFracBiomassPotSF())
    call SetSumKci(0._dp)

    ! 7. weed infestation and self-thinning of herbaceous perennial forage crops
    ! CC expansion due to weed infestation and/or CC decrease as a result of
    ! self-thinning
    ! 7.1 initialize
    call SetSimulation_RCadj(GetManagement_WeedRC())
    Cweed = 0_int8
    if (GetCrop_subkind() == subkind_Forage) then
        fi = MultiplierCCxSelfThinning(int(GetSimulation_YearSeason(),kind=int32), &
              int(GetCrop_YearCCx(),kind=int32), GetCrop_CCxRoot())
    else
        fi = 1._dp
    end if
    ! 7.2 fweed
    if (GetManagement_WeedRC() > 0_int8) then
        call SetfWeedNoS(CCmultiplierWeed(GetManagement_WeedRC(), &
              GetCrop_CCx(), GetManagement_WeedShape()))
        call SetCCxCropWeedsNoSFstress( roundc(((100._dp*GetCrop_CCx() &
                  * GetfWeedNoS()) + 0.49),mold=1)/100._dp) ! reference for plot with weed
        if (GetManagement_FertilityStress() > 0_int8) then
            fWeed = 1._dp
            if ((fi > 0._dp) .and. (GetCrop_subkind() == subkind_Forage)) then
                Cweed = 1_int8
                if (fi > 0.005_dp) then
                    ! calculate the adjusted weed cover
                    call SetSimulation_RCadj(roundc(GetManagement_WeedRC() &
                         + Cweed*(1._dp-fi)*GetCrop_CCx()*&
                           (1._dp-GetSimulation_EffectStress_RedCCX()/100._dp)*&
                           GetManagement_WeedAdj()/100._dp, mold=1_int8))
                    if (GetSimulation_RCadj() < (100._dp * (1._dp- fi/(fi + (1._dp-fi)*&
                          (GetManagement_WeedAdj()/100._dp))))) then
                        call SetSimulation_RCadj(roundc(100._dp * (1._dp- fi/(fi + &
                              (1._dp-fi)*(GetManagement_WeedAdj()/100._dp))),mold=1_int8))
                    end if
                    if (GetSimulation_RCadj() > 100_int8) then
                        call SetSimulation_RCadj(98_int8)
                    end if
                else
                    call SetSimulation_RCadj(100_int8)
                end if
            end if
        else
            if (GetCrop_subkind() == subkind_Forage) then
                RCadj_temp = GetSimulation_RCadj()
                fweed = CCmultiplierWeedAdjusted(GetManagement_WeedRC(), &
                          GetCrop_CCx(), GetManagement_WeedShape(), &
                          fi, GetSimulation_YearSeason(), &
                          GetManagement_WeedAdj(), &
                          RCadj_temp)
                call SetSimulation_RCadj(RCadj_temp)
            else
                fWeed = GetfWeedNoS()
            end if
        end if
    else
        call SetfWeedNoS(1._dp)
        fWeed = 1._dp
        call SetCCxCropWeedsNoSFstress(GetCrop_CCx())
    end if
    ! 7.3 CC total due to weed infestation
    call SetCCxTotal( fWeed * GetCrop_CCx() * (fi+Cweed*(1._dp-fi)*&
           GetManagement_WeedAdj()/100._dp))
    call SetCDCTotal( GetCrop_CDC() * (fWeed*GetCrop_CCx()*&
           (fi+Cweed*(1._dp-fi)*GetManagement_WeedAdj()/100._dp) + 2.29_dp)/ &
           (GetCrop_CCx()*(fi+Cweed*(1-fi)*GetManagement_WeedAdj()/100._dp) &
            + 2.29_dp))
    call SetGDDCDCTotal(GetCrop_GDDCDC() * (fWeed*GetCrop_CCx()*&
           (fi+Cweed*(1._dp-fi)*GetManagement_WeedAdj()/100._dp) + 2.29_dp)/ &
           (GetCrop_CCx()*(fi+Cweed*(1-fi)*GetManagement_WeedAdj()/100._dp) &
            + 2.29_dp))
    if (GetCrop_subkind() == subkind_Forage) then
        fi = MultiplierCCoSelfThinning(int(GetSimulation_YearSeason(),kind=int32), &
               int(GetCrop_YearCCx(),kind=int32), GetCrop_CCxRoot())
    else
        fi = 1._dp
    end if
    call SetCCoTotal(fWeed * GetCrop_CCo() * (fi+Cweed*(1._dp-fi)*&
            GetManagement_WeedAdj()/100._dp))

    ! 8. prepare output files
    ! Not applicable

    ! 9. first day
    call SetStartMode(.true.)
    bool_temp = (.not. GetSimulation_ResetIniSWC())
    call SetPreDay(bool_temp)
    call SetDayNri(GetSimulation_FromDayNr())
    call DetermineDate(GetSimulation_FromDayNr(), Day1, Month1, Year1) ! start simulation run
    call SetNoYear((Year1 == 1901));  ! for output file
end subroutine InitializeSimulationRunPart1


subroutine InitializeClimate()
    !! Creates the Climate SIM files and reads climate of first day

    ! 10. Climate
    ! create climate files
    call CreateDailyClimFiles(GetSimulation_FromDayNr(), &
               GetSimulation_ToDayNr())
    ! climatic data for first day
    call OpenClimFilesAndGetDataFirstDay(GetDayNri())
end subroutine InitializeClimate


subroutine InitializeSimulationRunPart2()
    !! Part2 (after reading the climate) of the initialization of a run
    !! Initializes parameters and states

    integer(int32) :: tHImax, Dayi, DayCC
    real(dp) :: SumGDDforDayCC
    real(dp) :: CCiniMin, CCiniMax, RatDGDD
    real(dp) :: ECe_temp, ECsw_temp, ECswFC_temp, KsSalt_temp
    real(dp) :: SumGDD_temp, SumGDDFromDay1_temp

    ! Sum of GDD before start of simulation
    call SetSimulation_SumGDD(0._dp)
    call SetSimulation_SumGDDfromDay1(0._dp)
    if ((GetCrop_ModeCycle() == modeCycle_GDDays) .and. &
        (GetCrop_Day1() < GetDayNri())) then
        SumGDD_temp = GetSimulation_SumGDD()
        SumGDDfromDay1_temp = GetSimulation_SumGDDfromDay1()
        call GetSumGDDBeforeSimulation(SumGDD_temp, SumGDDfromDay1_temp)
         ! GDDays before start of simulation
        call SetSimulation_SumGDD(SumGDD_temp)
        call SetSimulation_SumGDDFromDay1(SumGDDFromDay1_temp)
    end if
    call SetSumGDDPrev( GetSimulation_SumGDDfromDay1())

    ! Sum of GDD at end of first day
    call SetGDDayi(DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), GetTmin(), &
                   GetTmax(), GetSimulParam_GDDMethod()))
    if (GetDayNri() >= GetCrop_Day1()) then
        if (GetDayNri() == GetCrop_Day1()) then
            call SetSimulation_SumGDD(GetSimulation_SumGDD() + GetGDDayi())
        end if
        call SetSimulation_SumGDDfromDay1(GetSimulation_SumGDDfromDay1() + &
                 GetGDDayi())
    end if
    ! Reset cummulative sums of ETo and GDD for Run output
    call SetSumETo(0._dp)
    call SetSumGDD(0._dp)

    ! 11. Irrigation
    call SetIrriInterval(1)
    call SetGlobalIrriECw(.true.) ! In Versions < 3.2 - Irrigation water
                                  ! quality is not yet recorded on file
    call OpenIrrigationFile()

    ! 12. Adjusted time when starting as regrowth
    if (GetCrop_DaysToCCini() /= 0) then
        ! regrowth
        call SetGDDTadj(undef_int)
        call SetGDDayFraction(real(undef_int, kind=dp))
        if (GetCrop_DaysToCCini() == undef_int) then
            call SetTadj(GetCrop_DaysToFullCanopy() - &
                         GetCrop_DaysToGermination())
        else
            call SetTadj(GetCrop_DaysToCCini())
        end if
        call SetDayFraction(real((GetCrop_DaysToSenescence()- &
              GetCrop_DaysToFullCanopy()),kind=dp)/real(GetTadj() + &
              GetCrop_DaysToGermination() + &
              (GetCrop_DaysToSenescence()-GetCrop_DaysToFullCanopy()),kind=dp))
        if (GetCrop_ModeCycle() == modeCycle_GDDays) then
            if (GetCrop_GDDaysToCCini() == undef_int) then
                call SetGDDTadj(GetCrop_GDDaysToFullCanopy() - &
                       GetCrop_GDDaysToGermination())
            else
                call SetGDDTadj(GetCrop_GDDaysToCCini())
            end if
            call SetGDDayFraction(real(GetCrop_GDDaysToSenescence() - &
                  GetCrop_GDDaysToFullCanopy(),kind=dp)/ &
                  real(GetGDDTadj() + GetCrop_GDDaysToGermination() + &
                   (GetCrop_GDDaysToSenescence() -&
                    GetCrop_GDDaysToFullCanopy()),kind=dp))
        end if
    else
        ! sowing or transplanting
        call SetTadj(0)
        call SetGDDTadj(0)
        call SetDayFraction(real(undef_int, kind=dp))
        call SetGDDayFraction(real(undef_int, kind=dp))
    end if

    ! 13. Initial canopy cover
    ! 13.1 default value
    ! 13.1a RatDGDD for simulation of CanopyCoverNoStressSF (CCi with decline)
    RatDGDD = 1._dp
    if (GetCrop_ModeCycle() == modeCycle_GDDays) then
        if (GetCrop_GDDaysToFullCanopySF() < GetCrop_GDDaysToSenescence()) then
            RatDGDD = (GetCrop_DaysToSenescence() - &
                       GetCrop_DaysToFullCanopySF()) / &
                      real(GetCrop_GDDaysToSenescence() -&
                           GetCrop_GDDaysToFullCanopySF(), kind=dp)
        end if
    end if
    ! 13.1b DayCC for initial canopy cover
    Dayi = GetDayNri() - GetCrop_Day1()
    if (GetCrop_DaysToCCini() == 0) then
        ! sowing or transplant
        DayCC = Dayi
        call SetDayFraction(real(undef_int, kind=dp))
    else
        ! adjust time (calendar days) for regrowth
        DayCC = Dayi + GetTadj() + GetCrop_DaysToGermination() ! adjusted time scale
        if (DayCC > GetCrop_DaysToHarvest()) then
            DayCC = GetCrop_DaysToHarvest() ! special case where L123 > L1234
        end if
        if (DayCC > GetCrop_DaysToFullCanopy()) then
            if (Dayi <= GetCrop_DaysToSenescence()) then
                DayCC = GetCrop_DaysToFullCanopy()  + &
                         roundc(GetDayFraction() * &
                         (Dayi+GetTadj()+GetCrop_DaysToGermination() -&
                         GetCrop_DaysToFullCanopy()),mold=1) ! slow down
            else
                DayCC = Dayi ! switch time scale
            end if
        end if
    end if
    ! 13.1c SumGDDayCC for initial canopy cover
    SumGDDforDayCC = undef_int
    if (GetCrop_ModeCycle() == modeCycle_GDDays) then
        if (GetCrop_GDDaysToCCini() == 0) then
            SumGDDforDayCC = GetSimulation_SumGDDfromDay1() - GetGDDayi()
        else
            ! adjust time (Growing Degree Days) for regrowth
            SumGDDforDayCC = GetSimulation_SumGDDfromDay1() - GetGDDayi() + &
                             GetGDDTadj() + GetCrop_GDDaysToGermination()
            if (SumGDDforDayCC > GetCrop_GDDaysToHarvest()) then
                SumGDDforDayCC = GetCrop_GDDaysToHarvest()
                ! special case where L123 > L1234
            end if
            if (SumGDDforDayCC > GetCrop_GDDaysToFullCanopy()) then
                if (GetSimulation_SumGDDfromDay1() <= &
                    GetCrop_GDDaysToFullCanopy()) then
                    SumGDDforDayCC = GetCrop_GDDaysToFullCanopy() + &
                      roundc(GetGDDayFraction() * &
                       (GetSimulation_SumGDDfromDay1()+GetGDDTadj()+ &
                       GetCrop_GDDaysToGermination()- &
                       GetCrop_GDDaysToFullCanopy()),mold=1)
                    ! slow down
                else
                    SumGDDforDayCC = GetSimulation_SumGDDfromDay1() - &
                                      GetGDDayi()
                    ! switch time scale
                end if
            end if
        end if
    end if
    ! 13.1d CCi at start of day (is CCi at end of previous day)
    if (GetDayNri() <= GetCrop_Day1()) then
        if (GetCrop_DaysToCCini() /= 0) then
            ! regrowth which starts on 1st day
            if (GetDayNri() == GetCrop_Day1()) then
                call SetCCiPrev(CCiNoWaterStressSF(DayCC, &
                   GetCrop_DaysToGermination(), &
                   GetCrop_DaysToFullCanopySF(), &
                   GetCrop_DaysToSenescence(), GetCrop_DaysToHarvest(), &
                   GetCrop_GDDaysToGermination(), &
                   GetCrop_GDDaysToFullCanopySF(), &
                   GetCrop_GDDaysToSenescence(), GetCrop_GDDaysToHarvest(), &
                   GetCCoTotal(), GetCCxTotal(), GetCrop_CGC(), &
                   GetCrop_GDDCGC(), GetCDCTotal(), GetGDDCDCTotal(), &
                   SumGDDforDayCC, RatDGDD, &
                   GetSimulation_EffectStress_RedCGC(), &
                   GetSimulation_EffectStress_RedCCX(), &
                   GetSimulation_EffectStress_CDecline(), GetCrop_ModeCycle()))
            else
                call SetCCiPrev(0._dp)
            end if
        else
            ! sowing or transplanting
            call SetCCiPrev(0._dp)
            if (GetDayNri() == (GetCrop_Day1()+GetCrop_DaysToGermination())) then
                call SetCCiPrev(GetCCoTotal())
            end if
        end if
    else
        if (GetDayNri() > GetCrop_DayN()) then
            call SetCCiPrev(0._dp)  ! after cropping period
        else
            call SetCCiPrev(CCiNoWaterStressSF(DayCC, &
                GetCrop_DaysToGermination(), &
                GetCrop_DaysToFullCanopySF(), GetCrop_DaysToSenescence(), &
                GetCrop_DaysToHarvest(), GetCrop_GDDaysToGermination(), &
                GetCrop_GDDaysToFullCanopySF(), GetCrop_GDDaysToSenescence(), &
                GetCrop_GDDaysToHarvest(), GetCCoTotal(), GetCCxTotal(), &
                GetCrop_CGC(), GetCrop_GDDCGC(), GetCDCTotal(),&
                GetGDDCDCTotal(), SumGDDforDayCC, RatDGDD, &
                GetSimulation_EffectStress_RedCGC(), &
                GetSimulation_EffectStress_RedCCX(), &
                GetSimulation_EffectStress_CDecline(), GetCrop_ModeCycle()))
        end if
    end if
    ! 13.2 specified CCini (%)
    if ((GetSimulation_CCini() > 0._dp) .and. &
        (roundc(10000._dp*GetCCiPrev(), mold=1) > 0) .and. &
        (roundc(GetSimulation_CCini(), mold=1) /= &
            roundc(100._dp*GetCCiPrev(),mold=1))) then
        ! 13.2a Minimum CC
        CCiniMin = 100._dp * (GetCrop_SizeSeedling()/10000._dp)*&
                    (GetCrop_PlantingDens()/10000._dp)
        if (CCiniMin - roundc(CCiniMin*100._dp, mold=1)/100._dp >= 0.00001) then
            CCiniMin = roundc(CCiniMin*100._dp + 1._dp, mold=1)/100._dp
        else
            CCiniMin = roundc(CCiniMin*100._dp, mold=1)/100._dp
        end if
        ! 13.2b Maximum CC
        CCiniMax = 100._dp * GetCCiPrev()
        CCiniMax = roundc(CCiniMax*100._dp, mold=1)/100._dp
        ! 13.2c accept specified CCini
        if ((GetSimulation_CCini() >= CCiniMin) .and. &
            (GetSimulation_CCini() <= CCiniMax)) then
            call SetCCiPrev(GetSimulation_CCini()/100._dp)
        end if
    end if
    ! 13.3
    call SetCrop_CCxAdjusted(GetCCxTotal())
    call SetCrop_CCoAdjusted(GetCCoTotal())
    call SetTimeSenescence(0._dp)
    call SetCrop_CCxWithered(0._dp)
    call SetNoMoreCrop(.false.)
    call SetCCiActual(GetCCiPrev())

    ! 14. Biomass and re-setting of GlobalZero
    if (roundc(1000._dp*GetSimulation_Bini(), mold=1) > 0) then
        ! overwrite settings in GlobalZero (in Global)
        call SetSumWaBal_Biomass(GetSimulation_Bini())
        call SetSumWaBal_BiomassPot(GetSimulation_Bini())
        call SetSumWaBal_BiomassUnlim(GetSimulation_Bini())
        call SetSumWaBal_BiomassTot(GetSimulation_Bini())
    end if

    ! 15. Transfer of assimilates
    if ((GetCrop_subkind() == subkind_Forage)&
        ! only valid for perennial herbaceous forage crops
        .and. (trim(GetCropFileFull()) == &
               trim(GetSimulation_Storage_CropString()))&
        ! only for the same crop
        .and. (GetSimulation_YearSeason() > 1) &
        ! mobilization not possible in season 1
        .and. (GetSimulation_YearSeason() == &
               (GetSimulation_Storage_Season() + 1))) then
            ! season next to season in which storage took place
            ! mobilization of assimilates
            if (GetSimulation_YearSeason() == 2) then
                call SetTransfer_ToMobilize(GetSimulation_Storage_Btotal() *& 
                    0.2 * GetCrop_Assimilates_Mobilized()/100._dp)
            else
                call SetTransfer_ToMobilize(GetSimulation_Storage_Btotal() *&
                    GetCrop_Assimilates_Mobilized()/100._dp)
            endif
            if (roundc(1000._dp * GetTransfer_ToMobilize(),&
                   mold=1) > 0) then ! minimum 1 kg
                 call SetTransfer_Mobilize(.true.)
            else
                 call SetTransfer_Mobilize(.false.)
            end if
    else
        call SetSimulation_Storage_CropString(GetCropFileFull())
        ! no mobilization of assimilates
        call SetTransfer_ToMobilize(0._dp)
        call SetTransfer_Mobilize(.false.)
    end if
    ! Storage is off and zero at start of season
    call SetSimulation_Storage_Season(GetSimulation_YearSeason())
    call SetSimulation_Storage_Btotal(0._dp)
    call SetTransfer_Store(.false.)
    ! Nothing yet mobilized at start of season
    call SetTransfer_Bmobilized(0._dp)

    ! 16. Initial rooting depth
    ! 16.1 default value
    if (GetDayNri() <= GetCrop_Day1()) then
        call SetZiprev(real(undef_int, kind=dp))
    else
        if (GetDayNri() > GetCrop_DayN()) then
            call SetZiprev(real(undef_int, kind=dp))
        else
            call SetZiprev( ActualRootingDepth(GetDayNri()-GetCrop_Day1(),&
                  GetCrop_DaysToGermination(),&
                  GetCrop_DaysToMaxRooting(),&
                  GetCrop_DaysToHarvest(),&
                  GetCrop_GDDaysToGermination(),&
                  GetCrop_GDDaysToMaxRooting(),&
                  GetSumGDDPrev(),&
                  GetCrop_RootMin(),&
                  GetCrop_RootMax(),&
                  GetCrop_RootShape(),&
                  GetCrop_ModeCycle()) )
        end if
    end if
    ! 16.2 specified or default Zrini (m)
    if ((GetSimulation_Zrini() > 0._dp) .and. &
        (GetZiprev() > 0._dp) .and. &
        (GetSimulation_Zrini() <= GetZiprev())) then
        if ((GetSimulation_Zrini() >= GetCrop_RootMin()) .and. &
            (GetSimulation_Zrini() <= GetCrop_RootMax())) then
            call SetZiprev( GetSimulation_Zrini())
        else
            if (GetSimulation_Zrini() < GetCrop_RootMin()) then
                call SetZiprev( GetCrop_RootMin())
            else
                call SetZiprev( GetCrop_RootMax())
            end if
        end if
        if ((roundc(GetSoil_RootMax()*1000._dp, mold=1) < &
             roundc(GetCrop_RootMax()*1000._dp, mold=1)) &
            .and. (GetZiprev() > GetSoil_RootMax())) then
            call SetZiprev(real(GetSoil_RootMax(), kind=dp))
        end if
        call SetRootingDepth(GetZiprev())
        ! NOT NEEDED since RootingDepth is calculated in the RUN by considering
        ! Ziprev
    else
        call SetRootingDepth(ActualRootingDepth(GetDayNri()-GetCrop_Day1()+1, &
              GetCrop_DaysToGermination(), &
              GetCrop_DaysToMaxRooting(),&
              GetCrop_DaysToHarvest(),&
              GetCrop_GDDaysToGermination(),&
              GetCrop_GDDaysToMaxRooting(),&
              GetSumGDDPrev(),&
              GetCrop_RootMin(),&
              GetCrop_RootMax(),&
              GetCrop_RootShape(),&
              GetCrop_ModeCycle()))
    end if

    ! 17. Multiple cuttings
    call SetNrCut(0)
    call SetSumInterval(0)
    call SetSumGDDcuts(0._dp)
    call SetBprevSum(0._dp)
    call SetYprevSum(0._dp)
    call SetCutInfoRecord1_IntervalInfo(0)
    call SetCutInfoRecord2_IntervalInfo(0)
    call SetCutInfoRecord1_MassInfo(0._dp)
    call SetCutInfoRecord2_MassInfo(0._dp)
    call SetDayLastCut(0)
    call SetCGCref( GetCrop_CGC())
    call SetGDDCGCref(GetCrop_GDDCGC())
    if (GetManagement_Cuttings_Considered()) then
        call OpenHarvestInfo()
    end if

    ! 18. Tab sheets

    ! 19. Labels, Plots and displays
    if (GetManagement_BundHeight() < 0.01_dp) then
        call SetSurfaceStorage(0._dp)
        call SetECStorage(0._dp)
    end if
    if (GetRootingDepth() > 0._dp) then
        ! salinity in root zone
        ECe_temp = GetRootZoneSalt_ECe()
        ECsw_temp = GetRootZoneSalt_ECsw()
        ECswFC_temp = GetRootZoneSalt_ECswFC()
        KsSalt_temp = GetRootZoneSalt_KsSalt()
        call DetermineRootZoneSaltContent(GetRootingDepth(), ECe_temp,&
               ECsw_temp, ECswFC_temp, KsSalt_temp)
        call SetRootZoneSalt_ECe(ECe_temp)
        call SetRootZoneSalt_ECsw(ECsw_temp)
        call SetRootZoneSalt_ECswFC(ECswFC_temp)
        call SetRootZoneSalt_KsSalt(KsSalt_temp)
        call SetStressTot_Salt(((GetStressTot_NrD() - 1._dp)*GetStressTot_Salt() + &
              100._dp*(1._dp-GetRootZoneSalt_KsSalt()))/real(GetStressTot_NrD(), kind=dp))
    end if
    ! Harvest Index
    call SetSimulation_HIfinal(GetCrop_HI())
    call SetHItimesBEF(real(undef_int, kind=dp))
    call SetHItimesAT1(1._dp)
    call SetHItimesAT2(1._dp)
    call SetHItimesAT(1._dp)
    call SetalfaHI(real(undef_int, kind=dp))
    call SetalfaHIAdj(0._dp)
    if (GetSimulation_FromDayNr() <= (GetSimulation_DelayedDays() + &
        GetCrop_Day1() + GetCrop_DaysToFlowering())) then
        ! not yet flowering
        call SetScorAT1(0._dp)
        call SetScorAT2(0._dp)
    else
        ! water stress affecting leaf expansion
        ! NOTE: time to reach end determinancy  is tHImax (i.e. flowering/2 or
        ! senescence)
        if (GetCrop_DeterminancyLinked()) then
            tHImax = roundc(GetCrop_LengthFlowering()/2._dp, mold=1)
        else
            tHImax = (GetCrop_DaysToSenescence() - GetCrop_DaysToFlowering())
        end if
        if ((GetSimulation_FromDayNr() <= (GetSimulation_DelayedDays() + &
            GetCrop_Day1() + GetCrop_DaysToFlowering() + tHImax)) & ! not yet end period
            .and. (tHImax > 0)) then
            ! not yet end determinancy
            call SetScorAT1(1._dp/tHImax)
            call SetScorAT1(GetScorAT1() * (GetSimulation_FromDayNr() - &
                  (GetSimulation_DelayedDays() + GetCrop_Day1() + &
                   GetCrop_DaysToFlowering())))
            if (GetScorAT1() > 1._dp) then
                call SetScorAT1(1._dp)
            end if
        else
            call SetScorAT1(1._dp)  ! after period of effect
        end if
        ! water stress affecting stomatal closure
        ! period of effect is yield formation
        if (GetCrop_dHIdt() > 99._dp) then
            tHImax = 0
        else
            tHImax = roundc(GetCrop_HI()/GetCrop_dHIdt(), mold=1)
        end if
        if ((GetSimulation_FromDayNr() <= (GetSimulation_DelayedDays() + &
             GetCrop_Day1() + GetCrop_DaysToFlowering() + tHImax)) & ! not yet end period
             .and. (tHImax > 0)) then
            ! not yet end yield formation
            call SetScorAT2(1._dp/real(tHImax, kind=dp))
            call SetScorAT2(GetScorAT2() * (GetSimulation_FromDayNr() - &
                  (GetSimulation_DelayedDays() + GetCrop_Day1() + &
                   GetCrop_DaysToFlowering())))
            if (GetScorAT2() > 1._dp) then
                call SetScorAT2(1._dp)
            end if
        else
            call SetScorAT2(1._dp)  ! after period of effect
        end if
    end if

    if (GetOutDaily()) then
        call DetermineGrowthStage(GetDayNri(), GetCCiPrev())
    end if

    ! 20. Settings for start
    call SetStartMode(.true.)
    call SetStressLeaf(real(undef_int, kind=dp))
    call SetStressSenescence(real(undef_int, kind=dp))
end subroutine InitializeSimulationRunPart2


subroutine CreateEvalData(NrRun)
    integer(int8), intent(in) :: NrRun

    integer(int32) :: dayi, monthi, yeari, integer_temp
    character(len=:), allocatable :: TempString
    character(len=1025) :: StrNr, tempstring2

    ! open input file with field data
    call fObs_open(GetObservationsFilefull(), 'r') ! Observations recorded in File
    TempString = fObs_read() ! description
    TempString = fObs_read() ! AquaCrop Version number
    TempString = fObs_read()
    read(TempString, *) Zeval !  depth of sampled soil profile
    TempString = fObs_read()
    read(TempString, *) dayi
    TempString = fObs_read()
    read(TempString, *)monthi
    TempString = fObs_read()
    read(TempString, *) yeari
    integer_temp = GetDayNr1Eval()
    call DetermineDayNr(dayi, monthi, yeari, integer_temp)
    call SetDayNr1Eval(integer_temp)
    TempString = fObs_read() ! title
    TempString = fObs_read() ! title
    TempString = fObs_read() ! title
    TempString = fObs_read() ! title
    call SetLineNrEval(undef_int)
    TempString = fObs_read()
    if (.not. fObs_eof()) then
        call SetLineNrEval(11)
        read(TempString, *) integer_temp
        call SetDayNrEval(integer_temp)
        call SetDayNrEval(GetDayNr1Eval() + GetDayNrEval() -1)
        do while ((GetDayNrEval() < GetSimulation_FromDayNr()) &
                    .and. (GetLineNrEval() /= undef_int))
            TempString = fObs_read()
            if (fObs_eof()) then
                call SetLineNrEval(undef_int)
            else
                call SetLineNrEval(GetLineNrEval() + 1)
                read(TempString, *) integer_temp
                call SetDayNrEval(integer_temp)
                call SetDayNrEval(GetDayNr1Eval() + GetDayNrEval() -1)
            end if
        end do
    end if
    if (GetLineNrEval() == undef_int) then
        call fObs_close()
    end if
    ! open file with simulation results, field data
    if (GetSimulation_MultipleRun() .and. (GetSimulation_NrRuns() > 1)) then
        write(StrNr, '(i0)') NrRun
    else
        StrNr = ''
    end if
    call SetfEval_filename(GetPathNameSimul() // 'EvalData' // trim(StrNr) // '.OUT')

    call fEval_open(GetfEval_filename(), 'w')
    write(tempstring2, '(a)') GetAquaCropDescriptionWithTimeStamp()
    call fEval_write(trim(tempstring2))
    call fEval_write('Evaluation of simulation results - Data')
    write(TempString, '(f5.2)') GetZeval()
    call fEval_write('                                             ' // &
    '                                        for soil depth: ' // trim(TempString) // ' m')
    call fEval_write('   Day Month  Year   DAP Stage   CCsim   CCobs   CCstd    Bsim      ' // &
    'Bobs      Bstd   SWCsim  SWCobs   SWstd')
    call fEval_write('                                   %       %       %     ton/ha    ' // &
    'ton/ha    ton/ha    mm       mm      mm')
end subroutine CreateEvalData


subroutine OpenOutputRun(TheProjectType)
    integer(intEnum), intent(in) :: TheProjectType

    character(len=:), allocatable :: totalname
    character(len=1025) :: tempstring

    select case (TheProjectType)
    case(typeproject_TypePRO)
        totalname = GetPathNameOutp() // GetOutputName() // 'PROseason.OUT'
    case(typeproject_TypePRM)
        totalname = GetPathNameOutp() // GetOutputName() // 'PRMseason.OUT'
    end select

    call fRun_open(totalname, 'w')
    write(tempstring, '(a)') GetAquaCropDescriptionWithTimeStamp()
    call fRun_write(trim(tempstring))

    call fRun_write('')
    call fRun_write('    RunNr     Day1   Month1    Year1     Rain      ETo       GD     CO2' // &
        '      Irri   Infilt   Runoff    Drain   Upflow        E     E/Ex       Tr      TrW   Tr/Trx' // &
        '    SaltIn   SaltOut    SaltUp  SaltProf' // &
        '     Cycle   SaltStr  FertStr  WeedStr  TempStr   ExpStr   StoStr' // &
        '  BioMass  Brelative   HI    Y(dry)  Y(fresh)    WPet      Bin     Bout     DayN   MonthN    YearN')
    call fRun_write('                                           mm       mm  degC.day    ppm' // &
        '        mm       mm       mm       mm       mm       mm        %       mm       mm        %' // &
        '    ton/ha    ton/ha    ton/ha    ton/ha' // &
        '      days       %        %        %        %        %        %  ' // &
        '  ton/ha        %       %    ton/ha   ton/ha    kg/m3   ton/ha   ton/ha')
end subroutine OpenOutputRun


subroutine OpenOutputDaily(TheProjectType)
    integer(intEnum), intent(in) :: TheProjectType

    character(len=:), allocatable :: totalname
    character(len=1025) :: tempstring

    select case (TheProjectType)
    case(typeproject_TypePRO)
        totalname = GetPathNameOutp() // GetOutputName() // 'PROday.OUT'
    case(typeproject_TypePRM)
        totalname = GetPathNameOutp() // GetOutputName() // 'PRMday.OUT'
    end select

    call fDaily_open(totalname, 'w')
    write(tempstring, '(a)') GetAquaCropDescriptionWithTimeStamp()
    call fDaily_write(trim(tempstring))
end subroutine OpenOutputDaily


subroutine OpenPart1MultResults(TheProjectType)
    integer(intEnum), intent(in) :: TheProjectType

    character(len=:), allocatable :: totalname
    character(len=1025) :: tempstring

    select case (TheProjectType)
    case(typeproject_TypePRO)
        totalname = GetPathNameOutp() // GetOutputName() // 'PROharvests.OUT'
    case(typeproject_TypePRM)
        totalname = GetPathNameOutp() // GetOutputName() // 'PRMharvests.OUT'
    end select
    call SetfHarvest_filename(totalname)
    call fHarvest_open(GetfHarvest_filename(), 'w')
    write(tempstring, '(a)') GetAquaCropDescriptionWithTimeStamp()
    call fHarvest_write(trim(tempstring))
    call fHarvest_write('Biomass and Yield at Multiple cuttings')
end subroutine OpenPart1MultResults


subroutine CreateDailyClimFiles(FromSimDay, ToSimDay)
    integer(int32), intent(in) :: FromSimDay
    integer(int32), intent(in) :: ToSimDay

    character(len=:), allocatable :: totalname, totalnameOUT
    integer :: fETo, fRain, fTemp, fEToS, fRainS, fTempS, rc
    character(len=255) :: StringREAD
    integer(int32) :: i
    integer(int32) :: RunningDay
    real(dp) :: tmpRain
    real(dp) :: ETo_temp
    type(rep_DayEventDbl), dimension(31) :: TminDataSet_temp, TmaxDataSet_temp
    real(dp) :: Tmin_temp, Tmax_temp
    type(rep_DayEventDbl), dimension(31) :: EToDataSet_temp, RainDataSet_temp

    ! 1. ETo file
    if (GetEToFile() /= '(None)') then
        totalname = GetEToFilefull()
        if (FileExists(totalname)) then
            ! open file and find first day of simulation period
            select case (GetEToRecord_DataType())
            case(datatype_Daily)
                open(newunit=fETo, file=trim(totalname), status='old', &
                                                         action='read')
                read(fETo, *, iostat=rc) ! description
                read(fETo, *, iostat=rc) ! time step
                read(fETo, *, iostat=rc) ! day
                read(fETo, *, iostat=rc) ! month
                read(fETo, *, iostat=rc) ! year
                read(fETo, *, iostat=rc)
                read(fETo, *, iostat=rc)
                read(fETo, *, iostat=rc)
                do i = GetEToRecord_FromDayNr(), (FromSimDay - 1)
                    read(fETo, *, iostat=rc)
                end do
                read(fETo, *, iostat=rc) ETo_temp
                call SetETo(ETo_temp)
            case(datatype_decadely)
                EToDataSet_temp = GetEToDataSet()
                call GetDecadeEToDataSet(FromSimDay, EToDataSet_temp)
                call SetEToDataSet(EToDataSet_temp)
                i = 1
                do while (GetEToDataSet_DayNr(i) /= FromSimDay)
                    i = i+1
                end do
                call SetETo(GetEToDataSet_Param(i))
            case(datatype_Monthly)
                EToDataSet_temp = GetEToDataSet()
                call GetMonthlyEToDataSet(FromSimDay, EToDataSet_temp)
                call SetEToDataSet(EToDataSet_temp)
                i = 1
                do while (GetEToDataSet_DayNr(i) /= FromSimDay)
                    i = i+1
                end do
                call SetETo(GetEToDataSet_Param(i))
            end select

            ! create SIM file and record first day
            totalnameOUT = GetPathNameSimul() // 'EToData.SIM'
            open(newunit=fEToS, file=trim(totalnameOUT), status='replace', &
                                                         action='write')
            write(fEToS, '(f10.4)') GetETo()
            ! next days of simulation period
            do RunningDay = (FromSimDay + 1), ToSimDay
                select case (GetEToRecord_DataType())
                case(datatype_Daily)
                    if (rc == iostat_end) then
                        rewind(fETo)
                        read(fETo, *, iostat=rc) ! description
                        read(fETo, *, iostat=rc) ! time step
                        read(fETo, *, iostat=rc) ! day
                        read(fETo, *, iostat=rc) ! month
                        read(fETo, *, iostat=rc) ! year
                        read(fETo, *, iostat=rc)
                        read(fETo, *, iostat=rc)
                        read(fETo, *, iostat=rc)
                        read(fETo, *) ETo_temp
                        call SetETo(ETo_temp)
                    else
                        read(fETo, *) ETo_temp
                        call SetETo(ETo_temp)
                    end if
                case(datatype_Decadely)
                    if (RunningDay > GetEToDataSet_DayNr(31)) then
                        EToDataSet_temp = GetEToDataSet()
                        call GetDecadeEToDataSet(RunningDay, EToDataSet_temp)
                        call SetEToDataSet(EToDataSet_temp)
                    end if
                    i = 1
                    do while (GetEToDataSet_DayNr(i) /= RunningDay)
                        i = i+1
                    end do
                    call SetETo(GetEToDataSet_Param(i))
                case(datatype_Monthly)
                    if (RunningDay > GetEToDataSet_DayNr(31)) then
                        EToDataSet_temp = GetEToDataSet()
                        call GetMonthlyEToDataSet(RunningDay, EToDataSet_temp)
                        call SetEToDataSet(EToDataSet_temp)
                    end if
                    i = 1
                    do while (GetEToDataSet_DayNr(i) /= RunningDay)
                        i = i+1
                    end do
                    call SetETo(GetEToDataSet_Param(i))
                end select
                write(fEToS, '(f10.4)') GetETo()
            end do
            ! Close files
            if (GetEToRecord_DataType() == datatype_Daily) then
                close(fETo)
            end if
            close(fEToS)
        end if
    end if

    ! 2. Rain File
    if (GetRainFile() /= '(None)') then
        totalname = GetRainFilefull()
        if (FileExists(totalname)) then
            ! open file and find first day of simulation period
            select case (GetRainRecord_DataType())
            case(datatype_Daily)
                open(newunit=fRain, file=trim(totalname), status='old', &
                                                          action='read')
                read(fRain, *, iostat=rc) ! description
                read(fRain, *, iostat=rc) ! time step
                read(fRain, *, iostat=rc) ! day
                read(fRain, *, iostat=rc) ! month
                read(fRain, *, iostat=rc) ! year
                read(fRain, *, iostat=rc)
                read(fRain, *, iostat=rc)
                read(fRain, *, iostat=rc)
                do i = GetRainRecord_FromDayNr(), (FromSimDay - 1)
                    read(fRain, *, iostat=rc)
                end do
                read(fRain, *, iostat=rc) tmpRain
                call SetRain(tmpRain)
            case(datatype_Decadely)
                RainDataSet_temp = GetRainDataSet()
                call GetDecadeRainDataSet(FromSimDay, RainDataSet_temp)
                call SetRainDataSet(RainDataSet_temp)
                i = 1
                do while (GetRainDataSet_DayNr(i) /= FromSimDay)
                    i = i+1
                end do
                call SetRain(GetRainDataSet_Param(i))
            case(datatype_Monthly)
                RainDataSet_temp = GetRainDataSet()
                call GetMonthlyRainDataSet(FromSimDay, RainDataSet_temp)
                call SetRainDataSet(RainDataSet_temp)
                i = 1
                do while (GetRainDataSet_DayNr(i) /= FromSimDay)
                    i = i+1
                end do
                call SetRain(GetRainDataSet_Param(i))
            end select

            ! create SIM file and record first day
            totalnameOUT = GetPathNameSimul() // 'RainData.SIM'
            open(newunit=fRainS, file=trim(totalnameOUT), status='replace', &
                                                          action='write')
            write(fRainS, '(f10.4)') GetRain()
            ! next days of simulation period
            do RunningDay = (FromSimDay + 1), ToSimDay
                select case (GetRainRecord_DataType())
                case(datatype_daily)
                    if (rc == iostat_end) then
                        rewind(fRain)
                        read(fRain, *, iostat=rc) ! description
                        read(fRain, *, iostat=rc) ! time step
                        read(fRain, *, iostat=rc) ! day
                        read(fRain, *, iostat=rc) ! month
                        read(fRain, *, iostat=rc) ! year
                        read(fRain, *, iostat=rc)
                        read(fRain, *, iostat=rc)
                        read(fRain, *, iostat=rc)
                        read(fRain, *, iostat=rc) tmpRain
                        call SetRain(tmpRain)
                    else
                        read(fRain, *, iostat=rc) tmpRain
                        call SetRain(tmpRain)
                    end if
                case(datatype_Decadely)
                    if (RunningDay > GetRainDataSet_DayNr(31)) then
                        RainDataSet_temp = GetRainDataSet()
                        call GetDecadeRainDataSet(RunningDay, RainDataSet_temp)
                        call SetRainDataSet(RainDataSet_temp)
                    end if
                    i = 1
                    do while (GetRainDataSet_DayNr(i) /= RunningDay)
                        i = i+1
                    end do
                    call SetRain(GetRainDataSet_Param(i))
                case(datatype_monthly)
                    if (RunningDay > GetRainDataSet_DayNr(31)) then
                        RainDataSet_temp = GetRainDataSet()
                        call GetMonthlyRainDataSet(RunningDay, RainDataSet_temp)
                        call SetRainDataSet(RainDataSet_temp)
                    end if
                    i = 1
                    do while (GetRainDataSet_DayNr(i) /= RunningDay)
                        i = i+1
                    end do
                    call SetRain(GetRainDataSet_Param(i))
                end select
                write(fRainS, '(f10.4)') GetRain()
            end do
            ! Close files
            if (GetRainRecord_DataType() == datatype_Daily) then
                close(fRain)
            end if
            close(fRainS)
        end if
    end if

    ! 3. Temperature file
    if (GetTemperatureFile() /= '(None)') then
        totalname = GetTemperatureFilefull()
        if (FileExists(totalname)) then
            ! open file and find first day of simulation period
            select case (GetTemperatureRecord_DataType())
            case(datatype_daily)
                open(newunit=fTemp, file=trim(totalname), status='old', &
                                                          action='read')
                read(fTemp, *, iostat=rc) ! description
                read(fTemp, *, iostat=rc) ! time step
                read(fTemp, *, iostat=rc) ! day
                read(fTemp, *, iostat=rc) ! month
                read(fTemp, *, iostat=rc) ! year
                read(fTemp, *, iostat=rc)
                read(fTemp, *, iostat=rc)
                read(fTemp, *, iostat=rc)
                do i = GetTemperatureRecord_FromDayNr(), (FromSimDay - 1)
                    read(fTemp, *, iostat=rc)
                end do
                read(fTemp, '(a)', iostat=rc) StringREAD  ! i.e. DayNri
                Tmin_temp = GetTmin()
                Tmax_temp = GetTmax()
                call SplitStringInTwoParams(StringREAD, Tmin_temp, Tmax_temp)
                call SetTmin(Tmin_temp)
                call SetTmax(Tmax_temp)
            case(datatype_Decadely)
                TminDataSet_temp = GetTminDataSet()
                TmaxDataSet_temp = GetTmaxDataSet()
                call GetDecadeTemperatureDataSet(FromSimDay, TminDataSet_temp, &
                                                              TmaxDataSet_temp)
                call SetTminDataSet(TminDataSet_temp)
                call SetTmaxDataSet(TmaxDataSet_temp)
                i = 1
                do while (GetTminDataSet_DayNr(i) /= FromSimDay)
                    i = i+1
                end do
                call SetTmin(GetTminDataSet_Param(i))
                call SetTmax(GetTmaxDataSet_Param(i))
            case(datatype_Monthly)
                TminDataSet_temp = GetTminDataSet()
                TmaxDataSet_temp = GetTmaxDataSet()
                call GetMonthlyTemperatureDataSet(FromSimDay, TminDataSet_temp, &
                                                              TmaxDataSet_temp)
                call SetTminDataSet(TminDataSet_temp)
                call SetTmaxDataSet(TmaxDataSet_temp)
                i = 1
                do while (GetTminDataSet_DayNr(i) /= FromSimDay)
                    i = i+1
                end do
                call SetTmin(GetTminDataSet_Param(i))
                call SetTmax(GetTmaxDataSet_Param(i))
            end select

            ! create SIM file and record first day
            totalnameOUT = GetPathNameSimul() // 'TempData.SIM'
            open(newunit=fTempS, file=trim(totalnameOUT), status='replace', &
                                                          action='write')
            write(fTempS, '(2f10.4)') GetTmin(), GetTmax()
            ! next days of simulation period
            do RunningDay = (FromSimDay + 1), ToSimDay
                select case (GetTemperatureRecord_Datatype())
                case(datatype_Daily)
                    if (rc == iostat_end) then
                        rewind(fTemp)
                        read(fTemp, *, iostat=rc) ! description
                        read(fTemp, *, iostat=rc) ! time step
                        read(fTemp, *, iostat=rc) ! day
                        read(fTemp, *, iostat=rc) ! month
                        read(fTemp, *, iostat=rc) ! year
                        read(fTemp, *, iostat=rc)
                        read(fTemp, *, iostat=rc)
                        read(fTemp, *, iostat=rc)
                        read(fTemp, '(a)', iostat=rc) StringREAD
                        Tmin_temp = GetTmin()
                        Tmax_temp = GetTmax()
                        call SplitStringInTwoParams(StringREAD, Tmin_temp, &
                                                                Tmax_temp)
                        call SetTmin(Tmin_temp)
                        call SetTmax(Tmax_temp)
                    else
                        read(fTemp, *, iostat=rc) Tmin_temp, Tmax_temp
                        call SetTmin(Tmin_temp)
                        call SetTmax(Tmax_temp)
                    end if
                case(datatype_Decadely)
                    if (RunningDay > GetTminDataSet_DayNr(31)) then
                        TminDataSet_temp = GetTminDataSet()
                        TmaxDataSet_temp = GetTmaxDataSet()
                        call GetDecadeTemperatureDataSet(RunningDay, &
                                                         TminDataSet_temp, &
                                                         TmaxDataSet_temp)
                        call SetTminDataSet(TminDataSet_temp)
                        call SetTmaxDataSet(TmaxDataSet_temp)
                    end if
                    i = 1
                    do while (GetTminDataSet_DayNr(i) /= RunningDay)
                        i = i+1
                    end do
                    call SetTmin(GetTminDataSet_Param(i))
                    call SetTmax(GetTmaxDataSet_Param(i))
                case(datatype_Monthly)
                    if (RunningDay > GetTminDataSet_DayNr(31)) then
                        TminDataSet_temp = GetTminDataSet()
                        TmaxDataSet_temp = GetTmaxDataSet()
                        call GetMonthlyTemperatureDataSet(RunningDay, &
                                                          TminDataSet_temp, &
                                                          TmaxDataSet_temp)
                        call SetTminDataSet(TminDataSet_temp)
                        call SetTmaxDataSet(TmaxDataSet_temp)
                    end if
                    i = 1
                    do while (GetTminDataSet_DayNr(i) /= RunningDay)
                        i = i+1
                    end do
                    call SetTmin(GetTminDataSet_Param(i))
                    call SetTmax(GetTmaxDataSet_Param(i))
                end select
                write(fTempS, '(2f10.4)') GetTmin(), GetTmax()
            end do
            ! Close files
            if (GetTemperatureRecord_datatype() == datatype_Daily) then
                close(fTemp)
            end if
            close(fTempS)
        end if
    end if
end subroutine CreateDailyClimFiles


subroutine OpenHarvestInfo()
    character(len=:), allocatable :: totalname
    integer(int8) :: i
    character(len=:), allocatable :: TempString

    if (GetManFile() /= '(None)') then
        totalname = GetManFileFull()
    else
        totalname = trim(GetPathNameSimul())//'Cuttings.AqC'
    end if
    call fCuts_open(totalname, 'r')
    TempString = fCuts_read() ! description
    TempString = fCuts_read() ! AquaCrop version
    if (GetManFile() /= '(None)') then
        do i= 1, 10
            TempString = fCuts_read() ! management info
       end do
    end if
    do i = 1, 12
        TempString = fCuts_read()  ! cuttings info (already loaded)
    end do
    call GetNextHarvest
end subroutine OpenHarvestInfo


subroutine OpenClimFilesAndGetDataFirstDay(FirstDayNr)
    integer(int32), intent(in) :: FirstDayNr

    character(len=:), allocatable :: totalname
    integer(int32) :: i
    real(dp) :: tmpRain, ETo_temp
    real(dp) :: Tmin_temp, Tmax_temp
    character(len=1025) :: TempString

    ! ETo file
    if (GetEToFile() /= '(None)') then
        totalname = trim(GetPathNameSimul())//'EToData.SIM'
        call fEToSIM_open(totalname, 'r')
        if (FirstDayNr == GetSimulation_FromDayNr()) then
            TempString = fEToSIM_read()
            read(TempString, *) ETo_temp
            call SetETo(ETo_temp)
        else
            do i = GetSimulation_FromDayNr(), (FirstDayNr - 1)
                TempString = fEToSIM_read()
                read(TempString, *) ETo_temp
                call SetETo(ETo_temp)
            end do
            TempString = fEToSIM_read()
            read(TempString, *) ETo_temp
            call SetETo(ETo_temp)
        end if
    end if
    ! Rain file
    if (GetRainFile() /= '(None)') then
        totalname = trim(GetPathNameSimul())//'RainData.SIM'
        call fRainSIM_open(totalname, 'r')
        if (FirstDayNr == GetSimulation_FromDayNr()) then
            TempString = fRainSIM_read()
            read(TempString, *) tmpRain
            call SetRain(tmpRain)
        else
            do i = GetSimulation_FromDayNr(), (FirstDayNr - 1)
                TempString = fRainSIM_read()
                read(TempString, *) tmpRain
                call SetRain(tmpRain)
            end do
            TempString = fRainSIM_read()
            read(TempString, *) tmpRain
            call SetRain(tmpRain)
        end if
    end if
    ! Temperature file
    if (GetTemperatureFile() /= '(None)') then
        totalname = trim(GetPathNameSimul())//'TempData.SIM'
        call fTempSIM_open(totalname, 'r')
        if (FirstDayNr == GetSimulation_FromDayNr()) then
            TempString = fTempSIM_read()
            read(TempString, *) Tmin_temp, Tmax_temp
            call SetTmin(Tmin_temp)
            call SetTmax(Tmax_temp)
        else
            do i = GetSimulation_FromDayNr(), (FirstDayNr - 1)
                TempString = fTempSIM_read()
                read(TempString, *) Tmin_temp, Tmax_temp
                call SetTmin(Tmin_temp)
                call SetTmax(Tmax_temp)
            end do
            TempString = fTempSIM_read()
            read(TempString, *) Tmin_temp, Tmax_temp
            call SetTmin(Tmin_temp)
            call SetTmax(Tmax_temp)
        end if
    else
        call SetTmin(GetSimulParam_Tmin())
        call SetTmax(GetSimulParam_Tmax())
    end if
end subroutine OpenClimFilesAndGetDataFirstDay


subroutine InitializeSimulation(TheProjectFileStr, TheProjectType)
    character(len=*), intent(in) :: TheProjectFileStr
    integer(intenum), intent(in) :: TheProjectType

    call SetTheProjectFile(trim(TheProjectFileStr))
    call OpenOutputRun(TheProjectType) ! open seasonal results .out
    if (GetOutDaily()) then
        call OpenOutputDaily(TheProjectType)  ! Open Daily results .OUT
    end if
    if (GetPart1Mult()) then
        call OpenPart1MultResults(TheProjectType) ! Open Multiple harvests in season .OUT
    end if
end subroutine InitializeSimulation


subroutine FinalizeSimulation()

    call fRun_close() ! Close Run.out
    if (GetOutDaily()) then
        call fDaily_close()  ! Close Daily.OUT
    end if
    if (GetPart1Mult()) then
        call fHarvest_close()  ! Close Multiple harvests in season
    end if
end subroutine FinalizeSimulation


subroutine WriteSimPeriod(NrRun, TheProjectFile)
    integer(int8), intent(in) :: NrRun
    character(len=*), intent(in) :: TheProjectFile

    integer(int32) :: Day1, Month1, Year1, DayN, MonthN, YearN

    call DetermineDate(GetSimulation_FromDayNr(), Day1, Month1, Year1)
    ! Start simulation run
    call DetermineDate(GetSimulation_ToDayNr(), DayN, MonthN, YearN)
    ! End simulation run
    call WriteTheResults(NrRun, Day1, Month1, Year1, DayN, MonthN, YearN, &
                        GetSumWaBal_Rain(), GetSumETo(), GetSumGDD(), &
                        GetSumWaBal_Irrigation(), GetSumWaBal_Infiltrated(), &
                        GetSumWaBal_Runoff(), GetSumWaBal_Drain(), &
                        GetSumWaBal_CRwater(), GetSumWaBal_Eact(), &
                        GetSumWaBal_Epot(), GetSumWaBal_Tact(), &
                        GetSumWaBal_TrW(), GetSumWaBal_Tpot(), &
                        GetSumWaBal_SaltIn(), GetSumWaBal_SaltOut(), &
                        GetSumWaBal_CRsalt(), GetSumWaBal_Biomass(), &
                        GetSumWaBal_BiomassUnlim(), GetTransfer_Bmobilized(), &
                        GetSimulation_Storage_Btotal(), TheProjectFile)
end subroutine WriteSimPeriod


subroutine WriteIntermediatePeriod(TheProjectFile)
    character(len=*), intent(in) :: TheProjectFile

    integer(int32) :: Day1, Month1, Year1, DayN, MonthN, YearN
    real(dp) :: RPer, EToPer, GDDPer, IrriPer, InfiltPer, EPer, &
                ExPer, TrPer, TrWPer, TrxPer, DrainPer, BiomassPer, BUnlimPer
    real(dp) :: ROPer, CRwPer, SalInPer, SalOutPer, SalCRPer, BmobPer, BstoPer

    ! determine intermediate results
    call DetermineDate((GetPreviousDayNr()+1), Day1, Month1, Year1)
    call DetermineDate(GetDayNri(), DayN, MonthN, YearN)
    RPer = GetSumWaBal_Rain() - GetPreviousSum_Rain()
    EToPer = GetSumETo() - GetPreviousSumETo()
    GDDPer = GetSumGDD() - GetPreviousSumGDD()
    IrriPer = GetSumWaBal_Irrigation() - GetPreviousSum_Irrigation()
    InfiltPer = GetSumWaBal_Infiltrated() - GetPreviousSum_Infiltrated()
    EPer = GetSumWaBal_Eact() - GetPreviousSum_Eact()
    ExPer = GetSumWaBal_Epot() - GetPreviousSum_Epot()
    TrPer = GetSumWaBal_Tact() - GetPreviousSum_Tact()
    TrWPer = GetSumWaBal_TrW() - GetPreviousSum_TrW()
    TrxPer = GetSumWaBal_Tpot() - GetPreviousSum_Tpot()
    DrainPer = GetSumWaBal_Drain() - GetPreviousSum_Drain()
    BiomassPer = GetSumWaBal_Biomass() - GetPreviousSum_Biomass()
    BUnlimPer = GetSumWaBal_BiomassUnlim() - GetPreviousSum_BiomassUnlim()

    ROPer = GetSumWaBal_Runoff() - GetPreviousSum_Runoff()
    CRwPer = GetSumWaBal_CRwater() - GetPreviousSum_CRwater()
    SalInPer = GetSumWaBal_SaltIn() - GetPreviousSum_SaltIn()
    SalOutPer = GetSumWaBal_SaltOut() - GetPreviousSum_SaltOut()
    SalCRPer = GetSumWaBal_CRsalt() - GetPreviousSum_CRsalt()

    BmobPer = GetTransfer_Bmobilized() - GetPreviousBmob()
    BstoPer = GetSimulation_Storage_Btotal() - GetPreviousBsto()

    ! write
    call WriteTheResults(int(undef_int, kind=int8), Day1, Month1, Year1, DayN, &
                         MonthN, YearN, RPer, EToPer, GDDPer, IrriPer, InfiltPer, &
                         ROPer, DrainPer, CRwPer, EPer, ExPer, TrPer, TrWPer, &
                         TrxPer, SalInPer, SalOutPer, SalCRPer, BiomassPer, &
                         BUnlimPer, BmobPer, BstoPer, TheProjectFile)

    ! reset previous sums
    call SetPreviousDayNr(GetDayNri())
    call SetPreviousSum_Rain(GetSumWaBal_Rain())
    call SetPreviousSumETo(GetSumETo())
    call SetPreviousSumGDD(GetSumGDD())
    call SetPreviousSum_Irrigation(GetSumWaBal_Irrigation())
    call SetPreviousSum_Infiltrated(GetSumWaBal_Infiltrated())
    call SetPreviousSum_Eact(GetSumWaBal_Eact())
    call SetPreviousSum_Epot(GetSumWaBal_Epot())
    call SetPreviousSum_Tact(GetSumWaBal_Tact())
    call SetPreviousSum_TrW(GetSumWaBal_TrW())
    call SetPreviousSum_Tpot(GetSumWaBal_Tpot())
    call SetPreviousSum_Drain(GetSumWaBal_Drain())
    call SetPreviousSum_Biomass(GetSumWaBal_Biomass())
    call SetPreviousSum_BiomassPot(GetSumWaBal_BiomassPot())
    call SetPreviousSum_BiomassUnlim(GetSumWaBal_BiomassUnlim())

    call SetPreviousSum_Runoff(GetSumWaBal_Runoff())
    call SetPreviousSum_CRwater(GetSumWaBal_CRwater())
    call SetPreviousSum_SaltIn(GetSumWaBal_SaltIn())
    call SetPreviousSum_SaltOut(GetSumWaBal_SaltOut())
    call SetPreviousSum_CRsalt(GetSumWaBal_CRsalt())

    call SetPreviousBmob(GetTransfer_Bmobilized())
    call SetPreviousBsto(GetSimulation_Storage_Btotal())
end subroutine WriteIntermediatePeriod


!! ===BEGIN Subroutines and functions for AdvanceOneTimeStep ===

subroutine GetZandECgwt(ZiAqua, ECiAqua)
    integer(int32), intent(inout) :: ZiAqua
    real(dp), intent(inout) :: ECiAqua

    integer(int32) :: ZiIN
    type(CompartmentIndividual), dimension(max_No_compartments) :: Comp_temp

    ZiIN = ZiAqua
    if (GetGwTable_DNr1() == GetGwTable_DNr2()) then
        ZiAqua  = GetGwTable_Z1()
        ECiAqua = GetGwTable_EC1()
    else
        ZiAqua = GetGwTable_Z1() + &
                 roundc((GetDayNri() - GetGwTable_DNr1())* &
                 (GetGwTable_Z2() - GetGwTable_Z1())/ &
                 real(GetGwTable_DNr2() - GetGwTable_DNr1(), kind=dp), mold = 1)
        ECiAqua = GetGwTable_EC1() + &
                 (GetDayNri() - GetGwTable_DNr1())* &
                 (GetGwTable_EC2() - GetGwTable_EC1())/&
                 real(GetGwTable_DNr2() - GetGwTable_DNr1(), kind=dp)
    end if
    if (ZiAqua /= ZiIN) then
        Comp_temp = GetCompartment()
        call CalculateAdjustedFC((ZiAqua/100._dp), Comp_temp)
        call SetCompartment(Comp_temp)
    end if
end subroutine GetZandECgwt


integer(int32) function IrriOutSeason()
    integer(int32) :: DNr, Nri, i
    type(Rep_DayEventInt), dimension(5) :: IrriEvents
    logical :: TheEnd

    DNr = GetDayNri() - GetSimulation_FromDayNr() + 1
    do i = 1, 5
        IrriEvents(i) = GetIrriBeforeSeason_i(i)
    end do
    if (GetDayNri() > GetCrop_DayN()) then
        DNr = GetDayNri() - GetCrop_DayN()
        do i = 1, 5
            IrriEvents(i) = GetIrriAfterSeason_i(i)
        end do
    end if
    if (DNr < 1) then
        IrriOutSeason = 0
    else
        TheEnd = .false.
        Nri = 0
        loop: do
            Nri = Nri + 1
            if (IrriEvents(Nri)%DayNr == DNr) then
                IrriOutSeason = IrriEvents(Nri)%Param
                TheEnd = .true.
            else
                IrriOutSeason = 0
            end if
            if ((Nri == 5) .or. (IrriEvents(Nri)%DayNr == 0) &
              .or. (IrriEvents(Nri)%DayNr > DNr) &
              .or. (TheEnd)) exit loop
        end do loop
    end if
end function IrriOutSeason


integer(int32) function IrriManual()
    integer(int32) :: DNr
    character(len=:), allocatable :: StringREAD
    real(dp) :: Ir1, Ir2
    real(dp) :: IrriECw_temp

    if (GetIrriFirstDayNr() == undef_int) then
        DNr = GetDayNri() - GetCrop_Day1() + 1
    else
        DNr = GetDayNri() - GetIrriFirstDayNr() + 1
    end if
    if (GetIrriInfoRecord1_NoMoreInfo()) then
        IrriManual = 0
    else
        IrriManual = 0
        if (GetIrriInfoRecord1_TimeInfo() == DNr) then
            IrriManual = GetIrriInfoRecord1_DepthInfo()
            StringREAD = fIrri_read()
            if (fIrri_eof()) then
                call SetIrriInfoRecord1_NoMoreInfo(.true.)
            else
                call SetIrriInfoRecord1_NoMoreInfo(.false.)
                if (GetGlobalIrriECw()) then ! Versions before 3.2
                    call SplitStringInTwoParams(StringREAD, Ir1, Ir2)
                else
                    IrriECw_temp = GetSimulation_IrriECw()
                    call SplitStringInThreeParams(StringREAD, &
                            Ir1, Ir2, IrriECw_temp)
                    call SetSimulation_IrriECw(IrriECw_temp)
                end if
                call SetIrriInfoRecord1_TimeInfo(roundc(Ir1, mold=1))
                call SetIrriInfoRecord1_DepthInfo(roundc(Ir2, mold=1))
            end if
        end if
    end if
end function IrriManual


subroutine GetIrriParam(TargetTimeVal, TargetDepthVal)
    integer(int32), intent(inout) :: TargetTimeVal
    integer(int32), intent(inout) :: TargetDepthVal

    integer(int32) :: DayInSeason
    integer(int32) :: FromDay_temp, TimeInfo_temp, DepthInfo_temp
    real(dp)       :: IrriECw_temp
    character(len=:), allocatable :: TempString

    TargetTimeVal = -999
    TargetDepthVal = -999
    if ((GetDayNri() < GetCrop_Day1()) .or. &
        (GetDayNri() > GetCrop_DayN())) then
        call SetIrrigation(real(IrriOutSeason(), kind=dp))
    elseif (GetIrriMode() == IrriMode_Manual) then
        call SetIrrigation(real(IrriManual(), kind=dp))
    end if
    if ((GetIrriMode() == IrriMode_Generate) .and. &
        ((GetDayNri() >= GetCrop_Day1()) .and. &
         (GetDayNri() <= GetCrop_DayN()))) then
        ! read next line if required
        DayInSeason = GetDayNri() - GetCrop_Day1() + 1
        if (DayInSeason > GetIrriInfoRecord1_ToDay()) then
            ! read next line
            call SetIrriInfoRecord1(GetIrriInfoRecord2())

            TempString = fIrri_read()
            if (fIrri_eof()) then
                call SetIrriInfoRecord1_ToDay(GetCrop_DayN() - &
                        GetCrop_Day1() + 1)
            else
                call SetIrriInfoRecord2_NoMoreInfo(.false.)
                if (GetGlobalIrriECw()) then ! Versions before 3.2
                    read(TempString,*) FromDay_temp, &
                      TimeInfo_temp, DepthInfo_temp
                    call SetIrriInfoRecord2_FromDay(FromDay_temp)
                    call SetIrriInfoRecord2_TimeInfo(TimeInfo_temp)
                    call SetIrriInfoRecord2_DepthInfo(DepthInfo_temp)
                else
                    read(TempString,*) FromDay_temp, &
                      TimeInfo_temp, DepthInfo_temp, IrriEcw_temp
                    call SetIrriInfoRecord2_FromDay(FromDay_temp)
                    call SetIrriInfoRecord2_TimeInfo(TimeInfo_temp)
                    call SetIrriInfoRecord2_DepthInfo(DepthInfo_temp)
                    call SetSimulation_IrriEcw(IrriEcw_temp)
                end if
                call SetIrriInfoRecord1_ToDay(GetIrriInfoRecord2_FromDay() - 1)
            end if
        end if
        ! get TargetValues
        TargetDepthVal = GetIrriInfoRecord1_DepthInfo()
        select case (GetGenerateTimeMode())
        case (GenerateTimeMode_AllDepl)
            TargetTimeVal = GetIrriInfoRecord1_TimeInfo()
        case (GenerateTimeMode_AllRAW)
            TargetTimeVal = GetIrriInfoRecord1_TimeInfo()
        case (GenerateTimeMode_FixInt)
            TargetTimeVal = GetIrriInfoRecord1_TimeInfo()
            if (TargetTimeVal > GetIrriInterval()) then ! do not yet irrigate
                TargetTimeVal = 0
            elseif (TargetTimeVal == GetIrriInterval()) then ! irrigate
                TargetTimeVal = 1
            else
                ! still to solve
                TargetTimeVal = 1 ! temporary solution
            end if
            if ((TargetTimeVal == 1) .and. &
                (GetGenerateDepthMode() == GenerateDepthMode_FixDepth)) then
                call SetIrrigation(real(TargetDepthVal, kind=dp))
            end if
        case (GenerateTimeMode_WaterBetweenBunds)
            TargetTimeVal = GetIrriInfoRecord1_TimeInfo()
            if ( (GetManagement_BundHeight() >= 0.01_dp) &
                .and. (GetGenerateDepthMode() == GenerateDepthMode_FixDepth) &
                .and. (TargetTimeVal < (1000._dp * GetManagement_BundHeight())) &
                .and. (TargetTimeVal >= roundc(GetSurfaceStorage(), mold=1))) then
                call SetIrrigation(real(TargetDepthVal, kind=dp))
            else
                call SetIrrigation(0._dp)
            end if
            TargetTimeVal = -999 ! no need for check in SIMUL
       end select
    end if
end subroutine GetIrriParam


subroutine AdjustSWCRootZone(PreIrri)
    real(dp), intent(inout) :: PreIrri

    integer(int32) :: compi, layeri
    real(dp) :: SumDepth, ThetaPercRAW

    compi = 0
    SumDepth = 0
    PreIrri = 0._dp
    loop: do
        compi = compi + 1
        SumDepth = SumDepth + GetCompartment_Thickness(compi)
        layeri = GetCompartment_Layer(compi)
        ThetaPercRaw = GetSoilLayer_FC(layeri)/100._dp &
                        - GetSimulParam_PercRAW()/100._dp &
                        * GetCrop_pdef() &
                        * (GetSoilLayer_FC(layeri)/100._dp &
                                - GetSoilLayer_WP(layeri)/100._dp)
        if (GetCompartment_Theta(compi) < ThetaPercRaw) then
            PreIrri = PreIrri &
                      + (ThetaPercRaw - GetCompartment_Theta(compi)) &
                      *1000._dp*GetCompartment_Thickness(compi)
            call SetCompartment_Theta(compi, ThetaPercRaw)
        end if
        if ((SumDepth >= GetRootingDepth()) &
                .or. (compi == GetNrCompartments())) exit loop
    end do loop
end subroutine AdjustSWCRootZone


subroutine InitializeTransferAssimilates(Bin, Bout, AssimToMobilize, &
                                         AssimMobilized, FracAssim, &
                                         StorageON, MobilizationON, &
                                         HarvestNow)
    real(dp), intent(inout) :: Bin
    real(dp), intent(inout) :: Bout
    real(dp), intent(inout) :: AssimToMobilize
    real(dp), intent(inout) :: AssimMobilized
    real(dp), intent(inout) :: FracAssim
    logical, intent(inout) :: StorageON
    logical, intent(inout) :: MobilizationON
    logical, intent(in) :: HarvestNow

    real(dp) :: FracSto, tmob

    Bin = 0._dp
    Bout = 0._dp
    FracAssim = 0._dp
    if (GetCrop_subkind() == subkind_Forage) then
        ! only for perennial herbaceous forage crops
        FracAssim = 0._dp
        if (GetNoMoreCrop()) then
            StorageOn = .false.
            MobilizationOn = .false.
        else
            ! Start of storage period ?
            if ((GetDayNri() - GetSimulation_DelayedDays() - GetCrop_Day1() + 1) &
                == (GetCrop_DaysToHarvest() - GetCrop_Assimilates_Period() + 1)) then
                ! switch storage on
                StorageOn = .true.
                ! switch mobilization off
                if (MobilizationOn) then
                    AssimToMobilize = AssimMobilized
                end if
                MobilizationOn = .false.
            end if
            ! Fraction of assimilates transferred
            if (MobilizationOn) then
                tmob = (AssimToMobilize-AssimMobilized)/AssimToMobilize
                if (AssimToMobilize > AssimMobilized) then
                    FracAssim = (exp(-5._dp*tmob) - 1._dp)/(exp(-5._dp) - 1._dp)
                    if (GetCCiActual() > (0.9_dp*(GetCCxTotal() &
                        * (1._dp - real(GetSimulation_EffectStress_RedCCX(), kind=dp)/100._dp)))) then
                        FracAssim = FracAssim * ((GetCCxTotal()*(1._dp &
                        - real(GetSimulation_EffectStress_RedCCX(), kind=dp)/100._dp)) &
                        - GetCCiActual()) &
                        / (0.1*(GetCCxTotal()*(1._dp &
                        - real(GetSimulation_EffectStress_RedCCX(), kind=dp)/100._dp)))
                    end if
                    if (FracAssim < epsilon(0._dp)) then
                        FracAssim = 0._dp
                    end if
                else
                    ! everything is mobilized
                    FracAssim = 0._dp
                end if
            end if

            if ((StorageOn) .and. (GetCrop_Assimilates_Period() > 0)) then
                if (HarvestNow) then
                    FracSto = 0._dp
                else
                    if ((GetCCiActual() > GetManagement_Cuttings_CCcut()/100._dp) &
                        .and. (GetCCiActual() < (GetCCxTotal()*(1._dp &
                        - real(GetSimulation_EffectStress_RedCCX(),kind=dp)/100._dp)))) then
                        FracSto = (GetCCiActual() - GetManagement_Cuttings_CCcut()/100._dp) &
                            /((GetCCxTotal()*(1._dp - real(GetSimulation_EffectStress_RedCCX(),kind=dp)/100._dp)) &
                            - GetManagement_Cuttings_CCcut()/100._dp)
                    else 
                        FracSto = 1._dp
                    end if
                end if
                ! Use convex function
                FracAssim = FracSto * (GetCrop_Assimilates_Stored()/100._dp) &
                    * (1._dp-KsAny((((GetDayNri() - GetSimulation_DelayedDays() &
                    - GetCrop_Day1() + 1._dp) &
                    - (GetCrop_DaysToHarvest()-GetCrop_Assimilates_Period())) &
                    /real(GetCrop_Assimilates_Period(),kind=dp)),0._dp,1._dp,-5._dp))
            end if
            if (FracAssim < 0._dp) then
                FracAssim = 0._dp
            end if
            if (FracAssim > 1._dp) then
                FracAssim = 1._dp
            end if
        end if
    end if
end subroutine InitializeTransferAssimilates


subroutine GetPotValSF(DAP, SumGDDAdjCC, PotValSF)
    integer(int32), intent(in) :: DAP
    real(dp), intent(in) :: SumGDDAdjCC
    real(dp), intent(inout) :: PotValSF

    real(dp) :: RatDGDD

    RatDGDD = 1._dp
    if ((GetCrop_ModeCycle() == modecycle_GDDays) &
        .and. (GetCrop_GDDaysToFullCanopySF() < GetCrop_GDDaysToSenescence())) then
        RatDGDD = (GetCrop_DaysToSenescence()-GetCrop_DaysToFullCanopySF()) &
                    /(GetCrop_GDDaysToSenescence()-GetCrop_GDDaysToFullCanopySF())
    end if

    PotValSF = CCiNoWaterStressSF(DAP, GetCrop_DaysToGermination(), &
                    GetCrop_DaysToFullCanopySF(), GetCrop_DaysToSenescence(), &
                    GetCrop_DaysToHarvest(), GetCrop_GDDaysToGermination(), &
                    GetCrop_GDDaysToFullCanopySF(), GetCrop_GDDaysToSenescence(), &
                    GetCrop_GDDaysToHarvest(), GetCCoTotal(), GetCCxTotal(), &
                    GetCrop_CGC(), GetCrop_GDDCGC(), GetCDCTotal(), &
                    GetGDDCDCTotal(), SumGDDadjCC, RatDGDD, &
                    GetSimulation_EffectStress_RedCGC(), &
                    GetSimulation_EffectStress_RedCCX(), &
                    GetSimulation_EffectStress_CDecline(), GetCrop_ModeCycle())
    PotValSF = 100._dp * (1._dp/GetCCxCropWeedsNoSFstress()) * PotValSF
end subroutine GetPotValSF

!! ===END Subroutines and functions for AdvanceOneTimeStep ===


subroutine WriteEvaluationData(DAP)
    integer(int32), intent(in) :: DAP

    real(dp) :: SWCi, CCfield, CCstd, Bfield, Bstd, SWCfield, SWCstd
    integer(int32) :: Nr, Di, Mi, Yi
    character(len=1024) :: TempString
    integer(int32) :: DayNrEval_temp, DAP_temp

    ! 1. Prepare field data
    CCfield = real(undef_int, kind=dp)
    CCstd = real(undef_int, kind=dp)
    Bfield = real(undef_int, kind=dp)
    Bstd = real(undef_int, kind=dp)
    SWCfield = real(undef_int, kind=dp)
    SWCstd = real(undef_int, kind=dp)
    if ((GetLineNrEval() /= undef_int) .and. &
        (GetDayNrEval() == GetDayNri())) then
        ! read field data
        call fObs_rewind()
        do Nr = 1, (GetLineNrEval() -1)
            TempString = fObs_read()
        end do
        TempString = fObs_read()
        read(TempString, *) Nr, CCfield, CCstd, Bfield, Bstd, SWCfield, SWCstd
        ! get Day Nr for next field data
        TempString = fObs_read()
        if (fObs_eof()) then
            call SetLineNrEval(undef_int)
            call fObs_close()
        else
            call SetLineNrEval(GetLineNrEval() + 1)
            read(TempString, *) DayNrEval_temp
            call SetDayNrEval(DayNrEval_temp)
            call SetDayNrEval(GetDayNr1Eval() + GetDayNrEval() -1)
        end if
    end if
    ! 2. Date
    DAP_temp = DAP
    call DetermineDate(GetDayNri(), Di, Mi, Yi)
    if (GetClimRecord_FromY() == 1901) then
        Yi = Yi - 1901 + 1
    end if
    if (GetStageCode() == 0) then
        DAP_temp = undef_int ! before or after cropping
    end if
    ! 3. Write simulation results and field data
    SWCi = SWCZsoil(GetZeval())
    write(TempString, '(4i6, i5, 3f8.1, 3f10.3, 3f8.1)') &
           Di, Mi, Yi, DAP_temp, GetStageCode(), &
          (GetCCiActual()*100._dp), CCfield, CCstd, &
          GetSumWaBal_Biomass(), Bfield, Bstd, SWCi, SWCfield, SWCstd
    call fEval_write(trim(TempString))


    contains


    real(dp) function SWCZsoil(Zsoil)
        real(dp), intent(in) :: Zsoil

        integer(int32) :: compi
        real(dp) :: CumDepth, Factor, frac_value, SWCact

        CumDepth = 0._dp
        compi = 0
        SWCact = 0._dp
        loop : do
            compi = compi + 1
            CumDepth = CumDepth + GetCompartment_Thickness(compi)
            if (CumDepth <= Zsoil) then
                Factor = 1._dp
            else
                frac_value = Zsoil - (CumDepth - &
                             GetCompartment_Thickness(compi))
                if (frac_value > 0._dp) then
                    Factor = frac_value/GetCompartment_Thickness(compi)
                else
                     Factor = 0._dp
                end if
            end if
            SWCact = SWCact + Factor * 10._dp * &
                     (GetCompartment_Theta(compi)*100._dp) * &
                      GetCompartment_Thickness(compi)
            if ((roundc(100._dp*CumDepth, mold=1) >= &
                   roundc(100._dp*ZSoil, mold=1)) .or. &
                  (compi == GetNrCompartments())) exit loop
        end do loop
        SWCZsoil = SWCact
    end function SWCZsoil
end subroutine WriteEvaluationData


subroutine InitializeRunPart1(NrRun, TheProjectType)
    !! Part1 (before reading the climate) of the run initialization
    !! Loads the run input from the project file
    !! Initializes parameters and states
    !! Calls InitializeSimulationRunPart1
    integer(int8), intent(in) :: NrRun
    integer(intEnum), intent(in) :: TheProjectType

    type(rep_sum) :: SumWaBal_temp, PreviousSum_temp

    if (TheProjectType == typeproject_typenone) then
        ! Do nothing
        return
    end if

    call LoadSimulationRunProject(int(NrRun, kind=int32))

    call AdjustCompartments()
    SumWaBal_temp = GetSumWaBal()
    call GlobalZero(SumWaBal_temp)
    call SetSumWaBal(SumWaBal_temp)
    PreviousSum_temp = GetPreviousSum()
    call ResetPreviousSum(PreviousSum_temp)
    call SetPreviousSum(PreviousSum_temp)
    call InitializeSimulationRunPart1()

    contains


    subroutine AdjustCompartments()

        real(dp) :: TotDepth
        integer(int32) :: i
        type(CompartmentIndividual), &
                dimension(max_No_compartments) :: Comp_temp

        ! Adjust size of compartments if required
        TotDepth = 0._dp
        do i = 1, GetNrCompartments()
            TotDepth = TotDepth + GetCompartment_Thickness(i)
        end do
        if (GetSimulation_MultipleRunWithKeepSWC()) then
            ! Project with a sequence of simulation runs and KeepSWC
            if (roundc(GetSimulation_MultipleRunConstZrx()*1000._dp, mold=1) &
                > roundc(TotDepth*1000._dp, mold=1)) then
                call AdjustSizeCompartments(GetSimulation_MultipleRunConstZrx())
            end if
        else
            if (roundc(GetCrop_RootMax()*1000._dp, mold=1) &
                > roundc(TotDepth*1000._dp, mold=1)) then
                if (roundc(GetSoil_RootMax()*1000._dp, mold=1) &
                    == roundc(GetCrop_RootMax()*1000._dp, mold=1)) then
                    ! no restrictive soil layer
                    call AdjustSizeCompartments(GetCrop_RootMax())
                    ! adjust soil water content
                    Comp_temp = GetCompartment()
                    call CalculateAdjustedFC(GetZiAqua()/100._dp, Comp_temp)
                    call SetCompartment(Comp_temp)
                    if (GetSimulation_IniSWC_AtFC()) then
                        call ResetSWCToFC()
                    end if
                else
                    ! restrictive soil layer
                    if (roundc(GetSoil_RootMax()*1000._dp, mold=1) &
                        > roundc(TotDepth*1000._dp, mold=1)) then
                        call AdjustSizeCompartments(real(GetSoil_RootMax(), &
                                                            kind=dp))
                        ! adjust soil water content
                        Comp_temp = GetCompartment()
                        call CalculateAdjustedFC(GetZiAqua()/100._dp, Comp_temp)
                        call SetCompartment(Comp_temp)
                        if (GetSimulation_IniSWC_AtFC()) then
                            call ResetSWCToFC()
                        end if
                    end if
                end if
            end if
        end if
    end subroutine AdjustCompartments
end subroutine InitializeRunPart1


subroutine InitializeRunPart2(NrRun, TheProjectType)
    !! Part2 (after reading the climate) of the run initialization
    !! Calls InitializeSimulationRunPart2
    !! Initializes write out for the run
    integer(int8), intent(in) :: NrRun
    integer(intEnum), intent(in) :: TheProjectType

    call InitializeSimulationRunPart2()

    if (GetOutDaily()) then
        call WriteTitleDailyResults(TheProjectType, NrRun)
    end if

    if (GetPart1Mult()) then
        call WriteTitlePart1MultResults(TheProjectType, NrRun)
    end if

    if (GetPart2Eval() .and. (GetObservationsFile() /= '(None)')) then
        call CreateEvalData(NrRun)
    end if
end subroutine InitializeRunPart2


!--------duplicate nested in AdvanceOneTimeStep and FinalizeRun1----------!

subroutine RecordHarvest(NrCut, DayInSeason)
    integer(int32), intent(in) :: DayInSeason
    integer(int32), intent(in) :: NrCut

    integer(int32) :: Dayi, Monthi, Yeari
    logical :: NoYear
    character(len=1024) :: tempstring

    call fHarvest_open(GetfHarvest_filename(), 'a')
    call DetermineDate(GetCrop_Day1(), Dayi, Monthi, Yeari)
    NoYear = (Yeari == 1901)
    call DetermineDate(GetDayNri(), Dayi, Monthi, Yeari)
    if (NoYear) then
        Yeari = 9999
    end if

    if (NrCut == 9999) then
        ! last line at end of season
        write(tempstring, '(4i6, f34.3)') NrCut, Dayi, Monthi, Yeari, &
                                          GetSumWaBal_Biomass()
        call fHarvest_write(trim(tempstring), .false.)
        if (GetCrop_DryMatter() == undef_int) then
            write(tempstring, '(f20.3)') GetSumWaBal_YieldPart()
            call fHarvest_write(trim(tempstring))
        else
            write(tempstring, '(2f20.3)') GetSumWaBal_YieldPart(), &
                (GetSumWaBal_YieldPart()/(GetCrop_DryMatter()/100._dp))
            call fHarvest_write(trim(tempstring))
        end if
    else
        write(tempstring, '(6i6, f12.3, 2f10.3)') NrCut, Dayi, Monthi, Yeari, &
                DayInSeason, GetSumInterval(), (GetSumWaBal_Biomass()-GetBprevSum()), &
                GetSumWaBal_Biomass(), (GetSumWaBal_YieldPart()-GetYprevSum())
        call fHarvest_write(trim(tempstring), .false.)
        if (GetCrop_DryMatter() == undef_int) then
            write(tempstring, '(f10.3)') GetSumWaBal_YieldPart()
            call fHarvest_write(trim(tempstring))
        else
            write(tempstring, '(3f10.3)') GetSumWaBal_YieldPart(), &
                ((GetSumWaBal_YieldPart()-GetYprevSum()) &
                        /(GetCrop_DryMatter()/100._dp)), &
                (GetSumWaBal_YieldPart()/(GetCrop_DryMatter()/100._dp))
            call fHarvest_write(trim(tempstring))
        end if
    end if
end subroutine RecordHarvest

!--------end duplicate--------!


subroutine AdvanceOneTimeStep(WPi, HarvestNow)
    real(dp), intent(inout) :: WPi
    logical, intent(inout) :: HarvestNow

    real(dp) :: PotValSF, KsTr, TESTVALY, PreIrri, StressStomata, FracAssim
    integer(int32) :: VirtualTimeCC, DayInSeason
    real(dp) :: SumGDDadjCC, RatDGDD, &
                Biomass_temp, BiomassPot_temp, BiomassUnlim_temp, &
                BiomassTot_temp, YieldPart_temp, &
                ECe_temp, ECsw_temp, ECswFC_temp, KsSalt_temp
    type(rep_GwTable) :: GwTable_temp
    logical :: Store_temp, Mobilize_temp
    real(dp) :: ToMobilize_temp, Bmobilized_temp
    type(rep_EffectStress) :: EffectStress_temp
    logical :: SWCtopSOilConsidered_temp
    integer(int32) :: ZiAqua_temp
    real(dp) :: ECiAqua_temp, TactWeedInfested_temp,&
                Bin_temp, Bout_temp
    integer(int32) :: TargetTimeVal, TargetDepthVal
    integer(int8) :: PreviousStressLevel_temp, StressSFadjNEW_temp
    real(dp) :: CCxWitheredTpotNoS_temp, &
                StressLeaf_temp, StressSenescence_temp, TimeSenescence_temp, &
                SumKcTopStress_temp, SumKci_temp, WeedRCi_temp, &
                CCiActualWeedInfested_temp, HItimesBEF_temp, &
                ScorAT1_temp, ScorAT2_temp, HItimesAT1_temp, &
                HItimesAT2_temp, HItimesAT_temp, alfaHI_temp, &
                alfaHIAdj_temp, TESTVAL
    logical :: WaterTableInProfile_temp, NoMoreCrop_temp

    ! 1. Get ETo
    if (GetEToFile() == '(None)') then
        call SetETo(5.0_dp)
    end if

    ! 2. Get Rain
    if (GetRainFile() == '(None)') then
        call SetRain(0._dp)
    end if

    ! 3. Start mode
    if (GetStartMode()) then
        call SetStartMode(.false.)
    end if

    ! 4. Get depth and quality of the groundwater
    if (.not. GetSimulParam_ConstGwt()) then
        if (GetDayNri() > GetGwTable_DNr2()) then
            GwTable_temp = GetGwTable()
            call GetGwtSet(GetDayNri(), GwTable_temp)
            call SetGwTable(GwTable_temp)
        end if
        ZiAqua_temp = GetZiAqua()
        ECiAqua_temp = GetECiAqua()
        call GetZandECgwt(ZiAqua_temp, ECiAqua_temp)
        call SetZiAqua(ZiAqua_temp)
        call SetECiAqua(ECiAqua_temp)
        WaterTableInProfile_temp = GetWaterTableInProfile()
        call CheckForWaterTableInProfile((GetZiAqua()/100._dp), &
                      GetCompartment(), WaterTableInProfile_temp)
        call SetWaterTableInProfile(WaterTableInProfile_temp)
        if (GetWaterTableInProfile()) then
            call AdjustForWatertable()
        end if
    end if

    ! 5. Get Irrigation
    call SetIrrigation(0._dp)
    call GetIrriParam(TargetTimeVal, TargetDepthVal)

    ! 6. get virtual time for CC development
    SumGDDadjCC = real(undef_int, kind=dp)
    if (GetCrop_DaysToCCini() /= 0) then
        ! regrowth
        if (GetDayNri() >= GetCrop_Day1()) then
            ! time setting for canopy development
            VirtualTimeCC = (GetDayNri() - GetSimulation_DelayedDays() &
                             - GetCrop_Day1()) &
                             + GetTadj() + GetCrop_DaysToGermination()
            ! adjusted time scale
            if (VirtualTimeCC > GetCrop_DaysToHarvest()) then
                VirtualTimeCC = GetCrop_DaysToHarvest()
                ! special case where L123 > L1234
            end if
            if (VirtualTimeCC > GetCrop_DaysToFullCanopy()) then
                if ((GetDayNri() - GetSimulation_DelayedDays() - &
                     GetCrop_Day1()) <= GetCrop_DaysToSenescence()) then
                    VirtualTimeCC = GetCrop_DaysToFullCanopy() + &
                       roundc(GetDayFraction() * ((GetDayNri() - &
                              GetSimulation_DelayedDays() - &
                              GetCrop_Day1())+GetTadj()+ &
                              GetCrop_DaysToGermination() - &
                              GetCrop_DaysToFullCanopy()), mold=1) ! slow down
                else
                    VirtualTimeCC = GetDayNri() - &
                       GetSimulation_DelayedDays() - GetCrop_Day1() ! switch time scale
                end if
            end if
            if (GetCrop_ModeCycle() == modeCycle_GDDays) then
                SumGDDadjCC = GetSimulation_SumGDDfromDay1() + GetGDDTadj() + &
                              GetCrop_GDDaysToGermination()
                if (SumGDDadjCC > GetCrop_GDDaysToHarvest()) then
                    SumGDDadjCC = GetCrop_GDDaysToHarvest()
                    ! special case where L123 > L1234
                end if
                if (SumGDDadjCC > GetCrop_GDDaysToFullCanopy()) then
                    if (GetSimulation_SumGDDfromDay1() <= &
                        GetCrop_GDDaysToSenescence()) then
                        SumGDDadjCC = GetCrop_GDDaysToFullCanopy() &
                           + roundc(GetGDDayFraction() &
                             * (GetSimulation_SumGDDfromDay1() &
                             + GetGDDTadj()+GetCrop_GDDaysToGermination() &
                             - GetCrop_GDDaysToFullCanopy()), mold=1) ! slow down
                    else
                        SumGDDadjCC = GetSimulation_SumGDDfromDay1()
                        ! switch time scale
                    end if
                endif
            end if
            ! CC initial (at the end of previous day) when simulation starts
            ! before regrowth,
            if ((GetDayNri() == GetCrop_Day1()) .and. &
                (GetDayNri() > GetSimulation_FromDayNr())) then
                RatDGDD = 1._dp
                if ((GetCrop_ModeCycle() == modeCycle_GDDays) .and. &
                    (GetCrop_GDDaysToFullCanopySF() < &
                     GetCrop_GDDaysToSenescence())) then
                    RatDGDD = (GetCrop_DaysToSenescence() - &
                      GetCrop_DaysToFullCanopySF())/ &
                      real((GetCrop_GDDaysToSenescence() - &
                      GetCrop_GDDaysToFullCanopySF()), kind=dp)
                end if
                EffectStress_temp = GetSimulation_EffectStress()
                call CropStressParametersSoilFertility(&
                        GetCrop_StressResponse(), &
                        GetStressSFadjNEW(), EffectStress_temp)
                call SetSimulation_EffectStress(EffectStress_temp)
                call SetCCiPrev(CCiniTotalFromTimeToCCini(&
                        GetCrop_DaysToCCini(), &
                        GetCrop_GDDaysToCCini(), &
                        GetCrop_DaysToGermination(), &
                        GetCrop_DaysToFullCanopy(), &
                        GetCrop_DaysToFullCanopySF(), &
                        GetCrop_DaysToSenescence(), &
                        GetCrop_DaysToHarvest(), &
                        GetCrop_GDDaysToGermination(), &
                        GetCrop_GDDaysToFullCanopy(), &
                        GetCrop_GDDaysToFullCanopySF(), &
                        GetCrop_GDDaysToSenescence(), &
                        GetCrop_GDDaysToHarvest(), GetCrop_CCo(), &
                        GetCrop_CCx(), GetCrop_CGC(), &
                        GetCrop_GDDCGC(), GetCrop_CDC(), &
                        GetCrop_GDDCDC(), RatDGDD, &
                        GetSimulation_EffectStress_RedCGC(), &
                        GetSimulation_EffectStress_RedCCX(), &
                        GetSimulation_EffectStress_CDecline(), &
                        (GetCCxTotal()/GetCrop_CCx()), GetCrop_ModeCycle()))
                ! (CCxTotal/Crop.CCx) = fWeed
            end if
        else
            ! before start crop
            VirtualTimeCC = GetDayNri() - GetSimulation_DelayedDays() - &
                            GetCrop_Day1()
            if (GetCrop_ModeCycle() == modeCycle_GDDays) then
                SumGDDadjCC = GetSimulation_SumGDD()
            end if
        end if
    else
        ! sown or transplanted
        VirtualTimeCC = GetDayNri() - GetSimulation_DelayedDays() - &
                        GetCrop_Day1()
        if (GetCrop_ModeCycle() == modeCycle_GDDays) then
            SumGDDadjCC = GetSimulation_SumGDD()
        end if
        ! CC initial (at the end of previous day) when simulation starts
        ! before sowing/transplanting,
        if ((GetDayNri() == (GetCrop_Day1() + &
                             GetCrop_DaysToGermination())) &
            .and. (GetDayNri() > GetSimulation_FromDayNr())) then
            call SetCCiPrev(GetCCoTotal())
        end if
    end if

    ! 7. Rooting depth AND Inet day 1
    if (((GetCrop_ModeCycle() == modeCycle_CalendarDays) .and. &
        ((GetDayNri()-GetCrop_Day1()+1) < GetCrop_DaysToHarvest())) &
        .or. ((GetCrop_ModeCycle() == modeCycle_GDDays) &
          .and. (GetSimulation_SumGDD() < GetCrop_GDDaysToHarvest()))) then
        if (((GetDayNri()-GetSimulation_DelayedDays()) >= GetCrop_Day1()) .and. &
            ((GetDayNri()-GetSimulation_DelayedDays()) <= GetCrop_DayN())) then
            ! rooting depth at DAP (at Crop.Day1, DAP = 1)
            call SetRootingDepth(AdjustedRootingDepth(GetPlotVarCrop_ActVal(), &
                   GetPlotVarCrop_PotVal(), GetTpot(), GetTact(), &
                   GetStressLeaf(), GetStressSenescence(), &
                   (GetDayNri()-GetCrop_Day1()+1), &
                   GetCrop_DaysToGermination(), &
                   GetCrop_DaysToMaxRooting(), GetCrop_DaysToHarvest(), &
                   GetCrop_GDDaysToGermination(), &
                   GetCrop_GDDaysToMaxRooting(), &
                   GetCrop_GDDaysToHarvest(), GetSumGDDPrev(), &
                   (GetSimulation_SumGDD()), GetCrop_RootMin(), &
                   GetCrop_RootMax(), GetZiprev(), GetCrop_RootShape(), &
                   GetCrop_ModeCycle()))
            call SetZiprev(GetRootingDepth())
            ! IN CASE rootzone drops below groundwate table
            if ((GetZiAqua() >= 0._dp) .and. (GetRootingDepth() > &
                (GetZiAqua()/100._dp)) .and. (GetCrop_AnaeroPoint() > 0)) then
                call SetRootingDepth(GetZiAqua()/100._dp)
                if (GetRootingDepth() < GetCrop_RootMin()) then
                    call SetRootingDepth(GetCrop_RootMin())
                end if
             end if
        else
           call SetRootingDepth(0._dp)
        end if
    else
        call SetRootingDepth(GetZiprev())
    end if
    if ((GetRootingDepth() > 0._dp) .and. (GetDayNri() == GetCrop_Day1())) then
        ! initial root zone depletion day1 (for WRITE Output)
        SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
        call DetermineRootZoneWC(GetRootingDepth(), SWCtopSoilConsidered_temp)
        call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
        if (GetIrriMode() == IrriMode_Inet) then
           call AdjustSWCRootZone(PreIrri)  ! required to start germination
        end if
    end if

    ! 8. Transfer of Assimilates
    ToMobilize_temp = GetTransfer_ToMobilize()
    Bmobilized_temp = GetTransfer_Bmobilized()
    Store_temp = GetTransfer_Store()
    Mobilize_temp = GetTransfer_Mobilize()
    Bin_temp = GetBin()
    Bout_temp = GetBout()
    call InitializeTransferAssimilates(Bin_temp, Bout_temp, &
          ToMobilize_temp, Bmobilized_temp, FracAssim, Store_temp, &
          Mobilize_temp, HarvestNow)
    call SetTransfer_ToMobilize(ToMobilize_temp)
    call SetTransfer_Bmobilized(Bmobilized_temp)
    call SetTransfer_Store(Store_temp)
    call SetTransfer_Mobilize(Mobilize_temp)
    call SetBin(Bin_temp)
    call SetBout(Bout_temp)

    ! 9. RUN Soil water balance and actual Canopy Cover
    StressLeaf_temp = GetStressLeaf()
    StressSenescence_temp = GetStressSenescence()
    TimeSenescence_temp = GetTimeSenescence()
    NoMoreCrop_temp = GetNoMoreCrop()
    call BUDGET_module(GetDayNri(), TargetTimeVal, TargetDepthVal, &
           VirtualTimeCC, GetSumInterval(), GetDayLastCut(), &
           GetStressTot_NrD(), GetTadj(), GetGDDTadj(), GetGDDayi(), &
           GetCGCref(), GetGDDCGCref(), &
           GetCO2i(), GetCCxTotal(), GetCCoTotal(), GetCDCTotal(), &
           GetGDDCDCTotal(), SumGDDadjCC, &
           GetCoeffb0Salt(), GetCoeffb1Salt(), GetCoeffb2Salt(), &
           GetStressTot_Salt(), &
           GetDayFraction(), GetGDDayFraction(), FracAssim, &
           GetStressSFadjNEW(), GetTransfer_Store(), GetTransfer_Mobilize(), &
           StressLeaf_temp, StressSenescence_temp, TimeSenescence_temp, &
           NoMoreCrop_temp, TESTVAL)
    call SetStressLeaf(StressLeaf_temp)
    call SetStressSenescence(StressSenescence_temp)
    call SetTimeSenescence(TimeSenescence_temp)
    call SetNoMoreCrop(NoMoreCrop_temp)

    ! consider Pre-irrigation (6.) if IrriMode = Inet
    if ((GetRootingDepth() > 0._dp) .and. (GetDayNri() == GetCrop_Day1()) &
        .and. (GetIrriMode() == IrriMode_Inet)) then
         call SetIrrigation(GetIrrigation() + PreIrri)
         call SetSumWabal_Irrigation(GetSumWaBal_Irrigation() + PreIrri)
         PreIrri = 0._dp
     end if

     ! total number of days in the season
     if (GetCCiActual() > 0._dp) then
         if (GetStressTot_NrD() < 0) then
            call SetStressTot_NrD(1)
          else
            call SetStressTot_NrD(GetStressTot_NrD() + 1)
           end if
     end if

    ! 10. Potential biomass
    BiomassUnlim_temp = GetSumWaBal_BiomassUnlim()
    CCxWitheredTpotNoS_temp = GetCCxWitheredTpotNoS()
    call DeterminePotentialBiomass(VirtualTimeCC, SumGDDadjCC, &
         GetCO2i(), GetGDDayi(), CCxWitheredTpotNoS_temp, BiomassUnlim_temp)
    call SetCCxWitheredTpotNoS(CCxWitheredTpotNoS_temp)
    call SetSumWaBal_BiomassUnlim(BiomassUnlim_temp)

    ! 11. Biomass and yield
    if ((GetRootingDepth() > 0._dp) .and. (GetNoMoreCrop() .eqv. .false.)) then
        SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
        call DetermineRootZoneWC(GetRootingDepth(), SWCtopSoilConsidered_temp)
        call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
        ! temperature stress affecting crop transpiration
        if (GetCCiActual() <= 0.0000001_dp) then
             KsTr = 1._dp
        else
             KsTr = KsTemperature(0._dp, GetCrop_GDtranspLow(), GetGDDayi())
        end if
        call SetStressTot_Temp(((GetStressTot_NrD() - 1._dp)*GetStressTot_Temp() + &
                     100._dp*(1._dp-KsTr))/real(GetStressTot_NrD(), kind=dp))
        ! soil salinity stress
         ECe_temp = GetRootZoneSalt_ECe()
         ECsw_temp = GetRootZoneSalt_ECsw()
         ECswFC_temp = GetRootZoneSalt_ECswFC()
         KsSalt_temp = GetRootZoneSalt_KsSalt()
         call DetermineRootZoneSaltContent(GetRootingDepth(), ECe_temp, &
                 ECsw_temp, ECswFC_temp, KsSalt_temp)
         call SetRootZoneSalt_ECe(ECe_temp)
         call SetRootZoneSalt_ECsw(ECsw_temp)
         call SetRootZoneSalt_ECswFC(ECswFC_temp)
         call SetRootZoneSalt_KsSalt(KsSalt_temp)
         call SetStressTot_Salt(((GetStressTot_NrD() - 1._dp)*GetStressTot_Salt()&
                + 100._dp*(1._dp-GetRootZoneSalt_KsSalt()))/&
                   real(GetStressTot_NrD(), kind=dp))
         ! Biomass and yield
         Store_temp = GetTransfer_Store()
         Mobilize_temp = GetTransfer_Mobilize()
         ToMobilize_temp = GetTransfer_ToMobilize()
         Bmobilized_temp = GetTransfer_Bmobilized()
         Biomass_temp = GetSumWaBal_Biomass()
         BiomassPot_temp = GetSumWaBal_BiomassPot()
         BiomassUnlim_temp = GetSumWaBal_BiomassUnlim()
         BiomassTot_temp = GetSumWaBal_BiomassTot()
         YieldPart_temp = GetSumWaBal_YieldPart()
         TactWeedInfested_temp = GetTactWeedInfested()
         PreviousStressLevel_temp = GetPreviousStressLevel()
         StressSFadjNEW_temp = GetStressSFadjNEW()
         CCxWitheredTpotNoS_temp = GetCCxWitheredTpotNoS()
         Bin_temp = GetBin()
         Bout_temp = GetBout()
         SumKcTopStress_temp = GetSumKcTopStress()
         SumKci_temp = GetSumKci()
         WeedRCi_temp = GetWeedRCi()
         CCiActualWeedInfested_temp = GetCCiActualWeedInfested()
         HItimesBEF_temp = GetHItimesBEF()
         ScorAT1_temp = GetScorAT1()
         ScorAT2_temp = GetScorAT2()
         HItimesAT1_temp = GetHItimesAT1()
         HItimesAT2_temp = GetHItimesAT2()
         HItimesAT_temp = GetHItimesAT()
         alfaHI_temp = GetalfaHI()
         alfaHIAdj_temp = GetalfaHIAdj()
         call DetermineBiomassAndYield(GetDayNri(), GetETo(), GetTmin(), &
               GetTmax(), GetCO2i(), GetGDDayi(), GetTact(), GetSumKcTop(), &
               GetCGCref(), GetGDDCGCref(), GetCoeffb0(), GetCoeffb1(), &
               GetCoeffb2(), GetFracBiomassPotSF(), &
               GetCoeffb0Salt(), GetCoeffb1Salt(), GetCoeffb2Salt(), &
               GetStressTot_Salt(), SumGDDadjCC, GetCCiActual(), &
               FracAssim, VirtualTimeCC, GetSumInterval(), &
               Biomass_temp, BiomassPot_temp, BiomassUnlim_temp, &
               BiomassTot_temp, YieldPart_temp, WPi, &
               HItimesBEF_temp, ScorAT1_temp, ScorAT2_temp, &
               HItimesAT1_temp, HItimesAT2_temp, HItimesAT_temp, &
               alfaHI_temp, alfaHIAdj_temp, &
               SumKcTopStress_temp, SumKci_temp, &
               WeedRCi_temp, CCiActualWeedInfested_temp, &
               TactWeedInfested_temp, StressSFadjNEW_temp, &
               PreviousStressLevel_temp, Store_temp, &
               Mobilize_temp, ToMobilize_temp, &
               Bmobilized_temp, Bin_temp, Bout_temp, TESTVALY)
         call SetTransfer_Store(Store_temp)
         call SetTransfer_Mobilize(Mobilize_temp)
         call SetTransfer_ToMobilize(ToMobilize_temp)
         call SetTransfer_Bmobilized(Bmobilized_temp)
         call SetSumWaBal_Biomass(Biomass_temp)
         call SetSumWaBal_BiomassPot(BiomassPot_temp)
         call SetSumWaBal_BiomassUnlim(BiomassUnlim_temp)
         call SetSumWaBal_BiomassTot(BiomassTot_temp)
         call SetSumWaBal_YieldPart(YieldPart_temp)
         call SetTactWeedInfested(TactWeedInfested_temp)
         call SetBin(Bin_temp)
         call SetBout(Bout_temp)
         call SetPreviousStressLevel(int(PreviousStressLevel_temp, kind=int32))
         call SetStressSFadjNEW(int(StressSFadjNEW_temp, kind=int32))
         call SetCCxWitheredTpotNoS(CCxWitheredTpotNoS_temp)
         call SetSumKcTopStress(SumKcTopStress_temp)
         call SetSumKci(SumKci_temp)
         call SetWeedRCi(WeedRCi_temp)
         call SetCCiActualWeedInfested(CCiActualWeedInfested_temp)
         call SetHItimesBEF(HItimesBEF_temp)
         call SetScorAT1(ScorAT1_temp)
         call SetScorAT2(ScorAT2_temp)
         call SetHItimesAT1(HItimesAT1_temp)
         call SetHItimesAT2(HItimesAT2_temp)
         call SetHItimesAT(HItimesAT_temp)
         call SetalfaHI(alfaHI_temp)
         call SetalfaHIAdj(alfaHIAdj_temp)
    else
         !! SenStage = undef_int !GDL, 20220423, not used
         call SetWeedRCi(real(undef_int, kind=dp)) ! no crop and no weed infestation
         call SetCCiActualWeedInfested(0._dp) ! no crop
         call SetTactWeedInfested(0._dp) ! no crop
    end if

    ! 12. Reset after RUN
    if (GetPreDay() .eqv. .false.) then
        call SetPreviousDayNr(GetSimulation_FromDayNr() - 1)
    end if
    call SetPreDay(.true.)
    if (GetDayNri() >= GetCrop_Day1()) then
        call SetCCiPrev(GetCCiActual())
        if (GetZiprev() < GetRootingDepth()) then
            call SetZiprev(GetRootingDepth())
            ! IN CASE groundwater table does not affect root development
        end if
        call SetSumGDDPrev(GetSimulation_SumGDD())
    end if
    if (TargetTimeVal == 1) then
        call SetIrriInterval(0)
    end if

    ! 13. Cuttings
    if (GetManagement_Cuttings_Considered()) then
        HarvestNow = .false.
        DayInSeason = GetDayNri() - GetCrop_Day1() + 1
        call SetSumInterval(GetSumInterval() + 1)
        call SetSumGDDcuts( GetSumGDDcuts() + GetGDDayi())
        select case (GetManagement_Cuttings_Generate())
        case (.false.)
            if (GetManagement_Cuttings_FirstDayNr() /= undef_int) then
               ! adjust DayInSeason
                DayInSeason = GetDayNri() - &
                     GetManagement_Cuttings_FirstDayNr() + 1
            end if
            if ((DayInSeason >= GetCutInfoRecord1_FromDay()) .and. &
                (GetCutInfoRecord1_NoMoreInfo() .eqv. .false.)) then
                HarvestNow = .true.
                call GetNextHarvest()
            end if
            if (GetManagement_Cuttings_FirstDayNr() /= undef_int) then
               ! reset DayInSeason
                DayInSeason = GetDayNri() - GetCrop_Day1() + 1
            end if
        case (.true.)
            if ((DayInSeason > GetCutInfoRecord1_ToDay()) .and. &
                (GetCutInfoRecord1_NoMoreInfo() .eqv. .false.)) then
                call GetNextHarvest()
            end if
            select case (GetManagement_Cuttings_Criterion())
            case (TimeCuttings_IntDay)
                if ((GetSumInterval() >= GetCutInfoRecord1_IntervalInfo()) &
                     .and. (DayInSeason >= GetCutInfoRecord1_FromDay()) &
                     .and. (DayInSeason <= GetCutInfoRecord1_ToDay())) then
                    HarvestNow = .true.
                end if
            case (TimeCuttings_IntGDD)
                if ((GetSumGDDcuts() >= GetCutInfoRecord1_IntervalGDD()) &
                     .and. (DayInSeason >= GetCutInfoRecord1_FromDay()) &
                     .and. (DayInSeason <= GetCutInfoRecord1_ToDay())) then
                    HarvestNow = .true.
                end if
            case (TimeCuttings_DryB)
                if (((GetSumWabal_Biomass() - GetBprevSum()) >= &
                      GetCutInfoRecord1_MassInfo()) &
                     .and. (DayInSeason >= GetCutInfoRecord1_FromDay()) &
                     .and. (DayInSeason <= GetCutInfoRecord1_ToDay())) then
                    HarvestNow = .true.
                end if
            case (TimeCuttings_DryY)
                if (((GetSumWabal_YieldPart() - GetYprevSum()) >= &
                      GetCutInfoRecord1_MassInfo()) &
                    .and. (DayInSeason >= GetCutInfoRecord1_FromDay()) &
                    .and. (DayInSeason <= GetCutInfoRecord1_ToDay())) then
                    HarvestNow = .true.
                end if
            case (TimeCuttings_FreshY)
                ! OK if Crop.DryMatter = undef_int (not specified) HarvestNow
                ! remains false
                if ((((GetSumWaBal_YieldPart() - GetYprevSum())/&
                    (GetCrop_DryMatter()/100._dp)) >= &
                     GetCutInfoRecord1_MassInfo()) &
                    .and. (DayInSeason >= GetCutInfoRecord1_FromDay()) &
                    .and. (DayInSeason <= GetCutInfoRecord1_ToDay())) then
                    HarvestNow = .true.
                end if
            end select
        end select
        if (HarvestNow .eqv. .true.) then
            call SetNrCut(GetNrCut() + 1)
            call SetDayLastCut(DayInSeason)
            if (GetCCiPrev() > (GetManagement_Cuttings_CCcut()/100._dp)) then
                call SetCCiPrev(GetManagement_Cuttings_CCcut()/100._dp)
                ! ook nog CCwithered
                call SetCrop_CCxWithered(0._dp)  ! or CCiPrev ??
                call SetCCxWitheredTpotNoS(0._dp)
                   ! for calculation Maximum Biomass unlimited soil fertility
                call SetCrop_CCxAdjusted(GetCCiPrev()) ! new
            end if
            ! Record harvest
            if (GetPart1Mult()) then
                call RecordHarvest(GetNrCut(), DayInSeason)
            end if
            ! Reset
            call SetSumInterval(0)
            call SetSumGDDcuts(0._dp)
            call SetBprevSum(GetSumWaBal_Biomass())
            call SetYprevSum(GetSumWaBal_YieldPart())
        end if
    end if

    ! 14. Write results
    ! 14.a Summation
    call SetSumETo( GetSumETo() + GetETo())
    call SetSumGDD( GetSumGDD() + GetGDDayi())
    ! 14.b Stress totals
    if (GetCCiActual() > 0._dp) then
        ! leaf expansion growth
        if (GetStressLeaf() > - 0.000001_dp) then
            call SetStressTot_Exp(((GetStressTot_NrD() - 1._dp)*GetStressTot_Exp() &
                     + GetStressLeaf())/real(GetStressTot_NrD(), kind=dp))
        end if
        ! stomatal closure
        if (GetTpot() > 0._dp) then
            StressStomata = 100._dp *(1._dp - GetTact()/GetTpot())
            if (StressStomata > - 0.000001_dp) then
                call SetStressTot_Sto(((GetStressTot_NrD() - 1._dp) &
                    * GetStressTot_Sto() + StressStomata) / &
                    real(GetStressTot_NrD(), kind=dp))
            end if
        end if
    end if
    ! weed stress
    if (GetWeedRCi() > - 0.000001_dp) then
        call SetStressTot_Weed(((GetStressTot_NrD() - 1._dp)*GetStressTot_Weed() &
             + GetWeedRCi())/real(GetStressTot_NrD(), kind=dp))
    end if
    ! 14.c Assign crop parameters
    call SetPlotVarCrop_ActVal(GetCCiActual()/&
             GetCCxCropWeedsNoSFstress() * 100._dp)
    call SetPlotVarCrop_PotVal(100._dp * (1._dp/GetCCxCropWeedsNoSFstress()) * &
           CanopyCoverNoStressSF((VirtualTimeCC+GetSimulation_DelayedDays() &
             + 1), GetCrop_DaysToGermination(), GetCrop_DaysToSenescence(), &
             GetCrop_DaysToHarvest(), GetCrop_GDDaysToGermination(), &
             GetCrop_GDDaysToSenescence(), GetCrop_GDDaysToHarvest(), &
             (GetfWeedNoS()*GetCrop_CCo()), (GetfWeedNoS()*GetCrop_CCx()), &
             GetCGCref(), &
             (GetCrop_CDC()*(GetfWeedNoS()*GetCrop_CCx() + 2.29_dp)/&
             (GetCrop_CCx() + 2.29_dp)),&
             GetGDDCGCref(), &
             (GetCrop_GDDCDC()*(GetfWeedNoS()*GetCrop_CCx() + 2.29_dp)/&
             (GetCrop_CCx() + 2.29_dp)), &
             SumGDDadjCC, GetCrop_ModeCycle(), 0_int8, 0_int8))
    if ((VirtualTimeCC+GetSimulation_DelayedDays() + 1) <= &
         GetCrop_DaysToFullCanopySF()) then
        ! not yet canopy decline with soil fertility stress
        PotValSF = 100._dp * (1._dp/GetCCxCropWeedsNoSFstress()) * &
           CanopyCoverNoStressSF((VirtualTimeCC + &
            GetSimulation_DelayedDays() + 1), &
            GetCrop_DaysToGermination(), &
            GetCrop_DaysToSenescence(), GetCrop_DaysToHarvest(), &
            GetCrop_GDDaysToGermination(), GetCrop_GDDaysToSenescence(), &
            GetCrop_GDDaysToHarvest(),  GetCCoTotal(), GetCCxTotal(), &
            GetCrop_CGC(), GetCDCTotal(), &
            GetCrop_GDDCGC(), GetGDDCDCTotal(), &
            SumGDDadjCC, GetCrop_ModeCycle(), &
            GetSimulation_EffectStress_RedCGC(), &
            GetSimulation_EffectStress_RedCCX())
    else
        call GetPotValSF((VirtualTimeCC+GetSimulation_DelayedDays() + 1), &
               SumGDDAdjCC, PotValSF)
    end if
    ! 14.d Print ---------------------------------------
    if (GetOutputAggregate() > 0) then
        call CheckForPrint(GetTheProjectFile())
    end if
    if (GetOutDaily()) then
        call WriteDailyResults((GetDayNri() - GetSimulation_DelayedDays() &
                                - GetCrop_Day1()+1), WPi)
    end if
    if (GetPart2Eval() .and. (GetObservationsFile() /= '(None)')) then
        call WriteEvaluationData((GetDayNri()-GetSimulation_DelayedDays()-GetCrop_Day1()+1))
    end if

    ! 15. Prepare Next day
    ! 15.a Date
    call SetDayNri(GetDayNri() + 1)
    ! 15.b Irrigation
    if (GetDayNri() == GetCrop_Day1()) then
        call SetIrriInterval(1)
    else
        call SetIrriInterval(GetIrriInterval() + 1)
    end if
    ! 15.c Rooting depth
    ! 15.bis extra line for standalone
    if (GetOutDaily()) then
        call DetermineGrowthStage(GetDayNri(), GetCCiPrev())
    end if
    ! 15.extra - reset ageing of Kc at recovery after full senescence
    if (GetSimulation_SumEToStress() >= 0.1_dp) then
       call SetDayLastCut(GetDayNri())
    end if
end subroutine AdvanceOneTimeStep


subroutine ReadClimateNextDay()

    real(dp) :: ETo_tmp
    real(dp) :: tmpRain, Tmin_temp, Tmax_temp
    character(len=:), allocatable :: TempString

    ! Read Climate next day, Get GDDays and update SumGDDays
    if (GetDayNri() <= GetSimulation_ToDayNr()) then
        if (GetEToFile() /= '(None)') then
            TempString = fEToSIM_read()
            read(TempString,*) ETo_tmp
            call SetETo(ETo_tmp)
        end if
        if (GetRainFile() /= '(None)') then
            TempString = fRainSIM_read()
            read(TempString, *) tmpRain
            call SetRain(tmpRain)
        end if
        if (GetTemperatureFile() == '(None)') then
            call SetTmin(GetSimulParam_Tmin())
            call SetTmax(GetSimulParam_Tmax())
        else
            TempString = fTempSIM_read()
            read(TempString, *) Tmin_temp, Tmax_temp
            call SetTmin(Tmin_temp)
            call SetTmax(Tmax_temp)
        end if
    end if
end subroutine ReadClimateNextDay


subroutine SetGDDVariablesNextDay()

    if (GetDayNri() <= GetSimulation_ToDayNr()) then
        call SetGDDayi(DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), &
                GetTmin(), GetTmax(), GetSimulParam_GDDMethod()))
        if (GetDayNri() >= GetCrop_Day1()) then
            call SetSimulation_SumGDD(GetSimulation_SumGDD() + GetGDDayi())
            call SetSimulation_SumGDDfromDay1(GetSimulation_SumGDDfromDay1() &
                   + GetGDDayi())
        end if
    end if
end subroutine SetGDDVariablesNextDay


subroutine FinalizeRun1(NrRun, TheProjectFile, TheProjectType)
    integer(int8), intent(in) :: NrRun
    character(len=*), intent(in) :: TheProjectFile
    integer(intEnum), intent(in) :: TheProjectType

    ! 16. Finalise
    if ((GetDayNri()-1) == GetSimulation_ToDayNr()) then
        ! multiple cuttings
        if (GetPart1Mult()) then
            if (GetManagement_Cuttings_HarvestEnd() .eqv. .true.) then
                ! final harvest at crop maturity
                call SetNrCut(GetNrCut() + 1)
                call RecordHarvest(GetNrCut(), &
                                  (GetDayNri() - GetCrop_Day1()+1))
            end if
            call RecordHarvest((9999), &
                 (GetDayNri() - GetCrop_Day1()+1)) ! last line at end of season
        end if
        ! intermediate results
        if ((GetOutputAggregate() == 2) .or. (GetOutputAggregate() == 3) & ! 10-day and monthly results
            .and. ((GetDayNri()-1) > GetPreviousDayNr())) then
            call SetDayNri(GetDayNri()-1)
            call WriteIntermediatePeriod(TheProjectFile)
        end if
        call WriteSimPeriod(NrRun, TheProjectFile)
    end if
end subroutine FinalizeRun1


subroutine WriteDailyResults(DAP, WPi)
    integer(int32), intent(in) :: DAP
    real(dp), intent(in) :: WPi

    real(dp), parameter :: NoValD = undef_double
    integer(int32), parameter :: NoValI = undef_int
    integer(int32) :: Di, Mi, Yi, StrExp, StrSto, StrSalt, &
                      StrTr, StrW, Brel, Nr, DAP_loc
    integer(int32) :: Ratio1, Ratio2, Ratio3
    real(dp) :: KsTr, HI, KcVal, WPy, SaltVal, WPi_loc, tempreal
    logical :: SWCtopSoilConsidered_temp
    character(len=1025) :: tempstring

    DAP_loc = DAP
    WPi_loc = WPi

    call DetermineDate(GetDayNri(), Di, Mi, Yi)
    if (GetClimRecord_FromY() == 1901) then
        Yi = Yi - 1901 + 1
    end if
    if (GetStageCode() == 0) then
        DAP_loc = undef_int ! before or after cropping
    end if

    ! 0. info day
    write(tempstring, '(5i6)') Di, Mi, Yi, DAP_loc, GetStageCode()
    call fDaily_write(trim(tempstring), .false.)


    ! 1. Water balance
    if (GetOut1Wabal()) then
        if (GetZiAqua() == undef_int) then
            write(tempstring, '(f10.1, f8.1, f9.1, 3f7.1, 2f9.1, f8.2)') &
                    GetTotalWaterContent_EndDay(), GetRain(), GetIrrigation(), &
                    GetSurfaceStorage(), GetInfiltrated(), GetRunoff(), &
                    GetDrain(), GetCRwater(), undef_double
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(f10.1, f8.1, f9.1, 3f7.1, 2f9.1, f8.2)') &
                    GetTotalWaterContent_EndDay(), GetRain(), GetIrrigation(), &
                    GetSurfaceStorage(), GetInfiltrated(), GetRunoff(), &
                    GetDrain(), GetCRwater(), (GetZiAqua()/100._dp)
            call fDaily_write(trim(tempstring), .false.)
        end if
        if (GetTpot() > 0._dp) then
            Ratio1 = roundc(100._dp * GetTact()/GetTpot(), mold=1)
        else
            Ratio1 = 100
        end if

        if ((GetEpot()+GetTpot()) > 0._dp) then
            Ratio2 = roundc(100._dp * (GetEact()+GetTact())/(GetEpot()+GetTpot()), mold=1)
        else
            Ratio2 = 100
        end if

        if (GetEpot() > 0._dp) then
            Ratio3 = roundc(100._dp * GetEact()/GetEpot(), mold=1)
        else
            Ratio3 = 100
        end if

        if ((GetOut2Crop()) .or. (GetOut3Prof()) .or. (GetOut4Salt()) &
            .or. (GetOut5CompWC()) .or. (GetOut6CompEC()) &
            .or. (GetOut7Clim())) then
            write(tempstring, '(2f9.1, i7, 2f9.1, i6, f9.1, f8.1, i8)') &
                    GetEpot(), GetEact(), Ratio3, GetTpot(), GetTact(), Ratio1, &
                    (GetEpot()+GetTpot()), (GetEact()+GetTact()), Ratio2
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(2f9.1, i7, 2f9.1, i6, f9.1, f8.1, i8)') &
                    GetEpot(), GetEact(), Ratio3, GetTpot(), GetTact(), Ratio1, &
                    (GetEpot()+GetTpot()), (GetEact()+GetTact()), Ratio2
            call fDaily_write(trim(tempstring))
        end if
    end if


    ! 2. Crop development and yield
    if (GetOut2Crop()) then
        ! 1. relative transpiration
        if (GetTpot() > 0._dp) then
            Ratio1 = roundc(100._dp * GetTact()/GetTpot(), mold=1)
        else
            Ratio1 = 100
        end if
        ! 2. Water stresses
        if (GetStressLeaf() < 0._dp) then
            StrExp = undef_int
        else
            StrExp = roundc(GetStressLeaf(), mold=1)
        end if
        if (GetTpot() < epsilon(0._dp)) then
            StrSto = undef_int
        else
            StrSto = roundc(100._dp * (1._dp - GetTact()/GetTpot()), mold=1)
        end if

        ! 3. Salinity stress
        if (GetRootZoneSalt_KsSalt() < 0._dp) then
            StrSalt = undef_int
        else
            StrSalt = roundc(100._dp * (1._dp - GetRootZoneSalt_KsSalt()), &
                                                                    mold=1)
        end if

        ! 4. Air temperature stress
        if (GetCCiActual() <= 0.0000001_dp) then
            KsTr = 1._dp
        else
            KsTr = KsTemperature(0._dp, GetCrop_GDtranspLow(), GetGDDayi())
        end if

        if (KsTr < 1._dp) then
            StrTr = roundc((1._dp-KsTr)*100._dp, mold=1)
        else
            StrTr = 0._dp
        end if

        ! 5. Relative cover of weeds
        if (GetCCiActual() <= 0.0000001_dp) then
            StrW = undef_int
        else
            StrW = roundc(GetWeedRCi(), mold=1)
        end if

        ! 6. WPi adjustemnt
        if (GetSumWaBal_Biomass() <= 0.000001_dp) then
            WPi_loc = 0._dp
        end if

        ! 7. Harvest Index
        if ((GetSumWaBal_Biomass() > 0._dp) &
            .and. (GetSumWaBal_YieldPart() > 0._dp)) then
            HI = 100._dp * (GetSumWaBal_YieldPart())/(GetSumWaBal_Biomass())
        else
            HI = undef_double
        end if

        ! 8. Relative Biomass
        if ((GetSumWaBal_Biomass() > 0._dp) &
            .and. (GetSumWaBal_BiomassUnlim() > 0._dp)) then
            Brel = roundc(100._dp * GetSumWaBal_Biomass() &
                          /GetSumWaBal_BiomassUnlim(), mold=1)
            if (Brel > 100._dp) then
                Brel = 100._dp
            end if
        else
            Brel = undef_int
        end if

        ! 9. Kc coefficient
        if ((GetETo() > 0._dp) .and. (GetTpot() > 0._dp) &
            .and. (StrTr < 100)) then
            KcVal = GetTpot()/(GetETo()*KsTr)
        else
            KcVal = undef_int
        end if

        ! 10. Water Use Efficiency yield
        if (((GetSumWaBal_Tact() > 0._dp) .or. (GetSumWaBal_ECropCycle() > 0._dp)) &
            .and. (GetSumWaBal_YieldPart() > 0._dp)) then
            WPy = (GetSumWaBal_YieldPart()*1000._dp) &
                    /((GetSumWaBal_Tact()+GetSumWaBal_ECropCycle())*10._dp)
        else
            WPy = 0.0_dp
        end if

        ! write
        write(tempstring, '(f9.1, f8.2, 3i7, 2i7, 2f8.1, i7, f9.2, ' //&
                           '3f9.1, i6, f8.1, f10.3, f8.1, f9.3)') &
                GetGDDayi(), GetRootingDepth(), &
                StrExp, StrSto, &
                roundc(GetStressSenescence(), mold=1), StrSalt, StrW, &
                (GetCCiActual()*100._dp), (GetCCiActualWeedInfested()*100._dp), &
                StrTr, KcVal, GetTpot(), GetTact(), GetTactWeedInfested(), &
                Ratio1, (100._dp*WPi_loc), GetSumWaBal_Biomass(), HI, &
                GetSumWaBal_YieldPart()
        call fDaily_write(trim(tempstring), .false.)
        ! Fresh yield

        if ((GetCrop_DryMatter() == undef_int) &
            .or. (GetCrop_DryMatter() < epsilon(0._dp))) then
            write(tempstring, '(f9.3)') undef_double
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(f9.3)') GetSumWaBal_YieldPart() &
                                         /(GetCrop_DryMatter()/100._dp)
            call fDaily_write(trim(tempstring), .false.)
        end if

        ! finalize
        if ((GetOut3Prof()) .or. (GetOut4Salt()) .or. (GetOut5CompWC()) &
            .or. (GetOut6CompEC()) .or. (GetOut7Clim())) then
            write(tempstring, '(i8, f12.2, 2f9.3)') Brel, WPy, &
                                                    GetBin(), GetBout()
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(i8, f12.2, 2f9.3)') Brel, WPy, &
                                                    GetBin(), GetBout()
            call fDaily_write(trim(tempstring))
        end if
    end if

    ! 3. Profile/Root zone - Soil water content
    if (GetOut3Prof()) then
        write(tempstring, '(f10.1)') GetTotalWaterContent_EndDay()
        call fDaily_write(trim(tempstring), .false.)
        if (GetRootingDepth() < epsilon(0._dp)) then
            call SetRootZoneWC_Actual(undef_double)
        else
            if (roundc(GetSoil_RootMax()*1000._dp, mold=1) &
                == roundc(GetCrop_RootMax()*1000._dp, mold=1)) then
                SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
                call DetermineRootZoneWC(GetCrop_RootMax(), &
                                         SWCtopSoilConsidered_temp)
                call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
            else
                SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
                call DetermineRootZoneWC(real(GetSoil_RootMax(), kind=dp), SWCtopSoilConsidered_temp)
                call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
            end if
        end if
        write(tempstring, '(f9.1, f8.2)') GetRootZoneWC_actual(), GetRootingDepth()
        call fDaily_write(trim(tempstring), .false.)
        if (GetRootingDepth() < epsilon(0._dp)) then
            call SetRootZoneWC_Actual(undef_double)
            call SetRootZoneWC_FC(undef_double)
            call SetRootZoneWC_WP(undef_double)
            call SetRootZoneWC_SAT(undef_double)
            call SetRootZoneWC_Thresh(undef_double)
            call SetRootZoneWC_Leaf(undef_double)
            call SetRootZoneWC_Sen(undef_double)
        else
            SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
            call DetermineRootZoneWC(GetRootingDepth(), SWCtopSoilConsidered_temp)
            call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
        end if

        write(tempstring, '(f8.1, 5f10.1)') GetRootZoneWC_actual(), &
                GetRootZoneWC_SAT(), GetRootZoneWC_FC(), GetRootZoneWC_Leaf(), &
                GetRootZoneWC_Thresh(), GetRootZoneWC_Sen()
        call fDaily_write(trim(tempstring), .false.)
        if ((GetOut4Salt()) .or. (GetOut5CompWC()) .or. (GetOut6CompEC()) &
            .or. (GetOut7Clim())) then
            write(tempstring, '(f10.1)') GetRootZoneWC_WP()
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(f10.1)') GetRootZoneWC_WP()
            call fDaily_write(tempstring)
        end if
    end if

    ! 4. Profile/Root zone - soil salinity
    if (GetOut4Salt()) then
        write(tempstring, '(f9.3, 3f10.3)') GetSaltInfiltr(), &
                (GetDrain()*GetECdrain()*Equiv/100._dp), &
                (GetCRsalt()/100._dp), GetTotalSaltContent_EndDay()
        call fDaily_write(trim(tempstring), .false.)
        if (GetRootingDepth() < epsilon(0._dp)) then
            SaltVal = undef_int
            call SetRootZoneSalt_ECe(real(undef_int, kind=dp))
            call SetRootZoneSalt_ECsw(real(undef_int, kind=dp))
            call SetRootZoneSalt_KsSalt(1._dp)
        else
            SaltVal = (GetRootZoneWC_SAT()*GetRootZoneSalt_ECe()*Equiv)/100._dp
        end if
        if (GetZiAqua() == undef_int) then
            write(tempstring, '(f10.3, f8.2, f9.2, f8.2, i7, f8.2)') &
                    SaltVal, GetRootingDepth(), GetRootZoneSalt_ECe(), &
                    GetRootZoneSalt_ECsw(), &
                    roundc(100._dp*(1._dp-GetRootZoneSalt_KsSalt()), mold=1), &
                    undef_double
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(f10.3, f8.2, f9.2, f8.2, i7, f8.2)') &
                    SaltVal, GetRootingDepth(), GetRootZoneSalt_ECe(), &
                    GetRootZoneSalt_ECsw(), &
                    roundc(100._dp*(1._dp-GetRootZoneSalt_KsSalt()), mold=1), &
                    (GetZiAqua()/100._dp)
            call fDaily_write(trim(tempstring), .false.)
        end if
        if ((GetOut5CompWC()) .or. (GetOut6CompEC()) .or. (GetOut7Clim())) then
            write(tempstring, '(f8.2)') GetECiAqua()
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(f8.2)') GetECiAqua()
            call fDaily_write(trim(tempstring))
        end if
    end if

    ! 5. Compartments - Soil water content
    if (GetOut5CompWC()) then
        write(tempstring, '(f11.1)') (GetCompartment_Theta(1)*100._dp)
        call fDaily_write(trim(tempstring), .false.)
        do Nr = 2, (GetNrCompartments()-1)
            write(tempstring, '(f11.1)') &
                    (GetCompartment_Theta(Nr)*100._dp)
            call fDaily_write(trim(tempstring), .false.)
        end do
        if ((GetOut6CompEC()) .or. (GetOut7Clim())) then
            write(tempstring, '(f11.1)') &
                    (GetCompartment_Theta(GetNrCompartments())*100._dp)
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(f11.1)') &
                    (GetCompartment_Theta(GetNrCompartments())*100._dp)
            call fDaily_write(tempstring)
        end if
    end if

    ! 6. Compartmens - Electrical conductivity of the saturated soil-paste extract
    if (GetOut6CompEC()) then
        SaltVal = ECeComp(GetCompartment_i(1))
        write(tempstring, '(f11.1)') SaltVal
        call fDaily_write(trim(tempstring), .false.)
        do Nr = 2, (GetNrCompartments()-1)
            SaltVal = ECeComp(GetCompartment_i(Nr))
            write(tempstring, '(f11.1)') SaltVal
            call fDaily_write(trim(tempstring), .false.)
        end do
        SaltVal = ECeComp(GetCompartment_i(GetNrCompartments()))
        if (GetOut7Clim()) then
            write(tempstring, '(f11.1)') SaltVal
            call fDaily_write(trim(tempstring), .false.)
        else
            write(tempstring, '(f11.1)') SaltVal
            call fDaily_write(trim(tempstring))
        end if
    end if

    ! 7. Climate input parameters
    if (GetOut7Clim()) then
        tempreal = (GetTmin() + GetTmax())/2._dp
        write(tempstring, '(f9.1, 4f10.1, f10.2)') GetRain(), GetETo(), &
              GetTmin(), tempreal, GetTmax(), GetCO2i()
        call fDaily_write(trim(tempstring))
    end if
end subroutine WriteDailyResults


subroutine FileManagement()

    integer(int32) :: RepeatToDay
    real(dp) :: WPi
    logical :: HarvestNow

    WPi = 0._dp
    HarvestNow = .false.
    RepeatToDay = GetSimulation_ToDayNr()

    loop: do
        call AdvanceOneTimeStep(WPi, HarvestNow)
        call ReadClimateNextDay()
        call SetGDDVariablesNextDay()
        if ((GetDayNri()-1) == RepeatToDay) exit loop
    end do loop
end subroutine FileManagement


subroutine RunSimulation(TheProjectFile_, TheProjectType)
    character(len=*), intent(in) :: TheProjectFile_
    integer(intEnum), intent(in) :: TheProjectType

    integer(int8) :: NrRun
    integer(int32) :: NrRuns

    call InitializeSimulation(TheProjectFile_, TheProjectType)

    select case (TheProjectType)
    case(typeproject_TypePRO)
        NrRuns = 1
    case(typeproject_TypePRM)
        NrRuns = GetSimulation_NrRuns()
    end select

    do NrRun = 1, NrRuns
        call InitializeRunPart1(NrRun, TheProjectType)
        call InitializeClimate()
        call InitializeRunPart2(NrRun, TheProjectType)
        call FileManagement()
        call FinalizeRun1(NrRun, GetTheProjectFile(), TheProjectType)
        call FinalizeRun2(NrRun, TheProjectType)
    end do

    call FinalizeSimulation()
end subroutine RunSimulation

end module ac_run
