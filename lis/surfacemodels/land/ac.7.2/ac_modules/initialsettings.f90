module ac_initialsettings

use ac_defaultcropsoil, only :  ResetDefaultSoil, &
                                ResetDefaultCrop
use ac_global, only:    SetSimulParam_PercRAW, &
                        SetTnxReferenceFile, &
                        SetTnxReferenceFileFull, &
                        SetTnxReferenceYear, &
                        SetNrCompartments, &
                        SetSimulParam_CompDefThick, &
                        SetSimulParam_CropDay1, &
                        SetSimulParam_Tbase, &
                        SetSimulParam_Tupper, &
                        SetSimulParam_IrriFwInSeason, &
                        SetSimulParam_IrriFwOffSeason, &
                        SetSimulParam_RunOffDepth, &
                        SetSimulParam_CNcorrection, &
                        SetSimulParam_SaltDiff, &
                        SetSimulParam_SaltSolub, &
                        SetSimulParam_RootNrDF, &
                        SetSimulParam_IniAbstract, &
                        SetSimulParam_EffectiveRain_Method, &
                        EffectiveRainMethod_USDA, &
                        SetSimulParam_EvapDeclineFactor, &
                        SetSimulParam_KcWetBare, &
                        SetSimulParam_PercCCxHifinal, &
                        SetSimulParam_RootPercentZmin, &
                        SetSimulParam_MaxRootZoneExpansion, &
                        SetSimulParam_KsShapeFactorRoot, &
                        SetSimulParam_TAWGermination, &
                        SetSimulParam_pAdjFAO, &
                        SetSimulParam_DelayLowOxygen, &
                        SetSimulParam_ExpFsen, &
                        SetSimulParam_Beta, &
                        SetSimulParam_ThicknessTopSWC, &
                        SetSimulParam_EvapZmax, &
                        SetSimulParam_Tmin, &
                        SetSimulParam_Tmax, &
                        SetSimulParam_GDDMethod, &
                        SetPreDay, &
                        SetIniPercTAW, &
                        GetNrCompartments, &
                        max_No_compartments, &
                        SetGroundwaterFile, &
                        SetGroundwaterFilefull, &
                        GetGroundwaterFile, &
                        SetGroundwaterDescription, &
                        SetZiAqua, &
                        undef_int, &
                        SetECiAqua, &
                        SetSimulParam_ConstGwt, &
                        SetProfFile, &
                        SetProfFilefull, &
                        SetCrop_RootMin, &
                        SetCrop_RootMax, &
                        LoadProfile, &
                        CompleteProfileDescription, &
                        SetSimulation_CCini, &
                        SetSimulation_Bini, &
                        SetSimulation_Zrini, &
                        SetCropFile, &
                        SetCropFilefull, &
                        SetCrop_CCo, &
                        SetCrop_CCini, &
                        SetSoil_RootMax, &
                        GetPathNameSimul, &
                        GetCropFile, &
                        GetProfFile, &
                        GetCrop_PlantingDens, &
                        GetCrop_SizeSeedling, &
                        GetCrop_SizePlant, &
                        RootMaxInSoilProfile, &
                        GetCrop_RootMax, &
                        GetSoil_NrSoilLayers, &
                        GetSoilLayer, &
                        SetCrop_Day1, &
                        GetSimulParam_CropDay1, &
                        CompleteCropDescription, &
                        SetSimulation_YearSeason, &
                        NoCropCalendar, &
                        SetManFile, &
                        SetManFilefull, &
                        GetManFile, &
                        NoManagement, &
                        SetTemperatureFile, &
                        SetTemperatureFilefull, &
                        GetTemperatureFile, &
                        GetSimulParam_Tmin, &
                        GetSimulParam_Tmax, &
                        SetTemperatureDescription, &
                        SetTemperatureRecord_DataType, &
                        datatype_daily, &
                        SetTemperatureRecord_NrObs, &
                        SetTemperatureRecord_FromString, &
                        SetTemperatureRecord_ToString, &
                        SetTemperatureRecord_FromY, &
                        SetEToFile, &
                        SetEToFilefull, &
                        GetEToFile, &
                        SetEToDescription, &
                        SetEToRecord_DataType, &
                        SetEToRecord_NrObs, &
                        SetEToRecord_FromString, &
                        SetEToRecord_ToString, &
                        SetEToRecord_FromY, &
                        SetRainFile, &
                        SetRainFilefull, &
                        GetRainFile, &
                        SetRainDescription, &
                        SetRainRecord_DataType, &
                        SetRainRecord_NrObs, &
                        SetRainRecord_FromString, &
                        SetRainRecord_ToString, &
                        SetRainRecord_FromY, &
                        SetCO2File, &
                        SetCO2Filefull, &
                        GetCO2File, &
                        GetCO2Description, &
                        GetCO2Filefull, &
                        GenerateCO2Description, &
                        SetCO2Description, &
                        SetClimateFile, &
                        SetClimateFileFull, &
                        GetClimateFile, &
                        SetClimateDescription, &
                        SetClimData, &
                        SetSimulation_LinkCropToSimPeriod, &
                        GetCrop_Day1, &
                        GetCrop_DayN, &
                        AdjustCropYearToClimFile, &
                        SetCrop_DayN, &
                        GetClimFile, &
                        GetClimRecord_FromY, &
                        GetCLimRecord_NrObs, &
                        AdjustClimRecordTo, &
                        AdjustSimPeriod, &
                        SetIrriFile, &
                        SetIrriFilefull, &
                        GetIrriFile, &
                        NoIrrigation, &
                        SetOffSeasonFile, &
                        SetOffSeasonFilefull, &
                        GetOffSeasonFile, &
                        NoManagementOffSeason, &
                        SetProjectFile, &
                        SetProjectFilefull, &
                        GetProjectFile, &
                        SetSimulation_MultipleRun, &
                        SetSimulation_NrRuns, &
                        SetSimulation_MultipleRunWithKeepSWC, &
                        SetSimulation_MultipleRunCOnstZrx, &
                        SetMultipleProjectFile, &
                        SetMultipleProjectFilefull, &
                        GetProjectFileFull, &
                        SetObservationsFile, &
                        SetObservationsFilefull, &
                        GetObservationsFile, &
                        SetObservationsDescription, &
                        SetOutputName, &
                        SetOnset_Criterion, &
                        Criterion_RainPeriod, &
                        AirTCriterion_CumulGDD, &
                        SetOnset_AirTCriterion, &
                        AdjustOnsetSearchPeriod, &
                        SetETo, &
                        SetRain, &
                        SetIrrigation, &
                        SetSurfaceStorage, &
                        SetECstorage, &
                        SetDaySubmerged, &
                        GetSumWabal, &
                        GlobalZero, &
                        SetSumWaBal, &
                        SetDrain, &
                        SetRunoff, &
                        SetInfiltrated, &
                        SetCRwater, &
                        SetCRsalt, &
                        SetSimulation_ResetIniSWC, &
                        SetSimulation_EvapLimitON, &
                        SetMaxPlotNew, &
                        SetMaxPlotTr, &
                        SetSimulation_InitialStep, &
                        SetSimulation_LengthCuttingInterval, &
                        rep_sum, &
                        SetMultipleProjectDescription, &
                        GetProjectDescription, &
                        SetProjectDescription, &
                        GetSimulParam_CompDefThick, &
                        GetProfFileFull, &
                        SetCompartment_Thickness, &
                        SetSimulParam_EffectiveRain_RootNrEvap, &
                        SetSimulParam_EffectiveRain_ShowersInDecade, &
                        SetSimulParam_EffectiveRain_PercentEffRain
use ac_kinds, only: dp, &
                    int8, &
                    int32
implicit none


contains


subroutine InitializeSettings(use_default_soil_file,use_default_crop_file)
    logical, intent(in) :: use_default_soil_file
    logical, intent(in) :: use_default_crop_file
        !! Whether to make use of a 'DEFAULT.sol' soil file.

    character(len=1025) :: TempString1, TempString2, CO2descr
    integer(int32) :: Nri
    integer(int32) :: Crop_Day1_temp
    integer(int32) :: Crop_DayN_temp
    type(rep_sum) :: SumWaBal_temp

    ! 1. Program settings
    ! Settings of Program parameters
    ! 1a. General.PAR

    call SetSimulParam_PercRAW(50) ! Threshold [% of RAW]
                                   ! for determination of Inet
    call SetNrCompartments(12) ! Number of soil compartments (maximum is 12)
                          ! (not a program parameter)
    call SetSimulParam_CompDefThick(0.10_dp) ! Default thickness of
                                             ! soil compartments [m]
    call SetSimulParam_CropDay1(81) ! DayNumber of first day cropping period
                                    ! (1..365)
    call SetSimulParam_Tbase(10.0_dp)  ! Default base temperature (degC) below
                                       ! which no crop development
    call SetSimulParam_Tupper(30.0_dp) ! Default upper temperature threshold
                                       ! for crop development
    call SetSimulParam_IrriFwInSeason(100_int8) ! Percentage of soil surface
                                        ! wetted by irrigation in crop season
    call SetSimulParam_IrriFwOffSeason(100_int8) ! Percentage of soil surface
                                            ! wetted by irrigation off-season

    ! 1b. Soil.PAR - 6 parameters

    call SetSimulParam_RunoffDepth(0.30_dp) ! considered depth (m) of
              ! soil profile for calculation of mean soil water content
    call SetSimulParam_CNcorrection(.true.)
    call SetSimulParam_SaltDiff(20_int8) ! salt diffusion factor (%)
    call SetSimulParam_SaltSolub(100_int8) ! salt solubility (g/liter)
    call SetSimulParam_RootNrDF(16_int8) ! shape factor capillary rise factor
    call SetSimulParam_IniAbstract(5_int8) ! fixed in Version 5.0 cannot be
        ! changed since linked with equations for CN AMCII and CN converions

    ! 1c. Rainfall.PAR - 4 parameters
    call SetSimulParam_EffectiveRain_Method(EffectiveRainMethod_USDA)
    call SetSimulParam_EffectiveRain_PercentEffRain(70_int8) ! IF Method is
                                                             ! Percentage
    call SetSimulParam_EffectiveRain_ShowersInDecade(2_int8) ! For estimation
                                                          ! of surface run-off
    call SetSimulParam_EffectiveRain_RootNrEvap(5_int8) ! For reduction of
                                                        ! soil evaporation

    ! 1d. Crop.PAR  - 12 parameters
    call SetSimulParam_EvapDeclineFactor(4_int8) ! evaporation decline
                                                 ! factor in stage 2
    call SetSimulParam_KcWetBare(1.10_dp) ! Kc wet bare soil [-]
    call SetSimulParam_PercCCxHIfinal(5_int8) ! CC threshold below which HI no
                                              ! longer increase(% of 100)
    call SetSimulParam_RootPercentZmin(70) ! Starting depth of root sine
                                           ! function (% of Zmin)
    call SetSimulParam_MaxRootZoneExpansion(5.00_dp) ! fixed at 5 cm/day
    call SetSimulParam_KsShapeFactorRoot(-6_int8) ! Shape factor for effect
                                        ! water stress on rootzone expansion
    call SetSimulParam_TAWGermination(20_int8) ! Soil water content (% TAW)
                                 ! required at sowing depth for germination
    call SetSimulParam_pAdjFAO(1._dp) ! Adjustment factor for FAO-adjustment
                                    ! soil water depletion (p) for various ET
    call SetSimulParam_DelayLowOxygen(3) ! number of days for full effect of
                                         ! deficient aeration
    call SetSimulParam_ExpFsen(1.00_dp) ! exponent of senescence factor
                ! adjusting drop in photosynthetic activity of dying crop
    call SetSimulParam_Beta(12_int8) ! Decrease (percentage) of p(senescence)
                                  ! once early canopy senescence is triggered
    call SetSimulParam_ThicknessTopSWC(10_int8) ! Thickness top soil (cm) in
                            ! which soil water depletion has to be determined

    ! 1e. Field.PAR - 1 parameter
    call SetSimulParam_EvapZmax(30_int8) ! maximum water extraction depth by
                                         ! soil evaporation [cm]

    ! 1f. Temperature.PAR - 3 parameters
    call SetSimulParam_Tmin(12.0_dp) ! Default minimum temperature (degC) if no
                                     ! temperature file is specified
    call SetSimulParam_Tmax(28.0_dp) ! Default maximum temperature (degC) if no
                                     ! temperature file is specified
    call SetSimulParam_GDDMethod(3_int8) ! Default method for GDD calculations



    call SetPreDay(.false.)
    call SetIniPercTAW(50_int8) ! Default Value for Percentage TAW for Display
                                ! in Initial Soil Water Content Menu

    ! Default for soil compartments
    if (GetNrCompartments() > max_No_compartments) then
        ! Savety check of value in General.PAR;
        call SetNrCompartments(max_No_compartments)
    end if
    do Nri = 1, max_No_compartments
        ! required for formactivate ParamNew
        call SetCompartment_Thickness(Nri, GetSimulParam_CompDefThick())
    end do
    ! Default CropDay1 - Savety check of value in General.PAR
    do while (GetSimulParam_CropDay1() > 365)
        call SetSimulParam_CropDay1(GetSimulParam_CropDay1()-365)
    end do
    if (GetSimulParam_CropDay1() < 1) then
        call SetSimulParam_CropDay1(1)
    end if

    ! 2a. Ground water table
    call SetGroundWaterFile('(None)')
    call SetGroundWaterFilefull(GetGroundWaterFile())  ! no file
    call SetGroundWaterDescription('no shallow groundwater table')
    call SetZiAqua(undef_int)
    call SetECiAqua(real(undef_int, kind=dp))
    call SetSimulParam_ConstGwt(.true.)

    ! 2b. Soil profile and initial soil water content
    call ResetDefaultSoil(use_default_soil_file) ! Reset the soil profile to its default values
    ! required for SetSoil_RootMax(RootMaxInSoilProfile(GetCrop().RootMax,
                                 ! GetCrop().RootMin,GetSoil().NrSoilLayers,
                                 ! SoilLayer)) in LoadProfile
    call SetCrop_RootMin(0.30_dp) ! Minimum rooting depth (m)
    call SetCrop_RootMax(1.00_dp) ! Maximum rooting depth (m)

    if (use_default_soil_file) then
        call SetProfFile('DEFAULT.SOL')
        call SetProfFilefull(GetPathNameSimul() // GetProfFile())
        ! Crop.RootMin, RootMax, and Soil.RootMax are
        ! correctly calculated in LoadCrop
        call LoadProfile(GetProfFilefull())
    end if

    call CompleteProfileDescription ! Simulation.ResetIniSWC AND
                        ! specify_soil_layer whcih contains
                        ! PROCEDURE DeclareInitialCondAtFCandNoSalt,
                        ! in which SWCiniFile := '(None)', and settings
                        ! for Soil water and Salinity content

    ! 2c. Complete initial conditions (crop development)
    call SetSimulation_CCini(real(undef_int, kind=dp))
    call SetSimulation_Bini(0.000_dp)
    call SetSimulation_Zrini(real(undef_int, kind=dp))


    ! 3. Crop characteristics and cropping period
    call ResetDefaultCrop(use_default_crop_file) ! Reset the crop to its default values
    call SetCropFile('DEFAULT.CRO')
    call SetCropFilefull(GetPathNameSimul() // GetCropFile())
    ! LoadCrop ==============================
    call SetCrop_CCo((GetCrop_PlantingDens()/10000._dp) &
                        * (GetCrop_SizeSeedling()/10000._dp))
    call SetCrop_CCini((GetCrop_PlantingDens()/10000._dp) &
                        * (GetCrop_SizePlant()/10000._dp))
    ! maximum rooting depth in given soil profile
    call SetSoil_RootMax(RootMaxInSoilProfile(GetCrop_RootMax(), &
                                              GetSoil_NrSoilLayers(), &
                                              GetSoilLayer()))
    ! determine miscellaneous
    call SetCrop_Day1(GetSimulParam_CropDay1())
    call CompleteCropDescription()
    call SetSimulation_YearSeason(1_int8)
    call NoCropCalendar()

    ! 4. Field Management
    call SetManFile('(None)')
    call SetManFilefull(GetManFile())  ! no file
    call NoManagement()

    ! 5. Climate
    ! 5.1 Temperature
    call SetTemperatureFile('(None)')
    call SetTemperatureFilefull(GetTemperatureFile())  ! no file
    call SetTnxReferenceFile(GetTemperatureFile()) ! no file
    call SetTnxReferenceYear(2000) ! for refernce CO2 concentration
    write(TempString1, '(f8.1)') GetSimulParam_Tmin()
    write(TempString2, '(f8.1)') GetSimulParam_Tmax()
    call SetTemperatureDescription('')
    call SetTemperatureRecord_DataType(datatype_Daily)
    call SetTemperatureRecord_NrObs(0)
    call SetTemperatureRecord_FromString('any date')
    call SetTemperatureRecord_ToString('any date')
    call SetTemperatureRecord_FromY(1901)

    ! 5.2 ETo
    call SetEToFile('(None)')
    call SetEToFilefull(GetEToFile())  ! no file
    call SetEToDescription('')
    call SetEToRecord_DataType(datatype_Daily)
    call SetEToRecord_NrObs(0)
    call SetEToRecord_FromString('any date')
    call SetEToRecord_ToString('any date')
    call SetEToRecord_FromY(1901)

    ! 5.3 Rain
    call SetRainFile('(None)')
    call SetRainFilefull(GetRainFile())  ! no file
    call SetRainDescription('')
    call SetRainRecord_DataType(datatype_Daily)
    call SetRainRecord_NrObs(0)
    call SetRainRecord_FromString('any date')
    call SetRainRecord_ToString('any date')
    call SetRainRecord_FromY(1901)

    ! 5.4 CO2
    call SetCO2File('MaunaLoa.CO2')
    call SetCO2FileFull(GetPathNameSimul() // GetCO2File())
    CO2descr = GetCO2Description()
    call GenerateCO2Description(GetCO2FileFull(), CO2descr)
    call SetCO2Description(CO2descr)

    ! 5.5 Climate file
    call SetClimateFile('(None)')
    call SetClimateFileFull(GetClimateFile())
    call SetClimateDescription('')

    ! 5.6 Set Climate and Simulation Period
    call SetClimData
    call SetSimulation_LinkCropToSimPeriod(.true.)
    ! adjusting Crop.Day1 and Crop.DayN to ClimFile
    Crop_Day1_temp = GetCrop_Day1()
    Crop_DayN_temp = GetCrop_DayN()
    call AdjustCropYearToClimFile(Crop_Day1_temp, Crop_DayN_temp)
    call SetCrop_Day1(Crop_Day1_temp)
    call SetCrop_DayN(Crop_DayN_temp)
    ! adjusting ClimRecord.'TO' for undefined year with 365 days
    if ((GetClimFile() /= '(None)') .and. (GetClimRecord_FromY() == 1901) &
        .and. (GetClimRecord_NrObs() == 365)) then
        call AdjustClimRecordTo(GetCrop_DayN())
    end if
    ! adjusting simulation period
    call AdjustSimPeriod()

    ! 6. irrigation
    call SetIrriFile('(None)')
    call SetIrriFilefull(GetIrriFile())  ! no file
    call NoIrrigation()

    ! 7. Off-season
    call SetOffSeasonFile('(None)')
    call SetOffSeasonFileFull(GetOffSeasonFile())
    call NoManagementOffSeason()

    ! 8. Project and Multiple Project file
    call SetProjectFile('(None)')
    call SetProjectFileFull(GetProjectFile())
    call SetProjectDescription('No specific project')
    call SetSimulation_MultipleRun(.false.) ! No sequence of simulation
                                            ! runs in the project
    call SetSimulation_NrRuns(1)
    call SetSimulation_MultipleRunWithKeepSWC(.false.)
    call SetSimulation_MultipleRunConstZrx(real(undef_int, kind=dp))
    call SetMultipleProjectFile(GetProjectFile())
    call SetMultipleProjectFileFull(GetProjectFileFull())
    call SetMultipleProjectDescription(GetProjectDescription())

    ! 9. Observations file
    call SetObservationsFile('(None)')
    call SetObservationsFileFull(GetObservationsFile())
    call SetObservationsDescription('No field observations')

    ! 10. Output files
    call SetOutputName('Project')

    ! 11. Onset
    call SetOnset_Criterion(Criterion_RainPeriod)
    call SetOnset_AirTCriterion(AirTCriterion_CumulGDD)
    call AdjustOnsetSearchPeriod()

    ! 12. Simulation run
    call SetETo(5.0_dp)
    call SetRain(0._dp)
    call SetIrrigation(0._dp)
    call SetSurfaceStorage(0._dp)
    call SetECstorage(0.0_dp)
    call SetDaySubmerged(0)
    SumWaBal_temp = GetSumWaBal()
    call GlobalZero(SumWaBal_temp)
    call SetSumWaBal(SumWaBal_temp)
    call SetDrain(0.0_dp) ! added 4.0
    call SetRunoff(0.0_dp)! added 4.0
    call SetInfiltrated(0.0_dp) ! added 4.0
    call SetCRwater(0._dp) ! added 4.0
    call SetCRsalt(0._dp) ! added 4.0
    call SetSimulation_ResetIniSWC(.true.)
    call SetSimulation_EvapLimitON(.false.)
    call SetMaxPlotNew(50)
    call SetMaxPlotTr(10_int8)
    call SetSimulation_InitialStep(10) ! Length of period (days) for displaying
                                    ! intermediate results during simulation run
    call SetSimulation_LengthCuttingInterval(40) ! Default length of
                                                 ! cutting interval (days)
end subroutine InitializeSettings

end module ac_initialsettings
