!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: AC72_main
! \label{AC72_main}
!
! !REVISION HISTORY:
!   04 NOV 2024, Louise Busschaert; initial implementation
!
! !INTERFACE:
subroutine AC72_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod
    use LIS_constantsMod
    use AC72_lsmMod
    use ac72_prep_f
    use ESMF

    !!! MB_AC70
    use ac_global, only:    DegreesDay,&
                            GetCCiActual,&
                            GetCCiprev,&
                            GetCCiTopEarlySen,&
                            GetClimRecord,&
                            GetCompartment,&
                            GetCompartment_theta,&
                            GetCrop,&
                            GetCrop_Day1,&
                            GetCrop_DayN,&
                            GetCrop_DaysToCCini,&
                            GetCrop_DaysToFlowering,&
                            GetCrop_DaysToFullCanopy,&
                            GetCrop_DaysToFullCanopySF,&
                            GetCrop_DaysToGermination,&
                            GetCrop_DaysToHarvest,&
                            GetCrop_DaysToHIo,&
                            GetCrop_DaysToMaxRooting,&
                            GetCrop_DaysTosenescence,&
                            GetCrop_GDDaysToFlowering,&
                            GetCrop_GDDaysToGermination,&
                            GetCrop_GDDaysToGermination,&
                            GetCrop_GDDaysToHarvest,&
                            GetCrop_GDDaysToMaxRooting,&
                            GetCrop_GDDaysToSenescence,&
                            GetCrop_ModeCycle,&
                            GetCRsalt,&  
                            GetCRwater,& 
                            GetDaySubmerged,&
                            GetDrain,&  
                            GetEact,& 
                            GetECdrain,& 
                            GetECiAqua,& 
                            GetECstorage,& 
                            Geteffectiverain,&
                            GetEndSeason,&
                            GetEpot,& 
                            GetETo,&
                            GetEToRecord,&
                            GetEvapoEntireSoilSurface,&
                            GetGenerateDepthMode,&
                            GetGenerateTimeMode,&
                            GetInfiltrated,&
                            GetIniPercTAW,&
                            GetIrriAfterSeason,&
                            GetIrriBeforeSeason,&
                            GetIrriECw,&
                            GetIrriFirstDayNr,&
                            GetIrrigation,&
                            GetIrriMethod,&
                            GetIrriMode,&
                            GetManagement,&
                            GetManagement_Cuttings,&
                            GetMaxPlotNew,&
                            GetMaxPlotTr,&
                            GetNrCompartments,&
                            GetOnset,&
                            GetOut1Wabal,&
                            GetOut2Crop,&
                            GetOut3Prof,&
                            GetOut4Salt,&
                            GetOut5CompWC,&
                            GetOut6CompEC,&
                            GetOut7Clim,&
                            GetOut8Irri,&
                            GetOutDaily,&
                            GetOutputAggregate,&
                            GetPart1Mult,&
                            GetPart2Eval,&
                            GetPerennialPeriod,&
                            GetPreDay,&
                            GetRain,& 
                            GetRainRecord,&
                            GetRootingDepth,&
                            GetRootZoneSalt,&
                            GetRootZoneWC_Actual,&
                            GetRootZoneWC_FC,&
                            GetRootZoneWC_Leaf,&
                            GetRootZoneWC_SAT,&
                            GetRootZoneWC_Sen,&
                            GetRootZoneWC_Thresh,&
                            GetRootZoneWC_WP,&
                            GetRootZoneWC_ZtopAct,&
                            GetRootZoneWC_ZtopFC,&
                            GetRootZoneWC_ZtopThresh,&
                            GetRootZoneWC_ZtopWP,&
                            GetRunoff,& 
                            GetSaltInfiltr,&
                            GetSimulation,&
                            GetSimulation_EvapLimitON, &
                            GetSimulation_FromDayNr, &
                            GetSimulation_SumGDD,&
                            GetSimulation_SumGDD,&
                            GetSimulation_SumGDDfromDay1,&
                            GetSimulation_SumGDDfromDay1,&
                            GetSimulation_SWCtopSoilConsidered, &
                            GetSimulation_ToDayNr,&
                            GetSimulation_ToDayNr,&
                            GetSimulParam,&
                            GetSimulParam_GDDMethod,&
                            GetSimulParam_ThicknessTopSWC,&
                            GetSoil,&
                            GetSoilLayer,&
                            GetSumWaBal,&
                            GetSurf0,& 
                            GetSurfaceStorage,&
                            GetTact,&
                            GetTactWeedInfested,&
                            GetTemperatureRecord,&
                            GetTmax,& 
                            GetTmaxRun, &
                            GetTmaxRun_i, &
                            GetTmin,& 
                            GetTminRun, &
                            GetTminRun_i, &
                            GetTotalSaltContent,&
                            GetTotalWaterContent,&
                            GetTpot,&
                            GetZiAqua,&
                            IrriMode_Generate,&
                            IrriMode_Inet,&
                            IrriMode_Manual,&
                            IrriMode_NoIrri,&
                            ModeCycle_GDDays,&
                            SetCCiActual,&
                            SetCCiprev,&
                            SetCCiTopEarlySen,&
                            SetClimFile, &
                            SetClimRecord,&
                            SetClimRecord_DataType, &
                            SetClimRecord_fromd, &
                            SetClimRecord_fromdaynr, &
                            SetClimRecord_fromm, &
                            SetClimRecord_fromstring, &
                            SetClimRecord_fromy, &
                            SetClimRecord_NrObs, &
                            SetClimRecord_tod, &
                            SetClimRecord_todaynr, &
                            SetClimRecord_tom, &
                            SetClimRecord_tostring, &
                            SetClimRecord_toy,&
                            SetCompartment,&
                            SetCompartment_theta,&
                            SetCrop,&
                            SetCrop_DayN,&
                            SetCrop_DaysToCCini,&
                            SetCrop_DaysToFlowering,&
                            SetCrop_DaysToFullCanopy,&
                            SetCrop_DaysToFullCanopySF,&
                            SetCrop_DaysToGermination,&
                            SetCrop_DaysToHarvest,&
                            SetCrop_DaysToHIo,&
                            SetCrop_DaysToMaxRooting,&
                            SetCrop_DaysTosenescence,&
                            SetCRsalt,&  
                            SetCRwater,& 
                            SetDaySubmerged,&
                            SetDrain,&  
                            SetEact,& 
                            SetECdrain,& 
                            SetECiAqua,& 
                            SetECstorage,& 
                            Seteffectiverain,&
                            SetEndSeason,&
                            SetEpot,& 
                            SetETo,&
                            SetEToRecord,&
                            SetEvapoEntireSoilSurface,&
                            SetGenerateDepthMode,&
                            SetGenerateTimeMode,&
                            SetInfiltrated,&
                            SetIniPercTAW,&
                            SetIrriAfterSeason,&
                            SetIrriBeforeSeason,&
                            SetIrriECw,&
                            SetIrriFirstDayNr,&
                            SetIrrigation,&
                            SetIrriMethod,&
                            SetIrriMode,&
                            SetManagement,&
                            SetManagement_Cuttings,&
                            SetMaxPlotNew,&
                            SetMaxPlotTr,&
                            SetNrCompartments,&
                            SetOnset,&
                            SetOut1Wabal,&
                            SetOut2Crop,&
                            SetOut3Prof,&
                            SetOut4Salt,&
                            SetOut5CompWC,&
                            SetOut6CompEC,&
                            SetOut7Clim,&
                            SetOut8Irri,&
                            SetOutDaily,&
                            SetOutputAggregate,&
                            SetPart1Mult,&
                            SetPart2Eval,&
                            SetPerennialPeriod,&
                            SetPreDay,&
                            SetRain,& 
                            SetRainRecord,&
                            SetRootingDepth,&
                            SetRootZoneSalt,&
                            SetRootZoneWC_Actual,&
                            SetRootZoneWC_FC,&
                            SetRootZoneWC_Leaf,&
                            SetRootZoneWC_SAT,&
                            SetRootZoneWC_Sen,&
                            SetRootZoneWC_Thresh,&
                            SetRootZoneWC_WP,&
                            SetRootZoneWC_ZtopAct,&
                            SetRootZoneWC_ZtopFC,&
                            SetRootZoneWC_ZtopThresh,&
                            SetRootZoneWC_ZtopWP,&
                            SetRunoff,& 
                            SetSaltInfiltr,&
                            SetSimulation,&
                            SetSimulation_SumGDD, &
                            SetSimulation_SumGDDfromDay1,&
                            SetSimulation_ToDayNr, &
                            SetSimulParam,&
                            SetSimulParam_ThicknessTopSWC, &
                            SetSoil,&
                            SetSoilLayer,&
                            SetSumWaBal,&
                            SetSurf0,& 
                            SetSurfaceStorage,&
                            SetTact,&
                            SetTactWeedInfested,&
                            SetTemperatureRecord,&
                            SetTmax,& 
                            SetTmaxRun, &
                            SetTmaxRun_i, &
                            SetTmaxTnxReference12MonthsRun, &
                            SetTmin,&
                            SetTminRun, &
                            SetTminRun_i, &
                            SetTminTnxReference12MonthsRun, &
                            SetTnxReferenceFile, &
                            SetTnxReferenceYear,&
                            SetTotalSaltContent,&
                            SetTotalWaterContent,&
                            SetTpot,&
                            SetZiAqua,&
                            typeproject_typeprm, &
                            typeproject_typepro, &
                            undef_int

    use ac_project_input, only: ProjectInput 

    use ac_run, only:   AdvanceOneTimeStep, &
                        FinalizeRun1, &
                        FinalizeRun2, &
                        FinalizeSimulation, &
                        fIrri_close, &
                        fIrri_open,&
                        firri_read,&
                        GetalfaHI,&
                        GetalfaHIAdj,&
                        GetBin,&
                        GetBout,&
                        GetBprevSum,&
                        GetCCiActualWeedInfested,&
                        GetCCoTotal,&
                        GetCCxCropWeedsNoSFstress,&
                        GetCCxTotal,&
                        GetCCxWitheredTpotNoS,&
                        GetCDCTotal,&
                        GetCGCref,&
                        GetCO2i,&
                        GetCoeffb0,&
                        GetCoeffb0Salt,&
                        GetCoeffb1,&
                        GetCoeffb1Salt,&
                        GetCoeffb2,&
                        GetCoeffb2Salt,&
                        GetCrop_Tbase, &
                        GetCrop_Tupper, &
                        GetCutInfoRecord1,&
                        GetCutInfoRecord2,&
                        GetDayFraction,&
                        GetDayLastCut,&
                        GetDayNr1Eval,&
                        GetDayNrEval,&
                        GetDayNri,&
                        GetFracBiomassPotSF,&
                        GetfWeedNoS,&
                        GetGDDayFraction,&
                        GetGDDayi,&
                        GetGDDCDCTotal,&
                        GetGDDCGCref ,&
                        GetGDDTadj,&
                        GetGwTable,&
                        GetHItimesAT,&
                        GetHItimesAT1,&
                        GetHItimesAT2,&
                        GetHItimesBEF,&
                        GetIrriInfoRecord1, &
                        GetIrriInfoRecord1,&
                        GetIrriInfoRecord1_NoMoreInfo,&
                        GetIrriInfoRecord2, &
                        GetIrriInfoRecord2,&
                        GetIrriInterval,&
                        GetLineNrEval,&
                        GetNextSimFromDayNr ,&
                        GetNoMoreCrop,&
                        GetNoYear,&
                        GetNrCut,&
                        GetPlotVarCrop,&
                        GetPreviousBmob,&
                        GetPreviousBsto,&
                        GetPreviousDayNr,&
                        GetPreviousStressLevel,&
                        GetPreviousSum,&
                        GetPreviousSumETo,&
                        GetPreviousSumGDD,&
                        GetScorAT1,&
                        GetScorAT2,&
                        GetStageCode,&
                        GetStartMode,&
                        GetStressLeaf,&
                        GetStressSenescence ,&
                        GetStressSFadjNEW,&
                        GetStressTot,&
                        GetSumETo,&
                        GetSumGDD,&
                        GetSumGDDcuts,&
                        GetSumGDDPrev,&
                        GetSumInterval,&
                        GetSumKci,&
                        GetSumKcTop,&
                        GetSumKcTopStress,&
                        GetTadj,&
                        GetTheProjectFile,&
                        GetTimeSenescence ,&
                        GetTransfer,&
                        GetWaterTableInProfile,&
                        GetWeedRCi,&
                        GetYprevSum,&
                        GetZeval,&
                        GetZiprev,&
                        InitializeClimate,&
                        InitializeRunPart1, &
                        InitializeRunPart2, &
                        InitializeSimulation, &
                        InitializeSimulationRunPart2, &
                        ReadClimateNextDay, &
                        SetalfaHI,&
                        SetalfaHIAdj,&
                        SetBin,&
                        SetBout,&
                        SetBprevSum,&
                        SetCCiActualWeedInfested,&
                        SetCCoTotal,&
                        SetCCxCropWeedsNoSFstress,&
                        SetCCxTotal,&
                        SetCCxWitheredTpotNoS,&
                        SetCDCTotal,&
                        SetCGCref,&
                        SetCO2i,&
                        SetCoeffb0,&
                        SetCoeffb0Salt,&
                        SetCoeffb1,&
                        SetCoeffb1Salt,&
                        SetCoeffb2,&
                        SetCoeffb2Salt,&
                        SetCutInfoRecord1,&
                        SetCutInfoRecord2,&
                        SetDayFraction,&
                        SetDayLastCut,&
                        SetDayNr1Eval,&
                        SetDayNrEval,&
                        SetDayNri,&
                        SetFracBiomassPotSF,&
                        SetfWeedNoS,&
                        SetGDDayFraction,&
                        SetGDDayi,&
                        SetGDDCDCTotal,&
                        SetGDDCGCref ,&
                        SetGDDTadj,&
                        SetGDDVariablesNextDay, &
                        SetGwTable,&
                        SetHItimesAT,&
                        SetHItimesAT1,&
                        SetHItimesAT2,&
                        SetHItimesBEF,&
                        SetIrriInfoRecord1,&
                        SetIrriInfoRecord2,&
                        SetIrriInterval,&
                        SetLineNrEval,&
                        SetNextSimFromDayNr ,&
                        SetNoMoreCrop,&
                        SetNoYear,&
                        SetNrCut,&
                        SetPlotVarCrop,&
                        SetPreviousBmob,&
                        SetPreviousBsto,&
                        SetPreviousDayNr,&
                        SetPreviousStressLevel,&
                        SetPreviousSum,&
                        SetPreviousSumETo,&
                        SetPreviousSumGDD,&
                        SetScorAT1,&
                        SetScorAT2,&
                        SetStageCode,&
                        SetStartMode,&
                        SetStressLeaf,&
                        SetStressSenescence ,&
                        SetStressSFadjNEW,&
                        SetStressTot,&
                        SetSumETo,&
                        SetSumGDD,&
                        SetSumGDDcuts,&
                        SetSumGDDPrev,&
                        SetSumInterval,&
                        SetSumKci,&
                        SetSumKcTop,&
                        SetSumKcTopStress,&
                        SetTadj,&
                        SetTimeSenescence ,&
                        SetTransfer,&
                        SetWaterTableInProfile,&
                        SetWeedRCi,&
                        SetYprevSum,&
                        SetZeval,&
                        SetZiprev

! From ac_project_input
   use ac_project_input, only: ProjectInput, set_project_input
! From ac_kinds
    use ac_kinds, only: intEnum, &
                        int32, &
                        int8, &
                        sp

! From ac_startunit
    use ac_startunit, only:     FinalizeTheProgram, &
                                GetListProjectsFile, &
                                GetNumberOfProjects, &
                                GetProjectFileName, &
                                GetProjectType, &
                                GetSimulation_NrRuns, &
                                InitializeTheProgram, &
                                InitializeProject, &
                                WriteProjectsInfo

    implicit none
! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    integer              :: itemp, countertemp
    real                 :: dt
    real                 :: lat, lon
    real                 :: tmp_elev
    integer              :: row, col
    integer              :: year, month, day, hour, minute, second
    integer              :: simul_len
    logical              :: alarmCheck
    real                 :: TminRun_i, TmaxRun_i

    integer              :: status, ierr
    integer              :: c,r,l
    integer              :: ios, nid,rivid,fldid

    integer              :: tid

    integer              :: daynr, todaynr, iproject, nprojects
    logical              :: ListProjectFileExist
    character(len=:), allocatable :: ListProjectsFile, TheProjectFile
    integer              :: Crop_DaysToGermination, Crop_DaysToMaxRooting, Crop_DaysToFlowering
    integer              :: Crop_DaysToHarvest, Crop_DaysTosenescence, Crop_DaysToCCini
    integer              :: Crop_DaysToFullCanopy, Crop_DaysToFullCanopySF, Crop_DaysToHIo
    integer(int32) :: temp1    
    integer              :: irr_record_flag, DNr ! for irri file management
    character(250)       :: TempStr 

    real                 :: tmp_pres, tmp_precip, tmp_tmax, tmp_tmin   ! Weather Forcing
    real                 :: tmp_tdew, tmp_swrad, tmp_wind, tmp_eto     ! Weather Forcing

    ! For AdvanceOneTimeStep
    real                 :: tmp_wpi
!
! !DESCRIPTION:
!  This is the entry point for calling the AC72 physics.
!  This routine calls the {\AdvanceOneTimeStep} routine that runs AquaCrop

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for AC72
    ! check AC72 alarm. If alarm is ring, run model.
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "AC72 model alarm")
    if (alarmCheck) Then
        if (AC72_struc(n)%ac72(1)%read_Trecord.eq.1) then
        ! Read T record of next sim period
            call ac72_read_Trecord(n)
        endif
        do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            dt = LIS_rc%ts
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
            tmp_elev = LIS_domain(n)%tile(t)%elev

            !!------ This Block Is Where We Obtain Weather Forcing ------------------------------!!
            ! retrieve forcing data from AC72_struc(n)%ac72(t) and assign to local variables

            ! PRES: Daily average surface pressure (kPa)
            tmp_pres      = (AC72_struc(n)%ac72(t)%psurf / AC72_struc(n)%forc_count) / 1000

            ! PRECIP: Total daily precipitation (rain+snow) (mm)
            tmp_precip    = (AC72_struc(n)%ac72(t)%prcp / AC72_struc(n)%forc_count) * 3600. * 24. !Convert from kg/ms2 to mm/d

            ! TMAX: maximum daily air temperature (degC)
            tmp_tmax      = AC72_struc(n)%ac72(t)%tmax - LIS_CONST_TKFRZ !Convert from K to C

            ! TMIN: minimum daily air temperature (degC)
            tmp_tmin      = AC72_struc(n)%ac72(t)%tmin - LIS_CONST_TKFRZ !Convert from K to C 

            ! TDEW: average daily dewpoint temperature (degC)
            tmp_tdew      = (AC72_struc(n)%ac72(t)%tdew / AC72_struc(n)%forc_count) - LIS_CONST_TKFRZ !Convert from K to C

            ! SW_RAD: daily total incoming solar radiation (MJ/(m2d))
            tmp_swrad     = (AC72_struc(n)%ac72(t)%swdown / AC72_struc(n)%forc_count) * 0.0864 !Convert from W/m2 to MJ/(m2d)

            ! Wind: daily average wind speed (m/s)
            tmp_wind      = AC72_struc(n)%ac72(t)%wndspd / AC72_struc(n)%forc_count

            ! check validity of PRES
            if(tmp_pres .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable PRES in AC72"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! check validity of PRECIP
            if(tmp_precip .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable PRECIP in AC72"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! check validity of TMAX
            if(tmp_tmax .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable TMAX in AC72"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

           ! check validity of TMIN
            if(tmp_tmin .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable TMIN in AC72"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! check validity of TDEW
            if(tmp_tdew .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable TDEW in AC72"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! check validity of SW_RAD
            if(tmp_swrad .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable SW_RAD in AC72"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! check validity of WIND
            if(tmp_wind .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable WIND in AC72"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! Call ETo subroutine           
            call ac72_ETo_calc(tmp_pres, tmp_tmax, tmp_tmin, tmp_tdew, &
                               tmp_wind, tmp_swrad, &
                               tmp_elev, lat, tmp_eto)
            AC72_struc(n)%ac72(t)%eto = tmp_eto

            ! setting all global variables
            call SetRootZoneWC_Actual(AC72_struc(n)%ac72(t)%RootZoneWC_Actual)
            call SetRootZoneWC_FC(AC72_struc(n)%ac72(t)%RootZoneWC_FC)
            call SetRootZoneWC_WP(AC72_struc(n)%ac72(t)%RootZoneWC_WP)
            call SetRootZoneWC_SAT(AC72_struc(n)%ac72(t)%RootZoneWC_SAT)
            call SetRootZoneWC_Leaf(AC72_struc(n)%ac72(t)%RootZoneWC_Leaf)
            call SetRootZoneWC_Thresh(AC72_struc(n)%ac72(t)%RootZoneWC_Thresh)
            call SetRootZoneWC_Sen(AC72_struc(n)%ac72(t)%RootZoneWC_Sen)
            call SetRootZoneWC_ZtopAct(AC72_struc(n)%ac72(t)%RootZoneWC_ZtopAct)
            call SetRootZoneWC_ZtopFC(AC72_struc(n)%ac72(t)%RootZoneWC_ZtopFC)
            call SetRootZoneWC_ZtopWP(AC72_struc(n)%ac72(t)%RootZoneWC_ZtopWP)
            call SetRootZoneWC_ZtopThresh(AC72_struc(n)%ac72(t)%RootZoneWC_ZtopThresh)
            call SetCompartment(AC72_struc(n)%ac72(t)%Compartment)
            call SetTotalSaltContent(AC72_struc(n)%ac72(t)%TotalSaltContent)
            call SetTotalWaterContent(AC72_struc(n)%ac72(t)%TotalWaterContent)
            call Seteffectiverain(AC72_struc(n)%ac72(t)%effectiverain)
            call SetSumWaBal(AC72_struc(n)%ac72(t)%SumWaBal)
            call SetRootZoneSalt(AC72_struc(n)%ac72(t)%RootZoneSalt)
            call SetSimulation(AC72_struc(n)%ac72(t)%Simulation)
            call SetIrriInterval(AC72_struc(n)%ac72(t)%IrriInterval)
            call SetIrriInfoRecord1(AC72_struc(n)%ac72(t)%IrriInfoRecord1)
            call SetIrriInfoRecord2(AC72_struc(n)%ac72(t)%IrriInfoRecord2)
            call SetIrrigation(AC72_struc(n)%ac72(t)%Irrigation)
            do l=1, AC72_struc(n)%ac72(t)%NrCompartments
                 call SetCompartment_theta(l,AC72_struc(n)%ac72(t)%smc(l))
            enddo
            call SetIrriECw(AC72_struc(n)%ac72(t)%IrriECw) 
            call SetManagement(AC72_struc(n)%ac72(t)%Management) 
            call SetPerennialPeriod(AC72_struc(n)%ac72(t)%PerennialPeriod) 
            call SetSimulParam(AC72_struc(n)%ac72(t)%simulparam) 
            call SetManagement_Cuttings(AC72_struc(n)%ac72(t)%Cuttings) 
            call SetOnset(AC72_struc(n)%ac72(t)%onset) 
            call SetEndSeason(AC72_struc(n)%ac72(t)%endseason) 
            call SetCrop(AC72_struc(n)%ac72(t)%crop) 
            call SetSoil(AC72_struc(n)%ac72(t)%Soil) 
            call SetIrriBeforeSeason(AC72_struc(n)%ac72(t)%IrriBeforeSeason)
            call SetIrriAfterSeason(AC72_struc(n)%ac72(t)%IrriAfterSeason)
            call SetSoilLayer(AC72_struc(n)%ac72(t)%soillayer)
            call SetDayNri(AC72_struc(n)%ac72(t)%daynri)

            call SetGenerateTimeMode(AC72_struc(n)%ac72(t)%GenerateTimeMode) 
            call SetGenerateDepthMode(AC72_struc(n)%ac72(t)%GenerateDepthMode) 
            call SetIrriMode(AC72_struc(n)%ac72(t)%IrriMode) 
            call SetDaySubmerged(AC72_struc(n)%ac72(t)%DaySubmerged) 
            call SetMaxPlotNew(AC72_struc(n)%ac72(t)%MaxPlotNew) 
            call SetNrCompartments(AC72_struc(n)%ac72(t)%NrCompartments) 
            call SetIrriFirstDayNr(AC72_struc(n)%ac72(t)%IrriFirstDayNr) 
            call SetZiAqua(AC72_struc(n)%ac72(t)%ZiAqua) 
            call SetIniPercTAW(AC72_struc(n)%ac72(t)%IniPercTAW) 
            call SetMaxPlotTr(AC72_struc(n)%ac72(t)%MaxPlotTr)
            call SetOutputAggregate(AC72_struc(n)%ac72(t)%OutputAggregate) 

            call SetEvapoEntireSoilSurface(AC72_struc(n)%ac72(t)%EvapoEntireSoilSurface) 
            call SetPreDay(AC72_struc(n)%ac72(t)%PreDay) 
            call SetOutDaily(AC72_struc(n)%ac72(t)%OutDaily) 
            call SetOut1Wabal(AC72_struc(n)%ac72(t)%Out1Wabal) 
            call SetOut2Crop(AC72_struc(n)%ac72(t)%Out2Crop) 
            call SetOut3Prof(AC72_struc(n)%ac72(t)%Out3Prof) 
            call SetOut4Salt(AC72_struc(n)%ac72(t)%Out4Salt) 
            call SetOut5CompWC(AC72_struc(n)%ac72(t)%Out5CompWC) 
            call SetOut6CompEC(AC72_struc(n)%ac72(t)%Out6CompEC) 
            call SetOut7Clim(AC72_struc(n)%ac72(t)%Out7Clim)
            call SetOut8Irri(AC72_struc(n)%ac72(t)%Out8Irri) 
            call SetPart1Mult(AC72_struc(n)%ac72(t)%Part1Mult) 
            call SetPart2Eval(AC72_struc(n)%ac72(t)%Part2Eval) 

            !
            call SetCCiActual(AC72_struc(n)%ac72(t)%CCiActual)
            call SetCCiprev(AC72_struc(n)%ac72(t)%CCiprev)

            call SetCCiTopEarlySen(AC72_struc(n)%ac72(t)%CCiTopEarlySen)
            call SetCRsalt(AC72_struc(n)%ac72(t)%CRsalt)
            call SetCRwater(AC72_struc(n)%ac72(t)%CRwater)
            call SetECdrain(AC72_struc(n)%ac72(t)%ECdrain)
            call SetEciAqua(AC72_struc(n)%ac72(t)%ECiAqua)
            call SetECstorage(AC72_struc(n)%ac72(t)%ECstorage)
            call SetEact(AC72_struc(n)%ac72(t)%Eact)
            call SetEpot(AC72_struc(n)%ac72(t)%Epot)
            call SetDrain(AC72_struc(n)%ac72(t)%Drain)
            call SetInfiltrated(AC72_struc(n)%ac72(t)%Infiltrated)
            call SetRootingDepth(AC72_struc(n)%ac72(t)%RootingDepth)
            call SetRunoff(AC72_struc(n)%ac72(t)%Runoff)
            call SetSaltInfiltr(AC72_struc(n)%ac72(t)%SaltInfiltr)
            call SetSurf0(AC72_struc(n)%ac72(t)%Surf0)
            call SetSurfaceStorage(AC72_struc(n)%ac72(t)%SurfaceStorage)
            call SetTact(AC72_struc(n)%ac72(t)%Tact)
            call SetTpot(AC72_struc(n)%ac72(t)%Tpot)
            call SetTactWeedInfested(AC72_struc(n)%ac72(t)%TactWeedInfested)

            call SetGwTable(AC72_struc(n)%ac72(t)%GwTable)
            call SetPlotVarCrop(AC72_struc(n)%ac72(t)%PlotVarCrop)
            call SetStressTot(AC72_struc(n)%ac72(t)%StressTot)
            call SetCutInfoRecord1(AC72_struc(n)%ac72(t)%CutInfoRecord1)
            call SetCutInfoRecord2(AC72_struc(n)%ac72(t)%CutInfoRecord2)
            call SetTransfer(AC72_struc(n)%ac72(t)%Transfer)
            call SetPreviousSum(AC72_struc(n)%ac72(t)%PreviousSum)
            call SetTadj(AC72_struc(n)%ac72(t)%Tadj)
            call SetGDDTadj(AC72_struc(n)%ac72(t)%GDDTadj)
            call SetDayLastCut(AC72_struc(n)%ac72(t)%DayLastCut)
            call SetNrCut(AC72_struc(n)%ac72(t)%NrCut)
            call SetSumInterval(AC72_struc(n)%ac72(t)%SumInterval)
            call SetPreviousStressLevel(AC72_struc(n)%ac72(t)%PreviousStressLevel)
            call SetStressSFadjNEW(AC72_struc(n)%ac72(t)%StressSFadjNEW)
            call SetBin(AC72_struc(n)%ac72(t)%Bin)
            call SetBout(AC72_struc(n)%ac72(t)%Bout)
            call SetCO2i(AC72_struc(n)%ac72(t)%CO2i)
            call SetFracBiomassPotSF(AC72_struc(n)%ac72(t)%FracBiomassPotSF)
            call SetSumETo(AC72_struc(n)%ac72(t)%SumETo)
            call SetSumGDD(AC72_struc(n)%ac72(t)%SumGDD)
            call SetZiprev(AC72_struc(n)%ac72(t)%Ziprev)
            call SetSumGDDPrev(AC72_struc(n)%ac72(t)%SumGDDPrev)
            call SetCCxWitheredTpotNoS(AC72_struc(n)%ac72(t)%CCxWitheredTpotNoS)
            call SetCoeffb0(AC72_struc(n)%ac72(t)%Coeffb0)
            call SetCoeffb1(AC72_struc(n)%ac72(t)%Coeffb1)
            call SetCoeffb2(AC72_struc(n)%ac72(t)%Coeffb2)
            call SetCoeffb0Salt(AC72_struc(n)%ac72(t)%Coeffb0Salt)
            call SetCoeffb1Salt(AC72_struc(n)%ac72(t)%Coeffb1Salt)
            call SetCoeffb2Salt(AC72_struc(n)%ac72(t)%Coeffb2Salt)
            call SetStressLeaf(AC72_struc(n)%ac72(t)%StressLeaf)
            call SetStressSenescence(AC72_struc(n)%ac72(t)%StressSenescence)
            call SetDayFraction(AC72_struc(n)%ac72(t)%DayFraction)
            call SetGDDayFraction(AC72_struc(n)%ac72(t)%GDDayFraction)
            call SetCGCref(AC72_struc(n)%ac72(t)%CGCref)
            call SetGDDCGCref(AC72_struc(n)%ac72(t)%GDDCGCref)
            call SetTimeSenescence(AC72_struc(n)%ac72(t)%TimeSenescence)
            call SetSumKcTop(AC72_struc(n)%ac72(t)%SumKcTop)
            call SetSumKcTopStress(AC72_struc(n)%ac72(t)%SumKcTopStress)
            call SetSumKci(AC72_struc(n)%ac72(t)%SumKci)
            call SetCCoTotal(AC72_struc(n)%ac72(t)%CCoTotal)
            call SetCCxTotal(AC72_struc(n)%ac72(t)%CCxTotal)
            call SetCDCTotal(AC72_struc(n)%ac72(t)%CDCTotal)
            call SetGDDCDCTotal(AC72_struc(n)%ac72(t)%GDDCDCTotal)
            call SetCCxCropWeedsNoSFstress(AC72_struc(n)%ac72(t)%CCxCropWeedsNoSFstress)
            call SetWeedRCi(AC72_struc(n)%ac72(t)%WeedRCi)
            call SetCCiActualWeedInfested(AC72_struc(n)%ac72(t)%CCiActualWeedInfested)
            call SetfWeedNoS(AC72_struc(n)%ac72(t)%fWeedNoS)
            call SetZeval(AC72_struc(n)%ac72(t)%Zeval)
            call SetBprevSum(AC72_struc(n)%ac72(t)%BprevSum)
            call SetYprevSum(AC72_struc(n)%ac72(t)%YprevSum)
            call SetSumGDDcuts(AC72_struc(n)%ac72(t)%SumGDDcuts)
            call SetHItimesBEF(AC72_struc(n)%ac72(t)%HItimesBEF)
            call SetScorAT1(AC72_struc(n)%ac72(t)%ScorAT1)
            call SetScorAT2(AC72_struc(n)%ac72(t)%ScorAT2)
            call SetHItimesAT1(AC72_struc(n)%ac72(t)%HItimesAT1)
            call SetHItimesAT2(AC72_struc(n)%ac72(t)%HItimesAT2)
            call SetHItimesAT(AC72_struc(n)%ac72(t)%HItimesAT)
            call SetalfaHI(AC72_struc(n)%ac72(t)%alfaHI)
            call SetalfaHIAdj(AC72_struc(n)%ac72(t)%alfaHIAdj)
            call SetNextSimFromDayNr(AC72_struc(n)%ac72(t)%NextSimFromDayNr)
            call SetDayNr1Eval(AC72_struc(n)%ac72(t)%DayNr1Eval)
            call SetDayNrEval(AC72_struc(n)%ac72(t)%DayNrEval)
            call SetLineNrEval(int(AC72_struc(n)%ac72(t)%LineNrEval,kind=int32))
            call SetPreviousSumETo(AC72_struc(n)%ac72(t)%PreviousSumETo)
            call SetPreviousSumGDD(AC72_struc(n)%ac72(t)%PreviousSumGDD)
            call SetPreviousBmob(AC72_struc(n)%ac72(t)%PreviousBmob)
            call SetPreviousBsto(AC72_struc(n)%ac72(t)%PreviousBsto)
            call SetStageCode(AC72_struc(n)%ac72(t)%StageCode)
            call SetPreviousDayNr(AC72_struc(n)%ac72(t)%PreviousDayNr)
            call SetNoYear(AC72_struc(n)%ac72(t)%NoYear)
            call SetWaterTableInProfile(AC72_struc(n)%ac72(t)%WaterTableInProfile)
            call SetStartMode(AC72_struc(n)%ac72(t)%StartMode)
            call SetGDDayi(AC72_struc(n)%ac72(t)%GDDayi)
            call SetNoMoreCrop(AC72_struc(n)%AC72(t)%NoMoreCrop)

            ! Fixed var
            call SetOut3Prof(.true.) ! needed for correct rootzone sm
            call SetOutDaily(.true.)

            call SetTminRun(AC72_struc(n)%ac72(t)%Tmin_record)
            call SetTmaxRun(AC72_struc(n)%ac72(t)%Tmax_record)

            ! Set climate variables
            ! Round them to 4 digits after the comma as done in the AC standalone
            tmp_precip = anint(tmp_precip*10000)/10000
            call SetRain(tmp_precip)
            tmp_tmin = anint(tmp_tmin*10000)/10000
            call SetTmin(tmp_tmin)
            tmp_tmax = anint(tmp_tmax*10000)/10000
            call SetTmax(tmp_tmax)
            tmp_eto = anint(tmp_eto*10000)/10000
            call SetETo(tmp_eto)
            
            ! SumGDD calculation needed only for second day when not done within InitializeSimulationRunPart2
            if (GetDayNri()>GetSimulation_FromDayNr()) then
                ! Sum of GDD at end of first day
                call SetGDDayi(DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), GetTmin(), &
                        GetTmax(), GetSimulParam_GDDMethod()))
                if (GetDayNri() >= GetCrop_Day1()) then
                    call SetSimulation_SumGDD(GetSimulation_SumGDD() + GetGDDayi())
                    call SetSimulation_SumGDDfromDay1(GetSimulation_SumGDDfromDay1() + &
                        GetGDDayi())
                end if
            end if

            ! Start irrigation block
            call SetIrriMode(AC72_struc(n)%ac72(t)%IrriMode)
            irr_record_flag = 0
            if(AC72_struc(n)%ac72(t)%IrriMode.ne.IrriMode_NoIrri) then
                ! Irrigation ON for tile
                call SetIrriAfterSeason(AC72_struc(n)%ac72(t)%IrriAfterSeason)
                call SetIrriBeforeSeason(AC72_struc(n)%ac72(t)%IrriBeforeSeason)
                call SetIrriECw(AC72_struc(n)%ac72(t)%IrriECw) 
                call SetIrriInfoRecord1(AC72_struc(n)%ac72(t)%IrriInfoRecord1)
                call SetIrriInfoRecord2(AC72_struc(n)%ac72(t)%IrriInfoRecord2)
                call SetIrriMethod(AC72_struc(n)%ac72(t)%IrriMethod)
                call SetGenerateDepthMode(AC72_struc(n)%ac72(t)%GenerateDepthMode)
                call SetGenerateTimeMode(AC72_struc(n)%ac72(t)%GenerateTimeMode)

                ! For IrriMode_Generate and IrriMode_Manual
                ! Check if new record needs to be read
                if((AC72_struc(n)%ac72(t)%IrriMode.eq.IrriMode_Generate) &
                     .and. (GetDayNri() >= GetCrop_Day1()) &
                     .and. (GetDayNri() <= GetCrop_DayN())) then
                    if((AC72_struc(n)%ac72(t)%daynri-AC72_struc(n)%ac72(t)%Crop%Day1+1)&
                    .gt.AC72_struc(n)%ac72(t)%IrriInfoRecord1%ToDay) then
                        irr_record_flag = 1
                    endif
                elseif(AC72_struc(n)%ac72(t)%IrriMode.eq.IrriMode_Manual) then
                    ! Check start date of schedule
                    if(AC72_struc(n)%ac72(t)%IrriFirstDayNr.eq.undef_int)then
                        DNr = AC72_struc(n)%ac72(t)%daynri &
                              - AC72_struc(n)%ac72(t)%Crop%Day1 + 1
                    else
                        DNr = AC72_struc(n)%ac72(t)%daynri &
                              - AC72_struc(n)%ac72(t)%IrriFirstDayNr + 1
                    endif
                    if(AC72_struc(n)%ac72(t)%IrriInfoRecord1%TimeInfo.eq.DNr)then
                        irr_record_flag = 1
                    endif
                endif

                ! re-open irrigation file and read the previous records
                if(irr_record_flag.eq.1)then
                    call fIrri_open(trim(AC72_struc(n)%PathNameSimul)&
                                    //trim(AC72_struc(n)%Irrigation_Filename), 'r')
                        do i=1,AC72_struc(n)%ac72(t)%irri_lnr
                            TempStr = fIrri_read()
                        enddo
                    AC72_struc(n)%ac72(t)%irri_lnr = AC72_struc(n)%ac72(t)%irri_lnr + 1
                endif
            endif
            ! End irrigation block

            !!! initialize run (year)
            if (AC72_struc(n)%ac72(t)%InitializeRun.eq.1) then !make it flex
                call SetClimRecord_DataType(0_int8)
                call SetClimRecord_fromd(0)
                call SetClimRecord_fromdaynr(ProjectInput(1)%Simulation_DayNr1)
                call SetClimRecord_fromm(0)
                call SetClimRecord_fromstring('any date')
                call SetClimRecord_fromy(LIS_rc%syr)
                call SetClimRecord_NrObs(999)
                call SetClimRecord_tod(0)
                call SetClimRecord_todaynr(ProjectInput(GetSimulation_NrRuns())%Simulation_DayNrN)
                call SetClimRecord_tom(0)
                call SetClimRecord_tostring('any date')
                call SetClimRecord_toy(0)
                call SetClimFile('(External)')

                AC72_struc(n)%ac72(t)%WPi = 0.

                ! Set crop file (crop parameters are read when calling InitializeRunPart1)
                call set_project_input(AC72_struc(n)%ac72(t)%irun, &
                                       'Crop_Filename', &
                                        trim(AC72_struc(n)%ac72(t)%cropt)//'.CRO')

                ! Set Global variable to pass T record to AquaCrop
                call SetTminRun(AC72_struc(n)%ac72(t)%Tmin_record)
                call SetTmaxRun(AC72_struc(n)%ac72(t)%Tmax_record)

                ! Set Tmin and Tmax reference to compute the stress realtions
                call SetTminTnxReference12MonthsRun(AC72_struc(n)%ac72(t)%tmincli_monthly(:))
                call SetTmaxTnxReference12MonthsRun(AC72_struc(n)%ac72(t)%tmaxcli_monthly(:))

                ! Set reference year for CO2 for stress functions
                call SetTnxReferenceYear(AC72_struc(n)%tempcli_refyr)

                call SetTnxReferenceFile('(External)')

                ! InitializeRunPart
                call InitializeRunPart1(int(AC72_struc(n)%ac72(t)%irun,kind=int8), AC72_struc(n)%ac72(t)%TheProjectType)
                call InitializeSimulationRunPart2()
                AC72_struc(n)%ac72(t)%HarvestNow = .false. ! Initialize to false
                ! Check if enough GDDays to complete cycle
                if(GetCrop_ModeCycle().eq.ModeCycle_GDDays)then
                    if (((GetCrop_Day1()+GetCrop_DaysToHarvest()).gt.GetSimulation_ToDayNr()) &
                        .or.(GetCrop_DaysToHarvest()<1)) then
                        AC72_struc(n)%ac72(t)%cycle_complete = 0
                    else
                        AC72_struc(n)%ac72(t)%cycle_complete = 1
                    endif
                endif

                ! Overwrite the SMC to avoid problems when restarting
                do l=1, AC72_struc(n)%ac72(t)%NrCompartments
                    call SetCompartment_theta(l,AC72_struc(n)%ac72(t)%smc(l))
                enddo

                ! Irrigaton file management after InitializeRun
                if(GetIrriMode().ne.IrriMode_NoIrri) then
                    call fIrri_close()
                    ! Check for irrigation (irrigation file management)
                    if(GetIrriMode().eq.IrriMode_Manual)then
                        if(GetIrriInfoRecord1_NoMoreInfo())then
                            AC72_struc(n)%ac72(t)%irri_lnr = 9
                        else
                            AC72_struc(n)%ac72(t)%irri_lnr = 10
                        endif
                    elseif(GetIrriMode().eq.IrriMode_Generate)then
                        if(AC72_struc(n)%ac72(t)%IrriInfoRecord1%NoMoreInfo)then
                            AC72_struc(n)%ac72(t)%irri_lnr = 11
                        else
                            AC72_struc(n)%ac72(t)%irri_lnr = 12
                        endif
                    else ! no irrigation, set to 0
                        AC72_struc(n)%ac72(t)%irri_lnr = 0
                    endif
                endif
                ! End irrigation block
                AC72_struc(n)%ac72(t)%InitializeRun = 0 ! Initialization done
                AC72_struc(n)%ac72(t)%read_Trecord = 0
            end if

            ! Run AC
            tmp_wpi = AC72_struc(n)%ac72(t)%WPi
            call AdvanceOneTimeStep(tmp_wpi, AC72_struc(n)%ac72(t)%HarvestNow)
            AC72_struc(n)%ac72(t)%WPi = tmp_wpi

            ! Close irri file if opened
            if(irr_record_flag.eq.1)then
                call fIrri_close()
            endif

            ! Get all the ac72 variables and store in AC72_struc
            do l=1, AC72_struc(n)%ac72(t)%NrCompartments
                    AC72_struc(n)%ac72(t)%smc(l) = GetCompartment_theta(l)
            enddo
            
            AC72_struc(n)%ac72(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
            AC72_struc(n)%ac72(t)%RootZoneWC_FC = GetRootZoneWC_FC()
            AC72_struc(n)%ac72(t)%RootZoneWC_WP = GetRootZoneWC_WP()
            AC72_struc(n)%ac72(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
            AC72_struc(n)%ac72(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
            AC72_struc(n)%ac72(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
            AC72_struc(n)%ac72(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
            AC72_struc(n)%ac72(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
            AC72_struc(n)%ac72(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
            AC72_struc(n)%ac72(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
            AC72_struc(n)%ac72(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
            AC72_struc(n)%ac72(t)%Compartment = GetCompartment()
            AC72_struc(n)%ac72(t)%TotalSaltContent = GetTotalSaltContent()
            AC72_struc(n)%ac72(t)%TotalWaterContent = GetTotalWaterContent()
            AC72_struc(n)%ac72(t)%effectiverain = Geteffectiverain()
            AC72_struc(n)%ac72(t)%SumWaBal = GetSumWaBal()
            AC72_struc(n)%ac72(t)%RootZoneSalt = GetRootZoneSalt()
            AC72_struc(n)%ac72(t)%Simulation = GetSimulation()
            AC72_struc(n)%ac72(t)%IrriInterval = GetIrriInterval()
            AC72_struc(n)%ac72(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
            AC72_struc(n)%ac72(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
            AC72_struc(n)%ac72(t)%Irrigation = GetIrrigation()
            AC72_struc(n)%ac72(t)%IrriBeforeSeason = GetIrriBeforeSeason()
            AC72_struc(n)%ac72(t)%IrriAfterSeason = GetIrriAfterSeason()
            AC72_struc(n)%ac72(t)%SoilLayer = GetSoilLayer()
            AC72_struc(n)%ac72(t)%daynri = GetDayNri()
            do l=1, AC72_struc(n)%ac72(t)%NrCompartments
                AC72_struc(n)%ac72(t)%smc(l) = GetCompartment_theta(l)
            enddo
            AC72_struc(n)%ac72(t)%IrriECw = GetIrriECw()
            AC72_struc(n)%ac72(t)%Management = GetManagement()
            AC72_struc(n)%ac72(t)%PerennialPeriod = GetPerennialPeriod()
            AC72_struc(n)%ac72(t)%simulparam = GetSimulParam()
            AC72_struc(n)%ac72(t)%Cuttings = GetManagement_Cuttings()
            AC72_struc(n)%ac72(t)%onset = GetOnset()
            AC72_struc(n)%ac72(t)%endseason = GetEndSeason()
            AC72_struc(n)%ac72(t)%crop = GetCrop()
            AC72_struc(n)%ac72(t)%Soil = GetSoil()
            AC72_struc(n)%ac72(t)%TemperatureRecord = GetTemperatureRecord()
            AC72_struc(n)%ac72(t)%ClimRecord = GetClimRecord()
            AC72_struc(n)%ac72(t)%RainRecord = GetRainRecord()
            AC72_struc(n)%ac72(t)%EToRecord = GetEToRecord()

            AC72_struc(n)%ac72(t)%GenerateTimeMode = GetGenerateTimeMode()
            AC72_struc(n)%ac72(t)%GenerateDepthMode = GetGenerateDepthMode()
            AC72_struc(n)%ac72(t)%IrriMode = GetIrriMode()
            AC72_struc(n)%ac72(t)%IrriMethod = GetIrriMethod()
            AC72_struc(n)%ac72(t)%DaySubmerged = GetDaySubmerged()
            AC72_struc(n)%ac72(t)%MaxPlotNew = GetMaxPlotNew()
            AC72_struc(n)%ac72(t)%NrCompartments = GetNrCompartments()
            AC72_struc(n)%ac72(t)%IrriFirstDayNr = GetIrriFirstDayNr()
            AC72_struc(n)%ac72(t)%ZiAqua = GetZiAqua()
            AC72_struc(n)%ac72(t)%IniPercTAW = GetIniPercTAW()
            AC72_struc(n)%ac72(t)%MaxPlotTr = GetMaxPlotTr()
            AC72_struc(n)%ac72(t)%OutputAggregate = GetOutputAggregate()

            AC72_struc(n)%ac72(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
            AC72_struc(n)%ac72(t)%PreDay = GetPreDay()
            AC72_struc(n)%ac72(t)%OutDaily = GetOutDaily()
            AC72_struc(n)%ac72(t)%Out1Wabal = GetOut1Wabal()
            AC72_struc(n)%ac72(t)%Out2Crop = GetOut2Crop()
            AC72_struc(n)%ac72(t)%Out3Prof = GetOut3Prof()
            AC72_struc(n)%ac72(t)%Out4Salt = GetOut4Salt()
            AC72_struc(n)%ac72(t)%Out5CompWC = GetOut5CompWC()
            AC72_struc(n)%ac72(t)%Out6CompEC = GetOut6CompEC()
            AC72_struc(n)%ac72(t)%Out7Clim = GetOut7Clim()
            AC72_struc(n)%ac72(t)%Out8Irri = GetOut8Irri()
            AC72_struc(n)%ac72(t)%Part1Mult = GetPart1Mult()
            AC72_struc(n)%ac72(t)%Part2Eval = GetPart2Eval()

            !
            AC72_struc(n)%ac72(t)%CCiActual = GetCCiActual()
            AC72_struc(n)%ac72(t)%CCiprev = GetCCiprev()
            AC72_struc(n)%ac72(t)%CCiTopEarlySen = GetCCiTopEarlySen()
            AC72_struc(n)%ac72(t)%CRsalt = GetCRsalt ()
            AC72_struc(n)%ac72(t)%CRwater = GetCRwater()
            AC72_struc(n)%ac72(t)%ECdrain = GetECdrain()
            AC72_struc(n)%ac72(t)%ECiAqua = GetECiAqua()
            AC72_struc(n)%ac72(t)%ECstorage = GetECstorage()
            AC72_struc(n)%ac72(t)%Eact = GetEact()
            AC72_struc(n)%ac72(t)%Epot = GetEpot()
            AC72_struc(n)%ac72(t)%ETo = GetETo()
            AC72_struc(n)%ac72(t)%Drain = GetDrain()
            AC72_struc(n)%ac72(t)%Infiltrated = GetInfiltrated()
            AC72_struc(n)%ac72(t)%prcp = GetRain()
            AC72_struc(n)%ac72(t)%RootingDepth = GetRootingDepth()
            AC72_struc(n)%ac72(t)%Runoff = GetRunoff()
            AC72_struc(n)%ac72(t)%SaltInfiltr = GetSaltInfiltr()
            AC72_struc(n)%ac72(t)%Surf0 = GetSurf0()
            AC72_struc(n)%ac72(t)%SurfaceStorage = GetSurfaceStorage()
            AC72_struc(n)%ac72(t)%Tact = GetTact()
            AC72_struc(n)%ac72(t)%Tpot = GetTpot()
            AC72_struc(n)%ac72(t)%TactWeedInfested = GetTactWeedInfested()
            AC72_struc(n)%ac72(t)%tmax = GetTmax()
            AC72_struc(n)%ac72(t)%tmin =GetTmin()


            AC72_struc(n)%ac72(t)%GwTable = GetGwTable()
            AC72_struc(n)%ac72(t)%PlotVarCrop = GetPlotVarCrop()
            AC72_struc(n)%ac72(t)%StressTot = GetStressTot()
            AC72_struc(n)%ac72(t)%CutInfoRecord1 = GetCutInfoRecord1()
            AC72_struc(n)%ac72(t)%CutInfoRecord2 = GetCutInfoRecord2()
            AC72_struc(n)%ac72(t)%Transfer = GetTransfer()
            AC72_struc(n)%ac72(t)%PreviousSum = GetPreviousSum()
            AC72_struc(n)%ac72(t)%Tadj = GetTadj()
            AC72_struc(n)%ac72(t)%GDDTadj = GetGDDTadj()
            AC72_struc(n)%ac72(t)%DayLastCut = GetDayLastCut()
            AC72_struc(n)%ac72(t)%NrCut = GetNrCut()
            AC72_struc(n)%ac72(t)%SumInterval = GetSumInterval()
            AC72_struc(n)%ac72(t)%PreviousStressLevel = GetPreviousStressLevel()
            AC72_struc(n)%ac72(t)%StressSFadjNEW = GetStressSFadjNEW()
            AC72_struc(n)%ac72(t)%Bin = GetBin()
            AC72_struc(n)%ac72(t)%Bout = GetBout()
            AC72_struc(n)%ac72(t)%GDDayi = GetGDDayi()
            AC72_struc(n)%ac72(t)%CO2i = GetCO2i()
            AC72_struc(n)%ac72(t)%FracBiomassPotSF = GetFracBiomassPotSF()
            AC72_struc(n)%ac72(t)%SumETo = GetSumETo()
            AC72_struc(n)%ac72(t)%SumGDD = GetSumGDD()
            AC72_struc(n)%ac72(t)%Ziprev = GetZiprev()
            AC72_struc(n)%ac72(t)%SumGDDPrev = GetSumGDDPrev()
            AC72_struc(n)%ac72(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
            AC72_struc(n)%ac72(t)%Coeffb0 = GetCoeffb0()
            AC72_struc(n)%ac72(t)%Coeffb1 = GetCoeffb1()
            AC72_struc(n)%ac72(t)%Coeffb2 = GetCoeffb2()
            AC72_struc(n)%ac72(t)%Coeffb0Salt = GetCoeffb0Salt()
            AC72_struc(n)%ac72(t)%Coeffb1Salt = GetCoeffb1Salt()
            AC72_struc(n)%ac72(t)%Coeffb2Salt = GetCoeffb2Salt()
            AC72_struc(n)%ac72(t)%StressLeaf = GetStressLeaf()
            AC72_struc(n)%ac72(t)%StressSenescence = GetStressSenescence()
            AC72_struc(n)%ac72(t)%DayFraction = GetDayFraction()
            AC72_struc(n)%ac72(t)%GDDayFraction = GetGDDayFraction()
            AC72_struc(n)%ac72(t)%CGCref = GetCGCref()
            AC72_struc(n)%ac72(t)%GDDCGCref = GetGDDCGCref()
            AC72_struc(n)%ac72(t)%TimeSenescence = GetTimeSenescence()
            AC72_struc(n)%ac72(t)%SumKcTop = GetSumKcTop()
            AC72_struc(n)%ac72(t)%SumKcTopStress = GetSumKcTopStress()
            AC72_struc(n)%ac72(t)%SumKci = GetSumKci()
            AC72_struc(n)%ac72(t)%CCoTotal = GetCCoTotal()
            AC72_struc(n)%ac72(t)%CCxTotal = GetCCxTotal()
            AC72_struc(n)%ac72(t)%CDCTotal = GetCDCTotal()
            AC72_struc(n)%ac72(t)%GDDCDCTotal = GetGDDCDCTotal()
            AC72_struc(n)%ac72(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
            AC72_struc(n)%ac72(t)%WeedRCi = GetWeedRCi()
            AC72_struc(n)%ac72(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
            AC72_struc(n)%ac72(t)%fWeedNoS = GetfWeedNoS()
            AC72_struc(n)%ac72(t)%Zeval = GetZeval()
            AC72_struc(n)%ac72(t)%BprevSum = GetBprevSum()
            AC72_struc(n)%ac72(t)%YprevSum = GetYprevSum()
            AC72_struc(n)%ac72(t)%SumGDDcuts = GetSumGDDcuts()
            AC72_struc(n)%ac72(t)%HItimesBEF = GetHItimesBEF()
            AC72_struc(n)%ac72(t)%ScorAT1 = GetScorAT1()
            AC72_struc(n)%ac72(t)%ScorAT2 = GetScorAT2()
            AC72_struc(n)%ac72(t)%HItimesAT1 = GetHItimesAT1()
            AC72_struc(n)%ac72(t)%HItimesAT2 = GetHItimesAT2()
            AC72_struc(n)%ac72(t)%HItimesAT = GetHItimesAT()
            AC72_struc(n)%ac72(t)%alfaHI = GetalfaHI()
            AC72_struc(n)%ac72(t)%alfaHIAdj = GetalfaHIAdj()
            AC72_struc(n)%ac72(t)%NextSimFromDayNr = GetNextSimFromDayNr ()
            AC72_struc(n)%ac72(t)%DayNr1Eval = GetDayNr1Eval()
            AC72_struc(n)%ac72(t)%DayNrEval = GetDayNrEval()
            AC72_struc(n)%ac72(t)%LineNrEval = GetLineNrEval()
            AC72_struc(n)%ac72(t)%PreviousSumETo = GetPreviousSumETo()
            AC72_struc(n)%ac72(t)%PreviousSumGDD = GetPreviousSumGDD()
            AC72_struc(n)%ac72(t)%PreviousBmob = GetPreviousBmob()
            AC72_struc(n)%ac72(t)%PreviousBsto = GetPreviousBsto()
            AC72_struc(n)%ac72(t)%StageCode = GetStageCode()
            AC72_struc(n)%ac72(t)%PreviousDayNr = GetPreviousDayNr()
            AC72_struc(n)%ac72(t)%NoYear = GetNoYear()
            AC72_struc(n)%ac72(t)%WaterTableInProfile = GetWaterTableInProfile()
            AC72_struc(n)%ac72(t)%StartMode = GetStartMode()
            AC72_struc(n)%AC72(t)%NoMoreCrop = GetNoMoreCrop()

            ! Check for end of simulation period 
            ! (DayNri - 1 because DayNri is already for next day)
            if ((GetDayNri()-1) .eq. GetSimulation_ToDayNr()) then
                AC72_struc(n)%ac72(t)%InitializeRun = 1
                AC72_struc(n)%ac72(t)%read_Trecord = 1
                AC72_struc(n)%ac72(t)%irun = AC72_struc(n)%ac72(t)%irun + 1
            end if

            ! Diagnostic output variables
            ![ 1] output variable: smc (unit=m^3 m-3 ). ***  volumetric soil moisture
            do i=1, AC72_struc(n)%ac72(t)%NrCompartments
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = AC72_struc(n)%ac72(t)%smc(i), &
                                                    vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            ![ 2] output variable: biomass (unit=t/ha).  *** cummulative biomass
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Biomass, value = AC72_struc(n)%ac72(t)%SumWaBal%Biomass, &
                                                vlevel=1, unit="t ha-1", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 3] output variable: biomass (unit=mm).  *** actual rootzone water content
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_RootZoneWC_Actual, value = AC72_struc(n)%ac72(t)%RootZoneWC_Actual, &
                                                vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 4] output variable: RootZoneWC_WP (unit=mm).  *** rootzone water content at wilting point
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_RootZoneWC_WP, value = AC72_struc(n)%ac72(t)%RootZoneWC_WP, &
                                                vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 5] output variable: RootZoneWC_FC (unit=mm).  *** rootzone water content at field capacity
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_RootZoneWC_FC, value = AC72_struc(n)%ac72(t)%RootZoneWC_FC, &
                                                vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 6] output variable: Tact (unit=mm).  *** actual transpiration
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Tact, value = AC72_struc(n)%ac72(t)%Tact, &
                                                vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 7] output variable: Eact (unit=mm).  *** actual evaporation
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Eact, value = AC72_struc(n)%ac72(t)%Eact, &
                                                vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 8] output variable: AC72ETo (unit=mm).  *** reference evapotranspiration (Penman-Monteith)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_ETo, value = AC72_struc(n)%ac72(t)%eto, &
                                                vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 9] output variable: RootingDepth (unit=m).  *** rooting depth
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_RootingDepth, value = AC72_struc(n)%ac72(t)%RootingDepth, &
                                                vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 10] output variable: CCiActual (unit=-).  *** canopy cover
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_CCiActual, value = AC72_struc(n)%ac72(t)%CCiActual, &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 11] output variable: AC72Tmin (unit=deg C).  *** daily minimum temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Tmin, value = tmp_tmin, &
                                    vlevel=1, unit="degC", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 12] output variable: AC72Tmax (unit=deg C).  *** daily maximum temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Tmax, value = tmp_tmax, &
                                    vlevel=1, unit="degC", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 13] output variable: AC72Rain (unit=mm).  *** precipitation rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Rain, value = tmp_precip, &
                                    vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 13] output variable: yield (unit=t ha-1).  *** yield
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Yield, value = AC72_struc(n)%ac72(t)%SumWaBal%YieldPart, &
                                    vlevel=1, unit="t ha-1", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 14] output variable: irrigation (unit=mm).  *** irrigation
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_Irrigation, value = AC72_struc(n)%ac72(t)%Irrigation, &
                                    vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 15] output variable: StExp (unit=%).  *** expansion stress
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_StExp, value = GetStressLeaf(), &
                                    vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 16] output variable: StSen (unit=%).  *** senescence stress
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_StSen, value = GetStressSenescence(), &
                                    vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 17] output variable: cycle_complete (unit=binary).  *** Flag for completion of crop cycle within sim period
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC_cycle_complete, value = real(AC72_struc(n)%ac72(t)%cycle_complete), &
                                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            !  Reset forcings
            AC72_struc(n)%ac72(t)%tair = 0.0
            AC72_struc(n)%ac72(t)%tmax = 0.0
            AC72_struc(n)%ac72(t)%tmin = 0.0
            AC72_struc(n)%ac72(t)%tdew = 0.0
            AC72_struc(n)%ac72(t)%wndspd = 0.0            
            AC72_struc(n)%ac72(t)%psurf = 0.0
            AC72_struc(n)%ac72(t)%prcp = 0.0
            AC72_struc(n)%ac72(t)%eto = 0.0
            AC72_struc(n)%ac72(t)%swdown = 0.0
        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        AC72_struc(n)%forc_count = 0
    endif ! end of alarmCheck loop 

end subroutine AC72_main
