!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: Ac71_main
! \label{Ac71_main}
!
! !REVISION HISTORY:
!   18 JAN 2024, Louise Busschaert; initial implementation
!
! !INTERFACE:
subroutine Ac71_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod
    use LIS_constantsMod
    use Ac71_lsmMod
    use ac71_prep_f
    use ESMF

    !!! MB_AC70
    use ac_global, only: typeproject_typeprm, &
                        GetCrop_ModeCycle,&
                        ModeCycle_GDDays,&
                         typeproject_typepro, &
                         GetSimulParam_ThicknessTopSWC, &
                         GetRootZoneWC_Actual,&
                         GetRootZoneWC_FC,&
                         GetRootZoneWC_WP,&
                         GetRootZoneWC_SAT,&
                         GetRootZoneWC_Leaf,&
                         GetRootZoneWC_Thresh,&
                         GetRootZoneWC_Sen,&
                         GetRootZoneWC_ZtopAct,&
                         GetRootZoneWC_ZtopFC,&
                         GetRootZoneWC_ZtopWP,&
                         GetRootZoneWC_ZtopThresh,&
                         SetSimulParam_ThicknessTopSWC, &
                         SetRootZoneWC_Actual,&
                         SetRootZoneWC_FC,&
                         SetRootZoneWC_WP,&
                         SetRootZoneWC_SAT,&
                         SetCrop_DayN,&
                         SetRootZoneWC_Leaf,&
                         SetRootZoneWC_Thresh,&
                         SetRootZoneWC_Sen,&
                         SetRootZoneWC_ZtopAct,&
                         SetRootZoneWC_ZtopFC,&
                         SetRootZoneWC_ZtopWP,&
                         SetRootZoneWC_ZtopThresh,&
                         GetCompartment,&
                         SetCompartment,&
                         GetSoilLayer,&
                         SetSoilLayer,&
                         GetTotalSaltContent,&
                         GetTotalWaterContent,&
                         Geteffectiverain,&
                         GetSumWaBal,&
                         GetRootZoneSalt,&
                         GetSimulation,&
                         GetSimulation_SumGDD,&
                         GetSimulation_SumGDDfromDay1,&
                         SetTotalSaltContent,&
                         SetTotalWaterContent,&
                         Seteffectiverain,&
                         SetSumWaBal,&
                         SetRootZoneSalt,&
                         SetSimulation,&
                         GetIrrigation,&
                         SetIrrigation,&
                         GetCompartment_theta,&
                         SetCompartment_theta,&
                         GetIrriECw,&
                         GetManagement,&
                         GetPerennialPeriod,&
                         GetSimulParam,&
                         GetManagement_Cuttings,&
                         GetOnset,&
                         GetEndSeason,&
                         GetCrop,&
                         GetSoil,&
                         GetTemperatureRecord,&
                         GetClimRecord,&
                         GetRainRecord,&
                         GetEToRecord,&
                         SetIrriECw,&
                         SetManagement,&
                         SetPerennialPeriod,&
                         SetSimulParam,&
                         SetManagement_Cuttings,&
                         SetOnset,&
                         SetEndSeason,&
                         SetCrop,&
                         SetSoil,&
                         SetTemperatureRecord,&
                         SetClimRecord,&
                         SetRainRecord,&
                         SetEToRecord,&
                         GetSimulParam_GDDMethod,&

                         GetGenerateTimeMode,&
                        GetGenerateDepthMode,&
                        GetIrriMode,&
                        GetIrriMethod,&
                        GetDaySubmerged,&
                        GetMaxPlotNew,&
                        GetNrCompartments,&
                        GetIrriFirstDayNr,&
                        GetZiAqua,&
                        GetIniPercTAW,&
                        GetMaxPlotTr,&
                        GetOutputAggregate,&
                        GetEvapoEntireSoilSurface,&
                        GetPreDay,&
                        GetOutDaily,&
                        GetOut1Wabal,&
                        GetOut2Crop,&
                        GetOut3Prof,&
                        GetOut4Salt,&
                        GetOut5CompWC,&
                        GetOut6CompEC,&
                        GetOut7Clim,&
                        GetPart1Mult,&
                        GetPart2Eval,&
                        SetGenerateTimeMode,&
                        SetGenerateDepthMode,&
                        SetIrriMode,&
                        SetIrriMethod,&
                        SetDaySubmerged,&
                        SetMaxPlotNew,&
                        SetNrCompartments,&
                        SetIrriFirstDayNr,&
                        SetZiAqua,&
                        SetIniPercTAW,&
                        SetMaxPlotTr,&
                        SetOutputAggregate,&
                        SetEvapoEntireSoilSurface,&
                        SetPreDay,&
                        SetOutDaily,&
                        SetOut1Wabal,&
                        SetOut2Crop,&
                        SetOut3Prof,&
                        SetOut4Salt,&
                        SetOut5CompWC,&
                        SetOut6CompEC,&
                        SetOut7Clim,&
                        SetPart1Mult,&
                        SetPart2Eval,&


                        GetCCiActual,&
                        GetCCiprev,&
                        GetCCiTopEarlySen,&
                        GetCRsalt,&  
                        GetCRwater,& 
                        GetECdrain,& 
                        GetECiAqua,& 
                        GetECstorage,& 
                        GetEact,& 
                        GetEpot,& 
                        GetETo,&
                        GetDrain,&  
                        GetInfiltrated,&
                        GetRain,& 
                        GetRootingDepth,&
                        GetRunoff,& 
                        GetSaltInfiltr,&
                        GetSurf0,& 
                        GetSurfaceStorage,&
                        GetTact,&
                        GetTpot,&
                        GetTactWeedInfested,&
                        GetTmax,& 
                        GetTmin,& 
                        SetCCiActual,&
                        SetCCiprev,&
                        SetCCiTopEarlySen,&
                        SetCRsalt,&  
                        SetCRwater,& 
                        SetECdrain,& 
                        SetECiAqua,& 
                        SetECstorage,& 
                        SetEact,& 
                        SetEpot,& 
                        SetETo,&
                        SetDrain,&  
                        SetInfiltrated,&
                        SetRain,& 
                        SetRootingDepth,&
                        SetRunoff,& 
                        SetSaltInfiltr,&
                        SetSurf0,& 
                        SetSurfaceStorage,&
                        SetTact,&
                        SetTpot,&
                        SetTactWeedInfested,&
                        SetTmax,& 
                        SetTmin,&
                        GetIrriBeforeSeason,&
                        SetIrriBeforeSeason,&
                        GetIrriAfterSeason,&
                        SetIrriAfterSeason,&
                        GetCrop_Day1,&
                        DegreesDay,&
                         GetSimulation_ToDayNr, &
                         SetSimulation_ToDayNr, &
                        SetSimulation_SumGDDfromDay1,&
                        SetSimulation_SumGDD, &
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
                        SetClimFile, &
                        GetCrop_DaysToGermination,&
                        GetCrop_DaysToMaxRooting,&
                        GetCrop_DaysToFlowering,&
                        GetCrop_DaysToHarvest,&
                        GetCrop_DaysTosenescence,&
                        GetCrop_DaysToCCini,&
                        GetCrop_DaysToFullCanopy,&
                        GetCrop_DaysToFullCanopySF,&
                        GetCrop_DaysToHIo,&
                    GetSimulation_EvapLimitON, &
                    GetSimulation_SumGDD,&
                    GetSimulation_SumGDDfromDay1,&
                    GetSimulation_SWCtopSoilConsidered, &
                    GetSimulation_FromDayNr, &
                    GetSimulation_ToDayNr,&
                    IrriMode_Generate,&
                    IrriMode_Inet,&
                    IrriMode_Manual,&
                    IrriMode_NoIrri,&
                    undef_int, &
                    GetCrop_GDDaysToGermination,&
                    GetCrop_GDDaysToMaxRooting,&
                    GetCrop_GDDaysToHarvest,&
                    GetCrop_GDDaysToGermination,&
                    GetCrop_GDDaysToFlowering,&
                    GetCrop_GDDaysToSenescence,&
                    SetCrop_DaysToGermination,&
                    SetCrop_DaysToMaxRooting,&
                    SetCrop_DaysToFlowering,&
                    SetCrop_DaysToHarvest,&
                    SetCrop_DaysTosenescence,&
                    SetCrop_DaysToCCini,&
                    SetCrop_DaysToFullCanopy,&
                    SetCrop_DaysToFullCanopySF,&
                    SetCrop_DaysToHIo,&
                    SetTminRun_i, &
                    SetTminRun, &
                    SetTmaxRun_i, &
                    SetTmaxRun, &
                    GetTminRun_i, &
                    GetTminRun, &
                    GetTmaxRun_i, &
                    GetTmaxRun, &
                    GetCropFileSet, &
                    SetCropFileSet

           !!! MB_AC70
    !!! MB:
    use ac_project_input, only: ProjectInput 

    use ac_run, only:    SetDayNri,&
                         GetIrriInterval,&
                         GetIrriInfoRecord1,&
                         GetIrriInfoRecord2,&
                         SetIrriInterval,&
                         SetIrriInfoRecord1,&
                         SetIrriInfoRecord2,&
                         GetTheProjectFile,&
                    fIrri_close, &
                    fIrri_open,&
                    firri_read,&
                        GetGwTable,&
                        GetPlotVarCrop,&
                        GetStressTot,&
                        GetCutInfoRecord1,&
                        GetCutInfoRecord2,&
                        GetTransfer,&
                        GetPreviousSum,&
                        GetTadj,&
                        GetGDDTadj,&
                        GetDayLastCut,&
                        GetNrCut,&
                        GetSumInterval,&
                        GetPreviousStressLevel,&
                        GetStressSFadjNEW,&
                        GetBin,&
                        GetBout,&
                        GetGDDayi,&
                        GetCO2i,&
                        GetFracBiomassPotSF,&
                        GetSumETo,&
                        GetSumGDD,&
                        GetZiprev,&
                        GetSumGDDPrev,&
                        GetCCxWitheredTpotNoS,&
                        GetCoeffb0,&
                        GetCoeffb1,&
                        GetCoeffb2,&
                        GetCoeffb0Salt,&
                        GetCoeffb1Salt,&
                        GetCoeffb2Salt,&
                        GetStressLeaf,&
                        GetStressSenescence ,&
                        GetDayFraction,&
                        GetGDDayFraction,&
                        GetCGCref,&
                        GetGDDCGCref ,&
                        GetTimeSenescence ,&
                        GetSumKcTop,&
                        GetSumKcTopStress,&
                        GetSumKci,&
                        GetCCoTotal,&
                        GetCCxTotal,&
                        GetCDCTotal,&
                        GetGDDCDCTotal,&
                        GetCCxCropWeedsNoSFstress,&
                        GetWeedRCi,&
                        GetCCiActualWeedInfested,&
                        GetfWeedNoS,&
                        GetZeval,&
                        GetBprevSum,&
                        GetYprevSum,&
                        GetSumGDDcuts,&
                        GetHItimesBEF,&
                        GetScorAT1,&
                        GetScorAT2,&
                        GetHItimesAT1,&
                        GetHItimesAT2,&
                        GetHItimesAT,&
                        GetalfaHI,&
                        GetalfaHIAdj,&
                        GetNextSimFromDayNr ,&
                        GetDayNr1Eval,&
                        GetDayNrEval,&
                        GetLineNrEval,&
                        GetPreviousSumETo,&
                        GetPreviousSumGDD,&
                        GetPreviousBmob,&
                        GetPreviousBsto,&
                        GetStageCode,&
                        GetPreviousDayNr,&
                        GetNoYear,&
                        GetWaterTableInProfile,&
                        GetStartMode,&
                        GetNoMoreCrop,&
                        SetGwTable,&
                        SetPlotVarCrop,&
                        SetStressTot,&
                        SetCutInfoRecord1,&
                        SetCutInfoRecord2,&
                        SetTransfer,&
                        SetPreviousSum,&
                        SetTadj,&
                        SetGDDTadj,&
                        SetDayLastCut,&
                        SetNrCut,&
                        SetSumInterval,&
                        SetPreviousStressLevel,&
                        SetStressSFadjNEW,&
                        SetBin,&
                        SetBout,&
                        SetGDDayi,&
                        SetCO2i,&
                        SetFracBiomassPotSF,&
                        SetSumETo,&
                        SetSumGDD,&
                        SetZiprev,&
                        SetSumGDDPrev,&
                        SetCCxWitheredTpotNoS,&
                        SetCoeffb0,&
                        SetCoeffb1,&
                        SetCoeffb2,&
                        SetCoeffb0Salt,&
                        SetCoeffb1Salt,&
                        SetCoeffb2Salt,&
                        SetStressLeaf,&
                        SetStressSenescence ,&
                        SetDayFraction,&
                        SetGDDayFraction,&
                        SetCGCref,&
                        SetGDDCGCref ,&
                        SetTimeSenescence ,&
                        SetSumKcTop,&
                        SetSumKcTopStress,&
                        SetSumKci,&
                        SetCCoTotal,&
                        SetCCxTotal,&
                        SetCDCTotal,&
                        SetGDDCDCTotal,&
                        SetCCxCropWeedsNoSFstress,&
                        SetWeedRCi,&
                        SetCCiActualWeedInfested,&
                        SetfWeedNoS,&
                        SetZeval,&
                        SetBprevSum,&
                        SetYprevSum,&
                        SetSumGDDcuts,&
                        SetHItimesBEF,&
                        SetScorAT1,&
                        SetScorAT2,&
                        SetHItimesAT1,&
                        SetHItimesAT2,&
                        SetHItimesAT,&
                        SetalfaHI,&
                        SetalfaHIAdj,&
                        SetNextSimFromDayNr ,&
                        SetDayNr1Eval,&
                        SetDayNrEval,&
                        SetLineNrEval,&
                        SetPreviousSumETo,&
                        SetPreviousSumGDD,&
                        SetPreviousBmob,&
                        SetPreviousBsto,&
                        SetStageCode,&
                        SetPreviousDayNr,&
                        SetNoYear,&
                        SetWaterTableInProfile,&
                        SetStartMode,&
                        SetNoMoreCrop,&
                        AdvanceOneTimeStep, &
                        ReadClimateNextDay, &
                        SetGDDVariablesNextDay, &
                         FinalizeRun1, &
                         FinalizeRun2, &
                         GetDayNri,&
                         GetCrop_Tbase, &
                         GetCrop_Tupper, &
                         FinalizeSimulation, &
                         InitializeSimulation, &
                         InitializeRunPart1, &
                         InitializeRunPart2, &
                         InitializeSimulationRunPart2, &
                         InitializeClimate,&
                        GetIrriInfoRecord1, &
                        GetIrriInfoRecord1_NoMoreInfo,&
                        GetIrriInfoRecord2, &
                        fIrri_read

! From ac_project_input
   use ac_project_input, only: ProjectInput, set_project_input
! From ac_kinds
    use ac_kinds, only: intEnum, &
                        int32, &
                        int8, &
                        dp,&
                        sp

! From ac_startunit
    use ac_startunit, only:  FinalizeTheProgram, &
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

    !!! MB_AC71
    integer              :: daynr, todaynr, iproject, nprojects
    logical              :: ListProjectFileExist, phenological_stages_ensemble
    character(len=:), allocatable :: ListProjectsFile, TheProjectFile
    integer              :: Crop_DaysToGermination, Crop_DaysToMaxRooting, Crop_DaysToFlowering
    integer              :: Crop_DaysToHarvest, Crop_DaysTosenescence, Crop_DaysToCCini
    integer              :: Crop_DaysToFullCanopy, Crop_DaysToFullCanopySF, Crop_DaysToHIo
    integer(int32) :: temp1    

    !LB AC71
    integer              :: irr_record_flag, DNr ! for irri file management
    character(250)       :: TempStr 

    real                 :: tmp_pres, tmp_precip, tmp_tmax, tmp_tmin   ! Weather Forcing
    real                 :: tmp_tdew, tmp_swrad, tmp_wind, tmp_eto     ! Weather Forcing

    ! For type problem in AdvanceOneTimeStep
    real(dp)             :: tmp_wpi
!
! !DESCRIPTION:
!  This is the entry point for calling the Ac71 physics.
!  This routine calls the {\AdvanceOneTimeStep} routine that runs AquaCrop

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for Ac71
    ! check Ac71 alarm. If alarm is ring, run model.
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "Ac71 model alarm")
    if (alarmCheck) Then
        if (AC71_struc(n)%ac71(1)%InitializeRun.eq.1) then
            call ac71_read_Trecord(n)
        endif
        do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
            dt = LIS_rc%ts
            row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
            col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
            lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
            lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon
            tmp_elev = LIS_domain(n)%tile(t)%elev

            !!------ This Block Is Where We Obtain Weather Forcing ------------------------------!!
            ! retrieve forcing data from AC71_struc(n)%ac71(t) and assign to local variables

            ! PRES: Daily average surface pressure (kPa)
            tmp_pres      = (AC71_struc(n)%ac71(t)%psurf / Ac71_struc(n)%forc_count) / 1000

            ! PRECIP: Total daily precipitation (rain+snow) (mm)
            tmp_precip    = (AC71_struc(n)%ac71(t)%prcp / Ac71_struc(n)%forc_count) * 3600. * 24. !Convert from kg/ms2 to mm/d

            ! TMAX: maximum daily air temperature (degC)
            tmp_tmax      = AC71_struc(n)%ac71(t)%tmax - LIS_CONST_TKFRZ !Convert from K to C

            ! TMIN: minimum daily air temperature (degC)
            tmp_tmin      = AC71_struc(n)%ac71(t)%tmin - LIS_CONST_TKFRZ !Convert from K to C 

            ! TDEW: average daily dewpoint temperature (degC)
            tmp_tdew      = (AC71_struc(n)%ac71(t)%tdew / AC71_struc(n)%forc_count) - LIS_CONST_TKFRZ !Convert from K to C

            ! SW_RAD: daily total incoming solar radiation (MJ/(m2d))
            tmp_swrad     = (AC71_struc(n)%ac71(t)%swdown / AC71_struc(n)%forc_count) * 0.0864 !Convert from W/m2 to MJ/(m2d)

            ! Wind: daily average wind speed (m/s)
            tmp_wind      = AC71_struc(n)%ac71(t)%wndspd / AC71_struc(n)%forc_count

            ! check validity of PRES
            if(tmp_pres .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable PRES in Ac71"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! check validity of PRECIP
            if(tmp_precip .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable PRECIP in Ac71"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! check validity of TMAX
            if(tmp_tmax .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable TMAX in Ac71"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

           ! check validity of TMIN
            if(tmp_tmin .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable TMIN in Ac71"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! check validity of TDEW
            if(tmp_tdew .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable TDEW in Ac71"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! check validity of SW_RAD
            if(tmp_swrad .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable SW_RAD in Ac71"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif


            ! check validity of WIND
            if(tmp_wind .eq. LIS_rc%udef) then
                write(LIS_logunit, *) "undefined value found for forcing variable WIND in Ac71"
                write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
                call LIS_endrun()
            endif

            ! Call ETo subroutine           
            call ac71_ETo_calc(tmp_pres, tmp_tmax, tmp_tmin, tmp_tdew, &
                               tmp_wind, tmp_swrad, &
                               tmp_elev, lat, tmp_eto)
            AC71_struc(n)%ac71(t)%eto = tmp_eto

            ! setting all global variables
            call SetRootZoneWC_Actual(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_Actual,8))
            call SetRootZoneWC_FC(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_FC,8))
            call SetRootZoneWC_WP(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_WP,8))
            call SetRootZoneWC_SAT(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_SAT,8))
            call SetRootZoneWC_Leaf(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_Leaf,8))
            call SetRootZoneWC_Thresh(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_Thresh,8))
            call SetRootZoneWC_Sen(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_Sen,8))
            call SetRootZoneWC_ZtopAct(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopFC(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_ZtopFC,8))
            call SetRootZoneWC_ZtopWP(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_ZtopWP,8))
            call SetRootZoneWC_ZtopThresh(REAL(AC71_struc(n)%AC71(t)%RootZoneWC_ZtopThresh,8))
            call SetCompartment(AC71_struc(n)%AC71(t)%Compartment)
            call SetTotalSaltContent(AC71_struc(n)%AC71(t)%TotalSaltContent)
            call SetTotalWaterContent(AC71_struc(n)%AC71(t)%TotalWaterContent)
            call Seteffectiverain(AC71_struc(n)%AC71(t)%effectiverain)
            call SetSumWaBal(AC71_struc(n)%AC71(t)%SumWaBal)
            call SetRootZoneSalt(AC71_struc(n)%AC71(t)%RootZoneSalt)
            call SetSimulation(AC71_struc(n)%AC71(t)%Simulation)
            call SetIrriInterval(AC71_struc(n)%AC71(t)%IrriInterval)
            call SetIrriInfoRecord1(AC71_struc(n)%AC71(t)%IrriInfoRecord1)
            call SetIrriInfoRecord2(AC71_struc(n)%AC71(t)%IrriInfoRecord2)
            call SetIrrigation(REAL(AC71_struc(n)%AC71(t)%Irrigation,8))
            do l=1, AC71_struc(n)%AC71(t)%NrCompartments
                 call SetCompartment_theta(l,REAL(AC71_struc(n)%AC71(t)%smc(l),8))
            enddo
            call SetIrriECw(AC71_struc(n)%AC71(t)%IrriECw) 
            call SetManagement(AC71_struc(n)%AC71(t)%Management) 
            call SetPerennialPeriod(AC71_struc(n)%AC71(t)%PerennialPeriod) 
            call SetSimulParam(AC71_struc(n)%AC71(t)%simulparam) 
            call SetManagement_Cuttings(AC71_struc(n)%AC71(t)%Cuttings) 
            call SetOnset(AC71_struc(n)%AC71(t)%onset) 
            call SetEndSeason(AC71_struc(n)%AC71(t)%endseason) 
            call SetCrop(AC71_struc(n)%AC71(t)%crop) 
            call SetSoil(AC71_struc(n)%AC71(t)%Soil) 
            call SetIrriBeforeSeason(AC71_struc(n)%AC71(t)%IrriBeforeSeason)
            call SetIrriAfterSeason(AC71_struc(n)%AC71(t)%IrriAfterSeason)
            call SetSoilLayer(AC71_struc(n)%AC71(t)%soillayer)
            call SetDayNri(AC71_struc(n)%AC71(t)%daynri)

            call SetGenerateTimeMode(AC71_struc(n)%AC71(t)%GenerateTimeMode) 
            call SetGenerateDepthMode(AC71_struc(n)%AC71(t)%GenerateDepthMode) 
            call SetIrriMode(AC71_struc(n)%AC71(t)%IrriMode) 
            call SetDaySubmerged(AC71_struc(n)%AC71(t)%DaySubmerged) 
            call SetMaxPlotNew(AC71_struc(n)%AC71(t)%MaxPlotNew) 
            call SetNrCompartments(AC71_struc(n)%AC71(t)%NrCompartments) 
            call SetIrriFirstDayNr(AC71_struc(n)%AC71(t)%IrriFirstDayNr) 
            call SetZiAqua(AC71_struc(n)%AC71(t)%ZiAqua) 
            call SetIniPercTAW(AC71_struc(n)%AC71(t)%IniPercTAW) 
            call SetMaxPlotTr(AC71_struc(n)%AC71(t)%MaxPlotTr)
            call SetOutputAggregate(AC71_struc(n)%AC71(t)%OutputAggregate) 

            call SetEvapoEntireSoilSurface(AC71_struc(n)%AC71(t)%EvapoEntireSoilSurface) 
            call SetPreDay(AC71_struc(n)%AC71(t)%PreDay) 
            call SetOutDaily(AC71_struc(n)%AC71(t)%OutDaily) 
            call SetOut1Wabal(AC71_struc(n)%AC71(t)%Out1Wabal) 
            call SetOut2Crop(AC71_struc(n)%AC71(t)%Out2Crop) 
            call SetOut3Prof(AC71_struc(n)%AC71(t)%Out3Prof) 
            call SetOut4Salt(AC71_struc(n)%AC71(t)%Out4Salt) 
            call SetOut5CompWC(AC71_struc(n)%AC71(t)%Out5CompWC) 
            call SetOut6CompEC(AC71_struc(n)%AC71(t)%Out6CompEC) 
            call SetOut7Clim(AC71_struc(n)%AC71(t)%Out7Clim) 
            call SetPart1Mult(AC71_struc(n)%AC71(t)%Part1Mult) 
            call SetPart2Eval(AC71_struc(n)%AC71(t)%Part2Eval) 

            !
            call SetCCiActual(AC71_struc(n)%AC71(t)%CCiActual)
            call SetCCiprev(AC71_struc(n)%AC71(t)%CCiprev)

            call SetCCiTopEarlySen(AC71_struc(n)%AC71(t)%CCiTopEarlySen)
            call SetCRsalt(AC71_struc(n)%AC71(t)%CRsalt)
            call SetCRwater(AC71_struc(n)%AC71(t)%CRwater)
            call SetECdrain(AC71_struc(n)%AC71(t)%ECdrain)
            call SetEciAqua(AC71_struc(n)%AC71(t)%ECiAqua)
            call SetECstorage(AC71_struc(n)%AC71(t)%ECstorage)
            call SetEact(AC71_struc(n)%AC71(t)%Eact)
            call SetEpot(AC71_struc(n)%AC71(t)%Epot)
            call SetDrain(AC71_struc(n)%AC71(t)%Drain)
            call SetInfiltrated(AC71_struc(n)%AC71(t)%Infiltrated)
            call SetRootingDepth(AC71_struc(n)%AC71(t)%RootingDepth)
            call SetRunoff(AC71_struc(n)%AC71(t)%Runoff)
            call SetSaltInfiltr(AC71_struc(n)%AC71(t)%SaltInfiltr)
            call SetSurf0(AC71_struc(n)%AC71(t)%Surf0)
            call SetSurfaceStorage(AC71_struc(n)%AC71(t)%SurfaceStorage)
            call SetTact(AC71_struc(n)%AC71(t)%Tact)
            call SetTpot(AC71_struc(n)%AC71(t)%Tpot)
            call SetTactWeedInfested(AC71_struc(n)%AC71(t)%TactWeedInfested)

            call SetGwTable(AC71_struc(n)%AC71(t)%GwTable)
            call SetPlotVarCrop(AC71_struc(n)%AC71(t)%PlotVarCrop)
            call SetStressTot(AC71_struc(n)%AC71(t)%StressTot)
            call SetCutInfoRecord1(AC71_struc(n)%AC71(t)%CutInfoRecord1)
            call SetCutInfoRecord2(AC71_struc(n)%AC71(t)%CutInfoRecord2)
            call SetTransfer(AC71_struc(n)%AC71(t)%Transfer)
            call SetPreviousSum(AC71_struc(n)%AC71(t)%PreviousSum)
            call SetTadj(AC71_struc(n)%AC71(t)%Tadj)
            call SetGDDTadj(AC71_struc(n)%AC71(t)%GDDTadj)
            call SetDayLastCut(AC71_struc(n)%AC71(t)%DayLastCut)
            call SetNrCut(AC71_struc(n)%AC71(t)%NrCut)
            call SetSumInterval(AC71_struc(n)%AC71(t)%SumInterval)
            call SetPreviousStressLevel(int(AC71_struc(n)%AC71(t)%PreviousStressLevel,kind=int32))
            call SetStressSFadjNEW(int(AC71_struc(n)%AC71(t)%StressSFadjNEW,kind=int32))
            call SetBin(AC71_struc(n)%AC71(t)%Bin)
            call SetBout(AC71_struc(n)%AC71(t)%Bout)
            call SetCO2i(AC71_struc(n)%AC71(t)%CO2i)
            call SetFracBiomassPotSF(AC71_struc(n)%AC71(t)%FracBiomassPotSF)
            call SetSumETo(AC71_struc(n)%AC71(t)%SumETo)
            call SetSumGDD(AC71_struc(n)%AC71(t)%SumGDD)
            call SetZiprev(AC71_struc(n)%AC71(t)%Ziprev)
            call SetSumGDDPrev(AC71_struc(n)%AC71(t)%SumGDDPrev)
            call SetCCxWitheredTpotNoS(AC71_struc(n)%AC71(t)%CCxWitheredTpotNoS)
            call SetCoeffb0(AC71_struc(n)%AC71(t)%Coeffb0)
            call SetCoeffb1(AC71_struc(n)%AC71(t)%Coeffb1)
            call SetCoeffb2(AC71_struc(n)%AC71(t)%Coeffb2)
            call SetCoeffb0Salt(AC71_struc(n)%AC71(t)%Coeffb0Salt)
            call SetCoeffb1Salt(AC71_struc(n)%AC71(t)%Coeffb1Salt)
            call SetCoeffb2Salt(AC71_struc(n)%AC71(t)%Coeffb2Salt)
            call SetStressLeaf(AC71_struc(n)%AC71(t)%StressLeaf)
            call SetStressSenescence(AC71_struc(n)%AC71(t)%StressSenescence)
            call SetDayFraction(AC71_struc(n)%AC71(t)%DayFraction)
            call SetGDDayFraction(AC71_struc(n)%AC71(t)%GDDayFraction)
            call SetCGCref(AC71_struc(n)%AC71(t)%CGCref)
            call SetGDDCGCref(AC71_struc(n)%AC71(t)%GDDCGCref)
            call SetTimeSenescence(AC71_struc(n)%AC71(t)%TimeSenescence)
            call SetSumKcTop(AC71_struc(n)%AC71(t)%SumKcTop)
            call SetSumKcTopStress(AC71_struc(n)%AC71(t)%SumKcTopStress)
            call SetSumKci(AC71_struc(n)%AC71(t)%SumKci)
            call SetCCoTotal(AC71_struc(n)%AC71(t)%CCoTotal)
            call SetCCxTotal(AC71_struc(n)%AC71(t)%CCxTotal)
            call SetCDCTotal(AC71_struc(n)%AC71(t)%CDCTotal)
            call SetGDDCDCTotal(AC71_struc(n)%AC71(t)%GDDCDCTotal)
            call SetCCxCropWeedsNoSFstress(AC71_struc(n)%AC71(t)%CCxCropWeedsNoSFstress)
            call SetWeedRCi(AC71_struc(n)%AC71(t)%WeedRCi)
            call SetCCiActualWeedInfested(AC71_struc(n)%AC71(t)%CCiActualWeedInfested)
            call SetfWeedNoS(AC71_struc(n)%AC71(t)%fWeedNoS)
            call SetZeval(AC71_struc(n)%AC71(t)%Zeval)
            call SetBprevSum(AC71_struc(n)%AC71(t)%BprevSum)
            call SetYprevSum(AC71_struc(n)%AC71(t)%YprevSum)
            call SetSumGDDcuts(AC71_struc(n)%AC71(t)%SumGDDcuts)
            call SetHItimesBEF(AC71_struc(n)%AC71(t)%HItimesBEF)
            call SetScorAT1(AC71_struc(n)%AC71(t)%ScorAT1)
            call SetScorAT2(AC71_struc(n)%AC71(t)%ScorAT2)
            call SetHItimesAT1(AC71_struc(n)%AC71(t)%HItimesAT1)
            call SetHItimesAT2(AC71_struc(n)%AC71(t)%HItimesAT2)
            call SetHItimesAT(AC71_struc(n)%AC71(t)%HItimesAT)
            call SetalfaHI(AC71_struc(n)%AC71(t)%alfaHI)
            call SetalfaHIAdj(AC71_struc(n)%AC71(t)%alfaHIAdj)
            call SetNextSimFromDayNr(AC71_struc(n)%AC71(t)%NextSimFromDayNr)
            call SetDayNr1Eval(AC71_struc(n)%AC71(t)%DayNr1Eval)
            call SetDayNrEval(AC71_struc(n)%AC71(t)%DayNrEval)
            call SetLineNrEval(int(AC71_struc(n)%AC71(t)%LineNrEval,kind=int32))
            call SetPreviousSumETo(AC71_struc(n)%AC71(t)%PreviousSumETo)
            call SetPreviousSumGDD(AC71_struc(n)%AC71(t)%PreviousSumGDD)
            call SetPreviousBmob(AC71_struc(n)%AC71(t)%PreviousBmob)
            call SetPreviousBsto(AC71_struc(n)%AC71(t)%PreviousBsto)
            call SetStageCode(AC71_struc(n)%AC71(t)%StageCode)
            call SetPreviousDayNr(AC71_struc(n)%AC71(t)%PreviousDayNr)
            call SetNoYear(AC71_struc(n)%AC71(t)%NoYear)
            call SetWaterTableInProfile(AC71_struc(n)%AC71(t)%WaterTableInProfile)
            call SetStartMode(AC71_struc(n)%AC71(t)%StartMode)
            call SetNoMoreCrop(AC71_struc(n)%AC71(t)%NoMoreCrop)
            !call SetSimulation_ToDayNr(AC71_struc(n)%AC71(t)%Simulation%ToDayNr)
            call SetGDDayi(AC71_struc(n)%AC71(t)%GDDayi)

            ! Fixed var
            call SetOut3Prof(.true.) ! needed for correct rootzone sm
            call SetOutDaily(.true.)

            ! Test
            call SetCropFileSet(AC71_struc(n)%ac71(t)%CropFileSet)
            call SetTminRun(AC71_struc(n)%ac71(t)%Tmin_record)
            call SetTmaxRun(AC71_struc(n)%ac71(t)%Tmax_record)
            

            ! Set climate variables
            ! Round them to 4 digits after the comma as done in the AC standalone
            tmp_precip = anint(tmp_precip*10000)/10000
            call SetRain(real(tmp_precip, kind=dp))
            tmp_tmin = anint(tmp_tmin*10000)/10000
            call SetTmin(real(tmp_tmin, kind=dp))
            tmp_tmax = anint(tmp_tmax*10000)/10000
            call SetTmax(real(tmp_tmax, kind=dp))
            tmp_eto = anint(tmp_eto*10000)/10000
            call SetETo(real(tmp_eto, kind=dp))
            
            ! SumGDD calculation needed only for second day when not done within InitializeSimulationRunPart2
            if (GetDayNri()>GetSimulation_FromDayNr()) then
                ! Sum of GDD at end of first day ! Wait for GDD implementation from Michel
                call SetGDDayi(DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), GetTmin(), &
                        GetTmax(), GetSimulParam_GDDMethod()))
                if (GetDayNri() >= GetCrop_Day1()) then
                    call SetSimulation_SumGDD(GetSimulation_SumGDD() + GetGDDayi())
                    call SetSimulation_SumGDDfromDay1(GetSimulation_SumGDDfromDay1() + &
                        GetGDDayi())
                end if
            end if

            ! Irrigation management
            call SetIrriMode(AC71_struc(n)%ac71(t)%IrriMode)
            irr_record_flag = 0
            if(AC71_struc(n)%ac71(t)%IrriMode.ne.IrriMode_NoIrri) then
                ! Irrigation ON for tile
                call SetIrriAfterSeason(AC71_struc(n)%ac71(t)%IrriAfterSeason)
                call SetIrriBeforeSeason(AC71_struc(n)%ac71(t)%IrriBeforeSeason)
                call SetIrriECw(AC71_struc(n)%ac71(t)%IrriECw) 
                call SetIrriInfoRecord1(AC71_struc(n)%ac71(t)%IrriInfoRecord1)
                call SetIrriInfoRecord2(AC71_struc(n)%ac71(t)%IrriInfoRecord2)
                call SetIrriMethod(AC71_struc(n)%ac71(t)%IrriMethod)
                call SetGenerateDepthMode(AC71_struc(n)%ac71(t)%GenerateDepthMode)
                call SetGenerateTimeMode(AC71_struc(n)%ac71(t)%GenerateTimeMode)

                ! For IrriMode_Generate and IrriMode_Manual
                ! Check if new record needs to be read
                if(AC71_struc(n)%ac71(t)%IrriMode.eq.IrriMode_Generate) then
                    if((AC71_struc(n)%ac71(t)%daynri-AC71_struc(n)%ac71(t)%Crop%Day1+1)&
                    .gt.AC71_struc(n)%ac71(t)%IrriInfoRecord1%ToDay) then
                        irr_record_flag = 1
                    endif
                elseif(AC71_struc(n)%ac71(t)%IrriMode.eq.IrriMode_Manual) then
                    ! Check start date of schedule
                    if(AC71_struc(n)%ac71(t)%IrriFirstDayNr.eq.undef_int)then
                        DNr = AC71_struc(n)%ac71(t)%daynri &
                              - AC71_struc(n)%ac71(t)%Crop%Day1 + 1
                    else
                        DNr = AC71_struc(n)%ac71(t)%daynri &
                              - AC71_struc(n)%ac71(t)%IrriFirstDayNr + 1
                    endif
                    if(AC71_struc(n)%ac71(t)%IrriInfoRecord1%TimeInfo.eq.DNr)then
                        irr_record_flag = 1
                    endif
                endif

                ! re-open irrigation file and read the previous records
                if(irr_record_flag.eq.1)then
                    call fIrri_open(trim(AC71_struc(n)%PathNameSimul)&
                                    //trim(AC71_struc(n)%Irrigation_Filename), 'r')
                        do i=1,AC71_struc(n)%ac71(t)%irri_lnr
                            TempStr = fIrri_read()
                        enddo
                    AC71_struc(n)%ac71(t)%irri_lnr = AC71_struc(n)%ac71(t)%irri_lnr + 1
                endif
            endif
            ! End irrigation block

            !!! initialize run (year)
            if (AC71_struc(n)%ac71(t)%InitializeRun.eq.1) then !make it flex
                call SetClimRecord_DataType(0_int8)
                call SetClimRecord_fromd(0)
                call SetClimRecord_fromdaynr(ProjectInput(1)%Simulation_DayNr1)
                call SetClimRecord_fromm(0)
                call SetClimRecord_fromstring("")
                call SetClimRecord_fromy(LIS_rc%syr)
                call SetClimRecord_NrObs(999)
                call SetClimRecord_tod(0)
                call SetClimRecord_todaynr(ProjectInput(GetSimulation_NrRuns())%Simulation_DayNrN)
                call SetClimRecord_tom(0)
                call SetClimRecord_tostring("")
                call SetClimRecord_toy(0)
                call SetClimFile('(External)')

                AC71_struc(n)%ac71(t)%WPi = 0._dp

                ! Set crop file (crop parameters are read when calling InitializeRunPart1)
                call set_project_input(AC71_struc(n)%ac71(t)%irun, &
                                       'Crop_Filename', &
                                        trim(AC71_struc(n)%ac71(t)%cropt)//'.CRO')

                ! Set Global variable to pass T record to AquaCrop
                call SetTminRun(AC71_struc(n)%ac71(t)%Tmin_record)
                call SetTmaxRun(AC71_struc(n)%ac71(t)%Tmax_record)

                ! InitializeRunPart
                call InitializeRunPart1(int(AC71_struc(n)%ac71(t)%irun,kind=int8), AC71_struc(n)%ac71(t)%TheProjectType)
                call InitializeSimulationRunPart2()
                AC71_struc(n)%ac71(t)%HarvestNow = .false. ! Initialize to false
                ! Check if enough GDDays to complete cycle
                if(GetCrop_ModeCycle().eq.ModeCycle_GDDays)then
                    if (GetCrop_DaysToHarvest().eq.undef_int) then
                        AC71_struc(n)%ac71(t)%Crop%ModeCycle = -9
                    endif
                endif

                ! Irrigaton file management after InitializeRun
                if(GetIrriMode().ne.IrriMode_NoIrri) then
                    call fIrri_close()
                    ! Check for irrigation (irrigation file management)
                    if(GetIrriMode().eq.IrriMode_Manual)then
                        if(GetIrriInfoRecord1_NoMoreInfo())then
                            AC71_struc(n)%ac71(t)%irri_lnr = 9
                        else
                            AC71_struc(n)%ac71(t)%irri_lnr = 10
                        endif
                    elseif(GetIrriMode().eq.IrriMode_Generate)then
                        if(AC71_struc(n)%ac71(t)%IrriInfoRecord1%NoMoreInfo)then
                            AC71_struc(n)%ac71(t)%irri_lnr = 11
                        else
                            AC71_struc(n)%ac71(t)%irri_lnr = 12
                        endif
                    else ! no irrigation, set to 0
                        AC71_struc(n)%ac71(t)%irri_lnr = 0
                    endif
                endif
                ! End irrigation block
                AC71_struc(n)%ac71(t)%InitializeRun = 0 ! Initialization done
            end if


            ! Initialize for new crop cycle
            ! Reset Crop%DaysTo* to allow that members reach stages at different days
            phenological_stages_ensemble = .false.

            if (phenological_stages_ensemble) then
                if (GetDayNri() == GetCrop_Day1()) then
                AC71_struc(n)%ac71(t)%germ_reached = .false.
                AC71_struc(n)%ac71(t)%harv_reached = .false.
                AC71_struc(n)%ac71(t)%flowr_reached = .false.
                AC71_struc(n)%ac71(t)%MaxR_reached = .false.
                AC71_struc(n)%ac71(t)%Sene_reached = .false.

                !find calendar days for crop stages
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToGermination()) &
                    .and.  (.not. AC71_struc(n)%ac71(t)%germ_reached)) then ! from sow
                    AC71_struc(n)%ac71(t)%Crop%DaysToGermination = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%germ_reached = .true.
                end if
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToMaxRooting()) &
                    .and. (.not. AC71_struc(n)%ac71(t)%maxR_reached)) then ! from sowin
                    AC71_struc(n)%ac71(t)%Crop%DaysToMaxRooting = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%maxR_reached = .true.
                end if
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToFlowering()) &
                    .and.  (.not. AC71_struc(n)%ac71(t)%flowr_reached)) then ! from sowi
                    AC71_struc(n)%ac71(t)%Crop%DaysToFlowering = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%flowr_reached = .true.
                end if
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToSenescence()) &
                    .and. (.not. AC71_struc(n)%ac71(t)%sene_reached)) then ! from sowin
                    AC71_struc(n)%ac71(t)%Crop%DaysToSenescence = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%sene_reached = .true.
                end if
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToHarvest()) &
                    .and. (.not. AC71_struc(n)%ac71(t)%harv_reached)) then ! from sowing t
                    AC71_struc(n)%ac71(t)%Crop%DaysToHarvest = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%harv_reached = .true.
                end if
                end if
            end if ! Initialize crop stages done

            ! Run AC
            if (AC71_struc(n)%ac71(t)%Crop%ModeCycle.ne.-9) then ! do not run GDD if not enough
                tmp_wpi = REAL(AC71_struc(n)%ac71(t)%WPi,8)
                write(LIS_logunit,'(2f10.6)') lat, lon
                call AdvanceOneTimeStep(tmp_wpi, AC71_struc(n)%ac71(t)%HarvestNow)
                AC71_struc(n)%ac71(t)%WPi = tmp_wpi

                ! Close irri file if opened
                if(irr_record_flag.eq.1)then
                    call fIrri_close()
                endif

                ! Get all the ac71 variables and store in AC71_struc
                do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                        AC71_struc(n)%ac71(t)%smc(l) = GetCompartment_theta(l)
                enddo
                
                AC71_struc(n)%AC71(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
                AC71_struc(n)%AC71(t)%RootZoneWC_FC = GetRootZoneWC_FC()
                AC71_struc(n)%AC71(t)%RootZoneWC_WP = GetRootZoneWC_WP()
                AC71_struc(n)%AC71(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
                AC71_struc(n)%AC71(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
                AC71_struc(n)%AC71(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
                AC71_struc(n)%AC71(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
                AC71_struc(n)%AC71(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
                AC71_struc(n)%AC71(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
                AC71_struc(n)%AC71(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
                AC71_struc(n)%AC71(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
                AC71_struc(n)%AC71(t)%Compartment = GetCompartment()
                AC71_struc(n)%AC71(t)%TotalSaltContent = GetTotalSaltContent()
                AC71_struc(n)%AC71(t)%TotalWaterContent = GetTotalWaterContent()
                AC71_struc(n)%AC71(t)%effectiverain = Geteffectiverain()
                AC71_struc(n)%AC71(t)%SumWaBal = GetSumWaBal()
                AC71_struc(n)%AC71(t)%RootZoneSalt = GetRootZoneSalt()
                AC71_struc(n)%AC71(t)%Simulation = GetSimulation()
                AC71_struc(n)%AC71(t)%IrriInterval = GetIrriInterval()
                AC71_struc(n)%AC71(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
                AC71_struc(n)%AC71(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
                AC71_struc(n)%AC71(t)%Irrigation = GetIrrigation()
                AC71_struc(n)%AC71(t)%IrriBeforeSeason = GetIrriBeforeSeason()
                AC71_struc(n)%AC71(t)%IrriAfterSeason = GetIrriAfterSeason()
                AC71_struc(n)%AC71(t)%SoilLayer = GetSoilLayer()
                AC71_struc(n)%AC71(t)%daynri = GetDayNri()
                do l=1, AC71_struc(n)%AC71(t)%NrCompartments
                    AC71_struc(n)%AC71(t)%smc(l) = GetCompartment_theta(l)
                enddo
                !write(*,'(e23.15e3)') AC71_struc(n)%AC71(t)%AC71smc(1)
                AC71_struc(n)%AC71(t)%IrriECw = GetIrriECw()
                AC71_struc(n)%AC71(t)%Management = GetManagement()
                AC71_struc(n)%AC71(t)%PerennialPeriod = GetPerennialPeriod()
                AC71_struc(n)%AC71(t)%simulparam = GetSimulParam()
                AC71_struc(n)%AC71(t)%Cuttings = GetManagement_Cuttings()
                AC71_struc(n)%AC71(t)%onset = GetOnset()
                AC71_struc(n)%AC71(t)%endseason = GetEndSeason()
                AC71_struc(n)%AC71(t)%crop = GetCrop()
                AC71_struc(n)%AC71(t)%Soil = GetSoil()
                AC71_struc(n)%AC71(t)%TemperatureRecord = GetTemperatureRecord()
                AC71_struc(n)%AC71(t)%ClimRecord = GetClimRecord()
                AC71_struc(n)%AC71(t)%RainRecord = GetRainRecord()
                AC71_struc(n)%AC71(t)%EToRecord = GetEToRecord()

                AC71_struc(n)%AC71(t)%GenerateTimeMode = GetGenerateTimeMode()
                AC71_struc(n)%AC71(t)%GenerateDepthMode = GetGenerateDepthMode()
                AC71_struc(n)%AC71(t)%IrriMode = GetIrriMode()
                AC71_struc(n)%AC71(t)%IrriMethod = GetIrriMethod()
                AC71_struc(n)%AC71(t)%DaySubmerged = GetDaySubmerged()
                AC71_struc(n)%AC71(t)%MaxPlotNew = GetMaxPlotNew()
                AC71_struc(n)%AC71(t)%NrCompartments = GetNrCompartments()
                AC71_struc(n)%AC71(t)%IrriFirstDayNr = GetIrriFirstDayNr()
                AC71_struc(n)%AC71(t)%ZiAqua = GetZiAqua()
                AC71_struc(n)%AC71(t)%IniPercTAW = GetIniPercTAW()
                AC71_struc(n)%AC71(t)%MaxPlotTr = GetMaxPlotTr()
                AC71_struc(n)%AC71(t)%OutputAggregate = GetOutputAggregate()

                AC71_struc(n)%AC71(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
                AC71_struc(n)%AC71(t)%PreDay = GetPreDay()
                AC71_struc(n)%AC71(t)%OutDaily = GetOutDaily()
                AC71_struc(n)%AC71(t)%Out1Wabal = GetOut1Wabal()
                AC71_struc(n)%AC71(t)%Out2Crop = GetOut2Crop()
                AC71_struc(n)%AC71(t)%Out3Prof = GetOut3Prof()
                AC71_struc(n)%AC71(t)%Out4Salt = GetOut4Salt()
                AC71_struc(n)%AC71(t)%Out5CompWC = GetOut5CompWC()
                AC71_struc(n)%AC71(t)%Out6CompEC = GetOut6CompEC()
                AC71_struc(n)%AC71(t)%Out7Clim = GetOut7Clim()
                AC71_struc(n)%AC71(t)%Part1Mult = GetPart1Mult()
                AC71_struc(n)%AC71(t)%Part2Eval = GetPart2Eval()

                !
                AC71_struc(n)%AC71(t)%CCiActual = GetCCiActual()
                AC71_struc(n)%AC71(t)%CCiprev = GetCCiprev()
                AC71_struc(n)%AC71(t)%CCiTopEarlySen = GetCCiTopEarlySen()
                AC71_struc(n)%AC71(t)%CRsalt = GetCRsalt () ! gram/m2
                AC71_struc(n)%AC71(t)%CRwater = GetCRwater() ! mm/day
                AC71_struc(n)%AC71(t)%ECdrain = GetECdrain() ! EC drain water dS/m
                AC71_struc(n)%AC71(t)%ECiAqua = GetECiAqua() ! EC of the groundwater table in dS/m
                AC71_struc(n)%AC71(t)%ECstorage = GetECstorage() !EC surface storage dS/m
                AC71_struc(n)%AC71(t)%Eact = GetEact() ! mm/day
                AC71_struc(n)%AC71(t)%Epot = GetEpot() ! mm/day
                AC71_struc(n)%AC71(t)%ETo = GetETo() ! mm/day
                AC71_struc(n)%AC71(t)%Drain = GetDrain()  ! mm/day
                AC71_struc(n)%AC71(t)%Infiltrated = GetInfiltrated() ! mm/day
                AC71_struc(n)%AC71(t)%prcp = GetRain()  ! mm/day
                AC71_struc(n)%AC71(t)%RootingDepth = GetRootingDepth()
                AC71_struc(n)%AC71(t)%Runoff = GetRunoff()  ! mm/day
                AC71_struc(n)%AC71(t)%SaltInfiltr = GetSaltInfiltr() ! salt infiltrated in soil profile Mg/ha
                AC71_struc(n)%AC71(t)%Surf0 = GetSurf0()  ! surface water [mm] begin day
                AC71_struc(n)%AC71(t)%SurfaceStorage = GetSurfaceStorage() !mm/day
                AC71_struc(n)%AC71(t)%Tact = GetTact() ! mm/day
                AC71_struc(n)%AC71(t)%Tpot = GetTpot() ! mm/day
                AC71_struc(n)%AC71(t)%TactWeedInfested = GetTactWeedInfested() !mm/day
                AC71_struc(n)%AC71(t)%tmax = GetTmax() ! degC
                AC71_struc(n)%AC71(t)%tmin =GetTmin() ! degC


                AC71_struc(n)%AC71(t)%GwTable = GetGwTable()
                AC71_struc(n)%AC71(t)%PlotVarCrop = GetPlotVarCrop()
                AC71_struc(n)%AC71(t)%StressTot = GetStressTot()
                AC71_struc(n)%AC71(t)%CutInfoRecord1 = GetCutInfoRecord1()
                AC71_struc(n)%AC71(t)%CutInfoRecord2 = GetCutInfoRecord2()
                AC71_struc(n)%AC71(t)%Transfer = GetTransfer()
                AC71_struc(n)%AC71(t)%PreviousSum = GetPreviousSum()
                AC71_struc(n)%AC71(t)%Tadj = GetTadj()
                AC71_struc(n)%AC71(t)%GDDTadj = GetGDDTadj()
                AC71_struc(n)%AC71(t)%DayLastCut = GetDayLastCut()
                AC71_struc(n)%AC71(t)%NrCut = GetNrCut()
                AC71_struc(n)%AC71(t)%SumInterval = GetSumInterval()
                AC71_struc(n)%AC71(t)%PreviousStressLevel = GetPreviousStressLevel()
                AC71_struc(n)%AC71(t)%StressSFadjNEW = GetStressSFadjNEW()
                AC71_struc(n)%AC71(t)%Bin = GetBin()
                AC71_struc(n)%AC71(t)%Bout = GetBout()
                AC71_struc(n)%AC71(t)%GDDayi = GetGDDayi()
                AC71_struc(n)%AC71(t)%CO2i = GetCO2i()
                AC71_struc(n)%AC71(t)%FracBiomassPotSF = GetFracBiomassPotSF()
                AC71_struc(n)%AC71(t)%SumETo = GetSumETo()
                AC71_struc(n)%AC71(t)%SumGDD = GetSumGDD()
                AC71_struc(n)%AC71(t)%Ziprev = GetZiprev()
                AC71_struc(n)%AC71(t)%SumGDDPrev = GetSumGDDPrev()
                AC71_struc(n)%AC71(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
                AC71_struc(n)%AC71(t)%Coeffb0 = GetCoeffb0()
                AC71_struc(n)%AC71(t)%Coeffb1 = GetCoeffb1()
                AC71_struc(n)%AC71(t)%Coeffb2 = GetCoeffb2()
                AC71_struc(n)%AC71(t)%Coeffb0Salt = GetCoeffb0Salt()
                AC71_struc(n)%AC71(t)%Coeffb1Salt = GetCoeffb1Salt()
                AC71_struc(n)%AC71(t)%Coeffb2Salt = GetCoeffb2Salt()
                AC71_struc(n)%AC71(t)%StressLeaf = GetStressLeaf()
                AC71_struc(n)%AC71(t)%StressSenescence = GetStressSenescence()
                AC71_struc(n)%AC71(t)%DayFraction = GetDayFraction()
                AC71_struc(n)%AC71(t)%GDDayFraction = GetGDDayFraction()
                AC71_struc(n)%AC71(t)%CGCref = GetCGCref()
                AC71_struc(n)%AC71(t)%GDDCGCref = GetGDDCGCref()
                AC71_struc(n)%AC71(t)%TimeSenescence = GetTimeSenescence()
                AC71_struc(n)%AC71(t)%SumKcTop = GetSumKcTop()
                AC71_struc(n)%AC71(t)%SumKcTopStress = GetSumKcTopStress()
                AC71_struc(n)%AC71(t)%SumKci = GetSumKci()
                AC71_struc(n)%AC71(t)%CCoTotal = GetCCoTotal()
                AC71_struc(n)%AC71(t)%CCxTotal = GetCCxTotal()
                AC71_struc(n)%AC71(t)%CDCTotal = GetCDCTotal()
                AC71_struc(n)%AC71(t)%GDDCDCTotal = GetGDDCDCTotal()
                AC71_struc(n)%AC71(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
                AC71_struc(n)%AC71(t)%WeedRCi = GetWeedRCi()
                AC71_struc(n)%AC71(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
                AC71_struc(n)%AC71(t)%fWeedNoS = GetfWeedNoS()
                AC71_struc(n)%AC71(t)%Zeval = GetZeval()
                AC71_struc(n)%AC71(t)%BprevSum = GetBprevSum()
                AC71_struc(n)%AC71(t)%YprevSum = GetYprevSum()
                AC71_struc(n)%AC71(t)%SumGDDcuts = GetSumGDDcuts()
                AC71_struc(n)%AC71(t)%HItimesBEF = GetHItimesBEF()
                AC71_struc(n)%AC71(t)%ScorAT1 = GetScorAT1()
                AC71_struc(n)%AC71(t)%ScorAT2 = GetScorAT2()
                AC71_struc(n)%AC71(t)%HItimesAT1 = GetHItimesAT1()
                AC71_struc(n)%AC71(t)%HItimesAT2 = GetHItimesAT2()
                AC71_struc(n)%AC71(t)%HItimesAT = GetHItimesAT()
                AC71_struc(n)%AC71(t)%alfaHI = GetalfaHI()
                AC71_struc(n)%AC71(t)%alfaHIAdj = GetalfaHIAdj()
                AC71_struc(n)%AC71(t)%NextSimFromDayNr = GetNextSimFromDayNr ()
                AC71_struc(n)%AC71(t)%DayNr1Eval = GetDayNr1Eval()
                AC71_struc(n)%AC71(t)%DayNrEval = GetDayNrEval()
                AC71_struc(n)%AC71(t)%LineNrEval = GetLineNrEval()
                AC71_struc(n)%AC71(t)%PreviousSumETo = GetPreviousSumETo()
                AC71_struc(n)%AC71(t)%PreviousSumGDD = GetPreviousSumGDD()
                AC71_struc(n)%AC71(t)%PreviousBmob = GetPreviousBmob()
                AC71_struc(n)%AC71(t)%PreviousBsto = GetPreviousBsto()
                AC71_struc(n)%AC71(t)%StageCode = GetStageCode()
                AC71_struc(n)%AC71(t)%PreviousDayNr = GetPreviousDayNr()
                AC71_struc(n)%AC71(t)%NoYear = GetNoYear()
                AC71_struc(n)%AC71(t)%WaterTableInProfile = GetWaterTableInProfile()
                AC71_struc(n)%AC71(t)%StartMode = GetStartMode()
                AC71_struc(n)%AC71(t)%NoMoreCrop = GetNoMoreCrop()

                !Test
                AC71_struc(n)%ac71(t)%CropFileSet = GetCropFileSet()

                ! Check for end of simulation period 
                ! (DayNri - 1 because DayNri is already for next day)
                if ((GetDayNri()-1) .eq. GetSimulation_ToDayNr()) then
                    AC71_struc(n)%ac71(t)%InitializeRun = 1
                    AC71_struc(n)%ac71(t)%irun = AC71_struc(n)%ac71(t)%irun + 1
                end if


                ! Diagnostic output variables
                ![ 1] output variable: smc (unit=m^3 m-3 ). ***  volumetric soil moisture
                do i=1, AC71_struc(n)%ac71(t)%NrCompartments
                    call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = AC71_struc(n)%ac71(t)%smc(i), &
                                                        vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
                end do
                ![ 2] output variable: biomass (unit=t/ha).  *** cummulative biomass
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_biomass, value = real(AC71_struc(n)%ac71(t)%SumWaBal%Biomass,kind=sp), &
                                                    vlevel=1, unit="t ha-1", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 3] output variable: biomass (unit=mm).  *** actual rootzone water content
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RootZoneWC_Actual, value = real(AC71_struc(n)%ac71(t)%RootZoneWC_Actual,kind=sp), &
                                                    vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 4] output variable: RootZoneWC_WP (unit=mm).  *** rootzone water content at wilting point
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RootZoneWC_WP, value = real(AC71_struc(n)%ac71(t)%RootZoneWC_WP,kind=sp), &
                                                    vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 5] output variable: RootZoneWC_FC (unit=mm).  *** rootzone water content at field capacity
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RootZoneWC_FC, value = real(AC71_struc(n)%ac71(t)%RootZoneWC_FC,kind=sp), &
                                                    vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 6] output variable: Tact (unit=mm).  *** actual transpiration
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_Tact, value = real(AC71_struc(n)%ac71(t)%Tact,kind=sp), &
                                                    vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 7] output variable: Eact (unit=mm).  *** actual evaporation
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_Eact, value = real(AC71_struc(n)%ac71(t)%Eact,kind=sp), &
                                                    vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 8] output variable: AC71ETo (unit=mm).  *** reference evapotranspiration (Penman-Monteith)
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71ETo, value = real(AC71_struc(n)%ac71(t)%eto,kind=sp), &
                                                    vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 9] output variable: RootingDepth (unit=m).  *** rooting depth
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RootingDepth, value = real(AC71_struc(n)%ac71(t)%RootingDepth,kind=sp), &
                                                    vlevel=1, unit="m", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 10] output variable: CCiActual (unit=-).  *** canopy cover
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CCiActual, value = real(AC71_struc(n)%ac71(t)%CCiActual,kind=sp), &
                                                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 11] output variable: AC71Tmin (unit=deg C).  *** daily minimum temperature
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_Tmin, value = real(tmp_tmin,kind=sp), &
                                        vlevel=1, unit="degC", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 12] output variable: AC71Tmax (unit=deg C).  *** daily maximum temperature
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_Tmax, value = real(tmp_tmax,kind=sp), &
                                        vlevel=1, unit="degC", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 13] output variable: AC71Rain (unit=mm).  *** precipitation rate
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Rain, value = real(tmp_precip,kind=sp), &
                                        vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 13] output variable: GD (unit=degC).  *** growing degrees
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_GD, value = real(GetGDDayi(),kind=sp), &
                                        vlevel=1, unit="degC", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 14] output variable: yield (unit=t ha-1).  *** yield
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_Yield, value = real(AC71_struc(n)%ac71(t)%SumWaBal%YieldPart,kind=sp), &
                                        vlevel=1, unit="t ha-1", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 15] output variable: irrigation (unit=mm).  *** irrigation
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Irrigation, value = real(AC71_struc(n)%ac71(t)%Irrigation,kind=sp), &
                                        vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 16] output variable: StExp (unit=%).  *** expansion stress
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_StExp, value = real(GetStressLeaf(),kind=sp), &
                                        vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)
                ![ 17] output variable: StSen (unit=%).  *** senescence stress
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_StSen, value = real(GetStressSenescence(),kind=sp), &
                                        vlevel=1, unit="%", direction="-", surface_type = LIS_rc%lsm_index)
            endif

            !  Reset forcings
            AC71_struc(n)%ac71(t)%tair = 0.0
            AC71_struc(n)%ac71(t)%tmax = 0.0
            AC71_struc(n)%ac71(t)%tmin = 0.0
            AC71_struc(n)%ac71(t)%tdew = 0.0
            AC71_struc(n)%ac71(t)%wndspd = 0.0            
            AC71_struc(n)%ac71(t)%psurf = 0.0
            AC71_struc(n)%ac71(t)%prcp = 0.0
            AC71_struc(n)%ac71(t)%eto = 0.0
            AC71_struc(n)%ac71(t)%swdown = 0.0
        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        AC71_struc(n)%forc_count = 0
    endif ! end of alarmCheck loop 

end subroutine Ac71_main
