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
    use Ac71_lsmMod
    use ac71_prep_f, only : ac71_ETo_calc
   !use other modules
    use ESMF
    use LIS_routingMod, only : LIS_runoff_state

!   ! AC module imports
    use ac_global, only: &
                    DegreesDay,&
                    GetCCiActual,&
                    GetCCiTopEarlySen,&
                    GetCCiprev,&
                    GetCRsalt,&  
                    GetCRwater,& 
                    GetClimRecord,&
                    GetCompartment,&
                    GetCompartment_theta,&
                    GetCrop,&
                    GetCrop_Day1,&
                    GetDaySubmerged,&
                    GetDrain,&  
                    GetECdrain,& 
                    GetECiAqua,& 
                    GetECstorage,& 
                    GetETo,&
                    GetEToRecord,&
                    GetEact,& 
                    GetEndSeason,&
                    GetEpot,& 
                    GetEvapoEntireSoilSurface,&
                    GetGenerateDepthMode,&
                    GetGenerateTimeMode,&
                    GetInfiltrated,&
                    GetIniPercTAW,&
                    GetIrriAfterSeason,&
                    GetIrriBeforeSeason,&
                    GetIrriECw,&
                    GetIrriFirstDayNr,&
                    GetIrriMethod,&
                    GetIrriMode,&
                    GetIrrigation,&
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
                    GetOutDaily,&
                    GetOutputAggregate,&
                    GetPart1Mult,&
                    GetPart2Eval,&
                    GetPerennialPeriod,&
                    GetPreDay,&
                    GetRain,& 
                    GetRainRecord,&
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
                    GetRootingDepth,&
                    GetRunoff,& 
                    GetSaltInfiltr,&
                    GetSimulParam,&
                    GetSimulParam_GDDMethod,&
                    GetSimulation,&
                    GetSimulation_SumGDD,&
                    GetSimulation_SumGDDfromDay1,&
                    GetSoil,&
                    GetSoilLayer,&
                    GetSumWaBal,&
                    GetSurf0,& 
                    GetSurfaceStorage,&
                    GetTact,&
                    GetTactWeedInfested,&
                    GetTemperatureRecord,&
                    GetTmax,& 
                    GetTmin,& 
                    GetTotalSaltContent,&
                    GetTotalWaterContent,&
                    GetTpot,&
                    GetZiAqua,&
                    Geteffectiverain,&
                    SetCCiActual,&
                    SetCCiTopEarlySen,&
                    SetCCiprev,&
                    SetCRsalt,&  
                    SetCRwater,& 
                    SetClimFile, &
                    SetClimRecord,&
                    SetClimRecord_DataType, &
                    SetClimRecord_NrObs, &
                    SetClimRecord_fromd, &
                    SetClimRecord_fromdaynr, &
                    SetClimRecord_fromm, &
                    SetClimRecord_fromstring, &
                    SetClimRecord_fromy, &
                    SetClimRecord_tod, &
                    SetClimRecord_todaynr, &
                    SetClimRecord_tom, &
                    SetClimRecord_tostring, &
                    SetClimRecord_toy,&
                    SetCompartment,&
                    SetCompartment_theta,&
                    SetCrop,&
                    SetDaySubmerged,&
                    SetDrain,&  
                    SetECdrain,& 
                    SetECiAqua,& 
                    SetECstorage,& 
                    SetETo,&
                    SetEact,& 
                    SetEndSeason,&
                    SetEpot,& 
                    SetEvapoEntireSoilSurface,&
                    SetGenerateDepthMode,&
                    SetGenerateTimeMode,&
                    SetInfiltrated,&
                    SetIniPercTAW,&
                    SetIrriAfterSeason,&
                    SetIrriBeforeSeason,&
                    SetIrriECw,&
                    SetIrriFirstDayNr,&
                    SetIrriMode,&
                    SetIrrigation,&
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
                    SetOutDaily,&
                    SetOutputAggregate,&
                    SetPart1Mult,&
                    SetPart2Eval,&
                    SetPerennialPeriod,&
                    SetPreDay,&
                    SetRain,& 
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
                    SetRootingDepth,&
                    SetRunoff,& 
                    SetSaltInfiltr,&
                    SetSimulParam,&
                    SetSimulation,&
                    SetSimulation_SumGDD, &
                    SetSimulation_SumGDDfromDay1,&
                    SetSoil,&
                    SetSoilLayer,&
                    SetSumWaBal,&
                    SetSurf0,& 
                    SetSurfaceStorage,&
                    SetTact,&
                    SetTactWeedInfested,&
                    SetTemperatureRecord,&
                    SetTmax,& 
                    SetTmin,&
                    SetTotalSaltContent,&
                    SetTotalWaterContent,&
                    SetTpot,&
                    SetZiAqua,&
                    Seteffectiverain
    use ac_kinds, only: intEnum, &
                        int32, &
                        int8, &
                        dp,&
                        sp
    use ac_project_input, only: ProjectInput, set_project_input
    use ac_run, only: &
                    AdvanceOneTimeStep, &
                    FinalizeRun1, &
                    FinalizeRun2, &
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
                    GetGDDCDCTotal,&
                    GetGDDCGCref ,&
                    GetGDDTadj,&
                    GetGDDayFraction,&
                    GetGDDayi,&
                    GetGwTable,&
                    GetHItimesAT,&
                    GetHItimesAT1,&
                    GetHItimesAT2,&
                    GetHItimesBEF,&
                    GetIrriInfoRecord1,&
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
                    GetStressSFadjNEW,&
                    GetStressSenescence ,&
                    GetStressTot,&
                    GetSumETo,&
                    GetSumGDD,&
                    GetSumGDDPrev,&
                    GetSumGDDcuts,&
                    GetSumInterval,&
                    GetSumKcTop,&
                    GetSumKcTopStress,&
                    GetSumKci,&
                    GetTadj,&
                    GetTheProjectFile,&
                    GetTimeSenescence ,&
                    GetTransfer,&
                    GetWaterTableInProfile,&
                    GetWeedRCi,&
                    GetYprevSum,&
                    GetZeval,&
                    GetZiprev,&
                    GetalfaHI,&
                    GetalfaHIAdj,&
                    GetfWeedNoS,&
                    InitializeRunPart1, &
                    InitializeSimulation, &
                    InitializeSimulationRunPart2, &
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
                    SetGDDCDCTotal,&
                    SetGDDCGCref ,&
                    SetGDDTadj,&
                    SetGDDayFraction,&
                    SetGDDayi,&
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
                    SetStressSFadjNEW,&
                    SetStressSenescence ,&
                    SetStressTot,&
                    SetSumETo,&
                    SetSumGDD,&
                    SetSumGDDPrev,&
                    SetSumGDDcuts,&
                    SetSumInterval,&
                    SetSumKcTop,&
                    SetSumKcTopStress,&
                    SetSumKci,&
                    SetTadj,&
                    SetTimeSenescence ,&
                    SetTransfer,&
                    SetWaterTableInProfile,&
                    SetWeedRCi,&
                    SetYprevSum,&
                    SetZeval,&
                    SetZiprev,&
                    SetalfaHI,&
                    SetalfaHIAdj,&
                    SetfWeedNoS                     
    use ac_startunit, only: GetSimulation_NrRuns

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
    logical              :: alarmCheck

    integer              :: status
    integer              :: c,r,l
    integer              :: ios, nid,rivid,fldid

    integer              :: tid

    !!! MB_AC71
    integer              :: daynr, todaynr, iproject, nprojects
    logical              :: ListProjectFileExist
    character(len=:), allocatable :: ListProjectsFile, TheProjectFile

    !LB AC71

    real                 :: tmp_pres, tmp_precip, tmp_tmax, tmp_tmin   ! Weather Forcing
    real                 :: tmp_tdew, tmp_swrad, tmp_wind, tmp_eto     ! Weather Forcing
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
            tmp_tmax      = AC71_struc(n)%ac71(t)%tmax - 273.15 !Convert from K to C

            ! TMIN: minimum daily air temperature (degC)
            tmp_tmin      = AC71_struc(n)%ac71(t)%tmin - 273.15 !Convert from K to C 

            ! TDEW: average daily dewpoint temperature (degC)
            tmp_tdew      = (AC71_struc(n)%ac71(t)%tdew / AC71_struc(n)%forc_count) - 273.15 !Convert from K to C

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

            ! later --> add if statements to correct for 10m data (ERA5 10m T and 10m wind)
            ! Add in config file reference height for wind
             
            call ac71_ETo_calc(tmp_pres, tmp_tmax, tmp_tmin, tmp_tdew, &
                               tmp_wind, tmp_swrad, &
                               tmp_elev, lat, tmp_eto)
            AC71_struc(n)%ac71(t)%eto = tmp_eto

            ! setting all global variables
            do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                    call SetCompartment_theta(l,REAL(AC71_struc(n)%ac71(t)%smc(l),8))
            enddo
            call SetBin(AC71_struc(n)%ac71(t)%Bin)
            call SetBout(AC71_struc(n)%ac71(t)%Bout)
            call SetBprevSum(AC71_struc(n)%ac71(t)%BprevSum)
            call SetCCiActual(AC71_struc(n)%ac71(t)%CCiActual)
            call SetCCiActualWeedInfested(AC71_struc(n)%ac71(t)%CCiActualWeedInfested)
            call SetCCiTopEarlySen(AC71_struc(n)%ac71(t)%CCiTopEarlySen)
            call SetCCiprev(AC71_struc(n)%ac71(t)%CCiprev)
            call SetCCoTotal(AC71_struc(n)%ac71(t)%CCoTotal)
            call SetCCxCropWeedsNoSFstress(AC71_struc(n)%ac71(t)%CCxCropWeedsNoSFstress)
            call SetCCxTotal(AC71_struc(n)%ac71(t)%CCxTotal)
            call SetCCxWitheredTpotNoS(AC71_struc(n)%ac71(t)%CCxWitheredTpotNoS)
            call SetCDCTotal(AC71_struc(n)%ac71(t)%CDCTotal)
            call SetCGCref(AC71_struc(n)%ac71(t)%CGCref)
            call SetCO2i(AC71_struc(n)%ac71(t)%CO2i)
            call SetCRsalt(AC71_struc(n)%ac71(t)%CRsalt)
            call SetCRwater(AC71_struc(n)%ac71(t)%CRwater)
            call SetCoeffb0(AC71_struc(n)%ac71(t)%Coeffb0)
            call SetCoeffb0Salt(AC71_struc(n)%ac71(t)%Coeffb0Salt)
            call SetCoeffb1(AC71_struc(n)%ac71(t)%Coeffb1)
            call SetCoeffb1Salt(AC71_struc(n)%ac71(t)%Coeffb1Salt)
            call SetCoeffb2(AC71_struc(n)%ac71(t)%Coeffb2)
            call SetCoeffb2Salt(AC71_struc(n)%ac71(t)%Coeffb2Salt)
            call SetCompartment(AC71_struc(n)%ac71(t)%Compartment)
            call SetCrop(AC71_struc(n)%ac71(t)%crop) 
            call SetCutInfoRecord1(AC71_struc(n)%ac71(t)%CutInfoRecord1)
            call SetCutInfoRecord2(AC71_struc(n)%ac71(t)%CutInfoRecord2)
            call SetDayFraction(AC71_struc(n)%ac71(t)%DayFraction)
            call SetDayLastCut(AC71_struc(n)%ac71(t)%DayLastCut)
            call SetDayNr1Eval(AC71_struc(n)%ac71(t)%DayNr1Eval)
            call SetDayNrEval(AC71_struc(n)%ac71(t)%DayNrEval)
            call SetDayNri(AC71_struc(n)%ac71(t)%daynri)
            call SetDaySubmerged(AC71_struc(n)%ac71(t)%DaySubmerged) 
            call SetDrain(AC71_struc(n)%ac71(t)%Drain)
            call SetECdrain(AC71_struc(n)%ac71(t)%ECdrain)
            call SetECstorage(AC71_struc(n)%ac71(t)%ECstorage)
            call SetEact(AC71_struc(n)%ac71(t)%Eact)
            call SetEciAqua(AC71_struc(n)%ac71(t)%ECiAqua)
            call SetEndSeason(AC71_struc(n)%ac71(t)%endseason) 
            call SetEpot(AC71_struc(n)%ac71(t)%Epot)
            call SetEvapoEntireSoilSurface(AC71_struc(n)%ac71(t)%EvapoEntireSoilSurface) 
            call SetFracBiomassPotSF(AC71_struc(n)%ac71(t)%FracBiomassPotSF)
            call SetGDDCDCTotal(AC71_struc(n)%ac71(t)%GDDCDCTotal)
            call SetGDDCGCref(AC71_struc(n)%ac71(t)%GDDCGCref)
            call SetGDDTadj(AC71_struc(n)%ac71(t)%GDDTadj)
            call SetGDDayFraction(AC71_struc(n)%ac71(t)%GDDayFraction)
            call SetGDDayi(AC71_struc(n)%ac71(t)%GDDayi)
            call SetGenerateDepthMode(AC71_struc(n)%ac71(t)%GenerateDepthMode) 
            call SetGenerateTimeMode(AC71_struc(n)%ac71(t)%GenerateTimeMode) 
            call SetGwTable(AC71_struc(n)%ac71(t)%GwTable)
            call SetHItimesAT(AC71_struc(n)%ac71(t)%HItimesAT)
            call SetHItimesAT1(AC71_struc(n)%ac71(t)%HItimesAT1)
            call SetHItimesAT2(AC71_struc(n)%ac71(t)%HItimesAT2)
            call SetHItimesBEF(AC71_struc(n)%ac71(t)%HItimesBEF)
            call SetInfiltrated(AC71_struc(n)%ac71(t)%Infiltrated)
            call SetIniPercTAW(AC71_struc(n)%ac71(t)%IniPercTAW) 
            call SetIrriAfterSeason(AC71_struc(n)%ac71(t)%IrriAfterSeason)
            call SetIrriBeforeSeason(AC71_struc(n)%ac71(t)%IrriBeforeSeason)
            call SetIrriECw(AC71_struc(n)%ac71(t)%IrriECw) 
            call SetIrriFirstDayNr(AC71_struc(n)%ac71(t)%IrriFirstDayNr) 
            call SetIrriInfoRecord1(AC71_struc(n)%ac71(t)%IrriInfoRecord1)
            call SetIrriInfoRecord2(AC71_struc(n)%ac71(t)%IrriInfoRecord2)
            call SetIrriInterval(AC71_struc(n)%ac71(t)%IrriInterval)
            call SetIrriMode(AC71_struc(n)%ac71(t)%IrriMode) 
            call SetIrrigation(REAL(AC71_struc(n)%ac71(t)%Irrigation,8))
            call SetLineNrEval(int(AC71_struc(n)%ac71(t)%LineNrEval,kind=int32))
            call SetManagement(AC71_struc(n)%ac71(t)%Management) 
            call SetManagement_Cuttings(AC71_struc(n)%ac71(t)%Cuttings) 
            call SetMaxPlotNew(AC71_struc(n)%ac71(t)%MaxPlotNew) 
            call SetMaxPlotTr(AC71_struc(n)%ac71(t)%MaxPlotTr)
            call SetNextSimFromDayNr(AC71_struc(n)%ac71(t)%NextSimFromDayNr)
            call SetNoMoreCrop(AC71_struc(n)%ac71(t)%NoMoreCrop)
            call SetNoYear(AC71_struc(n)%ac71(t)%NoYear)
            call SetNrCompartments(AC71_struc(n)%ac71(t)%NrCompartments) 
            call SetNrCut(AC71_struc(n)%ac71(t)%NrCut)
            call SetOnset(AC71_struc(n)%ac71(t)%onset) 
            call SetOut1Wabal(AC71_struc(n)%ac71(t)%Out1Wabal) 
            call SetOut2Crop(AC71_struc(n)%ac71(t)%Out2Crop) 
            call SetOut3Prof(AC71_struc(n)%ac71(t)%Out3Prof) 
            call SetOut4Salt(AC71_struc(n)%ac71(t)%Out4Salt) 
            call SetOut5CompWC(AC71_struc(n)%ac71(t)%Out5CompWC) 
            call SetOut6CompEC(AC71_struc(n)%ac71(t)%Out6CompEC) 
            call SetOut7Clim(AC71_struc(n)%ac71(t)%Out7Clim) 
            call SetOutDaily(AC71_struc(n)%ac71(t)%OutDaily) 
            call SetOutputAggregate(AC71_struc(n)%ac71(t)%OutputAggregate) 
            call SetPart1Mult(AC71_struc(n)%ac71(t)%Part1Mult) 
            call SetPart2Eval(AC71_struc(n)%ac71(t)%Part2Eval) 
            call SetPerennialPeriod(AC71_struc(n)%ac71(t)%PerennialPeriod) 
            call SetPlotVarCrop(AC71_struc(n)%ac71(t)%PlotVarCrop)
            call SetPreDay(AC71_struc(n)%ac71(t)%PreDay) 
            call SetPreviousBmob(AC71_struc(n)%ac71(t)%PreviousBmob)
            call SetPreviousBsto(AC71_struc(n)%ac71(t)%PreviousBsto)
            call SetPreviousDayNr(AC71_struc(n)%ac71(t)%PreviousDayNr)
            call SetPreviousStressLevel(int(AC71_struc(n)%ac71(t)%PreviousStressLevel,kind=int32))
            call SetPreviousSum(AC71_struc(n)%ac71(t)%PreviousSum)
            call SetPreviousSumETo(AC71_struc(n)%ac71(t)%PreviousSumETo)
            call SetPreviousSumGDD(AC71_struc(n)%ac71(t)%PreviousSumGDD)
            call SetRootZoneSalt(AC71_struc(n)%ac71(t)%RootZoneSalt)
            call SetRootZoneWC_Actual(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_Actual,8))
            call SetRootZoneWC_FC(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_FC,8))
            call SetRootZoneWC_Leaf(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_Leaf,8))
            call SetRootZoneWC_SAT(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_SAT,8))
            call SetRootZoneWC_Sen(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_Thresh(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_WP(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_WP,8))
            call SetRootZoneWC_ZtopAct(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopFC(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopThresh(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopWP(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootingDepth(AC71_struc(n)%ac71(t)%RootingDepth)
            call SetRunoff(AC71_struc(n)%ac71(t)%Runoff)
            call SetSaltInfiltr(AC71_struc(n)%ac71(t)%SaltInfiltr)
            call SetScorAT1(AC71_struc(n)%ac71(t)%ScorAT1)
            call SetScorAT2(AC71_struc(n)%ac71(t)%ScorAT2)
            call SetSimulParam(AC71_struc(n)%ac71(t)%simulparam) 
            call SetSimulation(AC71_struc(n)%ac71(t)%Simulation)
            call SetSoil(AC71_struc(n)%ac71(t)%Soil) 
            call SetSoilLayer(AC71_struc(n)%ac71(t)%soillayer)
            call SetStageCode(AC71_struc(n)%ac71(t)%StageCode)
            call SetStartMode(AC71_struc(n)%ac71(t)%StartMode)
            call SetStressLeaf(AC71_struc(n)%ac71(t)%StressLeaf)
            call SetStressSFadjNEW(int(AC71_struc(n)%ac71(t)%StressSFadjNEW,kind=int32))
            call SetStressSenescence(AC71_struc(n)%ac71(t)%StressSenescence)
            call SetStressTot(AC71_struc(n)%ac71(t)%StressTot)
            call SetSumETo(AC71_struc(n)%ac71(t)%SumETo)
            call SetSumGDD(AC71_struc(n)%ac71(t)%SumGDD)
            call SetSumGDDPrev(AC71_struc(n)%ac71(t)%SumGDDPrev)
            call SetSumGDDcuts(AC71_struc(n)%ac71(t)%SumGDDcuts)
            call SetSumInterval(AC71_struc(n)%ac71(t)%SumInterval)
            call SetSumKcTop(AC71_struc(n)%ac71(t)%SumKcTop)
            call SetSumKcTopStress(AC71_struc(n)%ac71(t)%SumKcTopStress)
            call SetSumKci(AC71_struc(n)%ac71(t)%SumKci)
            call SetSumWaBal(AC71_struc(n)%ac71(t)%SumWaBal)
            call SetSurf0(AC71_struc(n)%ac71(t)%Surf0)
            call SetSurfaceStorage(AC71_struc(n)%ac71(t)%SurfaceStorage)
            call SetTact(AC71_struc(n)%ac71(t)%Tact)
            call SetTactWeedInfested(AC71_struc(n)%ac71(t)%TactWeedInfested)
            call SetTadj(AC71_struc(n)%ac71(t)%Tadj)
            call SetTimeSenescence(AC71_struc(n)%ac71(t)%TimeSenescence)
            call SetTotalSaltContent(AC71_struc(n)%ac71(t)%TotalSaltContent)
            call SetTotalWaterContent(AC71_struc(n)%ac71(t)%TotalWaterContent)
            call SetTpot(AC71_struc(n)%ac71(t)%Tpot)
            call SetTransfer(AC71_struc(n)%ac71(t)%Transfer)
            call SetWaterTableInProfile(AC71_struc(n)%ac71(t)%WaterTableInProfile)
            call SetWeedRCi(AC71_struc(n)%ac71(t)%WeedRCi)
            call SetYprevSum(AC71_struc(n)%ac71(t)%YprevSum)
            call SetZeval(AC71_struc(n)%ac71(t)%Zeval)
            call SetZiAqua(AC71_struc(n)%ac71(t)%ZiAqua) 
            call SetZiprev(AC71_struc(n)%ac71(t)%Ziprev)
            call SetalfaHI(AC71_struc(n)%ac71(t)%alfaHI)
            call SetalfaHIAdj(AC71_struc(n)%ac71(t)%alfaHIAdj)
            call Seteffectiverain(AC71_struc(n)%ac71(t)%effectiverain)
            call SetfWeedNoS(AC71_struc(n)%ac71(t)%fWeedNoS)

            !!! initialize run (year)

            if (AC71_struc(n)%ac71(t)%InitializeRun .eq. 1) then
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
                call set_project_input(int(AC71_struc(n)%ac71(t)%irun, kind=int32), &
                                       'Crop_Filename', &
                                        trim(AC71_struc(n)%ac71(t)%cropt)//'.CRO')

                call InitializeRunPart1(AC71_struc(n)%ac71(t)%irun, AC71_struc(n)%ac71(t)%TheProjectType)
                call InitializeSimulationRunPart2()
                AC71_struc(n)%ac71(t)%HarvestNow = .false.
                AC71_struc(n)%ac71(t)%InitializeRun = 0
            end if

            ! Set climate variables then advanceonetimestep
            call SetRain(real(tmp_precip,kind=dp))
            call SetTmin(real(tmp_tmin,kind=dp))
            call SetTmax(real(tmp_tmax,kind=dp))
            call SetETo(real(tmp_eto,kind=dp))
            ! Sum of GDD at end of first day ! Wait for GDD implementation from Michel
            call SetGDDayi(DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), GetTmin(), &
                    GetTmax(), GetSimulParam_GDDMethod()))
            if (GetDayNri() >= GetCrop_Day1()) then
                if (GetDayNri() == GetCrop_Day1()) then
                    call SetSimulation_SumGDD(GetSimulation_SumGDD() + GetGDDayi())
                end if
                call SetSimulation_SumGDDfromDay1(GetSimulation_SumGDDfromDay1() + &
                    GetGDDayi())
            end if

            ! Run AC
            call AdvanceOneTimeStep(AC71_struc(n)%ac71(t)%WPi, AC71_struc(n)%ac71(t)%HarvestNow)

            ! Get all the ac71 variables and store in AC71_struc
            do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                    AC71_struc(n)%ac71(t)%smc(l) = GetCompartment_theta(l)
            enddo
            AC71_struc(n)%ac71(t)%Bin = GetBin()
            AC71_struc(n)%ac71(t)%Bout = GetBout()
            AC71_struc(n)%ac71(t)%BprevSum = GetBprevSum()
            AC71_struc(n)%ac71(t)%CCiActual = GetCCiActual()
            AC71_struc(n)%ac71(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
            AC71_struc(n)%ac71(t)%CCiTopEarlySen = GetCCiTopEarlySen()
            AC71_struc(n)%ac71(t)%CCiprev = GetCCiprev()
            AC71_struc(n)%ac71(t)%CCoTotal = GetCCoTotal()
            AC71_struc(n)%ac71(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
            AC71_struc(n)%ac71(t)%CCxTotal = GetCCxTotal()
            AC71_struc(n)%ac71(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
            AC71_struc(n)%ac71(t)%CDCTotal = GetCDCTotal()
            AC71_struc(n)%ac71(t)%CGCref = GetCGCref()
            AC71_struc(n)%ac71(t)%CO2i = GetCO2i()
            AC71_struc(n)%ac71(t)%CRsalt = GetCRsalt ()
            AC71_struc(n)%ac71(t)%CRwater = GetCRwater()
            AC71_struc(n)%ac71(t)%ClimRecord = GetClimRecord()
            AC71_struc(n)%ac71(t)%Coeffb0 = GetCoeffb0()
            AC71_struc(n)%ac71(t)%Coeffb0Salt = GetCoeffb0Salt()
            AC71_struc(n)%ac71(t)%Coeffb1 = GetCoeffb1()
            AC71_struc(n)%ac71(t)%Coeffb1Salt = GetCoeffb1Salt()
            AC71_struc(n)%ac71(t)%Coeffb2 = GetCoeffb2()
            AC71_struc(n)%ac71(t)%Coeffb2Salt = GetCoeffb2Salt()
            AC71_struc(n)%ac71(t)%Compartment = GetCompartment()
            AC71_struc(n)%ac71(t)%CutInfoRecord1 = GetCutInfoRecord1()
            AC71_struc(n)%ac71(t)%CutInfoRecord2 = GetCutInfoRecord2()
            AC71_struc(n)%ac71(t)%Cuttings = GetManagement_Cuttings()
            AC71_struc(n)%ac71(t)%DayFraction = GetDayFraction()
            AC71_struc(n)%ac71(t)%DayLastCut = GetDayLastCut()
            AC71_struc(n)%ac71(t)%DayNr1Eval = GetDayNr1Eval()
            AC71_struc(n)%ac71(t)%DayNrEval = GetDayNrEval()
            AC71_struc(n)%ac71(t)%DaySubmerged = GetDaySubmerged()
            AC71_struc(n)%ac71(t)%Drain = GetDrain()
            AC71_struc(n)%ac71(t)%ECdrain = GetECdrain()
            AC71_struc(n)%ac71(t)%ECiAqua = GetECiAqua()
            AC71_struc(n)%ac71(t)%ECstorage = GetECstorage()
            AC71_struc(n)%ac71(t)%EToRecord = GetEToRecord()
            AC71_struc(n)%ac71(t)%Eact = GetEact()
            AC71_struc(n)%ac71(t)%Epot = GetEpot()
            AC71_struc(n)%ac71(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
            AC71_struc(n)%ac71(t)%FracBiomassPotSF = GetFracBiomassPotSF()
            AC71_struc(n)%ac71(t)%GDDCDCTotal = GetGDDCDCTotal()
            AC71_struc(n)%ac71(t)%GDDCGCref = GetGDDCGCref ()
            AC71_struc(n)%ac71(t)%GDDTadj = GetGDDTadj()
            AC71_struc(n)%ac71(t)%GDDayFraction = GetGDDayFraction()
            AC71_struc(n)%ac71(t)%GDDayi = GetGDDayi()
            AC71_struc(n)%ac71(t)%GenerateDepthMode = GetGenerateDepthMode()
            AC71_struc(n)%ac71(t)%GenerateTimeMode = GetGenerateTimeMode()
            AC71_struc(n)%ac71(t)%GwTable = GetGwTable()
            AC71_struc(n)%ac71(t)%HItimesAT = GetHItimesAT()
            AC71_struc(n)%ac71(t)%HItimesAT1 = GetHItimesAT1()
            AC71_struc(n)%ac71(t)%HItimesAT2 = GetHItimesAT2()
            AC71_struc(n)%ac71(t)%HItimesBEF = GetHItimesBEF()
            AC71_struc(n)%ac71(t)%Infiltrated = GetInfiltrated()
            AC71_struc(n)%ac71(t)%IniPercTAW = GetIniPercTAW()
            AC71_struc(n)%ac71(t)%IrriAfterSeason = GetIrriAfterSeason()
            AC71_struc(n)%ac71(t)%IrriBeforeSeason = GetIrriBeforeSeason()
            AC71_struc(n)%ac71(t)%IrriECw = GetIrriECw()
            AC71_struc(n)%ac71(t)%IrriFirstDayNr = GetIrriFirstDayNr()
            AC71_struc(n)%ac71(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
            AC71_struc(n)%ac71(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
            AC71_struc(n)%ac71(t)%IrriInterval = GetIrriInterval()
            AC71_struc(n)%ac71(t)%IrriMethod = GetIrriMethod()
            AC71_struc(n)%ac71(t)%IrriMode = GetIrriMode()
            AC71_struc(n)%ac71(t)%Irrigation = GetIrrigation()
            AC71_struc(n)%ac71(t)%LineNrEval = GetLineNrEval()
            AC71_struc(n)%ac71(t)%Management = GetManagement()
            AC71_struc(n)%ac71(t)%MaxPlotNew = GetMaxPlotNew()
            AC71_struc(n)%ac71(t)%MaxPlotTr = GetMaxPlotTr()
            AC71_struc(n)%ac71(t)%NextSimFromDayNr = GetNextSimFromDayNr ()
            AC71_struc(n)%ac71(t)%NoMoreCrop = GetNoMoreCrop()
            AC71_struc(n)%ac71(t)%NoYear = GetNoYear()
            AC71_struc(n)%ac71(t)%NrCompartments = GetNrCompartments()
            AC71_struc(n)%ac71(t)%NrCut = GetNrCut()
            AC71_struc(n)%ac71(t)%Out1Wabal = GetOut1Wabal()
            AC71_struc(n)%ac71(t)%Out2Crop = GetOut2Crop()
            AC71_struc(n)%ac71(t)%Out3Prof = GetOut3Prof()
            AC71_struc(n)%ac71(t)%Out4Salt = GetOut4Salt()
            AC71_struc(n)%ac71(t)%Out5CompWC = GetOut5CompWC()
            AC71_struc(n)%ac71(t)%Out6CompEC = GetOut6CompEC()
            AC71_struc(n)%ac71(t)%Out7Clim = GetOut7Clim()
            AC71_struc(n)%ac71(t)%OutDaily = GetOutDaily()
            AC71_struc(n)%ac71(t)%OutputAggregate = GetOutputAggregate()
            AC71_struc(n)%ac71(t)%Part1Mult = GetPart1Mult()
            AC71_struc(n)%ac71(t)%Part2Eval = GetPart2Eval()
            AC71_struc(n)%ac71(t)%PerennialPeriod = GetPerennialPeriod()
            AC71_struc(n)%ac71(t)%PlotVarCrop = GetPlotVarCrop()
            AC71_struc(n)%ac71(t)%PreDay = GetPreDay()
            AC71_struc(n)%ac71(t)%PreviousBmob = GetPreviousBmob()
            AC71_struc(n)%ac71(t)%PreviousBsto = GetPreviousBsto()
            AC71_struc(n)%ac71(t)%PreviousDayNr = GetPreviousDayNr()
            AC71_struc(n)%ac71(t)%PreviousStressLevel = GetPreviousStressLevel()
            AC71_struc(n)%ac71(t)%PreviousSum = GetPreviousSum()
            AC71_struc(n)%ac71(t)%PreviousSumETo = GetPreviousSumETo()
            AC71_struc(n)%ac71(t)%PreviousSumGDD = GetPreviousSumGDD()
            AC71_struc(n)%ac71(t)%RainRecord = GetRainRecord()
            AC71_struc(n)%ac71(t)%RootZoneSalt = GetRootZoneSalt()
            AC71_struc(n)%ac71(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
            AC71_struc(n)%ac71(t)%RootZoneWC_FC = GetRootZoneWC_FC()
            AC71_struc(n)%ac71(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
            AC71_struc(n)%ac71(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
            AC71_struc(n)%ac71(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
            AC71_struc(n)%ac71(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
            AC71_struc(n)%ac71(t)%RootZoneWC_WP = GetRootZoneWC_WP()
            AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
            AC71_struc(n)%ac71(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
            AC71_struc(n)%ac71(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
            AC71_struc(n)%ac71(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
            AC71_struc(n)%ac71(t)%RootingDepth = GetRootingDepth()
            AC71_struc(n)%ac71(t)%Runoff = GetRunoff()
            AC71_struc(n)%ac71(t)%SaltInfiltr = GetSaltInfiltr()
            AC71_struc(n)%ac71(t)%ScorAT1 = GetScorAT1()
            AC71_struc(n)%ac71(t)%ScorAT2 = GetScorAT2()
            AC71_struc(n)%ac71(t)%Simulation = GetSimulation()
            AC71_struc(n)%ac71(t)%Soil = GetSoil()
            AC71_struc(n)%ac71(t)%SoilLayer = GetSoilLayer()
            AC71_struc(n)%ac71(t)%StageCode = GetStageCode()
            AC71_struc(n)%ac71(t)%StartMode = GetStartMode()
            AC71_struc(n)%ac71(t)%StressLeaf = GetStressLeaf()
            AC71_struc(n)%ac71(t)%StressSFadjNEW = GetStressSFadjNEW()
            AC71_struc(n)%ac71(t)%StressSenescence = GetStressSenescence ()
            AC71_struc(n)%ac71(t)%StressTot = GetStressTot()
            AC71_struc(n)%ac71(t)%SumETo = GetSumETo()
            AC71_struc(n)%ac71(t)%SumGDD = GetSumGDD()
            AC71_struc(n)%ac71(t)%SumGDDPrev = GetSumGDDPrev()
            AC71_struc(n)%ac71(t)%SumGDDcuts = GetSumGDDcuts()
            AC71_struc(n)%ac71(t)%SumInterval = GetSumInterval()
            AC71_struc(n)%ac71(t)%SumKcTop = GetSumKcTop()
            AC71_struc(n)%ac71(t)%SumKcTopStress = GetSumKcTopStress()
            AC71_struc(n)%ac71(t)%SumKci = GetSumKci()
            AC71_struc(n)%ac71(t)%SumWaBal = GetSumWaBal()
            AC71_struc(n)%ac71(t)%Surf0 = GetSurf0()
            AC71_struc(n)%ac71(t)%SurfaceStorage = GetSurfaceStorage()
            AC71_struc(n)%ac71(t)%Tact = GetTact()
            AC71_struc(n)%ac71(t)%TactWeedInfested = GetTactWeedInfested()
            AC71_struc(n)%ac71(t)%Tadj = GetTadj()
            AC71_struc(n)%ac71(t)%TemperatureRecord = GetTemperatureRecord()
            AC71_struc(n)%ac71(t)%TimeSenescence = GetTimeSenescence ()
            AC71_struc(n)%ac71(t)%TotalSaltContent = GetTotalSaltContent()
            AC71_struc(n)%ac71(t)%TotalWaterContent = GetTotalWaterContent()
            AC71_struc(n)%ac71(t)%Tpot = GetTpot()
            AC71_struc(n)%ac71(t)%Transfer = GetTransfer()
            AC71_struc(n)%ac71(t)%WaterTableInProfile = GetWaterTableInProfile()
            AC71_struc(n)%ac71(t)%WeedRCi = GetWeedRCi()
            AC71_struc(n)%ac71(t)%YprevSum = GetYprevSum()
            AC71_struc(n)%ac71(t)%Zeval = GetZeval()
            AC71_struc(n)%ac71(t)%ZiAqua = GetZiAqua()
            AC71_struc(n)%ac71(t)%Ziprev = GetZiprev()
            AC71_struc(n)%ac71(t)%alfaHI = GetalfaHI()
            AC71_struc(n)%ac71(t)%alfaHIAdj = GetalfaHIAdj()
            AC71_struc(n)%ac71(t)%crop = GetCrop()
            AC71_struc(n)%ac71(t)%daynri = GetDayNri()
            AC71_struc(n)%ac71(t)%effectiverain = Geteffectiverain()
            AC71_struc(n)%ac71(t)%endseason = GetEndSeason()
            AC71_struc(n)%ac71(t)%fWeedNoS = GetfWeedNoS()
            AC71_struc(n)%ac71(t)%onset = GetOnset()
            AC71_struc(n)%ac71(t)%simulparam = GetSimulParam()

            !LB to change --> sim period
            if ((LIS_rc%mo .eq. 12) .AND. (LIS_rc%da .eq. 31)) then
                AC71_struc(n)%ac71(t)%InitializeRun = 1
                !call FinalizeRun1(AC71_struc(n)%ac71(t)%irun, GetTheProjectFile(), AC71_struc(n)%ac71(t)%TheProjectType)
                !call FinalizeRun2(AC71_struc(n)%ac71(t)%irun, AC71_struc(n)%ac71(t)%TheProjectType)
                AC71_struc(n)%ac71(t)%irun = AC71_struc(n)%ac71(t)%irun + 1
            end if

            ![ 1] output variable: smc (unit=m^3 m-3 ). ***  volumetric soil moisture, ice + liquid 
            do i=1, AC71_struc(n)%ac71(t)%NrCompartments
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = AC71_struc(n)%ac71(t)%smc(i),  &
                                                    vlevel=i, unit="kg m-2", direction="-", surface_type = LIS_rc%lsm_index)
                call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_SOILMOIST, value = AC71_struc(n)%ac71(t)%smc(i), &
                                                    vlevel=i, unit="m^3 m-3", direction="-", surface_type = LIS_rc%lsm_index)
            end do
            ![ 4] output variable: biomass (unit=t/ha). ***  leaf area index 
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71BIOMASS, value = real(AC71_struc(n)%ac71(t)%SumWaBal%Biomass,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            !call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71BIOMASS, value = real(AC71_struc(n)%ac71(t)%SumWaBal%Biomass,kind=sp), &
            !                                  vlevel=1, unit="t h-1", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_CCiPrev, value = real(AC71_struc(n)%ac71(t)%CCiPrev,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Irrigation, value = real(AC71_struc(n)%ac71(t)%Irrigation,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71RootZoneWC_Actual, value = real(AC71_struc(n)%ac71(t)%RootZoneWC_Actual,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71RootZoneWC_WP, value = real(AC71_struc(n)%ac71(t)%RootZoneWC_WP,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71RootZoneWC_FC, value = real(AC71_struc(n)%ac71(t)%RootZoneWC_FC,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Tact, value = real(AC71_struc(n)%ac71(t)%Tact,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Eact, value = real(AC71_struc(n)%ac71(t)%Eact,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71ETo, value = real(AC71_struc(n)%ac71(t)%eto,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71RootingDepth, value = real(AC71_struc(n)%ac71(t)%RootingDepth,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71CCiActual, value = real(AC71_struc(n)%ac71(t)%CCiActual,kind=sp), &
                                                vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Tmin, value = real(tmp_tmin,kind=sp), &
                                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Tmax, value = real(tmp_tmax,kind=sp), &
                                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Rain, value = real(tmp_precip,kind=sp), &
                                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

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
