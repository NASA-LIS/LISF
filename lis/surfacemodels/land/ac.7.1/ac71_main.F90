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

!   ! AC module imports
    use ac_global, only: &
                    DegreesDay,&
                    GenerateDepthMode_ToFC, &
                    GenerateTimeMode_AllRAW, &
                    GetCCiActual,&
                    GetCCiTopEarlySen,&
                    GetCCiprev,&
                    GetCompartment,&
                    GetCompartment_theta,&
                    GetCrop,&
                    GetCrop_Day1,&
                    GetCrop_DayN,&
                    SetCrop_Day1,&
                    SetCrop_DayN,&
                    SetCrop_DaysToGermination,&
                    SetCrop_DaysToMaxRooting,&
                    SetCrop_DaysToFlowering,&
                    SetCrop_DaysToHarvest,&
                    SetCrop_DaysTosenescence,&
                    SetCrop_DaysToCCini,&
                    SetCrop_DaysToFullCanopy,&
                    SetCrop_DaysToFullCanopySF,&
                    SetCrop_DaysToHIo,&
                    GetCrop_DaysToGermination,&
                    GetCrop_DaysToMaxRooting,&
                    GetCrop_DaysToFlowering,&
                    GetCrop_DaysToHarvest,&
                    GetCrop_DaysTosenescence,&
                    GetCrop_DaysToCCini,&
                    GetCrop_DaysToFullCanopy,&
                    GetCrop_DaysToFullCanopySF,&
                    GetCrop_DaysToHIo,&
                    GetCrop_CGC, &
                    GetCrop_GDDCGC,&
                    GetDaySubmerged,&
                    GetEact,&
                    GetECstorage,& 
                    GetETo,&
                    GetIrriAfterSeason, &
                    GetIrriBeforeSeason, &
                    GetIrriECw, &
                    GetIrrigation,&
                    GetIrriMode,&
                    GetManagement,&
                    GetPreDay,&
                    GetRootingDepth,&
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
                    GetSimulParam,&
                    GetSimulParam_GDDMethod,&
                    GetCrop_ModeCycle,&
                    GetSimulation_DelayedDays,&
                    GetSimulation,&
                    GetSimulation_EvapLimitON, &
                    GetSimulation_SumGDD,&
                    GetSimulation_SumGDDfromDay1,&
                    GetSimulation_SWCtopSoilConsidered, &
                    GetSimulation_FromDayNr, &
                    GetSimulation_ToDayNr, &
                    GetSoil,&
                    GetSoilLayer,&
                    GetSumWaBal,&
                    GetSurfaceStorage,&
                    GetTact,&
                    GetTactWeedInfested,&
                    GetTmax,& 
                    GetTmin,& 
                    GetTpot,&
                    GetCrop_GDDaysToGermination,&
                    GetCrop_GDDaysToMaxRooting,&
                    GetCrop_GDDaysToHarvest,&
                    GetCrop_GDDaysToGermination,&
                    GetCrop_GDDaysToFlowering,&
                    GetCrop_GDDaysToSenescence,&
                    IrriMode_Generate,&
                    IrriMode_Inet,&
                    IrriMode_Manual,&
                    IrriMode_NoIrri,&
                    SetCCiActual,&
                    SetCCiTopEarlySen,&
                    SetCCiprev,&
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
                    SetECstorage,& 
                    SetETo,&
                    SetIrriAfterSeason, &
                    SetIrriBeforeSeason, &
                    SetIrriECw, &
                    SetIrriMethod, &
                    SetIrriMode, &
                    SetGenerateDepthMode, &
                    SetGenerateTimeMode, &
                    SetManagement,&
                    SetOutputAggregate,&
                    SetOutDaily,&
                    SetOut3Prof,&
                    SetPart1Mult,&
                    SetPart2Eval,&
                    SetPreDay,&
                    SetRain,& 
                    SetRootingDepth,&
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
                    SetSimulation,&
                    SetSimulation_EvapLimitON, &
                    SetSimulation_SumGDD, &
                    SetSimulation_SumGDDfromDay1,&
                    SetSimulation_SWCtopSoilConsidered, &
                    SetSimulParam,&
                    SetSoil,&
                    SetSoilLayer,&
                    SetSumWaBal,&
                    SetSurfaceStorage,&
                    SetTact,&
                    SetTactWeedInfested,&
                    SetTmax,& 
                    SetTmin,&
                    SetTpot,&
                    SetTminRun_i, &
                    SetTminRun, &
                    SetTmaxRun_i, &
                    SetTmaxRun, &
                    GetTminRun_i, &
                    GetTminRun, &
                    GetTmaxRun_i, &
                    GetTmaxRun, &
                    undef_int
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
                    fIrri_close, &
                    fIrri_open,&
                    firri_read,&
                    GetBin,&
                    GetBout,&
                    GetCCiActualWeedInfested,&
                    GetCCoTotal,&
                    GetCCxCropWeedsNoSFstress,&
                    GetCCxTotal,&
                    GetCCxWitheredTpotNoS,&
                    GetCDCTotal,&
                    GetCoeffb0,&
                    GetCoeffb0Salt,&
                    GetCoeffb1,&
                    GetCoeffb1Salt,&
                    GetCoeffb2,&
                    GetCoeffb2Salt,&
                    GetCrop_Tbase, &
                    GetCrop_Tupper, &
                    GetDayFraction,&
                    GetDayNri,&
                    GetSimulation_FromDayNr,&
                    GetGDDCDCTotal,&
                    GetGDDTadj,&
                    GetGDDayFraction,&
                    GetGDDayi,&
                    GetHItimesAT,&
                    GetHItimesAT1,&
                    GetHItimesAT2,&
                    GetHItimesBEF,&
                    GetIrriInfoRecord1, &
                    GetIrriInfoRecord1_NoMoreInfo,&
                    GetIrriInfoRecord2, &
                    GetNoMoreCrop,&
                    GetPreviousStressLevel,&
                    GetScorAT1,&
                    GetScorAT2,&
                    GetStressLeaf,&
                    GetStressSFadjNEW,&
                    GetStressSenescence ,&
                    GetStressTot,&
                    GetSumGDD,&
                    GetSumGDDcuts,&
                    GetSumGDDPrev,&
                    GetSumInterval,&
                    GetSumKcTop,&
                    GetSumKcTopStress,&
                    GetSumKci,&
                    GetTadj,&
                    GetTheProjectFile,&
                    GetTimeSenescence ,&
                    GetWeedRCi,&
                    GetZiprev,&
                    GetalfaHI,&
                    GetalfaHIAdj,&
                    InitializeRunPart1, &
                    InitializeSimulation, &
                    InitializeSimulationRunPart2, &
                    SetBin,&
                    SetBout,&
                    SetCCiActualWeedInfested,&
                    SetCCoTotal,&
                    SetCCxCropWeedsNoSFstress,&
                    SetCCxTotal,&
                    SetCCxWitheredTpotNoS,&
                    SetCDCTotal,&
                    SetCGCref,&
                    SetCoeffb0,&
                    SetCoeffb0Salt,&
                    SetCoeffb1,&
                    SetCoeffb1Salt,&
                    SetCoeffb2,&
                    SetCoeffb2Salt,&
                    SetDayFraction,&
                    SetDayNri,&
                    SetGDDCDCTotal,&
                    SetGDDCGCref ,&
                    SetGDDTadj,&
                    SetGDDayFraction,&
                    SetGDDayi,&
                    SetHItimesAT,&
                    SetHItimesAT1,&
                    SetHItimesAT2,&
                    SetHItimesBEF,&
                    SetIrriInfoRecord1, &
                    SetIrriInfoRecord2, &
                    SetNextSimFromDayNr ,&
                    SetNoMoreCrop,&
                    SetNoYear,&
                    SetPreviousStressLevel,&
                    SetScorAT1,&
                    SetScorAT2,&
                    SetStartMode,&
                    SetStressLeaf,&
                    SetStressSFadjNEW,&
                    SetStressSenescence ,&
                    SetStressTot,&
                    SetSumGDD,&
                    SetSumGDDPrev,&
                    SetSumGDDcuts,&
                    SetSumInterval,&
                    SetSumKcTop,&
                    SetSumKcTopStress,&
                    SetSumKci,&
                    SetTadj,&
                    SetTimeSenescence ,&
                    SetWeedRCi,&
                    SetZiprev,&
                    SetalfaHI,&
                    SetalfaHIAdj,&
                    GetPlotVarCrop,&
                    SetPlotVarCrop,&
                    GetfWeedNoS,&
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
            ! Required vars for restart
            call SetalfaHI(REAL(AC71_struc(n)%ac71(t)%alfaHI, 8))
            call SetalfaHIAdj(REAL(AC71_struc(n)%ac71(t)%alfaHIAdj, 8))
            call SetBin(REAL(AC71_struc(n)%ac71(t)%Bin, 8))
            call SetBout(REAL(AC71_struc(n)%ac71(t)%Bout, 8))
            call SetCCiActual(REAL(AC71_struc(n)%ac71(t)%CCiActual, 8))
            call SetCCiActualWeedInfested(REAL(AC71_struc(n)%ac71(t)%CCiActualWeedInfested, 8))
            call SetCCiprev(REAL(AC71_struc(n)%ac71(t)%CCiprev, 8))
            call SetCCiTopEarlySen(REAL(AC71_struc(n)%ac71(t)%CCiTopEarlySen, 8))
            call SetCCxWitheredTpotNoS(REAL(AC71_struc(n)%ac71(t)%CCxWitheredTpotNoS, 8))
            call SetCompartment(AC71_struc(n)%ac71(t)%Compartment)
            do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                    call SetCompartment_theta(l,REAL(AC71_struc(n)%ac71(t)%smc(l),8))
            enddo
            call SetDayFraction(REAL(AC71_struc(n)%ac71(t)%DayFraction, 8))
            call SetDayNri(AC71_struc(n)%ac71(t)%daynri)
            call SetDaySubmerged(AC71_struc(n)%ac71(t)%DaySubmerged)
            call SetECstorage(REAL(AC71_struc(n)%ac71(t)%ECstorage, 8))
            call SetHItimesAT(REAL(AC71_struc(n)%ac71(t)%HItimesAT, 8))
            call SetHItimesAT1(REAL(AC71_struc(n)%ac71(t)%HItimesAT1, 8))
            call SetHItimesAT2(REAL(AC71_struc(n)%ac71(t)%HItimesAT2, 8))
            call SetHItimesBEF(REAL(AC71_struc(n)%ac71(t)%HItimesBEF, 8))
            call SetPlotVarCrop(AC71_struc(n)%ac71(t)%PlotVarCrop) !Not in restart yet
            call SetPreviousStressLevel(int(AC71_struc(n)%ac71(t)%PreviousStressLevel,kind=int32))
            call SetRootingDepth(REAL(AC71_struc(n)%ac71(t)%RootingDepth, 8)) ! Not in retstart, likely not needed
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
            call SetScorAT1(REAL(AC71_struc(n)%ac71(t)%ScorAT1, 8))
            call SetScorAT2(REAL(AC71_struc(n)%ac71(t)%ScorAT2, 8))
            call SetStressLeaf(REAL(AC71_struc(n)%ac71(t)%StressLeaf, 8))
            call SetStressSenescence(REAL(AC71_struc(n)%ac71(t)%StressSenescence, 8))
            call SetStressSFadjNEW(int(AC71_struc(n)%ac71(t)%StressSFadjNEW,kind=int32))
            call SetStressTot(AC71_struc(n)%ac71(t)%StressTot)
            call SetSumGDDcuts(REAL(AC71_struc(n)%ac71(t)%SumGDDcuts, 8))
            call SetSumGDD(REAL(AC71_struc(n)%ac71(t)%SumGDD, 8)) ! Not yet in restart
            call SetSumGDDPrev(REAL(AC71_struc(n)%ac71(t)%SumGDDPrev, 8))
            call SetSumInterval(AC71_struc(n)%ac71(t)%SumInterval)
            call SetSumKci(REAL(AC71_struc(n)%ac71(t)%SumKci, 8))
            call SetSumKcTopStress(REAL(AC71_struc(n)%ac71(t)%SumKcTopStress, 8))
            call SetSumWaBal(AC71_struc(n)%ac71(t)%SumWaBal)
            call SetSurfaceStorage(REAL(AC71_struc(n)%ac71(t)%SurfaceStorage, 8))
            call SetTact(REAL(AC71_struc(n)%ac71(t)%Tact, 8))
            call SetTactWeedInfested(REAL(AC71_struc(n)%ac71(t)%TactWeedInfested, 8))
            call SetTadj(AC71_struc(n)%ac71(t)%Tadj)
            call SetTimeSenescence(REAL(AC71_struc(n)%ac71(t)%TimeSenescence, 8))
            call SetTpot(REAL(AC71_struc(n)%ac71(t)%Tpot, 8))
            call SetWeedRCi(REAL(AC71_struc(n)%ac71(t)%WeedRCi, 8))
            call SetZiprev(REAL(AC71_struc(n)%ac71(t)%Ziprev, 8))
            call SetfWeedNoS(REAL(AC71_struc(n)%ac71(t)%fWeedNoS, 8)) ! not in restart yet

            ! Set logicals global variables required in restart
            if(AC71_struc(n)%ac71(t)%NoMoreCrop.eq.1)then
                call SetNoMoreCrop(.true.)
            else
                call SetNoMoreCrop(.false.)
            endif
            if(AC71_struc(n)%ac71(t)%Simulation%EvapLimitON.eq.1)then
                call SetSimulation_EvapLimitON(.true.)
            else
                call SetSimulation_EvapLimitON(.false.)
            endif
            if(AC71_struc(n)%ac71(t)%Simulation%SWCtopSoilConsidered.eq.1)then
                call SetSimulation_SWCtopSoilConsidered(.true.)
            else
                call SetSimulation_SWCtopSoilConsidered(.false.)
            endif

            ! Set in Initialize (not needed for restart)
            call SetCCoTotal(REAL(AC71_struc(n)%ac71(t)%CCoTotal, 8)) 
            call SetCCxCropWeedsNoSFstress(REAL(AC71_struc(n)%ac71(t)%CCxCropWeedsNoSFstress, 8)) 
            call SetCCxTotal(REAL(AC71_struc(n)%ac71(t)%CCxTotal, 8))
            call SetCDCTotal(REAL(AC71_struc(n)%ac71(t)%CDCTotal, 8))
            call SetCoeffb0(REAL(AC71_struc(n)%ac71(t)%Coeffb0, 8))
            call SetCoeffb0Salt(REAL(AC71_struc(n)%ac71(t)%Coeffb0Salt, 8))
            call SetCoeffb1(REAL(AC71_struc(n)%ac71(t)%Coeffb1, 8))
            call SetCoeffb1Salt(REAL(AC71_struc(n)%ac71(t)%Coeffb1Salt, 8))
            call SetCoeffb2(REAL(AC71_struc(n)%ac71(t)%Coeffb2, 8))
            call SetCoeffb2Salt(REAL(AC71_struc(n)%ac71(t)%Coeffb2Salt, 8))
            call SetCrop(AC71_struc(n)%ac71(t)%crop)
            call SetGDDayFraction(REAL(AC71_struc(n)%ac71(t)%GDDayFraction, 8))
            call SetGDDCDCTotal(REAL(AC71_struc(n)%ac71(t)%GDDCDCTotal, 8))
            call SetGDDTadj(AC71_struc(n)%ac71(t)%GDDTadj)
            call SetManagement(AC71_struc(n)%ac71(t)%Management)
            call SetSimulation(AC71_struc(n)%ac71(t)%Simulation)
            call SetSimulParam(AC71_struc(n)%ac71(t)%SimulParam)
            call SetSoil(AC71_struc(n)%ac71(t)%Soil)
            call SetSoilLayer(AC71_struc(n)%ac71(t)%SoilLayer) 
            call SetSumKcTop(REAL(AC71_struc(n)%ac71(t)%SumKcTop, 8))

            !! Fixed vars
            call SetCGCref(GetCrop_CGC()) ! Make sure crop is set before
            call SetGDDCGCref(GetCrop_GDDCGC()) ! Make sure crop is set before
            call SetNextSimFromDayNr(int(-9, kind=int32)) !Always undef_int in AquaCrop... Check src for v7.2
            call SetNoYear(.false.)
            call SetOutputAggregate(int(0,kind=int8))
            call SetOut3Prof(.true.)
            call SetOutDaily(.true.) ! AVoid writing daily output in console?
            call SetPart1Mult(.false.) 
            call SetPart2Eval(.false.)
            call SetPreDay(.true.) ! set to false in InitializeSettings
            call SetStartMode(.false.) ! Overwritten to .true. in InitalizeRunPart1
            

            ! Optional vars -> for later implementations
            ! Groundwater (not tested in LIS)
            !call SetSimulParam_ConstGwt(.true.) by default
            !call SetGwTable(AC71_struc(n)%ac71(t)%GwTable)
            !call SetSimulParam(AC71_struc(n)%ac71(t)%simulparam)
            !call SetZiAqua(AC71_struc(n)%ac71(t)%ZiAqua)
            !call SetEciAqua(AC71_struc(n)%ac71(t)%ECiAqua)
            !call SetWaterTableInProfile(AC71_struc(n)%ac71(t)%WaterTableInProfile)

            ! .OBS not used in LIS
            !call SetDayNr1Eval(AC71_struc(n)%ac71(t)%DayNr1Eval)
            !call SetDayNrEval(AC71_struc(n)%ac71(t)%DayNrEval)
            !call SetLineNrEval(int(AC71_struc(n)%ac71(t)%LineNrEval,kind=int32))

            ! Perennial mode
            !call SetCutInfoRecord1(AC71_struc(n)%ac71(t)%CutInfoRecord1)
            !call SetCutInfoRecord2(AC71_struc(n)%ac71(t)%CutInfoRecord2)
            !call SetDayLastCut(AC71_struc(n)%ac71(t)%DayLastCut)
            !call SetNrCut(AC71_struc(n)%ac71(t)%NrCut)
            !call SetPerennialPeriod(AC71_struc(n)%ac71(t)%PerennialPeriod)
            !call SetBprevSum(AC71_struc(n)%ac71(t)%BprevSum)
            !call SetTransfer(AC71_struc(n)%ac71(t)%Transfer)
            !call SetOnset(AC71_struc(n)%ac71(t)%onset)

            ! Set climate variables
            call SetRain(real(tmp_precip,kind=dp))
            call SetTmin(real(tmp_tmin,kind=dp))
            call SetTmax(real(tmp_tmax,kind=dp))
            call SetETo(real(tmp_eto,kind=dp))
            
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

                ! GDD mode: set end of stages to end of simulation period if undefined (-9)
                ! For calndar days, stages should already be < 365
                simul_len = GetSimulation_ToDayNr() - GetCrop_Day1() + 1
                ! Germination
                if (GetCrop_DaysToGermination().eq.undef_int) then
                    call SetCrop_DaysToGermination(simul_len)
                endif
                ! Flowering (only for Tuber and grain)
                if (GetCrop_DaysToFlowering().eq.undef_int) then
                    call SetCrop_DaysToFlowering(simul_len)
                endif
                ! MaxRooting
                if (GetCrop_DaysToMaxRooting().eq.undef_int) then
                    call SetCrop_DaysToMaxRooting(simul_len)
                endif
                ! Senescence
                if (GetCrop_DaysToSenescence().eq.undef_int) then
                    call SetCrop_DaysToSenescence(simul_len)
                endif
                ! Harvest
                if (GetCrop_DaysToHarvest().eq.undef_int) then
                    call SetCrop_DaysToHarvest(simul_len)
                    call SetCrop_DayN(GetCrop_Day1() + GetCrop_DaysToHarvest() - 1)
                endif
                ! CCini
                if (GetCrop_DaysToCCini().eq.undef_int) then
                    call SetCrop_DaysToCCini(simul_len)
                endif
                ! FullCanopy
                if (GetCrop_DaysToFullCanopy().eq.undef_int) then
                    call SetCrop_DaysToFullCanopy(simul_len)
                endif
                ! FullCanopySF
                if (GetCrop_DaysToFullCanopySF().eq.undef_int) then
                    call SetCrop_DaysToFullCanopySF(simul_len)
                endif
                ! HIo
                if (GetCrop_DaysToHIo().eq.undef_int) then
                    call SetCrop_DaysToHIo(simul_len)
                endif

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
            tmp_wpi = REAL(AC71_struc(n)%ac71(t)%WPi,8)
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
            
            AC71_struc(n)%ac71(t)%alfaHI = GetalfaHI()
            AC71_struc(n)%ac71(t)%alfaHIAdj = GetalfaHIAdj()
            AC71_struc(n)%ac71(t)%Bin = GetBin()
            AC71_struc(n)%ac71(t)%Bout = GetBout()
            AC71_struc(n)%ac71(t)%CCiActual = GetCCiActual()
            AC71_struc(n)%ac71(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
            AC71_struc(n)%ac71(t)%CCiprev = GetCCiprev()
            AC71_struc(n)%ac71(t)%CCiTopEarlySen = GetCCiTopEarlySen()
            AC71_struc(n)%ac71(t)%CCoTotal = GetCCoTotal()
            AC71_struc(n)%ac71(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
            AC71_struc(n)%ac71(t)%CCxTotal = GetCCxTotal()
            AC71_struc(n)%ac71(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
            AC71_struc(n)%ac71(t)%CDCTotal = GetCDCTotal()
            AC71_struc(n)%ac71(t)%Coeffb0 = GetCoeffb0()
            AC71_struc(n)%ac71(t)%Coeffb0Salt = GetCoeffb0Salt()
            AC71_struc(n)%ac71(t)%Coeffb1 = GetCoeffb1()
            AC71_struc(n)%ac71(t)%Coeffb1Salt = GetCoeffb1Salt()
            AC71_struc(n)%ac71(t)%Coeffb2 = GetCoeffb2()
            AC71_struc(n)%ac71(t)%Coeffb2Salt = GetCoeffb2Salt()
            AC71_struc(n)%ac71(t)%Compartment = GetCompartment()
            AC71_struc(n)%ac71(t)%Crop = GetCrop()
            AC71_struc(n)%ac71(t)%DayFraction = GetDayFraction()
            AC71_struc(n)%ac71(t)%daynri = GetDayNri()
            AC71_struc(n)%ac71(t)%DaySubmerged = GetDaySubmerged()
            AC71_struc(n)%ac71(t)%Eact = GetEact()
            AC71_struc(n)%ac71(t)%ECstorage = GetECstorage()
            AC71_struc(n)%ac71(t)%GDDayFraction = GetGDDayFraction()
            AC71_struc(n)%ac71(t)%GDDCDCTotal = GetGDDCDCTotal()
            AC71_struc(n)%ac71(t)%GDDTadj = GetGDDTadj()
            AC71_struc(n)%ac71(t)%HItimesAT = GetHItimesAT()
            AC71_struc(n)%ac71(t)%HItimesAT1 = GetHItimesAT1()
            AC71_struc(n)%ac71(t)%HItimesAT2 = GetHItimesAT2()
            AC71_struc(n)%ac71(t)%HItimesBEF = GetHItimesBEF()
            AC71_struc(n)%ac71(t)%Irrigation = GetIrrigation()
            AC71_struc(n)%ac71(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
            AC71_struc(n)%ac71(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
            AC71_struc(n)%ac71(t)%Management = GetManagement()
            AC71_struc(n)%ac71(t)%PlotVarCrop = GetPlotVarCrop()
            AC71_struc(n)%ac71(t)%PreDay = GetPreDay() ! be careful with restart
            AC71_struc(n)%ac71(t)%PreviousStressLevel = GetPreviousStressLevel()
            AC71_struc(n)%ac71(t)%RootingDepth = GetRootingDepth()
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
            AC71_struc(n)%ac71(t)%ScorAT1 = GetScorAT1()
            AC71_struc(n)%ac71(t)%ScorAT2 = GetScorAT2()
            AC71_struc(n)%ac71(t)%Simulation = GetSimulation()
            AC71_struc(n)%ac71(t)%SimulParam = GetSimulParam()
            AC71_struc(n)%ac71(t)%Soil = GetSoil()
            AC71_struc(n)%ac71(t)%SoilLayer = GetSoilLayer()
            AC71_struc(n)%ac71(t)%StressLeaf = GetStressLeaf()
            AC71_struc(n)%ac71(t)%fWeedNoS = GetfWeedNoS()
            AC71_struc(n)%ac71(t)%StressSenescence = GetStressSenescence()
            AC71_struc(n)%ac71(t)%StressSFadjNEW = GetStressSFadjNEW()
            AC71_struc(n)%ac71(t)%StressTot = GetStressTot()
            AC71_struc(n)%ac71(t)%SumGDDcuts = GetSumGDDcuts()
            AC71_struc(n)%ac71(t)%SumGDD = GetSumGDD()
            AC71_struc(n)%ac71(t)%SumGDDPrev = GetSumGDDPrev()
            AC71_struc(n)%ac71(t)%SumInterval = GetSumInterval()
            AC71_struc(n)%ac71(t)%SumKci = GetSumKci()
            AC71_struc(n)%ac71(t)%SumKcTop = GetSumKcTop()
            AC71_struc(n)%ac71(t)%SumKcTopStress = GetSumKcTopStress()
            AC71_struc(n)%ac71(t)%SumWaBal = GetSumWaBal()
            AC71_struc(n)%ac71(t)%SurfaceStorage = GetSurfaceStorage()
            AC71_struc(n)%ac71(t)%Tact = GetTact()
            AC71_struc(n)%ac71(t)%TactWeedInfested = GetTactWeedInfested()
            AC71_struc(n)%ac71(t)%Tadj = GetTadj()
            AC71_struc(n)%ac71(t)%TimeSenescence = GetTimeSenescence()
            AC71_struc(n)%ac71(t)%Tpot = GetTpot()
            AC71_struc(n)%ac71(t)%WeedRCi = GetWeedRCi()
            AC71_struc(n)%ac71(t)%Ziprev = GetZiprev()

            !logicals are stored as integers for restart file
            !NoMoreCrop
            if(GetNoMoreCrop()) then
                AC71_struc(n)%ac71(t)%NoMoreCrop = 1
            else
                AC71_struc(n)%ac71(t)%NoMoreCrop = 0
            endif
            !Simulation_EvapLimitON
            if(GetSimulation_EvapLimitON()) then
                AC71_struc(n)%ac71(t)%Simulation%EvapLimitON = 1
            else
                AC71_struc(n)%ac71(t)%Simulation%EvapLimitON = 0
            endif
            !Simulation_SWCtopSoilConsidered
            if(GetSimulation_SWCtopSoilConsidered()) then
                AC71_struc(n)%ac71(t)%Simulation%SWCtopSoilConsidered = 1
            else
                AC71_struc(n)%ac71(t)%Simulation%SWCtopSoilConsidered = 0
            endif

            ! Check for end of simulation period 
            ! (DayNri - 1 because DayNri is alreayd for next day)
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
            ![ 11] output variable: AC71Tmin (unit=K).  *** daily minimum temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_Tmin, value = real(tmp_tmin+LIS_CONST_TKFRZ,kind=sp), &
                                    vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 12] output variable: AC71Tmax (unit=K).  *** daily maximum temperature
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_Tmax, value = real(tmp_tmax+LIS_CONST_TKFRZ,kind=sp), &
                                    vlevel=1, unit="K", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 13] output variable: rainf (unit=kg m-2 s-1).  *** precipitation rate
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RAINF, value = real(tmp_precip/86400.,kind=sp), &
                                    vlevel=1, unit="kg m-2 s-1", direction="DN", surface_type = LIS_rc%lsm_index)
            ![ 14] output variable: yield (unit=t ha-1).  *** yield
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_Yield, value = real(AC71_struc(n)%ac71(t)%SumWaBal%YieldPart,kind=sp), &
                                    vlevel=1, unit="t ha-1", direction="-", surface_type = LIS_rc%lsm_index)
            ![ 15] output variable: irrigation (unit=mm).  *** irrigation
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Irrigation, value = real(AC71_struc(n)%ac71(t)%Irrigation,kind=sp), &
                                    vlevel=1, unit="mm", direction="-", surface_type = LIS_rc%lsm_index)

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
