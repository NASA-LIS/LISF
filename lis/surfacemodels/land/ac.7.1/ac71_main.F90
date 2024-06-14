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
   !use other modules
    use ESMF
    use LIS_routingMod, only : LIS_runoff_state

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
                    GetIrrigation, &
                    GetManagement,&
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
                    IrriMethod_MSprinkler, &
                    IrriMode_Generate, &
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
                    SetPart1Mult,&
                    SetPart2Eval,&
                    SetPreDay,&
                    SetRain,& 
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
                    GetTmaxRun

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
    logical              :: alarmCheck
    real                 :: TminRun_i, TmaxRun_i
    real                 :: sumGDD, sumGDD2

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
        !if ((AC71_struc(n)%ac71(1)%InitializeRun.eq.1).and.(AC71_struc(n)%GDD_Mode.eq.1)) then
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
            !Required vars to be set before simulation
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
            call SetPreviousStressLevel(int(AC71_struc(n)%ac71(t)%PreviousStressLevel,kind=int32))
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
            call SetSumGDD(REAL(AC71_struc(n)%ac71(t)%SumGDD, 8))
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
            call SetfWeedNoS(REAL(AC71_struc(n)%ac71(t)%fWeedNoS, 8))

            ! MB apparently needed for GDD runs
            call SetPlotVarCrop(AC71_struc(n)%ac71(t)%PlotVarCrop)

            if (.not. ((LIS_rc%mo .eq. AC71_struc(n)%Sim_AnnualStartMonth) &
                .AND. (LIS_rc%da .eq. AC71_struc(n)%Sim_AnnualStartDay))) then !make it flex
                ! Set logicals
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
                ! Can be false when sim is initialized
                call SetPreDay(.true.) ! set to false in InitializeSettings
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
            call SetSoil(AC71_struc(n)%ac71(t)%Soil)
            call SetSoilLayer(AC71_struc(n)%ac71(t)%SoilLayer) 
            call SetSumKcTop(REAL(AC71_struc(n)%ac71(t)%SumKcTop, 8))

            !! Fixed vars
            call SetCGCref(GetCrop_CGC()) ! Make sure crop is set before
            call SetGDDCGCref(GetCrop_GDDCGC()) ! Make sure crop is set before
            call SetNextSimFromDayNr(int(-9, kind=int32)) !Always undef_int in AquaCrop... Check src for v7.2
            call SetNoYear(.false.)
            call SetOutputAggregate(int(0,kind=int8)) ! Avoid writing out daily results in the console
            call SetPart1Mult(.false.) 
            call SetPart2Eval(.false.)
            call SetStartMode(.false.) ! Overwritten to .true. in InitalizeRunPart1
            !old, needed?: call SetSumGDDPrev(GetSimulation_SumGDD()) ! Make sure that Simulation is set before

            !! If irrigation ON
            if(LIS_rc%irrigation_type.ne."none") then
                if(trim(LIS_rc%irrigation_type).eq."Sprinkler") then
                    call SetIrriAfterSeason(AC71_struc(n)%ac71(t)%IrriAfterSeason)
                    call SetIrriBeforeSeason(AC71_struc(n)%ac71(t)%IrriBeforeSeason)
                    call SetIrriECw(AC71_struc(n)%ac71(t)%IrriECw) 
                    call SetIrriInfoRecord1(AC71_struc(n)%ac71(t)%IrriInfoRecord1)
                    call SetIrriInfoRecord2(AC71_struc(n)%ac71(t)%IrriInfoRecord2)
                    call SetIrriMethod(IrriMethod_MSprinkler)
                    call SetIrriMode(IrriMode_Generate)
                    call SetGenerateDepthMode(GenerateDepthMode_ToFC)
                    call SetGenerateTimeMode(GenerateTimeMode_AllRAW)
                else ! Other options can be implemented later
                    write(LIS_logunit, *) trim(LIS_rc%irrigation_type), " irrigation type not compatible with AquaCrop.7.1"
                    call LIS_endrun()
                endif
            endif

            call SetCrop_DaysToGermination(AC71_struc(n)%ac71(t)%Crop_DaysToGermination)
            call SetCrop_DaysToFlowering(AC71_struc(n)%ac71(t)%Crop_DaysToFlowering)
            call SetCrop_DaysToMaxRooting(AC71_struc(n)%ac71(t)%Crop_DaysToMaxRooting)
            call SetCrop_DaysToSenescence(AC71_struc(n)%ac71(t)%Crop_DaysToSenescence)
            call SetCrop_DaysToHarvest(AC71_struc(n)%ac71(t)%Crop_DaysToHarvest)
            call SetCrop_DaysToCCini(AC71_struc(n)%ac71(t)%Crop_DaysToCCini)
            call SetCrop_DaysToFullCanopy(AC71_struc(n)%ac71(t)%Crop_DaysToFullCanopy)
            call SetCrop_DaysToFullCanopySF(AC71_struc(n)%ac71(t)%Crop_DaysToFullCanopySF)
            call SetCrop_DaysToHIo(AC71_struc(n)%ac71(t)%Crop_DaysToHIo)

            ! Set climate variables
            call SetRain(real(tmp_precip,kind=dp))
            call SetTmin(real(tmp_tmin,kind=dp))
            call SetTmax(real(tmp_tmax,kind=dp))
            call SetETo(real(tmp_eto,kind=dp))
            
            ! SumGDD calculation needed only for second day when not done within InitializeSimulationRunPart2
            if (GetDayNri()>GetSimulation_FromDayNr()) then
            !if (AC71_struc(n)%ac71(t)%InitializeRun.eq.0) then !make it flex
            ! Sum of GDD at end of first day ! Wait for GDD implementation from Michel
            call SetGDDayi(DegreesDay(GetCrop_Tbase(), GetCrop_Tupper(), GetTmin(), &
                    GetTmax(), GetSimulParam_GDDMethod()))
            if (GetDayNri() >= GetCrop_Day1()) then
                call SetSimulation_SumGDD(GetSimulation_SumGDD() + GetGDDayi())
                call SetSimulation_SumGDDfromDay1(GetSimulation_SumGDDfromDay1() + &
                    GetGDDayi())
            end if
            end if
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
                call SetTminRun(AC71_struc(n)%Trecord(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index)%Tmin_record)
                call SetTmaxRun(AC71_struc(n)%Trecord(LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index)%Tmax_record)
                ! InitializeRunPart1
                call InitializeRunPart1(int(AC71_struc(n)%ac71(t)%irun,kind=int8), AC71_struc(n)%ac71(t)%TheProjectType)
                call InitializeSimulationRunPart2()
                AC71_struc(n)%ac71(t)%HarvestNow = .false.
                if(LIS_rc%irrigation_type.ne."none") then
                    call fIrri_close() !LB check if problem with specific irrigation file
                endif
                AC71_struc(n)%ac71(t)%InitializeRun = 0
            end if

            ! Initialize for new crop cycle

            ! Reset Crop_DaysTo* to allow that members reach stages at different days
            phenological_stages_ensemble = .false.

            if (phenological_stages_ensemble) then
                if (GetDayNri() == GetCrop_Day1()) then
                !call setCrop_DayN(GetSimulation_ToDayNr())
                AC71_struc(n)%ac71(t)%germ_reached = .false.
                AC71_struc(n)%ac71(t)%harv_reached = .false.
                AC71_struc(n)%ac71(t)%flowr_reached = .false.
                AC71_struc(n)%ac71(t)%MaxR_reached = .false.
                AC71_struc(n)%ac71(t)%Sene_reached = .false.
                ! Initialize to end of the year but with one day difference 
                ! due to internal AC if statements
                !AC71_struc(n)%ac71(t)%Crop_DaysToGermination = 361
                !AC71_struc(n)%ac71(t)%Crop_DaysToFlowering = 362
                !AC71_struc(n)%ac71(t)%Crop_DaysToMaxRooting = 363
                !AC71_struc(n)%ac71(t)%Crop_DaysToSenescence = 364
                !AC71_struc(n)%ac71(t)%Crop_DaysToHarvest= 365

                !find calendar days for crop stages
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToGermination()) .and.  (.not. AC71_struc(n)%ac71(t)%germ_reached)) then ! from sow
                    AC71_struc(n)%ac71(t)%Crop_DaysToGermination = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%germ_reached = .true.
                end if
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToMaxRooting()) .and. (.not. AC71_struc(n)%ac71(t)%maxR_reached)) then ! from sowin
                    AC71_struc(n)%ac71(t)%Crop_DaysToMaxRooting = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%maxR_reached = .true.
                end if
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToFlowering()) .and.  (.not. AC71_struc(n)%ac71(t)%flowr_reached)) then ! from sowi
                    AC71_struc(n)%ac71(t)%Crop_DaysToFlowering = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%flowr_reached = .true.
                end if
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToSenescence()) .and. (.not. AC71_struc(n)%ac71(t)%sene_reached)) then ! from sowin
                    AC71_struc(n)%ac71(t)%Crop_DaysToSenescence = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%sene_reached = .true.
                end if
                if ((GetSimulation_SumGDDfromDay1() >= GetCrop_GDDaysToHarvest()) .and. (.not. AC71_struc(n)%ac71(t)%harv_reached)) then ! from sowing t
                    AC71_struc(n)%ac71(t)%Crop_DaysToHarvest = GetDayNri() - GetCrop_Day1()
                    AC71_struc(n)%ac71(t)%harv_reached = .true.
                end if
                end if
            end if
           
            ! MB: tmp, just for debugging
            sumGDD = GetSimulation_SumGDD()
            sumGDD2 = GetSumGDD()
           Crop_DaysToGermination = GetCrop_DaysToGermination()
           Crop_DaysToMaxRooting =         GetCrop_DaysToMaxRooting()
           Crop_DaysToFlowering =        GetCrop_DaysToFlowering()
           Crop_DaysToHarvest =        GetCrop_DaysToHarvest()
           Crop_DaysTosenescence =        GetCrop_DaysTosenescence()
           Crop_DaysToCCini =         GetCrop_DaysToCCini()
           Crop_DaysToFullCanopy =         GetCrop_DaysToFullCanopy()
           Crop_DaysToFullCanopySF =         GetCrop_DaysToFullCanopySF()
           Crop_DaysToHIo =         GetCrop_DaysToHIo()
           temp1 = GetCrop_GDDaysToGermination()
           temp1 = GetCrop_GDDaysToSenescence()

            ! Run AC
            tmp_wpi = REAL(AC71_struc(n)%ac71(t)%WPi,8)
            call AdvanceOneTimeStep(tmp_wpi, AC71_struc(n)%ac71(t)%HarvestNow)
            AC71_struc(n)%ac71(t)%WPi = tmp_wpi
            
            ! MB: tmp, just for debugging
            sumGDD = GetSimulation_SumGDD()
            sumGDD2 = GetSumGDD()

            ! Get all the ac71 variables and store in AC71_struc
            do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                    AC71_struc(n)%ac71(t)%smc(l) = GetCompartment_theta(l)
            enddo
            
            ! MB likely not needed since done with AC71_struc(n)%ac71(t)%crop = GetCrop() ?
            AC71_struc(n)%ac71(t)%Crop_DaysToGermination = GetCrop_DaysToGermination()
            AC71_struc(n)%ac71(t)%Crop_DaysToFlowering = GetCrop_DaysToFlowering()
            AC71_struc(n)%ac71(t)%Crop_DaysToMaxRooting = GetCrop_DaysToMaxRooting()
            AC71_struc(n)%ac71(t)%Crop_DaysToSenescence = GetCrop_DaysToSenescence()
            AC71_struc(n)%ac71(t)%Crop_DaysToHarvest = GetCrop_DaysToHarvest()
            AC71_struc(n)%ac71(t)%Crop_DaysToCCini = GetCrop_DaysToCCini()
            AC71_struc(n)%ac71(t)%Crop_DaysToFullCanopy = GetCrop_DaysToFullCanopy()
            AC71_struc(n)%ac71(t)%Crop_DaysToFullCanopySF = GetCrop_DaysToFullCanopySF()
            AC71_struc(n)%ac71(t)%Crop_DaysToHIo = GetCrop_DaysToHIo()
            !
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
            AC71_struc(n)%ac71(t)%crop = GetCrop()
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
            AC71_struc(n)%ac71(t)%Management = GetManagement()
            AC71_struc(n)%ac71(t)%PreviousStressLevel = GetPreviousStressLevel()
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

            ! MB: Apparently needed for GDD
            AC71_struc(n)%ac71(t)%PlotVarCrop = GetPlotVarCrop()

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

            !! If irrigation ON
            if(LIS_rc%irrigation_type.ne."none") then
                AC71_struc(n)%ac71(t)%IrriAfterSeason = GetIrriAfterSeason()
                AC71_struc(n)%ac71(t)%IrriBeforeSeason = GetIrriBeforeSeason()
                AC71_struc(n)%ac71(t)%IrriECw = GetIrriECw()
                AC71_struc(n)%ac71(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
                AC71_struc(n)%ac71(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
            endif

            ! Check for end of simulation period
            if ((LIS_rc%mo .eq. AC71_struc(n)%Sim_AnnualEndMonth) &
                .and.(LIS_rc%da .eq. AC71_struc(n)%Sim_AnnualEndDay)) then
                AC71_struc(n)%ac71(t)%InitializeRun = 1
                AC71_struc(n)%ac71(t)%irun = AC71_struc(n)%ac71(t)%irun + 1
            end if

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
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_RootingDepth, value = real(AC71_struc(n)%ac71(t)%Ziprev,kind=sp), &
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
