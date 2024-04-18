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
                    GenerateDepthMode_ToFC, &
                    GenerateTimeMode_AllRAW, &
                    GetCCiActual,&
                    GetCCiTopEarlySen,&
                    GetCCiprev,&
                    GetCompartment,&
                    GetCompartment_theta,&
                    GetCrop,&
                    GetCrop_Day1,&
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
                    GetSimulation,&
                    GetSimulation_EvapLimitON, &
                    GetSimulation_SumGDD,&
                    GetSimulation_SumGDDfromDay1,&
                    GetSimulation_SWCtopSoilConsidered, &
                    GetSoil,&
                    GetSoilLayer,&
                    GetSumWaBal,&
                    GetSurfaceStorage,&
                    GetTact,&
                    GetTactWeedInfested,&
                    GetTmax,& 
                    GetTmin,& 
                    GetTpot,&
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
                    SetTpot
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
                    GetSumInterval,&
                    GetSumKcTop,&
                    GetSumKcTopStress,&
                    GetSumKci,&
                    GetTadj,&
                    GetTheProjectFile,&
                    GetTimeSenescence ,&
                    GetTransfer,&
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
                    SetTransfer,&
                    SetWeedRCi,&
                    SetZiprev,&
                    SetalfaHI,&
                    SetalfaHIAdj               
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
            call SetCCoTotal(REAL(AC71_struc(n)%ac71(t)%CCoTotal, 8))  ! Set in InitializeRunPart2 --> copy INIR2nto setup for restart
            call SetCCxCropWeedsNoSFstress(REAL(AC71_struc(n)%ac71(t)%CCxCropWeedsNoSFstress, 8))   ! Set in InitializeRunPart2 --> copy INIR2nto setup for restart
            call SetCCxTotal(REAL(AC71_struc(n)%ac71(t)%CCxTotal, 8))   ! Set in InitializeRunPart2 --> copy INIR2nto setup for restart
            call SetCDCTotal(REAL(AC71_struc(n)%ac71(t)%CDCTotal, 8))     ! Set in InitializeRunPart2 --> copy INIR2nto setup for restart
            call SetCoeffb0(REAL(AC71_struc(n)%ac71(t)%Coeffb0, 8))      ! Set in InitializeRunPart1 --> copy INIR2nto setup for restart
            call SetCoeffb0Salt(REAL(AC71_struc(n)%ac71(t)%Coeffb0Salt, 8))       ! Set in InitializeRunPart1 --> copy INIR2nto setup for restart
            call SetCoeffb1(REAL(AC71_struc(n)%ac71(t)%Coeffb1, 8))       ! Set in InitializeRunPart1 --> copy INIR2nto setup for restart
            call SetCoeffb1Salt(REAL(AC71_struc(n)%ac71(t)%Coeffb1Salt, 8))       ! Set in InitializeRunPart1 --> copy INIR2nto setup for restart
            call SetCoeffb2(REAL(AC71_struc(n)%ac71(t)%Coeffb2, 8))       ! Set in InitializeRunPart1 --> copy INIR2nto setup for restart
            call SetCoeffb2Salt(REAL(AC71_struc(n)%ac71(t)%Coeffb2Salt, 8))       ! Set in InitializeRunPart1 --> copy INIR2nto setup for restart
            call SetCrop(AC71_struc(n)%ac71(t)%crop) ! CCxWithered, CCxAdjusted, pActStom passed in retart
            call SetGDDayFraction(REAL(AC71_struc(n)%ac71(t)%GDDayFraction, 8)) ! Set in InitializeRunPart2 --> copy INIR2nto setup for restart
            call SetGDDCDCTotal(REAL(AC71_struc(n)%ac71(t)%GDDCDCTotal, 8)) ! Set in InitializeRunPart1 --> copy INIR2nto setup for restart
            call SetGDDTadj(AC71_struc(n)%ac71(t)%GDDTadj)  ! Set in InitializeRunPart2 --> copy INIR2nto setup for restart
            call SetManagement(AC71_struc(n)%ac71(t)%Management) ! WeedDeltaRC passed in restart
            call SetSimulation(AC71_struc(n)%ac71(t)%Simulation) ! A few vars are passed in restart (EffectStress)
            call SetSoil(AC71_struc(n)%ac71(t)%Soil)  !Not needed for restart, everything should be set in setup
            call SetSoilLayer(AC71_struc(n)%ac71(t)%SoilLayer) !Not needed for restart, everything should be set in setup
            call SetSumKcTop(REAL(AC71_struc(n)%ac71(t)%SumKcTop, 8))  ! Set in InitializeRunPart1 --> copy INIR2nto setup for restart

            !! Fixed vars
            call SetCGCref(GetCrop_CGC()) ! Make sure crop is set before
            call SetGDDCGCref(GetCrop_GDDCGC()) ! Make sure crop is set before
            call SetNextSimFromDayNr(int(-9, kind=int32)) !Always undef_int in AquaCrop... Check src for v7.2
            call SetNoYear(.false.)
            call SetOutputAggregate(int(0,kind=int8)) ! Avoid writing out daily results in the console
            call SetPart1Mult(.false.) 
            call SetPart2Eval(.false.)
            call SetStartMode(.false.) ! Overwritten to .true. in InitalizeRunPart1
            call SetSumGDDPrev(GetSimulation_SumGDD()) ! Make sure that Simulation is set before

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

                call InitializeRunPart1(int(AC71_struc(n)%ac71(t)%irun,kind=int8), AC71_struc(n)%ac71(t)%TheProjectType)
                call InitializeSimulationRunPart2()
                AC71_struc(n)%ac71(t)%HarvestNow = .false.
                if(LIS_rc%irrigation_type.ne."none") then
                    call fIrri_close() !LB check if problem with specific irrigation file
                endif
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
            tmp_wpi = REAL(AC71_struc(n)%ac71(t)%WPi,8)
            call AdvanceOneTimeStep(tmp_wpi, AC71_struc(n)%ac71(t)%HarvestNow)
            AC71_struc(n)%ac71(t)%WPi = tmp_wpi

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
            AC71_struc(n)%ac71(t)%StressSenescence = GetStressSenescence()
            AC71_struc(n)%ac71(t)%StressSFadjNEW = GetStressSFadjNEW()
            AC71_struc(n)%ac71(t)%StressTot = GetStressTot()
            AC71_struc(n)%ac71(t)%SumGDDcuts = GetSumGDDcuts()
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
            AC71_struc(n)%ac71(t)%Transfer = GetTransfer()
            AC71_struc(n)%ac71(t)%WeedRCi = GetWeedRCi()
            AC71_struc(n)%ac71(t)%Ziprev = GetZiprev()

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
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71Yield, value = real(AC71_struc(n)%ac71(t)%SumWaBal%YieldPart,kind=sp), &
                                    vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71ZiPrev, value = real(AC71_struc(n)%ac71(t)%ZiPrev,kind=sp), &
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
