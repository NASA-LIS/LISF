!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
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
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!
!   9/4/14: Shugong Wang; initial implementation for Ac71 with LIS-7
!   2/7/18: Soni Yatheendradas; code added for OPTUE to work
!   18 JAN 2024, Louise Busschaert; initial implementation for LIS 7 and AC71
!
! !INTERFACE:
subroutine Ac71_main(n)
! !USES:
    use LIS_coreMod
    use LIS_histDataMod
    use LIS_timeMgrMod, only : LIS_isAlarmRinging
    use LIS_constantsMod,  only : LIS_CONST_RHOFW
    use LIS_logMod, only     : LIS_logunit, LIS_endrun
    use LIS_FORC_AttributesMod 
    use Ac71_lsmMod
   !use other modules
    use ESMF
    use LIS_routingMod, only : LIS_runoff_state
    !!! MB_AC71
    use ac_global, only: typeproject_typeprm, &
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
                        GetCrop_CCini,&
                        GetCrop_CCx,&
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
                        SetClimFile

           !!! MB_AC71
    !!! MB:
    use ac_project_input, only: ProjectInput 

    use ac_run, only:    SetDayNri,&
                         GetIrriInterval,&
                         GetIrriInfoRecord1,&
                         GetIrriInfoRecord2,&
                         SetIrriInterval,&
                         SetIrriInfoRecord1,&
                         SetIrriInfoRecord1_TimeInfo,&
                         SetIrriInfoRecord2,&
                         GetTheProjectFile,&
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
                         InitializeClimate
                         
    use ac_startunit, only:  FinalizeTheProgram, &
                         GetListProjectsFile, &
                         GetNumberOfProjects, &
                         GetProjectFileName, &
                         GetProjectType, &
                         GetSimulation_NrRuns, &
                         InitializeTheProgram, &
                         InitializeProject, &
                         WriteProjectsInfo


    use ac_kinds, only: intEnum, &
                        int32, &
                        int8, &
                        dp,&
                        sp
    !!! MB_AC71

    implicit none

    !!! MB_AC71
    integer :: daynr, todaynr, iproject, nprojects
    logical :: ListProjectFileExist
    character(len=:), allocatable :: ListProjectsFile, TheProjectFile

    !!! MB_AC71

! !ARGUMENTS:
    integer, intent(in)  :: n
    integer              :: t
    integer              :: i
    integer              :: itemp, countertemp
    real                 :: dt
    real                 :: lat, lon
    real                 :: Tmin_movmean
    real                 :: Tmin_mplr
    integer              :: row, col
    integer              :: year, month, day, hour, minute, second
    logical              :: alarmCheck

    integer               :: status
    integer               :: c,r,l
    integer               :: ios, nid,rivid,fldid

    integer            :: tid
!
! !DESCRIPTION:
!  This is the entry point for calling the Ac71 physics.
!  This routine calls the {\tt ac\_driver\_71 } routine that performs the
!  land surface computations, to solve for water and energy equations.

!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!EOP

! define variables for Ac71
    
    !MB: AC71
    real                 :: tmp_PREC_ac        ! 
    real                 :: tmp_TMIN_ac        ! 
    real                 :: tmp_TMAX_ac        ! 
    real                 :: tmp_ETo_ac         ! 

    ! check Ac71 alarm. If alarm is ring, run model. 
    alarmCheck = LIS_isAlarmRinging(LIS_rc, "Ac71 model alarm")
    if (alarmCheck) Then

       do t = 1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          dt = LIS_rc%ts
          row = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%row
          col = LIS_surface(n, LIS_rc%lsm_index)%tile(t)%col
          lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
          lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

          ! PREC_ac
          tmp_PREC_ac      = AC71_struc(n)%ac71(t)%PREC_ac  / AC71_struc(n)%forc_count
          ! TMIN_ac
          tmp_TMIN_ac      = AC71_struc(n)%ac71(t)%TMIN_ac  / AC71_struc(n)%forc_count
          ! TMAX_ac
          tmp_TMAX_ac      = AC71_struc(n)%ac71(t)%TMAX_ac  / AC71_struc(n)%forc_count
          ! ETo_ac
          tmp_ETo_ac      = AC71_struc(n)%ac71(t)%ETo_ac  / AC71_struc(n)%forc_count

          ! check validity of PREC_ac
          if(tmp_PREC_ac .eq. LIS_rc%udef) then
              write(LIS_logunit, *) "undefined value found for forcing variable PREC_ac in Ac71"
              write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
              call LIS_endrun()
          endif
            
          ! check validity of TMIN
          if(tmp_TMIN_ac .eq. LIS_rc%udef) then
              write(LIS_logunit, *) "undefined value found for forcing variable TMIN in Ac71"
              write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
              call LIS_endrun()
          endif
          
          ! check validity of TMAX
          if(tmp_TMAX_ac .eq. LIS_rc%udef) then
              write(LIS_logunit, *) "undefined value found for forcing variable TMAX in Ac71"
              write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
              call LIS_endrun()
          endif
            
          ! check validity of ETo
          if(tmp_ETo_ac .eq. LIS_rc%udef) then
              write(LIS_logunit, *) "undefined value found for forcing variable ETo in Ac71"
              write(LIS_logunit, *) "for tile ", t, "latitude = ", lat, "longitude = ", lon
              call LIS_endrun()
          endif
            

            !!! MB_AC71

            ! setting all global variables
            call SetRootZoneWC_Actual(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_Actual,8))
            call SetRootZoneWC_FC(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_FC,8))
            call SetRootZoneWC_WP(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_WP,8))
            call SetRootZoneWC_SAT(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_SAT,8))
            call SetRootZoneWC_Leaf(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_Leaf,8))
            call SetRootZoneWC_Thresh(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_Sen(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopAct(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopFC(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopWP(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetRootZoneWC_ZtopThresh(REAL(AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct,8))
            call SetCompartment(AC71_struc(n)%ac71(t)%Compartment)
            call SetTotalSaltContent(AC71_struc(n)%ac71(t)%TotalSaltContent)
            call SetTotalWaterContent(AC71_struc(n)%ac71(t)%TotalWaterContent)
            call Seteffectiverain(AC71_struc(n)%ac71(t)%effectiverain)
            call SetSumWaBal(AC71_struc(n)%ac71(t)%SumWaBal)
            call SetRootZoneSalt(AC71_struc(n)%ac71(t)%RootZoneSalt)
            call SetSimulation(AC71_struc(n)%ac71(t)%Simulation)
            call SetIrriInterval(AC71_struc(n)%ac71(t)%IrriInterval)
            call SetIrriInfoRecord1(AC71_struc(n)%ac71(t)%IrriInfoRecord1)
            call SetIrriInfoRecord2(AC71_struc(n)%ac71(t)%IrriInfoRecord2)
            call SetIrrigation(REAL(AC71_struc(n)%ac71(t)%Irrigation,8))
            do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                 call SetCompartment_theta(l,REAL(AC71_struc(n)%ac71(t)%smc(l),8))
            enddo
            call SetIrriECw(AC71_struc(n)%ac71(t)%IrriECw) 
            call SetManagement(AC71_struc(n)%ac71(t)%Management) 
            call SetPerennialPeriod(AC71_struc(n)%ac71(t)%PerennialPeriod) 
            call SetSimulParam(AC71_struc(n)%ac71(t)%simulparam) 
            call SetManagement_Cuttings(AC71_struc(n)%ac71(t)%Cuttings) 
            call SetOnset(AC71_struc(n)%ac71(t)%onset) 
            call SetEndSeason(AC71_struc(n)%ac71(t)%endseason) 
            call SetCrop(AC71_struc(n)%ac71(t)%crop) 
            call SetSoil(AC71_struc(n)%ac71(t)%Soil) 
            call SetIrriBeforeSeason(AC71_struc(n)%ac71(t)%IrriBeforeSeason)
            call SetIrriAfterSeason(AC71_struc(n)%ac71(t)%IrriAfterSeason)
            call SetSoilLayer(AC71_struc(n)%ac71(t)%soillayer)
            call SetDayNri(AC71_struc(n)%ac71(t)%daynri)

            call SetGenerateTimeMode(AC71_struc(n)%ac71(t)%GenerateTimeMode) 
            call SetGenerateDepthMode(AC71_struc(n)%ac71(t)%GenerateDepthMode) 
            call SetIrriMode(AC71_struc(n)%ac71(t)%IrriMode) 
            call SetDaySubmerged(AC71_struc(n)%ac71(t)%DaySubmerged) 
            call SetMaxPlotNew(AC71_struc(n)%ac71(t)%MaxPlotNew) 
            call SetNrCompartments(AC71_struc(n)%ac71(t)%NrCompartments) 
            call SetIrriFirstDayNr(AC71_struc(n)%ac71(t)%IrriFirstDayNr) 
            call SetZiAqua(AC71_struc(n)%ac71(t)%ZiAqua) 
            call SetIniPercTAW(AC71_struc(n)%ac71(t)%IniPercTAW) 
            call SetMaxPlotTr(AC71_struc(n)%ac71(t)%MaxPlotTr)
            call SetOutputAggregate(AC71_struc(n)%ac71(t)%OutputAggregate) 

            call SetEvapoEntireSoilSurface(AC71_struc(n)%ac71(t)%EvapoEntireSoilSurface) 
            call SetPreDay(AC71_struc(n)%ac71(t)%PreDay) 
            call SetOutDaily(AC71_struc(n)%ac71(t)%OutDaily) 
            call SetOut1Wabal(AC71_struc(n)%ac71(t)%Out1Wabal) 
            call SetOut2Crop(AC71_struc(n)%ac71(t)%Out2Crop) 
            call SetOut3Prof(AC71_struc(n)%ac71(t)%Out3Prof) 
            call SetOut4Salt(AC71_struc(n)%ac71(t)%Out4Salt) 
            call SetOut5CompWC(AC71_struc(n)%ac71(t)%Out5CompWC) 
            call SetOut6CompEC(AC71_struc(n)%ac71(t)%Out6CompEC) 
            call SetOut7Clim(AC71_struc(n)%ac71(t)%Out7Clim) 
            call SetPart1Mult(AC71_struc(n)%ac71(t)%Part1Mult) 
            call SetPart2Eval(AC71_struc(n)%ac71(t)%Part2Eval) 

            !
            call SetCCiActual(AC71_struc(n)%ac71(t)%CCiActual)
            call SetCCiprev(AC71_struc(n)%ac71(t)%CCiprev)

            call SetCCiTopEarlySen(AC71_struc(n)%ac71(t)%CCiTopEarlySen)
            call SetCRsalt(AC71_struc(n)%ac71(t)%CRsalt)
            call SetCRwater(AC71_struc(n)%ac71(t)%CRwater)
            call SetECdrain(AC71_struc(n)%ac71(t)%ECdrain)
            call SetEciAqua(AC71_struc(n)%ac71(t)%ECiAqua)
            call SetECstorage(AC71_struc(n)%ac71(t)%ECstorage)
            call SetEact(AC71_struc(n)%ac71(t)%Eact)
            call SetEpot(AC71_struc(n)%ac71(t)%Epot)
            call SetDrain(AC71_struc(n)%ac71(t)%Drain)
            call SetInfiltrated(AC71_struc(n)%ac71(t)%Infiltrated)
            call SetRootingDepth(AC71_struc(n)%ac71(t)%RootingDepth)
            call SetRunoff(AC71_struc(n)%ac71(t)%Runoff)
            call SetSaltInfiltr(AC71_struc(n)%ac71(t)%SaltInfiltr)
            call SetSurf0(AC71_struc(n)%ac71(t)%Surf0)
            call SetSurfaceStorage(AC71_struc(n)%ac71(t)%SurfaceStorage)
            call SetTact(AC71_struc(n)%ac71(t)%Tact)
            call SetTpot(AC71_struc(n)%ac71(t)%Tpot)
            call SetTactWeedInfested(AC71_struc(n)%ac71(t)%TactWeedInfested)

            call SetGwTable(AC71_struc(n)%ac71(t)%GwTable)
            call SetPlotVarCrop(AC71_struc(n)%ac71(t)%PlotVarCrop)
            call SetStressTot(AC71_struc(n)%ac71(t)%StressTot)
            call SetCutInfoRecord1(AC71_struc(n)%ac71(t)%CutInfoRecord1)
            call SetCutInfoRecord2(AC71_struc(n)%ac71(t)%CutInfoRecord2)
            call SetTransfer(AC71_struc(n)%ac71(t)%Transfer)
            call SetPreviousSum(AC71_struc(n)%ac71(t)%PreviousSum)
            call SetTadj(AC71_struc(n)%ac71(t)%Tadj)
            call SetGDDTadj(AC71_struc(n)%ac71(t)%GDDTadj)
            call SetDayLastCut(AC71_struc(n)%ac71(t)%DayLastCut)
            call SetNrCut(AC71_struc(n)%ac71(t)%NrCut)
            call SetSumInterval(AC71_struc(n)%ac71(t)%SumInterval)
            call SetPreviousStressLevel(int(AC71_struc(n)%ac71(t)%PreviousStressLevel,kind=int32))
            call SetStressSFadjNEW(int(AC71_struc(n)%ac71(t)%StressSFadjNEW,kind=int32))
            call SetBin(AC71_struc(n)%ac71(t)%Bin)
            call SetBout(AC71_struc(n)%ac71(t)%Bout)
            call SetCO2i(AC71_struc(n)%ac71(t)%CO2i)
            call SetFracBiomassPotSF(AC71_struc(n)%ac71(t)%FracBiomassPotSF)
            call SetSumETo(AC71_struc(n)%ac71(t)%SumETo)
            call SetSumGDD(AC71_struc(n)%ac71(t)%SumGDD)
            call SetZiprev(AC71_struc(n)%ac71(t)%Ziprev)
            call SetSumGDDPrev(AC71_struc(n)%ac71(t)%SumGDDPrev)
            call SetCCxWitheredTpotNoS(AC71_struc(n)%ac71(t)%CCxWitheredTpotNoS)
            call SetCoeffb0(AC71_struc(n)%ac71(t)%Coeffb0)
            call SetCoeffb1(AC71_struc(n)%ac71(t)%Coeffb1)
            call SetCoeffb2(AC71_struc(n)%ac71(t)%Coeffb2)
            call SetCoeffb0Salt(AC71_struc(n)%ac71(t)%Coeffb0Salt)
            call SetCoeffb1Salt(AC71_struc(n)%ac71(t)%Coeffb1Salt)
            call SetCoeffb2Salt(AC71_struc(n)%ac71(t)%Coeffb2Salt)
            call SetStressLeaf(AC71_struc(n)%ac71(t)%StressLeaf)
            call SetStressSenescence(AC71_struc(n)%ac71(t)%StressSenescence)
            call SetDayFraction(AC71_struc(n)%ac71(t)%DayFraction)
            call SetGDDayFraction(AC71_struc(n)%ac71(t)%GDDayFraction)
            call SetCGCref(AC71_struc(n)%ac71(t)%CGCref)
            call SetGDDCGCref(AC71_struc(n)%ac71(t)%GDDCGCref)
            call SetTimeSenescence(AC71_struc(n)%ac71(t)%TimeSenescence)
            call SetSumKcTop(AC71_struc(n)%ac71(t)%SumKcTop)
            call SetSumKcTopStress(AC71_struc(n)%ac71(t)%SumKcTopStress)
            call SetSumKci(AC71_struc(n)%ac71(t)%SumKci)
            call SetCCoTotal(AC71_struc(n)%ac71(t)%CCoTotal)
            call SetCCxTotal(AC71_struc(n)%ac71(t)%CCxTotal)
            call SetCDCTotal(AC71_struc(n)%ac71(t)%CDCTotal)
            call SetGDDCDCTotal(AC71_struc(n)%ac71(t)%GDDCDCTotal)
            call SetCCxCropWeedsNoSFstress(AC71_struc(n)%ac71(t)%CCxCropWeedsNoSFstress)
            call SetWeedRCi(AC71_struc(n)%ac71(t)%WeedRCi)
            call SetCCiActualWeedInfested(AC71_struc(n)%ac71(t)%CCiActualWeedInfested)
            call SetfWeedNoS(AC71_struc(n)%ac71(t)%fWeedNoS)
            call SetZeval(AC71_struc(n)%ac71(t)%Zeval)
            call SetBprevSum(AC71_struc(n)%ac71(t)%BprevSum)
            call SetYprevSum(AC71_struc(n)%ac71(t)%YprevSum)
            call SetSumGDDcuts(AC71_struc(n)%ac71(t)%SumGDDcuts)
            call SetHItimesBEF(AC71_struc(n)%ac71(t)%HItimesBEF)
            call SetScorAT1(AC71_struc(n)%ac71(t)%ScorAT1)
            call SetScorAT2(AC71_struc(n)%ac71(t)%ScorAT2)
            call SetHItimesAT1(AC71_struc(n)%ac71(t)%HItimesAT1)
            call SetHItimesAT2(AC71_struc(n)%ac71(t)%HItimesAT2)
            call SetHItimesAT(AC71_struc(n)%ac71(t)%HItimesAT)
            call SetalfaHI(AC71_struc(n)%ac71(t)%alfaHI)
            call SetalfaHIAdj(AC71_struc(n)%ac71(t)%alfaHIAdj)
            call SetNextSimFromDayNr(AC71_struc(n)%ac71(t)%NextSimFromDayNr)
            call SetDayNr1Eval(AC71_struc(n)%ac71(t)%DayNr1Eval)
            call SetDayNrEval(AC71_struc(n)%ac71(t)%DayNrEval)
            call SetLineNrEval(int(AC71_struc(n)%ac71(t)%LineNrEval,kind=int32))
            call SetPreviousSumETo(AC71_struc(n)%ac71(t)%PreviousSumETo)
            call SetPreviousSumGDD(AC71_struc(n)%ac71(t)%PreviousSumGDD)
            call SetPreviousBmob(AC71_struc(n)%ac71(t)%PreviousBmob)
            call SetPreviousBsto(AC71_struc(n)%ac71(t)%PreviousBsto)
            call SetStageCode(AC71_struc(n)%ac71(t)%StageCode)
            call SetPreviousDayNr(AC71_struc(n)%ac71(t)%PreviousDayNr)
            call SetNoYear(AC71_struc(n)%ac71(t)%NoYear)
            call SetWaterTableInProfile(AC71_struc(n)%ac71(t)%WaterTableInProfile)
            call SetStartMode(AC71_struc(n)%ac71(t)%StartMode)
            call SetNoMoreCrop(AC71_struc(n)%ac71(t)%NoMoreCrop)
            !call SetSimulation_ToDayNr(AC71_struc(n)%ac71(t)%Simulation%ToDayNr)
            call SetGDDayi(AC71_struc(n)%ac71(t)%GDDayi)

            !!! initialize run (year)

             if (AC71_struc(n)%ac71(t)%InitializeRun .eq. 1) then
                ! Replaces LoadSimulationRunProject in LoadSimulationRunProject
                !call SetSimulation_YearSeason(YearSeason_temp)
                !call SetCrop_Day1(TempInt)
                !call SetCrop_DayN(TempInt)
                !SetCO2FileFull

                ! Run InitializeRunPart1 ( in future without LoadSimulationRunProject)
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
                        !call SetClimFile('EToRainTempFile')
                        !call SetClimDescription('Read ETo/RAIN/TEMP data set')

                AC71_struc(n)%ac71(t)%WPi = 0._dp

                call InitializeRunPart1(AC71_struc(n)%ac71(t)%irun, AC71_struc(n)%ac71(t)%TheProjectType);
                if (trim(LIS_rc%metforc(1)) == 'MERRA2_AC') then
                  call SetRain(real(tmp_PREC_ac,kind=dp))
                  call SetTmin(real(tmp_TMIN_ac,kind=dp))
                  call SetTmax(real(tmp_TMAX_ac,kind=dp))
                  call SetETo(real(tmp_ETo_ac,kind=dp))
                else ! read from AC input
                  call InitializeClimate();
                end if

                !call InitializeRunPart2(AC71_struc(n)%ac71(t)%irun, AC71_struc(n)%ac71(t)%TheProjectType);
                call InitializeSimulationRunPart2()
                !LB new for AC7.1
                AC71_struc(n)%ac71(t)%HarvestNow = .false.
                AC71_struc(n)%ac71(t)%InitializeRun = 0
            end if

            ! for MERRA2_AC --> first set climate variables then advanceonetimestep
            call SetRain(real(tmp_PREC_ac,kind=dp))
            call SetTmin(real(tmp_TMIN_ac,kind=dp))
            call SetTmax(real(tmp_TMAX_ac,kind=dp))
            call SetETo(real(tmp_ETo_ac,kind=dp))
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

            !write(LIS_logunit,*) &
            !                '[INFO] AdvanceOneTimeStep AquaCrop, day in month: ', LIS_rc
            call AdvanceOneTimeStep(AC71_struc(n)%ac71(t)%WPi, AC71_struc(n)%ac71(t)%HarvestNow)

            AC71_struc(n)%ac71(t)%RootZoneWC_Actual = GetRootZoneWC_Actual()
            AC71_struc(n)%ac71(t)%RootZoneWC_FC = GetRootZoneWC_FC()
            AC71_struc(n)%ac71(t)%RootZoneWC_WP = GetRootZoneWC_WP()
            AC71_struc(n)%ac71(t)%RootZoneWC_SAT = GetRootZoneWC_SAT()
            AC71_struc(n)%ac71(t)%RootZoneWC_Leaf = GetRootZoneWC_Leaf()
            AC71_struc(n)%ac71(t)%RootZoneWC_Thresh = GetRootZoneWC_Thresh()
            AC71_struc(n)%ac71(t)%RootZoneWC_Sen = GetRootZoneWC_Sen()
            AC71_struc(n)%ac71(t)%RootZoneWC_ZtopAct = GetRootZoneWC_ZtopAct()
            AC71_struc(n)%ac71(t)%RootZoneWC_ZtopFC = GetRootZoneWC_ZtopFC()
            AC71_struc(n)%ac71(t)%RootZoneWC_ZtopWP = GetRootZoneWC_ZtopWP()
            AC71_struc(n)%ac71(t)%RootZoneWC_ZtopThresh = GetRootZoneWC_ZtopThresh()
            AC71_struc(n)%ac71(t)%Compartment = GetCompartment()
            AC71_struc(n)%ac71(t)%TotalSaltContent = GetTotalSaltContent()
            AC71_struc(n)%ac71(t)%TotalWaterContent = GetTotalWaterContent()
            AC71_struc(n)%ac71(t)%effectiverain = Geteffectiverain()
            AC71_struc(n)%ac71(t)%SumWaBal = GetSumWaBal()
            AC71_struc(n)%ac71(t)%RootZoneSalt = GetRootZoneSalt()
            AC71_struc(n)%ac71(t)%Simulation = GetSimulation()
            AC71_struc(n)%ac71(t)%IrriInterval = GetIrriInterval()
            AC71_struc(n)%ac71(t)%IrriInfoRecord1 = GetIrriInfoRecord1()
            AC71_struc(n)%ac71(t)%IrriInfoRecord2 = GetIrriInfoRecord2()
            AC71_struc(n)%ac71(t)%Irrigation = GetIrrigation()
            AC71_struc(n)%ac71(t)%IrriBeforeSeason = GetIrriBeforeSeason()
            AC71_struc(n)%ac71(t)%IrriAfterSeason = GetIrriAfterSeason()
            AC71_struc(n)%ac71(t)%SoilLayer = GetSoilLayer()
            AC71_struc(n)%ac71(t)%daynri = GetDayNri()
            do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                 AC71_struc(n)%ac71(t)%smc(l) = GetCompartment_theta(l)
            enddo
            !write(*,'(e23.15e3)') AC71_struc(n)%ac71(t)%ac71smc(1)
            AC71_struc(n)%ac71(t)%IrriECw = GetIrriECw()
            AC71_struc(n)%ac71(t)%Management = GetManagement()
            AC71_struc(n)%ac71(t)%PerennialPeriod = GetPerennialPeriod()
            AC71_struc(n)%ac71(t)%simulparam = GetSimulParam()
            AC71_struc(n)%ac71(t)%Cuttings = GetManagement_Cuttings()
            AC71_struc(n)%ac71(t)%onset = GetOnset()
            AC71_struc(n)%ac71(t)%endseason = GetEndSeason()
            AC71_struc(n)%ac71(t)%crop = GetCrop()
            AC71_struc(n)%ac71(t)%Soil = GetSoil()
            AC71_struc(n)%ac71(t)%TemperatureRecord = GetTemperatureRecord()
            AC71_struc(n)%ac71(t)%ClimRecord = GetClimRecord()
            AC71_struc(n)%ac71(t)%RainRecord = GetRainRecord()
            AC71_struc(n)%ac71(t)%EToRecord = GetEToRecord()

            AC71_struc(n)%ac71(t)%GenerateTimeMode = GetGenerateTimeMode()
            AC71_struc(n)%ac71(t)%GenerateDepthMode = GetGenerateDepthMode()
            AC71_struc(n)%ac71(t)%IrriMode = GetIrriMode()
            AC71_struc(n)%ac71(t)%IrriMethod = GetIrriMethod()
            AC71_struc(n)%ac71(t)%DaySubmerged = GetDaySubmerged()
            AC71_struc(n)%ac71(t)%MaxPlotNew = GetMaxPlotNew()
            AC71_struc(n)%ac71(t)%NrCompartments = GetNrCompartments()
            AC71_struc(n)%ac71(t)%IrriFirstDayNr = GetIrriFirstDayNr()
            AC71_struc(n)%ac71(t)%ZiAqua = GetZiAqua()
            AC71_struc(n)%ac71(t)%IniPercTAW = GetIniPercTAW()
            AC71_struc(n)%ac71(t)%MaxPlotTr = GetMaxPlotTr()
            AC71_struc(n)%ac71(t)%OutputAggregate = GetOutputAggregate()

            AC71_struc(n)%ac71(t)%EvapoEntireSoilSurface = GetEvapoEntireSoilSurface()
            AC71_struc(n)%ac71(t)%PreDay = GetPreDay()
            AC71_struc(n)%ac71(t)%OutDaily = GetOutDaily()
            AC71_struc(n)%ac71(t)%Out1Wabal = GetOut1Wabal()
            AC71_struc(n)%ac71(t)%Out2Crop = GetOut2Crop()
            AC71_struc(n)%ac71(t)%Out3Prof = GetOut3Prof()
            AC71_struc(n)%ac71(t)%Out4Salt = GetOut4Salt()
            AC71_struc(n)%ac71(t)%Out5CompWC = GetOut5CompWC()
            AC71_struc(n)%ac71(t)%Out6CompEC = GetOut6CompEC()
            AC71_struc(n)%ac71(t)%Out7Clim = GetOut7Clim()
            AC71_struc(n)%ac71(t)%Part1Mult = GetPart1Mult()
            AC71_struc(n)%ac71(t)%Part2Eval = GetPart2Eval()

            !
            AC71_struc(n)%ac71(t)%CCiActual = GetCCiActual()
            AC71_struc(n)%ac71(t)%CCiprev = GetCCiprev()
            AC71_struc(n)%ac71(t)%CCiTopEarlySen = GetCCiTopEarlySen()
            AC71_struc(n)%ac71(t)%CRsalt = GetCRsalt () ! gram/m2
            AC71_struc(n)%ac71(t)%CRwater = GetCRwater() ! mm/day
            AC71_struc(n)%ac71(t)%ECdrain = GetECdrain() ! EC drain water dS/m
            AC71_struc(n)%ac71(t)%ECiAqua = GetECiAqua() ! EC of the groundwater table in dS/m
            AC71_struc(n)%ac71(t)%ECstorage = GetECstorage() !EC surface storage dS/m
            AC71_struc(n)%ac71(t)%Eact = GetEact() ! mm/day
            AC71_struc(n)%ac71(t)%Epot = GetEpot() ! mm/day
            AC71_struc(n)%ac71(t)%ETo_ac = GetETo() ! mm/day
            AC71_struc(n)%ac71(t)%Drain = GetDrain()  ! mm/day
            AC71_struc(n)%ac71(t)%Infiltrated = GetInfiltrated() ! mm/day
            AC71_struc(n)%ac71(t)%PREC_ac = GetRain()  ! mm/day
            AC71_struc(n)%ac71(t)%RootingDepth = GetRootingDepth()
            AC71_struc(n)%ac71(t)%Runoff = GetRunoff()  ! mm/day
            AC71_struc(n)%ac71(t)%SaltInfiltr = GetSaltInfiltr() ! salt infiltrated in soil profile Mg/ha
            AC71_struc(n)%ac71(t)%Surf0 = GetSurf0()  ! surface water [mm] begin day
            AC71_struc(n)%ac71(t)%SurfaceStorage = GetSurfaceStorage() !mm/day
            AC71_struc(n)%ac71(t)%Tact = GetTact() ! mm/day
            AC71_struc(n)%ac71(t)%Tpot = GetTpot() ! mm/day
            AC71_struc(n)%ac71(t)%TactWeedInfested = GetTactWeedInfested() !mm/day
            AC71_struc(n)%ac71(t)%Tmax_ac = GetTmax() ! degC
            AC71_struc(n)%ac71(t)%Tmin_ac =GetTmin() ! degC


            AC71_struc(n)%ac71(t)%GwTable = GetGwTable()
            AC71_struc(n)%ac71(t)%PlotVarCrop = GetPlotVarCrop()
            AC71_struc(n)%ac71(t)%StressTot = GetStressTot()
            AC71_struc(n)%ac71(t)%CutInfoRecord1 = GetCutInfoRecord1()
            AC71_struc(n)%ac71(t)%CutInfoRecord2 = GetCutInfoRecord2()
            AC71_struc(n)%ac71(t)%Transfer = GetTransfer()
            AC71_struc(n)%ac71(t)%PreviousSum = GetPreviousSum()
            AC71_struc(n)%ac71(t)%Tadj = GetTadj()
            AC71_struc(n)%ac71(t)%GDDTadj = GetGDDTadj()
            AC71_struc(n)%ac71(t)%DayLastCut = GetDayLastCut()
            AC71_struc(n)%ac71(t)%NrCut = GetNrCut()
            AC71_struc(n)%ac71(t)%SumInterval = GetSumInterval()
            AC71_struc(n)%ac71(t)%PreviousStressLevel = GetPreviousStressLevel()
            AC71_struc(n)%ac71(t)%StressSFadjNEW = GetStressSFadjNEW()
            AC71_struc(n)%ac71(t)%Bin = GetBin()
            AC71_struc(n)%ac71(t)%Bout = GetBout()
            AC71_struc(n)%ac71(t)%GDDayi = GetGDDayi()
            AC71_struc(n)%ac71(t)%CO2i = GetCO2i()
            AC71_struc(n)%ac71(t)%FracBiomassPotSF = GetFracBiomassPotSF()
            AC71_struc(n)%ac71(t)%SumETo = GetSumETo()
            AC71_struc(n)%ac71(t)%SumGDD = GetSumGDD()
            AC71_struc(n)%ac71(t)%Ziprev = GetZiprev()
            AC71_struc(n)%ac71(t)%SumGDDPrev = GetSumGDDPrev()
            AC71_struc(n)%ac71(t)%CCxWitheredTpotNoS = GetCCxWitheredTpotNoS()
            AC71_struc(n)%ac71(t)%Coeffb0 = GetCoeffb0()
            AC71_struc(n)%ac71(t)%Coeffb1 = GetCoeffb1()
            AC71_struc(n)%ac71(t)%Coeffb2 = GetCoeffb2()
            AC71_struc(n)%ac71(t)%Coeffb0Salt = GetCoeffb0Salt()
            AC71_struc(n)%ac71(t)%Coeffb1Salt = GetCoeffb1Salt()
            AC71_struc(n)%ac71(t)%Coeffb2Salt = GetCoeffb2Salt()
            AC71_struc(n)%ac71(t)%StressLeaf = GetStressLeaf()
            AC71_struc(n)%ac71(t)%StressSenescence = GetStressSenescence ()
            AC71_struc(n)%ac71(t)%DayFraction = GetDayFraction()
            AC71_struc(n)%ac71(t)%GDDayFraction = GetGDDayFraction()
            AC71_struc(n)%ac71(t)%CGCref = GetCGCref()
            AC71_struc(n)%ac71(t)%GDDCGCref = GetGDDCGCref ()
            AC71_struc(n)%ac71(t)%TimeSenescence = GetTimeSenescence ()
            AC71_struc(n)%ac71(t)%SumKcTop = GetSumKcTop()
            AC71_struc(n)%ac71(t)%SumKcTopStress = GetSumKcTopStress()
            AC71_struc(n)%ac71(t)%SumKci = GetSumKci()
            AC71_struc(n)%ac71(t)%CCoTotal = GetCCoTotal()
            AC71_struc(n)%ac71(t)%CCxTotal = GetCCxTotal()
            AC71_struc(n)%ac71(t)%CDCTotal = GetCDCTotal()
            AC71_struc(n)%ac71(t)%GDDCDCTotal = GetGDDCDCTotal()
            AC71_struc(n)%ac71(t)%CCxCropWeedsNoSFstress = GetCCxCropWeedsNoSFstress()
            AC71_struc(n)%ac71(t)%WeedRCi = GetWeedRCi()
            AC71_struc(n)%ac71(t)%CCiActualWeedInfested = GetCCiActualWeedInfested()
            AC71_struc(n)%ac71(t)%fWeedNoS = GetfWeedNoS()
            AC71_struc(n)%ac71(t)%Zeval = GetZeval()
            AC71_struc(n)%ac71(t)%BprevSum = GetBprevSum()
            AC71_struc(n)%ac71(t)%YprevSum = GetYprevSum()
            AC71_struc(n)%ac71(t)%SumGDDcuts = GetSumGDDcuts()
            AC71_struc(n)%ac71(t)%HItimesBEF = GetHItimesBEF()
            AC71_struc(n)%ac71(t)%ScorAT1 = GetScorAT1()
            AC71_struc(n)%ac71(t)%ScorAT2 = GetScorAT2()
            AC71_struc(n)%ac71(t)%HItimesAT1 = GetHItimesAT1()
            AC71_struc(n)%ac71(t)%HItimesAT2 = GetHItimesAT2()
            AC71_struc(n)%ac71(t)%HItimesAT = GetHItimesAT()
            AC71_struc(n)%ac71(t)%alfaHI = GetalfaHI()
            AC71_struc(n)%ac71(t)%alfaHIAdj = GetalfaHIAdj()
            AC71_struc(n)%ac71(t)%NextSimFromDayNr = GetNextSimFromDayNr ()
            AC71_struc(n)%ac71(t)%DayNr1Eval = GetDayNr1Eval()
            AC71_struc(n)%ac71(t)%DayNrEval = GetDayNrEval()
            AC71_struc(n)%ac71(t)%LineNrEval = GetLineNrEval()
            AC71_struc(n)%ac71(t)%PreviousSumETo = GetPreviousSumETo()
            AC71_struc(n)%ac71(t)%PreviousSumGDD = GetPreviousSumGDD()
            AC71_struc(n)%ac71(t)%PreviousBmob = GetPreviousBmob()
            AC71_struc(n)%ac71(t)%PreviousBsto = GetPreviousBsto()
            AC71_struc(n)%ac71(t)%StageCode = GetStageCode()
            AC71_struc(n)%ac71(t)%PreviousDayNr = GetPreviousDayNr()
            AC71_struc(n)%ac71(t)%NoYear = GetNoYear()
            AC71_struc(n)%ac71(t)%WaterTableInProfile = GetWaterTableInProfile()
            AC71_struc(n)%ac71(t)%StartMode = GetStartMode()
            AC71_struc(n)%ac71(t)%NoMoreCrop = GetNoMoreCrop()
          

            if ((LIS_rc%mo .eq. 12) .AND. (LIS_rc%da .eq. 31)) then
                AC71_struc(n)%ac71(t)%InitializeRun = 1
                !call FinalizeRun1(AC71_struc(n)%ac71(t)%irun, GetTheProjectFile(), AC71_struc(n)%ac71(t)%TheProjectType)
                !call FinalizeRun2(AC71_struc(n)%ac71(t)%irun, AC71_struc(n)%ac71(t)%TheProjectType)
                AC71_struc(n)%ac71(t)%irun = AC71_struc(n)%ac71(t)%irun + 1
            end if
            !!! MB_AC71


            ! MB: AC71
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
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71RootingDepth, value = real(AC71_struc(n)%ac71(t)%RootingDepth,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)
            call LIS_diagnoseSurfaceOutputVar(n, t, LIS_MOC_AC71CCiActual, value = real(AC71_struc(n)%ac71(t)%CCiActual,kind=sp), &
                                              vlevel=1, unit="-", direction="-", surface_type = LIS_rc%lsm_index)

            AC71_struc(n)%ac71(t)%PREC_ac = 0.0
            AC71_struc(n)%ac71(t)%TMIN_ac = 0.0
            AC71_struc(n)%ac71(t)%TMAX_ac = 0.0
            AC71_struc(n)%ac71(t)%ETo_ac = 0.0
        enddo ! end of tile (t) loop
        ! reset forcing counter to be zero
        AC71_struc(n)%forc_count = 0 
    endif ! end of alarmCheck loop 

end subroutine Ac71_main
