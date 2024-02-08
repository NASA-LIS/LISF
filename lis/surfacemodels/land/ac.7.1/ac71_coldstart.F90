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
! !ROUTINE: Ac71_coldstart
! \label{Ac71_coldstart}
!
! !REVISION HISTORY:
!  This subroutine is generated with the Model Implementation Toolkit developed
!  by Shugong Wang for the NASA Land Information System Version 7. The initial 
!  specification of the subroutine is defined by Sujay Kumar. 
!   18 JAN 2024, Louise Busschaert; initial implementation for LIS 7 and AC71
!
! !INTERFACE:
subroutine Ac71_coldstart(mtype)
! !USES:
    use LIS_coreMod, only: LIS_rc
    use LIS_logMod, only: LIS_logunit
    use LIS_timeMgrMod, only: LIS_date2time
    use Ac71_lsmMod
 
!!! MB: AC71
    use ac_global, only: GetSimulParam_ThicknessTopSWC,&
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
                         GetTotalSaltContent,&
                         GetTotalWaterContent,&
                         Geteffectiverain,&
                         GetSumWaBal,&
                         GetRootZoneSalt,&
                         GetSimulation,&
                         GetCompartment,&
                         GetCompartment_theta,&
                         GetSoilLayer,&
                         GetIrrigation,&
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
                         GetIrriBeforeSeason,&
                         GetIrriAfterSeason


    use ac_run, only: GetIrriInterval,&
                      GetIrriInfoRecord1,&
                      GetIrriInfoRecord2,&

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
                    GetDayNri
              
              
    use ac_kinds, only: dp
!
! !DESCRIPTION:
!
!  This routine initializes the AC71 state variables with some
!  predefined values constantly for the entire domain. 
!
!EOP
 
    implicit none
    integer :: mtype
    integer :: t, l, n, i
    integer :: c, r
    

    do n=1, LIS_rc%nnest
        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: Ac71_coldstart -- cold-starting Ac71"
            do t=1, LIS_rc%npatch(n,mtype)
                do l=1, AC71_struc(n)%nsoil
                    AC71_struc(n)%ac71(t)%sh2o(l) = AC71_struc(n)%init_sh2o(l)
                enddo
                AC71_struc(n)%ac71(t)%irun = 1
                AC71_struc(n)%daynrinextclimaterecord = 1

                AC71_struc(n)%ac71(t)%InitializeRun = 1 ! gets 1 at end of year 

                do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                    AC71_struc(n)%ac71(t)%smc(l) = GetCompartment_theta(l)
                enddo              
            enddo
        endif
    
        LIS_rc%yr = LIS_rc%syr
        LIS_rc%mo = LIS_rc%smo
        LIS_rc%da = LIS_rc%sda
        LIS_rc%hr = LIS_rc%shr
        LIS_rc%mn = LIS_rc%smn
        LIS_rc%ss = LIS_rc%sss
        
        call LIS_date2time(LIS_rc%time, LIS_rc%doy, LIS_rc%gmt, LIS_rc%yr,      &
                           LIS_rc%mo, LIS_rc%da, LIS_rc%hr, LIS_rc%mn, LIS_rc%ss)
        write(LIS_logunit,*) "MSG: Ac71_coldstart -- ",     &
                             "Using the specified start time ", LIS_rc%time
    enddo
end subroutine Ac71_coldstart
