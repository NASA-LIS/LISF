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
    
    ! added by Shugong Wang
    integer :: isnow
    real, allocatable, dimension(:) :: zsnso 
    real, allocatable, dimension(:) :: tsno
    real, allocatable, dimension(:) :: snice
    real, allocatable, dimension(:) :: snliq 
    real, allocatable, dimension(:) :: zsoil 
    ! end add

    !EMK...Temporary arrays.
    real :: tmp_swe(1,1),tmp_tgxy(1,1),tmp_snodep(1,1)
    integer :: tmp_isnowxy(1,1)

    do n=1, LIS_rc%nnest
       isnow = -AC71_struc(n)%nsnow
        ! added by shugong
        allocate(zsnso(-AC71_struc(n)%nsnow+1:AC71_struc(n)%nsoil))
        allocate(tsno(-AC71_struc(n)%nsnow+1:0))
        allocate(snice(-AC71_struc(n)%nsnow+1:0))
        allocate(snliq(-AC71_struc(n)%nsnow+1:0))
        allocate(zsoil(AC71_struc(n)%nsoil))
        zsoil(1) = -AC71_struc(n)%sldpth(1)
        do l=2, AC71_struc(n)%nsoil
          zsoil(l) = zsoil(l-1) - AC71_struc(n)%sldpth(l) 
        enddo
        ! end add 

        if (trim(LIS_rc%startcode) .eq. "coldstart") then
            write(LIS_logunit,*) "MSG: Ac71_coldstart -- cold-starting Ac71"
            do t=1, LIS_rc%npatch(n,mtype)
                AC71_struc(n)%ac71(t)%albold = AC71_struc(n)%init_albold
                AC71_struc(n)%ac71(t)%sneqvo = AC71_struc(n)%init_sneqvo
                ! only soil temperature is intialized, snow temperature is calculated by snow_init 
                do l=1, AC71_struc(n)%nsoil
                    AC71_struc(n)%ac71(t)%sstc(AC71_struc(n)%nsnow+l) = AC71_struc(n)%init_stc(l)
                enddo
                do l=1, AC71_struc(n)%nsoil
                    AC71_struc(n)%ac71(t)%sh2o(l) = AC71_struc(n)%init_sh2o(l)
                enddo
                AC71_struc(n)%ac71(t)%tah = AC71_struc(n)%init_tah
                AC71_struc(n)%ac71(t)%eah = AC71_struc(n)%init_eah
                AC71_struc(n)%ac71(t)%fwet = AC71_struc(n)%init_fwet
                AC71_struc(n)%ac71(t)%canliq = AC71_struc(n)%init_canliq
                AC71_struc(n)%ac71(t)%canice = AC71_struc(n)%init_canice
                AC71_struc(n)%ac71(t)%tv = AC71_struc(n)%init_tv
                AC71_struc(n)%ac71(t)%tg = AC71_struc(n)%init_tg
                AC71_struc(n)%ac71(t)%qsnow = AC71_struc(n)%init_qsnow
                AC71_struc(n)%ac71(t)%snowh = AC71_struc(n)%init_snowh
                AC71_struc(n)%ac71(t)%sneqv = AC71_struc(n)%init_sneqv
                AC71_struc(n)%ac71(t)%zwt = AC71_struc(n)%init_zwt
                AC71_struc(n)%ac71(t)%wa = AC71_struc(n)%init_wa
                AC71_struc(n)%ac71(t)%wt = AC71_struc(n)%init_wt
                AC71_struc(n)%ac71(t)%wslake = AC71_struc(n)%init_wslake
                AC71_struc(n)%ac71(t)%lfmass = AC71_struc(n)%init_lfmass
                AC71_struc(n)%ac71(t)%rtmass = AC71_struc(n)%init_rtmass
                AC71_struc(n)%ac71(t)%stmass = AC71_struc(n)%init_stmass
                AC71_struc(n)%ac71(t)%wood = AC71_struc(n)%init_wood
                AC71_struc(n)%ac71(t)%stblcp = AC71_struc(n)%init_stblcp
                AC71_struc(n)%ac71(t)%fastcp = AC71_struc(n)%init_fastcp
                AC71_struc(n)%ac71(t)%lai = AC71_struc(n)%init_lai
                AC71_struc(n)%ac71(t)%sai = AC71_struc(n)%init_sai
                AC71_struc(n)%ac71(t)%cm = AC71_struc(n)%init_cm
                AC71_struc(n)%ac71(t)%ch = AC71_struc(n)%init_ch
                AC71_struc(n)%ac71(t)%tauss = AC71_struc(n)%init_tauss
                AC71_struc(n)%ac71(t)%smcwtd = AC71_struc(n)%init_smcwtd
                AC71_struc(n)%ac71(t)%deeprech = AC71_struc(n)%init_deeprech
                AC71_struc(n)%ac71(t)%rech = AC71_struc(n)%init_rech
                AC71_struc(n)%ac71(t)%zlvl = AC71_struc(n)%init_zlvl 
                AC71_struc(n)%ac71(t)%irun = 1
                AC71_struc(n)%daynrinextclimaterecord = 1

                AC71_struc(n)%ac71(t)%InitializeRun = 1 ! gets 1 at end of year 

                do l=1, AC71_struc(n)%ac71(t)%NrCompartments
                    AC71_struc(n)%ac71(t)%smc(l) = GetCompartment_theta(l)
                enddo

                do l=1, 40
                    AC71_struc(n)%ac71(t)%Tmin_ac_antecedent(l) = 0.0
                enddo
                ! added by shugong 
                zsnso = 0.0 
                !EMK...snow_init_71 is expecting several arrays which
                !are being passed as scalars.  Although no memory corruption
                !occurs here because of the declared array dimensions (all 1),
                !this is still technically a syntax error.  So, we will
                !copy the required fields to temporary arrays with the
                !correct declarations and pass those instead.
                !call snow_init_71(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                !                  AC71_struc(n)%nsnow,          & !input 
                !                  AC71_struc(n)%nsoil,          & !input 
                !                  zsoil,                            & !input
                !                  AC71_struc(n)%init_sneqv,     & !input
                !                  AC71_struc(n)%init_tg,        & !input
                !                  AC71_struc(n)%init_snowh,     & !input
                !                  zsnso, tsno, snice, snliq, isnow) ! output 
                tmp_swe(1,1) = AC71_struc(n)%init_sneqv
                tmp_tgxy(1,1) = AC71_struc(n)%init_tg
                tmp_snodep(1,1) = AC71_struc(n)%init_snowh
                call snow_init_71(1, 1, 1, 1, 1, 1, 1, 1,           & !input 
                                  AC71_struc(n)%nsnow,          & !input 
                                  AC71_struc(n)%nsoil,          & !input 
                                  zsoil,                            & !input
                                  tmp_swe,         & !input
                                  tmp_tgxy,        & !input
                                  tmp_snodep,      & !input
                                  zsnso, tsno, snice, snliq, &
                                  tmp_isnowxy) ! output
                isnow = tmp_isnowxy(1,1)
                AC71_struc(n)%ac71(t)%snowice(1:AC71_struc(n)%nsnow) = snice(-AC71_struc(n)%nsnow+1:0)
                AC71_struc(n)%ac71(t)%snowliq(1:AC71_struc(n)%nsnow) = snliq(-AC71_struc(n)%nsnow+1:0)
                AC71_struc(n)%ac71(t)%zss(1:AC71_struc(n)%nsnow+AC71_struc(n)%nsoil) = zsnso(-AC71_struc(n)%nsnow+1:AC71_struc(n)%nsoil) 
                AC71_struc(n)%ac71(t)%sstc(AC71_struc(n)%nsnow+isnow+1:AC71_struc(n)%nsnow) = tsno(isnow+1:0) 
                AC71_struc(n)%ac71(t)%isnow = isnow                

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
        deallocate(zsnso)
        deallocate(tsno)
        deallocate(snice)
        deallocate(snliq)
        deallocate(zsoil) 
    enddo
end subroutine Ac71_coldstart
