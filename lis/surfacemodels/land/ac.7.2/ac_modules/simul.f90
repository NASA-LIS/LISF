module ac_simul

use ac_global, only: ActiveCells, &
                     ac_zero_threshold, &
                     adjustedksstotoecsw, &
                     BMRange, &
                     CalculateAdjustedFC, &
                     CalculateETpot, &
                     CanopyCoverNoStressSF, &
                     ccinowaterstresssf, &
                     CheckForWaterTableInProfile, &
                     CO2Ref, &
                     CompartmentIndividual, &
                     CropStressParametersSoilFertility, &
                     datatype_daily, &
                     datatype_decadely, &
                     datatype_monthly, &
                     DetermineCNIandIII, &
                     DetermineDate, &
                     DetermineRootZoneSaltContent, &
                     DetermineRootzoneWC, &
                     DetermineSaltContent, &
                     ECeComp, &
                     ECswComp, &
                     EffectiveRainMethod_Percentage, &
                     EffectiveRainMethod_USDA, &
                     Equiv, &
                     EvapZmin, &
                     fAdjustedForCO2, &
                     GenerateDepthMode_FixDepth, &
                     GenerateTimeMode_AllDepl, &
                     GenerateTimeMode_AllRAW, &
                     GetCCiActual, &
                     GetCCiPrev, &
                     GetCCiTopEarlySen, &
                     GetCompartment, &
                     getcompartment_dayanaero, &
                     GetCompartment_Depo, &
                     GetCompartment_FCadj, &
                     GetCompartment_fluxout, &
                     GetCompartment_i, &
                     GetCompartment_Layer, &
                     GetCompartment_Salt, &
                     GetCompartment_Smax, &
                     GetCompartment_theta, &
                     GetCompartment_Thickness, &
                     GetCompartment_WFactor, &
                     GetCrop, &
                     GetCrop, &
                     GetCrop_aCoeff, &
                     GetCrop_AdaptedToCO2, &
                     GetCrop_AnaeroPoint, &
                     getcrop_bcoeff, &
                     GetCrop_CCEffectEvapLate, &
                     GetCrop_CCo, &
                     GetCrop_CCoAdjusted, &
                     GetCrop_CCsaltDistortion, &
                     GetCrop_CCx, &
                     GetCrop_CCxAdjusted, &
                     GetCrop_CCxRoot, &
                     GetCrop_CCxWithered, &
                     GetCrop_CDC, &
                     GetCrop_CGC, &
                     GetCrop_Day1, &
                     GetCrop_DayN, &
                     GetCrop_DaysToCCini, &
                     GetCrop_DaysToFlowering, &
                     GetCrop_DaysToFullCanopy, &
                     GetCrop_DaysToFullCanopySF, &
                     GetCrop_DaysToGermination, &
                     GetCrop_DaysToHarvest, &
                     GetCrop_DaysToSenescence, &
                     GetCrop_DeterminancyLinked, &
                     GetCrop_dHIdt, &
                     getcrop_dhimax, &
                     GetCrop_ECemax, &
                     GetCrop_ECemin, &
                     getcrop_fexcess, &
                     GetCrop_GDDaysToFlowering, &
                     GetCrop_GDDaysToFullCanopy, &
                     GetCrop_GDDaysToFullCanopySF, &
                     GetCrop_GDDaysToGermination, &
                     GetCrop_GDDaysToHarvest, &
                     GetCrop_GDDaysToSenescence, &
                     GetCrop_GDDCDC, &
                     GetCrop_GDDCGC, &
                     GetCrop_GDDLengthFlowering, &
                     GetCrop_GDtranspLow, &
                     GetCrop_HI, &
                     getcrop_hiincrease, &
                     GetCrop_KcDecline, &
                     GetCrop_KcTop, &
                     GetCrop_KsShapeFactorLeaf, &
                     GetCrop_KsShapeFactorSenescence, &
                     GetCrop_KsShapeFactorStomata, &
                     GetCrop_LengthFlowering, &
                     GetCrop_ModeCycle, &
                     GetCrop_pActStom, &
                     GetCrop_pdef, &
                     GetCrop_Planting, &
                     GetCrop_pLeafAct, &
                     GetCrop_pLeafDefLL, &
                     GetCrop_pLeafDefUL, &
                     GetCrop_pMethod, &
                     GetCrop_pPollination, &
                     GetCrop_pSenAct, &
                     GetCrop_pSenescence, &
                     GetCrop_ResponseECsw, &
                     GetCrop_RootMin, &
                     GetCrop_RootMin, &
                     Getcrop_smaxbot, &
                     GetCrop_SmaxTop, &
                     GetCrop_StressResponse, &
                     getcrop_stressresponse_calibrated, &
                     GetCrop_subkind, &
                     GetCrop_SumEToDelaySenescence, &
                     GetCrop_Tbase, &
                     GetCrop_Tcold, &
                     getcrop_theat, &
                     GetCrop_Tupper, &
                     GetCrop_WP, &
                     GetCrop_WPy, &
                     GetCrop_YearCCx, &
                     GetCRsalt, &
                     GetCRwater, &
                     GetDaySubmerged, &
                     GetDrain, &
                     GetEact, &
                     GetECdrain, &
                     GetECiAqua, &
                     GetECstorage, &
                     GetEpot, &
                     GetETo, &
                     GetEvapoEntireSoilSurface, &
                     GetGenerateDepthMode, &
                     GetGenerateTimeMode, &
                     GetInfiltrated, &
                     GetInfiltrated, &
                     GetIrriECw_PostSeason, &
                     GetIrriECw_PreSeason, &
                     GetIrrigation, &
                     GetIrriMode, &
                     GetManagement_BundHeight, &
                     GetManagement_CNcorrection, &
                     GetManagement_EffectMulchInS, &
                     GetManagement_EffectMulchOffS, &
                     GetManagement_FertilityStress, &
                     GetManagement_Mulch, &
                     GetManagement_RunoffON, &
                     GetManagement_SoilCoverAfter, &
                     GetManagement_SoilCoverBefore, &
                     GetManagement_WeedAdj, &
                     GetManagement_WeedDeltaRC, &
                     GetNrCompartments, &
                     GetPreDay, &
                     GetRain, &
                     GetRainRecord_DataType, &
                     GetRootingDepth, &
                     GetRootZoneSalt_ECe, &
                     GetRootZoneSalt_ECsw, &
                     GetRootZoneSalt_ECswFC, &
                     GetRootZoneSalt_KsSalt, &
                     GetRootZoneWC_Actual, &
                     GetRootZoneWC_FC, &
                     GetRootZoneWC_SAT, &
                     GetRootZoneWC_Thresh, &
                     GetRootZoneWC_WP, &
                     GetRootZoneWC_ZtopAct, &
                     GetRootZoneWC_ZtopFC, &
                     GetRootZoneWC_ZtopThresh, &
                     GetRootZoneWC_ZtopWP, &
                     GetRunoff, &
                     getsimulation_dayanaero, &
                     GetSimulation_DelayedDays, &
                     GetSimulation_EffectStress, &
                     GetSimulation_EffectStress_CDecline, &
                     GetSimulation_EffectStress_RedCCX, &
                     GetSimulation_EffectStress_RedCGC, &
                     GetSimulation_EffectStress_RedKsSto, &
                     GetSimulation_EffectStress_RedWP, &
                     GetSimulation_EvapStartStg2, &
                     GetSimulation_EvapWCsurf, &
                     GetSimulation_EvapZ, &
                     GetSimulation_Germinate, &
                     GetSimulation_HIfinal, &
                     GetSimulation_IrriECw, &
                     GetSimulation_ProtectedSeedling, &
                     GetSimulation_RCadj, &
                     GetSimulation_SalinityConsidered, &
                     GetSimulation_SCor, &
                     GetSimulation_Storage_Btotal, &
                     GetSimulation_SumEToStress, &
                     GetSimulation_SumGDD, &
                     GetSimulation_SWCtopSoilConsidered, &
                     GetSimulation_YearSeason, &
                     GetSimulParam_Beta, &
                     GetSimulParam_CNcorrection, &
                     GetSimulParam_DelayLowOxygen, &
                     GetSimulParam_EffectiveRain_Method, &
                     GetSimulParam_EffectiveRain_PercentEffRain, &
                     GetSimulParam_EffectiveRain_RootNrEvap, &
                     GetSimulParam_EffectiveRain_ShowersInDecade, &
                     GetSimulParam_EvapDeclineFactor, &
                     GetSimulParam_EvapZmax, &
                     GetSimulParam_IniAbstract, &
                     GetSimulParam_IrriFwInSeason, &
                     GetSimulParam_IrriFwOffSeason, &
                     GetSimulParam_pAdjFAO, &
                     GetSimulParam_PercCCxHIfinal, &
                     GetSimulParam_PercRAW, &
                     GetSimulParam_RootNrDF, &
                     GetSimulParam_RunoffDepth, &
                     GetSimulParam_SaltSolub, &
                     GetSimulParam_TAWGermination, &
                     GetSimulParam_Tmax, &
                     GetSimulParam_Tmin, &
                     GetSoil, &
                     GetSoil_CNvalue, &
                     GetSoil_NrSoilLayers, &
                     GetSoil_REW, &
                     GetSoilLayer_CRa, &
                     GetSoilLayer_CRb, &
                     GetSoilLayer_Dx, &
                     GetSoilLayer_FC, &
                     GetSoilLayer_GravelVol, &
                     GetSoilLayer_InfRate, &
                     GetSoilLayer_SaltMobility_i, &
                     GetSoilLayer_SAT, &
                     GetSoilLayer_SC, &
                     GetSoilLayer_SCP1, &
                     GetSoilLayer_tau, &
                     GetSoilLayer_Thickness, &
                     GetSoilLayer_UL, &
                     GetSoilLayer_WaterContent, &
                     GetSoilLayer_WP, &
                     GetSumWaBal_Biomass, &
                     GetSumWaBal_CRsalt, &
                     GetSumWaBal_CRwater, &
                     GetSumWaBal_Drain, &
                     GetSumWaBal_Eact, &
                     GetSumWaBal_ECropCycle, &
                     GetSumWaBal_Epot, &
                     GetSumWaBal_Infiltrated, &
                     GetSumWaBal_Irrigation, &
                     GetSumWaBal_Rain, &
                     GetSumWaBal_Runoff, &
                     GetSumWaBal_SaltIn, &
                     GetSumWaBal_SaltOut, &
                     GetSumWaBal_Tact, &
                     GetSumWaBal_Tpot, &
                     GetSumWaBal_TrW, &
                     GetSurf0, &
                     GetSurfaceStorage, &
                     GetTact, &
                     GetTactWeedInfested, &
                     GetTotalSaltContent_BeginDay, &
                     GetTotalSaltContent_EndDay, &
                     GetTotalWaterContent_BeginDay, &
                     GetTotalWaterContent_EndDay, &
                     GetTpot, &
                     GetWeedRC, &
                     GetZiAqua, &
                     HarvestIndexDay, &
                     HImultiplier, &
                     IrriMode_Generate, &
                     IrriMode_Inet, &
                     KsAny, &
                     KsAny, &
                     KsSalinity, &
                     KsTemperature, &
                     LengthCanopyDecline, &
                     max_No_compartments, &
                     MaxCRatDepth, &
                     modeCycle_CalendarDays, &
                     modeCycle_GDDays, &
                     MultiplierCCxSelfThinning, &
                     plant_regrowth, &
                     plant_seed, &
                     pMethod_FAOCorrection, &
                     pMethod_NoCorrection, &
                     rep_Crop, &
                     rep_EffectStress, &
                     rep_Soil, &
                     SaltSolutionDeposit, &
                     SaltSolutionDeposit, &
                     SetCCiActual, &
                     SetCCiPrev, &
                     SetCCiTopEarlySen, &
                     SetCompartment, &
                     SetCompartment_DayAnaero, &
                     SetCompartment_Depo, &
                     SetCompartment_fluxout, &
                     SetCompartment_i, &
                     SetCompartment_Salt, &
                     SetCompartment_theta, &
                     SetCompartment_WFactor, &
                     SetCrop_CCoAdjusted, &
                     SetCrop_CCxAdjusted, &
                     SetCrop_CCxWithered, &
                     SetCrop_DaysTOFullCanopySF, &
                     SetCrop_DaysToFullCanopySF, &
                     SetCrop_GDDaysToFullCanopySF, &
                     SetCrop_GDDaysToFullCanopySF, &
                     SetCrop_pActStom, &
                     SetCrop_pLeafAct, &
                     SetCrop_pPollination, &
                     SetCrop_pSenAct, &
                     SetCRsalt, &
                     SetCRwater, &
                     SetDaySubmerged, &
                     SetDrain, &
                     SetEact, &
                     SetECdrain, &
                     SetECstorage, &
                     SetEpot, &
                     SetEvapoEntireSoilSurface, &
                     SetInfiltrated, &
                     SetIrrigation, &
                     SetManagement_WeedDeltaRC, &
                     SetRootingdepth, &
                     SetRootZoneSalt_ECe, &
                     SetRootZoneSalt_ECsw, &
                     SetRootZoneSalt_ECswFC, &
                     SetRootZoneSalt_KsSalt, &
                     SetRunoff, &
                     SetRunoff, &
                     SetSaltInfiltr, &
                     setsimulation_dayanaero, &
                     SetSimulation_DelayedDays, &
                     SetSimulation_EffectStress, &
                     SetSimulation_EffectStress_CDecline, &
                     SetSimulation_EffectStress_RedCCx, &
                     SetSimulation_EffectStress_RedCGC, &
                     SetSimulation_EffectStress_RedKsSto, &
                     SetSimulation_EffectStress_RedWP, &
                     SetSimulation_EvapLimitON, &
                     SetSimulation_EvapStartStg2, &
                     SetSimulation_EvapWCsurf, &
                     SetSimulation_EvapZ, &
                     SetSimulation_Germinate, &
                     SetSimulation_HIfinal, &
                     SetSimulation_ProtectedSeedling, &
                     SetSimulation_Storage_Btotal, &
                     SetSimulation_SumEToStress, &
                     setsimulation_sumgdd, &
                     SetSimulation_SWCtopSoilConsidered, &
                     SetSoilLayer_WaterContent, &
                     SetSumWaBal_CRsalt, &
                     SetSumWaBal_CRwater, &
                     SetSumWaBal_Drain, &
                     SetSumWaBal_Eact, &
                     SetSumWaBal_ECropCycle, &
                     SetSumWaBal_Epot, &
                     SetSumWaBal_Infiltrated, &
                     SetSumWaBal_Irrigation, &
                     SetSumWaBal_Rain, &
                     SetSumWaBal_Runoff, &
                     SetSumWaBal_SaltIn, &
                     SetSumWaBal_SaltOut, &
                     SetSumWaBal_Tact, &
                     SetSumWaBal_Tpot, &
                     SetSumWaBal_TrW, &
                     SetSurf0, &
                     SetSurfaceStorage, &
                     SetTact, &
                     SetTact, &
                     SetTotalSaltContent_BeginDay, &
                     SetTotalSaltContent_EndDay, &
                     SetTotalSaltContent_ErrorDay, &
                     SetTotalWaterContent_BeginDay, &
                     SetTotalWaterContent_EndDay, &
                     SetTotalWaterContent_ErrorDay, &
                     SetTpot, &
                     SoilEvaporationReductionCoefficient, &
                     Subkind_Forage, &
                     subkind_Grain,  &
                     subkind_Grain, &
                     subkind_Tuber, &
                     subkind_Vegetative, &
                     TimeToMaxCanopySF, &
                     undef_double, &
                     undef_int
use ac_kinds, only:  dp, &
                     int8, &
                     int32, &
                     intEnum
use ac_tempprocessing, only: CropStressParametersSoilSalinity, &
                             GrowingDegreeDays, &
                             SumCalendarDays
use ac_utils, only: roundc
implicit none


integer(intEnum), parameter :: whichtheta_AtSat = 0
    !! index of AtSat in whichtheta enumerated type
integer(intEnum), parameter :: whichtheta_AtFC = 1
    !! index of AtFC in whichtheta enumerated type
integer(intEnum), parameter :: whichtheta_AtWP = 2
    !! index of AtWP in whichtheta enumerated type
integer(intEnum), parameter :: whichtheta_AtAct = 3
    !! index of AtAct in whichtheta enumerated type


integer(intEnum), parameter :: control_begin_day = 0
    !! index of beginday in control enumerated type
integer(intEnum), parameter :: control_end_day = 1
    !! index of endday in control enumerated type


contains


real(dp) function GetCDCadjustedNoStressNew(CCx, CDC, CCxAdjusted)
    real(dp), intent(in) :: CCx
    real(dp), intent(in) :: CDC
    real(dp), intent(in) :: CCxAdjusted

    real(dp) :: CDCadjusted

    CDCadjusted = CDC * ((CCxadjusted+2.29_dp)/(CCx+2.29_dp))
    GetCDCadjustedNoStressNew = CDCadjusted
end function GetCDCadjustedNoStressNew


subroutine AdjustpLeafToETo(EToMean, pLeafULAct, pLeafLLAct)
    real(dp), intent(in) :: EToMean
    real(dp), intent(inout) :: pLeafULAct
    real(dp), intent(inout) :: pLeafLLAct


    pLeafLLAct = GetCrop_pLeafDefLL()
    pLeafULAct = GetCrop_pLeafDefUL()
    if (GetCrop_pMethod() == pMethod_FAOCorrection) then
        pLeafLLAct = GetCrop_pLeafDefLL() + GetSimulParam_pAdjFAO()* 0.04_dp &
                        *(5._dp-EToMean)*log10(10._dp-9._dp*GetCrop_pLeafDefLL())
        if (pLeafLLAct > 1.0) then
            pLeafLLAct = 1.0_dp
        end if
        if (pLeafLLAct < 0) then
            pLeafLLAct = 0._dp
        end if
        pLeafULAct = GetCrop_pLeafDefUL() + GetSimulParam_pAdjFAO()* 0.04_dp &
                        *(5._dp-EToMean)*log10(10._dp-9._dp*GetCrop_pLeafDefUL())
        if (pLeafULAct > 1.0) then
            pLeafULAct = 1.0_dp
        end if
        if (pLeafULAct < 0) then
            pLeafULAct = 0._dp
        end if
    end if
end subroutine AdjustpLeafToETo


subroutine DeterminePotentialBiomass(VirtualTimeCC, SumGDDadjCC, CO2i, GDDayi, &
                                              CCxWitheredTpotNoS, BiomassUnlim)
    integer(int32), intent(in) :: VirtualTimeCC
    real(dp), intent(in) :: SumGDDadjCC
    real(dp), intent(in) :: CO2i
    real(dp), intent(in) :: GDDayi
    real(dp), intent(inout) :: CCxWitheredTpotNoS
    real(dp), intent(inout) :: BiomassUnlim

    real(dp) :: CCiPot,  WPi, fSwitch, TpotForB, EpotTotForB
    integer(int32) :: DAP, DaysYieldFormation, DayiAfterFlowering
    real(dp) :: Tmin_local, Tmax_local

    ! potential biomass - unlimited soil fertiltiy
    ! 1. - CCi
    CCiPot = CanopyCoverNoStressSF((VirtualTimeCC + GetSimulation_DelayedDays() + 1), &
                                  GetCrop_DaysToGermination(), GetCrop_DaysToSenescence(), &
                                  GetCrop_DaysToHarvest(), GetCrop_GDDaysToGermination(), &
                                  GetCrop_GDDaysToSenescence(), GetCrop_GDDaysToHarvest(), &
                                  GetCrop_CCo(), GetCrop_CCx(), GetCrop_CGC(), &
                                  GetCrop_CDC(), GetCrop_GDDCGC(), GetCrop_GDDCDC(), &
                                  SumGDDadjCC, GetCrop_ModeCycle(), 0_int8, 0_int8)
    if (CCiPot < 0._dp) then
        CCiPot = 0._dp
    end if
    if (CCiPot > CCxWitheredTpotNoS) then
        CCxWitheredTpotNoS = CCiPot
    end if

    ! 2. - Calculation of Tpot
    if (GetCrop_ModeCycle() == modeCycle_CalendarDays) then
        DAP = VirtualTimeCC
    else
        ! growing degree days
        Tmin_local = GetSimulParam_Tmin()
        Tmax_local = GetSimulParam_Tmax()
        DAP = SumCalendarDays(roundc(SumGDDadjCC, mold=1), GetCrop_Day1(), GetCrop_Tbase(), &
                    GetCrop_Tupper(), Tmin_local, Tmax_local)
        DAP = DAP + GetSimulation_DelayedDays() ! are not considered when working with GDDays
    end if
    call CalculateETpot(DAP, GetCrop_DaysToGermination(), GetCrop_DaysToFullCanopy(), &
                   GetCrop_DaysToSenescence(), GetCrop_DaysToHarvest(), 0, CCiPot, &
                   GetETo(), GetCrop_KcTop(), GetCrop_KcDecline(), GetCrop_CCx(), &
                   CCxWitheredTpotNoS, real(GetCrop_CCEffectEvapLate(), kind=dp), CO2i, GDDayi, &
                   GetCrop_GDtranspLow(), TpotForB, EpotTotForB)

    ! 3. - WPi for that day
    ! 3a - given WPi
    WPi = (GetCrop_WP()/100._dp)
    ! 3b - WPi decline in reproductive stage  (works with calendar days)
    if (((GetCrop_subkind() == subkind_Grain) .or. (GetCrop_subkind() == subkind_Tuber)) &
        .and. (GetCrop_WPy() < 100._dp) .and. (GetCrop_dHIdt() > 0._dp) &
        .and. (VirtualTimeCC >= GetCrop_DaysToFlowering())) then
        ! WPi in reproductive stage
        fSwitch = 1._dp
        DaysYieldFormation = roundc(GetCrop_HI()/GetCrop_dHIdt(), mold=1)
        DayiAfterFlowering = VirtualTimeCC - GetCrop_DaysToFlowering()
        if ((DaysYieldFormation > 0) .and. (DayiAfterFlowering < &
                                              (DaysYieldFormation/3._dp))) then
            fSwitch = DayiAfterFlowering/(DaysYieldFormation/3._dp)
        end if
        WPi =  WPi * (1._dp - (1._dp-GetCrop_WPy()/100._dp)*fSwitch)
    end if
    ! 3c - adjustment WPi for CO2
    if (roundc(100._dp*CO2i, mold=1) /= roundc(100._dp*CO2Ref, mold=1)) then
        WPi = WPi * fAdjustedForCO2(CO2i, GetCrop_WP(), GetCrop_AdaptedToCO2())
    end if

    ! 4. - Potential Biomass
    if (GetETo() > 0._dp) then
        BiomassUnlim = BiomassUnlim + WPi * TpotForB/real(GetETo(), kind=dp) ! ton/ha
    end if
end subroutine DeterminePotentialBiomass


subroutine DetermineBiomassAndYield(dayi, ETo, TminOnDay, TmaxOnDay, CO2i, &
                                    GDDayi, Tact, SumKcTop, CGCref, GDDCGCref, &
                                    Coeffb0, Coeffb1, Coeffb2, FracBiomassPotSF, &
                                    Coeffb0Salt, Coeffb1Salt, Coeffb2Salt, &
                                    AverageSaltStress, SumGDDadjCC, CCtot, &
                                    FracAssim, VirtualTimeCC, SumInterval, &
                                    Biomass, BiomassPot, BiomassUnlim, &
                                    BiomassTot, YieldPart, WPi, HItimesBEF, &
                                    ScorAT1, ScorAT2, HItimesAT1, HItimesAT2, &
                                    HItimesAT, alfa, alfaMax, SumKcTopStress, &
                                    SumKci, &
                                    WeedRCi, CCw, Trw, StressSFadjNEW, &
                                    PreviousStressLevel, StoreAssimilates, &
                                    MobilizeAssimilates, AssimToMobilize, &
                                    AssimMobilized, Bin, Bout, TESTVAL)

    integer(int32), intent(in) :: dayi
    real(dp), intent(in) :: ETo
    real(dp), intent(in) :: TminOnDay
    real(dp), intent(in) :: TmaxOnDay
    real(dp), intent(in) :: CO2i
    real(dp), intent(in) :: GDDayi
    real(dp), intent(in) :: Tact
    real(dp), intent(in) :: SumKcTop
    real(dp), intent(in) :: CGCref
    real(dp), intent(in) :: GDDCGCref
    real(dp), intent(in) :: Coeffb0
    real(dp), intent(in) :: Coeffb1
    real(dp), intent(in) :: Coeffb2
    real(dp), intent(in) :: FracBiomassPotSF
    real(dp), intent(in) :: Coeffb0Salt
    real(dp), intent(in) :: Coeffb1Salt
    real(dp), intent(in) :: Coeffb2Salt
    real(dp), intent(in) :: AverageSaltStress
    real(dp), intent(in) :: SumGDDadjCC
    real(dp), intent(in) :: CCtot
    real(dp), intent(inout) :: FracAssim
    integer(int32), intent(in) :: VirtualTimeCC
    integer(int32), intent(in) :: SumInterval
    real(dp), intent(inout) :: Biomass
    real(dp), intent(inout) :: BiomassPot
    real(dp), intent(inout) :: BiomassUnlim
    real(dp), intent(inout) :: BiomassTot
    real(dp), intent(inout) :: YieldPart
    real(dp), intent(inout) :: WPi
    real(dp), intent(inout) :: HItimesBEF
    real(dp), intent(inout) :: ScorAT1
    real(dp), intent(inout) :: ScorAT2
    real(dp), intent(inout) :: HItimesAT1
    real(dp), intent(inout) :: HItimesAT2
    real(dp), intent(inout) :: HItimesAT
    real(dp), intent(inout) :: alfa
    real(dp), intent(inout) :: alfaMax
    real(dp), intent(inout) :: SumKcTopStress
    real(dp), intent(inout) :: SumKci
    real(dp), intent(inout) :: WeedRCi
    real(dp), intent(inout) :: CCw
    real(dp), intent(inout) :: Trw
    integer(int32), intent(inout) :: StressSFadjNEW
    integer(int32), intent(inout) :: PreviousStressLevel
    logical, intent(inout) :: StoreAssimilates
    logical, intent(inout) :: MobilizeAssimilates
    real(dp), intent(inout) :: AssimToMobilize
    real(dp), intent(inout) :: AssimMobilized
    real(dp), intent(inout) :: Bin
    real(dp), intent(inout) :: Bout
    real(dp), intent(inout) :: TESTVAL

    real(dp), parameter :: TempRange = 5._dp
    real(dp), parameter :: k = 2._dp

    real(dp) :: RatioBM, RBM, HItimesTotal, pLeafULAct, pLeafLLAct, &
                pStomatULAct, pLL, Ksleaf, Ksstomatal, KsPolWS, KsPolCs, &
                KsPolHs, KsPol, Wrel, Dcor, fFlor, fSwitch, fCCx,WPsf, WPunlim, &
                BioAdj,CCtotStar, CCwStar, croppol_temp
    integer(int32) :: tmax1, tmax2, DayCor, DayiAfterFlowering, &
                      DaysYieldFormation, wdrc_temp, HIfinal_temp
    integer(int8) :: PercentLagPhase
    logical :: SWCtopSoilConsidered_temp

    TESTVAL = undef_int

    ! 0. Reference HarvestIndex for that day (alfa in percentage) + Information on PercentLagPhase (for estimate WPi)
    if ((GetCrop_subkind() == Subkind_Tuber) .or. (GetCrop_Subkind() == Subkind_grain) &
        .or. (GetCrop_Subkind() == Subkind_Vegetative) .or. &
                                    (GetCrop_Subkind() == Subkind_Forage)) then
        ! DaysToFlowering corresponds with Tuberformation
        if ((GetCrop_Subkind() == Subkind_Vegetative) .and. &
                                (GetCrop_Planting() == plant_Regrowth) .or. &
                                (GetCrop_Subkind() == Subkind_Forage) .and. &
                                (GetCrop_Planting() == plant_Regrowth)) then
            alfa = GetCrop_HI()
        else
            HIfinal_temp = GetSimulation_HIfinal()
            alfa = HarvestIndexDay((dayi-GetCrop_Day1()), GetCrop_DaysToFlowering(), &
                                   GetCrop_HI(), GetCrop_dHIdt(), GetCCiactual(), &
                                   GetCrop_CCxAdjusted(), GetCrop_CCxWithered(), GetSimulParam_PercCCxHIfinal(), &
                                   GetCrop_Planting(), PercentLagPhase, HIfinal_temp)
            call SetSimulation_HIfinal(HIfinal_temp)
        end if
    end if


    WPi = (GetCrop_WP()/100._dp)

    ! 1. biomass
    if (ETo > 0._dp) then
        ! 1.1 WPi for that day
        ! 1.1a - given WPi
        WPi = (GetCrop_WP()/100._dp)
        ! 1.1b - adjustment WPi for reproductive stage (works with calendar days)
        if (((GetCrop_subkind() == Subkind_Tuber) .or. &
                    (GetCrop_Subkind() == Subkind_grain)) .and. (alfa > 0._dp)) then
            ! WPi switch to WP for reproductive stage
            fSwitch = 1._dp
            DaysYieldFormation = roundc(GetCrop_HI()/GetCrop_dHIdt(), mold=1)
            if (DaysYieldFormation > 0) then
                if (GetCrop_DeterminancyLinked()) then
                    fSwitch = PercentLagPhase/100._dp
                else
                    DayiAfterFlowering = dayi - GetSimulation_DelayedDays() - &
                                      GetCrop_Day1() - GetCrop_DaysToFlowering()
                    if (DayiAfterFlowering < (DaysYieldFormation/3._dp)) then
                        fSwitch = DayiAfterFlowering/(DaysYieldFormation/3._dp)
                    end if
                end if
            end if
            WPi =  WPi * (1._dp - (1._dp-GetCrop_WPy()/100._dp)*fSwitch)  ! switch in Lag Phase
        end if


        ! 1.1c - adjustment WPi for CO2
        if (roundc(100._dp*CO2i, mold=1) /= roundc(100._dp*CO2Ref, mold=1)) then
            WPi = WPi * fAdjustedForCO2(CO2i, GetCrop_WP(), GetCrop_AdaptedToCO2())
        end if


        ! 1.1d - adjustment WPi for Soil Fertility
        WPsf = WPi          ! no water stress, but fertility stress
        WPunlim = WPi       ! no water stress, no fertiltiy stress
        if (GetSimulation_EffectStress_RedWP() > 0._dp) then ! Reductions are zero if no fertility stress
            ! water stress and fertility stress
            if ((SumKci/real(SumKcTopStress, dp)) < 1._dp) then
                if (ETo > 0._dp) then
                    SumKci = SumKci + Tact/ETo
                end if
                if (SumKci > 0._dp) then
                    WPi = WPi * (1._dp - (GetSimulation_EffectStress_RedWP()/100._dp) &
                                            * exp(k*log(SumKci/SumKcTopStress)) )
                end if
            else
                WPi = WPi * (1._dp - GetSimulation_EffectStress_RedWP()/100._dp)
            end if
        elseif (ETo > 0._dp) then
            SumKci = SumKci + Tact/ETo
        end if


        ! 1.2 actual biomass
        if ((GetSimulation_RCadj() > 0._dp) .and. (roundc(CCtot*10000._dp, mold=1) > 0._dp)) then
            ! weed infestation
            ! green canopy cover of the crop in weed-infested field
            if (GetManagement_WeedDeltaRC() /= 0) then
                if (GetCrop_subkind() == Subkind_Forage) then
                    fCCx = MultiplierCCxSelfThinning(int(GetSimulation_YearSeason(), kind=int32), &
                                           int(GetCrop_YearCCx(), kind=int32), GetCrop_CCxRoot())
                else
                    fCCx = 1._dp
                end if
                wdrc_temp = GetManagement_WeedDeltaRC()
                WeedRCi = GetWeedRC(VirtualTimeCC, SumGDDadjCC, fCCx, &
                GetSimulation_RCadj(), GetManagement_WeedAdj(), wdrc_temp, &
                GetCrop_DaysToFullCanopySF(), GetCrop_DaysToSenescence(), &
                GetCrop_GDDaysToFullCanopySF(), GetCrop_GDDaysToSenescence(), &
                GetCrop_ModeCycle())
                call SetManagement_WeedDeltaRC(wdrc_temp)
            else
                WeedRCi = GetSimulation_RCadj()
            end if
            CCw = CCtot * (1._dp-WeedRCi/100._dp)
            ! correction for micro-advection
            CCtotStar = 1.72_dp*CCtot - 1._dp*(CCtot*CCtot) + 0.30_dp*(CCtot*CCtot*CCtot)
            if (CCtotStar < 0._dp) then
                CCtotStar = 0._dp
            end if
            if (CCtotStar > 1) then
                CCtotStar = 1._dp
            end if
            if (CCw > 0.0001_dp) then
                CCwStar = CCw + (CCtotStar - CCtot)
            else
                CCwStar = 0._dp
            end if
            ! crop transpiration in weed-infested field
            if (CCtotStar <= 0.0001_dp) then
                TrW = 0._dp
            else
                TrW = Tact * (CCwStar/CCtotStar)
            end if
            ! crop biomass in weed-infested field
            Biomass = Biomass + WPi *(TrW/ETo)  ! ton/ha
        else
            WeedRCi = 0.0_dp
            CCw = CCtot
            TrW = Tact
            Biomass = Biomass + WPi *(Tact/ETo)  ! ton/ha
        end if

        ! Transfer of assimilates
        if (GetCrop_subkind() == subkind_Forage) then
            ! only for perennial herbaceous forage crops
            ! 1. Mobilize assimilates at start of season
            if (MobilizeAssimilates .eqv. .true.) then
                ! mass to mobilize
                Bin = FracAssim * WPi *(TrW/ETo)  ! ton/ha
                if ((AssimMobilized + Bin) > AssimToMobilize) then
                    Bin = AssimToMobilize - AssimMobilized
                end if
                ! cumulative mass mobilized
                AssimMobilized = AssimMobilized + Bin
                ! switch mobilize off when all mass is transfered
                if (roundc(1000._dp*AssimToMobilize, mold=1) &
                   <= roundc(1000._dp *AssimMobilized, mold=1)) then
                    MobilizeAssimilates = .false.
                end if
            end if
            ! 2. Store assimilates at end of season
            if (StoreAssimilates .eqv. .true.) then
                ! mass to store
                Bout = FracAssim * WPi *(TrW/ETo)  ! ton/ha
                ! cumulative mass stored
                call SetSimulation_Storage_Btotal(GetSimulation_Storage_Btotal() + Bout)
            end if
            TESTVAL = FracAssim
        end if

        Biomass = Biomass + Bin - Bout  ! ton/ha ! correction for transferred assimilates

        ! actual total biomass (crop and weeds)
        BiomassTot = BiomassTot + WPi *(Tact/ETo)  ! ton/ha  for dynamic adjustment of soil fertility stress
        BiomassTot = BiomassTot + Bin - Bout ! correction for transferred assimilates

        ! 1.3 potential biomass - unlimited soil fertiltiy
        BiomassUnlim = BiomassUnlim + Bin - Bout ! correction for transferred assimilates

    end if

    ! 1.4 potential biomass for given soil fertility
    BiomassPot =  FracBiomassPotSF * BiomassUnlim ! ton/ha

    ! 2. yield
    tmax1 = undef_int
    if ((GetCrop_subkind() == subkind_Tuber) .or. (GetCrop_Subkind() == subkind_Grain)) then
        ! DaysToFlowering corresponds with Tuberformation
        if (dayi > (GetSimulation_DelayedDays() + GetCrop_Day1() + GetCrop_DaysToFlowering())) then
            ! calculation starts when flowering has started

            ! 2.2 determine HImultiplier at the start of flowering
            ! effect of water stress before flowering (HItimesBEF)
            if (HItimesBEF < - 0.1_dp) then
                ! i.e. undefined at the start of flowering
                if (BiomassPot < 0.0001_dp) then
                    HItimesBEF = 1._dp
                else
                    RatioBM = Biomass/BiomassPot
                    ! Not correct if weed infestation and no fertility stress
                    ! for that case BiomassPot might be larger (but cannot be calculated since WP is unknown)
                    if (RatioBM > 1._dp) then
                        RatioBM = 1_dp
                    end if
                    RBM = BMRange(int(GetCrop_HIincrease(), kind=int32))
                    HItimesBEF = HImultiplier(RatioBM, RBM, GetCrop_HIincrease())
                end if
                if (GetCCiActual() <= 0.01_dp) then
                    if ((GetCrop_CCxWithered() > 0._dp) &
                        .and. (GetCCiActual() < GetCrop_CCxWithered())) then
                        HItimesBEF = 0._dp ! no green canopy cover left at start of flowering;
                    else
                        HItimesBEF = 1._dp
                    end if
                end if
            end if

            ! 2.3 Relative water content for that day
            SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
            call DetermineRootZoneWC(GetRootingDepth(), SWCtopSoilConsidered_temp)
            call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
            if (GetSimulation_SWCtopSoilConsidered() .eqv. .true.) then ! top soil is relative wetter than total root zone
                Wrel = (GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopAct())/ &
                       (GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopWP()) ! top soil
            else
                Wrel = (GetRootZoneWC_FC() - GetRootZoneWC_Actual())/ &
                       (GetRootZoneWC_FC() - GetRootZoneWC_WP()) ! total root zone
            end if

            ! 2.4 Failure of Pollination during flowering (alfaMax in percentage)
            if (GetCrop_Subkind() == Subkind_grain) then ! - only valid for fruit/grain crops (flowers)
                if ((dayi <= (GetSimulation_DelayedDays() + GetCrop_Day1() + &
                   GetCrop_DaysToFlowering() + GetCrop_LengthFlowering())) & ! calculation limited to flowering period
                    .and. ((GetCCiactual()*100._dp) > GetSimulParam_PercCCxHIfinal())) then
                    ! sufficient green canopy remains
                    ! 2.4a - Fraction of flowers which are flowering on day  (fFlor)
                    fFlor = FractionFlowering(dayi)
                    ! 2.4b - Ks(pollination) water stress
                    pLL = 1._dp
                    croppol_temp = GetCrop_pPollination()
                    KsPolWS = KsAny(Wrel, croppol_temp, pLL, 0._dp)
                    ! 2.4c - Ks(pollination) cold stress
                    KsPolCS = KsTemperature((GetCrop_Tcold()-TempRange), real(GetCrop_Tcold(), kind=dp), TminOnDay)
                    ! 2.4d - Ks(pollination) heat stress
                    KsPolHS = KsTemperature((GetCrop_Theat()+TempRange), real(GetCrop_Theat(), kind=dp), TmaxOnDay)
                    ! 2.4e - Adjust alfa
                    KsPol = KsPolWS
                    if (KsPol > KsPolCS) then
                        KsPol = KsPolCS
                    end if
                    if (KsPol > KsPolHS) then
                        KsPol = KsPolHS
                    end if
                    alfaMax = alfaMax + (KsPol * (1 + GetCrop_fExcess()/100._dp) * fFlor * GetCrop_HI())
                    if (alfaMax > GetCrop_HI()) then
                        alfaMax = GetCrop_HI()
                    end if
                end if
            else
                alfaMax = GetCrop_HI() ! for Tuber crops (no flowering)
            end if

            ! 2.5 determine effect of water stress affecting leaf expansion after flowering
            ! from start flowering till end of determinancy
            if (GetCrop_DeterminancyLinked()) then
                tmax1 = roundc(GetCrop_LengthFlowering()/2._dp, mold=1)
            else
                tmax1 = (GetCrop_DaysToSenescence() - GetCrop_DaysToFlowering())
            end if
            if ((HItimesBEF > 0.99_dp) & ! there is green canopy cover at start of flowering;
                .and. (dayi <= (GetSimulation_DelayedDays() + GetCrop_Day1() &
                      + GetCrop_DaysToFlowering()+ tmax1)) & ! and not yet end period
                .and. (tmax1 > 0) & ! otherwise no effect
                .and. (roundc(GetCrop_aCoeff(), mold=1) /= undef_int) & ! otherwise no effect
                ! possible precision issue in pascal code
                ! added -epsilon(0._dp) for zero-diff with pascal version output
                .and. (GetCCiactual() > (0.001_dp-epsilon(0._dp)))) then ! and as long as green canopy cover remains (for correction to stresses)
                ! determine KsLeaf
                call AdjustpLeafToETo(ETo, pLeafULAct, pLeafLLAct)
                Ksleaf = KsAny(Wrel, pLeafULAct, pLeafLLAct, GetCrop_KsShapeFactorLeaf())
                ! daily correction
                Dcor = (1._dp + (1._dp-Ksleaf)/GetCrop_aCoeff())
                ! weighted correction
                ScorAT1 = ScorAT1 + Dcor/tmax1
                DayCor = dayi - (GetSimulation_DelayedDays() + GetCrop_Day1() + GetCrop_DaysToFlowering())
                HItimesAT1  = (tmax1*1._dp/DayCor) * ScorAT1
            end if

            ! 2.6 determine effect of water stress affecting stomatal closure after flowering
            ! during yield formation
            if (GetCrop_dHIdt() > 99._dp) then
                tmax2 = 0
            else
                tmax2 = roundc(GetCrop_HI()/GetCrop_dHIdt(), mold=1)
            end if
            if ((HItimesBEF > 0.99_dp) & ! there is green canopy cover at start of flowering;
                .and. (dayi <= (GetSimulation_DelayedDays() + GetCrop_Day1() &
                      + GetCrop_DaysToFlowering() + tmax2)) & ! and not yet end period
                .and. (tmax2 > 0) & ! otherwise no effect
                .and. (roundc(GetCrop_bCoeff(), mold=1) /= undef_int) & ! otherwise no effect
                ! possible precision issue in pascal code
                ! added -epsilon(0._dp) for zero-diff with pascal version output
                .and. (GetCCiactual() > (0.001_dp-epsilon(0._dp)))) then ! and as long as green canopy cover remains (for correction to stresses)
                ! determine KsStomatal
                call AdjustpStomatalToETo(ETo, pStomatULAct)
                pLL = 1._dp
                Ksstomatal = KsAny(Wrel, pStomatULAct, pLL, GetCrop_KsShapeFactorStomata())
                ! daily correction
                if (Ksstomatal > 0.001_dp) then
                    Dcor = (exp(0.10_dp*log(Ksstomatal))) * (1._dp-(1._dp-Ksstomatal)/GetCrop_bCoeff())
                else
                    Dcor = 0._dp
                end if
                ! weighted correction
                ScorAT2 = ScorAT2 + Dcor/tmax2
                DayCor = dayi - (GetSimulation_DelayedDays() + GetCrop_Day1() + GetCrop_DaysToFlowering())
                HItimesAT2  = (tmax2*1._dp/DayCor) * ScorAT2
            end if

            ! 2.7 total multiplier after flowering
            if ((tmax2 == 0) .and. (tmax1 == 0)) then
                HItimesAT = 1._dp
            else
                if (tmax2 == 0) then
                    HItimesAT = HItimesAT1
                else
                    if (tmax1 == 0) then
                        HItimesAT = HItimesAT2
                    elseif (tmax1 <= tmax2) then
                        HItimesAT = HItimesAT2 * ((tmax1*HItimesAT1 + (tmax2-tmax1))/tmax2)
                        if (roundc(GetCrop_bCoeff(), mold=1) == undef_int) then
                            HItimesAT = HItimesAT1
                        end if
                        if (roundc(GetCrop_aCoeff(), mold=1) == undef_int) then
                            HItimesAT = HItimesAT2
                        end if
                    else
                        HItimesAT = HItimesAT1 * ((tmax2*HItimesAT2 + (tmax1-tmax2))/tmax1)
                        if (roundc(GetCrop_bCoeff(), mold=1) == undef_int) then
                            HItimesAT = HItimesAT1
                        end if
                        if (roundc(GetCrop_aCoeff(), mold=1) == undef_int) then
                            HItimesAT = HItimesAT2
                        end if
                    end if
                end if
            end if

            ! 2.8 Limit HI to allowable maximum increase
            HItimesTotal = HItimesBEF * HItimesAT
            if (HItimesTotal > (1._dp +(GetCrop_DHImax()/100._dp))) then
                HItimesTotal = 1._dp +(GetCrop_DHImax()/100._dp)
            end if

            ! 2.9 Yield
            if (alfaMax >= alfa) then
                YieldPart = Biomass * HItimesTotal*(alfa/100._dp)
            else
                YieldPart = Biomass * HItimesTotal*(alfaMax/100._dp)
            end if
        end if
    end if

    ! 2bis. yield leafy vegetable crops and forage crops
    if ((GetCrop_subkind() == subkind_Vegetative) .or. (GetCrop_subkind() == subkind_Forage)) then
        if (dayi >= (GetSimulation_DelayedDays() + GetCrop_Day1() + GetCrop_DaysToFlowering())) then
            ! calculation starts at crop day 1 (since days to flowering is 0)
            if (roundc(100._dp*ETo, mold=1)> 0._dp) then
                ! with correction for transferred assimilates
                if (GetSimulation_RCadj() > 0._dp) then
                    YieldPart = YieldPart + (WPi*(TrW/ETo) + Bin - Bout) * (alfa/100._dp)
                else
                    YieldPart = YieldPart + (WPi*(Tact/ETo) + Bin - Bout) * (alfa/100._dp)
                end if
            end if
        end if
    end if


    ! 3. Dynamic adjustment of soil fertility stress
    if ((GetManagement_FertilityStress() > 0) .and. (BiomassUnlim > 0.001_dp) &
                            .and. GetCrop_StressResponse_Calibrated()) then
        BioAdj = 100._dp * (FracBiomassPotSF + (FracBiomassPotSF - BiomassTot/BiomassUnlim))
        if (BioAdj >= 100._dp) then
            StressSFadjNEW = 0
        else
            if (BioAdj <= epsilon(1._dp)) then
                StressSFadjNEW = 80
            else
                StressSFadjNEW = roundc(Coeffb0 + Coeffb1*BioAdj + Coeffb2*BioAdj*BioAdj, &
                                        mold=1_int8)
                if (StressSFadjNEW < 0) then
                    StressSFadjNEW = GetManagement_FertilityStress()
                end if
                if (StressSFadjNEW > 80) then
                    StressSFadjNEW = 80
                end if
            end if
            if (StressSFadjNEW > GetManagement_FertilityStress()) then
                StressSFadjNEW = GetManagement_FertilityStress()
            end if
        end if
        if ((GetCrop_Subkind() == Subkind_grain) .and. GetCrop_DeterminancyLinked() &
            .and. (dayi > (GetSimulation_DelayedDays() + GetCrop_Day1() &
                                + GetCrop_DaysToFlowering() + tmax1))) then
            ! potential vegetation period is exceeded
            if (StressSFadjNEW < PreviousStressLevel) then
                StressSFadjNEW = PreviousStressLevel
            end if
            if (StressSFadjNEW > GetManagement_FertilityStress()) then
                StressSFadjNEW = GetManagement_FertilityStress()
            end if
        end if
    else
        if (GetManagement_FertilityStress() == 0 .or. &
            .not. GetCrop_StressResponse_Calibrated()) then
            ! no (calibrated) soil fertility stress
            StressSFadjNEW = 0
        else
            ! BiomassUnlim is too small
            StressSFadjNEW = GetManagement_FertilityStress()
        end if
    end if

    PreviousStressLevel = StressSFadjNEW
    SumKcTopStress = (1._dp - StressSFadjNEW/100._dp) * SumKcTop


    contains


    real(dp) function FractionFlowering(Dayi)
      integer(int32), intent(in) :: Dayi

      real(dp) :: f1, f2, F
      integer(int32) :: DiFlor

      if (GetCrop_LengthFlowering() <= 1) then
          F = 1._dp
      else
          DiFlor = dayi - (GetSimulation_DelayedDays() + &
                            GetCrop_Day1() + GetCrop_DaysToFlowering())
          f2 = FractionPeriod(DiFlor)
          DiFlor = (dayi-1) - (GetSimulation_DelayedDays() + &
                            GetCrop_Day1() + GetCrop_DaysToFlowering())
          f1 = FractionPeriod(DiFlor)
          if (abs(f1-f2) < ac_zero_threshold) then
              F = 0._dp
          else
              F = (100._dp * ((f1+f2)/2._dp)/GetCrop_LengthFlowering())
          end if
      end if
      FractionFlowering = F
    end function FractionFlowering


    real(dp) function FractionPeriod(DiFlor)
        integer(int32), intent(in) :: DiFlor

        real(dp) :: fi, TimePerc

        if (DiFlor <= epsilon(1._dp)) then
            fi = 0._dp
        else
            TimePerc = 100._dp * (DiFlor * 1._dp/GetCrop_LengthFlowering())
            if (TimePerc > 100._dp) then
                fi = 1._dp
            else
                fi = 0.00558_dp * exp(0.63_dp*log(TimePerc)) - &
                     0.000969_dp * TimePerc - 0.00383_dp
                if (fi < 0._dp) then
                    fi = 0._dp
                end if
            end if
        end if
        FractionPeriod = fi
    end function FractionPeriod


    integer(int32) function YearWeighingFactor(CropFirstDayNr)
        integer(int32), intent(in) :: CropFirstDayNr

        integer(int32) :: Dayi, Monthi, Yeari

        call DetermineDate(CropFirstDayNr, Dayi, Monthi, Yeari)
        YearWeighingFactor = Yeari
    end function YearWeighingFactor
end subroutine DetermineBiomassAndYield


subroutine AdjustpStomatalToETo(MeanETo, pStomatULAct)
    real(dp), intent(in) :: MeanETo
    real(dp), intent(inout) :: pStomatULAct


    select case (GetCrop_pMethod())
        case (pMethod_NoCorrection)
            pStomatULAct = GetCrop_pdef()

        case (pMethod_FAOCorrection)
             pStomatULAct = GetCrop_pdef() + GetSimulParam_pAdjFAO() * &
             (0.04_dp *(5._dp-MeanETo))*log10(10._dp-9._dp*GetCrop_pdef())
    end select
    if (pStomatULAct > 1) then
        pStomatULAct = 1._dp
    end if
    if (pStomatULAct < 0._dp) then
        pStomatULAct = 0._dp
    end if
end subroutine AdjustpStomatalToETo


subroutine AdjustpSenescenceToETo(EToMean, TimeSenescence, WithBeta, pSenAct)
    real(dp), intent(in) :: EToMean
    real(dp), intent(in) :: TimeSenescence
    logical, intent(in) :: WithBeta
    real(dp), intent(inout) :: pSenAct

    pSenAct = GetCrop_pSenescence()
    if (GetCrop_pMethod() == pMethod_FAOCorrection) then
        pSenAct = GetCrop_pSenescence() + GetSimulParam_pAdjFAO() &
                            * 0.04_dp*(5._dp-EToMean) &
                            * log10(10._dp-9._dp*GetCrop_pSenescence())
        if ((TimeSenescence > 0.0001_dp) .and. WithBeta) then
            pSenAct = pSenAct * (1._dp-GetSimulParam_Beta()/100._dp)
        end if
        if (pSenAct < 0._dp) then
            pSenAct = 0._dp
        end if
        if (pSenAct >= 1.0_dp) then
            pSenAct = 0.98_dp ! otherwise senescence is not possible at WP
        end if
    end if
end subroutine AdjustpSenescenceToETo


subroutine CheckGermination()

    real(dp) :: Zroot, WCGermination
    logical :: SWCtopSoilConsidered_temp

    ! total root zone is considered
    Zroot = GetCrop_RootMin()
    SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
    call DetermineRootZoneWC(Zroot, SWCtopSoilConsidered_temp)
    call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
    WCGermination = GetRootZoneWC_WP() + (GetRootZoneWC_FC() - &
                    GetRootZoneWC_WP()) * (GetSimulParam_TAWGermination()/100._dp)
    if (GetRootZoneWC_Actual() < WCGermination) then
        call SetSimulation_DelayedDays(GetSimulation_DelayedDays() + 1)
        call SetSimulation_SumGDD(0._dp)
    else
        call SetSimulation_Germinate(.true.)
        if (GetCrop_Planting() == plant_Seed) then
            call SetSimulation_ProtectedSeedling(.true.)
        else
            call SetSimulation_ProtectedSeedling(.false.)
        end if
    end if
end subroutine CheckGermination


subroutine calculate_transpiration(Tpot, Coeffb0Salt, Coeffb1Salt, Coeffb2Salt)
    real(dp), intent(in) :: Tpot
    real(dp), intent(in) :: Coeffb0Salt
    real(dp), intent(in) :: Coeffb1Salt
    real(dp), intent(in) :: Coeffb2Salt

    real(dp) :: WtoExtract, theta_critical, alfa, sinkMM
    integer(int32) :: compi, layeri, pre_layer
    real(dp) :: DeltaWC, InetThreshold
    real(dp) :: TpotMAX, RedFact, RedFactECsw
    real(dp) :: Wrel, WrelSalt, pStomatLLAct, crop_pActStom_tmp
    real(dp) :: CompiECe, CompiECsw, CompiECswFC
    logical :: SWCtopSoilConsidered_temp
    type(CompartmentIndividual), dimension(max_No_compartments) :: Comp_temp
    type(CompartmentIndividual) :: Compi_temp


    call SetTact(0.0_dp)

    if (Tpot > 0._dp) then
        ! 1. maximum transpiration in actual root zone
        if (GetIrriMode() == IrriMode_Inet) then
            ! salinity stress not considered
            TpotMAX = Tpot
        else
            SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
            call DetermineRootZoneWC(GetRootingDepth(), SWCtopSoilConsidered_temp)
            call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)

            ! --- 1. Effect of water stress and ECe (total rootzone)
            WrelSalt = (GetRootZoneWC_FC()-GetRootZoneWC_Actual())/ &
                       (GetRootZoneWC_FC()-GetRootZoneWC_WP())

            ! --- 2. Effect of water stress
            pStomatLLAct = 1._dp
            if (GetSimulation_SWCtopSoilConsidered() .eqv. .true.) then
                ! top soil is relative wetter than total root zone
                if (GetRootZoneWC_ZtopAct() < &
                                  (0.999_dp * GetRootZoneWC_ZtopThresh())) then
                    Wrel = (GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopAct())/ &
                           (GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopWP())
                    crop_pActStom_tmp = GetCrop_pActStom()
                    RedFact = (1._dp - GetSimulation_EffectStress_RedKsSto()/100._dp) &
                              * KsAny(Wrel, crop_pActStom_tmp, pStomatLLAct, 0.0_dp) ! where (0.0) is linear
                else
                    RedFact = (1._dp - GetSimulation_EffectStress_RedKsSto()/100._dp)
                end if
            else
                ! total root zone
                if (GetRootZoneWC_Actual() < (0.999_dp * GetRootZoneWC_Thresh())) then
                    Wrel = (GetRootZoneWC_FC()-GetRootZoneWC_Actual())/ &
                            (GetRootZoneWC_FC()-GetRootZoneWC_WP())
                    crop_pActStom_tmp = GetCrop_pActStom()
                    RedFact = (1._dp - GetSimulation_EffectStress_RedKsSto()/100._dp) &
                              * KsAny(Wrel, crop_pActStom_tmp, pStomatLLAct, 0.0_dp) ! where (0.0) is linear
                    call SetCrop_pActStom(crop_pActStom_tmp)
                else
                    RedFact = (1._dp - GetSimulation_EffectStress_RedKsSto()/100._dp)
                end if
            end if

            if (RedFact < 0._dp) then
                RedFact = 0._dp
            end if
            if (RedFact > 1._dp) then
                RedFact = 1._dp
            end if

            ! --- 3. Extra effect of ECsw (salt in total root zone is considered)
            if (GetSimulation_SalinityConsidered()) then
                RedFactECsw = AdjustedKsStoToECsw(GetCrop_ECemin(), &
                              GetCrop_ECemax(), GetCrop_ResponseECsw(), &
                              GetRootZoneSalt_ECe(), GetRootZoneSalt_ECsw(), &
                              GetRootZoneSalt_ECswFC(), WrelSalt, Coeffb0Salt, &
                              Coeffb1Salt, Coeffb2Salt, RedFact)
            else
                RedFactECsw = RedFact
            end if

            ! --- 4. Conclusion (adjustment of TpotMAX considering Water and Salt stress)
            TpotMAX = RedFactECsw * Tpot

            ! 1.b anaerobic conditions in root zone (total root zone is considered)
            call DetermineRootZoneAnaeroConditions(GetRootZoneWC_SAT(), &
                                              GetRootZoneWC_Actual(), &
                                              real(GetCrop_AnaeroPoint(), kind=dp), &
                                              GetRootingDepth(), RedFact)
            TpotMAX = RedFact * TpotMax
        end if

        ! 2. extraction of TpotMax out of the compartments
        ! 2.a initial settings
        Comp_temp = GetCompartment()
        call calculate_rootfraction_compartment(GetRootingDepth(), Comp_temp)
        call calculate_sink_values(TpotMAX, GetRootingDepth(), Comp_temp, GetCrop())
        call SetCompartment(Comp_temp)
        compi = 0
        pre_layer = 0
        loop: do
            compi = compi + 1
            layeri = GetCompartment_Layer(compi)
            if (layeri > pre_layer) then
                call calculate_theta_critical(layeri, theta_critical)
                pre_layer = layeri
            end if
            ! 2.b calculate alfa
            if (GetIrriMode() == IrriMode_Inet) then
                alfa = 1._dp
            else
                ! effect of water stress and ECe
                if (GetCompartment_theta(compi) >= theta_critical) then
                    alfa = (1._dp - GetSimulation_EffectStress_RedKsSto()/100._dp)
                elseif (GetCompartment_theta(compi) > (GetSoilLayer_WP(layeri)/100._dp)) then
                    if (theta_critical > (GetSoilLayer_WP(layeri)/100._dp)) then
                        Wrel = (GetSoilLayer_FC(layeri)/100._dp &
                        - GetCompartment_theta(compi)) &
                            /(GetSoilLayer_FC(layeri)/100._dp &
                                - GetSoilLayer_WP(layeri)/100._dp)
                        pStomatLLAct = 1._dp
                        crop_pActStom_tmp = GetCrop_pActStom()
                        alfa = (1._dp - GetSimulation_EffectStress_RedKsSto()/100._dp) &
                               * KsAny(Wrel, crop_pActStom_tmp, pStomatLLAct, &
                                                 GetCrop_KsShapeFactorStomata())
                        call SetCrop_pActStom(crop_pActStom_tmp)
                    else
                        alfa = (1._dp - GetSimulation_EffectStress_RedKsSto()/100._dp)
                    end if
                else
                    alfa = 0._dp
                end if
                ! extra effect of ECsw
                if (GetSimulation_SalinityConsidered()) then
                    WrelSalt = (GetSoilLayer_FC(layeri)/100._dp &
                                - GetCompartment_theta(compi)) &
                                    /(GetSoilLayer_FC(layeri)/100._dp &
                                        - GetSoilLayer_WP(layeri)/100._dp)
                    CompiECe = ECeComp(GetCompartment_i(compi))
                    CompiECsw = ECswComp(GetCompartment_i(compi), .false.)
                    CompiECswFC = ECswComp(GetCompartment_i(compi), .true.)
                    RedFactECsw = AdjustedKsStoToECsw(GetCrop_ECemin(), &
                                  GetCrop_ECemax(), GetCrop_ResponseECsw(), &
                                  CompiECe, CompiECsw, CompiECswFC, WrelSalt, &
                                  Coeffb0Salt, Coeffb1Salt, Coeffb2Salt, alfa)
                else
                    RedFactECsw = alfa
                end if
                alfa = RedFactECsw
            end if
            if (GetCrop_AnaeroPoint() > 0._dp) then
                Compi_temp = GetCompartment_i(compi)
                call Correction_Anaeroby(Compi_temp, alfa)
                call SetCompartment_i(compi, Compi_temp)
            end if
            ! 2.c extract water
            sinkMM = 1000._dp * (alfa * GetCompartment_WFactor(compi) * &
                    GetCompartment_Smax(compi)) * GetCompartment_Thickness(compi)
            WtoExtract = TpotMAX-GetTact()
            if (WtoExtract < sinkMM) then
                sinkMM = WtoExtract
            end if
            call SetCompartment_theta(compi, GetCompartment_theta(compi) &
                 - sinkMM/(1000._dp*GetCompartment_Thickness(compi)* &
                 (1._dp - GetSoilLayer_GravelVol(layeri)/100._dp)))
            WtoExtract = WtoExtract - sinkMM
            call SetTact(GetTact() + sinkMM)
            if ((WtoExtract < epsilon(1._dp) .or. &
                                    (compi == GetNrCompartments()))) exit loop
        end do loop


        ! 3. add net irrigation water requirement
        if (GetIrriMode() == IrriMode_Inet) then
            ! total root zone is considered
            SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
            call DetermineRootZoneWC(GetRootingDepth(), SWCtopSoilConsidered_temp)
            call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
            InetThreshold = GetRootZoneWC_FC() - GetSimulParam_PercRAW()/100._dp &
                            *(GetRootZoneWC_FC() - GetRootZoneWC_Thresh())
            if (GetRootZoneWC_Actual() < InetThreshold) then
                pre_layer = 0
                do compi = 1, GetNrCompartments()
                    layeri = GetCompartment_Layer(compi)
                    if (layeri > pre_layer) then
                        call calculate_theta_critical(layeri, theta_critical)
                        InetThreshold = GetSoilLayer_FC(layeri)/100._dp - &
                                        GetSimulParam_PercRAW()/100._dp* &
                                        (GetSoilLayer_FC(layeri)/100._dp - theta_critical)
                        pre_layer = layeri
                    end if
                    DeltaWC = GetCompartment_WFactor(compi) * &
                              (InetThreshold - GetCompartment_Theta(compi)) &
                              *1000._dp*GetCompartment_Thickness(compi)* &
                              (1._dp - GetSoilLayer_GravelVol(layeri)/100._dp)
                    call SetCompartment_Theta(compi, GetCompartment_theta(compi) &
                         + DeltaWC/(1000._dp*GetCompartment_Thickness(compi)* &
                        (1._dp - GetSoilLayer_GravelVol(layeri)/100._dp)))
                    call SetIrrigation(GetIrrigation() + DeltaWC)
                end do
            end if
        end if
    end if


    contains


    subroutine calculate_theta_critical(layeri, theta_critical)
        integer(int32), intent(in) :: layeri
        real(dp), intent(inout) :: theta_critical

        real(dp) :: theta_TAW

        theta_TAW = GetSoilLayer_FC(layeri)/100._dp &
                                             - GetSoilLayer_WP(layeri)/100._dp
        theta_critical = GetSoilLayer_FC(layeri)/100._dp - theta_TAW &
                                                           * GetCrop_pActStom()
    end subroutine calculate_theta_critical


    subroutine calculate_rootfraction_compartment(RootingDepth, Compartment)
        real(dp), intent(in) :: RootingDepth
        type(CompartmentIndividual), dimension(max_No_compartments), &
                                    intent(inout) :: Compartment

        real(dp) ::  frac_value, cumdepth
        integer(int32) :: compi, i

        cumdepth = 0._dp
        compi = 0
        loop: do
            compi = compi + 1
            cumdepth = cumdepth + Compartment(compi)%Thickness
            if (cumdepth <= RootingDepth) then
                Compartment(compi)%WFactor = 1
            else
                frac_value = RootingDepth - (cumdepth - Compartment(compi)%Thickness)
                if (frac_value > 0._dp) then
                    Compartment(compi)%WFactor = frac_value/Compartment(compi)%Thickness
                else
                    Compartment(compi)%WFactor = 0._dp
                end if
            end if
            if ((cumdepth >= RootingDepth) .or. &
                                      (compi == GetNrCompartments())) exit loop
        end do loop
        do i = compi+1, GetNrCompartments()
            Compartment(i)%WFactor = 0._dp
        end do
    end subroutine calculate_rootfraction_compartment


    subroutine calculate_sink_values(Tpot, RootingDepth, Compartment, Crop)
        real(dp), intent(in) :: Tpot
        real(dp), intent(in) :: RootingDepth
        type(CompartmentIndividual), dimension(max_No_compartments), &
                                    intent(inout) :: Compartment
        type(rep_crop), intent(in) :: Crop


        real(dp) ::       sink_value, StopComp, SbotComp, cumdepth
        integer(int32) :: compi, i

        if (GetIrriMode() == IrriMode_Inet) then
            sink_value = (GetCrop_SmaxTop() + GetCrop_SmaxBot())/2._dp
            do compi = 1, GetNrCompartments()
                Compartment(compi)%Smax = sink_value
            end do
        else
            cumdepth = 0._dp
            compi = 0
            SbotComp = GetCrop_SmaxTop()
            loop: do
                compi = compi + 1
                StopComp = SbotComp
                cumdepth = cumdepth + Compartment(compi)%Thickness
                if (cumdepth <= RootingDepth) then
                    SbotComp = GetCrop_SmaxBot() * GetSimulation_SCor() &
                               + (GetCrop_SmaxTop()- GetCrop_SmaxBot() &
                                    *GetSimulation_SCor()) &
                                       * (RootingDepth - cumdepth)/RootingDepth
                else
                    SbotComp = GetCrop_SmaxBot()*GetSimulation_SCor()
                end if
                Compartment(compi)%Smax = ((StopComp + SbotComp)/2._dp)
                if (Compartment(compi)%Smax > 0.06_dp) then
                    Compartment(compi)%Smax = 0.06_dp
                end if
                if ((cumdepth >= RootingDepth) .or. &
                                      (compi == GetNrCompartments())) exit loop
            end do loop
            do i = (compi + 1), GetNrCompartments()
                Compartment(i)%Smax = 0._dp
            end do
        end if
    end subroutine calculate_sink_values


    subroutine Correction_Anaeroby(Comp, alfa)
        type(CompartmentIndividual), intent(inout) :: Comp
        real(dp), intent(inout) :: alfa

        real(dp) :: alfaAN
        integer(int32) :: ini
        if ((GetDaySubmerged() >= GetSimulParam_DelayLowOxygen()) .and. &
                                         (GetCrop_AnaeroPoint() > 0._dp)) then
            alfaAN = 0._dp
        elseif (Comp%theta > (GetSoilLayer_SAT(Comp%Layer) &
                        - GetCrop_AnaeroPoint())/100._dp) then
            Comp%DayAnaero = Comp%DayAnaero + 1
            if (Comp%DayAnaero >= GetSimulParam_DelayLowOxygen()) then
                ini = 0
                Comp%DayAnaero = GetSimulParam_DelayLowOxygen()
            else
                ini = 1
            end if
            alfaAN = (GetSoilLayer_SAT(Comp%Layer)/100._dp &
                          - Comp%theta)/(GetCrop_AnaeroPoint()/100._dp)
            if (alfaAN < 0._dp) then
                alfaAN = 0._dp
            end if
            if (GetSimulParam_DelayLowOxygen() > 1._dp) then
                alfaAN = (ini+(Comp%DayAnaero-1._dp)*alfaAN) &
                            /(ini+Comp%DayAnaero-1._dp)
            end if
        else
            alfaAN = 1._dp
            Comp%DayAnaero = 0._dp
        end if
        if (alfa > alfaAN) then
            alfa = alfaAN
        end if
    end subroutine Correction_Anaeroby


    subroutine DetermineRootZoneAnaeroConditions(Wsat, Wact, AnaeVol, Zr, RedFact)
        real(dp), intent(in) :: Wsat
        real(dp), intent(in) :: Wact
        real(dp), intent(in) :: AnaeVol
        real(dp), intent(in) :: Zr
        real(dp), intent(inout) :: RedFact

        real(dp) :: SATVol, ACTVol
        RedFact = 1
        if ((AnaeVol > 0._dp) .and. (Zr > 0._dp)) then
            SATVol = Wsat/(10._dp*Zr)
            ACTVol = Wact/(10._dp*Zr)
            if (ACTVol > SATVol) then
                ACTVol = SATVol
            end if
            if (ActVol > (SatVol-AnaeVol)) then
                call SetSimulation_DayAnaero(GetSimulation_DayAnaero() + 1_int8)
                if (GetSimulation_DayAnaero() > GetSimulParam_DelayLowOxygen()) then
                    call SetSimulation_DayAnaero(int(GetSimulParam_DelayLowOxygen(), kind=int8))
                end if
                RedFact = 1._dp - (1._dp-((SATVol - ACTVol)/AnaeVol)) &
                                * (GetSimulation_DayAnaero()*1._dp &
                                    /GetSimulParam_DelayLowOxygen())
            else
                call SetSimulation_DayAnaero(0_int8)
            end if
        else
            call SetSimulation_DayAnaero(0_int8)
        end if
    end subroutine DetermineRootZoneAnaeroConditions
end subroutine calculate_transpiration


subroutine surface_transpiration(Coeffb0Salt, Coeffb1Salt, Coeffb2Salt)
    real(dp), intent(in) :: Coeffb0Salt
    real(dp), intent(in) :: Coeffb1Salt
    real(dp), intent(in) :: Coeffb2Salt

    real(dp) :: Part
    integer(int32) :: compi
    real(dp) :: KsReduction, SaltSurface
    real(dp) :: Tact_temp

    call SetDaySubmerged(GetDaySubmerged() + 1)
    do compi = 1, GetNrCompartments()
        call SetCompartment_DayAnaero(compi, GetCompartment_DayAnaero(compi) + 1)
        if (GetCompartment_DayAnaero(compi) > GetSimulParam_DelayLowOxygen()) then
            call SetCompartment_DayAnaero(compi, GetSimulParam_DelayLowOxygen())
        end if
    end do
    if (GetCrop_AnaeroPoint() > 0._dp) then
        Part = (1._dp-GetDaySubmerged()/real(GetSimulParam_DelayLowOxygen(), kind=dp))
    else
        Part = 1._dp
    end if
    KsReduction = KsSalinity(GetSimulation_SalinityConsidered(), GetCrop_ECemin(), &
                  GetCrop_ECemax(), GetECstorage(), 0.0_dp)
    SaltSurface = GetSurfaceStorage()*GetECstorage()*Equiv
    if (GetSurfaceStorage() > KsReduction*Part*GetTpot()) then
        call SetSurfaceStorage(GetSurfaceStorage() - KsReduction*Part*GetTpot())
        call SetTact(KsReduction*Part*GetTpot())
        ! salinisation of surface storage layer
        call SetECstorage(SaltSurface/(GetSurfaceStorage()*Equiv))
    else
        call SetTact(GetSurfaceStorage() -0.1_dp)
        call SetSurfaceStorage(0.1_dp) ! zero give error in already updated salt balance
    end if
    if (GetTact() < KsReduction*Part*GetTpot()) then
        Tact_temp = GetTact()   !(*Protect Tact from changes in the next routine*)
        call calculate_transpiration((KsReduction*Part*GetTpot()-GetTact()), &
                                         Coeffb0Salt, Coeffb1Salt, Coeffb2Salt)
        call SetTact(Tact_temp + GetTact())
    end if
end subroutine surface_transpiration


!-----------------------------------------------------------------------------
! BUDGET_module
!-----------------------------------------------------------------------------

real(dp) function calculate_delta_theta(theta_in, thetaAdjFC, NrLayer)
    real(dp), intent(in) :: theta_in
    real(dp), intent(in) :: thetaAdjFC
    integer(int32), intent(in) :: NrLayer

    real(dp) :: DeltaX, theta, theta_sat, theta_fc

    theta = theta_in
    theta_sat = GetSoilLayer_SAT(NrLayer) / 100.0_dp
    theta_fc = GetSoilLayer_FC(NrLayer) / 100.0_dp
    if (theta > theta_sat) then
        theta = theta_sat
    end if
    if (theta <= thetaAdjFC/100.0_dp) then
        DeltaX = 0.0_dp
    else
        DeltaX = GetSoilLayer_tau(NrLayer)&
                 * (theta_sat - theta_fc)&
                 * (exp(theta - theta_fc) - 1.0_dp)&
                 / (exp(theta_sat - theta_fc) - 1.0_dp)
        if ((theta - DeltaX) < thetaAdjFC) then
            DeltaX = theta - thetaAdjFC
        end if
    end if
    calculate_delta_theta = DeltaX
end function calculate_delta_theta


real(dp) function calculate_theta(delta_theta, thetaAdjFC, NrLayer)
    real(dp), intent(in) :: delta_theta
    real(dp), intent(in) :: thetaAdjFC
    integer(int32), intent(in) :: NrLayer

    real(dp) :: ThetaX, theta_sat, theta_fc, tau

    theta_sat = GetSoilLayer_SAT(NrLayer) / 100.0_dp
    theta_fc = GetSoilLayer_FC(NrLayer) / 100.0_dp
    tau = GetSoilLayer_tau(NrLayer)
    if (delta_theta <= epsilon(0.0_dp)) then
        ThetaX = thetaAdjFC
    elseif (tau > 0.0_dp) then
        ThetaX = theta_fc&
            + log(1.0_dp&
                  + delta_theta&
                  * (exp(theta_sat - theta_fc) - 1.0_dp)&
                  / (tau * (theta_sat - theta_fc)))
        if (ThetaX < thetaAdjFC) then
            ThetaX = thetaAdjFC
        end if
    else
        ! to stop draining
        ThetaX = theta_sat + 0.1_dp
    end if
    calculate_theta = ThetaX
end function calculate_theta


subroutine calculate_drainage()
    integer(int32) ::  i, compi, layeri, pre_nr
    real(dp) :: drainsum, delta_theta, drain_comp, drainmax, theta_x, excess
    real(dp) :: pre_thick
    logical :: drainability

    drainsum = 0.0_dp
    do compi=1, GetNrCompartments()
        ! 1. Calculate drainage of compartment
        ! ====================================
        layeri = GetCompartment_Layer(compi)
        if (GetCompartment_theta(compi) &
                 > GetCompartment_FCadj(compi)/100.0_dp) then
            delta_theta = calculate_delta_theta(GetCompartment_theta(compi), &
                (GetCompartment_FCadj(compi)/100.0_dp), layeri)
        else
            delta_theta = 0.0_dp
        end if
        drain_comp = delta_theta * 1000.0_dp * GetCompartment_Thickness(compi) &
                      * (1 - GetSoilLayer_GravelVol(layeri)/100.0_dp)


        ! 2. Check drainability
        ! =====================
        excess = 0.0_dp
        pre_thick = 0.0_dp
        do i = 1, (compi-1)
            pre_thick = pre_thick + GetCompartment_Thickness(i)
        end do
        drainmax = delta_theta * 1000.0_dp * pre_thick &
                   * (1 - GetSoilLayer_GravelVol(layeri)/100.0_dp)
        if (drainsum <= drainmax) then
            drainability = .true.
        else
            drainability = .false.
        end if

        ! 3. Drain compartment
        ! ====================
        if (drainability) then
            call SetCompartment_theta(compi, &
                     GetCompartment_theta(compi)-delta_theta)
            drainsum = drainsum + drain_comp
            call CheckDrainsum(layeri, drainsum, excess)
        else  ! drainability == .false.
            delta_theta = drainsum/(1000.0_dp * pre_thick&
                                    *(1-GetSoilLayer_GravelVol(layeri)/100.0_dp))
            theta_x = calculate_theta(delta_theta, &
                (GetCompartment_FCadj(compi)/100.0_dp), layeri)

            if (theta_x <= GetSoilLayer_SAT(layeri)/100.0_dp) then
                call SetCompartment_theta(compi, &
                         GetCompartment_theta(compi) &
                         + drainsum/(1000.0_dp*GetCompartment_Thickness(compi) &
                                     *(1-GetSoilLayer_GravelVol(layeri)/100.0_dp)))
                if (GetCompartment_theta(compi) > theta_x) then
                    drainsum = (GetCompartment_theta(compi) - theta_x) &
                               * 1000.0_dp * GetCompartment_Thickness(compi) &
                               * (1 - GetSoilLayer_GravelVol(layeri)/100.0_dp)
                    delta_theta = calculate_delta_theta(theta_x, &
                        (GetCompartment_FCadj(compi)/100.0_dp), layeri)
                    drainsum = drainsum +  delta_theta * 1000.0_dp &
                                           * GetCompartment_Thickness(compi) &
                                           * (1 - GetSoilLayer_GravelVol(layeri)&
                                                  /100.0_dp)
                    call CheckDrainsum(layeri, drainsum, excess)
                    call SetCompartment_theta(compi, theta_x - delta_theta)
                elseif (GetCompartment_theta(compi) &
                         > GetCompartment_FCadj(compi)/100.0_dp) then
                    delta_theta = calculate_delta_theta(&
                        GetCompartment_theta(compi), &
                        (GetCompartment_FCadj(compi)/100.0_dp), &
                        layeri)
                    call SetCompartment_theta(compi, &
                             GetCompartment_theta(compi) - delta_theta)
                    drainsum = delta_theta * 1000.0_dp &
                               * GetCompartment_Thickness(compi) &
                               * (1 - GetSoilLayer_GravelVol(layeri)/100.0_dp)
                    call CheckDrainsum(layeri, drainsum, excess)
                else
                    drainsum = 0.0_dp
                end if
            end if ! theta_x <= SoilLayer[layeri].SAT/100

            if (theta_x > GetSoilLayer_SAT(layeri)/100.0_dp) then
                call SetCompartment_theta(compi,&
                         GetCompartment_theta(compi)&
                         + drainsum/(1000.0_dp*GetCompartment_Thickness(compi) &
                                     *(1-GetSoilLayer_GravelVol(layeri)/100.0_dp)))
                if (GetCompartment_theta(compi) &
                         <= GetSoilLayer_SAT(layeri)/100.0_dp) then
                    if (GetCompartment_theta(compi) &
                            > GetCompartment_FCadj(compi)/100.0_dp) then
                        delta_theta = calculate_delta_theta(&
                            GetCompartment_theta(compi), &
                            (GetCompartment_FCadj(compi)/100.0_dp),&
                            layeri)
                        call SetCompartment_theta(compi, &
                                 GetCompartment_theta(compi) - delta_theta)
                        drainsum = delta_theta * 1000.0_dp &
                                   * GetCompartment_Thickness(compi) &
                                   *(1-GetSoilLayer_GravelVol(layeri)/100.0_dp)
                        call CheckDrainsum(layeri, drainsum, excess)
                    else
                        drainsum = 0.0_dp
                    end if
                end if
                if (GetCompartment_theta(compi)&
                        > GetSoilLayer_SAT(layeri)/100.0_dp) then
                    excess = (GetCompartment_theta(compi)&
                               - (GetSoilLayer_SAT(layeri)/100.0_dp)) &
                             * 1000.0_dp * GetCompartment_Thickness(compi) &
                             * (1 - GetSoilLayer_GravelVol(layeri)/100.0_dp)
                    delta_theta = calculate_delta_theta(&
                         GetCompartment_theta(compi), &
                         (GetCompartment_FCadj(compi)/100),&
                         layeri)
                    call SetCompartment_theta(compi, &
                             GetSoilLayer_SAT(layeri)/100.0_dp - delta_theta)
                    drain_comp = delta_theta * 1000.0_dp&
                                 * GetCompartment_Thickness(compi)&
                                 * (1 - GetSoilLayer_GravelVol(layeri)/100.0_dp)
                    drainmax = delta_theta * 1000.0_dp * pre_thick&
                               * (1 - GetSoilLayer_GravelVol(layeri)/100.0_dp)
                    if (drainmax > excess) then
                        drainmax = excess
                    end if
                    excess = excess - drainmax
                    drainsum = drainmax + drain_comp
                    call CheckDrainsum(layeri, drainsum, excess)
                end if
            end if ! theta_x > SoilLayer[layeri].SAT/100
        end if ! drainability = false

        call SetCompartment_fluxout(compi, drainsum)


        ! 4. Redistribute excess
        ! ======================
        if (excess > 0.0_dp) then
            pre_nr = compi + 1
            loop: do
                pre_nr = pre_nr - 1
                layeri = GetCompartment_Layer(pre_nr)
                if (pre_nr < compi) then
                    call SetCompartment_fluxout(pre_nr,&
                             GetCompartment_fluxout(pre_nr) - excess)
                end if
                call SetCompartment_theta(pre_nr,&
                    GetCompartment_theta(pre_nr)&
                    + excess&
                    / (1000.0_dp*GetCompartment_Thickness(pre_nr)&
                       *(1-GetSoilLayer_GravelVol(GetCompartment_Layer(pre_nr))&
                           /100.0_dp)))
                if (GetCompartment_theta(pre_nr) &
                        > GetSoilLayer_SAT(layeri)/100.0_dp) then
                    excess = (GetCompartment_theta(pre_nr) &
                              - GetSoilLayer_SAT(layeri)/100) &
                             * 1000.0_dp * GetCompartment_Thickness(pre_nr) &
                             * (1-GetSoilLayer_GravelVol(GetCompartment_Layer(pre_nr))&
                                    /100.0_dp)
                    call SetCompartment_theta(pre_nr,&
                             GetSoilLayer_SAT(layeri)/100.0_dp)
                else
                    excess = 0.0_dp
                end if
                if ((abs(excess) < epsilon(0._dp)) .or. (pre_nr == 1)) exit loop
            end do loop
            ! redistribute excess
        end if

    !Do-loop
    end do
    call SetDrain(drainsum)


    contains


    subroutine CheckDrainsum(layeri, drainsum, excess)
        integer(int32), intent(in) :: layeri
        real(dp), intent(inout) :: drainsum
        real(dp), intent(inout) :: excess

        if (drainsum > GetSoilLayer_InfRate(layeri)) then
            excess = excess + drainsum - GetSoilLayer_InfRate(layeri)
            drainsum = GetSoilLayer_InfRate(layeri)
        end if
    end subroutine CheckDrainsum
end subroutine calculate_drainage


subroutine calculate_weighting_factors(Depth, Compartment)
    real(dp), intent(in) :: Depth
    type(CompartmentIndividual), dimension(max_No_compartments), intent(inout) :: Compartment

    integer(int32) :: i, compi
    real(dp) :: CumDepth, xx, wx

    CumDepth = 0.0_dp
    xx = 0.0_dp
    compi = 0
    loop: do
        compi = compi + 1
        CumDepth = CumDepth + Compartment(compi)%Thickness
        if (CumDepth > Depth) then
            CumDepth = Depth
        end if
        wx = 1.016_dp * (1.0_dp - EXP(-4.16_dp * CumDepth/Depth))
        Compartment(compi)%WFactor = wx - xx
        if (Compartment(compi)%WFactor > 1.0_dp) then
            Compartment(compi)%WFactor = 1.0_dp
        end if
        if (Compartment(compi)%WFactor < 0.0_dp) then
            Compartment(compi)%WFactor = 0.0_dp
        end if
        xx = wx
        if ((CumDepth >= Depth) .or. (compi == GetNrCompartments())) exit loop
    enddo loop
    do i = (compi + 1), GetNrCompartments()
        Compartment(i)%WFactor = 0.0_dp
    end do
end subroutine calculate_weighting_factors


subroutine calculate_runoff(MaxDepth)
    real(dp), intent(in) :: MaxDepth

    real(dp) :: SUM, CNA, Shower, term, S
    integer(int8) :: CN2, CN1, CN3

    CN2 = roundc(GetSoil_CNvalue()&
                 * (100 + GetManagement_CNcorrection())/100.0_dp,&
                 mold=1_int8)
    if (GetRainRecord_DataType() == datatype_daily) then
        if (GetSimulParam_CNcorrection()) then
            call calculate_relative_wetness_topsoil(SUM)
            call DetermineCNIandIII(CN2, CN1, CN3)
            CNA = real(roundc(CN1+(CN3-CN1)*SUM, mold=1_int32), kind=dp)
        else
            CNA = real(CN2, kind=dp)
        end if
        Shower = GetRain()
    else
        CNA = real(CN2, kind=dp)
        Shower = (GetRain()*10.0_dp) &
                 / GetSimulParam_EffectiveRain_ShowersInDecade()
    end if
    S = 254.0_dp * (100.0_dp/CNA - 1.0_dp)
    term = Shower - (GetSimulParam_IniAbstract()/100.0_dp) * S
    if (term <= epsilon(0.0_dp)) then
        call SetRunoff(0.0_dp);
    else
        call SetRunoff(term**2&
             / (Shower + (1.0_dp - (GetSimulParam_IniAbstract()/100.0_dp)) * S))
    end if
    if ((GetRunoff() > 0.0_dp) .and. ( &
            (GetRainRecord_DataType() == datatype_decadely)&
             .or. (GetRainRecord_DataType() == datatype_monthly))) then
        if (GetRunoff() >= Shower) then
            call SetRunoff(GetRain())
        else
            call SetRunoff(GetRunoff() &
                   * (GetSimulParam_EffectiveRain_ShowersInDecade()/10.14_dp))
            if (GetRunoff() > GetRain()) then
                call SetRunoff(GetRain())
            end if
        end if
    end if


    contains


    subroutine calculate_relative_wetness_topsoil(SUM)
        real(dp), intent(inout) :: SUM

        real(dp) :: CumDepth, theta
        integer(int32) :: compi, layeri
        type(CompartmentIndividual), dimension(max_No_compartments) :: Compartment_temp

        Compartment_temp = GetCompartment()
        call calculate_weighting_factors(MaxDepth, Compartment_temp)
        call SetCompartment(Compartment_temp)
        SUM = 0.0_dp
        compi = 0
        CumDepth = 0.0_dp

        loop : do
            compi = compi + 1
            layeri = GetCompartment_Layer(compi)
            CumDepth = CumDepth + GetCompartment_Thickness(compi)
            if (GetCompartment_theta(compi) < GetSoilLayer_WP(layeri)/100.0_dp) then
                theta = GetSoilLayer_WP(layeri)/100.0_dp
            else
                theta = GetCompartment_theta(compi)
            end if
            SUM = SUM + GetCompartment_WFactor(compi) &
                 * (theta-GetSoilLayer_WP(layeri)/100.0_dp) &
                 / (GetSoilLayer_FC(layeri)/100.0_dp - GetSoilLayer_WP(layeri)/100.0_dp)
            if ((CumDepth >= MaxDepth) .or. (compi == GetNrCompartments())) exit loop
        end do loop

        if (SUM < 0.0_dp) then
            SUM = 0.0_dp
        end if
        if (SUM > 1.0_dp) then
            SUM = 1.0_dp
        end if
    end subroutine calculate_relative_wetness_topsoil
end subroutine calculate_runoff


subroutine Calculate_irrigation(SubDrain, TargetTimeVal, TargetDepthVal)
    real(dp), intent(inout) :: SubDrain
    integer(int32), intent(inout) :: TargetTimeVal
    integer(int32), intent(in) :: TargetDepthVal

    real(dp) :: ZrWC, RAWi
    logical :: SWCtopSoilConsidered_temp

    ! total root zone is considered
    SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
    call DetermineRootZoneWC(GetRootingDepth(), SWCtopSoilConsidered_temp)
    call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
    ZrWC = GetRootZoneWC_Actual() - GetEpot() - GetTpot() &
           + GetRain() - GetRunoff() - SubDrain
    if (GetGenerateTimeMode() == GenerateTimeMode_AllDepl) then
        if ((GetRootZoneWC_FC() - ZrWC) >= TargetTimeVal) then
            TargetTimeVal = 1
        else
            TargetTimeVal = 0
        end if
    end if
    if (GetGenerateTimeMode() == GenerateTimeMode_AllRAW) then
        RAWi = TargetTimeVal/100._dp &
                * (GetRootZoneWC_FC() - GetRootZoneWC_Thresh())
        if ((GetRootZoneWC_FC() - ZrWC) >= RAWi) then
            TargetTimeVal = 1
        else
            TargetTimeVal = 0
        end if
    end if
    if (TargetTimeVal == 1) then
        if (GetGenerateDepthMode() == GenerateDepthMode_FixDepth) then
            call SetIrrigation(real(TargetDepthVal, kind=dp))
        else
            call SetIrrigation((GetRootZoneWC_FC() - ZrWc) &
                                            + TargetDepthVal)
            if (GetIrrigation() < 0._dp) then
                call SetIrrigation(0._dp)
            end if
        end if
    else
        call SetIrrigation(0._dp)
    end if
end subroutine Calculate_irrigation


subroutine CalculateEffectiveRainfall(SubDrain)
    real(dp), intent(inout) :: SubDrain

    real(dp) :: EffecRain, ETcropMonth, RainMonth, &
                DrainMax, Zr, depthi, DTheta, RestTheta
    integer(int32) :: compi

    if (GetRain() > 0._dp) then
        ! 1. Effective Rainfall
        EffecRain = (GetRain()-GetRunoff())
        select case (GetSimulParam_EffectiveRain_Method())
        case (EffectiveRainMethod_Percentage)
            EffecRain = (GetSimulParam_EffectiveRain_PercentEffRain()/100._dp) &
                                * (GetRain()-GetRunoff())
        case (EffectiveRainMethod_USDA)
            ETcropMonth = ((GetEpot()+GetTpot())*30._dp)/25.4_dp ! inch/month
            RainMonth = ((GetRain()-GetRunoff())*30._dp)/25.4_dp ! inch/Month
            if (RainMonth > 0.1_dp) then
                EffecRain = (0.70917_dp*exp(0.82416_dp*log(RainMonth))-0.11556_dp) &
                                * (exp(0.02426_dp*ETcropMonth*log(10._dp)))
                                                                 ! inch/month
            else
                EffecRain = RainMonth
            end if
            EffecRain = EffecRain*(25.4_dp/30._dp) ! mm/day
        end select
    end if
    if (EffecRain < 0._dp) then
        EffecRain = 0._dp
    end if
    if (EffecRain > (GetRain()-GetRunoff())) then
        EffecRain = (GetRain()-GetRunoff())
    end if
    SubDrain = (GetRain()-GetRunoff()) - EffecRain

    ! 2. Verify Possibility of SubDrain
    if (SubDrain > 0._dp) then
        DrainMax = GetSoilLayer_InfRate(1)
        if (GetSurfaceStorage() > 0._dp) then
            DrainMax = 0._dp
        else
            Zr = GetRootingDepth()
            if (Zr <= epsilon(0._dp)) then
                Zr = (GetSimulParam_EvapZmax()/100._dp)
            end if
            compi = 0
            depthi = 0._dp
            DTheta = (EffecRain/Zr)/1000._dp
            loop: do
                compi = compi + 1
                depthi = depthi + GetCompartment_Thickness(compi)
                RestTheta = GetSoilLayer_SAT(GetCompartment_Layer(compi)) &
                            /100._dp - (GetCompartment_theta(compi) + DTheta)
                if (RestTheta <= epsilon(0._dp)) then
                    DrainMax = 0._dp
                end if
                if (GetSoilLayer_InfRate(GetCompartment_Layer(compi)) &
                                                        < DrainMax) then
                    DrainMax = GetSoilLayer_InfRate(GetCompartment_Layer(compi))
                end if
                if ((depthi >= Zr) &
                    .or. (compi >= GetNrCompartments())) exit loop
            end do loop
        end if
        if (SubDrain > DrainMax) then
            if (GetManagement_Bundheight() < 0.001_dp) then
                call SetRunoff(GetRunoff() + (SubDrain-DrainMax))
            end if
            SubDrain = DrainMax
        end if
    end if
end subroutine CalculateEffectiveRainfall


subroutine calculate_CapillaryRise(CRwater, CRsalt)
    real(dp), intent(inout) :: CRwater
    real(dp), intent(inout) :: CRsalt

    real(dp) :: Zbottom, MaxMM, DThetaMax, DTheta, LimitMM, &
                CRcomp, SaltCRi, DrivingForce, ZtopNextLayer, &
                Krel, ThetaThreshold
    integer(int32) :: compi, SCellAct, layeri

    Zbottom = 0._dp
    do compi = 1, GetNrCompartments()
        Zbottom = Zbottom + GetCompartment_Thickness(compi)
    end do

    ! start at the bottom of the soil profile
    compi = GetNrCompartments()
    MaxMM = MaxCRatDepth(GetSoilLayer_CRa(GetCompartment_Layer(compi)), &
                         GetSoilLayer_CRb(GetCompartment_Layer(compi)), &
                         GetSoilLayer_InfRate(GetCompartment_Layer(compi)), &
                         (Zbottom - GetCompartment_Thickness(compi)/2._dp), &
                         (GetZiAqua()/100._dp))

    ! check restrictions on CR from soil layers below
    ZtopNextLayer = 0._dp
    do layeri = 1, GetCompartment_Layer(GetNrCompartments())
        ZtopNextLayer = ZtopNextLayer + GetSoilLayer_Thickness(layeri)
    end do
    layeri = GetCompartment_Layer(GetNrCompartments())
    do while ((ZtopNextLayer < (GetZiAqua()/100._dp)) &
                .and. (layeri < GetSoil_NrSoilLayers()))
        layeri = layeri + 1
        LimitMM = MaxCRatDepth(GetSoilLayer_CRa(layeri), &
                               GetSoilLayer_CRb(layeri), &
                               GetSoilLayer_InfRate(layeri), &
                               ZtopNextLayer, (GetZiAqua()/100._dp))
        if (MaxMM > LimitMM) then
            MaxMM = LimitMM
        end if
        ZtopNextLayer = ZtopNextLayer + GetSoilLayer_Thickness(layeri)
    end do

    loop: do while ((roundc(MaxMM*1000._dp, mold=1) > 0) &
            .and. (compi > 0) &
            .and. (roundc(GetCompartment_fluxout(compi)*1000._dp, mold=1) == 0))
        ! Driving force
        if ((GetCompartment_theta(compi) &
                >= GetSoilLayer_WP(GetCompartment_Layer(compi))/100._dp) &
            .and. (GetSimulParam_RootNrDF() > 0_int8)) then
            DrivingForce = 1._dp &
                          - (exp(GetSimulParam_RootNrDF() &
                            * log(GetCompartment_theta(compi) &
                                - GetSoilLayer_WP(GetCompartment_Layer(compi)) &
                                                                    /100._dp)) &
                          /exp(GetSimulParam_RootNrDF() &
                            *log(GetCompartment_FCadj(compi)/100._dp &
                      - GetSoilLayer_WP(GetCompartment_Layer(compi))/100._dp)))
        else
            DrivingForce = 1._dp
        end if
        ! relative hydraulic conductivity
        ThetaThreshold = (GetSoilLayer_WP(GetCompartment_Layer(compi))/100._dp &
                          + GetSoilLayer_FC(GetCompartment_Layer(compi)) &
                                                                /100._dp)/2._dp
        if (GetCompartment_Theta(compi) < ThetaThreshold) then
            if ((GetCompartment_Theta(compi) &
                <= GetSoilLayer_WP(GetCompartment_Layer(compi))/100._dp) &
              .or. (ThetaThreshold &
                <= GetSoilLayer_WP(GetCompartment_Layer(compi))/100._dp)) then
                Krel = 0._dp
            else
                Krel = (GetCompartment_Theta(compi) &
                        - GetSoilLayer_WP(GetCompartment_Layer(compi))/100._dp) &
                      /(ThetaThreshold &
                        - GetSoilLayer_WP(GetCompartment_Layer(compi))/100._dp)
            end if
        else
            Krel = 1._dp
        end if

        ! room available to store water
        DTheta = GetCompartment_FCadj(compi)/100._dp &
                - GetCompartment_Theta(compi)
        if ((DTheta > 0._dp) &
            .and. ((Zbottom - GetCompartment_Thickness(compi)/2._dp) &
                    < (GetZiAqua()/100._dp))) then
            ! water stored
            DThetaMax = Krel * DrivingForce &
                        * MaxMM/(1000._dp*GetCompartment_Thickness(compi))
            if (DTheta >= DThetaMax) then
                call SetCompartment_Theta(compi, &
                                          GetCompartment_Theta(compi) &
                                                         + DThetaMax)
                CRcomp = DThetaMax*1000._dp*GetCompartment_Thickness(compi) &
                         * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
                MaxMM = 0._dp
            else
                call SetCompartment_Theta(compi, &
                                         GetCompartment_FCadj(compi)/100._dp)
                CRcomp = DTheta*1000._dp*GetCompartment_Thickness(compi) &
                         * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
                MaxMM = Krel * MaxMM - CRcomp
            end if
            CRwater = CRwater + CRcomp
            ! salt stored
            SCellAct = ActiveCells(GetCompartment_i(compi))
            SaltCRi = Equiv * CRcomp * GetECiAqua() ! gram/m2
            call SetCompartment_Salt(compi, SCellAct, &
                                     GetCompartment_Salt(compi, SCellAct) &
                                                               + SaltCRi)
            CRsalt = CRsalt + SaltCRi
        end if
        Zbottom = Zbottom - GetCompartment_Thickness(compi)
        compi = compi - 1
        if (compi < 1) exit loop
        LimitMM = MaxCRatDepth(&
                    GetSoilLayer_CRa(GetCompartment_Layer(compi)), &
                    GetSoilLayer_CRb(GetCompartment_Layer(compi)), &
                    GetSoilLayer_InfRate(GetCompartment_Layer(compi)), &
                    (Zbottom - GetCompartment_Thickness(compi)/2._dp), &
                    (GetZiAqua()/100._dp))
        if (MaxMM > LimitMM) then
            MaxMM = LimitMM
        end if
    end do loop
end subroutine calculate_CapillaryRise


subroutine CheckWaterSaltBalance(dayi,&
              InfiltratedRain,  &
              control, InfiltratedIrrigation,&
              InfiltratedStorage, Surf0, ECInfilt, ECdrain, &
              HorizontalWaterFlow, HorizontalSaltFlow, SubDrain)
    integer(int32), intent(in) :: dayi
    real(dp), intent(in) :: InfiltratedRain
    integer(intEnum), intent(in) :: control
    real(dp), intent(in) :: InfiltratedIrrigation
    real(dp), intent(in) :: InfiltratedStorage
    real(dp), intent(inout) :: Surf0
    real(dp), intent(inout) :: ECInfilt
    real(dp), intent(inout) :: ECdrain
    real(dp), intent(inout) :: HorizontalWaterFlow
    real(dp), intent(inout) :: HorizontalSaltFlow
    real(dp), intent(inout) :: SubDrain

    integer(int32) :: compi, layeri, celli
    real(dp) :: Surf1, ECw

    select case (control)
    case (control_begin_day)
        call SetTotalWaterContent_BeginDay(0._dp) ! mm
        Surf0 = GetSurfaceStorage() ! mm
        call SetTotalSaltContent_BeginDay(0._dp) ! Mg/ha
        do compi =1, GetNrCompartments()
            call SetTotalWaterContent_BeginDay(GetTotalWaterContent_BeginDay() &
               + GetCompartment_theta(compi)*1000._dp* &
                 GetCompartment_Thickness(compi) &
               * (1._dp - &
                  GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._dp))
            call SetCompartment_fluxout(compi, 0._dp)
            do celli = 1, GetSoilLayer_SCP1(GetCompartment_Layer(compi))
                call SetTotalSaltContent_BeginDay(&
                       GetTotalSaltContent_BeginDay() &
                       + (GetCompartment_Salt(compi, celli) + &
                          GetCompartment_Depo(compi, celli))/100._dp) ! Mg/ha
            end do
        end do
        call SetDrain(0._dp)
        call SetRunoff(0._dp)
        ! Eact is set to 0 at the beginning of the evaporation process
        call SetTact(0._dp)
        call SetInfiltrated(0._dp)
        ECinfilt = 0._dp
        SubDrain = 0._dp
        ECdrain = 0._dp
        HorizontalWaterFlow = 0._dp
        HorizontalSaltFlow = 0._dp
        call SetCRwater(0._dp)
        call SetCRsalt(0._dp)

    case (control_end_day)
        call SetInfiltrated(InfiltratedRain+InfiltratedIrrigation &
                            +InfiltratedStorage)
        do layeri = 1, GetSoil_NrSoilLayers()
            call SetSoilLayer_WaterContent(layeri, 0._dp)
        end do
        call SetTotalWaterContent_EndDay(0._dp)
        Surf1 = GetSurfaceStorage()
        call SetTotalSaltContent_EndDay(0._dp)

        ! quality of irrigation water
        if (dayi < GetCrop_Day1()) then
            ECw = GetIrriECw_PreSeason()
        else
            ECw = GetSimulation_IrriECw()
            if (dayi > GetCrop_DayN()) then
                ECw = GetIrriECw_PostSeason()
            end if
        end if

        do compi = 1, GetNrCompartments()
            call SetTotalWaterContent_EndDay(GetTotalWaterContent_EndDay() &
               + GetCompartment_theta(compi)*1000._dp*&
                 GetCompartment_Thickness(compi) &
               * (1._dp -&
                  GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._dp))
            call SetSoilLayer_WaterContent(GetCompartment_Layer(compi), &
                    GetSoilLayer_WaterContent(GetCompartment_Layer(compi)) &
                    + GetCompartment_theta(compi)*1000._dp*&
                          GetCompartment_theta(compi) &
                    * (1._dp - &
                       GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._dp))
            do celli = 1, GetSoilLayer_SCP1(GetCompartment_Layer(compi))
                call SetTotalSaltContent_EndDay(GetTotalSaltContent_EndDay() &
                   + (GetCompartment_Salt(compi, celli) + &
                      GetCompartment_Depo(compi, celli))/100._dp) ! Mg/ha
            end do
        end do
        call SetTotalWaterContent_ErrorDay(GetTotalWaterContent_BeginDay() &
                + Surf0 &
                - (GetTotalWaterContent_EndDay()+GetDrain()+GetRunoff()+GetEact()&
                + GetTact()+Surf1-GetRain()-GetIrrigation()-GetCRwater()-HorizontalWaterFlow))
        call SetTotalSaltContent_ErrorDay(GetTotalSaltContent_BeginDay() &
                - GetTotalSaltContent_EndDay() & ! Mg/ha
                + InfiltratedIrrigation*ECw*Equiv/100._dp &
                + InfiltratedStorage*ECinfilt*Equiv/100._dp &
                - GetDrain()*ECdrain*Equiv/100._dp &
                + GetCRsalt()/100._dp &
                + HorizontalSaltFlow)
        call SetSumWaBal_Epot(GetSumWaBal_Epot() + GetEpot())
        call SetSumWaBal_Tpot(GetSumWaBal_Tpot() + GetTpot())
        call SetSumWaBal_Rain(GetSumWaBal_Rain() + GetRain())
        call SetSumWaBal_Irrigation(GetSumWaBal_Irrigation() + GetIrrigation())
        call SetSumWaBal_Infiltrated(GetSumWaBal_Infiltrated() + &
                  GetInfiltrated())
        call SetSumWaBal_Runoff(GetSumWaBal_Runoff() + GetRunoff())
        call SetSumWaBal_Drain(GetSumWaBal_Drain() + GetDrain())
        call SetSumWaBal_Eact(GetSumWaBal_Eact() + GetEact())
        call SetSumWaBal_Tact(GetSumWaBal_Tact() + GetTact())
        call SetSumWaBal_TrW(GetSumWaBal_TrW() + GetTactWeedInfested())
        call SetSumWaBal_CRwater(GetSumWaBal_CRwater() + GetCRwater())

        if (((dayi-GetSimulation_DelayedDays()) >= GetCrop_Day1() ) &
            .and. ((dayi-GetSimulation_DelayedDays()) <= GetCrop_DayN())) then
            ! in growing cycle
            if (GetSumWaBal_Biomass() > 0._dp) then
                ! biomass was already produced (i.e. CC present)
                ! and still canopy cover
                if (GetCCiActual() > 0._dp) then
                    call SetSumWaBal_ECropCycle(GetSumWaBal_ECropCycle() &
                           + GetEact())
                end if
            else
                call SetSumWaBal_ECropCycle(GetSumWaBal_ECropCycle() &
                           + GetEact()) ! before germination
            end if
        end if
        call SetSumWaBal_CRsalt(GetSumWaBal_CRsalt() + GetCRsalt()/100._dp)
        call SetSumWaBal_SaltIn(GetSumWaBal_SaltIn() + &
               (InfiltratedIrrigation*ECw+InfiltratedStorage*ECinfilt)*Equiv/100._dp)
        call SetSumWaBal_SaltOut(GetSumWaBal_SaltOut() + &
                GetDrain()*ECdrain*Equiv/100._dp)
    end select
end subroutine CheckWaterSaltBalance


subroutine calculate_saltcontent(InfiltratedRain, InfiltratedIrrigation, &
                                 InfiltratedStorage, SubDrain, dayi)
    real(dp), intent(in) :: InfiltratedRain
    real(dp), intent(in) :: InfiltratedIrrigation
    real(dp), intent(in) :: InfiltratedStorage
    integer(int32), intent(in) :: dayi
    real(dp), intent(in) :: SubDrain

    real(dp) ::   SaltIN, SaltOUT, mmIN, DeltaTheta, Theta, SAT, &
                  mm1, mm2, Dx, limit, Dif, UL
    real(dp) :: Zr, depthi, ECsubdrain, ECcel, DeltaZ, ECsw1, ECsw2, &
                ECsw, SM1, SM2, DS1, DS2, DS

    integer(int32) :: compi, celi, celiM1, Ni
    real(dp) :: ECw
    real(dp) :: Salt_temp, Salt2_temp, Depo_temp, Depo2_temp
    type(CompartmentIndividual) :: Compi_temp


    mmIN = InfiltratedRain + InfiltratedIrrigation + InfiltratedStorage

    ! quality of irrigation water
    if (dayi < GetCrop_Day1()) then
        ECw = GetIrriECw_PreSeason()
    else
        ECw = GetSimulation_IrriECw()
        if (dayi > GetCrop_DayN()) then
            ECw = GetIrriECw_PostSeason()
        end if
    end if

    ! initialise salt balance
    SaltIN = InfiltratedIrrigation*ECw*Equiv &
            + InfiltratedStorage*GetECstorage()*Equiv
    call SetSaltInfiltr(SaltIN/100._dp)
                ! salt infiltrated in soil profile kg/ha
    SaltOut= 0._dp


    do compi = 1, GetNrCompartments()
        ! 0. Set compartment parameters
        SAT = (GetSoilLayer_SAT(GetCompartment_Layer(compi)))/100._dp  ! m3/m3
        UL = GetSoilLayer_UL(GetCompartment_Layer(compi)) ! m3/m3
                                      ! Upper limit of SC salt cel
        Dx = GetSoilLayer_Dx(GetCompartment_Layer(compi)) ! m3/m3
                            ! Size of salts cel (expect last one)

        ! 1. Initial situation before drain and infiltration
        DeltaTheta = mmIN &
           /(1000._dp*GetCompartment_Thickness(compi) &
               *(1._dp - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                  /100._dp))
        Theta = GetCompartment_theta(compi) - DeltaTheta &
                + GetCompartment_fluxout(compi) &
                        /(1000._dp*GetCompartment_Thickness(compi))

        ! 2. Determine active SaltCels and Add IN
        Theta = Theta + DeltaTheta
        if (Theta <= UL) then
            celi = 0
            do while (Theta > Dx*celi)
                celi = celi + 1
            end do
        else
            celi = GetSoilLayer_SCP1(GetCompartment_Layer(compi))
        end if
        if (celi == 0) then
            celi = 1  ! XXX would be best to avoid celi=0 to begin with
        end if
        if (DeltaTheta > 0._dp) then
            call SetCompartment_Salt(compi, celi, &
                                     GetCompartment_Salt(compi, celi) &
                                        + SaltIN)
        end if

        ! 3. Mixing
        if (celi > 1) then
            do Ni = 1, (celi-1)
                mm1 = Dx*1000._dp*GetCompartment_Thickness(compi) &
                        * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
                if (Ni < GetSoilLayer_SC(GetCompartment_Layer(compi))) then
                    mm2 = mm1
                elseif (Theta > SAT) then
                    mm2 = (Theta-UL)*1000._dp*GetCompartment_Thickness(compi) &
                            * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
                else
                    mm2 = (SAT-UL)*1000._dp*GetCompartment_Thickness(compi) &
                            * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
                end if
                Dif = GetSoilLayer_SaltMobility_i(GetCompartment_Layer(compi), &
                                                    Ni)
                Salt_temp = GetCompartment_Salt(compi, Ni)
                Salt2_temp = GetCompartment_Salt(compi, Ni+1)
                Depo_temp = GetCompartment_Depo(compi, Ni)
                Depo2_temp = GetCompartment_Depo(compi, Ni+1)
                call Mixing(Dif, mm1, mm2, Salt_temp, Salt2_temp, &
                                           Depo_temp, Depo2_temp)
                call SetCompartment_Salt(compi, Ni, Salt_temp)
                call SetCompartment_Salt(compi, Ni+1, Salt2_temp)
                call SetCompartment_Depo(compi, Ni, Depo_temp)
                call SetCompartment_Depo(compi, Ni+1, Depo2_temp)
            end do
        end if

        ! 4. Drain
        SaltOut = 0._dp
        if (GetCompartment_fluxout(compi) > 0._dp) then
            DeltaTheta = GetCompartment_fluxout(compi) &
                         /(1000._dp*GetCompartment_Thickness(compi) &
                            * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                     /100._dp))
            do while (DeltaTheta > 0._dp)
                if (celi < GetSoilLayer_SCP1(GetCompartment_Layer(compi))) then
                    limit = (celi-1._dp)*Dx
                else
                    limit = UL
                end if
                if ((Theta - DeltaTheta) < limit) then
                    SaltOut = SaltOut + GetCompartment_Salt(compi, celi) &
                              + GetCompartment_Depo(compi, celi)
                    call SetCompartment_Salt(compi, celi, 0._dp)
                    mm1 = (Theta - limit)*1000._dp &
                           * GetCompartment_Thickness(compi) * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
                    if (SaltOut > (GetSimulParam_SaltSolub() * mm1)) then
                        call SetCompartment_Depo(compi, celi, &
                                   SaltOut - (GetSimulParam_SaltSolub() * mm1))
                        SaltOut = (GetSimulParam_SaltSolub() * mm1)
                    else
                        call SetCompartment_Depo(compi, celi, 0._dp)
                    end if
                    DeltaTheta = DeltaTheta - (Theta-limit)
                    Theta = limit
                    celi = celi - 1
                else
                    SaltOut = SaltOut &
                              + (GetCompartment_Salt(compi, celi) &
                                + GetCompartment_Depo(compi, celi)) &
                              * (DeltaTheta/(Theta-limit))
                    call SetCompartment_Salt(compi, celi, &
                                             GetCompartment_Salt(compi, celi) &
                                             * (1._dp-DeltaTheta/(Theta-limit)))
                    call SetCompartment_Depo(compi, celi, &
                                             GetCompartment_Depo(compi, celi) &
                                             * (1._dp-DeltaTheta/(Theta-limit)))
                    mm1 = DeltaTheta*1000._dp*GetCompartment_Thickness(compi) &
                            * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
                    if (SaltOut > (GetSimulParam_SaltSolub() * mm1)) then
                        call SetCompartment_Depo(&
                                compi, celi, &
                                GetCompartment_Depo(compi, celi) &
                                    + (SaltOut - GetSimulParam_SaltSolub() &
                                                                    * mm1))
                        SaltOut = (GetSimulParam_SaltSolub() * mm1)
                    end if
                    DeltaTheta = 0._dp
                    mm1 = GetSoilLayer_Dx(GetCompartment_Layer(compi)) &
                            *1000._dp*GetCompartment_Thickness(compi) &
                            * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
                    if (celi == GetSoilLayer_SCP1(GetCompartment_Layer(compi))) then
                        mm1 = 2._dp*mm1
                    end if
                    Salt_temp = GetCompartment_Salt(compi, celi)
                    Depo_temp = GetCompartment_Depo(compi, celi)
                    call SaltSolutionDeposit(mm1, Salt_temp, Depo_temp)
                    call SetCompartment_Salt(compi, celi, Salt_temp)
                    call SetCompartment_Depo(compi, celi, Depo_temp)
                end if
            end do
        end if

        mmIN = GetCompartment_fluxout(compi)
        SaltIN = SaltOUT
    end do

    if (GetDrain() > 0.001_dp) then
        call SetECdrain(SaltOUT/(GetDrain()*Equiv))
    end if

    ! 5. vertical salt diffusion
    celi = ActiveCells(GetCompartment_i(1))
    SM2 = GetSoilLayer_SaltMobility_i(GetCompartment_Layer(1), celi)/4._dp
    ECsw2 = ECswComp(GetCompartment_i(1), .false.) ! not at FC
    mm2 = GetCompartment_Theta(1)*1000._dp*GetCompartment_Thickness(1) &
            * (1._dp - GetSoilLayer_GravelVol(GetCompartment_Layer(1))/100._dp)
    do compi = 2, GetNrCompartments()
        celiM1 = celi
        SM1 = SM2
        ECsw1 = ECsw2
        mm1 = mm2
        celi =  ActiveCells(GetCompartment_i(compi))
        SM2 = GetSoilLayer_SaltMobility_i(GetCompartment_Layer(compi), &
                                          celi)/4._dp
        ECsw2 = ECswComp(GetCompartment_i(compi), .false.) ! not at FC
        mm2 = GetCompartment_Theta(compi)*1000._dp &
                * GetCompartment_Thickness(compi) &
                * (1._dp &
                    - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                /100._dp)
        ECsw = (ECsw1*mm1+ECsw2*mm2)/(mm1+mm2)
        DS1 = (ECsw1 - (ECsw1+(ECsw-ECsw1)*SM1))*mm1*Equiv
        DS2 = (ECsw2 - (ECsw2+(ECsw-ECsw2)*SM2))*mm2*Equiv
        if (abs(DS2) < abs(DS1)) then
            DS = abs(DS2)
        else
            DS = abs(DS1)
        end if
        if (DS > 0._dp) then
            if (ECsw1 > ECsw) then
                DS = DS*(-1._dp)
            end if
            Compi_temp = GetCompartment_i(compi-1)
            call MoveSaltTo(Compi_temp, celiM1, DS)
            call SetCompartment_i(compi-1, Compi_temp)
            DS = DS*(-1._dp)
            Compi_temp = GetCompartment_i(compi)
            call MoveSaltTo(Compi_temp, celi, DS)
            call SetCompartment_i(compi, Compi_temp)
        end if
    end do

    ! 6. Internal salt movement as a result of SubDrain
    ! SubDrain part of non-effective rainfall (10-day & monthly input)
    if (SubDrain > 0._dp) then
        Zr = GetRootingDepth()
        if (Zr >= epsilon(0._dp)) then
            Zr = (GetSimulParam_EvapZmax()/100._dp) ! in meter
        end if
        compi = 0
        depthi = 0._dp
        ECsubdrain = 0._dp

        ! extract
        loop: do
            compi = compi + 1
            depthi = depthi + GetCompartment_Thickness(compi)
            if (depthi <= Zr) then
                DeltaZ = GetCompartment_Thickness(compi)
            else
                DeltaZ = GetCompartment_Thickness(compi) - (depthi-Zr)
            end if
            celi = ActiveCells(GetCompartment_i(compi))
            if (celi < GetSoilLayer_SCP1(GetCompartment_Layer(compi))) then
                mm1 = GetSoilLayer_Dx(GetCompartment_Layer(compi))*1000._dp &
                        * GetCompartment_Thickness(compi) &
                        * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
            else
                mm1 = 2._dp*GetSoilLayer_Dx(GetCompartment_Layer(compi)) &
                        *1000._dp*GetCompartment_Thickness(compi) &
                        * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
            end if
            ECcel = GetCompartment_Salt(compi, celi)/(mm1*Equiv)
            ECsubdrain = (ECcel*mm1*(DeltaZ/GetCompartment_Thickness(compi)) &
                                                       +ECsubdrain*SubDrain) &
                         /(mm1*(DeltaZ/GetCompartment_Thickness(compi)) &
                                                              +SubDrain)
            call SetCompartment_Salt(&
                    compi, celi, &
                    (1._dp - (DeltaZ/GetCompartment_Thickness(compi))) &
                        * GetCompartment_Salt(compi, celi) &
                    + (DeltaZ/GetCompartment_Thickness(compi)) &
                        *ECsubdrain*mm1*Equiv)
            Salt_temp = GetCompartment_Salt(compi, celi)
            Depo_temp = GetCompartment_Depo(compi, celi)
            call SaltSolutionDeposit(mm1, Salt_temp, Depo_temp)
            call SetCompartment_Salt(compi, celi, Salt_temp)
            call SetCompartment_Depo(compi, celi, Depo_temp)
            if ((depthi >= Zr) .or. (compi >= GetNrCompartments())) exit loop
        end do loop

        ! dump
        if (compi >= GetNrCompartments()) then
            SaltOUT = GetECdrain()*(GetDrain()*Equiv) + ECsubdrain*SubDrain*Equiv
            call SetECdrain(SaltOUT/(GetDrain()*Equiv))
        else
            compi = compi + 1
            celi = ActiveCells(GetCompartment_i(compi))
            if (celi < GetSoilLayer_SCP1(GetCompartment_Layer(compi))) then
                mm1 = GetSoilLayer_Dx(GetCompartment_Layer(compi))*1000._dp &
                        * GetCompartment_Thickness(compi) &
                        * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
            else
                mm1 = 2._dp*GetSoilLayer_Dx(GetCompartment_Layer(compi)) &
                        *1000._dp*GetCompartment_Thickness(compi) &
                        * (1._dp &
                          - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                      /100._dp)
            end if
            call SetCompartment_Salt(compi, celi, &
                                     GetCompartment_Salt(compi, celi) &
                                         + ECsubdrain*SubDrain*Equiv)
            Salt_temp = GetCompartment_Salt(compi, celi)
            Depo_temp = GetCompartment_Depo(compi, celi)
            call SaltSolutionDeposit(mm1, Salt_temp, Depo_temp)
            call SetCompartment_Salt(compi, celi, Salt_temp)
            call SetCompartment_Depo(compi, celi, Depo_temp)
        end if
    end if


    contains


    subroutine Mixing(Dif, mm1, mm2, Salt1, Salt2, Depo1, Depo2)
        real(dp), intent(in) :: Dif
        real(dp), intent(in) :: mm1
        real(dp), intent(in) :: mm2
        real(dp), intent(inout) :: Salt1
        real(dp), intent(inout) :: Salt2
        real(dp), intent(inout) :: Depo1
        real(dp), intent(inout) :: Depo2

        real(dp) :: EC1, EC2, ECmix

        call SaltSolutionDeposit(mm1, Salt1, Depo1)
        EC1 = Salt1/(mm1*Equiv)
        call SaltSolutionDeposit(mm2, Salt2, Depo2)
        EC2 = Salt2/(mm2*Equiv)
        ECmix = (EC1*mm1+EC2*mm2)/(mm1+mm2)
        EC1 = EC1 + (ECmix-EC1)*Dif
        EC2 = EC2 + (ECmix-EC2)*Dif
        Salt1 = EC1*mm1*Equiv
        call SaltSolutionDeposit(mm1, Salt1, Depo1)
        Salt2 = EC2*mm2*Equiv
        call SaltSolutionDeposit(mm2, Salt2, Depo2)
    end subroutine Mixing


    subroutine MoveSaltTo(Compx, celx, DS)
        type(CompartmentIndividual), intent(inout) :: Compx
        integer(int32), intent(in) :: celx
        real(dp), intent(in) :: DS

        real(dp) :: mmx
        integer(int32) :: celx_local

        celx_local = celx
        if (DS >= epsilon(0._dp)) then
            Compx%Salt(celx_local) = Compx%Salt(celx_local) + DS
            mmx = GetSoilLayer_Dx(Compx%Layer)*1000._dp*Compx%Thickness &
                    * (1._dp - GetSoilLayer_GravelVol(Compx%Layer)/100._dp)
            if (celx_local == GetSoilLayer_SCP1(Compx%Layer)) then
                mmx = 2._dp*mmx
            end if
            call SaltSolutionDeposit(mmx, Compx%Salt(celx_local), &
                                     Compx%Depo(celx_local))
        else
            celx_local = GetSoilLayer_SCP1(Compx%Layer)
            Compx%Salt(celx_local) = Compx%Salt(celx_local) + DS
            mmx = 2._dp*GetSoilLayer_Dx(Compx%Layer)*1000._dp*Compx%Thickness &
                    * (1._dp - GetSoilLayer_GravelVol(Compx%Layer)/100._dp)
            call SaltSolutionDeposit(mmx, Compx%Salt(celx_local), &
                                     Compx%Depo(celx_local))
            mmx = mmx/2._dp
            do while (Compx%Salt(celx_local) < 0._dp)
                ! index zero problem
                ! TO DO: likely also happened with original Pascal code, 
                !        but Pascal code tolerates it
                if (celx_local == 1_int32) then 
                    celx_local = size(Compx%Salt)
                end if
                Compx%Salt(celx_local-1) = Compx%Salt(celx_local-1) &
                                           + Compx%Salt(celx_local)
                Compx%Salt(celx_local) = 0._dp
                celx_local = celx_local - 1
                call SaltSolutionDeposit(mmx, Compx%Salt(celx_local), &
                                         Compx%Depo(celx_local))
            end do
        end if
    end subroutine MoveSaltTo
end subroutine calculate_saltcontent


subroutine calculate_Extra_runoff(InfiltratedRain, InfiltratedIrrigation, &
                                  InfiltratedStorage, SubDrain)
    real(dp), intent(inout) :: InfiltratedRain
    real(dp), intent(inout) :: InfiltratedIrrigation
    real(dp), intent(inout) :: InfiltratedStorage
    real(dp), intent(inout) :: SubDrain

    real(dp) :: FracSubDrain

    InfiltratedStorage = 0._dp
    InfiltratedRain = GetRain() - GetRunoff()
    if (InfiltratedRain > 0._dp) then
        FracSubDrain = SubDrain/InfiltratedRain
    else
        FracSubDrain = 0._dp
    end if
    if ((GetIrrigation()+InfiltratedRain) &
            > GetSoilLayer_InfRate(GetCompartment_Layer(1))) then
        if (GetIrrigation() > GetSoilLayer_InfRate(GetCompartment_Layer(1))) then
            InfiltratedIrrigation = GetSoilLayer_InfRate(GetCompartment_Layer(1))
            call SetRunoff(GetRain() + (GetIrrigation()-InfiltratedIrrigation))
            InfiltratedRain = 0._dp
            SubDrain = 0._dp
        else
            InfiltratedIrrigation = GetIrrigation()
            InfiltratedRain = GetSoilLayer_InfRate(GetCompartment_Layer(1)) &
                                - InfiltratedIrrigation
            SubDrain = FracSubDrain*InfiltratedRain
            call SetRunoff(GetRain() - InfiltratedRain)
        end if
    else
        InfiltratedIrrigation = GetIrrigation()
    end if
end subroutine calculate_Extra_runoff


subroutine calculate_surfacestorage(InfiltratedRain, InfiltratedIrrigation, &
                                    InfiltratedStorage, ECinfilt, SubDrain, &
                                    dayi)
    real(dp), intent(inout) :: InfiltratedRain
    real(dp), intent(inout) :: InfiltratedIrrigation
    real(dp), intent(inout) :: InfiltratedStorage
    real(dp), intent(inout) :: ECinfilt
    real(dp), intent(in) :: SubDrain
    integer(int32), intent(in) :: dayi

    real(dp) :: Sum
    real(dp) :: ECw

    InfiltratedRain = 0._dp
    InfiltratedIrrigation = 0._dp
    if (GetRainRecord_DataType() == datatype_Daily) then
        Sum = GetSurfaceStorage() + GetIrrigation() + GetRain()
    else
        Sum = GetSurfaceStorage() + GetIrrigation() + GetRain() &
              - GetRunoff() - SubDrain
    end if
    if (Sum > 0._dp) then
        ! quality of irrigation water
        if (dayi < GetCrop_Day1()) then
            ECw = GetIrriECw_PreSeason()
        else
            ECw = GetSimulation_IrriECw()
            if (dayi > GetCrop_DayN()) then
                ECw = GetIrriECw_PostSeason()
            end if
        end if
        ! quality of stored surface water
        call SetECstorage((GetECstorage() &
                            * GetSurfaceStorage() &
                            + ECw*GetIrrigation()) /Sum)
        ! quality of infiltrated water (rain and/or irrigation and/or stored surface water)
        ECinfilt = GetECstorage()
        ! surface storage
        if (Sum > GetSoilLayer_InfRate(GetCompartment_Layer(1))) then
            InfiltratedStorage = GetSoilLayer_InfRate(GetCompartment_Layer(1))
            call SetSurfaceStorage(Sum - InfiltratedStorage)
        else
            if (GetRainRecord_DataType() == datatype_Daily) then
                InfiltratedStorage = Sum
            else
                InfiltratedStorage = GetSurfaceStorage() + GetIrrigation()
                InfiltratedRain = GetRain() - GetRunoff()
            end if
            call SetSurfaceStorage(0._dp)
        end if
        ! extra run-off
        if (GetSurfaceStorage() > (GetManagement_BundHeight()*1000._dp)) then
            call SetRunoff(GetRunoff() &
                            + (GetSurfaceStorage() &
                            - GetManagement_BundHeight()*1000._dp))
            call SetSurfaceStorage(GetManagement_BundHeight()*1000._dp)
        end if
    else
        InfiltratedStorage = 0._dp
        call SetECstorage(0._dp)
    end if
end subroutine calculate_surfacestorage


subroutine calculate_infiltration(InfiltratedRain, InfiltratedIrrigation, &
                                  InfiltratedStorage, SubDrain)
    real(dp), intent(inout) :: InfiltratedRain
    real(dp), intent(inout) :: InfiltratedIrrigation
    real(dp), intent(inout) :: InfiltratedStorage
    real(dp), intent(inout) :: SubDrain

    integer(int32) :: compi, layeri, pre_comp
    real(dp) :: RunoffIni, amount_still_to_store, factor, &
                delta_theta_nul, delta_theta_SAT, theta_nul, &
                drain_max, diff, excess
    real(dp) :: EffecRain, Zr, depthi, DeltaZ, StorableMM


    ! calculate_infiltration
    ! A -  INFILTRATION versus STORAGE in Rootzone (= EffecRain)
    if (GetRainRecord_DataType() == datatype_Daily) then
        amount_still_to_store = InfiltratedRain + InfiltratedIrrigation &
                                + InfiltratedStorage
        EffecRain = 0._dp
    else
        amount_still_to_store = InfiltratedIrrigation + InfiltratedStorage
        EffecRain = InfiltratedRain - SubDrain
    end if

    ! B - INFILTRATION through TOP soil surface
    if (amount_still_to_store > 0._dp) then
        RunoffIni = GetRunoff()
        compi = 0

        loop1: do
            compi = compi + 1
            layeri = GetCompartment_Layer(compi)

            !1. Calculate multiplication factor
            !====================================
            factor = calculate_factor(layeri, compi)

            !2. Calculate theta nul
            !========================
            delta_theta_nul = amount_still_to_store &
                             /(1000._dp * GetCompartment_Thickness(compi) &
                              * (1._dp-GetSoilLayer_GravelVol(layeri)/100._dp))
            delta_theta_SAT = calculate_delta_theta(&
                                GetSoilLayer_SAT(layeri)/100._dp, &
                                GetSoilLayer_FC(layeri)/100._dp, &
                                layeri)

            if (delta_theta_nul < delta_theta_SAT) then
                theta_nul = calculate_theta(delta_theta_nul, &
                                            GetSoilLayer_FC(layeri)/100._dp, &
                                            layeri)
                if (theta_nul <= (GetCompartment_FCadj(compi)/100._dp)) then
                    theta_nul = GetCompartment_FCadj(compi)/100._dp
                    delta_theta_nul = calculate_delta_theta( &
                                        theta_nul, &
                                        GetSoilLayer_FC(layeri)/100._dp, &
                                        layeri)
                end if
                if (theta_nul > GetSoilLayer_SAT(layeri)/100._dp) then
                    theta_nul = GetSoilLayer_SAT(layeri)/100._dp
                end if
            else
                theta_nul = GetSoilLayer_SAT(layeri)/100._dp
                delta_theta_nul = delta_theta_SAT
            end if

            !3. Calculate drain max
            !========================
            drain_max = factor * delta_theta_nul * 1000._dp &
                            * GetCompartment_Thickness(compi) &
                            * (1._dp-GetSoilLayer_GravelVol(layeri)/100._dp)
            if ((GetCompartment_fluxout(compi) + drain_max) &
                        > GetSoilLayer_InfRate(layeri)) then
                drain_max = GetSoilLayer_InfRate(layeri) &
                            - GetCompartment_fluxout(compi)
            end if

            !4. Store water
            !================
            diff = theta_nul - GetCompartment_theta(compi)
            if (diff > 0._dp) then
                call SetCompartment_theta(compi, GetCompartment_theta(compi) &
                           + amount_still_to_store &
                             /(1000._dp * GetCompartment_Thickness(compi) &
                                 * (1._dp &
                                    - GetSoilLayer_GravelVol(layeri)/100._dp)))
                if (GetCompartment_theta(compi) > theta_nul) then
                    amount_still_to_store = (GetCompartment_theta(compi) &
                                                - theta_nul) &
                               * 1000._dp &
                               * GetCompartment_Thickness(compi) &
                               * (1._dp-GetSoilLayer_GravelVol(layeri)/100._dp)
                    call SetCompartment_theta(compi, theta_nul)
                else
                    amount_still_to_store = 0.0_dp
                end if
            end if
            call SetCompartment_fluxout(compi, GetCompartment_fluxout(compi) &
                                               + amount_still_to_store)

            !5. Redistribute excess
            !========================
            excess = amount_still_to_store - drain_max
            if (excess < 0._dp) then
                excess = 0._dp
            end if
            amount_still_to_store = amount_still_to_store - excess

            if (excess > 0._dp) then
                pre_comp = compi + 1
                loop2: do
                    pre_comp = pre_comp - 1
                    layeri = GetCompartment_Layer(pre_comp)
                    call SetCompartment_fluxout(pre_comp, &
                                    GetCompartment_fluxout(pre_comp) &
                                                           - excess)
                    call SetCompartment_theta(&
                            pre_comp, GetCompartment_theta(pre_comp) &
                               + excess/(1000._dp &
                                  * GetCompartment_Thickness(pre_comp) &
                                  * (1._dp &
                        - GetSoilLayer_GravelVol(GetCompartment_Layer(pre_comp))&
                                                                    /100._dp)))
                    if (GetCompartment_theta(pre_comp) &
                            > GetSoilLayer_SAT(layeri)/100._dp) then
                        excess = (GetCompartment_theta(pre_comp) &
                            - GetSoilLayer_SAT(layeri)/100._dp) * 1000._dp &
                                * GetCompartment_Thickness(pre_comp) &
                                * (1._dp &
                        -GetSoilLayer_GravelVol(GetCompartment_Layer(pre_comp))&
                                                                      /100._dp)
                        call SetCompartment_theta( &
                                pre_comp, &
                                GetSoilLayer_SAT(layeri)/100._dp)
                    else
                        excess = 0.0_dp
                    end if
                    if ((excess < epsilon(0._dp)) .or. (pre_comp == 1)) exit loop2
                end do loop2
                if (excess > 0._dp) then
                    call SetRunoff(GetRunoff() + excess)
                end if
            end if

            if ((amount_still_to_store <= epsilon(0._dp)) &
                    .or. (compi == GetNrCompartments())) exit loop1
        end do loop1
        if (amount_still_to_store > 0._dp) then
            call SetDrain(GetDrain() + amount_still_to_store)
        end if

        !6. Adjust infiltrated water
        !=============================
        if (GetRunoff() > RunoffIni) then
            if (GetManagement_Bundheight() >= 0.01_dp) then
                call SetSurfaceStorage(GetSurfaceStorage() &
                                        + (GetRunoff() &
                                        - RunoffIni))
                InfiltratedStorage = InfiltratedStorage &
                                     - (GetRunoff()-RunoffIni)
                if (GetSurfaceStorage() &
                            > GetManagement_BundHeight()*1000._dp) then
                    call SetRunoff(RunoffIni &
                                   + (GetSurfaceStorage() &
                                        - GetManagement_BundHeight()*1000._dp))
                    call SetSurfaceStorage(GetManagement_BundHeight()*1000._dp)
                else
                    call SetRunoff(RunoffIni)
                end if
            else
                InfiltratedRain = InfiltratedRain - (GetRunoff()-RunoffIni)
                if (InfiltratedRain < 0._dp) then
                    InfiltratedIrrigation = InfiltratedIrrigation &
                                            + InfiltratedRain
                    InfiltratedRain = 0._dp
                end if
            end if

            ! INFILTRATION through TOP soil surface
        end if
    end if

    ! C - STORAGE in Subsoil (= SubDrain)
    if (SubDrain > 0._dp) then
        amount_still_to_store = SubDrain

        ! Where to store
        Zr = GetRootingDepth()
        if (Zr <= 0._dp) then
            Zr = GetSimulParam_EvapZmax()/100._dp
        end if
        compi = 0
        depthi = 0._dp
        loop3: do
            compi = compi + 1
            depthi = depthi + GetCompartment_Thickness(compi)
            if ((depthi >= Zr) &
                .or. (compi >= GetNrCompartments())) exit loop3
        end do loop3
        if (depthi > Zr) then
            DeltaZ = (depthi - Zr)
        else
            DeltaZ = 0._dp
        end if

        ! Store
        do while((amount_still_to_store > 0._dp) &
                .and. ((compi < GetNrCompartments()) &
                    .or. (DeltaZ > 0._dp)))
            if (abs(DeltaZ) < epsilon(0._dp)) then
                compi = compi + 1
                DeltaZ = GetCompartment_Thickness(compi)
            end if
            StorableMM = (GetSoilLayer_SAT(GetCompartment_Layer(compi))&
                                                            /100._dp &
                            - GetCompartment_Theta(compi)) * 1000._dp &
                                 * DeltaZ * (1._dp &
                     - GetSoilLayer_GravelVol(GetCompartment_Layer(compi))&
                                                                 /100._dp)
            if (StorableMM > amount_still_to_store) then
               call SetCompartment_theta(&
                      compi, &
                      GetCompartment_Theta(compi) &
                       + (amount_still_to_store)&
                          /(1000._dp*GetCompartment_Thickness(compi) &
                            * (1._dp &
                    - GetSoilLayer_GravelVol(GetCompartment_Layer(compi))&
                                                                /100._dp)))
                amount_still_to_store = 0._dp
            else
                amount_still_to_store = amount_still_to_store - StorableMM
                call SetCompartment_theta(&
                        compi, &
                        GetCompartment_Theta(compi) &
                        + (StorableMM)/(1000._dp &
                            * GetCompartment_Thickness(compi) &
                            * (1._dp &
                     - GetSoilLayer_GravelVol(GetCompartment_Layer(compi))&
                                                                /100._dp)))
            end if
            DeltaZ = 0._dp
            if (amount_still_to_store &
                  > GetSoilLayer_InfRate(GetCompartment_Layer(compi))) then
                SubDrain = SubDrain &
                            - (amount_still_to_store &
                        - GetSoilLayer_InfRate(GetCompartment_Layer(compi)))
                EffecRain = EffecRain &
                            + (amount_still_to_store &
                        - GetSoilLayer_InfRate(GetCompartment_Layer(compi)))
                amount_still_to_store = GetSoilLayer_InfRate(&
                                            GetCompartment_Layer(compi))
            end if
        end do

        ! excess
        if (amount_still_to_store > 0._dp) then
            call SetDrain(GetDrain() + amount_still_to_store)
        end if
        ! STORAGE in Subsoil (= SubDrain)
    end if

    ! D - STORAGE in Rootzone (= EffecRain)
    if (EffecRain > 0._dp) then
        Zr = GetRootingDepth()
        if (Zr <= epsilon(0._dp)) then
            Zr = GetSimulParam_EvapZmax()/100._dp
        end if
        amount_still_to_store = EffecRain

        ! Store
        ! step 1 fill to FC (from top to bottom)
        compi = 0
        depthi = 0._dp
        loop4: do
            compi = compi + 1
            depthi = depthi + GetCompartment_Thickness(compi)
            if (depthi <= Zr) then
                DeltaZ = GetCompartment_Thickness(compi)
            else
                DeltaZ = GetCompartment_Thickness(compi) &
                         - (depthi-Zr)
            end if
            StorableMM = (GetCompartment_FCadj(compi)/100._dp &
                            - GetCompartment_Theta(compi))*1000._dp*DeltaZ &
                            * (1._dp &
                      - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                  /100._dp)
            if (StorableMM < 0._dp) then
                StorableMM = 0._dp
            end if
            if (StorableMM > amount_still_to_store) then
                call SetCompartment_theta(&
                        compi, &
                        GetCompartment_Theta(compi) &
                        + amount_still_to_store &
                            /(1000._dp*GetCompartment_Thickness(compi) &
                              *(1._dp &
                       - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                 /100._dp)))
                amount_still_to_store = 0._dp
            elseif (StorableMM > 0._dp) then
                call SetCompartment_theta(&
                        compi, &
                        GetCompartment_Theta(compi) &
                        + StorableMM &
                            /(1000._dp*GetCompartment_Thickness(compi) &
                                * (1._dp &
                    - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                              /100._dp)))
                amount_still_to_store = amount_still_to_store - StorableMM
            end if
            if ((depthi >= Zr) &
                    .or. (compi >= GetNrCompartments()) &
                    .or. (amount_still_to_store <= epsilon(0._dp))) &
                            exit loop4
        end do loop4

        ! step 2 fill to SATURATION (from bottom to top)
        if (amount_still_to_store > 0._dp) then
            loop5: do
                if (depthi > Zr) then
                    DeltaZ = GetCompartment_Thickness(compi) - (depthi-Zr)
                else
                    DeltaZ = GetCompartment_Thickness(compi)
                end if
                StorableMM = (GetSoilLayer_SAT(GetCompartment_Layer(compi)) &
                                                                   /100._dp &
                             - GetCompartment_Theta(compi))*1000._dp*DeltaZ &
                                * (1._dp &
                      - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                  /100._dp)
                if (StorableMM < 0._dp) then
                    StorableMM = 0._dp
                end if
                if (StorableMM > amount_still_to_store) then
                    call SetCompartment_theta(&
                            compi, &
                            GetCompartment_theta(compi) &
                            + amount_still_to_store &
                                /(1000._dp*GetCompartment_Thickness(compi) &
                                    * (1._dp &
                      - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                /100._dp)))
                    amount_still_to_store = 0._dp
                elseif (StorableMM > 0._dp) then
                    call SetCompartment_theta(&
                            compi, &
                            GetCompartment_Theta(compi) &
                            + StorableMM &
                                /(1000._dp*GetCompartment_Thickness(compi) &
                                    *(1._dp &
                      - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                /100._dp)))
                    amount_still_to_store = amount_still_to_store &
                                           - StorableMM
                end if
                compi = compi - 1
                depthi = depthi - GetCompartment_Thickness(compi)
                if ((compi == 0) &
                    .or. (amount_still_to_store <= epsilon(0._dp))) &
                            exit loop5
            end do loop5
        end if

        ! excess
        if (amount_still_to_store > 0._dp) then
            if (InfiltratedRain > 0._dp) then
                InfiltratedRain = InfiltratedRain - amount_still_to_store
            end if
            if (GetManagement_Bundheight() >= 0.01_dp) then
                call SetSurfaceStorage(GetSurfaceStorage() &
                                      + amount_still_to_store)
                if (GetSurfaceStorage() &
                        > (GetManagement_BundHeight()*1000._dp)) then
                    call SetRunoff(GetRunoff() &
                                   + (GetSurfaceStorage() &
                                   - GetManagement_BundHeight()*1000._dp))
                    call SetSurfaceStorage(GetManagement_BundHeight() &
                                                            *1000._dp)
                end if
            else
                call SetRunoff(GetRunoff() + amount_still_to_store)
            end if
        end if
        ! STORAGE in Rootzone (= EffecRain)
    end if


    contains


    real(dp) function Calculate_factor(layeri, compi)
        integer(int32), intent(in) :: layeri
        integer(int32), intent(in) :: compi

        real(dp) :: delta_theta_SAT

        delta_theta_SAT = calculate_delta_theta(GetSoilLayer_SAT(layeri)/100._dp, &
                                                GetSoilLayer_FC(layeri)/100._dp, &
                                                layeri)
        if (delta_theta_SAT > 0._dp) then
            Calculate_factor = GetSoilLayer_InfRate(layeri)&
                                /(delta_theta_SAT * 1000._dp &
                                    * GetCompartment_Thickness(compi) &
                                    * (1._dp-GetSoilLayer_GravelVol(layeri) &
                                                                    /100._dp))
        else
            Calculate_factor = 1._dp
        end if
    end function Calculate_factor
end subroutine calculate_infiltration


subroutine DetermineCCiGDD(CCxTotal, CCoTotal, &
                           StressLeaf, FracAssim, MobilizationON, &
                           StorageON, SumGDDAdjCC, VirtualTimeCC, &
                           StressSenescence, TimeSenescence, NoMoreCrop, &
                           CDCTotal, GDDayFraction, &
                           GDDayi, GDDCDCTotal, GDDTadj)
    real(dp), intent(in) :: CCxTotal
    real(dp), intent(in) :: CCoTotal
    real(dp), intent(inout) :: StressLeaf
    real(dp), intent(in) :: FracAssim
    logical, intent(in) :: MobilizationON
    logical, intent(in) :: StorageON
    real(dp), intent(in) :: SumGDDAdjCC
    integer(int32), intent(in) :: VirtualTimeCC
    real(dp), intent(inout) :: StressSenescence
    real(dp), intent(inout) :: TimeSenescence
    logical, intent(inout) :: NoMoreCrop
    real(dp), intent(in) :: CDCTotal
    real(dp), intent(in) :: GDDayFraction
    real(dp), intent(in) :: GDDayi
    real(dp), intent(in) :: GDDCDCTotal
    integer(int32), intent(in) :: GDDTadj

    real(dp), parameter :: CCdormant = 0.05_dp

    real(dp) :: pLeafLLAct , GDDCGCadjusted, GDDCDCadjusted, &
                CCiSen, GDDtTemp, CCxSF, CGCGDDSF, CCxSFCD, &
                RatDGDD, KsRED, CCibis
    integer(int32) :: GDDtFinalCCx
    logical :: WithBeta
    logical :: TheSenescenceON

    !! test Version 6.2
    real(dp) :: KsSen
    real(dp) :: Crop_pLeafAct_temp
    real(dp) :: Crop_pSenAct_temp
    real(dp) :: Crop_CCxAdjusted_temp


    if ((SumGDDadjCC <= GetCrop_GDDaysToGermination()) &
          .or. (roundc(SumGDDadjCC, mold=1) > GetCrop_GDDaysToHarvest())) then
        call SetCCiActual(0._dp)
    else
        ! growing season (once germinated)
        ! 1. find some parameters
        CGCGDDSF = GetCrop_GDDCGC() &
                    * (1._dp - GetSimulation_EffectStress_RedCGC()/100._dp)
        GDDCGCadjusted = CGCGDDSF

        RatDGDD = 1._dp
        if (GetCrop_GDDaysToFullCanopySF() < GetCrop_GDDaysToSenescence()) then
            RatDGDD = (GetCrop_DaysToSenescence() &
                        - GetCrop_DaysToFullCanopySF()) &
                      /real(GetCrop_GDDaysToSenescence() &
                        - GetCrop_GDDaysToFullCanopySF(), kind=dp)
        end if

        CCxSF = CCxTotal*(1._dp - GetSimulation_EffectStress_RedCCX()/100._dp)
        ! maximum canopy cover than can be reached
        ! (considering soil fertility/salinity, weed stress)
        if (SumGDDadjCC <= GetCrop_GDDaysToFullCanopySF()) then
            CCxSFCD = CCxSF ! no canopy decline before max canopy can be reached
        else
            ! canopy decline due to soil fertility
            if (SumGDDadjCC < GetCrop_GDDaysToSenescence()) then
                CCxSFCD = CCiNoWaterStressSF(&
                            (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                            GetCrop_DaysToGermination(), &
                            GetCrop_DaysToFullCanopySF(), &
                            GetCrop_DaysToSenescence(), &
                            GetCrop_DaysToHarvest(), &
                            GetCrop_GDDaysToGermination(), &
                            GetCrop_GDDaysToFullCanopySF(), &
                            GetCrop_GDDaysToSenescence(), &
                            GetCrop_GDDaysToHarvest(), &
                            CCoTotal, CCxTotal, GetCrop_CGC(), &
                            GetCrop_GDDCGC(), CDCTotal, GDDCDCTotal, &
                            SumGDDadjCC, RatDGDD, &
                            GetSimulation_EffectStress_RedCGC(), &
                            GetSimulation_EffectStress_RedCCX(), &
                            GetSimulation_EffectStress_CDecline(), &
                            GetCrop_ModeCycle())
            else
                CCxSFCD = CCxSF &
                          - (RatDGDD &
                                * GetSimulation_EffectStress_CDecline()/100._dp) &
                          * (GetCrop_GDDaysToSenescence() &
                                - GetCrop_GDDaysToFullCanopySF())
            end if
            if (CCxSFCD < 0._dp) then
                CCxSFCD = 0._dp
            end if
        end if
        StressLeaf = undef_int
        if ((abs(SumGDDadjCC - GetCrop_GDDaysToGermination()) < epsilon(0._dp)) &
                .and. (GetCrop_DaysToCCini() == 0)) then
            call SetCCiPrev(CCoTotal)
        end if

        ! time of potential vegetative growth
        GDDtFinalCCx = GetCrop_GDDaysToSenescence() ! non determinant crop
        if ((GetCrop_subkind() == subkind_Grain) &
                .and. (GetCrop_DeterminancyLinked())) then
            ! determinancy
            ! reduce GDDtFinalCCx in f(determinancy of crop)
            if (GetCrop_DaysToCCini() /= 0) then
                ! regrowth
                GDDtFinalCCx = GetCrop_GDDaysToFullCanopy() &
                               + roundc(GDDayFraction &
                                        * (GetCrop_GDDaysToFlowering() &
                                            + (GetCrop_GDDLengthFlowering() &
                                                                    /2._dp) &
                                            + GDDTadj &
                                            + GetCrop_GDDaysToGermination() &
                                            - GetCrop_GDDaysToFullCanopy()), &
                                mold=1) ! slow down
            else
                ! sown or transplant
                GDDtFinalCCx = GetCrop_GDDaysToFlowering() &
                               + roundc(GetCrop_GDDLengthFlowering()/2._dp, &
                                                                    mold=1)
            end if
            if (GDDtFinalCCx > GetCrop_GDDaysToSenescence()) then
                GDDtFinalCCx = GetCrop_GDDaysToSenescence()
            end if
        end if

        ! Crop.pLeafAct and Crop.pSenAct for plotting root zone depletion in RUN
        Crop_pLeafAct_temp = GetCrop_pLeafAct()
        call AdjustpLeafToETo(GetETo(), Crop_pLeafAct_temp, pLeafLLAct)
        call SetCrop_pLeafAct(Crop_pLeafAct_temp)
        WithBeta = .true.
        Crop_pSenAct_temp = GetCrop_pSenAct()
        call AdjustpSenescenceToETo(GetETo(), TimeSenescence, WithBeta, Crop_pSenAct_temp)
        call SetCrop_pSenAct(Crop_pSenAct_temp)

        ! 2. Canopy can still develop (stretched to GDDtFinalCCx)
        if (SumGDDadjCC < GDDtFinalCCx) then
            ! Canopy can stil develop (stretched to GDDtFinalCCx)
            if ((GetCCiPrev() <= GetCrop_CCoAdjusted()) &
                .or. (SumGDDadjCC <= GDDayi) &
                .or. ((GetSimulation_ProtectedSeedling()) &
                        .and. (GetCCiPrev() <= (1.25_dp * CCoTotal)))) then
                ! 2.a First day or very small CC as a result of senescence
                ! (no adjustment for leaf stress)
                if (GetSimulation_ProtectedSeedling()) then
                    call SetCCiActual(CanopyCoverNoStressSF(&
                                (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                                 GetCrop_DaysToGermination(), &
                                 GetCrop_DaysToSenescence(), &
                                 GetCrop_DaysToHarvest(), &
                                 GetCrop_GDDaysToGermination(), &
                                 GetCrop_GDDaysToSenescence(), &
                                 GetCrop_GDDaysToHarvest(), &
                                 CCoTotal, CCxTotal, GetCrop_CGC(), &
                                 CDCTotal, GetCrop_GDDCGC(), &
                                 GDDCDCadjusted, SumGDDadjCC, &
                                 GetCrop_ModeCycle(), &
                                 GetSimulation_EffectStress_RedCGC(), &
                                 GetSimulation_EffectStress_RedCCX()))
                    if (GetCCiActual() > (1.25_dp * CCoTotal)) then
                        call SetSimulation_ProtectedSeedling(.false.)
                    end if
                else
                    call SetCCiActual(GetCrop_CCoAdjusted() * exp(CGCGDDSF * GDDayi))
                end if
            ! 2.b CC > CCo
            else
                if (GetCCiPrev() < (0.97999_dp*CCxSF)) then
                    call DetermineGDDCGCadjusted(GDDCGCadjusted)
                    if (GDDCGCadjusted > ac_zero_threshold) then
                        ! Crop.GDDCGC or GDDCGCadjusted > 0
                        Crop_CCxAdjusted_temp = GetCrop_CCxAdjusted()
                        call DetermineCCxAdjusted(Crop_CCxAdjusted_temp)
                        call SetCrop_CCxAdjusted(Crop_CCxAdjusted_temp)
                        if (GetCrop_CCxAdjusted() < 0._dp) then
                            call SetCCiActual(GetCCiPrev())
                        elseif (abs(GetCCiPrev() - 0.97999_dp*CCxSF) &
                                    < 0.001_dp) then
                            call SetCCiActual(CanopyCoverNoStressSF(&
                                (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                                GetCrop_DaysToGermination(), &
                                GetCrop_DaysToSenescence(), &
                                GetCrop_DaysToHarvest(), &
                                GetCrop_GDDaysToGermination(), &
                                GetCrop_GDDaysToSenescence(), &
                                GetCrop_GDDaysToHarvest(), &
                                CCoTotal, CCxTotal, GetCrop_CGC(), &
                                CDCTotal, GetCrop_GDDCGC(), GDDCDCadjusted, &
                                SumGDDadjCC, GetCrop_ModeCycle(), &
                                GetSimulation_EffectStress_RedCGC(), &
                                GetSimulation_EffectStress_RedCCX()))
                        else
                            GDDtTemp = RequiredGDD(GetCCiprev(), &
                                                   GetCrop_CCoAdjusted(), &
                                                   GetCrop_CCxAdjusted(), &
                                                   GDDCGCadjusted)
                            if (GDDtTemp < 0._dp) then
                                call SetCCiActual(GetCCiPrev())
                            else
                                GDDtTemp = GDDtTemp + GDDayi
                                call SetCCiActual(CCatGDDTime(GDDtTemp, &
                                                        GetCrop_CCoAdjusted(), &
                                                        GDDCGCadjusted, &
                                                        GetCrop_CCxAdjusted()))
                            end if
                        end if
                    else
                        ! GDDCGCadjusted = 0 - too dry for leaf expansion
                        call SetCCiActual(GetCCiPrev())
                        if (GetCCiActual() > GetCrop_CCoAdjusted()) then
                            call SetCrop_CCoAdjusted(CCoTotal)
                        else
                            call SetCrop_CCoAdjusted(GetCCiActual())
                        end if
                    end if
                else
                    call SetCCiActual(CanopyCoverNoStressSF(&
                            (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                            GetCrop_DaysToGermination(), &
                            GetCrop_DaysToSenescence(), &
                            GetCrop_DaysToHarvest(), &
                            GetCrop_GDDaysToGermination(), &
                            GetCrop_GDDaysToSenescence(), &
                            GetCrop_GDDaysToHarvest(), &
                            CCoTotal, CCxTotal, GetCrop_CGC(), CDCTotal, &
                            GetCrop_GDDCGC(), GDDCDCadjusted, SumGDDadjCC, &
                            GetCrop_ModeCycle(), &
                            GetSimulation_EffectStress_RedCGC(), &
                            GetSimulation_EffectStress_RedCCX()))
                    call SetCrop_CCoAdjusted(CCoTotal)
                    StressLeaf = -33._dp ! maximum canopy is reached;
                end if
                if (GetCCiActual() > CCxSFCD) then
                    call SetCCiActual(CCxSFCD)
                    StressLeaf = -33._dp ! maximum canopy is reached;
                end if
            end if
            call SetCrop_CCxAdjusted(GetCCiActual())

        ! 3. Canopy can no longer develop
        ! (Mid-season (from tFinalCCx) or Late season stage)
        else
            StressLeaf = -33._dp ! maximum canopy is reached;
            if (GetCrop_CCxAdjusted() < 0._dp) then
                call SetCrop_CCxAdjusted(GetCCiPrev())
            end if

            if (SumGDDadjCC < GetCrop_GDDaysToSenescence()) then ! mid-season
                if (GetCrop_CCxAdjusted() > 0.97999_dp*CCxSF) then
                    call SetCCiActual(CanopyCoverNoStressSF(&
                                (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                                GetCrop_DaysToGermination(), &
                                GetCrop_DaysToSenescence(), &
                                GetCrop_DaysToHarvest(), &
                                GetCrop_GDDaysToGermination(), &
                                GetCrop_GDDaysToSenescence(), &
                                GetCrop_GDDaysToHarvest(), &
                                CCoTotal, CCxTotal, GetCrop_CGC(), &
                                CDCTotal, GetCrop_GDDCGC(), &
                                GDDCDCadjusted, SumGDDadjCC, &
                                GetCrop_ModeCycle(), &
                                GetSimulation_EffectStress_RedCGC(), &
                                GetSimulation_EffectStress_RedCCX()))
                    call SetCrop_CCxAdjusted(GetCCiActual())
                else
                    call SetCCiActual(CanopyCoverNoStressSF(&
                                (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                                GetCrop_DaysToGermination(), &
                                GetCrop_DaysToSenescence(), &
                                GetCrop_DaysToHarvest(), &
                                GetCrop_GDDaysToGermination(), &
                                GetCrop_GDDaysToSenescence(), &
                                GetCrop_GDDaysToHarvest(), &
                                CCoTotal, &
                                (GetCrop_CCxAdjusted() &
                                    /(1._dp &
                                       - GetSimulation_EffectStress_RedCCx() &
                                                                   /100._dp)), &
                                GetCrop_CGC(), CDCTotal, GetCrop_GDDCGC(), &
                                GDDCDCadjusted, SumGDDadjCC, &
                                GetCrop_ModeCycle(), &
                                GetSimulation_EffectStress_RedCGC(), &
                                GetSimulation_EffectStress_RedCCX()))
                end if
                if (GetCCiActual() > CCxSFCD) then
                    call SetCCiActual(CCxSFCD)
                end if
            ! late season
            else
                StressSenescence = undef_int ! to avoid display of zero stress
                                             ! in late season
                if (GetCrop_CCxAdjusted() > CCxSFCD) then
                    call SetCrop_CCxAdjusted(CCxSFCD)
                end if
                if (GetCrop_CCxAdjusted() < 0.01_dp) then
                    call SetCCiActual(0._dp)
                else
                    ! calculate CC in late season
                    ! CCibis = CC which canopy declines
                    ! (soil fertility/salinity stress) further in late season
                    if (GetCrop_GDDaysToSenescence() <= GetCrop_GDDaysToFullCanopySF()) then
                        CCibis = GetCCiActual()
                    else
                        CCibis = CCxSF &
                                - (RatDGDD*GetSimulation_EffectStress_CDecline() &
                                                                       /100._dp) &
                                * (exp(2._dp &
                                      * log(SumGDDadjCC &
                                            - GetCrop_GDDaysToFullCanopySF())) &
                                    /(GetCrop_GDDaysToSenescence() &
                                        - GetCrop_GDDaysToFullCanopySF()))
                    end if
                    if (CCibis < 0._dp) then
                        call SetCCiActual(0._dp)
                    else
                        ! CCiActual = CC with natural senescence in late season
                        GDDCDCadjusted = GetGDDCDCadjustedNoStress(&
                                                    CCxTotal, &
                                                    GDDCDCTotal, &
                                                    GetCrop_CCxAdjusted())
                        if (SumGDDadjCC &
                                < (GetCrop_GDDaysToSenescence() &
                                    + LengthCanopyDecline(&
                                                    GetCrop_CCxAdjusted(), &
                                                    GDDCDCadjusted))) then
                            call SetCCiActual(GetCrop_CCxAdjusted() &
                                 * (1._dp - 0.05_dp &
                                    * (exp((SumGDDadjCC &
                                            - GetCrop_GDDaysToSenescence()) &
                                          * 3.33_dp &
                                          * GDDCDCadjusted &
                                                /(GetCrop_CCxAdjusted() &
                                                            + 2.29_dp)) &
                                        - 1._dp)))
                            ! CCiActual becomes CCibis, when canopy decline is more severe
                            if (CCibis < GetCCiActual()) then
                                call SetCCiActual(CCibis)
                            end if
                        else
                            call SetCCiActual(0._dp)
                        end if
                    end if
                end if
                ! late season
            end if
            ! 3. Canopy can no longer develop (Mid-season (from tFinalCCx)
            ! or Late season stage)
        end if

        ! 4. Canopy senescence due to water stress ?
        if ((SumGDDadjCC < GetCrop_GDDaysToSenescence()) &
                            ! not yet late season stage
            .or. (TimeSenescence > 0._dp)) then
            ! in late season with ongoing early senesence
            ! (TimeSenescence in GDD)
            StressSenescence = 0._dp
            WithBeta = .true.
            Crop_pSenAct_temp = GetCrop_pSenAct()
            call AdjustpSenescenceToETo(GetETo(), TimeSenescence, &
                                        WithBeta, Crop_pSenAct_temp)
            call SetCrop_pSenAct(Crop_pSenAct_temp)
            KsRED = 1._dp ! effect of soil salinity
                          ! on the threshold for senescence
            if (GetSimulation_SWCtopSoilConsidered()) then
                ! top soil is relative wetter than total root zone
                if ((GetRootZoneWC_ZtopAct() &
                        < (GetRootZoneWC_ZtopFC() &
                            - GetCrop_pSenAct() * KsRED &
                                * (GetRootZoneWC_ZtopFC() &
                                    - GetRootZoneWC_ZtopWP()))) &
                    .and. (GetSimulation_ProtectedSeedling() .eqv. .false.)) then
                    TheSenescenceON = .true.
                else
                    TheSenescenceON = .false.
                end if
            else
                if ((GetRootZoneWC_Actual() &
                        < (GetRootZoneWC_FC() &
                            - GetCrop_pSenAct() * KsRED &
                                * (GetRootZoneWC_FC() - GetRootZoneWC_WP()))) &
                    .and. (GetSimulation_ProtectedSeedling() .eqv. .false.)) then
                    TheSenescenceON = .true.
                else
                    TheSenescenceON = .false.
                end if
            end if

            if (TheSenescenceON) then
                ! CanopySenescence
                call SetSimulation_EvapLimitON(.true.)
                ! consider withered crop when not yet in late season
                if (abs(TimeSenescence) < epsilon(0._dp)) then
                    call SetCCiTopEarlySen(GetCCiActual()) ! CC before canopy decline
                end if
                TimeSenescence = TimeSenescence + GDDayi
                call DetermineGDDCDCadjustedWaterStress(GDDCDCadjusted, KsSen)
                if (GetCCiTopEarlySen() < 0.001_dp) then
                    if ((GetSimulation_SumEToStress() &
                            > GetCrop_SumEToDelaySenescence()) &
                      .or. (abs(GetCrop_SumEToDelaySenescence()) < epsilon(0._dp))) then
                        CCiSen = 0._dp ! no crop anymore
                    else
                        if (CCdormant > GetCrop_CCo()) then
                            CCiSen = GetCrop_CCo() &
                                    + (1._dp &
                                       - GetSimulation_SumEToStress() &
                                        / GetCrop_SumEToDelaySenescence()) &
                                    * (CCdormant - GetCrop_CCo())
                        else
                            CCiSen = GetCrop_CCo()
                        end if
                    end if
                else
                    if (((TimeSenescence*GDDCDCadjusted*3.33_dp) &
                                /(GetCCiTopEarlySen()+2.29_dp) > 100._dp) &
                        ! e power too large and in any case CCisen << 0
                        .or. (GetCCiprev() >= 1.05_dp * GetCCiTopEarlySen())) then
                        ! Ln of negative or zero value
                        if ((GetSimulation_SumEToStress() &
                            > GetCrop_SumEToDelaySenescence()) &
                            .or. (abs(GetCrop_SumEToDelaySenescence()) < epsilon(0._dp))) then
                            CCiSen = 0._dp ! no crop anymore
                        else
                            if (CCdormant > GetCrop_CCo()) then
                                CCiSen = GetCrop_CCo() &
                                         + (1._dp &
                                            - GetSimulation_SumEToStress() &
                                                /GetCrop_SumEToDelaySenescence()) &
                                         * (CCdormant - GetCrop_CCo())
                            else
                                CCiSen = GetCrop_CCo()
                            end if
                        end if
                    else
                        ! GDDCDC is adjusted to degree of stress
                        ! time required to reach CCiprev with GDDCDCadjusted
                        GDDtTemp = (log(1._dp &
                                        + (1._dp - GetCCiprev()/GetCCiTopEarlySen()) &
                                                                     /0.05_dp)) &
                                    /(GDDCDCadjusted &
                                        * 3.33_dp/(GetCCiTopEarlySen() + 2.29_dp))
                        ! add 1 day to tTemp and calculate CCiSen with CDCadjusted
                        CCiSen = GetCCiTopEarlySen() &
                                * (1._dp &
                                    - 0.05_dp &
                                        * (exp((GDDtTemp+GDDayi) &
                                                * GDDCDCadjusted &
                                                * 3.33_dp &
                                                /(GetCCiTopEarlySen()+2.29_dp)) &
                                           -1._dp))
                    end if
                    if (CCiSen < 0._dp) then
                        CCiSen = 0._dp
                    end if
                    if ((GetCrop_SumEToDelaySenescence() > 0._dp) &
                        .and. (GetSimulation_SumEToStress() &
                                <= GetCrop_SumEToDelaySenescence())) then
                        if ((CCiSen < GetCrop_CCo()) &
                            .or. (CCiSen < CCdormant)) then
                            if (CCdormant > GetCrop_CCo()) then
                                CCiSen = GetCrop_CCo() &
                                         + (1._dp &
                                            - GetSimulation_SumEToStress() &
                                               /GetCrop_SumEToDelaySenescence()) &
                                            * (CCdormant - GetCrop_CCo())
                            else
                                CCiSen = GetCrop_CCo()
                            end if
                        end if
                    end if
                end if
                if (SumGDDadjCC < GetCrop_GDDaysToSenescence()) then
                    ! before late season
                    if (CCiSen > CCxSFCD) then
                        CCiSen = CCxSFCD
                    end if
                    call SetCCiActual(CCiSen)
                    if (GetCCiActual() > GetCCiPrev()) then
                        call SetCCiActual(GetCCiPrev()) ! to avoid jump in CC
                    end if
                    ! when GDDCGCadjusted increases as a result of watering
                    call SetCrop_CCxAdjusted(GetCCiActual())
                    if (GetCCiActual() < CCoTotal) then
                        call SetCrop_CCoAdjusted(GetCCiActual())
                    else
                        call SetCrop_CCoAdjusted(CCoTotal)
                    end if
                else
                    ! in late season
                    if (CCiSen < GetCCiActual()) then
                        call SetCCiActual(CCiSen)
                    end if
                end if

                if ((roundc(10000._dp*CCiSen, mold=1) <= (10000._dp*CCdormant)) &
                    .or. (roundc(10000._dp*CCiSen, mold=1) &
                            <= roundc(10000._dp*GetCrop_CCo(), mold=1))) then
                    call SetSimulation_SumEToStress(GetSimulation_SumEToStress() &
                                                    + GetETo())
                end if
            else
                ! no water stress, resulting in canopy senescence
                if ((TimeSenescence > 0._dp) &
                    .and. (SumGDDadjCC > GetCrop_GDDaysToSenescence())) then
                    ! rewatering in late season of an early declining canopy
                    Crop_CCxAdjusted_temp = GetCrop_CCxAdjusted()
                    call GetNewCCxandGDDCDC(GetCCiprev(), GDDCDCTotal, &
                                            CCxSF, Crop_CCxAdjusted_temp, &
                                            GDDCDCadjusted)
                    call SetCrop_CCxAdjusted(Crop_CCxAdjusted_temp)
                    call SetCCiActual(CanopyCoverNoStressSF(&
                                (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                                GetCrop_DaysToGermination(), &
                                GetCrop_DaysToSenescence(), &
                                GetCrop_DaysToHarvest(), &
                                GetCrop_GDDaysToGermination(), &
                                GetCrop_GDDaysToSenescence(), &
                                GetCrop_GDDaysToHarvest(), &
                                CCoTotal, &
                                (GetCrop_CCxAdjusted() &
                                    /(1._dp - GetSimulation_EffectStress_RedCCx() &
                                                                       /100._dp)), &
                                GetCrop_CGC(), CDCTotal, GetCrop_GDDCGC(), &
                                GDDCDCadjusted,SumGDDadjCC, &
                                GetCrop_ModeCycle(), &
                                GetSimulation_EffectStress_RedCGC(), &
                                GetSimulation_EffectStress_RedCCX()))
                end if
                TimeSenescence = 0._dp  ! No early senescence or back to normal
                StressSenescence = 0._dp
                call SetSimulation_SumEToStress(0._dp)
            end if
        end if

        ! 5. Adjust Crop.CCxWithered - required for correction
        ! of Transpiration of dying green canopy
        if (GetCCiActual() > GetCrop_CCxWithered()) then
            call SetCrop_CCxWithered(GetCCiActual())
        end if

        ! 6. correction for late-season stage for rounding off errors
        if (SumGDDadjCC > GetCrop_GDDaysToSenescence()) then
            if (GetCCiActual() > GetCCiprev()) then
                call SetCCiActual(GetCCiprev())
            end if
        end if

        ! 7. no crop as a result of fertiltiy and/or water stress
        if (roundc(1000._dp*GetCCiActual(), mold=1) <= 0) then
            NoMoreCrop = .true.
        end if
    end if


    contains


    subroutine DetermineGDDCGCadjusted(GDDCGCadjusted)
        real(dp), intent(inout) :: GDDCGCadjusted

        real(dp) :: Wrelative
        real(dp) :: KsLeaf
        real(dp) :: SWCeffectiveRootZone, FCeffectiveRootZone, &
                    WPeffectiveRootZone

        ! determine FC and PWP
        if (GetSimulation_SWCtopSoilConsidered()) then
            ! top soil is relative wetter than total root zone
            SWCeffectiveRootZone = GetRootZoneWC_ZtopAct()
            Wrelative = (GetRootZoneWC_ZtopFC() &
                            - GetRootZoneWC_ZtopAct()) &
                        /(GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopWP())
                                                                ! top soil
            FCeffectiveRootZone = GetRootZoneWC_ZtopFC()
            WPeffectiveRootZone = GetRootZoneWC_ZtopWP()
        else
            SWCeffectiveRootZone = GetRootZoneWC_Actual()
            Wrelative = (GetRootZoneWC_FC() - GetRootZoneWC_Actual()) &
                            /(GetRootZoneWC_FC() - GetRootZoneWC_WP())
                                                        ! total root zone
            FCeffectiveRootZone = GetRootZoneWC_FC()
            WPeffectiveRootZone = GetRootZoneWC_WP()
        end if

        ! Canopy stress and effect of water stress on CGCGDD
        if (SWCeffectiveRootZone >= FCeffectiveRootZone) then
            GDDCGCadjusted = CGCGDDSF
            StressLeaf = 0._dp
        else
            if (SWCeffectiveRootZone <= WPeffectiveRootZone) then
                GDDCGCadjusted = 0._dp
                StressLeaf = 100._dp
            else
                if (Wrelative <= GetCrop_pLeafAct()) then
                    GDDCGCadjusted = CGCGDDSF
                    StressLeaf = 0._dp
                elseif (Wrelative >= pLeafLLAct) then
                    GDDCGCadjusted = 0._dp
                    StressLeaf = 100._dp
                else
                    KsLeaf = KsAny(Wrelative, GetCrop_pLeafAct(), &
                                   pLeafLLAct, GetCrop_KsShapeFactorLeaf())
                    GDDCGCadjusted = CGCGDDSF * KsLeaf
                    StressLeaf = 100._dp * (1._dp - KsLeaf)
                end if
            end if
        end if
    end subroutine DetermineGDDCGCadjusted


    real(dp) function RequiredGDD(CCiToFind, CCo, CCx, GDDCGCadjusted)
        real(dp), intent(in) :: CCiToFind
        real(dp), intent(in) :: CCo
        real(dp), intent(in) :: CCx
        real(dp), intent(in) :: GDDCGCadjusted

        real(dp) :: GDDCGCx

        ! Only when SumGDDadj > GDDayi
        ! and CCx < CCiToFind
        ! 1. GDDCGCx to reach CCiToFind on previous day (= SumGDDadj - GDDayi )
        if (CCiToFind <= CCx/2._dp) then
            GDDCGCx = (log(CCiToFind/CCo))/(SumGDDadjCC-GDDayi)
        else
            GDDCGCx = (log((0.25_dp*CCx*CCx/CCo) &
                                /(CCx-CCiToFind)))/(SumGDDadjCC-GDDayi)
        end if
        ! 2. Required GDD
        RequiredGDD = (SumGDDadjCC-GDDayi) * GDDCGCx/GDDCGCadjusted
    end function RequiredGDD


    real(dp) function CCatGDDTime(GDDtfictive, CCoGiven, GDDCGCGiven, CCxGiven)
        real(dp), intent(in) :: GDDtfictive
        real(dp), intent(in) :: CCoGiven
        real(dp), intent(in) :: GDDCGCGiven
        real(dp), intent(in) :: CCxGiven

        real(dp) :: CCi

        CCi = CCoGiven * exp(GDDCGCGiven * GDDtfictive)
        if (CCi > CCxGiven/2._dp) then
            CCi = CCxGiven &
                  - 0.25_dp * (CCxGiven/CCoGiven) &
                            * CCxGiven &
                            * exp(-GDDCGCGiven*GDDtfictive)
        end if
        CCatGDDTime = CCi
    end function CCatGDDTime


    subroutine DetermineCCxAdjusted(CCxAdjusted)
        real(dp), intent(inout) :: CCxAdjusted

        real(dp) :: GDDtfictive

        ! 1. find time (GDDtfictive) required to reach CCiPrev
        ! (CCi of previous day) with GDDCGCadjusted
        GDDtfictive = RequiredGDD(GetCCiprev(), GetCrop_CCoAdjusted(), &
                                  CCxSF, GDDCGCadjusted)

        ! 2. Get CCxadjusted (reached at end of stretched crop development)
        if (GDDtfictive > 0._dp) then
            GDDtfictive = GDDtfictive &
                          + (GDDtFinalCCx - SumGDDadjCC) &
                          + GDDayi
            CCxAdjusted = CCatGDDTime(GDDtfictive, GetCrop_CCoAdjusted(), &
                                      GDDCGCadjusted, CCxSF)
        else
            CCxAdjusted = undef_double ! this means CCiActual := CCiPrev
        end if
    end subroutine DetermineCCxAdjusted


    real(dp) function GetGDDCDCadjustedNoStress(CCx, GDDCDC, CCxAdjusted)
        real(dp), intent(in) :: CCx
        real(dp), intent(in) :: GDDCDC
        real(dp), intent(in) :: CCxAdjusted

        real(dp) :: GDDCDCadjusted

        GDDCDCadjusted = GDDCDC * ((CCxadjusted+2.29_dp)/(CCx+2.29_dp))
        GetGDDCDCadjustedNoStress = GDDCDCadjusted
    end function GetGDDCDCadjustedNoStress


    subroutine DetermineGDDCDCadjustedWaterStress(GDDCDCadjusted, KsSen)
        real(dp), intent(inout) :: GDDCDCadjusted
        real(dp), intent(inout) :: KsSen

        real(dp) :: Wrelative
        real(dp) :: pSenLL
        real(dp) :: pSenAct
        logical :: WithBeta

        pSenLL = 0.999_dp ! WP
        if (GetSimulation_SWCtopSoilConsidered()) then
        ! top soil is relative wetter than total root zone
            Wrelative = (GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopAct()) &
                        /(GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopWP())
                                                                ! top soil
        else
            Wrelative = (GetRootZoneWC_FC() - GetRootZoneWC_Actual()) &
                        /(GetRootZoneWC_FC() - GetRootZoneWC_WP())
                                                ! total root zone
        end if

        WithBeta = .false.
        call AdjustpSenescenceToETo(GetETo(), TimeSenescence, &
                                    WithBeta, pSenAct)
        if (Wrelative <= pSenAct) then
            GDDCDCadjusted = 0.0001_dp ! extreme small decline
            StressSenescence = 0._dp
            KsSen = 1._dp
        elseif (Wrelative >= pSenLL) then
            GDDCDCadjusted = GDDCDCTotal &
                            * ((CCxSFCD+2.29_dp)/(CCxTotal+2.29_dp))
                                                        ! full speed
            StressSenescence = 100._dp
            KsSen = 0._dp
        else
            KsSen = KsAny(Wrelative, pSenAct, pSenLL, &
                          GetCrop_KsShapeFactorSenescence())
            if (KsSen > ac_zero_threshold) then
                GDDCDCadjusted = GDDCDCTotal &
                                 * ((CCxSFCD+2.29_dp)/(CCxTotal+2.29_dp)) &
                                 * (1._dp - exp(8._dp*log(KsSen)))
                StressSenescence = 100._dp * (1._dp - KsSen)
            else
                GDDCDCadjusted = 0.0001_dp ! extreme small decline
                StressSenescence = 0._dp
            end if
        end if
    end subroutine DetermineGDDCDCadjustedWaterStress


    subroutine GetNewCCxandGDDCDC(CCiPrev, GDDCDC, CCx, CCxAdjusted, &
                                 GDDCDCadjusted)
        real(dp), intent(in) :: CCiPrev
        real(dp), intent(in) :: GDDCDC
        real(dp), intent(in) :: CCx
        real(dp), intent(inout) :: CCxAdjusted
        real(dp), intent(inout) :: GDDCDCadjusted

        CCxAdjusted = CCiPrev &
                        /(1._dp - 0.05_dp &
                                *(exp((SumGDDadjCC - GDDayi &
                                       - GetCrop_GDDaysToSenescence()) &
                                      * GDDCDC * 3.33_dp/(CCX+2.29_dp))-1._dp))
        GDDCDCadjusted = GDDCDC * (CCxAdjusted+2.29_dp)/(CCx+2.29_dp)
    end subroutine GetNewCCxandGDDCDC
end subroutine DetermineCCiGDD


subroutine EffectSoilFertilitySalinityStress(StressSFadjNEW, Coeffb0Salt, &
                                             Coeffb1Salt, Coeffb2Salt, &
                                             NrDayGrow, StressTotSaltPrev, &
                                             VirtualTimeCC)
    integer(int32), intent(inout) :: StressSFadjNEW
    real(dp), intent(in) :: Coeffb0Salt, Coeffb1Salt, Coeffb2Salt
    integer(int32), intent(in) :: NrDayGrow
    real(dp), intent(in) :: StressTotSaltPrev
    integer(int32), intent(in) :: VirtualTimeCC

    type(rep_EffectStress) :: FertilityEffectStress, SalinityEffectStress
    real(dp) :: SaltStress, CCxRedD
    integer(int8) :: CCxRed
    real(dp) :: ECe_temp, ECsw_temp, ECswFC_temp, KsSalt_temp
    integer(int8) :: RedCGC_temp, RedCCX_temp
    integer(int32) :: Crop_DaysToFullCanopySF_temp
    type(rep_EffectStress) :: EffectStress_temp

    if (GetSimulation_SalinityConsidered()) then
        ECe_temp = GetRootZoneSalt_ECe()
        ECsw_temp = GetRootZoneSalt_ECsw()
        ECswFC_temp = GetRootZoneSalt_ECswFC()
        KsSalt_temp = GetRootZoneSalt_KsSalt()
        call DetermineRootZoneSaltContent(GetRootingDepth(), &
                                          ECe_temp, ECsw_temp, &
                                          ECswFC_temp, KsSalt_temp)
        call SetRootZoneSalt_ECe(ECe_temp)
        call SetRootZoneSalt_ECsw(ECsw_temp)
        call SetRootZoneSalt_ECswFC(ECswFC_temp)
        call SetRootZoneSalt_KsSalt(KsSalt_temp)
        SaltStress = (NrDayGrow*StressTotSaltPrev + 100._dp &
                            *(1._dp-GetRootZoneSalt_KsSalt())) &
                     /(NrDayGrow+1._dp)
    else
        SaltStress = 0._dp
    end if
    if ((VirtualTimeCC < GetCrop_DaysToGermination()) &
            .or. (VirtualTimeCC > (GetCrop_DayN()-GetCrop_Day1())) &
            .or. (GetSimulation_Germinate() .eqv. .false.) &
            .or. ((StressSFAdjNEW == 0) .and. (SaltStress <= 0.1_dp))) then
        ! no soil fertility and salinity stress
        EffectStress_temp = GetSimulation_EffectStress()
        call NoEffectStress(EffectStress_temp)
        call SetSimulation_EffectStress(EffectStress_temp)
        call SetCrop_DaysToFullCanopySF(GetCrop_DaysToFullCanopy())
        if (GetCrop_ModeCycle() == modeCycle_GDDays) then
            call SetCrop_GDDaysToFullCanopySF(GetCrop_GDDaysToFullCanopy())
        end if
    else
        ! Soil fertility
        if (StressSFAdjNEW == 0) then
            call NoEffectStress(FertilityEffectStress)
        else
            call CropStressParametersSoilFertility(GetCrop_StressResponse(), &
                                                   StressSFAdjNEW, &
                                                   FertilityEffectStress)
        end if
        ! Soil Salinity
        CCxRedD = real(roundc(Coeffb0Salt + Coeffb1Salt * SaltStress &
                              + Coeffb2Salt * SaltStress * SaltStress, &
                                                      mold=1), kind=dp)
        if ((CCxRedD < 0._dp) &
                .or. (SaltStress <= 0.1_dp) &
                .or. (GetSimulation_SalinityConsidered() .eqv. .false.)) then
            call NoEffectStress(SalinityEffectStress)
        else
            if ((CCxRedD > 100._dp) .or. (SaltStress >= 99.9_dp)) then
                CCxRed = 100_int8
            else
                CCxRed = roundc(CCxRedD, mold=1_int8)
            end if
            call CropStressParametersSoilSalinity(CCxRed, &
                                                  GetCrop_CCsaltDistortion(), &
                                                  GetCrop_CCo(), &
                                                  GetCrop_CCx(), &
                                                  GetCrop_CGC(), &
                                                  GetCrop_GDDCGC(), &
                                                  GetCrop_DeterminancyLinked(), &
                                                  GetCrop_DaysToFullCanopy(), &
                                                  GetCrop_DaysToFlowering(), &
                                                  GetCrop_LengthFlowering(), &
                                                  GetCrop_DaysToHarvest(), &
                                                  GetCrop_GDDaysToFullCanopy(), &
                                                  GetCrop_GDDaysToFlowering(), &
                                                  GetCrop_GDDLengthFlowering(), &
                                                  GetCrop_GDDaysToHarvest(), &
                                                  GetCrop_ModeCycle(), &
                                                  SalinityEffectStress)
        end if
        ! Assign integrated effect of the stresses
        call SetSimulation_EffectSTress_RedWP(FertilityEffectStress%RedWP)
        call SetSimulation_EffectSTress_RedKsSto(SalinityEffectStress%RedKsSto)
        if (FertilityEffectStress%RedCGC > SalinityEffectStress%RedCGC) then
            call SetSimulation_EffectSTress_RedCGC(FertilityEffectStress%RedCGC)
        else
            call SetSimulation_EffectSTress_RedCGC(SalinityEffectStress%RedCGC)
        end if
        if (FertilityEffectStress%RedCCX > SalinityEffectStress%RedCCX) then
            call SetSimulation_EffectSTress_RedCCX(FertilityEffectStress%RedCCX)
        else
            call SetSimulation_EffectSTress_RedCCX(SalinityEffectStress%RedCCX)
        end if
        if (FertilityEffectStress%CDecline > SalinityEffectStress%CDecline) then
            call SetSimulation_EffectSTress_CDecline(FertilityEffectStress%CDecline)
        else
            call SetSimulation_EffectSTress_CDecline(SalinityEffectStress%CDecline)
        end if
        ! adjust time to maximum canopy cover
        RedCGC_temp = GetSimulation_EffectStress_RedCGC()
        RedCCX_temp = GetSimulation_EffectStress_RedCCX()
        Crop_DaysToFullCanopySF_temp = GetCrop_DaysToFullCanopySF()
        call TimeToMaxCanopySF(GetCrop_CCo(), GetCrop_CGC(), GetCrop_CCx(), &
                               GetCrop_DaysToGermination(), &
                               GetCrop_DaysToFullCanopy(), &
                               GetCrop_DaysToSenescence(), &
                               GetCrop_DaysToFlowering(), &
                               GetCrop_LengthFlowering(), &
                               GetCrop_DeterminancyLinked(), &
                               Crop_DaysToFullCanopySF_temp, RedCGC_temp, &
                               RedCCX_temp, StressSFAdjNEW)
        call SetSimulation_EffectStress_RedCGC(RedCGC_temp)
        call SetSimulation_EffectStress_RedCCX(RedCCX_temp)
        call SetCrop_DaysToFullCanopySF(Crop_DaysToFullCanopySF_temp)
        if (GetCrop_ModeCycle() == modeCycle_GDDays) then
            if ((abs(GetManagement_FertilityStress()) > epsilon(0._dp)) &
                    .or. (abs(SaltStress) > epsilon(0._dp))) then
                call SetCrop_GDDaysToFullCanopySF(&
                             GrowingDegreeDays(GetCrop_DaysToFullCanopySF(), &
                                               GetCrop_Day1(), &
                                               GetCrop_Tbase(), &
                                               GetCrop_Tupper(), &
                                               GetSimulParam_Tmin(), &
                                               GetSimulParam_Tmax()))
            else
                call SetCrop_GDDaysToFullCanopySF(GetCrop_GDDaysToFullCanopy())
            end if
        end if
    end if


    contains


    subroutine NoEffectStress(TheEffectStress)
        type(rep_EffectStress), intent(inout) :: TheEffectStress

        TheEffectStress%RedCGC = 0._dp
        TheEffectStress%RedCCX = 0._dp
        TheEffectStress%RedWP = 0._dp
        TheEffectStress%CDecline = 0._dp
        TheEffectStress%RedKsSto = 0._dp
    end subroutine NoEffectStress
end subroutine EffectSoilFertilitySalinityStress


subroutine PrepareStage1()
    type(rep_Soil) :: Soil_temp

    Soil_temp = GetSoil()

    if (GetSurfaceStorage() > ac_zero_threshold) then
        call SetSimulation_EvapWCsurf(Soil_temp%REW*1._dp)
    else
        call SetSimulation_EvapWCsurf(GetRain() + GetIrrigation() - GetRunOff())
        if (GetSimulation_EvapWCsurf() > Soil_temp%REW) then
            call SetSimulation_EvapWCsurf(Soil_temp%REW*1._dp)
        end if
    end if
    call SetSimulation_EvapStartStg2(int(undef_int,kind=int8))
    call SetSimulation_EvapZ(EvapZmin/100._dp)
end subroutine PrepareStage1


real(dp) function WCEvapLayer(Zlayer, AtTheta)
    real(dp), intent(in) :: Zlayer
    integer(intEnum), intent(in) :: AtTheta

    real(dp) :: Ztot, Wx, fracZ
    integer(int32) :: compi

    Wx = 0.0_dp
    Ztot = 0.0_dp
    compi = 0
    do while ((abs(Zlayer-Ztot) > 0.0001_dp) &
            .and. (compi < GetNrCompartments()))
        compi = compi + 1
        if ((Ztot + GetCompartment_Thickness(compi)) > Zlayer) then
            fracZ = (Zlayer - Ztot)/(GetCompartment_Thickness(compi))
        else
            fracZ = 1._dp
        end if
        select case (AtTheta)
            case(whichtheta_AtSAT)
            Wx = Wx + 10._dp &
                    * GetSoilLayer_SAT(GetCompartment_Layer(compi)) &
                    * fracZ * GetCompartment_Thickness(compi) &
                    * (1._dp &
                        - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                    /100._dp)
            case (whichtheta_AtFC)
            Wx = Wx + 10._dp &
                    * GetSoilLayer_FC(GetCompartment_Layer(compi)) &
                    * fracZ * GetCompartment_Thickness(compi) &
                    * (1._dp &
                        - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                    /100._dp)
            case (whichtheta_AtWP)
            Wx = Wx + 10._dp &
                    * GetSoilLayer_WP(GetCompartment_Layer(compi)) &
                    * fracZ * GetCompartment_Thickness(compi) &
                    * (1._dp &
                        - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                    /100._dp)
            case default
                Wx = Wx + 1000._dp &
                        * GetCompartment_Theta(compi) * fracZ &
                        * GetCompartment_Thickness(compi) &
                        * (1._dp &
                            - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                        /100._dp)
        end select
        Ztot = Ztot + fracZ * GetCompartment_Thickness(compi)
    end do
    WCEvapLayer = Wx
end function WCEvapLayer


subroutine PrepareStage2()

    integer(intEnum) :: AtTheta
    integer(int8) :: EvapStartStg2
    real(dp) :: WSAT, WFC, Wact

    call SetSimulation_EvapZ(EvapZmin/100)
    AtTheta = whichtheta_AtSat
    WSAT = WCEvapLayer(GetSimulation_EvapZ(), AtTheta)
    AtTheta = whichtheta_AtFC
    WFC = WCEvapLayer(GetSimulation_EvapZ(), AtTheta)
    AtTheta = whichtheta_AtAct
    Wact = WCEvapLayer(GetSimulation_EvapZ(), AtTheta)

    if ((Wact - (WFC-GetSoil_REW())) <= epsilon(0._dp)) then
        EvapStartStg2 = 0_int8
    else
        EvapStartStg2 = roundc(100._dp * (Wact - (WFC-GetSoil_REW())) &
                               / (WSAT - (WFC-GetSoil_REW())), mold=1_int8)
    end if
    call SetSimulation_EvapStartStg2(EvapStartStg2)
end subroutine PrepareStage2


subroutine CalculateEvaporationSurfaceWater()

    real(dp) :: SaltSurface

    if (GetSurfaceStorage() > GetEpot()) then
        SaltSurface = GetSurfaceStorage()*GetECstorage()*Equiv
        call SetEact(GetEpot())
        call SetSurfaceStorage(GetSurfaceStorage() - GetEact())
        call SetECstorage(SaltSurface/(GetSurfaceStorage()*Equiv))
            ! salinisation of surface storage layer
    else
        call SetEact(GetSurfaceStorage())
        call SetSurfaceStorage(0._dp)
        call SetSimulation_EvapWCsurf(real(GetSoil_REW(), kind=dp))
        call SetSimulation_EvapZ(EvapZmin/100._dp)
        if (GetSimulation_EvapWCsurf() < 0.0001_dp) then
            call PrepareStage2()
        else
            call SetSimulation_EvapStartStg2(int(undef_int, kind=int8))
        end if
    end if
end subroutine CalculateEvaporationSurfaceWater


subroutine AdjustEpotMulchWettedSurface(dayi, EpotTot, Epot, EvapWCsurface)
    integer(int32), intent(in) :: dayi
    real(dp), intent(in) :: EpotTot
    real(dp), intent(inout) :: Epot
    real(dp), intent(inout) :: EvapWCsurface

    real(dp) :: EpotIrri

    ! 1. Mulches (reduction of EpotTot to Epot)
    if (GetSurfaceStorage() <= ac_zero_threshold) then
        if (dayi < GetCrop_Day1()) then ! before season
            Epot = EpotTot &
                    * (1._dp - (GetManagement_EffectMulchOffS()/100._dp) &
                               *(GetManagement_SoilCoverBefore()/100._dp))
        else
            if (dayi < GetCrop_Day1()+GetCrop_DaysToHarvest()) then ! in season
                Epot = EpotTot &
                        * (1._dp &
                            - (GetManagement_EffectMulchInS()/100._dp) &
                             * (GetManagement_Mulch()/100._dp))
            else
                Epot = EpotTot &
                        * (1._dp &
                            - (GetManagement_EffectMulchOffS()/100._dp) &
                             * (GetManagement_SoilCoverAfter()/100._dp))
            end if
        end if
    else
        Epot = EpotTot ! flooded soil surface
    end if

    ! 2a. Entire soil surface wetted ?
    if (GetIrrigation() > 0._dp) then
        ! before season
        if ((dayi < GetCrop_Day1()) &
            .and. (GetSimulParam_IrriFwOffSeason() < 100)) then
            call SetEvapoEntireSoilSurface(.false.)
        end if
        ! in season
        if ((dayi >= GetCrop_Day1()) &
            .and. (dayi < GetCrop_Day1()+GetCrop_DaysToHarvest()) &
            .and. (GetSimulParam_IrriFwInSeason() < 100)) then
            call SetEvapoEntireSoilSurface(.false.)
        end if
        ! after season
        if ((dayi >= GetCrop_Day1()+GetCrop_DaysToHarvest()) &
            .and.(GetSimulParam_IrriFwOffSeason() < 100)) then
            call SetEvapoEntireSoilSurface(.false.)
        end if
    end if
    if ((GetRain() > 1._dp) .or. (GetSurfaceStorage() > 0._dp)) then
        call SetEvapoEntireSoilSurface(.true.)
    end if
    if ((dayi >= GetCrop_Day1()) &
        .and. (dayi < GetCrop_Day1()+GetCrop_DaysToHarvest()) &
        .and. (GetIrriMode() == IrriMode_Inet)) then
        call SetEvapoEntireSoilSurface(.true.)
    end if

    ! 2b. Correction for Wetted surface by Irrigation
    if (.not.GetEvapoEntireSoilSurface()) then
        if ((dayi >= GetCrop_Day1()) &
            .and. (dayi < GetCrop_Day1()+GetCrop_DaysToHarvest())) then
            ! in season
            EvapWCsurface = EvapWCsurface &
                            * (GetSimulParam_IrriFwInSeason()/100._dp)
            EpotIrri = EpotTot * (GetSimulParam_IrriFwInSeason()/100._dp)
        else
            ! off-season
            EvapWCsurface = EvapWCsurface &
                            * (GetSimulParam_IrriFwOffSeason()/100._dp)
            EpotIrri = EpotTot * (GetSimulParam_IrriFwOffSeason()/100._dp)
        end if
        if (GetEact() > EpotIrri) then
            EpotIrri = GetEact()  ! Eact refers to the previous day
        end if
        if (EpotIrri < Epot) then
            Epot = Epotirri
        end if
    end if
end subroutine AdjustEpotMulchWettedSurface


subroutine ConcentrateSalts()

    integer(int32) :: compi, celWet, celi
    real(dp) :: SaltTot, mm
    real(dp) :: Salt_temp, Depo_temp

    do compi = 1, GetNrCompartments()
        SaltTot = 0.0_dp
        celWet = ActiveCells(GetCompartment_i(compi))
        if (celWet < GetSoilLayer_SCP1(GetCompartment_Layer(compi))) then
            do celi = (celWet+1), GetSoilLayer_SCP1(GetCompartment_Layer(compi))
                SaltTot = SaltTot + GetCompartment_Salt(compi, celi)&
                          + GetCompartment_Depo(compi, celi)
                call SetCompartment_Salt(compi, celi, 0.0_dp)
                call SetCompartment_Depo(compi, celi, 0.0_dp)
            end do
        end if
        if (SaltTot > 0.0_dp) then
            call SetCompartment_Salt(compi, celWet, &
                GetCompartment_Salt(compi, celWet) + SaltTot)
            mm = GetSoilLayer_Dx(GetCompartment_Layer(compi))*1000.0_dp&
                 * GetCompartment_Thickness(compi)&
                 * (1 - GetSoilLayer_GravelVol(GetCompartment_Layer(compi))&
                        / 100.0_dp)
            Salt_temp = GetCompartment_Salt(compi, celWet)
            Depo_temp = GetCompartment_Depo(compi, celWet)
            call SaltSolutionDeposit(mm, Salt_temp, Depo_temp)
            call SetCompartment_Salt(compi, celWet, Salt_temp)
            call SetCompartment_Depo(compi, celWet, Depo_temp)
        end if
    end do
end subroutine ConcentrateSalts


subroutine ExtractWaterFromEvapLayer(EvapToLose, Zact, Stg1)
    real(dp), intent(in) :: EvapToLose
    real(dp), intent(in) :: Zact
    logical, intent(in) :: Stg1

    real(dp) :: EvapLost, Wx, Wairdry, AvailableW, &
                Ztot, fracZ, StillToExtract
    integer(int32) :: compi

    EvapLost = 0._dp
    compi = 0
    Ztot = 0._dp
    loop: do
        compi = compi + 1
        if ((Ztot + GetCompartment_Thickness(compi)) > Zact) then
            fracZ = (Zact-Ztot)/GetCompartment_Thickness(compi)
        else
            fracZ = 1._dp
        end if
        Wairdry = 10._dp &
                  * GetSoilLayer_WP(GetCompartment_Layer(compi))/2._dp &
                  * GetCompartment_Thickness(compi) &
                  * (1._dp &
                     - GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._dp)
        Wx = 1000._dp * GetCompartment_Theta(compi) &
             * GetCompartment_Thickness(compi) &
             * (1._dp &
                  - GetSoilLayer_GravelVol(GetCompartment_Layer(compi))/100._dp)
        AvailableW = (Wx-Wairdry)*fracZ
        StillToExtract = (EvapToLose-EvapLost)
        if (AvailableW > 0._dp) then
            if (AvailableW > StillToExtract) then
                call SetEact(GetEact() + StillToExtract)
                EvapLost = EvapLost + StillToExtract
                Wx = Wx - StillToExtract
            else
                call SetEact(GetEact() + AvailableW)
                EvapLost = EvapLost + AvailableW
                Wx = Wx - AvailableW
            end if
            call SetCompartment_Theta(&
                    compi, &
                    Wx/(1000._dp * GetCompartment_Thickness(compi) &
                        * (1._dp &
                           - GetSoilLayer_GravelVol(GetCompartment_Layer(compi))&
                                                                     /100._dp)))
        end if
        Ztot = Ztot + fracZ * (GetCompartment_Thickness(compi))
        if ((Compi >= GetNrCompartments()) &
            .or. (abs(StillToExtract) < ac_zero_threshold) &
            .or. (Ztot >= 0.999999_dp*Zact)) exit loop
    end do loop
    if (Stg1) then
        call SetSimulation_EvapWCsurf(GetSimulation_EvapWCsurf() - EvapLost)
        if (abs(EvapToLose-EvapLost) > 0.0001_dp) then
            ! not enough water left in the compartment to store WCsurf
            call SetSimulation_EvapWCsurf(0._dp)
        end if
    end if
end subroutine ExtractWaterFromEvapLayer


subroutine CalculateSoilEvaporationStage1()

    real(dp) :: Eremaining
    logical :: Stg1

    Stg1 = .true.
    Eremaining = GetEpot() - GetEact()
    if (GetSimulation_EvapWCsurf() > Eremaining) then
        call ExtractWaterFromEvapLayer(Eremaining, EvapZmin, Stg1)
    else
        call ExtractWaterFromEvapLayer(GetSimulation_EvapWCsurf(), EvapZmin, &
                                                                       Stg1)
    end if
    if (GetSimulation_EvapWCsurf() < ac_zero_threshold) then
        call PrepareStage2()
    end if
end subroutine CalculateSoilEvaporationStage1


subroutine CalculateSoilEvaporationStage2()

    integer(int32), parameter :: NrOfStepsInDay = 20
    real(dp), parameter :: FractionWtoExpandZ = 0.4_dp
    integer(intEnum) :: AtTheta
    real(dp) :: Wupper, Wlower, Wact, Eremaining, Wrel, Kr, &
                Elost, MaxSaltExDepth, SX, Zi, SaltDisplaced, UL, DeltaX
    integer(int32) :: i, compi, SCell1, SCellEnd
    logical :: Stg1, BoolCell
    real(dp), dimension(11) :: ThetaIniEvap
    integer(int32), dimension(11) :: SCellIniEvap

    ! Step 1. Conditions before soil evaporation
    compi = 1
    MaxSaltExDepth = GetCompartment_Thickness(1)
    do while ((MaxSaltExDepth < GetSimulParam_EvapZmax()) &
                .and. (compi < GetNrCompartments()))
        compi = compi + 1
        ThetaIniEvap(compi-1) = GetCompartment_Theta(compi)
        SCellIniEvap(compi-1) = ActiveCells(GetCompartment_i(compi))
        MaxSaltExDepth = MaxSaltExDepth + GetCompartment_Thickness(compi)
    end do

    ! Step 2. Soil evaporation
    Stg1 = .false.
    Eremaining = GetEpot() - GetEact()
    call GetLimitsEvapLayer(real(GetSimulation_EvapStartStg2(), kind=dp), &
                            Wupper, Wlower)
    do i = 1, NrOfStepsInDay
        AtTheta = whichtheta_AtAct
        Wact = WCEvapLayer(GetSimulation_EvapZ(), AtTheta)
        Wrel = (Wact-Wlower)/(Wupper-Wlower)
        if (GetSimulParam_EvapZmax() > EvapZmin) then
            do while ((Wrel < (FractionWtoExpandZ &
                * (GetSimulParam_EvapZmax() &
                    -(100._dp*GetSimulation_EvapZ())) &
                        /(GetSimulParam_EvapZmax()-EvapZmin))) &
                .and. (GetSimulation_EvapZ() &
                            < GetSimulParam_EvapZmax()/100._dp))
                call SetSimulation_EvapZ(GetSimulation_EvapZ() + 0.001_dp)
                                                                ! add 1 mm
                call GetLimitsEvapLayer(real(GetSimulation_EvapStartStg2(), kind=dp), &
                                        Wupper, Wlower)
                AtTheta = whichtheta_AtAct
                Wact = WCEvapLayer(GetSimulation_EvapZ(), AtTheta)
                Wrel = (Wact-Wlower)/(Wupper-Wlower)
            end do
            Kr = SoilEvaporationReductionCoefficient(Wrel, &
                               real(GetSimulParam_EvapDeclineFactor(), kind=dp))
        end if
        if (abs(GetETo() - 5._dp) > 0.01_dp) then
            ! correction for evaporative demand
            ! adjustment of Kr (not considered yet)
        end if
        Elost = Kr * (Eremaining/NrOfStepsInDay)
        call ExtractWaterFromEvapLayer(Elost, GetSimulation_EvapZ(), Stg1)
    end do

    ! Step 3. Upward salt transport
    SX = SaltTransportFactor(Getcompartment_Theta(1))
    if (SX > 0.01_dp) then
        SCell1 = ActiveCells(GetCompartment_i(1))
        compi = 2
        Zi = GetCompartment_Thickness(1) + GetCompartment_Thickness(2)
        do while ((roundc(Zi*100._dp, mold=1) &
                    <= roundc(MaxSaltExDepth*100._dp, mold=1)) &
            .and. (compi <= GetNrCompartments()) &
            .and. (roundc(ThetaIniEvap(compi-1)*100000._dp, mold=1) &
                    /= roundc(GetCompartment_theta(compi)*100000._dp, mold=1)))
            ! move salt to compartment 1
            SCellEnd = ActiveCells(GetCompartment_i(compi))
            BoolCell = .false.
            UL = GetSoilLayer_UL(GetCompartment_Layer(compi))
            DeltaX = GetSoilLayer_Dx(GetCompartment_Layer(compi))
            loop: do
                if (SCellEnd < SCellIniEvap(compi-1)) then
                    SaltDisplaced = SX &
                         * GetCompartment_Salt(compi, SCellIniEvap(compi-1))
                    call SetCompartment_Salt(compi, SCellIniEvap(compi-1), &
                                             GetCompartment_Salt(compi, &
                                             SCellIniEvap(compi-1)) - SaltDisplaced)
                    SCellIniEvap(compi-1) = SCellIniEvap(compi-1) - 1
                    ThetaIniEvap(compi-1) = DeltaX * SCellIniEvap(compi-1)
                else
                    BoolCell = .true.
                    if (SCellEnd == GetSoilLayer_SCP1(GetCompartment_Layer(compi))) then
                        SaltDisplaced = SX &
                            * GetCompartment_Salt(compi, SCellIniEvap(compi)) &
                            * (ThetaIniEvap(compi-1) &
                                - GetCompartment_theta(compi)) &
                                            /(ThetaIniEvap(compi-1)-UL)
                    else
                        SaltDisplaced = SX &
                            * GetCompartment_Salt(compi, SCellIniEvap(compi-1)) &
                            * (ThetaIniEvap(compi-1) &
                                - GetCompartment_theta(compi)) &
                                  /(ThetaIniEvap(compi-1)-(DeltaX*(SCellEnd-1)))
                    end if
                    call SetCompartment_Salt(compi, SCellIniEvap(compi-1), &
                              GetCompartment_Salt(compi, SCellIniEvap(compi-1)) &
                                            - SaltDisplaced)
                end if
                call SetCompartment_Salt(1, SCell1, &
                               GetCompartment_Salt(1, SCell1) + SaltDisplaced)
                if (BoolCell) exit loop
            end do loop
            compi = compi + 1
            if (compi <= GetNrCompartments()) then
                Zi = Zi + GetCompartment_Thickness(compi)
            end if
        end do
    end if


    contains


    subroutine GetLimitsEvapLayer(xProc, Wupper, Wlower)
        real(dp), intent(in) :: xProc
        real(dp), intent(inout) :: Wupper
        real(dp), intent(inout) :: Wlower

        integer(intEnum) :: AtTheta
        real(dp) :: WSAT, WFC

        AtTheta = whichtheta_AtSat
        WSAT = WCEvapLayer(GetSimulation_EvapZ(), AtTheta)
        AtTheta = whichtheta_AtFC
        WFC = WCEvapLayer(GetSimulation_EvapZ(), AtTheta)
        Wupper = (xProc/100._dp) &
                 * (WSAT - (WFC-GetSoil_REW())) &
                 + (WFC-GetSoil_REW())
        AtTheta = whichtheta_AtWP
        Wlower = WCEvapLayer(GetSimulation_EvapZ(), AtTheta)/2
    end subroutine GetLimitsEvapLayer


    real(dp) function SaltTransportFactor(theta)
        real(dp), intent(in) :: theta

        real(dp) :: x

        if (theta <= GetSoilLayer_WP(1)/200._dp) then
            SaltTransportFactor = 0._dp
        else
            x = (theta*100._dp - GetSoilLayer_WP(1)/2._dp) &
                /(GetSoilLayer_SAT(1) - GetSoilLayer_WP(1)/2._dp)
            SaltTransportFactor = exp(x*log(10._dp)+log(x/10._dp))
        end if
    end function SaltTransportFactor
end subroutine CalculateSoilEvaporationStage2


subroutine DetermineCCi(CCxTotal, CCoTotal, StressLeaf, FracAssim, &
                        MobilizationON, StorageON, Tadj, VirtualTimeCC, &
                        StressSenescence, TimeSenescence, NoMoreCrop, &
                        CDCTotal, DayFraction, &
                        GDDCDCTotal, TESTVAL)
    real(dp), intent(in) :: CCxTotal
    real(dp), intent(in) :: CCoTotal
    real(dp), intent(inout) :: StressLeaf
    real(dp), intent(in) :: FracAssim
    logical, intent(in) :: MobilizationON
    logical, intent(in) :: StorageON
    integer(int32), intent(in) :: Tadj
    integer(int32), intent(in) :: VirtualTimeCC
    real(dp), intent(inout) :: StressSenescence
    real(dp), intent(inout) :: TimeSenescence
    logical, intent(inout) :: NoMoreCrop
    real(dp), intent(in) :: CDCTotal
    real(dp), intent(in) :: DayFraction
    real(dp), intent(in) :: GDDCDCTotal
    real(dp), intent(inout) :: TESTVAL

    real(dp), parameter :: CCdormant = 0.05_dp
    real(dp) :: pLeafLLAct , CGCadjusted, CDCadjusted, &
                CCiSen, tTemp, CCxSF, CGCSF, CCxSFCD, KsRED, CCibis
    integer(int32) :: tFinalCCx
    logical :: WithBeta
    logical :: TheSenescenceON
    real(dp) :: KsSen
    !! test Version 6.2
    real(dp) :: Crop_pLeafAct_temp
    real(dp) :: Crop_pSenAct_temp
    real(dp) :: Crop_CCxAdjusted_temp

    ! DetermineCCi
    if ((VirtualTimeCC < GetCrop_DaysToGermination()) &
        .or. (VirtualTimeCC > (GetCrop_DayN()-GetCrop_Day1()))) then
        call SetCCiActual(0._dp)
    else
        ! growing season (once germinated)
        ! 1. find some parameters
        CGCSF = GetCrop_CGC() &
                * (1._dp - GetSimulation_EffectStress_RedCGC()/100._dp)
        CGCadjusted = CGCSF
        CCxSF = CCxTotal &
                * (1._dp - GetSimulation_EffectStress_RedCCX()/100._dp)

        ! maximum canopy cover than can be reached
        ! (considering soil fertility/salinity, weed stress)
        if (VirtualTimeCC <= GetCrop_DaysToFullCanopySF()) then
            CCxSFCD = CCxSF ! no correction before maximum canopy is reached
        else
            if (VirtualTimeCC < GetCrop_DaysToSenescence()) then
                CCxSFCD = CCiNoWaterStressSF(&
                            (VirtualTimeCC + GetSimulation_DelayedDays()+1), &
                            GetCrop_DaysToGermination(), &
                            GetCrop_DaysToFullCanopySF(), &
                            GetCrop_DaysToSenescence(), &
                            GetCrop_DaysToHarvest(), &
                            GetCrop_GDDaysToGermination(), &
                            GetCrop_GDDaysToFullCanopySF(), &
                            GetCrop_GDDaysToSenescence(), &
                            GetCrop_GDDaysToHarvest(), &
                            CCoTotal, CCxTotal, GetCrop_CGC(), &
                            GetCrop_GDDCGC(), CDCTotal, GDDCDCTotal, &
                            GetSimulation_SumGDD(), 1._dp, &
                            GetSimulation_EffectStress_RedCGC(), &
                            GetSimulation_EffectStress_RedCCX(), &
                            GetSimulation_EffectStress_CDecline(), &
                            GetCrop_ModeCycle())
            else
                CCxSFCD = CCxSF &
                          - (GetSimulation_EffectStress_CDecline()/100._dp) &
                          * (GetCrop_DaysToSenescence() &
                                - GetCrop_DaysToFullCanopySF())
            end if
            if (CCxSFCD < 0._dp) then
                CCxSFCD = 0._dp
            end if
        end if
        StressLeaf = undef_int
        if (VirtualTimeCC == GetCrop_DaysToGermination()) then
            call SetCCiPrev(CCoTotal)
        end if

        ! time of potentional vegetative growth
        tFinalCCx = GetCrop_DaysToSenescence() ! undeterminant crop
        if ((GetCrop_subkind() == subkind_Grain) &
                .and. (GetCrop_DeterminancyLinked())) then
            ! determinant crop
            ! reduce tFinalCC in f(determinancy of crop)
            if (GetCrop_DaysToCCini() /= 0) then
                ! regrowth  (adjust to slower time)
                tFinalCCx = GetCrop_DaysToFullCanopy() &
                            + roundc(DayFraction &
                                * ((GetCrop_DaysToFlowering() &
                                    + (GetCrop_LengthFlowering()/2._dp) &
                                    - GetSimulation_DelayedDays()) &
                                    + Tadj + GetCrop_DaysToGermination() &
                                    - GetCrop_DaysToFullCanopy()), mold=1)
            else
                ! sown or transplant
                tFinalCCx = GetCrop_DaysToFlowering() &
                            + roundc(GetCrop_LengthFlowering()/2._dp, mold=1)
            end if
            if (tFinalCCx > GetCrop_DaysToSenescence()) then
                tFinalCCx = GetCrop_DaysToSenescence()
            end if
        end if

        ! Crop.pLeafAct and Crop.pSenAct for
        ! plotting root zone depletion in RUN
        Crop_pLeafAct_temp = GetCrop_pLeafAct()
        call AdjustpLeafToETo(GetETo(), Crop_pLeafAct_temp, pLeafLLAct)
        call SetCrop_pLeafAct(Crop_pLeafAct_temp)
        WithBeta = .true.
        Crop_pSenAct_temp = GetCrop_pSenAct()
        call AdjustpSenescenceToETo(GetETo(), TimeSenescence, WithBeta, &
                                    Crop_pSenAct_temp)
        call SetCrop_pSenAct(Crop_pSenAct_temp)

        ! 2. Canopy can still develop (stretched to tFinalCCx)
        if (VirtualTimeCC < tFinalCCx) then
            ! Canopy can stil develop (stretched to tFinalCCx)
            if ((GetCCiPrev() <= GetCrop_CCoAdjusted()) &
                .or. (VirtualTimeCC <= 1) &
                .or. ((GetSimulation_ProtectedSeedling()) &
                    .and. (GetCCiPrev() <= (1.25_dp * CCoTotal)))) then
                ! 2.a first day or very small CC as a result of senescence
                ! (no adjustment for leaf stress)
                if (GetSimulation_ProtectedSeedling()) then
                    call SetCCiActual(CanopyCoverNoStressSF(&
                            (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                            GetCrop_DaysToGermination(), &
                            GetCrop_DaysToSenescence(), &
                            GetCrop_DaysToHarvest(), &
                            GetCrop_GDDaysToGermination(), &
                            GetCrop_GDDaysToSenescence(), &
                            GetCrop_GDDaysToHarvest(), &
                            CCoTotal, CCxTotal, GetCrop_CGC(), &
                            CDCTotal, GetCrop_GDDCGC(), GDDCDCTotal, &
                            GetSimulation_SumGDD(), GetCrop_ModeCycle(), &
                            GetSimulation_EffectStress_RedCGC(), &
                            GetSimulation_EffectStress_RedCCX()))
                    if (GetCCiActual() > (1.25_dp * CCoTotal)) then
                        call SetSimulation_ProtectedSeedling(.false.)
                    end if
                else
                    ! this results in CC increase when during senescence CC
                    ! becomes smaller than CCini)
                    if (VirtualTimeCC == 1) then
                        call SetCCiActual(GetCrop_CCoAdjusted() &
                                            * exp(CGCSF*2._dp))
                    else
                        call SetCCiActual(GetCrop_CCoAdjusted() &
                                            * exp(CGCSF*1._dp))
                    end if
                end if

                ! 2.b CC > CCo
            else
                if (GetCCiPrev() < 0.97999_dp*CCxSF) then
                    call DetermineCGCadjusted(CGCadjusted)
                    if (CGCadjusted > ac_zero_threshold) then
                        ! CGCSF or CGCadjusted > 0
                        Crop_CCxAdjusted_temp = GetCrop_CCxAdjusted()
                        call DetermineCCxAdjusted(Crop_CCxAdjusted_temp)
                        call SetCrop_CCxAdjusted(Crop_CCxAdjusted_temp)
                        if (GetCrop_CCxAdjusted() < 0) then
                            call SetCCiActual(GetCCiPrev())
                        elseif (abs(GetCCiPrev() - 0.97999_dp*CCxSF) < 0.001_dp) then
                            call SetCCiActual(CanopyCoverNoStressSF(&
                                (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                                GetCrop_DaysToGermination(), &
                                GetCrop_DaysToSenescence(), &
                                GetCrop_DaysToHarvest(), &
                                GetCrop_GDDaysToGermination(), &
                                GetCrop_GDDaysToSenescence(), &
                                GetCrop_GDDaysToHarvest(), &
                                CCoTotal, CCxTotal, GetCrop_CGC(), &
                                CDCTotal, GetCrop_GDDCGC(), GDDCDCTotal, &
                                GetSimulation_SumGDD(), GetCrop_ModeCycle(), &
                                GetSimulation_EffectStress_RedCGC(), &
                                GetSimulation_EffectStress_RedCCX()))
                        else
                            tTemp = RequiredTimeNew(GetCCiPrev(), &
                                                    GetCrop_CCoAdjusted(), &
                                                    GetCrop_CCxAdjusted(), &
                                                    CGCadjusted)
                            if (tTemp < 0._dp) then
                                call SetCCiActual(GetCCiPrev())
                            else
                                tTemp = tTemp + 1._dp
                                call SetCCiActual(CCatTime(&
                                        tTemp, GetCrop_CCoAdjusted(), CGCadjusted, &
                                        GetCrop_CCxAdjusted()))
                            end if
                        end if
                    else
                        ! CGCadjusted = 0 - too dry for leaf expansion
                        call SetCCiActual(GetCCiPrev())
                        if (GetCCiActual() > GetCrop_CCoAdjusted()) then
                            call SetCrop_CCoAdjusted(CCoTotal)
                        else
                            call SetCrop_CCoAdjusted(GetCCiActual())
                        end if
                    end if
                else
                    call SetCCiActual(CanopyCoverNoStressSF(&
                            (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                            GetCrop_DaysToGermination(), &
                            GetCrop_DaysToSenescence(), &
                            GetCrop_DaysToHarvest(), &
                            GetCrop_GDDaysToGermination(), &
                            GetCrop_GDDaysToSenescence(), &
                            GetCrop_GDDaysToHarvest(), &
                            CCoTotal, CCxTotal, GetCrop_CGC(), CDCTotal, &
                            GetCrop_GDDCGC(), GDDCDCTotal, &
                            GetSimulation_SumGDD(), GetCrop_ModeCycle(), &
                            GetSimulation_EffectStress_RedCGC(), &
                            GetSimulation_EffectStress_RedCCX()))
                    call SetCrop_CCoAdjusted(CCoTotal)
                    StressLeaf = -33._dp ! maximum canopy is reached;
                    ! no increase anymore of CGC after cutting
                end if
                if (GetCCiActual() > CCxSFCD) then
                    call SetCCiActual(CCxSFCD)
                    StressLeaf = -33._dp ! maximum canopy is reached;
                    ! no increase anymore of CGC after cutting
                end if
            end if
            call SetCrop_CCxAdjusted(GetCCiActual())

            ! 3. Canopy can no longer develop (Mid-season (from tFinalCCx) or Late season stage)
        else
            StressLeaf = -33._dp ! maximum canopy is reached;
            if (GetCrop_CCxAdjusted() < 0._dp) then
                call SetCrop_CCxAdjusted(GetCCiPrev())
            end if

            if (VirtualTimeCC < GetCrop_DaysToSenescence()) then ! mid-season
                if (GetCrop_CCxAdjusted() > 0.97999_dp*CCxSF) then
                    call SetCCiActual(CanopyCoverNoStressSF(&
                            (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                            GetCrop_DaysToGermination(), &
                            GetCrop_DaysToSenescence(), &
                            GetCrop_DaysToHarvest(), &
                            GetCrop_GDDaysToGermination(), &
                            GetCrop_GDDaysToSenescence(), &
                            GetCrop_GDDaysToHarvest(), &
                            CCoTotal, CCxTotal, GetCrop_CGC(), &
                            CDCTotal, GetCrop_GDDCGC(), GDDCDCTotal, &
                            GetSimulation_SumGDD(), GetCrop_ModeCycle(), &
                            GetSimulation_EffectStress_RedCGC(), &
                            GetSimulation_EffectStress_RedCCX()))
                    call SetCrop_CCxAdjusted(GetCCiActual())
                else
                    call SetCCiActual(CanopyCoverNoStressSF(&
                            (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                            GetCrop_DaysToGermination(), &
                            GetCrop_DaysToSenescence(), &
                            GetCrop_DaysToHarvest(), &
                            GetCrop_GDDaysToGermination(), &
                            GetCrop_GDDaysToSenescence(), &
                            GetCrop_GDDaysToHarvest(), &
                            CCoTotal, &
                            (GetCrop_CCxAdjusted() &
                                /(1._dp - GetSimulation_EffectStress_RedCCx() &
                                                                  /100._dp)), &
                            GetCrop_CGC(), CDCTotal, GetCrop_GDDCGC(), &
                            GDDCDCTotal, GetSimulation_SumGDD(), &
                            GetCrop_ModeCycle(), &
                            GetSimulation_EffectStress_RedCGC(), &
                            GetSimulation_EffectStress_RedCCX()))
                end if
                if (GetCCiActual() > CCxSFCD) then
                    call SetCCiActual(CCxSFCD)
                end if
                ! late season
            else
                StressSenescence = undef_int
                ! to avoid display of zero stress in late season
                if (GetCrop_CCxAdjusted() > CCxSFCD) then
                    call SetCrop_CCxAdjusted(CCxSFCD)
                end if
                if (GetCrop_CCxAdjusted() < 0.01_dp) then
                    call SetCCiActual(0._dp)
                else
                    ! calculate CC in late season
                    ! CCibis = CC which canopy declines
                    ! (soil fertility/salinity stress) further in late season
                    CCibis = CCxSF - (GetSimulation_EffectStress_CDecline() &
                                                                  /100._dp) &
                            * (exp(2._dp * log( &
                                (VirtualTimeCC+GetSimulation_DelayedDays()+1._dp) &
                                    - GetCrop_DaysToFullCanopySF())) &
                                /(GetCrop_DaysToSenescence() &
                                  - GetCrop_DaysToFullCanopySF()))
                    if (CCibis < 0._dp) then
                        call SetCCiActual(0._dp)
                    else
                        ! CCiActual = CC with natural senescence in late season
                        Crop_CCxAdjusted_temp = GetCrop_CCxAdjusted()
                        CDCadjusted = GetCDCadjustedNoStressNew(&
                                            CCxTotal,CDCTotal, &
                                            Crop_CCxAdjusted_temp)
                        call SetCrop_CCxAdjusted(Crop_CCxAdjusted_temp)
                        if ((VirtualTimeCC+GetSimulation_DelayedDays()+1) &
                            < (GetCrop_DaysToSenescence() &
                                + LengthCanopyDecline(GetCrop_CCxAdjusted(), &
                                                      CDCadjusted))) then
                            call SetCCiActual(GetCrop_CCxAdjusted() &
                                     * (1._dp - 0.05_dp * (exp(&
                                ((VirtualTimeCC+GetSimulation_DelayedDays()+1) &
                                    - GetCrop_DaysToSenescence()) &
                                * 3.33_dp &
                                * CDCadjusted/(GetCrop_CCxAdjusted() + 2.29_dp)) &
                                                                     - 1._dp)))
                            ! CCiActual becomes CCibis, when canopy decline is more severe
                            if (CCibis < GetCCiActual()) then
                                call SetCCiActual(CCibis)
                            end if
                        else
                            call SetCCiActual(0._dp)
                        end if
                    end if
                    ! late season
                end if
                ! 3. Canopy can no longer develop
                ! (Mid-season (from tFinalCCx) or Late season stage)
            end if
        end if

        ! 4. Canopy senescence due to water stress ?
        if ((VirtualTimeCC < GetCrop_DaysToSenescence()) &
                                ! not yet late season stage
            .or. (TimeSenescence > 0._dp)) then
            ! in late season with ongoing early senesence
            ! (TimeSenescence in days)
            StressSenescence = 0._dp
            WithBeta = .true.
            Crop_pSenAct_temp = GetCrop_pSenAct()
            call AdjustpSenescenceToETo(GetETo(), TimeSenescence, &
                                        WithBeta, Crop_pSenAct_temp)
            call SetCrop_pSenAct(Crop_pSenAct_temp)
            KsRED = 1._dp  ! effect of soil salinity on the
                           ! threshold for senescence
            if (GetSimulation_SWCtopSoilConsidered()) then
                ! top soil is relative wetter than total root zone
                if ((GetRootZoneWC_ZtopAct() &
                        < (GetRootZoneWC_ZtopFC() &
                            - GetCrop_pSenAct()*KsRED &
                            * (GetRootZoneWC_ZtopFC() &
                                - GetRootZoneWC_ZtopWP()))) &
                    .and. (.not.GetSimulation_ProtectedSeedling())) then
                    TheSenescenceON = .true.
                else
                    TheSenescenceON = .false.
                end if
            else
                if ((GetRootZoneWC_Actual() &
                    < (GetRootZoneWC_FC() &
                        - GetCrop_pSenAct()*KsRED &
                        * (GetRootZoneWC_FC() - GetRootZoneWC_WP()))) &
                    .and. (.not. GetSimulation_ProtectedSeedling())) then
                    TheSenescenceON = .true.
                else
                    TheSenescenceON = .false.
                end if
            end if

            if (TheSenescenceON) then
                ! CanopySenescence
                call SetSimulation_EvapLimitON(.true.)
                ! consider withered crop when not yet in late season
                if (abs(TimeSenescence) < epsilon(0._dp)) then
                    call SetCCiTopEarlySen(GetCCiActual())
                    ! CC before canopy decline
                end if
                TimeSenescence = TimeSenescence + 1._dp  ! add 1 day
                call DetermineCDCadjustedWaterStress(CDCadjusted, KsSen)
                if (GetCCiTopEarlySen() < 0.001_dp) then
                    if ((GetSimulation_SumEToStress() &
                        > GetCrop_SumEToDelaySenescence()) &
                        .or. (abs(GetCrop_SumEToDelaySenescence()) &
                              < epsilon(0._dp))) then
                        CCiSen = 0._dp ! no crop anymore
                    else
                        if (CCdormant > GetCrop_CCo()) then
                            CCiSen = GetCrop_CCo() &
                                + (1._dp &
                                  - GetSimulation_SumEToStress() &
                                        /GetCrop_SumEToDelaySenescence()) &
                                * (CCdormant - GetCrop_CCo())
                        else
                            CCiSen = GetCrop_CCo()
                        end if
                    end if
                else
                    if (((TimeSenescence*CDCTotal*3.33_dp) &
                            /(GetCCiTopEarlySen()+2.29_dp) > 100._dp) &
                            ! e power too large and in any case CCisen << 0
                        .or. (GetCCiPrev() &
                                >= 1.05_dp * GetCCiTopEarlySen())) then
                                ! Ln of negative or zero value
                        if ((GetSimulation_SumEToStress() &
                                > GetCrop_SumEToDelaySenescence()) &
                            .or. (abs(GetCrop_SumEToDelaySenescence()) &
                                    < epsilon(0._dp))) then
                            CCiSen = 0._dp ! no crop anymore
                        else
                            if (CCdormant > GetCrop_CCo()) then
                                CCiSen = GetCrop_CCo() &
                                    + (1._dp &
                                        - GetSimulation_SumEToStress() &
                                            /GetCrop_SumEToDelaySenescence()) &
                                    * (CCdormant - GetCrop_CCo())
                            else
                                CCiSen = GetCrop_CCo()
                            end if
                        end if
                    else
                        ! CDC is adjusted to degree of stress
                        ! time required to reach CCiprev with CDCadjusted
                        if (abs(GetCCiTopEarlySen()) < epsilon(0._dp)) then
                            call SetCCiTopEarlySen(epsilon(1._dp))
                        end if
                        if (abs(CDCadjusted) < epsilon(0._dp)) then
                            CDCadjusted = epsilon(1._dp)
                        end if
                        tTemp = (log(1._dp &
                                     + (1._dp &
                                        - GetCCiPrev()/GetCCiTopEarlySen()) &
                                                                    /0.05_dp)) &
                                /(CDCadjusted*3.33_dp &
                                    /(GetCCiTopEarlySen()+2.29_dp))
                        ! add 1 day to tTemp and calculate CCiSen
                        ! with CDCadjusted
                        CCiSen = GetCCiTopEarlySen() &
                                 * (1._dp - 0.05_dp &
                                            * (exp((tTemp+1._dp) &
                                                    *CDCadjusted &
                                                    *3.33_dp &
                                                    /(GetCCiTopEarlySen() &
                                                                  +2.29)) &
                                                -1))
                    end if

                    if (CCiSen < 0._dp) then
                        CCiSen = 0._dp
                    end if
                    if ((GetCrop_SumEToDelaySenescence() > 0._dp) &
                        .and. (GetSimulation_SumEToStress() &
                                <= GetCrop_SumEToDelaySenescence())) then
                        if ((CCiSen < GetCrop_CCo()) &
                            .or. (CCiSen < CCdormant)) then
                            if (CCdormant > GetCrop_CCo()) then
                                CCiSen = GetCrop_CCo() &
                                        + (1._dp &
                                            - GetSimulation_SumEToStress() &
                                              /GetCrop_SumEToDelaySenescence()) &
                                        * (CCdormant - GetCrop_CCo())
                            else
                                CCiSen = GetCrop_CCo()
                            end if
                        end if
                    end if
                end if
                if (VirtualTimeCC < GetCrop_DaysToSenescence()) then
                    ! before late season
                    if (CCiSen > CCxSFCD) then
                        CCiSen = CCxSFCD
                    end if
                    call SetCCiActual(CCiSen)
                    if (GetCCiActual() > GetCCiPrev()) then
                        call SetCCiActual(GetCCiPrev())
                        ! to avoid jump in CC
                    end if
                    ! when CGCadjusted increases as a result of watering
                    call SetCrop_CCxAdjusted(GetCCiActual())
                    if (GetCCiActual() < CCoTotal) then
                        call SetCrop_CCoAdjusted(GetCCiActual())
                    else
                        call SetCrop_CCoAdjusted(CCoTotal)
                    end if
                else
                    ! in late season
                    if (CCiSen < GetCCiActual()) then
                        call SetCCiActual(CCiSen)
                    end if
                end if

                if ((roundc(10000._dp*CCiSen, mold=1) &
                    <= (10000._dp*CCdormant)) &
                    .or. (roundc(10000._dp*CCiSen, mold=1) &
                            <= roundc(10000._dp*GetCrop_CCo(), mold=1))) then
                    call SetSimulation_SumEToStress(&
                            GetSimulation_SumEToStress() + GetETo())
                end if
            else
                ! no water stress, resulting in canopy senescence
                TimeSenescence = 0._dp
                ! No early senescence or back to normal
                StressSenescence = 0._dp
                call SetSimulation_SumEToStress(0._dp)
                if ((VirtualTimeCC > GetCrop_DaysToSenescence()) &
                    .and. (GetCCiActual() > GetCCiPrev())) then
                    ! result of a rewatering in late season of
                    ! an early declining canopy
                    Crop_CCxAdjusted_temp = GetCrop_CCxAdjusted()
                    call GetNewCCxandCDC(GetCCiPrev(), CDCTotal, &
                                         CCxSF, Crop_CCxAdjusted_temp, &
                                         CDCadjusted)
                    call SetCrop_CCxAdjusted(Crop_CCxAdjusted_temp)
                    call SetCCiActual(CanopyCoverNoStressSF(&
                            (VirtualTimeCC+GetSimulation_DelayedDays()+1), &
                            GetCrop_DaysToGermination(), &
                            GetCrop_DaysToSenescence(), &
                            GetCrop_DaysToHarvest(), &
                            GetCrop_GDDaysToGermination(), &
                            GetCrop_GDDaysToSenescence(), &
                            GetCrop_GDDaysToHarvest(), &
                            CCoTotal, &
                            (GetCrop_CCxAdjusted() &
                                /(1._dp - GetSimulation_EffectStress_RedCCx() &
                                                                    /100._dp)), &
                            GetCrop_CGC(), CDCadjusted, &
                            GetCrop_GDDCGC(), GDDCDCTotal, &
                            GetSimulation_SumGDD(), GetCrop_ModeCycle(), &
                            GetSimulation_EffectStress_RedCGC(), &
                            GetSimulation_EffectStress_RedCCX()))
                end if
            end if
        end if

        ! 5. Adjust GetCrop().CCxWithered - required for correction
        ! of Transpiration of dying green canopy
        if (GetCCiActual() > GetCrop_CCxWithered()) then
            call SetCrop_CCxWithered(GetCCiActual())
        end if

        ! 6. correction for late-season stage for rounding off errors
        if (VirtualTimeCC > GetCrop_DaysToSenescence()) then
            if (GetCCiActual() > GetCCiPrev()) then
                call SetCCiActual(GetCCiPrev())
            end if
        end if

        ! 7. no crop as a result of fertiltiy and/or water stress
        if (roundc(1000._dp*GetCCiActual(), mold=1) <= 0) then
            NoMoreCrop = .true.
        end if

        ! test
        TESTVAL = CGCadjusted
    end if


    contains


    subroutine DetermineCGCadjusted(CGCadjusted)
        real(dp), intent(inout) :: CGCadjusted

        real(dp) :: Wrelative
        real(dp) :: KsLeaf
        real(dp) :: SWCeffectiveRootZone, FCeffectiveRootZone, &
                    WPeffectiveRootZone

        ! determine FC and PWP
        if (GetSimulation_SWCtopSoilConsidered()) then
            ! top soil is relative wetter than total root zone
            SWCeffectiveRootZone = GetRootZoneWC_ZtopAct()
            Wrelative = (GetRootZoneWC_ZtopFC() &
                         - GetRootZoneWC_ZtopAct()) &
                            /(GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopWP())
            FCeffectiveRootZone = GetRootZoneWC_ZtopFC()
            WPeffectiveRootZone = GetRootZoneWC_ZtopWP()
        else
            ! total rootzone is wetter than top soil
            SWCeffectiveRootZone = GetRootZoneWC_Actual()
            Wrelative = (GetRootZoneWC_FC() - GetRootZoneWC_Actual()) &
                            /(GetRootZoneWC_FC() - GetRootZoneWC_WP())
            FCeffectiveRootZone = GetRootZoneWC_FC()
            WPeffectiveRootZone = GetRootZoneWC_WP()
        end if

        ! Canopy stress and effect of soil water stress on CGC
        if (SWCeffectiveRootZone >= FCeffectiveRootZone) then
            CGCadjusted = CGCSF
            StressLeaf = 0._dp
        elseif (SWCeffectiveRootZone <= WPeffectiveRootZone) then
            CGCadjusted = 0._dp
            StressLeaf = 100._dp
        else
            if (Wrelative <= GetCrop_pLeafAct()) then
                CGCadjusted = CGCSF
                StressLeaf = 0._dp
            elseif (Wrelative >= pLeafLLAct) then
                CGCadjusted = 0._dp
                StressLeaf = 100._dp
            else
                KsLeaf = KsAny(Wrelative, GetCrop_pLeafAct(), &
                               pLeafLLAct, GetCrop_KsShapeFactorLeaf())
                CGCadjusted = CGCSF * KsLeaf
                StressLeaf = 100._dp * (1._dp - KsLeaf)
            end if
        end if

    end subroutine DetermineCGCadjusted


    subroutine DetermineCDCadjustedWaterStress(CDCadjusted, KsSen)
        real(dp), intent(inout) :: CDCadjusted
        real(dp), intent(inout) :: KsSen

        real(dp) :: Wrelative
        real(dp) :: pSenLL
        real(dp) :: pSenAct
        logical :: WithBeta

        pSenLL = 0.999_dp ! WP
        if (GetSimulation_SWCtopSoilConsidered()) then
        ! top soil is relative wetter than total root zone
            Wrelative = (GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopAct()) &
                        /(GetRootZoneWC_ZtopFC() - GetRootZoneWC_ZtopWP())
                                                                ! top soil
        else
            Wrelative = (GetRootZoneWC_FC() - GetRootZoneWC_Actual()) &
                        /(GetRootZoneWC_FC() - GetRootZoneWC_WP())
                                                 ! total root zone
        end if
        WithBeta = .false.
        call AdjustpSenescenceToETo(GetETo(), TimeSenescence, &
                                    WithBeta, pSenAct)
        if (Wrelative <= pSenAct) then
            CDCadjusted = 0.001_dp ! extreme small decline
            StressSenescence = 0._dp
            KsSen = 1._dp
        elseif (Wrelative >= pSenLL) then
            CDCadjusted = CDCTotal * (CCxSFCD+2.29_dp)/(CCxTotal+2.29_dp)
                                                            ! full speed
            StressSenescence = 100._dp
            KsSen = 0._dp
        else
            KsSen = KsAny(Wrelative, pSenAct, pSenLL, &
                          GetCrop_KsShapeFactorSenescence())
            if (KsSen > ac_zero_threshold) then
                CDCadjusted = CDCTotal &
                              * ((CCxSFCD+2.29_dp)/(CCxTotal+2.29_dp)) &
                                    * (1._dp - exp(8._dp*log(KsSen)))
                StressSenescence = 100._dp * (1._dp - KsSen)
            else
                CDCadjusted = 0._dp
                StressSenescence = 0._dp
            end if
        end if
    end subroutine DetermineCDCadjustedWaterStress


    real(dp) function RequiredTimeNew(CCiToFind, CCo, CCx, CGCadjusted)
        real(dp), intent(in) :: CCiToFind
        real(dp), intent(in) :: CCo
        real(dp), intent(in) :: CCx
        real(dp), intent(in) :: CGCadjusted

        real(dp) :: CGCx

        ! Only when VirtualTime > 1
        ! and CCx < CCiToFind
        ! 1. CGCx to reach CCiToFind on previous day (= VirtualTime -1 )
        if (CCiToFind <= CCx/2._dp) then
            CGCx = (log(CCiToFind/CCo))/VirtualTimeCC
        else
            CGCx = (log((0.25_dp*CCx*CCx/CCo)/(CCx-CCiToFind)))/VirtualTimeCC
        end if
        ! 2. Required time
        RequiredTimeNew = VirtualTimeCC * CGCx/CGCadjusted
    end function RequiredTimeNew


    real(dp) function CCatTime(tfictive, CCoGiven, CGCGiven, CCxGiven)
        real(dp), intent(in) :: tfictive
        real(dp), intent(in) :: CCoGiven
        real(dp), intent(in) :: CGCGiven
        real(dp), intent(in) :: CCxGiven

        real(dp) :: CCi

        CCi = CCoGiven * exp(CGCGiven * tfictive)
        if (CCi > CCxGiven/2._dp) then
            CCi = CCxGiven - 0.25_dp &
                             * (CCxGiven/CCoGiven) &
                             * CCxGiven &
                             * exp(-CGCGiven*tfictive)
        end if
        CCatTime = CCi
    end function CCatTime


    subroutine DetermineCCxAdjusted(CCxAdjusted)
        real(dp), intent(inout) :: CCxAdjusted

        real(dp) :: tfictive

        ! 1. find time (tfictive) required to reach CCiPrev
        !    (CCi of previous day) with CGCadjusted
        tfictive = RequiredTimeNew(GetCCiPrev(), GetCrop_CCoAdjusted(), &
                                   CCxSF, CGCadjusted)

        ! 2. Get CCxadjusted (reached at end of stretched crop development)
        if (tfictive > 0._dp) then
            tfictive = tfictive + (tFinalCCx - VirtualTimeCC)
            CCxAdjusted = CCatTime(tfictive, GetCrop_CCoAdjusted(), &
                                   CGCadjusted, CCxSF)
        else
            CCxAdjusted = undef_double ! this means CCiActual := CCiPrev
        end if
    end subroutine DetermineCCxAdjusted


    subroutine GetNewCCxandCDC(CCiPrev, CDC, CCx, CCxAdjusted, CDCadjusted)
        real(dp), intent(in) :: CCiPrev
        real(dp), intent(in) :: CDC
        real(dp), intent(in) :: CCx
        real(dp), intent(inout) :: CCxAdjusted
        real(dp), intent(inout) :: CDCadjusted

        CCxAdjusted = CCiPrev &
                      /(1._dp - 0.05_dp &
                            * (exp((VirtualTimeCC-GetCrop_DaysToSenescence()) &
                                          *CDC*3.33_dp/(CCX+2.29_dp))-1._dp))
        ! CDCadjusted := CDC * CCxAdjusted/CCx;
        CDCadjusted = CDC * (CCxAdjusted+2.29_dp)/(CCx+2.29_dp)
    end subroutine GetNewCCxandCDC
end subroutine DetermineCCi


subroutine FeedbackCC()

    if (((GetCCiActual() - GetCCiPrev()) > 0.005_dp) &
        ! canopy is still developing
        .and. (GetTact() < epsilon(0._dp))) then
        ! due to aeration stress or ETo = 0
        call SetCCiActual(GetCCiPrev())
        ! no transpiration, no crop developmentc
    end if
end subroutine FeedbackCC


subroutine HorizontalInflowGWTable(DepthGWTmeter, HorizontalSaltFlow, &
                                   HorizontalWaterFlow)
    real(dp), intent(in) :: DepthGWTmeter
    real(dp), intent(inout) :: HorizontalSaltFlow
    real(dp), intent(inout) :: HorizontalWaterFlow

    real(dp) :: Ztot, Zi, DeltaTheta, SaltAct, SaltAdj
    integer(int32) :: compi, celli
    type(CompartmentIndividual) :: Compi_temp

    Ztot = 0._dp
    do compi = 1, GetNrCompartments()
        Ztot = Ztot + GetCompartment_Thickness(compi)
        Zi = Ztot - GetCompartment_Thickness(compi)/2._dp
        if (Zi >= DepthGWTmeter) then
            ! soil water content is at saturation
            if (GetCompartment_Theta(compi) &
                < GetSoilLayer_SAT(GetCompartment_Layer(compi))/100._dp) then
                DeltaTheta = GetSoilLayer_SAT(GetCompartment_Layer(compi))/100._dp &
                            - GetCompartment_Theta(compi)
                call SetCompartment_theta(&
                        compi, &
                        GetSoilLayer_SAT(GetCompartment_Layer(compi))/100._dp)
                HorizontalWaterFlow = HorizontalWaterFlow &
                                      + 1000._dp * DeltaTheta &
                                            * GetCompartment_Thickness(compi) &
                                            * (1._dp &
                      - GetSoilLayer_GravelVol(GetCompartment_Layer(compi)) &
                                                                  /100._dp)
            end if
            ! ECe is equal to the EC of the groundwater table
            if (abs(ECeComp(GetCompartment_i(compi)) - GetECiAqua()) &
                        > 0.0001_dp) then
                SaltAct = 0._dp
                do celli = 1, GetSoilLayer_SCP1(GetCompartment_Layer(compi))
                    SaltAct = SaltAct &
                              + (GetCompartment_Salt(compi, celli) &
                                    + GetCompartment_Depo(compi, celli))&
                                /100._dp ! Mg/ha
                end do
                Compi_temp = GetCompartment_i(compi)
                call DetermineSaltContent(GetECiAqua(), Compi_temp)
                call SetCompartment_i(compi, Compi_temp)
                SaltAdj = 0._dp
                do celli = 1, GetSoilLayer_SCP1(GetCompartment_Layer(compi))
                    SaltAdj = SaltAdj &
                              + (GetCompartment_Salt(compi, celli) &
                                    + GetCompartment_Depo(compi, celli))&
                                 /100._dp ! Mg/ha
                end do
                HorizontalSaltFlow = HorizontalSaltFlow + (SaltAdj - SaltAct)
            end if
        end if
    end do
end subroutine HorizontalInflowGWTable


subroutine BUDGET_module(dayi, TargetTimeVal, TargetDepthVal, VirtualTimeCC, &
                         SumInterval, DayLastCut, NrDayGrow, Tadj, GDDTadj, &
                         GDDayi, CGCref, GDDCGCref, CO2i, CCxTotal, CCoTotal, &
                         CDCTotal, GDDCDCTotal, SumGDDadjCC, Coeffb0Salt, &
                         Coeffb1Salt, Coeffb2Salt, StressTotSaltPrev, &
                         DayFraction, GDDayFraction, FracAssim, &
                         StressSFadjNEW, StorageON, MobilizationON, &
                         StressLeaf, StressSenescence, TimeSenescence, &
                         NoMoreCrop, TESTVAL)
    integer(int32), intent(in) :: dayi
    integer(int32), intent(in) :: TargetTimeVal
    integer(int32), intent(in) :: TargetDepthVal
    integer(int32), intent(in) :: VirtualTimeCC
    integer(int32), intent(in) :: SumInterval
    integer(int32), intent(in) :: DayLastCut
    integer(int32), intent(in) :: NrDayGrow
    integer(int32), intent(in) :: Tadj
    integer(int32), intent(in) :: GDDTadj
    real(dp), intent(in) :: GDDayi
    real(dp), intent(in) :: CGCref
    real(dp), intent(in) :: GDDCGCref
    real(dp), intent(in) :: CO2i
    real(dp), intent(in) :: CCxTotal
    real(dp), intent(in) :: CCoTotal
    real(dp), intent(in) :: CDCTotal
    real(dp), intent(in) :: GDDCDCTotal
    real(dp), intent(in) :: SumGDDadjCC
    real(dp), intent(in) :: Coeffb0Salt
    real(dp), intent(in) :: Coeffb1Salt
    real(dp), intent(in) :: Coeffb2Salt
    real(dp), intent(in) :: StressTotSaltPrev
    real(dp), intent(in) :: DayFraction
    real(dp), intent(in) :: GDDayFraction
    real(dp), intent(in) :: FracAssim
    integer(int32), intent(in) :: StressSFadjNEW
    logical, intent(in) :: StorageON
    logical, intent(in) :: MobilizationON
    real(dp), intent(inout) :: StressLeaf
    real(dp), intent(inout) :: StressSenescence
    real(dp), intent(inout) :: TimeSenescence
    logical, intent(inout) :: NoMoreCrop
    real(dp), intent(inout) :: TESTVAL


    integer(intEnum) ::  control
    real(dp) :: InfiltratedRain, InfiltratedIrrigation, &
                InfiltratedStorage, EpotTot, SubDrain
    integer(int32) :: DAP
    real(dp) :: ECInfilt
        !! EC of the infiltrated water (surface storage)
    logical :: WaterTableInProfile
    real(dp) :: HorizontalWaterFlow, HorizontalSaltFlow
    logical :: SWCtopSoilConsidered_temp
    real(dp) :: EvapWCsurf_temp, CRwater_temp, Tpot_temp, Epot_temp
    type(CompartmentIndividual), dimension(max_No_compartments) :: Comp_temp
    real(dp) :: Crop_pActStom_temp
    real(dp) :: CRsalt_temp, ECdrain_temp, Surf0_temp
    integer(int32) :: TargetTimeVal_loc
    integer(int32) :: StressSFadjNEW_loc

    TargetTimeVal_loc = TargetTimeVal
    StressSFadjNEW_loc = StressSFadjNEW

    ! 1. Soil water balance
    control = control_begin_day
    ECdrain_temp = GetECdrain()
    Surf0_temp = GetSurf0()
    call CheckWaterSaltBalance(dayi, InfiltratedRain, control, &
                               InfiltratedIrrigation, InfiltratedStorage, &
                               Surf0_temp, ECInfilt, ECdrain_temp, &
                               HorizontalWaterFlow, HorizontalSaltFlow, &
                               SubDrain)
    call SetECdrain(ECdrain_temp)
    call SetSurf0(Surf0_temp)

    ! 2. Adjustments in presence of Groundwater table
    call CheckForWaterTableInProfile(GetZiAqua()/100._dp, GetCompartment(), &
                                     WaterTableInProfile)
    Comp_temp = GetCompartment()
    call CalculateAdjustedFC(GetZiAqua()/100._dp, Comp_temp)
    call SetCompartment(Comp_temp)

    ! 3. Drainage
    call calculate_drainage()

    ! 4. Runoff
    if (GetManagement_Bundheight() < 0.001_dp) then
        call SetDaySubmerged(0)
        if ((GetManagement_RunoffON()) .and. (GetRain() > 0.1_dp)) then
            call calculate_runoff(GetSimulParam_RunoffDepth())
        end if
    end if

    ! 5. Infiltration (Rain and Irrigation)
    if ((GetRainRecord_DataType() == datatype_decadely) &
            .or. (GetRainRecord_DataType() == datatype_monthly)) then
        call CalculateEffectiveRainfall(SubDrain)
    end if
    if (((GetIrriMode() == IrriMode_Generate) &
        .and. (GetIrrigation() < epsilon(0._dp))) &
            .and. (TargetTimeVal_loc /= -999)) then
        call Calculate_irrigation(SubDrain, TargetTimeVal_loc, TargetDepthVal)
    end if
    if (GetManagement_Bundheight() >= 0.01_dp) then
        call calculate_surfacestorage(InfiltratedRain, InfiltratedIrrigation, &
                                      InfiltratedStorage, ECinfilt, SubDrain, &
                                      dayi)
    else
        call calculate_Extra_runoff(InfiltratedRain, InfiltratedIrrigation, &
                                    InfiltratedStorage, SubDrain)
    end if
    call calculate_infiltration(InfiltratedRain, InfiltratedIrrigation, &
                                InfiltratedStorage, SubDrain)

    ! 6. Capillary Rise
    CRwater_temp = GetCRwater()
    CRsalt_temp = GetCRsalt()
    call calculate_CapillaryRise(CRwater_temp, CRsalt_temp)
    call SetCRwater(CRwater_temp)
    call SetCRsalt(CRsalt_temp)

    ! 7. Salt balance
    call calculate_saltcontent(InfiltratedRain, InfiltratedIrrigation, &
                               InfiltratedStorage, SubDrain, dayi)


    ! 8. Check Germination
    if ((.not. GetSimulation_Germinate()) .and. (dayi >=GetCrop_Day1())) then
        call CheckGermination()
    end if

    ! 9. Determine effect of soil fertiltiy and soil salinity stress
    if (.not. NoMoreCrop) then
        call EffectSoilFertilitySalinityStress(StressSFadjNEW_loc, Coeffb0Salt, &
                                               Coeffb1Salt, Coeffb2Salt, &
                                               NrDayGrow, StressTotSaltPrev, &
                                               VirtualTimeCC)
    end if


    ! 10. Canopy Cover (CC)
    if (.not. NoMoreCrop) then
        ! determine water stresses affecting canopy cover
        SWCtopSoilConsidered_temp = GetSimulation_SWCtopSoilConsidered()
        call DetermineRootZoneWC(GetRootingDepth(), SWCtopSoilConsidered_temp)
        call SetSimulation_SWCtopSoilConsidered(SWCtopSoilConsidered_temp)
        ! determine canopy cover
        select case (GetCrop_ModeCycle())
            case(modecycle_GDDays)
            call DetermineCCiGDD(CCxTotal, CCoTotal, StressLeaf, FracAssim, &
                                 MobilizationON, StorageON, SumGDDAdjCC, &
                                 VirtualTimeCC, StressSenescence, &
                                 TimeSenescence, NoMoreCrop, CDCTotal, &
                                 GDDayFraction, &
                                 GDDayi, GDDCDCTotal, GDDTadj)
            case default
            call DetermineCCi(CCxTotal, CCoTotal, StressLeaf, FracAssim, &
                              MobilizationON, StorageON, Tadj, VirtualTimeCC, &
                              StressSenescence, TimeSenescence, NoMoreCrop, &
                              CDCTotal, &
                              DayFraction, GDDCDCTotal, TESTVAL)
        end select
    end if

    ! 11. Determine Tpot and Epot
    ! 11.1 Days after Planting
    if (GetCrop_ModeCycle() == modecycle_Calendardays) then
        DAP = VirtualTimeCC
    else
        ! growing degree days - to position correctly where in cycle
        DAP = SumCalendarDays(roundc(SumGDDadjCC, mold=1), GetCrop_Day1(), &
                              GetCrop_Tbase(), GetCrop_Tupper(), &
                              GetSimulParam_Tmin(), GetSimulParam_Tmax())
        DAP = DAP + GetSimulation_DelayedDays()
            ! are not considered when working with GDDays
    end if

    ! 11.2 Calculation
    Tpot_temp = GetTpot()
    call CalculateETpot(DAP, GetCrop_DaysToGermination(), &
                        GetCrop_DaysToFullCanopy(), GetCrop_DaysToSenescence(), &
                        GetCrop_DaysToHarvest(), DayLastCut, GetCCiActual(), &
                        GetETo(), GetCrop_KcTop(), GetCrop_KcDecline(), &
                        GetCrop_CCxAdjusted(), GetCrop_CCxWithered(), &
                        real(GetCrop_CCEffectEvapLate(), kind=dp), CO2i, &
                        GDDayi, GetCrop_GDtranspLow(), Tpot_temp, EpotTot)
    call SetTpot(Tpot_temp)
    call SetEpot(EpotTot)
        ! adjustment Epot for mulch and partial wetting in next step
    Crop_pActStom_temp = GetCrop_pActStom()
    call AdjustpStomatalToETo(GetETo(), Crop_pActStom_temp)
    call SetCrop_pActStom(Crop_pActStom_temp)

    ! 12. Evaporation
    if (.not. GetPreDay()) then
        call PrepareStage2()
            ! Initialize Simulation.EvapstartStg2 (REW is gone)
    end if
    if ((GetRain() > 0._dp) &
        .or. ((GetIrrigation() > 0._dp) &
            .and. (GetIrriMode() /= IrriMode_Inet))) then
        call PrepareStage1()
    end if
    EvapWCsurf_temp = GetSimulation_EvapWCsurf()
    Epot_temp = GetEpot()
    call AdjustEpotMulchWettedSurface(dayi, EpotTot, Epot_temp, EvapWCsurf_temp)
    call SetEpot(Epot_temp)
    call SetSimulation_EvapWCsurf(EvapWCsurf_temp)
    if (((GetRainRecord_DataType() == datatype_Decadely) &
            .or. (GetRainRecord_DataType() == datatype_Monthly)) &
        .and. (GetSimulParam_EffectiveRain_RootNrEvap() > 0)) then
        ! reduction soil evaporation
        call SetEpot(GetEpot() &
                    * (exp((1._dp/GetSimulParam_EffectiveRain_RootNrEvap())&
                            *log((GetSoil_REW()+1._dp)/20._dp))))
    end if
    ! actual evaporation
    call SetEact(0._dp)
    if (GetEpot() > 0._dp) then
        ! surface water
        if (GetSurfaceStorage() > 0._dp) then
            call CalculateEvaporationSurfaceWater
        end if
        ! stage 1 evaporation
        if ((abs(GetEpot() - GetEact()) > ac_zero_threshold) &
            .and. (GetSimulation_EvapWCsurf() > 0._dp)) then
            call CalculateSoilEvaporationStage1()
        end if
        ! stage 2 evaporation
        if (abs(GetEpot() - GetEact()) > ac_zero_threshold) then
            call CalculateSoilEvaporationStage2()
        end if
    end if
    ! Reset redcution Epot for 10-day or monthly rainfall data
    if (((GetRainRecord_DataType() == datatype_Decadely) &
            .or. (GetRainRecord_DataType() == datatype_Monthly)) &
        .and. (GetSimulParam_EffectiveRain_RootNrEvap() > 0._dp)) then
        call SetEpot(GetEpot()&
                    /(exp((1._dp/GetSimulParam_EffectiveRain_RootNrEvap()) &
                           *log((GetSoil_REW()+1._dp)/20._dp))))
    end if


    ! 13. Transpiration
    if ((.not. NoMoreCrop) .and. (GetRootingDepth() > 0.0001_dp)) then
        if ((GetSurfaceStorage() > 0._dp) &
            .and. ((GetCrop_AnaeroPoint() == 0) &
                  .or. (GetDaySubmerged() < GetSimulParam_DelayLowOxygen()))) then
            call surface_transpiration(Coeffb0Salt, Coeffb1Salt, Coeffb2Salt)
        else
            call calculate_transpiration(GetTpot(), Coeffb0Salt, Coeffb1Salt, &
                                         Coeffb2Salt)
        end if
    end if
    if (GetSurfaceStorage() < epsilon(0._dp)) then
        call SetDaySubmerged(0)
    end if
    call FeedbackCC()

    ! 14. Adjustment to groundwater table
    if (WaterTableInProfile) then
        call HorizontalInflowGWTable(GetZiAqua()/100._dp, HorizontalSaltFlow, &
                                     HorizontalWaterFlow)
    end if

    ! 15. Salt concentration
    call ConcentrateSalts()

    ! 16. Soil water balance
    control = control_end_day
    ECdrain_temp = GetECdrain()
    Surf0_temp = GetSurf0()
    call CheckWaterSaltBalance(dayi, InfiltratedRain, control, &
                               InfiltratedIrrigation, InfiltratedStorage, &
                               Surf0_temp, ECInfilt, ECdrain_temp, &
                               HorizontalWaterFlow, HorizontalSaltFlow, &
                               SubDrain)
    call SetECdrain(ECdrain_temp)
    call SetSurf0(Surf0_temp)
end subroutine BUDGET_module

!-----------------------------------------------------------------------------
! end BUDGET_module
!-----------------------------------------------------------------------------

end module ac_simul
