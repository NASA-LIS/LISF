!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !MODULE: LVT_metric_pluginMod
!  \label(LVT_metric_pluginMod)
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  27 Jun 2011;   Sujay Kumar  Initial Specification
! 
!EOP
module LVT_metric_pluginMod

  implicit none
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  PUBLIC :: LVT_metric_plugin  
contains
!BOP
! !ROUTINE: LVT_metric_plugin
!  \label{LVT_metric_plugin}
!
! !DESCRIPTION:
!
!
! !INTERFACE:
  subroutine LVT_metric_plugin

    use LVT_pluginIndices
    use LVT_MEANMod, only : LVT_initMEAN, LVT_diagnoseMEAN, LVT_computeMEAN,&
         LVT_writeMetric_MEAN, LVT_resetMetric_MEAN, LVT_writerestart_MEAN, &
         LVT_readrestart_MEAN

    use LVT_MinMod, only : LVT_initMin, LVT_diagnoseMin, LVT_computeMin,&
         LVT_writeMetric_Min, LVT_resetMetric_Min, LVT_writerestart_Min, &
         LVT_readrestart_Min

    use LVT_MaxMod, only : LVT_initMax, LVT_diagnoseMax, LVT_computeMax,&
         LVT_writeMetric_Max, LVT_resetMetric_Max, LVT_writerestart_Max, &
         LVT_readrestart_Max

    use LVT_SumMod, only : LVT_initSum, LVT_diagnoseSum, LVT_computeSum,&
         LVT_writeMetric_Sum, LVT_resetMetric_Sum, LVT_writerestart_Sum, &
         LVT_readrestart_Sum

    use LVT_StdevMod, only : LVT_initStdev, LVT_diagnoseStdev, LVT_computeStdev,&
         LVT_writeMetric_Stdev, LVT_resetMetric_Stdev, LVT_writerestart_Stdev, &
         LVT_readrestart_Stdev

    use LVT_VarianceMod, only : LVT_initVariance, LVT_diagnoseVariance, LVT_computeVariance,&
         LVT_writeMetric_Variance, LVT_resetMetric_Variance, LVT_writerestart_Variance, &
         LVT_readrestart_Variance

    use LVT_AnomalyMod, only : LVT_initAnomaly, LVT_diagnoseAnomaly, &
         LVT_computeAnomaly, LVT_writeMetric_Anomaly, &
         LVT_resetMetric_Anomaly, LVT_writerestart_Anomaly, &
         LVT_readrestart_Anomaly

    use LVT_RMSEMod, only : LVT_initRMSE, LVT_diagnoseRMSE, LVT_computeRMSE,&
         LVT_writeMetric_RMSE, LVT_resetMetric_RMSE, LVT_writerestart_RMSE, &
         LVT_readrestart_RMSE

    use LVT_BIASMod, only : LVT_initBIAS, LVT_diagnoseBIAS, LVT_computeBIAS,&
         LVT_writeMetric_BIAS, LVT_resetMetric_BIAS, LVT_writerestart_BIAS, &
         LVT_readrestart_BIAS

    use LVT_MAEMod, only : LVT_initMAE, LVT_diagnoseMAE, LVT_computeMAE,&
         LVT_writeMetric_MAE, LVT_resetMetric_MAE, LVT_writerestart_MAE,&
         LVT_readrestart_MAE

    use LVT_PODYMod, only : LVT_initPODY, LVT_diagnosePODY, LVT_computePODY,&
         LVT_writeMetric_PODY, LVT_resetMetric_PODY, LVT_writerestart_PODY,&
         LVT_readrestart_PODY

    use LVT_PODNMod, only : LVT_initPODN, LVT_diagnosePODN, LVT_computePODN,&
         LVT_writeMetric_PODN, LVT_resetMetric_PODN, LVT_writerestart_PODN,&
         LVT_readrestart_PODN

    use LVT_FARMod, only : LVT_initFAR, LVT_diagnoseFAR, LVT_computeFAR,&
         LVT_writeMetric_FAR, LVT_resetMetric_FAR, LVT_writerestart_FAR,&
         LVT_readrestart_FAR

    ! EMK
    use LVT_DFRMod, only : LVT_initDFR, LVT_diagnoseDFR, LVT_computeDFR,&
         LVT_writeMetric_DFR, LVT_resetMetric_DFR, LVT_writerestart_DFR,&
         LVT_readrestart_DFR

    use LVT_POFDMod, only : LVT_initPOFD, LVT_diagnosePOFD, LVT_computePOFD,&
         LVT_writeMetric_POFD, LVT_resetMetric_POFD, LVT_writerestart_POFD,&
         LVT_readrestart_POFD

    use LVT_CSIMod, only : LVT_initCSI, LVT_diagnoseCSI, LVT_computeCSI,&
         LVT_writeMetric_CSI, LVT_resetMetric_CSI, LVT_writerestart_CSI,&
         LVT_readrestart_CSI

    use LVT_ACCMod, only : LVT_initACC, LVT_diagnoseACC, LVT_computeACC,&
         LVT_writeMetric_ACC, LVT_resetMetric_ACC, LVT_writerestart_ACC,&
         LVT_readrestart_ACC

    use LVT_FBIASMod, only : LVT_initFBIAS, LVT_diagnoseFBIAS, &
         LVT_computeFBIAS, LVT_writeMetric_FBIAS, LVT_resetMetric_FBIAS,&
         LVT_writerestart_FBIAS, LVT_readrestart_FBIAS

    ! EMK
    use LVT_EFMod, only : LVT_initEF, LVT_diagnoseEF, &
         LVT_computeEF, LVT_writeMetric_EF, LVT_resetMetric_EF,&
         LVT_writerestart_EF, LVT_readrestart_EF

    ! EMK
    use LVT_FFMod, only : LVT_initFF, LVT_diagnoseFF, &
         LVT_computeFF, LVT_writeMetric_FF, LVT_resetMetric_FF,&
         LVT_writerestart_FF, LVT_readrestart_FF

    ! EMK
    use LVT_HSSMod, only : LVT_initHSS, LVT_diagnoseHSS, &
         LVT_computeHSS, LVT_writeMetric_HSS, LVT_resetMetric_HSS,&
         LVT_writerestart_HSS, LVT_readrestart_HSS

    ! EMK
    use LVT_PSSMod, only : LVT_initPSS, LVT_diagnosePSS, &
         LVT_computePSS, LVT_writeMetric_PSS, LVT_resetMetric_PSS,&
         LVT_writerestart_PSS, LVT_readrestart_PSS

    ! EMK
    use LVT_CSSMod, only : LVT_initCSS, LVT_diagnoseCSS, &
         LVT_computeCSS, LVT_writeMetric_CSS, LVT_resetMetric_CSS,&
         LVT_writerestart_CSS, LVT_readrestart_CSS

    use LVT_ETSMod, only : LVT_initETS, LVT_diagnoseETS, LVT_computeETS,&
         LVT_writeMetric_ETS, LVT_resetMetric_ETS,&
         LVT_writerestart_ETS, LVT_readrestart_ETS

    use LVT_RawCorrMod, only : LVT_initRawCorr, LVT_diagnoseRawCorr, &
         LVT_computeRawCorr, LVT_writeMetric_RawCorr, LVT_resetMetric_RawCorr,&
         LVT_writerestart_RawCorr, LVT_readrestart_RawCorr

    use LVT_RnkCorrMod, only : LVT_initRnkCorr, LVT_diagnoseRnkCorr, &
         LVT_computeRnkCorr, LVT_writeMetric_RnkCorr, LVT_resetMetric_RnkCorr,&
         LVT_writerestart_RnkCorr, LVT_readrestart_RnkCorr

    use LVT_AnomalyRnkCorrMod, only : LVT_initAnomalyRnkCorr, LVT_diagnoseAnomalyRnkCorr, &
         LVT_computeAnomalyRnkCorr, LVT_writeMetric_AnomalyRnkCorr, LVT_resetMetric_AnomalyRnkCorr,&
         LVT_writerestart_AnomalyRnkCorr, LVT_readrestart_AnomalyRnkCorr

    use LVT_AnomalyCorrMod, only : LVT_initAnomalyCorr, &
         LVT_diagnoseAnomalyCorr, LVT_computeAnomalyCorr,&
         LVT_writeMetric_AnomalyCorr, LVT_resetMetric_AnomalyCorr,&
         LVT_writerestart_AnomalyCorr, LVT_readrestart_AnomalyCorr

    use LVT_AnomalyRMSEMod, only : LVT_initAnomalyRMSE, &
         LVT_diagnoseAnomalyRMSE, LVT_computeAnomalyRMSE,&
         LVT_writeMetric_AnomalyRMSE, LVT_resetMetric_AnomalyRMSE, &
         LVT_writerestart_AnomalyRMSE, LVT_readrestart_AnomalyRMSE

    use LVT_NSEMod, only : LVT_initNSE, LVT_diagnoseNSE, LVT_computeNSE,&
         LVT_writeMetric_NSE, LVT_resetMetric_NSE,&
         LVT_writerestart_NSE, LVT_readrestart_NSE


    use LVT_ubRMSEMod, only : LVT_initubRMSE, LVT_diagnoseubRMSE, &
         LVT_computeubRMSE, LVT_writeMetric_ubRMSE, LVT_resetMetric_ubRMSE, &
         LVT_writerestart_ubRMSE, LVT_readrestart_ubRMSE

    use LVT_AreaMetricMod, only : LVT_initAreaMetric, &
         LVT_diagnoseAreaMetric, LVT_computeAreaMetric,&
         LVT_writeMetric_AreaMetric, LVT_resetMetric_AreaMetric,&
         LVT_writerestart_AreaMetric, LVT_readrestart_AreaMetric

    use LVT_waveletStatMod, only : LVT_initwaveletStat, LVT_diagnosewaveletStat, &
         LVT_computewaveletStat, LVT_writeMetric_waveletStat, &
         LVT_resetMetric_waveletStat, LVT_writerestart_waveletStat,&
         LVT_readrestart_waveletStat


    use LVT_HNMod, only : LVT_initHN, LVT_diagnoseHN, &
         LVT_computeHN, LVT_writeMetric_HN, &
         LVT_resetMetric_HN, LVT_writerestart_HN,&
         LVT_readrestart_HN

    use LVT_SPIMod, only : LVT_initSPI, LVT_diagnoseSPI, &
         LVT_computeSPI, LVT_writeMetric_SPI, &
         LVT_resetMetric_SPI, LVT_writerestart_SPI,&
         LVT_readrestart_SPI

    use LVT_SdSIMod, only : LVT_initSdSI, LVT_diagnoseSdSI, &
         LVT_computeSdSI, LVT_writeMetric_SdSI, &
         LVT_resetMetric_SdSI, LVT_writerestart_SdSI,&
         LVT_readrestart_SdSI

    use LVT_SRIMod, only : LVT_initSRI, LVT_diagnoseSRI, &
         LVT_computeSRI, LVT_writeMetric_SRI, &
         LVT_resetMetric_SRI, LVT_writerestart_SRI,&
         LVT_readrestart_SRI

    use LVT_SSWIMod, only : LVT_initSSWI, LVT_diagnoseSSWI, &
         LVT_computeSSWI, LVT_writeMetric_SSWI, &
         LVT_resetMetric_SSWI, LVT_writerestart_SSWI,&
         LVT_readrestart_SSWI

    use LVT_SGWIMod, only : LVT_initSGWI, LVT_diagnoseSGWI, &
         LVT_computeSGWI, LVT_writeMetric_SGWI, &
         LVT_resetMetric_SGWI, LVT_writerestart_SGWI,&
         LVT_readrestart_SGWI

    use LVT_percentileMod, only : LVT_initpercentile, LVT_diagnosepercentile, &
         LVT_computepercentile, LVT_writeMetric_percentile, &
         LVT_resetMetric_percentile, LVT_writerestart_percentile,&
         LVT_readrestart_percentile

    use LVT_RFVMod, only : LVT_initRFV, LVT_diagnoseRFV, &
         LVT_computeRFV, LVT_writeMetric_RFV, &
         LVT_resetMetric_RFV, LVT_writerestart_RFV,&
         LVT_readrestart_RFV

    use LVT_MaxTimeMod, only : LVT_initMaxTime, LVT_diagnoseMaxTime, &
         LVT_computeMaxTime, LVT_writeMetric_MaxTime, LVT_resetMetric_MaxTime, &
         LVT_writerestart_MaxTime, LVT_readrestart_MaxTime

    use LVT_MinTimeMod, only : LVT_initMinTime, LVT_diagnoseMinTime, LVT_computeMinTime,&
         LVT_writeMetric_MinTime, LVT_resetMetric_MinTime, LVT_writerestart_MinTime, &
         LVT_readrestart_MinTime


    use LVT_TendencyMod, only : LVT_initTendency, LVT_diagnoseTendency,&
         LVT_computeTendency,&
         LVT_writeMetric_Tendency, LVT_resetMetric_Tendency, &
         LVT_writerestart_Tendency, LVT_readrestart_Tendency

    use LVT_TendencyCorrMod, only : LVT_initTendencyCorr, LVT_diagnoseTendencyCorr,&
         LVT_computeTendencyCorr,&
         LVT_writeMetric_TendencyCorr, LVT_resetMetric_TendencyCorr,&
         LVT_writerestart_TendencyCorr, LVT_readrestart_TendencyCorr
    
    use LVT_ZscoreMod, only : LVT_initZscore, LVT_diagnoseZscore, LVT_computeZscore,&
         LVT_writeMetric_Zscore, LVT_resetMetric_Zscore, LVT_writerestart_Zscore, &
         LVT_readrestart_Zscore

    use LVT_TrendMod, only : LVT_initTrend, LVT_diagnoseTrend, LVT_computeTrend,&
         LVT_writeMetric_Trend, LVT_resetMetric_Trend, LVT_writerestart_Trend, &
         LVT_readrestart_Trend

    use LVT_MetricEntropyMod, only : LVT_initMetricEntropy, &
         LVT_diagnoseMetricEntropy, LVT_computeMetricEntropy, &
         LVT_writeMetric_MetricEntropy, LVT_resetMetric_MetricEntropy,&
         LVT_writerestart_MetricEntropy, LVT_readrestart_MetricEntropy
         
    use LVT_InformationGainMod, only : LVT_initInformationGain, &
         LVT_diagnoseInformationGain, LVT_computeInformationGain, &
         LVT_writeMetric_InformationGain, LVT_resetMetric_InformationGain,&
         LVT_writerestart_InformationGain, LVT_readrestart_InformationGain

    use LVT_FluctuationComplexityMod, only : LVT_initFluctuationComplexity, &
         LVT_diagnoseFluctuationComplexity, LVT_computeFluctuationComplexity, &
         LVT_writeMetric_FluctuationComplexity, &
         LVT_resetMetric_FluctuationComplexity, &
         LVT_writerestart_FluctuationComplexity, &
         LVT_readrestart_FluctuationComplexity

    use LVT_EffectiveComplexityMod, only : LVT_initEffectiveComplexity, &
         LVT_diagnoseEffectiveComplexity, &
         LVT_computeEffectiveComplexity, &
         LVT_writeMetric_EffectiveComplexity, &
         LVT_resetMetric_EffectiveComplexity, &
         LVT_writerestart_EffectiveComplexity, &
         LVT_readrestart_EffectiveComplexity


    use LVT_TCMod, only : LVT_initTC, LVT_diagnoseTC, &
         LVT_computeTC, LVT_writeMetric_TC, LVT_resetMetric_TC,&
         LVT_writerestart_TC, LVT_readrestart_TC

    use LVT_RELMod, only : LVT_initREL, LVT_diagnoseREL, LVT_computeREL,&
         LVT_writeMetric_REL, LVT_resetMetric_REL, LVT_writerestart_REL, &
         LVT_readrestart_REL

    use LVT_RESMod, only : LVT_initRES, LVT_diagnoseRES, LVT_computeRES,&
         LVT_writeMetric_RES, LVT_resetMetric_RES, LVT_writerestart_RES, &
         LVT_readrestart_RES

    use LVT_VULMod, only : LVT_initVUL, LVT_diagnoseVUL, LVT_computeVUL,&
         LVT_writeMetric_VUL, LVT_resetMetric_VUL, LVT_writerestart_VUL, &
         LVT_readrestart_VUL

    use LVT_KMEANSMod, only : LVT_initKMEANS, LVT_diagnoseKMEANS, LVT_computeKMEANS,&
         LVT_writeMetric_KMEANS, LVT_resetMetric_KMEANS, LVT_writerestart_KMEANS, &
         LVT_readrestart_KMEANS

    ! Tian decomposition of bias...EMK
    use LVT_THBMod, only : LVT_initTHB, LVT_diagnoseTHB, LVT_computeTHB,&
         LVT_writeMetric_THB, LVT_resetMetric_THB, LVT_writerestart_THB, &
         LVT_readrestart_THB

    use LVT_TMBMod, only : LVT_initTMB, LVT_diagnoseTMB, LVT_computeTMB,&
         LVT_writeMetric_TMB, LVT_resetMetric_TMB, LVT_writerestart_TMB, &
         LVT_readrestart_TMB

    use LVT_TFBMod, only : LVT_initTFB, LVT_diagnoseTFB, LVT_computeTFB,&
         LVT_writeMetric_TFB, LVT_resetMetric_TFB, LVT_writerestart_TFB, &
         LVT_readrestart_TFB

    use LVT_InformationEntropyMod, only : LVT_initInformationEntropy, &
         LVT_diagnoseInformationEntropy, LVT_computeInformationEntropy,&
         LVT_writeMetric_InformationEntropy, &
         LVT_resetMetric_InformationEntropy, &
         LVT_writerestart_InformationEntropy, &
         LVT_readrestart_InformationEntropy

    use LVT_ConditionalEntropyMod, only : LVT_initConditionalEntropy, &
         LVT_diagnoseConditionalEntropy, LVT_computeConditionalEntropy,&
         LVT_writeMetric_ConditionalEntropy, &
         LVT_resetMetric_ConditionalEntropy, &
         LVT_writerestart_ConditionalEntropy, &
         LVT_readrestart_ConditionalEntropy


    use LVT_RelativeEntropyMod, only : LVT_initRelativeEntropy, &
         LVT_diagnoseRelativeEntropy, LVT_computeRelativeEntropy,&
         LVT_writeMetric_RelativeEntropy, &
         LVT_resetMetric_RelativeEntropy, &
         LVT_writerestart_RelativeEntropy, &
         LVT_readrestart_RelativeEntropy

    use LVT_JointEntropyMod, only : LVT_initJointEntropy, &
         LVT_diagnoseJointEntropy, LVT_computeJointEntropy,&
         LVT_writeMetric_JointEntropy, &
         LVT_resetMetric_JointEntropy, &
         LVT_writerestart_JointEntropy, &
         LVT_readrestart_JointEntropy


    use LVT_MutualInformationMod, only : LVT_initMutualInformation, &
         LVT_diagnoseMutualInformation, LVT_computeMutualInformation,&
         LVT_writeMetric_MutualInformation, LVT_resetMetric_MutualInformation,&
         LVT_writerestart_MutualInformation, LVT_readrestart_MutualInformation

#if 0 
    use LVT_ensMEANMod, only : LVT_initensMEAN, LVT_diagnoseensMEAN, &
         LVT_computeensMEAN, LVT_writeMetric_ensMEAN, &
         LVT_resetMetric_ensMEAN, LVT_writerestart_ensMEAN,&
         LVT_readrestart_ensMEAN

    use LVT_ensStdevMod, only : LVT_initensStdev, LVT_diagnoseensStdev, &
         LVT_computeensStdev, LVT_writeMetric_ensStdev, &
         LVT_resetMetric_ensStdev, LVT_writerestart_ensStdev,&
         LVT_readrestart_ensStdev

    use LVT_ensSpreadMod, only : LVT_initensSpread, LVT_diagnoseensSpread, &
         LVT_computeensSpread, LVT_writeMetric_ensSpread, &
         LVT_resetMetric_ensSpread, LVT_writerestart_ensSpread,&
         LVT_readrestart_ensSpread

    use LVT_ensLLMod, only : LVT_initensLL, LVT_diagnoseensLL, &
         LVT_computeensLL, LVT_writeMetric_ensLL, &
         LVT_resetMetric_ensLL, LVT_writerestart_ensLL, &
         LVT_readrestart_ensLL

    use LVT_ensXcorrMod, only : LVT_initensXcorr, LVT_diagnoseensXcorr, &
         LVT_computeensXcorr, LVT_writeMetric_ensXcorr, &
         LVT_resetMetric_ensXcorr, LVT_writerestart_ensXcorr, &
         LVT_readrestart_ensXcorr

    use LVT_ensSkillMod, only : LVT_initensSkill, LVT_diagnoseensSkill, &
         LVT_computeensSkill, LVT_writeMetric_ensSkill, &
         LVT_resetMetric_ensSkill, LVT_writerestart_ensSkill,&
         LVT_readrestart_ensSkill

    use LVT_ensMEMod, only : LVT_initensME, LVT_diagnoseensME, &
         LVT_computeensME, LVT_writeMetric_ensME, &
         LVT_resetMetric_ensME, LVT_writerestart_ensME,&
         LVT_readrestart_ensME

    use LVT_ensMeanBiasMod, only : LVT_initensMeanBias, LVT_diagnoseensMeanBias, &
         LVT_computeensMeanBias, LVT_writeMetric_ensMeanBias, &
         LVT_resetMetric_ensMeanBias, LVT_writerestart_ensMeanBias,&
         LVT_readrestart_ensMeanBias

    use LVT_ensPercentileMod, only : LVT_initensPercentile, LVT_diagnoseensPercentile, &
         LVT_computeensPercentile, LVT_writeMetric_ensPercentile, &
         LVT_resetMetric_ensPercentile, LVT_writerestart_ensPercentile,&
         LVT_readrestart_ensPercentile

    use LVT_PSDMod, only : LVT_initPSD, &
         LVT_diagnosePSD, LVT_computePSD, &
         LVT_writeMetric_PSD, LVT_resetMetric_PSD,&
         LVT_writerestart_PSD, LVT_readrestart_PSD

    use LVT_KStestMod, only : LVT_initKStest, &
         LVT_diagnoseKStest, LVT_computeKStest,&
         LVT_writeMetric_KStest, LVT_resetMetric_KStest,&
         LVT_writerestart_KStest, LVT_readrestart_KStest


#endif
!EOP

    ! EMK Start of regular metrics

    call registermetricinit(LVT_MEANid,LVT_initMEAN)
    call registermetricdiagnose(LVT_MEANid, LVT_diagnoseMEAN)
    call registermetriccompute(LVT_MEANid, LVT_computeMEAN)
    call registermetricwriteentry(LVT_MEANid,&
         LVT_writeMetric_MEAN)
    call registermetricreset(LVT_MEANid,LVT_resetMetric_MEAN)
    call registermetricwriterestart(LVT_MEANid,LVT_writerestart_MEAN)
    call registermetricreadrestart(LVT_MEANid,LVT_readrestart_MEAN)

    call registermetricinit(LVT_Minid,LVT_initMin)
    call registermetricdiagnose(LVT_Minid, LVT_diagnoseMin)
    call registermetriccompute(LVT_Minid, LVT_computeMin)
    call registermetricwriteentry(LVT_Minid,&
         LVT_writeMetric_Min)
    call registermetricreset(LVT_Minid,LVT_resetMetric_Min)
    call registermetricwriterestart(LVT_Minid,LVT_writerestart_Min)
    call registermetricreadrestart(LVT_Minid,LVT_readrestart_Min)

    call registermetricinit(LVT_Maxid,LVT_initMax)
    call registermetricdiagnose(LVT_Maxid, LVT_diagnoseMax)
    call registermetriccompute(LVT_Maxid, LVT_computeMax)
    call registermetricwriteentry(LVT_Maxid,&
         LVT_writeMetric_Max)
    call registermetricreset(LVT_Maxid,LVT_resetMetric_Max)
    call registermetricwriterestart(LVT_Maxid,LVT_writerestart_Max)
    call registermetricreadrestart(LVT_Maxid,LVT_readrestart_Max)

    call registermetricinit(LVT_Sumid,LVT_initSum)
    call registermetricdiagnose(LVT_Sumid, LVT_diagnoseSum)
    call registermetriccompute(LVT_Sumid, LVT_computeSum)
    call registermetricwriteentry(LVT_Sumid,&
         LVT_writeMetric_Sum)
    call registermetricreset(LVT_Sumid,LVT_resetMetric_Sum)
    call registermetricwriterestart(LVT_Sumid,LVT_writerestart_Sum)
    call registermetricreadrestart(LVT_Sumid,LVT_readrestart_Sum)

    call registermetricinit(LVT_Stdevid,LVT_initStdev)
    call registermetricdiagnose(LVT_Stdevid, LVT_diagnoseStdev)
    call registermetriccompute(LVT_Stdevid, LVT_computeStdev)
    call registermetricwriteentry(LVT_Stdevid,&
         LVT_writeMetric_Stdev)
    call registermetricreset(LVT_Stdevid,LVT_resetMetric_Stdev)
    call registermetricwriterestart(LVT_Stdevid,LVT_writerestart_Stdev)
    call registermetricreadrestart(LVT_Stdevid,LVT_readrestart_Stdev)

    call registermetricinit(LVT_Varianceid,LVT_initVariance)
    call registermetricdiagnose(LVT_Varianceid, LVT_diagnoseVariance)
    call registermetriccompute(LVT_Varianceid, LVT_computeVariance)
    call registermetricwriteentry(LVT_Varianceid,&
         LVT_writeMetric_Variance)
    call registermetricreset(LVT_Varianceid,LVT_resetMetric_Variance)
    call registermetricwriterestart(LVT_Varianceid,LVT_writerestart_Variance)
    call registermetricreadrestart(LVT_Varianceid,LVT_readrestart_Variance)

    call registermetricinit(LVT_Anomalyid,LVT_initAnomaly)
    call registermetricdiagnose(LVT_Anomalyid, LVT_diagnoseAnomaly)
    call registermetriccompute(LVT_Anomalyid, LVT_computeAnomaly)
    call registermetricwriteentry(LVT_Anomalyid,&
         LVT_writeMetric_Anomaly)
    call registermetricreset(LVT_Anomalyid,LVT_resetMetric_Anomaly)
    call registermetricwriterestart(LVT_Anomalyid,LVT_writerestart_Anomaly)
    call registermetricreadrestart(LVT_Anomalyid,LVT_readrestart_Anomaly)

    call registermetricinit(LVT_RMSEid,LVT_initRMSE)
    call registermetricdiagnose(LVT_RMSEid, LVT_diagnoseRMSE)
    call registermetriccompute(LVT_RMSEid, LVT_computeRMSE)
    call registermetricwriteentry(LVT_RMSEid,&
         LVT_writeMetric_RMSE)
    call registermetricreset(LVT_RMSEid,LVT_resetMetric_RMSE)
    call registermetricwriterestart(LVT_RMSEid,LVT_writerestart_RMSE)
    call registermetricreadrestart(LVT_RMSEid,LVT_readrestart_RMSE)

    call registermetricinit(LVT_BIASid,LVT_initBIAS)
    call registermetricdiagnose(LVT_BIASid, LVT_diagnoseBIAS)
    call registermetriccompute(LVT_BIASid, LVT_computeBIAS)
    call registermetricwriteentry(LVT_BIASid,&
         LVT_writeMetric_BIAS)
    call registermetricreset(LVT_BIASid,LVT_resetMetric_BIAS)
    call registermetricwriterestart(LVT_BIASid,LVT_writerestart_BIAS)
    call registermetricreadrestart(LVT_BIASid,LVT_readrestart_BIAS)

    call registermetricinit(LVT_MAEid,LVT_initMAE)
    call registermetricdiagnose(LVT_MAEid, LVT_diagnoseMAE)
    call registermetriccompute(LVT_MAEid, LVT_computeMAE)
    call registermetricwriteentry(LVT_MAEid,&
         LVT_writeMetric_MAE)
    call registermetricreset(LVT_MAEid,LVT_resetMetric_MAE)
    call registermetricwriterestart(LVT_MAEid,LVT_writerestart_MAE)
    call registermetricreadrestart(LVT_MAEid,LVT_readrestart_MAE)

    call registermetricinit(LVT_PODYid,LVT_initPODY)
    call registermetricdiagnose(LVT_PODYid, LVT_diagnosePODY)
    call registermetriccompute(LVT_PODYid, LVT_computePODY)
    call registermetricwriteentry(LVT_PODYid,&
         LVT_writeMetric_PODY)
    call registermetricreset(LVT_PODYid,LVT_resetMetric_PODY)
    call registermetricwriterestart(LVT_PODYid,LVT_writerestart_PODY)
    call registermetricreadrestart(LVT_PODYid,LVT_readrestart_PODY)

    call registermetricinit(LVT_PODNid,LVT_initPODN)
    call registermetricdiagnose(LVT_PODNid, LVT_diagnosePODN)
    call registermetriccompute(LVT_PODNid, LVT_computePODN)
    call registermetricwriteentry(LVT_PODNid,&
         LVT_writeMetric_PODN)
    call registermetricreset(LVT_PODNid,LVT_resetMetric_PODN)
    call registermetricwriterestart(LVT_PODNid,LVT_writerestart_PODN)
    call registermetricreadrestart(LVT_PODNid,LVT_readrestart_PODN)

    call registermetricinit(LVT_FARid,LVT_initFAR)
    call registermetricdiagnose(LVT_FARid, LVT_diagnoseFAR)
    call registermetriccompute(LVT_FARid, LVT_computeFAR)
    call registermetricwriteentry(LVT_FARid,&
         LVT_writeMetric_FAR)
    call registermetricreset(LVT_FARid,LVT_resetMetric_FAR)
    call registermetricwriterestart(LVT_FARid,LVT_writerestart_FAR)
    call registermetricreadrestart(LVT_FARid,LVT_readrestart_FAR)

    call registermetricinit(LVT_POFDid,LVT_initPOFD)
    call registermetricdiagnose(LVT_POFDid, LVT_diagnosePOFD)
    call registermetriccompute(LVT_POFDid, LVT_computePOFD)
    call registermetricwriteentry(LVT_POFDid,&
         LVT_writeMetric_POFD)
    call registermetricreset(LVT_POFDid,LVT_resetMetric_POFD)
    call registermetricwriterestart(LVT_POFDid,LVT_writerestart_POFD)
    call registermetricreadrestart(LVT_POFDid,LVT_readrestart_POFD)

    call registermetricinit(LVT_CSIid,LVT_initCSI)
    call registermetricdiagnose(LVT_CSIid, LVT_diagnoseCSI)
    call registermetriccompute(LVT_CSIid, LVT_computeCSI)
    call registermetricwriteentry(LVT_CSIid,&
         LVT_writeMetric_CSI)
    call registermetricreset(LVT_CSIid,LVT_resetMetric_CSI)
    call registermetricwriterestart(LVT_CSIid,LVT_writerestart_CSI)
    call registermetricreadrestart(LVT_CSIid,LVT_readrestart_CSI)

    call registermetricinit(LVT_ACCid,LVT_initACC)
    call registermetricdiagnose(LVT_ACCid, LVT_diagnoseACC)
    call registermetriccompute(LVT_ACCid, LVT_computeACC)
    call registermetricwriteentry(LVT_ACCid,&
         LVT_writeMetric_ACC)
    call registermetricreset(LVT_ACCid,LVT_resetMetric_ACC)
    call registermetricwriterestart(LVT_ACCid,LVT_writerestart_ACC)
    call registermetricreadrestart(LVT_ACCid,LVT_readrestart_ACC)

    call registermetricinit(LVT_FBIASid,LVT_initFBIAS)
    call registermetricdiagnose(LVT_FBIASid, LVT_diagnoseFBIAS)
    call registermetriccompute(LVT_FBIASid, LVT_computeFBIAS)
    call registermetricwriteentry(LVT_FBIASid,&
         LVT_writeMetric_FBIAS)
    call registermetricreset(LVT_FBIASid,LVT_resetMetric_FBIAS)
    call registermetricwriterestart(LVT_FBIASid,LVT_writerestart_FBIAS)
    call registermetricreadrestart(LVT_FBIASid,LVT_readrestart_FBIAS)

    call registermetricinit(LVT_ETSid,LVT_initETS)
    call registermetricdiagnose(LVT_ETSid, LVT_diagnoseETS)
    call registermetriccompute(LVT_ETSid, LVT_computeETS)
    call registermetricwriteentry(LVT_ETSid,&
         LVT_writeMetric_ETS)
    call registermetricreset(LVT_ETSid,LVT_resetMetric_ETS)
    call registermetricwriterestart(LVT_ETSid,LVT_writerestart_ETS)
    call registermetricreadrestart(LVT_ETSid,LVT_readrestart_ETS)

    call registermetricinit(LVT_RCORRid,LVT_initRawCorr)
    call registermetricdiagnose(LVT_RCORRid, LVT_diagnoseRawCorr)
    call registermetriccompute(LVT_RCORRid, LVT_computeRawCorr)
    call registermetricwriteentry(LVT_RCORRid,&
         LVT_writeMetric_RawCorr)
    call registermetricreset(LVT_RCORRid,LVT_resetMetric_RawCorr)
    call registermetricwriterestart(LVT_RCORRid,LVT_writerestart_RawCorr)
    call registermetricreadrestart(LVT_RCORRid,LVT_readrestart_RawCorr)

    call registermetricinit(LVT_RNKCORRid,LVT_initRnkCorr)
    call registermetricdiagnose(LVT_RNKCORRid, LVT_diagnoseRnkCorr)
    call registermetriccompute(LVT_RNKCORRid, LVT_computeRnkCorr)
    call registermetricwriteentry(LVT_RNKCORRid,&
         LVT_writeMetric_RnkCorr)
    call registermetricreset(LVT_RNKCORRid,LVT_resetMetric_RnkCorr)
    call registermetricwriterestart(LVT_RNKCORRid,LVT_writerestart_RnkCorr)
    call registermetricreadrestart(LVT_RNKCORRid,LVT_readrestart_RnkCorr)

    call registermetricinit(LVT_ARNKCORRid,LVT_initAnomalyRnkCorr)
    call registermetricdiagnose(LVT_ARNKCORRid, LVT_diagnoseAnomalyRnkCorr)
    call registermetriccompute(LVT_ARNKCORRid, LVT_computeAnomalyRnkCorr)
    call registermetricwriteentry(LVT_ARNKCORRid,&
         LVT_writeMetric_AnomalyRnkCorr)
    call registermetricreset(LVT_ARNKCORRid,LVT_resetMetric_AnomalyRnkCorr)
    call registermetricwriterestart(LVT_ARNKCORRid,LVT_writerestart_AnomalyRnkCorr)
    call registermetricreadrestart(LVT_ARNKCORRid,LVT_readrestart_AnomalyRnkCorr)

    call registermetricinit(LVT_ACORRid,LVT_initAnomalyCorr)
    call registermetricdiagnose(LVT_ACORRid, LVT_diagnoseAnomalyCorr)
    call registermetriccompute(LVT_ACORRid, LVT_computeAnomalyCorr)
    call registermetricwriteentry(LVT_ACORRid,&
         LVT_writeMetric_AnomalyCorr)
    call registermetricreset(LVT_ACORRid,LVT_resetMetric_AnomalyCorr)
    call registermetricwriterestart(LVT_ACORRid,LVT_writerestart_AnomalyCorr)
    call registermetricreadrestart(LVT_ACORRid,LVT_readrestart_AnomalyCorr)

    call registermetricinit(LVT_ARMSEid,LVT_initAnomalyRMSE)
    call registermetricdiagnose(LVT_ARMSEid, LVT_diagnoseAnomalyRMSE)
    call registermetriccompute(LVT_ARMSEid, LVT_computeAnomalyRMSE)
    call registermetricwriteentry(LVT_ARMSEid,&
         LVT_writeMetric_AnomalyRMSE)
    call registermetricreset(LVT_ARMSEid,LVT_resetMetric_AnomalyRMSE)
    call registermetricwriterestart(LVT_ARMSEid,LVT_writerestart_AnomalyRMSE)
    call registermetricreadrestart(LVT_ARMSEid,LVT_readrestart_AnomalyRMSE)

    call registermetricinit(LVT_NSEid,LVT_initNSE)
    call registermetricdiagnose(LVT_NSEid, LVT_diagnoseNSE)
    call registermetriccompute(LVT_NSEid, LVT_computeNSE)
    call registermetricwriteentry(LVT_NSEid,&
         LVT_writeMetric_NSE)
    call registermetricreset(LVT_NSEid,LVT_resetMetric_NSE)
    call registermetricwriterestart(LVT_NSEid,LVT_writerestart_NSE)
    call registermetricreadrestart(LVT_NSEid,LVT_readrestart_NSE)

    call registermetricinit(LVT_ubRMSEid,LVT_initubRMSE)
    call registermetricdiagnose(LVT_ubRMSEid, LVT_diagnoseubRMSE)
    call registermetriccompute(LVT_ubRMSEid, LVT_computeubRMSE)
    call registermetricwriteentry(LVT_ubRMSEid,&
         LVT_writeMetric_ubRMSE)
    call registermetricreset(LVT_ubRMSEid,LVT_resetMetric_ubRMSE)
    call registermetricwriterestart(LVT_ubRMSEid,LVT_writerestart_ubRMSE)
    call registermetricreadrestart(LVT_ubRMSEid,LVT_readrestart_ubRMSE)

    call registermetricinit(LVT_AREAid,LVT_initAreaMetric)
    call registermetricdiagnose(LVT_AREAid, LVT_diagnoseAreaMetric)
    call registermetriccompute(LVT_AREAid, LVT_computeAreaMetric)
    call registermetricwriteentry(LVT_AREAid,&
         LVT_writeMetric_AreaMetric)
    call registermetricreset(LVT_AREAid,LVT_resetMetric_AreaMetric)
    call registermetricwriterestart(LVT_AREAid,LVT_writerestart_AreaMetric)
    call registermetricreadrestart(LVT_AREAid,LVT_readrestart_AreaMetric)

    call registermetricinit(LVT_waveletstatid,LVT_initwaveletStat)
    call registermetricdiagnose(LVT_waveletstatid, LVT_diagnosewaveletStat)
    call registermetriccompute(LVT_waveletstatid, LVT_computewaveletStat)
    call registermetricwriteentry(LVT_waveletstatid,&
         LVT_writeMetric_waveletStat)
    call registermetricreset(LVT_waveletstatid,LVT_resetMetric_waveletStat)
    call registermetricwriterestart(LVT_waveletstatid,&
         LVT_writerestart_waveletStat)
    call registermetricreadrestart(LVT_waveletstatid,&
         LVT_readrestart_waveletStat)

    call registermetricinit(LVT_hnid,LVT_initHN)
    call registermetricdiagnose(LVT_hnid, LVT_diagnoseHN)
    call registermetriccompute(LVT_hnid, LVT_computeHN)
    call registermetricwriteentry(LVT_hnid, LVT_writeMetric_HN)
    call registermetricreset(LVT_hnid,LVT_resetMetric_HN)
    call registermetricwriterestart(LVT_hnid,LVT_writerestart_HN)
    call registermetricreadrestart(LVT_hnid, LVT_readrestart_HN)

    call registermetricinit(LVT_spiid,LVT_initSPI)
    call registermetricdiagnose(LVT_spiid, LVT_diagnoseSPI)
    call registermetriccompute(LVT_spiid, LVT_computeSPI)
    call registermetricwriteentry(LVT_spiid,&
         LVT_writeMetric_SPI)
    call registermetricreset(LVT_spiid,LVT_resetMetric_SPI)
    call registermetricwriterestart(LVT_spiid,LVT_writerestart_SPI)
    call registermetricreadrestart(LVT_spiid, LVT_readrestart_SPI)

    call registermetricinit(LVT_sriid,LVT_initSRI)
    call registermetricdiagnose(LVT_sriid, LVT_diagnoseSRI)
    call registermetriccompute(LVT_sriid, LVT_computeSRI)
    call registermetricwriteentry(LVT_sriid,&
         LVT_writeMetric_SRI)
    call registermetricreset(LVT_sriid,LVT_resetMetric_SRI)
    call registermetricwriterestart(LVT_sriid,LVT_writerestart_SRI)
    call registermetricreadrestart(LVT_sriid, LVT_readrestart_SRI)

    call registermetricinit(LVT_sswiid,LVT_initSSWI)
    call registermetricdiagnose(LVT_sswiid, LVT_diagnoseSSWI)
    call registermetriccompute(LVT_sswiid, LVT_computeSSWI)
    call registermetricwriteentry(LVT_sswiid,&
         LVT_writeMetric_SSWI)
    call registermetricreset(LVT_sswiid,LVT_resetMetric_SSWI)
    call registermetricwriterestart(LVT_sswiid,LVT_writerestart_SSWI)
    call registermetricreadrestart(LVT_sswiid, LVT_readrestart_SSWI)

    call registermetricinit(LVT_sgwiid,LVT_initSGWI)
    call registermetricdiagnose(LVT_sgwiid, LVT_diagnoseSGWI)
    call registermetriccompute(LVT_sgwiid, LVT_computeSGWI)
    call registermetricwriteentry(LVT_sgwiid,&
         LVT_writeMetric_SGWI)
    call registermetricreset(LVT_sgwiid,LVT_resetMetric_SGWI)
    call registermetricwriterestart(LVT_sgwiid,LVT_writerestart_SGWI)
    call registermetricreadrestart(LVT_sgwiid, LVT_readrestart_SGWI)

    call registermetricinit(LVT_percentileid,LVT_initpercentile)
    call registermetricdiagnose(LVT_percentileid, LVT_diagnosepercentile)
    call registermetriccompute(LVT_percentileid, LVT_computepercentile)
    call registermetricwriteentry(LVT_percentileid,&
         LVT_writeMetric_percentile)
    call registermetricreset(LVT_percentileid,LVT_resetMetric_percentile)
    call registermetricwriterestart(LVT_percentileid,LVT_writerestart_percentile)
    call registermetricreadrestart(LVT_percentileid, LVT_readrestart_percentile)

    call registermetricinit(LVT_RFVid,LVT_initRFV)
    call registermetricdiagnose(LVT_RFVid, LVT_diagnoseRFV)
    call registermetriccompute(LVT_RFVid, LVT_computeRFV)
    call registermetricwriteentry(LVT_RFVid,&
         LVT_writeMetric_RFV)
    call registermetricreset(LVT_RFVid,LVT_resetMetric_RFV)
    call registermetricwriterestart(LVT_RFVid,LVT_writerestart_RFV)
    call registermetricreadrestart(LVT_RFVid, LVT_readrestart_RFV)

    call registermetricinit(LVT_MinTimeid,LVT_initMinTime)
    call registermetricdiagnose(LVT_MinTimeid, LVT_diagnoseMinTime)
    call registermetriccompute(LVT_MinTimeid, LVT_computeMinTime)
    call registermetricwriteentry(LVT_MinTimeid,&
         LVT_writeMetric_MinTime)
    call registermetricreset(LVT_MinTimeid,LVT_resetMetric_MinTime)
    call registermetricwriterestart(LVT_MinTimeid,LVT_writerestart_MinTime)
    call registermetricreadrestart(LVT_MinTimeid,LVT_readrestart_MinTime)

    call registermetricinit(LVT_MaxTimeid,LVT_initMaxTime)
    call registermetricdiagnose(LVT_MaxTimeid, LVT_diagnoseMaxTime)
    call registermetriccompute(LVT_MaxTimeid, LVT_computeMaxTime)
    call registermetricwriteentry(LVT_MaxTimeid,&
         LVT_writeMetric_MaxTime)
    call registermetricreset(LVT_MaxTimeid,LVT_resetMetric_MaxTime)
    call registermetricwriterestart(LVT_MaxTimeid,LVT_writerestart_MaxTime)
    call registermetricreadrestart(LVT_MaxTimeid,LVT_readrestart_MaxTime)

    call registermetricinit(LVT_Tendencyid,LVT_initTendency)
    call registermetricdiagnose(LVT_Tendencyid, LVT_diagnoseTendency)
    call registermetriccompute(LVT_Tendencyid, LVT_computeTendency)
    call registermetricwriteentry(LVT_Tendencyid,&
         LVT_writeMetric_Tendency)
    call registermetricreset(LVT_Tendencyid,LVT_resetMetric_Tendency)
    call registermetricwriterestart(LVT_Tendencyid,LVT_writerestart_Tendency)
    call registermetricreadrestart(LVT_Tendencyid,LVT_readrestart_Tendency)

    call registermetricinit(LVT_TendencyCorrid,LVT_initTendencyCorr)
    call registermetricdiagnose(LVT_TendencyCorrid, LVT_diagnoseTendencyCorr)
    call registermetriccompute(LVT_TendencyCorrid, LVT_computeTendencyCorr)
    call registermetricwriteentry(LVT_TendencyCorrid,&
         LVT_writeMetric_TendencyCorr)
    call registermetricreset(LVT_TendencyCorrid,LVT_resetMetric_TendencyCorr)
    call registermetricwriterestart(LVT_TendencyCorrid,LVT_writerestart_TendencyCorr)
    call registermetricreadrestart(LVT_TendencyCorrid,LVT_readrestart_TendencyCorr)

    call registermetricinit(LVT_Zscoreid,LVT_initZscore)
    call registermetricdiagnose(LVT_Zscoreid, LVT_diagnoseZscore)
    call registermetriccompute(LVT_Zscoreid, LVT_computeZscore)
    call registermetricwriteentry(LVT_Zscoreid,&
         LVT_writeMetric_Zscore)
    call registermetricreset(LVT_Zscoreid,LVT_resetMetric_Zscore)
    call registermetricwriterestart(LVT_Zscoreid,LVT_writerestart_Zscore)
    call registermetricreadrestart(LVT_Zscoreid,LVT_readrestart_Zscore)

    call registermetricinit(LVT_Trendid,LVT_initTrend)
    call registermetricdiagnose(LVT_Trendid, LVT_diagnoseTrend)
    call registermetriccompute(LVT_Trendid, LVT_computeTrend)
    call registermetricwriteentry(LVT_Trendid,&
         LVT_writeMetric_Trend)
    call registermetricreset(LVT_Trendid,LVT_resetMetric_Trend)
    call registermetricwriterestart(LVT_Trendid,LVT_writerestart_Trend)
    call registermetricreadrestart(LVT_Trendid,LVT_readrestart_Trend)

    call registermetricinit(LVT_SdSIid,LVT_initSdSI)
    call registermetricdiagnose(LVT_SdSIid, LVT_diagnoseSdSI)
    call registermetriccompute(LVT_SdSIid, LVT_computeSdSI)
    call registermetricwriteentry(LVT_SdSIid,&
         LVT_writeMetric_SdSI)
    call registermetricreset(LVT_SdSIid,LVT_resetMetric_SdSI)
    call registermetricwriterestart(LVT_SdSIid,LVT_writerestart_SdSI)
    call registermetricreadrestart(LVT_SdSIid, LVT_readrestart_SdSI)

    call registermetricinit(LVT_TCid,LVT_initTC)
    call registermetricdiagnose(LVT_TCid, LVT_diagnoseTC)
    call registermetriccompute(LVT_TCid, LVT_computeTC)
    call registermetricwriteentry(LVT_TCid,&
         LVT_writeMetric_TC)
    call registermetricreset(LVT_TCid,LVT_resetMetric_TC)
    call registermetricwriterestart(LVT_TCid,LVT_writerestart_TC)
    call registermetricreadrestart(LVT_TCid,LVT_readrestart_TC)

    ! EMK
    call registermetricinit(LVT_DFRid,LVT_initDFR)
    call registermetricdiagnose(LVT_DFRid, LVT_diagnoseDFR)
    call registermetriccompute(LVT_DFRid, LVT_computeDFR)
    call registermetricwriteentry(LVT_DFRid,LVT_writeMetric_DFR)
    call registermetricreset(LVT_DFRid,LVT_resetMetric_DFR)
    call registermetricwriterestart(LVT_DFRid,LVT_writerestart_DFR)
    call registermetricreadrestart(LVT_DFRid,LVT_readrestart_DFR)

    ! EMK
    call registermetricinit(LVT_EFid,LVT_initEF)
    call registermetricdiagnose(LVT_EFid, LVT_diagnoseEF)
    call registermetriccompute(LVT_EFid, LVT_computeEF)
    call registermetricwriteentry(LVT_EFid,LVT_writeMetric_EF)
    call registermetricreset(LVT_EFid,LVT_resetMetric_EF)
    call registermetricwriterestart(LVT_EFid,LVT_writerestart_EF)
    call registermetricreadrestart(LVT_EFid,LVT_readrestart_EF)

    ! EMK
    call registermetricinit(LVT_FFid,LVT_initFF)
    call registermetricdiagnose(LVT_FFid, LVT_diagnoseFF)
    call registermetriccompute(LVT_FFid, LVT_computeFF)
    call registermetricwriteentry(LVT_FFid,LVT_writeMetric_FF)
    call registermetricreset(LVT_FFid,LVT_resetMetric_FF)
    call registermetricwriterestart(LVT_FFid,LVT_writerestart_FF)
    call registermetricreadrestart(LVT_FFid,LVT_readrestart_FF)

    ! EMK
    call registermetricinit(LVT_HSSid,LVT_initHSS)
    call registermetricdiagnose(LVT_HSSid, LVT_diagnoseHSS)
    call registermetriccompute(LVT_HSSid, LVT_computeHSS)
    call registermetricwriteentry(LVT_HSSid,LVT_writeMetric_HSS)
    call registermetricreset(LVT_HSSid,LVT_resetMetric_HSS)
    call registermetricwriterestart(LVT_HSSid,LVT_writerestart_HSS)
    call registermetricreadrestart(LVT_HSSid,LVT_readrestart_HSS)

    ! EMK
    call registermetricinit(LVT_PSSid,LVT_initPSS)
    call registermetricdiagnose(LVT_PSSid, LVT_diagnosePSS)
    call registermetriccompute(LVT_PSSid, LVT_computePSS)
    call registermetricwriteentry(LVT_PSSid,LVT_writeMetric_PSS)
    call registermetricreset(LVT_PSSid,LVT_resetMetric_PSS)
    call registermetricwriterestart(LVT_PSSid,LVT_writerestart_PSS)
    call registermetricreadrestart(LVT_PSSid,LVT_readrestart_PSS)

    ! EMK
    call registermetricinit(LVT_CSSid,LVT_initCSS)
    call registermetricdiagnose(LVT_CSSid, LVT_diagnoseCSS)
    call registermetriccompute(LVT_CSSid, LVT_computeCSS)
    call registermetricwriteentry(LVT_CSSid,LVT_writeMetric_CSS)
    call registermetricreset(LVT_CSSid,LVT_resetMetric_CSS)
    call registermetricwriterestart(LVT_CSSid,LVT_writerestart_CSS)
    call registermetricreadrestart(LVT_CSSid,LVT_readrestart_CSS)

    call registermetricinit(LVT_RELid,LVT_initREL)
    call registermetricdiagnose(LVT_RELid, LVT_diagnoseREL)
    call registermetriccompute(LVT_RELid, LVT_computeREL)
    call registermetricwriteentry(LVT_RELid,&
         LVT_writeMetric_REL)
    call registermetricreset(LVT_RELid,LVT_resetMetric_REL)
    call registermetricwriterestart(LVT_RELid,LVT_writerestart_REL)
    call registermetricreadrestart(LVT_RELid,LVT_readrestart_REL)

    call registermetricinit(LVT_RESid,LVT_initRES)
    call registermetricdiagnose(LVT_RESid, LVT_diagnoseRES)
    call registermetriccompute(LVT_RESid, LVT_computeRES)
    call registermetricwriteentry(LVT_RESid,&
         LVT_writeMetric_RES)
    call registermetricreset(LVT_RESid,LVT_resetMetric_RES)
    call registermetricwriterestart(LVT_RESid,LVT_writerestart_RES)
    call registermetricreadrestart(LVT_RESid,LVT_readrestart_RES)

    call registermetricinit(LVT_VULid,LVT_initVUL)
    call registermetricdiagnose(LVT_VULid, LVT_diagnoseVUL)
    call registermetriccompute(LVT_VULid, LVT_computeVUL)
    call registermetricwriteentry(LVT_VULid,&
         LVT_writeMetric_VUL)
    call registermetricreset(LVT_VULid,LVT_resetMetric_VUL)
    call registermetricwriterestart(LVT_VULid,LVT_writerestart_VUL)
    call registermetricreadrestart(LVT_VULid,LVT_readrestart_VUL)

    call registermetricinit(LVT_KMEANSid,LVT_initKMEANS)
    call registermetricdiagnose(LVT_KMEANSid, LVT_diagnoseKMEANS)
    call registermetriccompute(LVT_KMEANSid, LVT_computeKMEANS)
    call registermetricwriteentry(LVT_KMEANSid,&
         LVT_writeMetric_KMEANS)
    call registermetricreset(LVT_KMEANSid,LVT_resetMetric_KMEANS)
    call registermetricwriterestart(LVT_KMEANSid,LVT_writerestart_KMEANS)
    call registermetricreadrestart(LVT_KMEANSid,LVT_readrestart_KMEANS)

    ! Tian decomposition of bias...EMK
    call registermetricinit(LVT_THBid,LVT_initTHB)
    call registermetricdiagnose(LVT_THBid, LVT_diagnoseTHB)
    call registermetriccompute(LVT_THBid, LVT_computeTHB)
    call registermetricwriteentry(LVT_THBid,&
         LVT_writeMetric_THB)
    call registermetricreset(LVT_THBid,LVT_resetMetric_THB)
    call registermetricwriterestart(LVT_THBid,LVT_writerestart_THB)
    call registermetricreadrestart(LVT_THBid,LVT_readrestart_THB)

    call registermetricinit(LVT_TMBid,LVT_initTMB)
    call registermetricdiagnose(LVT_TMBid, LVT_diagnoseTMB)
    call registermetriccompute(LVT_TMBid, LVT_computeTMB)
    call registermetricwriteentry(LVT_TMBid,&
         LVT_writeMetric_TMB)
    call registermetricreset(LVT_TMBid,LVT_resetMetric_TMB)
    call registermetricwriterestart(LVT_TMBid,LVT_writerestart_TMB)
    call registermetricreadrestart(LVT_TMBid,LVT_readrestart_TMB)

    call registermetricinit(LVT_TFBid,LVT_initTFB)
    call registermetricdiagnose(LVT_TFBid, LVT_diagnoseTFB)
    call registermetriccompute(LVT_TFBid, LVT_computeTFB)
    call registermetricwriteentry(LVT_TFBid,&
         LVT_writeMetric_TFB)
    call registermetricreset(LVT_TFBid,LVT_resetMetric_TFB)
    call registermetricwriterestart(LVT_TFBid,LVT_writerestart_TFB)
    call registermetricreadrestart(LVT_TFBid,LVT_readrestart_TFB)

    ! EMK End of regular metrics

    ! EMK Start of information content metrics

    call registermetricinit(LVT_mentropyid,LVT_initMetricEntropy)
    call registermetricdiagnose(LVT_mentropyid, LVT_diagnoseMetricEntropy)
    call registermetriccompute(LVT_mentropyid, LVT_computeMetricEntropy)
    call registermetricwriteentry(LVT_mentropyid,&
         LVT_writeMetric_MetricEntropy)
    call registermetricreset(LVT_mentropyid,LVT_resetMetric_MetricEntropy)
    call registermetricwriterestart(LVT_mentropyid,LVT_writerestart_MetricEntropy)
    call registermetricreadrestart(LVT_mentropyid,LVT_readrestart_MetricEntropy)

    call registermetricinit(LVT_igainid,LVT_initInformationGain)
    call registermetricdiagnose(LVT_igainid, LVT_diagnoseInformationGain)
    call registermetriccompute(LVT_igainid, LVT_computeInformationGain)
    call registermetricwriteentry(LVT_igainid,&
         LVT_writeMetric_InformationGain)
    call registermetricreset(LVT_igainid,LVT_resetMetric_InformationGain)
    call registermetricwriterestart(LVT_igainid,LVT_writerestart_InformationGain)
    call registermetricreadrestart(LVT_igainid,LVT_readrestart_InformationGain)

    call registermetricinit(LVT_fcomplexityid,LVT_initFluctuationComplexity)
    call registermetricdiagnose(LVT_fcomplexityid, &
         LVT_diagnoseFluctuationComplexity)
    call registermetriccompute(LVT_fcomplexityid, &
         LVT_computeFluctuationComplexity)
    call registermetricwriteentry(LVT_fcomplexityid,&
         LVT_writeMetric_FluctuationComplexity)
    call registermetricreset(LVT_fcomplexityid, &
         LVT_resetMetric_FluctuationComplexity)
    call registermetricwriterestart(LVT_fcomplexityid,&
         LVT_writerestart_FluctuationComplexity)
    call registermetricreadrestart(LVT_fcomplexityid,&
         LVT_readrestart_FluctuationComplexity)

    call registermetricinit(LVT_ecomplexityid,LVT_initEffectiveComplexity)
    call registermetricdiagnose(LVT_ecomplexityid, LVT_diagnoseEffectiveComplexity)
    call registermetriccompute(LVT_ecomplexityid, LVT_computeEffectiveComplexity)
    call registermetricwriteentry(LVT_ecomplexityid,&
         LVT_writeMetric_EffectiveComplexity)
    call registermetricreset(LVT_ecomplexityid,LVT_resetMetric_EffectiveComplexity)
    call registermetricwriterestart(LVT_ecomplexityid,&
         LVT_writerestart_EffectiveComplexity)
    call registermetricreadrestart(LVT_ecomplexityid,&
         LVT_readrestart_EffectiveComplexity)

    call registermetricinit(LVT_IEid,LVT_initInformationEntropy)
    call registermetricdiagnose(LVT_IEid, LVT_diagnoseInformationEntropy)
    call registermetriccompute(LVT_IEid, LVT_computeInformationEntropy)
    call registermetricwriteentry(LVT_IEid,&
         LVT_writeMetric_InformationEntropy)
    call registermetricreset(LVT_IEid,LVT_resetMetric_InformationEntropy)
    call registermetricwriterestart(LVT_IEid,LVT_writerestart_InformationEntropy)
    call registermetricreadrestart(LVT_IEid,LVT_readrestart_InformationEntropy)

    call registermetricinit(LVT_CEid,LVT_initConditionalEntropy)
    call registermetricdiagnose(LVT_CEid, LVT_diagnoseConditionalEntropy)
    call registermetriccompute(LVT_CEid, LVT_computeConditionalEntropy)
    call registermetricwriteentry(LVT_CEid,&
         LVT_writeMetric_ConditionalEntropy)
    call registermetricreset(LVT_CEid,LVT_resetMetric_ConditionalEntropy)
    call registermetricwriterestart(LVT_CEid,LVT_writerestart_ConditionalEntropy)
    call registermetricreadrestart(LVT_CEid,LVT_readrestart_ConditionalEntropy)


    call registermetricinit(LVT_REid,LVT_initRelativeEntropy)
    call registermetricdiagnose(LVT_REid, LVT_diagnoseRelativeEntropy)
    call registermetriccompute(LVT_REid, LVT_computeRelativeEntropy)
    call registermetricwriteentry(LVT_REid,&
         LVT_writeMetric_RelativeEntropy)
    call registermetricreset(LVT_REid,LVT_resetMetric_RelativeEntropy)
    call registermetricwriterestart(LVT_REid,LVT_writerestart_RelativeEntropy)
    call registermetricreadrestart(LVT_REid,LVT_readrestart_RelativeEntropy)

    call registermetricinit(LVT_JEid,LVT_initJointEntropy)
    call registermetricdiagnose(LVT_JEid, LVT_diagnoseJointEntropy)
    call registermetriccompute(LVT_JEid, LVT_computeJointEntropy)
    call registermetricwriteentry(LVT_JEid,&
         LVT_writeMetric_JointEntropy)
    call registermetricreset(LVT_JEid,LVT_resetMetric_JointEntropy)
    call registermetricwriterestart(LVT_JEid,LVT_writerestart_JointEntropy)
    call registermetricreadrestart(LVT_JEid,LVT_readrestart_JointEntropy)


    call registermetricinit(LVT_miid,LVT_initMutualInformation)
    call registermetricdiagnose(LVT_miid, LVT_diagnoseMutualInformation)
    call registermetriccompute(LVT_miid, LVT_computeMutualInformation)
    call registermetricwriteentry(LVT_miid,&
         LVT_writeMetric_MutualInformation)
    call registermetricreset(LVT_miid,LVT_resetMetric_MutualInformation)
    call registermetricwriterestart(LVT_miid,&
         LVT_writerestart_MutualInformation)
    call registermetricreadrestart(LVT_miid,&
         LVT_readrestart_MutualInformation)

    ! EMK End of information content metrics
  end subroutine LVT_metric_plugin
end module LVT_metric_pluginMod
