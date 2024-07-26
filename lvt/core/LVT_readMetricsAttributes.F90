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
! !ROUTINE: LVT_readMetricsAttributes
!  \label{LVT_readMetricsAttributes}
!
! !INTERFACE: 
subroutine LVT_readMetricsAttributes(attribFile)
! 
! !USES: 
  use ESMF
  use LVT_statsDataMod, only : LVT_metrics
  use LVT_logMod,       only : LVT_verify, LVT_abort

  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in) :: attribFile
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the attributes of different analysis
!  metrics within LVT. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  type(ESMF_Config)            :: attribConfig
  integer                      :: rc
  character*100                  :: messages ( 20 )
  
  messages(:) = ''

  attribConfig = ESMF_ConfigCreate(rc=rc)
  call ESMF_ConfigLoadFile(attribConfig,trim(attribFile),rc=rc)
  call LVT_verify(rc,'loading file '//trim(attribFile)//' failed')

  call ESMF_ConfigFindLabel(attribConfig,"Mean:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%mean,"MEAN",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Min:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%min,"MIN",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Max:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%max,"MAX",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Sum:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%sum,"SUM",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standard deviation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%stdev,"STDEV",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Variance:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%variance,"VARIANCE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Anomaly:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%anomaly,"Anomaly",rc)

  call ESMF_ConfigFindLabel(attribConfig,"RMSE:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%rmse,"RMSE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Bias:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%bias,"BIAS",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Mean absolute error:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%mae,"MAE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Probability of detection (PODy):",&
       rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%pody,"PODY",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Probability of detection (PODn):",&
       rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%podn,"PODN",rc)

  call ESMF_ConfigFindLabel(attribConfig,"False alarm ratio (FAR):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%far,"FAR",rc)

  call ESMF_ConfigFindLabel(attribConfig, &
       "Probability of false detection (POFD):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%pofd,"POFD",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Critical success index (CSI):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%csi,"CSI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Accuracy measure (ACC):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%acc,"ACC",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Frequency bias (FBIAS):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%fbias,"FBIAS",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Equitable threat score (ETS):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ets,"ETS",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Raw correlation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%rcorr,"RCORR",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Rank correlation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%rnkcorr,"RNKCORR",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Anomaly rank correlation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%arnkcorr,"ARNKCORR",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Anomaly correlation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%acorr,"ACORR",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Anomaly RMSE:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%armse,"ARMSE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Nash sutcliffe efficiency:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%nse,"NSE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"ubRMSE:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ubrmse,"ubRMSE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Area metric:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%area,"AREA",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Wavelet stat:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%waveletStat, &
       "Waveletstat",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Hausdorff norm:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%hn,"Hnorm",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standard precipitation index:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%spi,"SPI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standard runoff index:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%sri,"SRI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standardized soil water index:", &
       rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%sswi,"SSWI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standardized ground water index:", &
       rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%sgwi,"SGWI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Percentile:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%percentile,&
       "Percentile",rc)

  call ESMF_ConfigFindLabel(attribConfig,"River flow variate:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%rfv,&
       "RFV",rc)

  call ESMF_ConfigFindLabel(attribConfig,"MinTime:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%mintime,"MINTIME",rc)

  call ESMF_ConfigFindLabel(attribConfig,"MaxTime:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%maxtime,"MAXTIME",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Tendency:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%tendency,"TENDENCY",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Tendency correlation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%tendencycorr,&
       "TENDENCYCORR",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Z score:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%zscore,"ZSCORE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Trend:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%trend,"TREND",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Standard dS index:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%SdSI,"SdSI",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Triple collocation:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%tc,"TC",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Detection failure ratio (DFR):", &
       rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%dfr,"DFR",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Event frequency (EF):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ef,"EF",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Forecast frequency (FF):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ff,"FF",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Heidke skill score (HSS):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%hss,"HSS",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Peirce skill score (PSS):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%pss,"PSS",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Clayton skill score (CSS):",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%css,"CSS",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Reliability:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%rel,"REL",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Resilience:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%res,"RES",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Vulnerability:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%vul,"VUL",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Kmeans:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%kmeans,"KMEANS",rc)

  ! Tian bias decomposition...EMK
  call ESMF_ConfigFindLabel(attribConfig,"Tian Hit Bias:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%thb,"THB",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Tian Miss Bias:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%tmb,"TMB",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Tian False Alarm Bias:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%tfb,"TFB",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Metric entropy:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%mentropy,"Mentropy",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Information gain:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%igain,"Igain",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Fluctuation complexity:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%fcomplexity, &
       "Fcomplexity",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Effective complexity:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ecomplexity, &
       "Ecomplexity",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Information entropy:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ie, &
       "IE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Conditional entropy:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%ce, &
       "CE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Relative entropy:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%re, &
       "RE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Joint entropy:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%je, &
       "JE",rc)

  call ESMF_ConfigFindLabel(attribConfig,"Mutual information:",rc=rc)
  call get_metric_attributes(attribConfig, LVT_metrics%mi, &
       "MI",rc)


end subroutine LVT_readMetricsAttributes

!BOP
! 
! !ROUTINE: get_metric_attributes
! \label{get_metric_attributes}
!
! !INTERFACE:
subroutine get_metric_attributes(attribConfig,attribEntry,&
     short_name, status)
! 
! !USES: 
  use ESMF
  use LVT_statsDataMod
  use LVT_logMod, only : LVT_verify

  implicit none
!
! !INPUT PARAMETERS: 
  type(ESMF_Config),      intent(inout) :: attribConfig
  type(LVT_metricEntry), intent(inout) :: attribEntry
  character(len=*),        intent(in)   :: short_name
  integer,                 intent(in)   :: status
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine reads the run time specification for a specific
!  analysis metric. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  integer                               :: rc

  if(status.eq.0) then 
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%selectOpt,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading select option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%timeOpt,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading in-time option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%writeTS,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading write time series files option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%extractTS,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading extract time series for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%minthreshold,&
          default=0.0,rc=rc)
     call LVT_verify(rc,'Error reading min threshold option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%maxthreshold,&
          default=0.0,rc=rc)
     call LVT_verify(rc,'Error reading max threshold option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%computeSC,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading compute seasonal cycle option for '//trim(short_name))
     call ESMF_ConfigGetAttribute(attribConfig,attribEntry%computeADC,&
          default=0,rc=rc)
     call LVT_verify(rc,'Error reading compute average diurnal cycle option for '//trim(short_name))
     attribEntry%short_name  = (short_name)
  else

     attribEntry%selectOpt = 0
     attribEntry%timeOpt = 0 
     attribEntry%short_name = (short_name)

  endif
end subroutine get_metric_attributes
