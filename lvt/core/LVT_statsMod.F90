!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: LVT_statsMod
! \label{LVT_statsMod}
!
! !INTERFACE:
module LVT_statsMod
!
! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif
  use grib_api
  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_TSMod
  use LVT_CIMod
  use LVT_logMod
  use LVT_pluginIndices

  implicit none 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  The code in this file controls the flow of various statistics computations
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  02 Oct 2008: Sujay Kumar; Initial version
! 
!EOP
!BOP
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_statsInit
  public :: LVT_diagnoseStats
  public :: LVT_computeStats
  public :: LVT_checkTavgSpecs
  public :: LVT_writeSummaryStats
  public :: LVT_writeSummaryStats2
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------

!EOP

contains

!BOP
! 
! !ROUTINE: LVT_statsInit
! \label{LVT_statsInit}
!
! !INTERFACE:   
  subroutine LVT_statsInit
! 
! !USES: 
    use ESMF
    use LVT_timeMgrMod,         only : LVT_clock, LVT_calendar, LVT_seconds2time
    use LVT_InformationContentMod, only : LVT_initInformationContent
    use LVT_StratStatsMod,      only : LVT_initStratStats
    use LVT_CIMod,              only : LVT_initCI

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine initializes data structures required for various statistics
!  computations. It also initializes all the time series files, masking 
!  routines and variables, stratification routines and variables. Finally
!  it reads the restart files if the LVT run is kicked off of a previous
!  analysis
! 
!   The routines invoked are: 
!    \begin{description}
!    \item[LVT\_TSinit] (\ref{LVT_TSinit}) \newline
!      initializes the data structures required for time series
!      calculations and outputs
!    \item[LVT\_initInformationContent] (\ref{LVT_initInformationContent}) \newline
!      initializes the information theory based metrics calculations
!    \item[registerMetricEntry] (\ref{registerMetricEntry}) \newline
!      registers the specified metric 
!    \item[initMetricFiles] (\ref{initMetricFiles}) \newline
!      initializes the output file handles related to a metric
!    \item[LVT\_initCI] (\ref{LVT_initCI}) \newline
!      initializes the confidence interval calculations
!    \item[LVT\_initStratStats] (\ref{LVT_initStratStats}) \newline
!      initializes the data stratification routines and data structures
!    \item[LVT\_getDataStream1Ptr] (\ref{LVT_getDataStream1Ptr}) \newline
!      extracts the object at the head of the data stream 1 linked list
!    \item[LVT\_getDataStream2Ptr] (\ref{LVT_getDataStream2Ptr}) \newline
!      extracts the object at the head of the data stream 2 linked list
!    \item[LVT\_getS] (\ref{LVT_getstatsEntryPtr}) \newline
!      extracts the object at the head of the stats linked list
!    \item[] (\ref{initStatsEntry}) \newline
!      initializes the stats entry (stats are specified for each selected
!      variable pair)
!    \item[] (\ref{readmetricrestart}) \newline
!      reads the restart files for a specified metric. 
!    \end{description}       
!EOP

    integer                 :: ftn
    character*2             :: fadc
    integer                 :: da, hr, mn, ss
    real(ESMF_KIND_R8)      :: nts
    integer                 :: i,m,k,rc
    real                    :: nc,nr
    character*500           :: metricsAttribFile
    integer                 :: exp_v1, exp_v2
    character*50            :: anomalyTwindow
    type(ESMF_Time)         :: startTime, stopTime
    type(ESMF_TimeInterval) :: timestep
    integer                 :: selectLevels(LVT_rc%nDataStreams)

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_metadataEntry), pointer :: ds3
    type(LVT_statsEntry),    pointer :: stats


    call ESMF_ConfigGetAttribute(LVT_config,metricsAttribFile,&
         label="Metrics attributes file:",&
         rc=rc)
    call LVT_verify(rc,'Metrics attributes file: not defined')
!----------------------------------------------------------------------------
! read in the attributes of the metrics to be computed
!----------------------------------------------------------------------------
    
    call LVT_readMetricsAttributes(metricsAttribFile)

! interval for anomaly metrics calculations
    if((LVT_metrics%anomaly%selectOpt.gt.0).or.&
         (LVT_metrics%acorr%selectOpt.gt.0).or.&
         (LVT_metrics%armse%selectOpt.gt.0)) then 
       call ESMF_ConfigGetAttribute(LVT_config,&
            anomalyTwindow,&
            label="Averaging window for computing mean values in anomaly calculations:",rc=rc)
       if(rc.ne.0) then 
          write(LVT_logunit,*) &
               "[ERR] Averaging window for computing mean values in anomaly calculations: not defined"
          write(LVT_logunit,*) "[ERR] Supported options are 'monthly' or 'yearly'"
          call LVT_endrun()
       endif

       if(anomalyTwindow.eq."monthly") then 
          LVT_rc%anomalyTlength = 12
       else
          LVT_rc%anomalyTlength = 1
       endif
    endif
! variable based stratification, if enabled there are three levels 
! 1 - without any stratification
! 2 - above the stratification threshold, 
! 3 - below the stratification threshold, 
!
    LVT_rc%strat_nlevels = 1
    if(LVT_rc%var_based_strat .gt. 0) then 
       LVT_rc%strat_nlevels = 3 
    endif

    call LVT_seconds2time(LVT_rc%tavgInterval, da,hr,mn,ss)          

    LVT_rc%timeAvgOpt = 0        
    LVT_rc%prev_mo_rst = -1
    LVT_rc%monthCount = 0 
    LVT_rc%dayCount = 0 
! TSinit call must come before initMetricFiles because this routine
! sets the location information for ASCII time series files

    call LVT_TSinit()

!for seasonal cycle computations
    LVT_rc%nasc = 0 
    if(LVT_rc%scInterval.eq.1) then 
       LVT_rc%nasc = 12
       allocate(LVT_rc%scname(LVT_rc%nasc))
       LVT_rc%scname(1) = 'JAN'
       LVT_rc%scname(2) = 'FEB'
       LVT_rc%scname(3) = 'MAR'
       LVT_rc%scname(4) = 'APR'      
       LVT_rc%scname(5) = 'MAY'
       LVT_rc%scname(6) = 'JUN'
       LVT_rc%scname(7) = 'JUL'
       LVT_rc%scname(8) = 'AUG'      
       LVT_rc%scname(9) = 'SEP'
       LVT_rc%scname(10) = 'OCT'
       LVT_rc%scname(11) = 'NOV'
       LVT_rc%scname(12) = 'DEC'      
    elseif(LVT_rc%scInterval.eq.2) then 
       LVT_rc%nasc = 4 ! =12/3
       allocate(LVT_rc%scname(LVT_rc%nasc))
       LVT_rc%scname(1) = 'DJF'
       LVT_rc%scname(2) = 'MAM'
       LVT_rc%scname(3) = 'JJA'
       LVT_rc%scname(4) = 'SON'      
    elseif(LVT_rc%scInterval.eq.6) then 
       LVT_rc%nasc = 2 ! =12/6
       allocate(LVT_rc%scname(LVT_rc%nasc))
       LVT_rc%scname(1) = 'HF1'
       LVT_rc%scname(2) = 'HF2'
    elseif(LVT_rc%scInterval.eq.12) then 
       LVT_rc%nasc = 1
       allocate(LVT_rc%scname(LVT_rc%nasc))
       LVT_rc%scname(2) = 'YYR'
    endif
! The average diurnal cycle will be resolved at the stats
! output frequency
! 
    LVT_rc%nadc = nint(86400.0/LVT_rc%statswriteint)
    allocate(LVT_rc%adcname(LVT_rc%nadc))
!TODO: need to assign adcnames    
    do k=1,LVT_rc%nadc
       write(fadc,'(i2.2)') k
       LVT_rc%adcname(k) ='TINDEX_'//trim(fadc) 
    enddo

! Total number of LVT output times
    call ESMF_ClockGet(LVT_clock,runtimeStepCount=nts,&
         rc=rc)
    call LVT_verify(rc,&
         'Error in clockGet in LVT_statsinit')
!subtracting 1 because the ticktime advances the clock for
!the first analysis step. So the actual analysis period is
!from starting time + dt to ending time. 
    LVT_rc%nts = (nint(nts) + 1)-1

!Total number of LVT temporal averaging calculations
    
    call ESMF_ClockGet(LVT_clock, startTime = startTime, &
         stopTime = stopTime, rc=rc)
    call ESMF_TimeIntervalSet(timeStep, s = LVT_rc%tavgInterval, &
         rc=rc)
    LVT_rc%ntavgs = nint((stopTime-startTime)/timestep) + 1
    
    if(LVT_rc%computeICmetrics.eq.1) then 
       call LVT_initInformationContent
    endif
    

    if(LVT_rc%computeEnsMetrics.eq.1) then 
       LVT_rc%metric_sindex = LVT_ENSMETRIC_SINDEX
       LVT_rc%metric_eindex = LVT_ENSMETRIC_EINDEX
    elseif(LVT_rc%computeICmetrics.eq.1) then 
       call registerMetricEntry(LVT_mentropyid,LVT_metrics%mentropy)
       call registerMetricEntry(LVT_igainid,LVT_metrics%igain)
       call registerMetricEntry(LVT_fcomplexityid,LVT_metrics%fcomplexity)
       call registerMetricEntry(LVT_ecomplexityid,LVT_metrics%ecomplexity)
       LVT_rc%metric_sindex = LVT_ICMETRIC_SINDEX
       LVT_rc%metric_eindex = LVT_ICMETRIC_EINDEX     
    else
       call registerMetricEntry(LVT_MEANId,LVT_metrics%mean)
       call registerMetricEntry(LVT_MINId,LVT_metrics%min)
       call registerMetricEntry(LVT_MAXId,LVT_metrics%max)
       call registerMetricEntry(LVT_SUMId,LVT_metrics%sum)
       call registerMetricEntry(LVT_StdevId,LVT_metrics%stdev)
       call registerMetricEntry(LVT_VarianceId,LVT_metrics%variance)
       call registerMetricEntry(LVT_AnomalyId,LVT_metrics%anomaly)
       call registerMetricEntry(LVT_RMSEId,LVT_metrics%rmse)
       call registerMetricEntry(LVT_BIASId,LVT_metrics%bias)
       call registerMetricEntry(LVT_MAEId,LVT_metrics%mae)
       call registerMetricEntry(LVT_PODYId,LVT_metrics%pody)
       call registerMetricEntry(LVT_PODNId,LVT_metrics%podn)
       call registerMetricEntry(LVT_FARId,LVT_metrics%far)
       call registerMetricEntry(LVT_POFDId,LVT_metrics%POFD)
       call registerMetricEntry(LVT_CSIId,LVT_metrics%CSI)
       call registerMetricEntry(LVT_ACCId,LVT_metrics%ACC)
       call registerMetricEntry(LVT_FBIASId,LVT_metrics%FBIAS)
       call registerMetricEntry(LVT_ETSId,LVT_metrics%ETS)
       call registerMetricEntry(LVT_RCORRId,LVT_metrics%RCORR)
       call registerMetricEntry(LVT_RNKCORRId,LVT_metrics%RNKCORR)
       call registerMetricEntry(LVT_ACORRId,LVT_metrics%ACORR)
       call registerMetricEntry(LVT_ARMSEId,LVT_metrics%ARMSE)
       call registerMetricEntry(LVT_NSEId,LVT_metrics%NSE)
       call registerMetricEntry(LVT_ubRMSEId,LVT_metrics%ubRMSE)       
       call registerMetricEntry(LVT_AREAId,LVT_metrics%AREA)
       call registerMetricEntry(LVT_waveletStatId,LVT_metrics%waveletStat)
       call registerMetricEntry(LVT_hnId, LVT_metrics%hn)
       call registerMetricEntry(LVT_spiId, LVT_metrics%spi)
       call registerMetricEntry(LVT_sriId, LVT_metrics%sri)
       call registerMetricEntry(LVT_sswiId, LVT_metrics%sswi)
       call registerMetricEntry(LVT_sgwiId, LVT_metrics%sgwi)
       call registerMetricEntry(LVT_percentileId, LVT_metrics%percentile)
       call registerMetricEntry(LVT_rfvId, LVT_metrics%rfv)
       call registerMetricEntry(LVT_MINTIMEId,LVT_metrics%mintime)
       call registerMetricEntry(LVT_MAXTIMEId,LVT_metrics%maxtime)
       call registerMetricEntry(LVT_TendencyId,LVT_metrics%tendency)
       call registerMetricEntry(LVT_TendencyCorrId,LVT_metrics%tendencyCorr)
       call registerMetricEntry(LVT_ZscoreId,LVT_metrics%zscore)
       call registerMetricEntry(LVT_TrendId,LVT_metrics%trend)
       call registerMetricEntry(LVT_SdSIId, LVT_metrics%SdSI)
       call registerMetricEntry(LVT_TCId,LVT_metrics%tc)
       call registerMetricEntry(LVT_DFRId,LVT_metrics%dfr) ! EMK
       call registerMetricEntry(LVT_EFId,LVT_metrics%ef) ! EMK
       call registerMetricEntry(LVT_FFId,LVT_metrics%ff) ! EMK
       call registerMetricEntry(LVT_HSSId,LVT_metrics%hss) ! EMK
       call registerMetricEntry(LVT_PSSId,LVT_metrics%pss) ! EMK
       call registerMetricEntry(LVT_CSSId,LVT_metrics%css) ! EMK
       call registerMetricEntry(LVT_RELId,LVT_metrics%rel) 
       call registerMetricEntry(LVT_RESId,LVT_metrics%res) 
       call registerMetricEntry(LVT_VULId,LVT_metrics%vul) 
       call registerMetricEntry(LVT_KMEANSId,LVT_metrics%kmeans)
       ! Tian bias decomposition...EMK
       call registerMetricEntry(LVT_THBId,LVT_metrics%thb)
       call registerMetricEntry(LVT_TMBId,LVT_metrics%tmb)
       call registerMetricEntry(LVT_TFBId,LVT_metrics%tfb)

       LVT_rc%metric_sindex = LVT_METRIC_SINDEX
       LVT_rc%metric_eindex = LVT_METRIC_EINDEX
    endif   

    if(LVT_metrics%waveletStat%selectOpt.eq.1) then 
       
       !padding to the nearest 2 power dimension. 
       
       nc = 2**ceiling(log(float(LVT_rc%lnc))/log(2.0))
       nr = 2**ceiling(log(float(LVT_rc%lnr))/log(2.0))
       
       nc = max(nc,nr)
       nr = max(nc,nr)
       
       exp_v1 = nint(log(nc)/log(2.0))
       exp_v2 = nint(log(nr)/log(2.0))      
       
       LVT_rc%nscales = exp_v1
       
       allocate(LVT_rc%lnc_sc(LVT_rc%nscales))
       allocate(LVT_rc%lnr_sc(LVT_rc%nscales))
       
       do i=1,LVT_rc%nscales 
          LVT_rc%lnc_sc(i) = 2**exp_v1
          LVT_rc%lnr_sc(i) = 2**exp_v2
          
          exp_v1 = exp_v1 - 1
          exp_v2 = exp_v2 - 1
          
       enddo
    endif
       
    do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
       call initMetricFiles(LVT_metricsPtr(m)%metricEntryPtr) 
    enddo
    
    call LVT_initCI()
    call LVT_initStratStats()

!    if(LVT_rc%runmode.ne.LVT_dastatId) then 
    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    if(LVT_rc%nDataStreams.gt.2) then 
       call LVT_getDataStream3Ptr(ds3)
    endif

    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model)) 
       selectLevels(1) = model%selectNlevs
       selectLevels(2) = obs%selectNlevs
       if(LVT_rc%nDataStreams.gt.2) then 
          selectLevels(3) = ds3%selectNlevs
       endif
    
       call initStatsEntry(stats, selectLevels)
       
       model => model%next
       obs   => obs%next
       stats => stats%next
       
    enddo
!---------------------------------------------------------------------
! determine if time series files are required (both ascii and 
! domain files) and the number of passes through the data
!---------------------------------------------------------------------
    LVT_rc%wtsout = 0
    LVT_rc%pass = 0
    LVT_rc%extractTS = 0 
    LVT_rc%computeErrSC = 0
    do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
       if(LVT_metricsPtr(m)%metricEntryPtr%timeOpt.eq.1) then 
          LVT_rc%wtsout = 1
       endif
       if(LVT_metricsPtr(m)%metricEntryPtr%timeOpt.eq.1.and.&
            LVT_metricsPtr(m)%metricEntryPtr%extractTS.eq.1) then 
          LVT_rc%extractTS = 1
       endif
       if(LVT_metricsPtr(m)%metricEntryPtr%timeOpt.eq.1.and.&
            LVT_metricsPtr(m)%metricEntryPtr%computeSC.eq.1) then 
          LVT_rc%computeErrSC = 1
          if(LVT_rc%nasc.eq.0) then 
             write(LVT_logunit,*) & 
                  '[ERR] The seasonal cycle interval must be set'
             write(LVT_logunit,*) &
                  '[ERR] when seasonal cycle is computed'
             call LVT_endrun()
          endif
       endif
       
       if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.eq.1.or.&
            LVT_metricsPtr(m)%metricEntryPtr%timeOpt.eq.1) then 
          LVT_rc%pass = max(LVT_rc%pass,&
               LVT_metricsPtr(m)%metricEntryPtr%npass)
       endif
    enddo
    if(LVT_rc%pass.eq.0) then 
       write(LVT_logunit,*) &
            '[ERR] LVT determines the number of passes through the data to be zero'
       write(LVT_logunit,*) &
            '[ERR] No analysis will be conducted. Please check the settings in the metric attributes file'
       call LVT_endrun()
    endif

    if(LVT_rc%statswriteint.ge.86400.and.LVT_rc%computeADC.eq.1) then  
       write(LVT_logunit,*) '[ERR] Please set the stats output interval to less than 86400 (day)' 
       write(LVT_logunit,*) '[ERR] when computing average diurnal cycle of error metrics '
       write(LVT_logunit,*) '[ERR] option is enabled'
       call LVT_endrun()
    endif
    if((LVT_rc%startmode).eq."restart") then 

       if(.not. LVT_metrics%percentile%selectOpt.eq.1) then 
!for percentile calculations, they are handled differently
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=(LVT_rc%rstfile),form='unformatted')
          write(LVT_logunit,*) '[INFO] Reading restart file ',(LVT_rc%rstfile)
          
!          read(ftn) LVT_rc%curr_pass
          do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
             if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
                call readmetricrestart(m,ftn)
             endif
          enddo
          call LVT_releaseUnitNumber(ftn)
       endif
    endif

    if(LVT_rc%var_strat_index.gt.0) then 
       model => LVT_histData%ptr_into_ds1_list(&
            LVT_rc%var_strat_index)%dataEntryPtr
       allocate(LVT_stats%strat_var(LVT_rc%ngrid,&
            LVT_rc%nensem,model%vlevels))
    endif
    
!special checks to determine the number of passes for derived variables
    call LVT_checkForDerivedVariableDependendency()

  end subroutine LVT_statsInit

!BOP
! 
! !ROUTINE: initMetricFiles
! \label{initMetricFiles}
!
! !INTERFACE: 
  subroutine initMetricFiles(metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine creates the filenames (METADATA, Summary stats
!  and the time series) associated with each analysis metric
! 
!  In a restart mode, the timeseries files already written are copied
!  and parsed to the start of the restart run. The subsequent 
!  analysis is then appended to the existing time series files. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:     
    type(LVT_metricEntry) :: metric
!EOP  
    integer                 :: ftn1,ftn2
    integer                 :: i,m
    integer                 :: ios
    type(ESMF_Time)         :: ftime,currTime
    character*4             :: fens
    character*500           :: filename
    character*500           :: meta_output_file
    character*500           :: summ_output_file
    integer                 :: yr,mo,da,hr,mn
    character*5000          :: c_line
    integer                 :: status

    if(metric%timeopt.eq.1.and.metric%extractTS.eq.1) then 

       if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
          allocate(metric%ftn_ts_loc(LVT_rc%ntslocs,LVT_rc%nensem))

          call system("mkdir -p "//trim(LVT_rc%statsodir))       
          
          do m=1,LVT_rc%nensem
             do i=1,LVT_rc%ntslocs
                
                if(LVT_rc%nensem.gt.1) then 
                   write(fens,fmt='(i4.4)') m
                   filename = trim(LVT_rc%statsodir)//'/'//&
                        trim(metric%short_name)//'_'//&
                        trim(LVT_TSobj(i)%tslocname)//&
                        '_'//trim(fens)//&
                        '.dat'
                else
                   filename = trim(LVT_rc%statsodir)//'/'//&
                        trim(metric%short_name)//'_'//&
                        trim(LVT_TSobj(i)%tslocname)//&
                        '.dat'
                endif
                metric%ftn_ts_loc(i,m) = LVT_getNextUnitNumber()
                
                if((LVT_rc%startmode).eq."coldstart") then 
                   open(metric%ftn_ts_loc(i,m),file=(filename),&
                        form='formatted')
                elseif((LVT_rc%startmode).eq."restart") then
                   !copy the files over before
                   !             call system('cp '//(filename)//' temp')
                   ftn1 = LVT_getNextUnitNumber()
                   ftn2 = LVT_getNextUnitNumber()
                   open(ftn1,file='temp',form='formatted')
                   open(ftn2,file=(filename),form='formatted')
                   ios = 0 
                   do while(ios.eq.0) 
                      read(ftn2,'(a)',iostat=ios) c_line
                      !get the time information:
                      if(ios.ne.0) exit
                      read(c_line(1:4),*) yr
                      read(c_line(6:7),*) mo
                      read(c_line(9:10),*) da
                      read(c_line(12:13),*) hr
                      read(c_line(15:16),*) mn
                      
                      call ESMF_TimeSet(fTime,  yy=yr, &
                           mm = mo, &
                           dd = da, &
                           h = hr, &
                           m = mn, &
                           calendar = LVT_calendar, &
                           rc=status)
                      call LVT_verify(status, 'error in initMetricFiles')    
                      call ESMF_TimeSet(currTime,  yy=LVT_rc%syr, &
                           mm = LVT_rc%smo, &
                           dd = LVT_rc%sda, &
                           h = LVT_rc%shr, &
                           m = LVT_rc%smn, &
                           calendar = LVT_calendar, &
                           rc=status)
                      call LVT_verify(status, 'error in initMetricFiles')    
                      if(ftime.gt.currTime) then 
                         ios = -1
                         exit
                      endif
                      write(ftn1,'(a)') trim(c_line)
                   enddo
                   call LVT_releaseUnitNumber(ftn1)
                   call LVT_releaseUnitNumber(ftn2)
!after parsing, move the temp file back. 
                   call system('mv temp '//trim(filename))

                   open(metric%ftn_ts_loc(i,m),file=(filename),&
                        ACCESS = 'APPEND',form='formatted')
                endif
             enddo
          end do
       else
          allocate(metric%ftn_ts_loc(LVT_rc%ntslocs,1))

          call system("mkdir -p "//trim(LVT_rc%statsodir))       
          
          do i=1,LVT_rc%ntslocs            
             filename = trim(LVT_rc%statsodir)//'/'//&
                  trim(metric%short_name)//'_'//&
                  trim(LVT_TSobj(i)%tslocname)//&
                  '.dat'
             metric%ftn_ts_loc(i,1) = LVT_getNextUnitNumber()
             
             if((LVT_rc%startmode).eq."coldstart") then 
                open(metric%ftn_ts_loc(i,1),file=(filename),&
                     form='formatted')
             elseif((LVT_rc%startmode).eq."restart") then
                !copy the files over before
                   !             call system('cp '//(filename)//' temp')
                ftn1 = LVT_getNextUnitNumber()
                ftn2 = LVT_getNextUnitNumber()
                open(ftn1,file='temp',form='formatted')
                open(ftn2,file=(filename),form='formatted')
                ios = 0 
                do while(ios.eq.0) 
                   read(ftn2,'(a)',iostat=ios) c_line
                   !get the time information:
                   if(ios.ne.0) exit
                   read(c_line(1:4),*) yr
                   read(c_line(6:7),*) mo
                   read(c_line(9:10),*) da
                   read(c_line(12:13),*) hr
                   read(c_line(15:16),*) mn
                   
                   call ESMF_TimeSet(fTime,  yy=yr, &
                        mm = mo, &
                        dd = da, &
                        h = hr, &
                        m = mn, &
                        calendar = LVT_calendar, &
                        rc=status)
                   call LVT_verify(status, 'error in initMetricFiles')    
                   call ESMF_TimeSet(currTime,  yy=LVT_rc%syr, &
                        mm = LVT_rc%smo, &
                        dd = LVT_rc%sda, &
                        h = LVT_rc%shr, &
                        m = LVT_rc%smn, &
                        calendar = LVT_calendar, &
                        rc=status)
                   call LVT_verify(status, 'error in initMetricFiles')    
                   if(ftime.gt.currTime) then 
                      ios = -1
                      exit
                   endif
                   write(ftn1,'(a)') trim(c_line)
                enddo
                call LVT_releaseUnitNumber(ftn1)
                call LVT_releaseUnitNumber(ftn2)
!after parsing, move the temp file back. 
                call system('mv temp '//trim(filename))

                open(metric%ftn_ts_loc(i,1),file=(filename),&
                     ACCESS = 'APPEND',form='formatted')
             end if
          enddo
       end if
    endif
    if((LVT_rc%lvt_out_format).eq."binary".and.metric%selectOpt.eq.1) then 
       call system("mkdir -p "//trim(LVT_rc%statsodir))
       meta_output_file = trim(LVT_rc%statsodir)//'/'//&
            trim(metric%short_name)//'_'//& 
            'METADATA.dat'
       metric%ftn_meta_out = LVT_getNextUnitNumber()
       open(metric%ftn_meta_out,file=(meta_output_file))

    endif

    if(metric%selectOpt.eq.1) then 
       call system("mkdir -p "//trim(LVT_rc%statsodir))
       summ_output_file = trim(LVT_rc%statsodir)//'/'//&
            trim(metric%short_name)//'_'//& 
            'SUMMARY_STATS.dat'
       metric%ftn_summ = LVT_getNextUnitNumber()
       open(metric%ftn_summ,file=(summ_output_file))

    endif

  end subroutine initMetricFiles


!BOP
! 
! !ROUTINE: initStatsEntry
! \label{initStatsEntry}
!
! !INTERFACE: 
  subroutine initStatsEntry(stats, selectNlevs)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine initializes the objects to hold different statistics 
!  computations
!
!   The arguments are: 
!   \begin{description}
!    \item[model] object to hold model variable information
!    \item[stats] object to hold statistics computations
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:     
    type(LVT_statsEntry) :: stats
    integer              :: selectNlevs(LVT_rc%nDataStreams)
!EOP
    integer            :: m
    
    if(stats%selectOpt.eq.1.and.selectNlevs(1).ge.1) then 
! Hardcoded the second index to 4 to account for multiple fields within each
! metric (Assuming that we don't need more than 2 fields within a metric.
! each field will need to save space for 'model' and 'obs'
!
       allocate(stats%vid_total(LVT_NMETRICS,20))
       allocate(stats%vid_count_total(LVT_NMETRICS,20))
       allocate(stats%vid_stdev_total(LVT_NMETRICS,20))
       allocate(stats%vid_count_stdev_total(LVT_NMETRICS,20))
       allocate(stats%vid_ts(LVT_NMETRICS,20))
       allocate(stats%vid_count_ts(LVT_NMETRICS,20))
       allocate(stats%vid_sc_total(LVT_rc%nasc,LVT_NMETRICS,20))
       allocate(stats%vid_adc_total(LVT_rc%nadc,LVT_NMETRICS,20))

       do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
          if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
             
             LVT_metricsPtr(m)%metricEntryPtr%nLevs = 1
             LVT_metricsPtr(m)%metricEntryPtr%customNames = .false.

             call initmetric(m,selectNlevs,stats,&
                  LVT_metricsPtr(m)%metricEntryPtr)
          endif
       enddo
    endif
    
  end subroutine initStatsEntry

!BOP
! 
! !ROUTINE: registerMetricEntry
! \label{registerMetricEntry}
!
! !INTERFACE: 
  subroutine registerMetricEntry(metric_index, metric)
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
! 
!EOP
!BOP
! 
    
    integer                  :: metric_index
    type(LVT_metricEntry), target :: metric

    LVT_MetricsPtr(metric_index)%metricEntryPtr => metric

  end subroutine registerMetricEntry

  subroutine LVT_checkTavgSpecs()

    type(LVT_metaDataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats
    
    logical :: tavg_check_status
    logical :: inst_check_status
    logical :: accum_check_status

    call LVT_getDataStream1Ptr(model)

    tavg_check_status = .false. 
    inst_check_status = .false.
    accum_check_status = .false. 

    do while(associated(model))

       ! EMK...Skip checking Tair_f_min
       if (trim(model%short_name) .eq. "RHMin") then
          model => model%next
       end if

       if(model%selectNlevs.ge.1.and. model%timeAvgOpt.eq.0) then 
          inst_check_status = .true. 
       elseif(model%selectNlevs.ge.1.and. model%timeAvgOpt.eq.1) then 
          tavg_check_status = .true. 
       elseif(model%selectNlevs.ge.1.and. model%timeAvgOpt.eq.3) then 
          accum_check_status = .true. 
       endif

       model => model%next
    enddo
    
    if(inst_check_status.and.tavg_check_status) then 
       write(LVT_logunit,*) '[ERR] LVT does not support the simultaneous use of '
       write(LVT_logunit,*) '[ERR] both time averaged and instantaneous variables '
       write(LVT_logunit,*) '[ERR] in analysis. Please create separate instances '
       write(LVT_logunit,*) '[ERR] for analyzing such variables.'
       call LVT_endrun()
    endif
    if(inst_check_status.and.accum_check_status) then 
       write(LVT_logunit,*) '[ERR] LVT does not support the simultaneous use of '
       write(LVT_logunit,*) '[ERR] both instantaneous and accumulated variables '
       write(LVT_logunit,*) '[ERR] in analysis. Please create separate instances '
       write(LVT_logunit,*) '[ERR] for analyzing such variables.'
       call LVT_endrun()
    endif
    LVT_rc%timeAvgOpt = -1
    if(inst_check_status) then 
       LVT_rc%timeAvgOpt = 0
    endif
    !EMK...Add separate entry for accumulations
!    if(tavg_check_status.or.accum_check_status) then 
!       LVT_rc%timeAvgOpt = 1
!    endif    
    if(tavg_check_status) then 
       LVT_rc%timeAvgOpt = 1
    endif    
    if(accum_check_status) then 
       LVT_rc%timeAvgOpt = 3
    endif    

  end subroutine LVT_checkTavgSpecs
    

!BOP
! 
! !ROUTINE: LVT_diagnoseStats
! \label{LVT_diagnoseStats}
!
! !INTERFACE:
  subroutine LVT_diagnoseStats(pass)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This routine invokes the methods for computing the specified statistics
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer       :: m
    integer       :: pass

    if(LVT_rc%computeFlag) then     
       do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
          if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
             call diagnosemetric(m,pass)
          endif
       enddo
    endif
    
  end subroutine LVT_diagnoseStats

!BOP
! 
! !ROUTINE: LVT_computeStats
! \label{LVT_computeStats}
!
! !INTERFACE: 
  subroutine LVT_computeStats(pass)
! 
! !USES:   
    use LVT_timeMgrMod,      only : LVT_calendar
    use LVT_DataStreamsMod  

    implicit none
! !ARGUMENTS: 
    integer :: pass
!EOP
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine issues the calls to compute the specified set of statistics 
!  There are two different alarms at play here. Based on the time averaging
!  interval, the 'computeFlag' variable (logical) is set. If the value of
!  'computeFlag' is true, then it triggers the invocation of the compute
!  part of each metric that is enabled (For example, the RMSE is calculated
!  at this point, whereas in all previous timesteps, the values are simply
!  logged).
! 
!  The second alarm has to do with the stats writing interval, which 
!  could be greater than or equal to the time averaging interval. When this
!  is invoked, the average of all metric values (calculated at the
!  time averaging intervals) is written out. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    logical                 :: alarmCheck,alarmCheck_rst
    integer                 :: m,nfrac
    integer                 :: mfactor
    integer                 :: yr, mo, da, hr, mn, ss
    integer                 :: nyr, nmo, nda, nhr, nmn, nss
    type(ESMF_Time)         :: currTime
    type(ESMF_TimeInterval) :: ts
    integer                 :: ftn
    character*500           :: rstfile
    character(len=12)       :: cdate
    character(len=4)        :: cdate1
    character*500           :: dir_string
    integer                 :: status
    logical                 :: alarmCheck_total

    yr = LVT_rc%yr
    mo = LVT_rc%mo
    da = LVT_rc%da
    hr = LVT_rc%hr
    mn = LVT_rc%mn
    ss = LVT_rc%ss   

    alarmCheck = .false. 

    if(LVT_rc%tsconv.eq."dekad") then 
       alarmCheck = .true.
    else          
       if(LVT_rc%timeAvgOpt.eq.0) then !instantaneous variables
          if(mod(LVT_rc%statswriteint,31536000).eq.0) then !yearly alarm
             if(LVT_rc%nmo.eq.LVT_rc%use_shift_mo) then 
                if(LVT_rc%nyr.ne.LVT_rc%prev_yr_sout) then 
                   LVT_rc%prev_yr_sout = LVT_rc%nyr
                   alarmCheck = .true. 
                endif
             endif
          elseif(mod(LVT_rc%statswriteint,15552000).eq.0)then !6 monthly 
             if(LVT_rc%nmo.ne.LVT_rc%prev_mo_sout) then 
                LVT_rc%prev_mo_sout = LVT_rc%nmo
                LVT_rc%monthCount_sout = LVT_rc%monthCount_sout + 1
                if(LVT_rc%monthCount_sout.eq.6) then 
                   alarmCheck = .true. 
                   LVT_rc%monthCount_sout = 0 
                endif
             endif
          elseif(mod(LVT_rc%statswriteint,7776000).eq.0)then !3 monthly 
             if(LVT_rc%nmo.ne.LVT_rc%prev_mo_sout) then 
                LVT_rc%prev_mo_sout = LVT_rc%nmo
                LVT_rc%monthCount_sout = LVT_rc%monthCount_sout + 1
                if(LVT_rc%monthCount_sout.eq.3) then 
                   alarmCheck = .true. 
                   LVT_rc%monthCount_sout = 0 
                endif
             endif
             
          elseif(mod(LVT_rc%statswriteint,2592000).eq.0) then 
             if(LVT_rc%nmo.ne.LVT_rc%prev_mo_sout) then 
                LVT_rc%prev_mo_sout = LVT_rc%nmo
                alarmCheck = .true. 
             endif
          elseif(mod(LVT_rc%statswriteint,604800).eq.0) then 
             if(mod(real(LVT_rc%nhr)*3600+60*real(LVT_rc%nmn) + & 
                  float(LVT_rc%nss),&
                  real(LVT_rc%statswriteint)).eq.0) then 
                LVT_rc%dayCount_sout = LVT_rc%dayCount_sout + 1
                if(LVT_rc%dayCount_sout.eq.7) then 
                   alarmCheck = .true. 
                   LVT_rc%dayCount_sout = 0
                endif
             endif
          elseif(LVT_rc%statswriteint.le.86400) then 
             if(mod(real(LVT_rc%nhr)*3600+60*real(LVT_rc%nmn)+&
                  float(LVT_rc%nss),&
                  real(LVT_rc%statswriteint)).eq.0) then        
                alarmCheck = .true.
             endif
          else
             write(LVT_logunit,*) '[ERR] The support for the specified stats output'
             write(LVT_logunit,*) '[ERR] interval is not supported currently. '  
             write(LVT_logunit,*) '[ERR] Please contact the LVT development team.'
             call LVT_endrun()
          endif
       elseif(LVT_rc%timeAvgOpt.eq.1.or.LVT_rc%timeAvgOpt.eq.3) then !time averaged variables
          if(mod(LVT_rc%statswriteint,31536000).eq.0) then !yearly alarm
             if(LVT_rc%mo.eq.LVT_rc%use_shift_mo) then 
                if(LVT_rc%yr.ne.LVT_rc%prev_yr_sout) then 
                   LVT_rc%prev_yr_sout = LVT_rc%yr
                   alarmCheck = .true. 
                endif
             endif
          elseif(mod(LVT_rc%statswriteint,15552000).eq.0)then !6 monthly 
             if(LVT_rc%mo.ne.LVT_rc%prev_mo_sout) then 
                LVT_rc%prev_mo_sout = LVT_rc%mo
                LVT_rc%monthCount_sout = LVT_rc%monthCount_sout + 1
                if(LVT_rc%monthCount_sout.eq.6) then 
                   alarmCheck = .true. 
                   LVT_rc%monthCount_sout = 0 
                endif
             endif
          elseif(mod(LVT_rc%statswriteint,7776000).eq.0)then !3 monthly 
             if(LVT_rc%mo.ne.LVT_rc%prev_mo_sout) then 
                LVT_rc%prev_mo_sout = LVT_rc%mo
                LVT_rc%monthCount_sout = LVT_rc%monthCount_sout + 1
                if(LVT_rc%monthCount_sout.eq.3) then 
                   alarmCheck = .true. 
                   LVT_rc%monthCount_sout = 0 
                endif
             endif
             
          elseif(mod(LVT_rc%statswriteint,2592000).eq.0) then 
             if(LVT_rc%mo.ne.LVT_rc%prev_mo_sout) then 
                LVT_rc%prev_mo_sout = LVT_rc%mo
                alarmCheck = .true. 
             endif
          elseif(mod(LVT_rc%statswriteint,604800).eq.0) then 
             if(mod(real(LVT_rc%hr)*3600+60*real(LVT_rc%mn) + & 
                  float(LVT_rc%ss),&
                  real(LVT_rc%statswriteint)).eq.0) then 
                LVT_rc%dayCount_sout = LVT_rc%dayCount_sout + 1
                if(LVT_rc%dayCount_sout.eq.7) then 
                   alarmCheck = .true. 
                   LVT_rc%dayCount_sout = 0
                endif
             endif
          elseif(LVT_rc%statswriteint.le.86400) then 
             if(mod(real(LVT_rc%hr)*3600+60*real(LVT_rc%mn)+&
                  float(LVT_rc%ss),&
                  real(LVT_rc%statswriteint)).eq.0) then        
                alarmCheck = .true.
             endif
          else
             write(LVT_logunit,*) '[ERR] The support for the specified stats output'
             write(LVT_logunit,*) '[ERR] interval is not supported currently. '  
             write(LVT_logunit,*) '[ERR] Please contact the LVT development team.'
             call LVT_endrun()
          endif
       endif
          
    endif

    if(LVT_rc%computeFlag) then 

       if(alarmCheck.and.LVT_rc%wtsout.eq.1) then 
          call system("mkdir -p "//trim(LVT_rc%statsodir))
         
          do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
             if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
                call createTSfiles(pass, LVT_metricsPtr(m)%metricEntryPtr)
             endif
          enddo
       endif

!---------------------------------------------------------------------------
! If the endtime is reached, then the metrics are computed. 
! The alarmCheck_total is passed to the metric calculations, which 
! then figures out the temporal averaging of the metrics within 
! a stats writing interval
!---------------------------------------------------------------------------
       alarmCheck_total = alarmCheck.or.LVT_rc%endtime.eq.1
       
       do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
          if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
             call computemetric(m,pass,alarmCheck_total)
          endif
       enddo

       
       if((alarmCheck).and.LVT_rc%wtsout.eq.1) then 
          call outputTimeSeriesStats(pass)

          do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
             if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
                call finalizeTSfile(pass, LVT_metricsPtr(m)%metricEntryPtr)
             endif
          enddo
       endif

       do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
          if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
             call resetmetric(m, alarmCheck_total)
          endif
       enddo
       
       if(LVT_rc%endtime.eq.1) then              
          do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
             if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
                call createOutputFile(LVT_metricsPtr(m)%metricEntryPtr,pass)
             endif
          enddo
          
          call outputFinalStats(pass)
          
          do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
             if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
                call finalizeOutputFile(LVT_metricsPtr(m)%metricEntryPtr,pass)
             endif
          enddo
       endif
    
       LVT_stats%computeFlag = .false. 

!       call LVT_resetDataStreams
       
    endif

    alarmCheck_rst = .false.
    if(LVT_rc%wrst.ne.0) then 
       if(mod(LVT_rc%statswriteint,31536000).eq.0) then !yearly alarm
          if(yr.ne.LVT_rc%nyr) then 
             alarmCheck_rst = .true. 
          else
             alarmCheck_rst = .false. 
          endif
       elseif(mod(LVT_rc%restartInterval,2592000).eq.0) then !monthly alarm
          if(mo.ne.LVT_rc%nmo) then 
             alarmCheck_rst = .true. 
          else
             alarmCheck_rst = .false. 
          endif
       elseif(mod(real(hr)*3600+60*&
            real(mn)+real(ss),&
            real(LVT_rc%restartInterval)).eq.0.0) then 
          alarmcheck_rst = .true. 
       endif
       if(alarmCheck_rst) then 
          dir_string ="mkdir -p "//trim(LVT_rc%statsodir)//'/RST' 
          call system(dir_string)
          
          write(unit=cdate,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2)') &
               yr, mo, da, hr, mn
          
          rstfile = trim(LVT_rc%statsodir)//'/RST/LVT.'//cdate//'.rst'
          
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=(rstfile),form='unformatted')
          write(LVT_logunit,*) "[INFO] Writing restart file ",trim(rstfile)
          
          !       write(ftn) pass
          do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
             if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.gt.0) then
                call writemetricrestart(m,ftn,pass)
             endif
          enddo
          
          call LVT_releaseUnitNumber(ftn)
       else
          alarmcheck_rst = .false. 
       endif
    endif
  end subroutine LVT_computeStats

!BOP
! 
! !ROUTINE: LVT_writeSummaryStats
! \label{LVT_writeSummaryStats}
!
! !INTERFACE:
  subroutine LVT_writeSummaryStats(ftn, l, metricname, nsize, metric, npts, &
       varname, ci)
! 
! !USES:   
   implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine outputs a summary of various statistics computed 
!  during the LVT analysis. The statistics are computed for the whole
!  analysis domain and for each time series location specified. 
!
!  The arguments are: 
!  \begin{description}
!  \item[ftn]   unit number of the file
!  \item[metricname]  name of the metric 
!  \item[metric]    array containing values of the computed metric
!  \item[npts] array containing the number of counts used in computing 
!              the metric
!  \item[varname]  name of the variable being written out
!  \item[ci]  confidence interval associated with the computed metric values
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !LVT_writeSummaryStats
! !ARGUMENTS:   
   integer             :: ftn
   integer             :: l
   character(len=*)    :: metricname
   integer             :: nsize
   real                :: metric(nsize, LVT_rc%nensem)
   integer             :: npts(nsize, LVT_rc%nensem)
   character(len=*)    :: varname
   real                :: ci(LVT_rc%nensem)
!EOP  
   

   integer                 :: i,c,r,t,tid,m
   real                    :: sum_v
   integer                 :: nsum_v
   real                    :: wsum_v
   integer                 :: wnsum_v
   real, allocatable       :: sum_sd_v(:)
   integer, allocatable    :: nsum_sd_v(:)
   real, allocatable       :: wsum_sd_v(:)
   integer, allocatable    :: wnsum_sd_v(:)
!   real,    allocatable    :: ts_v(LVT_LIS_rc(1)%ntiles)
   real,    allocatable    :: ts_v(:)
   integer                 :: nts_v
   real,   allocatable     :: ci_val(:)
   type(ESMF_Time)         :: startTime,stopTime
   type(ESMF_TimeInterval) :: timeStep
   integer                 :: nTimeSpan
   integer                 :: status
   !percentage of data coverage
   real                    :: pct_coverage(nsize, LVT_rc%nensem)
   real                    :: pct_cov_avg
   integer                 :: npct_cov_avg
   real, allocatable       :: pct_cov_ts(:)
   integer, allocatable    :: npct_cov_ts(:)

   !total number of timesteps in the LVT analysis time period
   call ESMF_ClockGet(LVT_clock, startTime = startTime, &
        stopTime = stopTime, rc=status)
   call LVT_verify(status, 'ESMF_ClockGet failed in LVT_writeSummaryStats')
   
   call ESMF_TimeIntervalSet(timeStep, s = LVT_rc%tavgInterval, &
        rc=status)
   call LVT_verify(status, &
        'ESMF_TimeIntervalSet failed in LVT_writeSummaryStats')
   
   nTimeSpan = nint((stopTime - startTime)/timestep) + 1
   do t=1,nsize
      do m=1,LVT_rc%nensem
         pct_coverage(t,m) = (float(npts(t,m))/float(nTimeSpan))*100.0
      enddo
   enddo

   allocate(sum_sd_v(LVT_rc%ntslocs))
   allocate(nsum_sd_v(LVT_rc%ntslocs))
   allocate(ci_val(LVT_rc%ntslocs))
   
   allocate(pct_cov_ts(LVT_rc%ntslocs))
   allocate(npct_cov_ts(LVT_rc%ntslocs))
   
   pct_cov_ts = 0
   npct_cov_ts = 0 

   if(LVT_rc%computeEnsMetrics.eq.1.and.&
        nsize.eq.LVT_LIS_rc(1)%ntiles) then 
      allocate(ts_v(LVT_LIS_rc(1)%ntiles))
   else
      allocate(ts_v(nsize*LVT_rc%nensem))
   endif
   
   sum_v = 0 
   nsum_v = 0
   sum_sd_v = 0 
   nsum_sd_v = 0 
   ci_val = 0 
   pct_cov_avg = 0 
   npct_cov_avg = 0 

   do t=1,nsize
      do m=1,LVT_rc%nensem
         if(metric(t,m).ne.LVT_rc%udef) then 
            sum_v = sum_v + metric(t,m)
            nsum_v = nsum_v + 1
            pct_cov_avg = pct_cov_avg + pct_coverage(t,m)
            npct_cov_avg = npct_cov_avg + 1
         endif
      enddo
   enddo
   
   if(nsum_v.gt.0) then 
      sum_v = sum_v/nsum_v
   else
      sum_v = LVT_rc%udef
   endif
   if(npct_cov_avg.gt.0) then 
      pct_cov_avg = pct_cov_avg/npct_cov_avg
   else
      pct_cov_avg = 0.0
   endif

   !subdomain stats
   do i=1,LVT_rc%ntslocs   
      sum_sd_v(i) = 0 
      nsum_sd_v(i) = 0 
      ts_v = 0 
      nts_v = 0

      if(LVT_rc%computeEnsMetrics.eq.1.and.&
           nsize.eq.LVT_LIS_rc(1)%ntiles) then 
         do t=1,LVT_LIS_rc(1)%ntiles
            do m=1,LVT_rc%nensem
               c = LVT_LIS_domain(1)%tile(t)%col
               r = LVT_LIS_domain(1)%tile(t)%row
               if(c.ge.LVT_TSobj(i)%ts_cindex1.and.&
                    c.le.LVT_TSobj(i)%ts_cindex2.and.&
                    r.ge.LVT_TSobj(i)%ts_rindex1.and.&
                    r.le.LVT_TSobj(i)%ts_rindex2) then 
                  if(metric(t,m).ne.LVT_rc%udef) then 
                     sum_sd_v(i) = sum_sd_v(i) + metric(t,m)
                     nsum_sd_v(i) = nsum_sd_v(i) + 1
                     nts_v = nts_v + 1
                     ts_v(nts_v) = metric(t,m) 
                     pct_cov_ts(i) = pct_cov_ts(i) + pct_coverage(t,m)
                     npct_cov_ts(i) = npct_cov_ts(i) + 1
                  endif
               endif
            enddo
         enddo
      else
         if(LVT_rc%tsspecstyle.eq.1.or.&
              LVT_rc%tsspecstyle.eq.2.or.&
              LVT_rc%tsspecstyle.eq.3) then 
            do t=1,LVT_rc%ngrid
               c = LVT_domain%grid(t)%col
               r = LVT_domain%grid(t)%row
               do m=1,LVT_rc%nensem
                  if(c.ge.LVT_TSobj(i)%ts_cindex1.and.&
                       c.le.LVT_TSobj(i)%ts_cindex2.and.&
                       r.ge.LVT_TSobj(i)%ts_rindex1.and.&
                       r.le.LVT_TSobj(i)%ts_rindex2) then 
                     if(metric(t,m).ne.LVT_rc%udef) then 
                        sum_sd_v(i) = sum_sd_v(i) + metric(t,m)
                        nsum_sd_v(i) = nsum_sd_v(i) + 1
                        nts_v = nts_v + 1
                        ts_v(nts_v) = metric(t,m)    
                        pct_cov_ts(i) = pct_cov_ts(i) + pct_coverage(t,m)
                        npct_cov_ts(i) = npct_cov_ts(i) + 1
                     endif
                  endif
               enddo
            enddo
         else
            do t=1,LVT_TSobj(i)%npts
               tid = LVT_TSobj(i)%ts_tindex(t)
               do m=1,LVT_rc%nensem
                  if(tid.ne.-1) then 
                     if(metric(tid,m).ne.LVT_rc%udef) then 
                        sum_sd_v(i) = sum_sd_v(i)+metric(tid,m)
                        nsum_sd_v(i) = nsum_sd_v(i) + 1
                        nts_v = nts_v + 1
                        ts_v(nts_v) = metric(tid,m)  
                        pct_cov_ts(i) = pct_cov_ts(i) + pct_coverage(tid,m)
                        npct_cov_ts(i) = npct_cov_ts(i) + 1
                     endif
                  endif
               enddo
            enddo
         endif
      end if
      if(nts_v.gt.1) then 
         call LVT_computeCI(ts_v(1:nts_v),nts_v,&
              LVT_rc%pval_CI,ci_val(i))
      else
         ci_val(i) = LVT_rc%udef
      endif
      if(nsum_sd_v(i).gt.0) then 
         sum_sd_v(i) = sum_sd_v(i)/nsum_sd_v(i)
      else
         sum_sd_v(i) = LVT_rc%udef
      endif
      if(npct_cov_ts(i).gt.0) then
         pct_cov_ts(i) = pct_cov_ts(i)/npct_cov_ts(i)
      else
         pct_cov_ts(i) = 0.0
      endif
   enddo
   write(ftn,*) '---------------------------------------------------------'
   write(ftn,*) 'VAR: ',trim(varname)
   if(l.eq.2) then 
      write(ftn,*) 'Stratified using ',trim(LVT_rc%vname_strat), &
           ' above ',LVT_rc%strat_var_threshold
   elseif(l.eq.3) then 
      write(ftn,*) 'Stratified using ',trim(LVT_rc%vname_strat), &
           ' below ',LVT_rc%strat_var_threshold
   endif
   write(ftn,*) '---------------------------------------------------------'
   write(ftn,fmt='(a10)',advance='no') 'ALL: '
   write(ftn,fmt='(E14.5,a5,E14.5,E14.5)') sum_v, &
        ' +/- ',sum(ci(:))/LVT_rc%nensem, pct_cov_avg
   do i=1,LVT_rc%ntslocs
      write(ftn, fmt='(a10)',advance='no') trim(LVT_TSobj(i)%tslocname)//':'
      if(ci_val(i).ne.LVT_rc%udef) then 
         write(ftn, '(E14.5,a5,E14.5,E14.5)') sum_sd_v(i), &
              ' +/- ',ci_val(i), pct_cov_ts(i)
      else
         write(ftn, '(E14.5,a5,a14,E14.5)') sum_sd_v(i), &
              ' +/- ','     -       ', pct_cov_ts(i)
      endif
   enddo
   deallocate(sum_sd_v)
   deallocate(nsum_sd_v)
   deallocate(ci_val)
   deallocate(ts_v)
   deallocate(pct_cov_ts)
   deallocate(npct_cov_ts)

  end subroutine LVT_writeSummaryStats

!BOP
! 
! !ROUTINE: LVT_writeSummaryStats2
! \label{LVT_writeSummaryStats2}
!
! !LVT_writeSummaryStats
  subroutine LVT_writeSummaryStats2(ftn,metricname, nsize, metric, npts, &
       varname, ci)

   implicit none
! !ARGUMENTS:   
   integer             :: ftn
   character(len=*)    :: metricname
   integer             :: nsize
   real                :: metric(nsize, LVT_rc%nensem)
   integer             :: npts(nsize, LVT_rc%nensem)
   character(len=*)    :: varname
   real                :: ci(LVT_rc%nensem)
   integer             :: dummy
! 
! !DESCRIPTION: 
!  This routine outputs a summary of various statistics computed 
!  during the LVT analysis. The statistics are computed for the whole
!  analysis domain and for each time series location specified. 
!
!  The arguments are: 
!  \begin{description}
!  \item[ftn]   unit number of the file
!  \item[metricname]  name of the metric 
!  \item[metric]    array containing values of the computed metric
!  \item[npts] array containing the number of counts used in computing 
!              the metric
!  \item[varname]  name of the variable being written out
!  \item[ci]  confidence interval associated with the computed metric values
!  \end{description}
!EOP  
   
   integer             :: i,c,r,m,t
   real                :: sum_v
   integer             :: nsum_v
   real                :: wsum_v
   integer             :: wnsum_v
   real, allocatable       :: sum_sd_v(:)
   integer, allocatable    :: nsum_sd_v(:)
   real, allocatable       :: wsum_sd_v(:)
   integer, allocatable    :: wnsum_sd_v(:)

   real, allocatable       :: ts_v(:)
   integer                 :: nts_v
   real,   allocatable     :: ci_val(:)


   allocate(sum_sd_v(LVT_rc%ntslocs))
   allocate(nsum_sd_v(LVT_rc%ntslocs))
   allocate(ci_val(LVT_rc%ntslocs))
   if(LVT_rc%computeEnsMetrics.eq.1.and.&
        nsize.eq.LVT_LIS_rc(1)%ntiles) then 
      allocate(ts_v(LVT_LIS_rc(1)%ntiles))
   else
      allocate(ts_v(nsize))
   endif

   sum_v = 0 
   nsum_v = 0
   sum_sd_v = 0 
   nsum_sd_v = 0 
   ci_val = 0 
   
   do t=1,nsize
      do m=1,LVT_rc%nensem
         if(metric(t,m).ne.LVT_rc%udef) then 
            sum_v = sum_v + metric(t,m)
            nsum_v = nsum_v + 1
         endif
      enddo
   enddo

   if(nsum_v.gt.0) then 
      sum_v = sum_v/nsum_v
   else
      sum_v = LVT_rc%udef
   endif
   !subdomain stats
   do i=1,LVT_rc%ntslocs   
      sum_sd_v(i) = 0 
      nsum_sd_v(i) = 0 
      ts_v = 0 
      nts_v = 0
      if(LVT_rc%computeEnsMetrics.eq.1.and.&
           nsize.eq.LVT_LIS_rc(1)%ntiles) then 
         do t=1,LVT_LIS_rc(1)%ntiles
            do m=1,LVT_rc%nensem
               c = LVT_LIS_domain(1)%tile(t)%col
               r = LVT_LIS_domain(1)%tile(t)%row
               if(c.ge.LVT_TSobj(i)%ts_cindex1.and.&
                    c.le.LVT_TSobj(i)%ts_cindex2.and.&
                    r.ge.LVT_TSobj(i)%ts_rindex1.and.&
                    r.le.LVT_TSobj(i)%ts_rindex2) then 
                  if(metric(t,m).ne.LVT_rc%udef) then 
                     sum_sd_v(i) = sum_sd_v(i) + metric(t,1)
                     nsum_sd_v(i) = nsum_sd_v(i) + 1
                     nts_v = nts_v + 1
                     ts_v(nts_v) = metric(t,m)                   
                  endif
               endif
            enddo
         enddo
      else
         do t=1,LVT_rc%ngrid
            c = LVT_domain%grid(t)%col
            r = LVT_domain%grid(t)%row
            if(c.ge.LVT_TSobj(i)%ts_cindex1.and.&
                 c.le.LVT_TSobj(i)%ts_cindex2.and.&
                 r.ge.LVT_TSobj(i)%ts_rindex1.and.&
                 r.le.LVT_TSobj(i)%ts_rindex2) then 
               do m=1,LVT_rc%nensem
                  if(metric(t,m).ne.LVT_rc%udef) then 
                     sum_sd_v(i) = sum_sd_v(i) + metric(t,m)
                     nsum_sd_v(i) = nsum_sd_v(i) + 1
                     nts_v = nts_v + 1
                     ts_v(nts_v) = metric(t,m)                   
                  endif
               enddo
            endif
         enddo
      endif
      if(nts_v.gt.1) then 
         call LVT_computeCI(ts_v(1:nts_v),nts_v,&
              LVT_rc%pval_CI,ci_val(i))
      else
         ci_val(i) = LVT_rc%udef
      endif
      if(nsum_sd_v(i).gt.0) then 
         sum_sd_v(i) = sum_sd_v(i)/nsum_sd_v(i)
      else
         sum_sd_v(i) = LVT_rc%udef
      endif
   enddo
   write(ftn,*) '---------------------------------------------------------'
   write(ftn,*) 'VAR: ',trim(varname)
   write(ftn,*) '---------------------------------------------------------'
   write(ftn,fmt='(a10)',advance='no') 'ALL: '
   write(ftn,fmt='(E14.5,a5,E14.5,I14)') sum_v, &
        ' +/- ',ci, nsum_v
   do i=1,LVT_rc%ntslocs
      write(ftn, fmt='(a10)',advance='no') trim(LVT_TSobj(i)%tslocname)//':'
      if(ci_val(i).ne.LVT_rc%udef) then 
         write(ftn, '(E14.5,a5,E14.5,I14)') sum_sd_v(i), &
              ' +/- ',ci_val(i), nsum_sd_v(i)
      else
         write(ftn, '(E14.5,a5,a14,I14)') sum_sd_v(i), &
              ' +/- ','     -       ', nsum_sd_v(i)
      endif
   enddo
   deallocate(sum_sd_v)
   deallocate(nsum_sd_v)
   deallocate(ci_val)
   deallocate(ts_v)

 end subroutine LVT_writeSummaryStats2


!BOP
! 
! !ROUTINE: createOutputFile
! \label(createOutputFile)
!
! !INTERFACE:
  subroutine createOutputFile(metric,pass)
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
! 
!EOP

    type(LVT_metricEntry)      :: metric
    integer                    :: pass
    
    character(len=12)  :: cdate
    character(len=4)   :: cdate1
    integer            :: iret
    character*500      :: fname_total

    if(pass.eq.metric%npass.and.metric%selectOpt.eq.1) then 
       call system("mkdir -p "//trim(LVT_rc%statsodir))
       write(unit=cdate,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2)') LVT_rc%yr, LVT_rc%mo, &
            LVT_rc%da, LVT_rc%hr, LVT_rc%mn 
       
       write(unit=cdate1,fmt='(a2,i2.2)') '.d',LVT_rc%nnest
       if((LVT_rc%lvt_out_format).eq."binary") then 
          fname_total = trim(LVT_rc%statsodir)//&
               '/LVT_'//trim(metric%short_name)//'_'//&
               'FINAL.'//cdate//cdate1//'.gs4r'
       elseif((LVT_rc%lvt_out_format).eq."netcdf") then 
          fname_total = trim(LVT_rc%statsodir)//&
               '/LVT_'//trim(metric%short_name)//'_'//&
               'FINAL.'//cdate//cdate1//'.nc'
       endif
       
       if(metric%selectOpt.eq.1) then 
          if((LVT_rc%lvt_out_format).eq."binary") then 
             metric%ftn_total = LVT_getNextUnitNumber()
             open(metric%ftn_total,file=trim(fname_total),&
                  form='unformatted')
          elseif((LVT_rc%lvt_out_format).eq."netcdf") then 
#if (defined USE_NETCDF4)
             iret = nf90_create(path=trim(fname_total), &
                  cmode =nf90_hdf5, &
                  ncid = metric%ftn_total)
#endif
#if (defined USE_NETCDF3)
             iret = nf90_create(path=trim(fname_total), &
                  cmode =nf90_clobber, &
                  ncid = metric%ftn_total)
#endif
          endif
       endif
    endif
  end subroutine createOutputFile

!BOP
! 
! !ROUTINE: createTSfiles
! \label(createTSfiles)
!
! !INTERFACE:
  subroutine createTSfiles(pass, metric)
! 
! !USES:   
    use LVT_timeMgrMod,  only : LVT_tick

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
! 
!EOP

    integer                 :: pass
    type(LVT_metricEntry)   :: metric

    character*500           :: fname_ts
    character(len=12)       :: cdate
    character(len=4)        :: cdate1
    integer                 :: status,iret
    type(ESMF_Time)         :: currTime
    type(ESMF_TimeInterval) :: ts
    integer                 :: yr, mo, da, hr, mn, ss,doy
    real*8                  :: time
    real                    :: gmt


    if(metric%timeOpt.eq.1.and.metric%writeTS.eq.1.and.&
         metric%npass.eq.pass) then 

       yr = LVT_rc%yr
       mo = LVT_rc%mo
       da = LVT_rc%da
       hr = LVT_rc%hr
       mn = LVT_rc%mn
       ss = LVT_rc%ss

       write(unit=cdate,fmt='(i4.4,i2.2,i2.2,i2.2,i2.2)') &
            yr, mo, da, hr, mn
              
       write(unit=cdate1,fmt='(a2,i2.2)') '.d',LVT_rc%nnest
       
       if((LVT_rc%lvt_out_format).eq."binary") then 
          fname_ts = trim(LVT_rc%statsodir)//&
               '/'//trim(metric%short_name)//'_'//&
               'TS.'//cdate//cdate1//'.gs4r'
          metric%ftn_ts = LVT_getNextUnitNumber()             
          open(metric%ftn_ts,file=trim(fname_ts),form='unformatted')
       elseif((LVT_rc%lvt_out_format).eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
          fname_ts = trim(LVT_rc%statsodir)//&
               '/'//trim(metric%short_name)//'_'//&
               'TS.'//cdate//cdate1//'.nc'
#if(defined USE_NETCDF4) 
          iret = nf90_create(path=(fname_ts),cmode =nf90_hdf5, &
               ncid = metric%ftn_ts)
#endif

#if(defined USE_NETCDF3) 
          iret = nf90_create(path=(fname_ts),cmode=nf90_clobber, &
               ncid = metric%ftn_ts)
#endif
#endif
       endif
    endif
  end subroutine CreateTSfiles

!BOP
! 
! !ROUTINE: finalizeTSfile
! \label(finalizeTSfile)
!
! !INTERFACE:
  subroutine finalizeTSfile(pass, metric)
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
! 
!EOP

    integer               :: pass
    type(LVT_metricEntry) :: metric
    integer               :: iret

    if(metric%timeOpt.eq.1.and.metric%writeTS.eq.1.and.&
         pass.eq.metric%npass) then 
       if((LVT_rc%lvt_out_format).eq."binary") then 
          call LVT_releaseUnitNumber(metric%ftn_ts)
       elseif((LVT_rc%lvt_out_format).eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
          iret = nf90_close(metric%ftn_ts)
#endif
       endif
    endif
    
  end subroutine FinalizeTSfile


!BOP
! 
! !ROUTINE: outputTimeSeriesStats
! \label{outputTimeSeriesStats}
!
! !INTERFACE: 
  subroutine outputTimeSeriesStats(pass)
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
! 
!EOP
!BOP
! 
  

    implicit none
    integer :: pass
! 
! !DESCRIPTION: 
!   This subroutine invokes the calls to write specified temporal statistics to 
!   a file on disk. 
!EOP
    integer :: count
    integer :: index

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

    count = 1

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)

    do while(associated(model))
       call writeSingleHeaderEntry(pass, count,&
            model, obs, stats)
       model  => model%next
       stats => stats%next
       obs => obs%next
    enddo

    call closeTSHeaderEntries(pass)

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)

    do while(associated(model))
       call writeSingleEntry(pass,&
            model, obs, stats)

       model    => model%next
       obs      => obs%next
       stats    => stats%next
    enddo

  end subroutine OutputTimeSeriesStats


!BOP
! 
! !ROUTINE:
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
! 
!EOP
!BOP
! 
! !ROUTINE: writeSingleHeaderEntry
! \label{writeSingleHeaderEntry}
! 
! !INTERFACE:   
  subroutine writeSingleHeaderEntry(pass, count, model, obs, stats)
! !USES: 
    use LVT_historyMod, only  : LVT_writevar_data_header, LVT_writevar_gridded
    implicit none

! !ARGUMENTS: 
    integer                 :: count
    integer                 :: pass
    type(LVT_metadataEntry) :: model
    type(LVT_metadataEntry) :: obs
    type(LVT_statsEntry)    :: stats

! 
! !DESCRIPTION:
!  This routine writes the header information into each statistics file. 
!
!EOP

    integer             :: k 
    integer             :: m
    integer             :: t   
    character*500       :: short_name_ds1, short_name_ds2
    character*500       :: long_name, standard_name, units

    if(stats%selectOpt.eq.1) then 
       do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
          if(LVT_metricsPtr(m)%metricEntryPtr%timeOpt.eq.1.and.&
               LVT_metricsPtr(m)%metricEntryPtr%writeTS.eq.1.and.&
               pass.eq.LVT_metricsPtr(m)%metricEntryPtr%npass) then           

             if(LVT_metricsPtr(m)%metricEntryPtr%obsData) then 
                short_name_ds1    = trim(model%short_name)//"_from_"//&
                     trim(stats%short_name)//"_ds1"
                short_name_ds2    = trim(obs%short_name)//"_from_"//&
                     trim(stats%short_name)//"_ds2"
                long_name     = model%long_name
                standard_name = model%standard_name
                units         = model%units
             else
                short_name_ds1    = stats%short_name
                long_name     = stats%long_name
                standard_name = stats%standard_name
                units         = stats%units
             endif
             
             call LVT_writevar_data_header(&
                  LVT_metricsPtr(m)%metricEntryPtr%ftn_ts, &
                  LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                  short_name_ds1,&
                  standard_name,&
                  long_name,&
                  units,&
                  stats%vid_ts(m,1), model%selectNlevs,&
                  LVT_metricsPtr(m)%metricEntryPtr%nLevs,count)
             
             call LVT_writevar_data_header(&
                  LVT_metricsPtr(m)%metricEntryPtr%ftn_ts, &
                  LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                  "COUNT_"//trim(short_name_ds1), &
                  "COUNT_"//trim(standard_name),&
                  "Number of points of "//trim(long_name),&
                  "-", stats%vid_count_ts(m,1),&
                  model%selectNlevs, LVT_metricsPtr(m)%metricEntryPtr%nLevs,count+1)  

             if(LVT_metricsPtr(m)%metricEntryPtr%obsdata) then 
                call LVT_writevar_data_header(&
                     LVT_metricsPtr(m)%metricEntryPtr%ftn_ts, &
                     LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                     trim(short_name_ds2),&
                     trim(obs%standard_name),&
                     "Observations of "//trim(obs%long_name), & 
                     trim(stats%units),&
                     stats%vid_ts(m,2), model%selectNlevs,&
                     LVT_metricsPtr(m)%metricEntryPtr%nLevs,count+1)
                call LVT_writevar_data_header(&
                     LVT_metricsPtr(m)%metricEntryPtr%ftn_ts, &
                     LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                     "COUNT_"//trim(short_name_ds2), &
                     "COUNT_"//trim(obs%standard_name), &
                     "Number of observation points of "//trim(obs%long_name), &
                     "-",stats%vid_count_ts(m,2),&
                     model%selectNlevs, LVT_metricsPtr(m)%metricEntryPtr%nLevs,&
                     count+1)         
             endif             
          endif
       enddo
       count = count+1
    endif
    
  end subroutine writeSingleHeaderEntry

!BOP
! 
! !ROUTINE:
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
! 
!EOP
!BOP
! !ROUTINE: closeTSHeaderEntries
! \label{closeTSHeaderEntries}
!
! !INTERFACE: 
  subroutine closeTSHeaderEntries(pass)
! !USES: 
    use LVT_historyMod, only : LVT_close_data_header

    integer              :: pass
    integer              :: m

    do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
       if(LVT_metricsPtr(m)%metricEntryPtr%timeOpt.eq.1.and.&
            LVT_metricsPtr(m)%metricEntryPtr%writeTS.eq.1.and.&
            pass.eq.LVT_metricsPtr(m)%metricEntryPtr%npass) then 
          call LVT_close_data_header(LVT_metricsPtr(m)%metricEntryPtr%ftn_ts)
       endif
    enddo
    
  end subroutine closeTSHeaderEntries

!BOP
! 
! !ROUTINE:
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
! 
!EOP
!BOP
! 
! !ROUTINE: writeSingleEntry
! \label{writeSingleEntry}
! 
! !INTERFACE:   
  subroutine writeSingleEntry(pass, model,obs,stats)
! !USES: 
    use LVT_historyMod, only  : LVT_writevar_gridded

    implicit none

! !ARGUMENTS: 
    integer                 :: pass
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    integer                 :: vlevels
! 
! !DESCRIPTION:
!  This routine writes the specified set of temporal statistics for a 
!  single variable. 
!EOP

    integer             :: k 
    integer             :: t
    integer             :: l
    integer             :: m

    do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
       if(LVT_metricsPtr(m)%metricEntryPtr%timeOpt.eq.1.and.&
            LVT_metricsPtr(m)%metricEntryPtr%writeTS.eq.1.and.&
            pass.eq.LVT_metricsPtr(m)%metricEntryPtr%npass) then 
          call writemetricentry(m,pass,0,model%selectNlevs,stats,obs)
       endif
    enddo

  end subroutine writeSingleEntry

!BOP
! 
! !ROUTINE:
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
! 
!EOP
!BOP
! !ROUTINE: outputFinalStats
! \label{outputFinalStats}
! 
! !INTERFACE: 
  subroutine outputFinalStats(pass)

    implicit none
! !ARGUMENTS: 
    integer          :: pass
! 
! !DESCRIPTION: 
!  This routine writes the final set of statistics at the end of the analysis. 
!
!EOP
    integer                          :: count
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

    count = 1

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)

    do while(associated(model))
       call writeFinalSingleHeaderEntry(pass,count,&
            model, obs, stats)
       model  => model%next
       obs => obs%next
       stats => stats%next
    enddo
       
    call closeFinalHeaderEntries(pass)

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)

    do while(associated(model))

       call writeFinalSingleEntry(pass,model, obs, stats)

       model  => model%next
       obs    => obs%next
       stats  => stats%next
    enddo

  end subroutine OutputFinalStats


!BOP
! 
! !ROUTINE:
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
! 
!EOP
!BOP
! !ROUTINE: writeFinalSingleHeaderEntry
! \label{writeFinalSingleHeaderEntry}
! 
! !INTERFACE:   
  subroutine writeFinalSingleHeaderEntry(pass,count, model, obs, stats)
! !USES: 

    use LVT_historyMod, only  : LVT_writevar_data_header, LVT_writevar_gridded

    implicit none
! !ARGUMENTS: 
    integer                 :: pass
    integer                 :: count
    type(LVT_metadataEntry) :: model
    type(LVT_metadataEntry) :: obs
    type(LVT_statsEntry)    :: stats

!
! !DESCRIPTION: 
!  This routine writes the specified set of statistics for a 
!  single variable at the end of the analysis. 
!EOP

    integer             :: k,m,i
    integer             :: kk
    character*100       :: short_name_ds1, short_name_ds2
    character*100       :: long_name, standard_name, units

    kk = count
    if(stats%selectOpt.eq.1) then 
       do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
          if(pass.eq.LVT_metricsPtr(m)%metricEntryPtr%npass) then 
             if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.eq.1) then
                if(LVT_metricsPtr(m)%metricEntryPtr%customNames) then 

                   do i=1,LVT_metricsPtr(m)%metricEntryPtr%nFields
                      call LVT_writevar_data_header(&
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                           trim(stats%short_name)//'_'//&
                           trim(LVT_metricsPtr(m)%metricEntryPtr%mName(i)),&
                           (stats%standard_name),&
                           (stats%long_name),&
                           (stats%units),&
                           stats%vid_total(m,(i-1)*2+1), model%selectNlevs, &
                           LVT_metricsPtr(m)%metricEntryPtr%nLevs, &
                           count+(i-1))
                      call LVT_writevar_data_header(&
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                           'COUNT_'//trim(stats%short_name)//'_'//&
                           trim(LVT_metricsPtr(m)%metricEntryPtr%mName(i)), &
                           trim(stats%standard_name),&
                           trim(stats%long_name),&
                           stats%units, stats%vid_count_total(m,(i-1)*2+1),&
                           model%selectNlevs, &
                           LVT_metricsPtr(m)%metricEntryPtr%nLevs,&
                           count+1+(i-1))        
                      if(LVT_metricsPtr(m)%metricEntryPtr%obsData) then 
                         call LVT_writevar_data_header(&
                              LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                              LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                              "OBS_"//trim(stats%short_name)//'_'//&
                              trim(LVT_metricsPtr(m)%metricEntryPtr%mName(i)),&
                              trim(obs%standard_name), &
                              "Observations of "//trim(obs%long_name), &
                              trim(stats%units),&
                              stats%vid_total(m,(i-1)*2+2), model%selectNlevs,&
                              LVT_metricsPtr(m)%metricEntryPtr%nLevs,&
                              count+2+(i-1))
                         call LVT_writevar_data_header(&
                              LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                              LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                              "COUNT_OBS_"//trim(stats%short_name)//'_'//&
                              trim(LVT_metricsPtr(m)%metricEntryPtr%mName(i)), &
                              "COUNT_OBS_"//trim(stats%standard_name),&
                              "Number of observation points of "//&
                              trim(stats%long_name),&
                              stats%units, stats%vid_count_total(m,(i-1)*2+2),&
                              model%selectNlevs, &
                              LVT_metricsPtr(m)%metricEntryPtr%nLevs,&
                              count+3+(i-1))        
                         
                      endif
                   enddo
                else
                   if(LVT_metricsPtr(m)%metricEntryPtr%obsData) then 
                      short_name_ds1  = trim(model%short_name)//"_from_"//&
                           trim(stats%short_name)//"_ds1"
                      short_name_ds2  = trim(obs%short_name)//"_from_"//&
                           trim(stats%short_name)//"_ds2"
                      long_name     = model%long_name
                      standard_name = model%standard_name
                      units         = model%units
                   else
                      short_name_ds1    = stats%short_name
                      long_name     = stats%long_name
                      standard_name = stats%standard_name
                      units         = stats%units
                   endif
                   
                   call LVT_writevar_data_header(&
                        LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                        LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                        trim(short_name_ds1), &
                        trim(standard_name), &
                        long_name, &
                        units,&
                        stats%vid_total(m,1), model%selectNlevs, &
                        LVT_metricsPtr(m)%metricEntryPtr%nLevs,&
                        count)
                   kk= kk+1
                   call LVT_writevar_data_header(&
                        LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                        LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                        "COUNT_"//trim(short_name_ds1),&
                        "COUNT_"//trim(standard_name),&
                        "Number of points in "//trim(long_name),&
                        "-",stats%vid_count_total(m,1), model%selectNlevs,&
                        LVT_metricsPtr(m)%metricEntryPtr%nLevs,kk)
                   kk= kk+1
                   if(LVT_metricsPtr(m)%metricEntryPtr%stdevFlag) then 
                      call LVT_writevar_data_header(&
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                           "SD_"//(short_name_ds1),&
                           "SD_"//(standard_name),&
                           "standard deviation "//(long_name),&
                           (units),&
                           stats%vid_stdev_total(m,1), model%selectNlevs,&
                           LVT_metricsPtr(m)%metricEntryPtr%nLevs,kk)
                      kk= kk+1
                      call LVT_writevar_data_header(&
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                           "COUNT_SD_"//trim(short_name_ds1), &
                           "COUNT_SD_"//trim(standard_name),&
                           "Number of points (SD) "//trim(long_name),&
                           "-", stats%vid_count_stdev_total(m,1),&
                           model%selectNlevs, &
                           LVT_metricsPtr(m)%metricEntryPtr%nLevs,kk)         
                   endif
                   if(LVT_metricsPtr(m)%metricEntryPtr%obsdata) then 
                      kk = kk + 1
                      call LVT_writevar_data_header(&
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                           trim(short_name_ds2), &
                           trim(obs%standard_name), &
                           "Observations of "//trim(obs%long_name), &
                           trim(stats%units),stats%vid_total(m,2), model%selectNlevs,&
                           LVT_metricsPtr(m)%metricEntryPtr%nLevs,kk)
                      kk= kk+1
                      call LVT_writevar_data_header(&
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                           LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                           "COUNT_"//trim(short_name_ds2),&
                           "COUNT_"//trim(obs%standard_name),&
                           "Number of observation points of "//trim(obs%long_name),&
                           "-",stats%vid_count_total(m,2), model%selectNlevs,&
                           LVT_metricsPtr(m)%metricEntryPtr%nLevs,kk)        
                   endif
                   
                   if(LVT_metricsPtr(m)%metricEntryPtr%computeSC.eq.1) then 
                      do k=1, LVT_rc%nasc
                         kk = kk+k
                         call LVT_writevar_data_header(&
                              LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                              LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                              trim(short_name_ds1)//'_'//&
                              trim(LVT_rc%scname(k)),&
                              trim(standard_name)//'_'//&
                              trim(LVT_rc%scname(k)),&
                              'Seasonal cycle '//trim(long_name)//'_'//&
                              trim(LVT_rc%scname(k)),&
                              trim(stats%units),stats%vid_sc_total(k,m,1), &
                              model%selectNlevs,&
                              LVT_metricsPtr(m)%metricEntryPtr%nLevs,kk)
                      enddo
                      if(LVT_metricsPtr(m)%metricEntryPtr%obsdata) then 
                         do k=1, LVT_rc%nasc
                            kk = kk+k
                            call LVT_writevar_data_header(&
                                 LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                                 LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                                 trim(short_name_ds2)//"_"//&
                                 trim(LVT_rc%scname(k)),&
                                 "DS2_"//trim(standard_name)//'_'//&
                                 trim(LVT_rc%scname(k)),&
                                 "DS2_"//trim(long_name)//'_'//&
                                 trim(LVT_rc%scname(k)),&
                                 "-",stats%vid_sc_total(k,m,2),&
                                 model%selectNlevs,&
                                 LVT_metricsPtr(m)%metricEntryPtr%nLevs,kk)
                         enddo
                      endif
                   endif
                   if(LVT_metricsPtr(m)%metricEntryPtr%computeADC.eq.1) then 
                      do k=1, LVT_rc%nadc
                         kk = kk+k
                         call LVT_writevar_data_header(&
                              LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                              LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                              trim(short_name_ds1)//"_"//&
                              trim(LVT_rc%adcname(k)),&
                              trim(standard_name)//'_'//&
                              trim(LVT_rc%adcname(k)),&
                              trim(long_name)//'_'//&
                              trim(LVT_rc%adcname(k)),&
                              trim(stats%units),stats%vid_adc_total(k,m,1), &
                              1,LVT_metricsPtr(m)%metricEntryPtr%nLevs,kk)
                      enddo
                      if(LVT_metricsPtr(m)%metricEntryPtr%obsdata) then 
                         do k=1, LVT_rc%nadc
                            kk = kk+k
                            call LVT_writevar_data_header(&
                                 LVT_metricsPtr(m)%metricEntryPtr%ftn_total, &
                                 LVT_metricsPtr(m)%metricEntryPtr%ftn_meta_out,&
                                 trim(short_name_ds2)//'_'//&
                                 trim(LVT_rc%adcname(k)),&
                                 "DS2_"//trim(stats%standard_name)//'_'//&
                                 trim(LVT_rc%adcname(k)),&
                                 "DS2_"//trim(stats%long_name)//'_'//&
                                 trim(LVT_rc%adcname(k)),&
                                 "-",stats%vid_adc_total(k,m,2), &
                                 1,LVT_metricsPtr(m)%metricEntryPtr%nLevs,kk)
                         enddo
                      endif
                   endif
                endif
             endif
          endif           
       enddo
       count = count + 1
    endif
    
  end subroutine writeFinalSingleHeaderEntry

!BOP
! 
! !ROUTINE:
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
! 
!EOP
!BOP
! !ROUTINE: writeFinalSingleEntry
! \label{writeFinalSingleEntry}
! 
! !INTERFACE:   
  subroutine writeFinalSingleEntry(pass,model,obs,stats)
! !USES: 
    use LVT_historyMod, only  : LVT_writevar_gridded

    implicit none
! !ARGUMENTS: 
    integer             :: pass
    type(LVT_statsEntry)    :: stats
    type(LVT_metadataEntry) :: obs
    type(LVT_metadataEntry) :: model
!
! !DESCRIPTION: 
!  This routine writes the specified set of statistics for a 
!  single variable at the end of the analysis. 
!EOP

    integer             :: k 
    integer             :: c,r,m
    integer             :: t,i
    
    if(stats%selectOpt.eq.1) then 
       do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
          call writemetricentry(m,pass,1,model%selectNlevs,stats,obs)
       enddo
    endif
    
  end subroutine writeFinalSingleEntry

  
!BOP
! 
! !ROUTINE:
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
! 
!EOP
!BOP
! !ROUTINE: closeFinalHeaderEntries
! \label{closeFinalHeaderEntries}
! 
! !INTERFACE:   
  subroutine closeFinalHeaderEntries(pass)
! !USES: 
    use LVT_historyMod, only  : LVT_close_data_header

    implicit none
! !ARGUMENTS: 
    integer,   intent(in)   :: pass
!
! !DESCRIPTION: 
!  This routine writes the specified set of statistics for a 
!  single variable at the end of the analysis. 
!EOP
    integer                 :: m

    do m=LVT_rc%metric_sindex,LVT_rc%metric_eindex
       if(pass.eq.LVT_metricsPtr(m)%metricEntryPtr%npass) then     
          if(LVT_metricsPtr(m)%metricEntryPtr%selectOpt.eq.1) then
             call LVT_close_data_header(LVT_metricsPtr(m)%metricEntryPtr%ftn_total)
          endif
       endif
    enddo
    
  end subroutine closeFinalHeaderEntries

!BOP
! 
! !ROUTINE: finalizeOutputfile
! \label(finalizeOutputfile)
!
! !INTERFACE:
  subroutine finalizeOutputfile(metric,pass)
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
! 
!EOP

    type(LVT_metricEntry) :: metric
    integer               :: pass
    integer               :: iret

    if(pass.eq.metric%npass) then 
       if(metric%selectOpt.eq.1) then 
          if((LVT_rc%lvt_out_format).eq."binary") then 
             call LVT_releaseUnitNumber(metric%ftn_total)
          elseif((LVT_rc%lvt_out_format).eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
             iret = nf90_close(metric%ftn_total)
#endif
          endif
       endif
    endif
  end subroutine FinalizeOutputfile


  subroutine LVT_checkForDerivedVariableDependendency()
    
    use LVT_LISoutputHandlerMod

    integer :: source
    integer :: npass

    npass = 0
    do source =1,2

       if (trim(LVT_rc%obssource(source)).eq."LIS output") then 
          if ((LVT_MOC_RELSMC(source).gt.0).and.&
               (LVT_LIS_MOC_RELSMC(source).lt.0)) then
             npass = 2
          endif
       else !non-LIS data 
          if ((LVT_MOC_RELSMC(source).gt.0).and.&
               (trim(LVT_rc%obssource(source)).ne."none")) then 
             npass = 2
          endif
       endif
    enddo

    LVT_rc%pass = max(npass, LVT_rc%pass)

  end subroutine LVT_checkForDerivedVariableDependendency

end module LVT_statsMod
