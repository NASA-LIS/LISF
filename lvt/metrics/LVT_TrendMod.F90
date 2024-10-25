!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! 
! !MODULE: LVT_TrendMod
! \label(LVT_TrendMod)
!
! !INTERFACE:
module LVT_TrendMod
! 
! !USES: 
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_timeMgrMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod

  implicit none

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the computations required to compute trends 
!  of desired variables. The output is the slope of the trendline
!  across the time averaged datastream values (at the metric computation
!  frequency). A mann kendall test is applied to filter
!  out stastistically insignificant values.  
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  1 June 2017    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initTrend
  public :: LVT_diagnoseTrend
  public :: LVT_computeTrend
  public :: LVT_writeMetric_Trend
  public :: LVT_resetMetric_Trend
  public :: LVT_writerestart_Trend
  public :: LVT_readrestart_Trend
!EOP
  
  private

contains
  subroutine LVT_initTrend(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This routine allocates and initializes the datastructures required for
!  the trend calculations.  
! 
!EOP
    integer                 :: m

    allocate(stats%trend(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%trend(m)%model_value_total(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%ntavgs))
          stats%trend(m)%model_value_total = LVT_rc%udef
          allocate(stats%trend(m)%count_model_value_total(LVT_rc%ngrid, selectNlevs(1)))
          stats%trend(m)%count_model_value_total = 0.0

          allocate(stats%trend(m)%model_slope(LVT_rc%ngrid, selectNlevs(1)))
          stats%trend(m)%model_slope = LVT_rc%udef

          allocate(stats%trend(m)%model_sig(LVT_rc%ngrid, selectNlevs(1)))
          stats%trend(m)%model_slope = LVT_rc%udef

          allocate(stats%trend(m)%model_value_ci(selectNlevs(1)))
          stats%trend(m)%model_value_ci = LVT_rc%udef
          
          if(LVT_rc%obssource(2).ne."none") then 
             if(selectNlevs(2).ge.1) then 
                allocate(stats%trend(m)%obs_value_total(LVT_rc%ngrid, selectNlevs(2), &
                     LVT_rc%ntavgs))
                stats%trend(m)%obs_value_total = LVT_rc%udef
                allocate(stats%trend(m)%count_obs_value_total(LVT_rc%ngrid, selectNlevs(2)))
                stats%trend(m)%count_obs_value_total = 0.0

                allocate(stats%trend(m)%obs_slope(LVT_rc%ngrid, selectNlevs(2)))
                stats%trend(m)%obs_slope = LVT_rc%udef

                allocate(stats%trend(m)%obs_sig(LVT_rc%ngrid, selectNlevs(2)))
                stats%trend(m)%obs_sig = LVT_rc%udef

                allocate(stats%trend(m)%obs_value_ci(selectNlevs(1)))
                stats%trend(m)%obs_value_ci = LVT_rc%udef
             endif
          endif
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    if(LVT_rc%obssource(2).ne."none") then 
       metric%obsData = .true. 
    else
       metric%obsData = .false. 
    endif
    metric%stdevFlag = .false. 

    if(LVT_metrics%trend%timeOpt.eq.1) then 
       write(LVT_logunit,*) '[ERR] Calculations of trends are only done across the entire'
       write(LVT_logunit,*) '[ERR] analysis time period. Please turn off the flags in '
       write(LVT_logunit,*) '[ERR] metrics attributes table that enables the temporal'
       write(LVT_logunit,*) '[ERR] calculations'
       call LVT_endrun()
    endif

  end subroutine LVT_initTrend

!BOP
! 
! !ROUTINE: LVT_diagnoseTrend
! \label{LVT_diagnoseTrend}
!
! !INTERFACE: 
  subroutine LVT_diagnoseTrend(pass)
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the computations for 
!   calculating the trend of desired variables. This routine will be
!   called at each metric computation frequency for each set of 
!   variables that are selected. 
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelTrend](\ref{diagnoseSingleModelTrend})
!     updates the trend computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer       :: pass
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.LVT_metrics%trend%npass) then
       if(LVT_metrics%trend%selectOpt.eq.1.or.&
            LVT_metrics%trend%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleTrend(model,obs,stats,&
                  LVT_metrics%trend)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseTrend

!BOP
! 
! !ROUTINE: diagnoseSingleTrend
! \label{diagnoseSingleTrend}
!
! !INTERFACE: 
  subroutine diagnoseSingleTrend(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the trend computation of the 
!   specified variable. At each metric computation frequency,
!   the time averaged variables are stored in memory.  
!
!  The arguments are: 
!
!  \begin{description}
!   \item[model] model variable object
!   \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metadataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP
    integer    :: t,k,n,m,m_k,o_k,tind, tindex
    
    type(ESMF_Time)         :: currTime, startTime
    type(ESMF_TimeInterval) :: timestep
    integer                 :: rc

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
!compute tindex 
           
       call ESMF_ClockGet(LVT_clock, startTime = startTime, &
            currTime = currTime, rc=rc)
       call ESMF_TimeIntervalSet(timeStep, s = LVT_rc%tavgInterval, &
            rc=rc)
       tindex = nint((currTime-startTime)/timestep) + 1

       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             do m=1,LVT_rc%nensem
                m_k = k+model%startNlevs -1
                if(model%count(t,m,m_k).gt.0) then
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         stats%trend(m)%model_value_total(t,k,tindex) = &
                              model%value(t,m,m_k)
                         stats%trend(m)%count_model_value_total(t,k) = &
                              stats%trend(m)%count_model_value_total(t,k) + 1
                      endif
                   endif
                
                endif
             end do
          enddo
          do k=1,obs%selectNlevs
             do m=1,LVT_rc%nensem
                o_k = k+obs%startNlevs -1
                if(LVT_rc%obssource(2).ne."none") then
                   if(obs%selectNlevs.ge.1) then
                      if(obs%count(t,m,o_k).gt.0) then 
                         if(metric%selectOpt.eq.1.and.&
                              obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            stats%trend(m)%obs_value_total(t,k,tindex) = &
                                 obs%value(t,m,o_k)
                            stats%trend(m)%count_obs_value_total(t,k) = &
                              stats%trend(m)%count_obs_value_total(t,k) + 1
                            
                         endif

                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine diagnoseSingleTrend

!BOP
! 
! !ROUTINE: LVT_computeTrend
! \label{LVT_computeTrend}
!
! !INTERFACE: 
  subroutine LVT_computeTrend(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the trend values for 
!   desired variables. Since the trend metric is computed only 
!   across the entire analysis time period, this routine is 
!   invoked at the end of the analysis, for each variable. 
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelTrend](\ref{computeSingleModelTrend})
!     computes the trend values for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer               :: pass
    logical               :: alarm
!EOP
    integer     :: i,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    
    if(pass.eq.LVT_metrics%trend%npass) then
       if(LVT_metrics%trend%selectOpt.eq.1.or.&
            LVT_metrics%trend%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%trend%timeOpt.eq.1.and.&
                  LVT_metrics%trend%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%trend%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%trend%ftn_ts_loc(i,1),200,advance='no') &
                           LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                           LVT_rc%hr,'',LVT_rc%mn, '' 
                   enddo
                endif
             endif
          endif
             
200       format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             
             call computeSingleTrend(alarm,model,obs,stats,&
                  LVT_metrics%trend)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%trend%timeOpt.eq.1.and.&
                  LVT_metrics%trend%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%trend%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%trend%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif
             endif
          endif
       endif
    endif
  end subroutine LVT_computeTrend
  

!BOP
! 
! !ROUTINE: computeSingleTrend
! \label{computeSingleTrend}
!
! !INTERFACE: 
  subroutine computeSingleTrend(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the trend values at the end of the 
!  analysis. A mann kendall test is applied to filter out
!  statistically insignificant values. 
! 
!  The arguments are: 
!
!  \begin{description}
!    \item[model] model variable object
!    \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP

    integer  :: t,l,k,m
    
    real,    allocatable :: model_value_avg(:,:,:)
    real,    allocatable :: obs_value_avg(:,:,:)   

       
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(stats%trend(m)%count_model_value_total(t,k).gt.&
                        LVT_rc%obsCountThreshold) then 
                      
                      call mann_kendall(&
                           stats%trend(m)%model_value_total(t,k,:),&
                           LVT_rc%ntavgs,&
                           LVT_rc%udef, &
                           metric%minthreshold, &
                           stats%trend(m)%model_slope(t,k),&
                           stats%trend(m)%model_sig(t,k))
                      
                   else
                      stats%trend(m)%model_slope(t,k) = LVT_rc%udef
                      stats%trend(m)%model_sig(t,k) = LVT_rc%udef
                   endif
                enddo

                do k=1,obs%selectNlevs
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then
                         if(stats%trend(m)%count_obs_value_total(t,k).gt.&
                              LVT_rc%obsCountThreshold) then 

                            call mann_kendall(&
                                 stats%trend(m)%obs_value_total(t,k,:),&
                                 LVT_rc%ntavgs,&
                                 LVT_rc%udef, &
                                 metric%minthreshold, &
                                 stats%trend(m)%obs_slope(t,k),&
                                 stats%trend(m)%obs_sig(t,k))
                         else
                            stats%trend(m)%obs_slope(t,k) = LVT_rc%udef
                            stats%trend(m)%obs_sig(t,k) = LVT_rc%udef
                         endif
                      end if
                   endif
                enddo
             enddo
          enddo

          do m=1,LVT_rc%nensem          
             do k=1,model%selectNlevs                
                call LVT_computeCI(stats%trend(m)%model_slope(:,k),&
                     LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%trend(m)%model_value_ci(k))
             enddo
             if(LVT_rc%obssource(2).ne."none") then 
                do k=1,obs%selectNlevs
                   call LVT_computeCI(stats%trend(m)%obs_slope(:,k),&
                        LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%trend(m)%obs_value_ci(k))
                enddo
             endif
          enddo
#if 0 
!----------------------------------------------------------------------
!  External data based stratification 
!----------------------------------------------------------------------
          if(LVT_rc%data_based_strat.eq.1) then 

             allocate(model_value_avg(LVT_rc%ngrid,&
                  model%vlevels,&
                  LVT_rc%strat_nlevels))
             model_value_avg = 0.0

             do m=1,LVT_rc%nensem        
                do t=1,LVT_rc%ngrid                      
                   do k=1,model%selectNlevs                
                      do l=1, LVT_rc%strat_nlevels
                         model_value_avg(t,k,l) = &
                              model_value_avg(t,k,l) + & 
                              stats%trend(m)%model_value_total(t,k,l)
                      enddo
                   enddo
                enddo
             enddo
             
             do t=1,LVT_rc%ngrid                      
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%strat_nlevels
                      model_value_avg(t,k,l) = &
                           model_value_avg(t,k,l)/LVT_rc%nensem
                   enddo
                enddo
             enddo
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   
                allocate(obs_value_avg(LVT_rc%ngrid,&
                     obs%vlevels,&
                     LVT_rc%strat_nlevels))
                obs_value_avg = 0.0

                do m=1,LVT_rc%nensem        
                   do t=1,LVT_rc%ngrid
                      do k=1,obs%selectNlevs                
                         do l=1, LVT_rc%strat_nlevels
                            obs_value_avg(t,k,l) = &
                                 obs_value_avg(t,k,l) + & 
                                 stats%trend(m)%obs_value_total(t,k,l)
                         enddo
                      enddo
                   enddo
                enddo
                do t=1,LVT_rc%ngrid
                   do k=1,obs%selectNlevs                
                      do l=1, LVT_rc%strat_nlevels                
                         obs_value_avg(t,k,l) = &
                              obs_value_avg(t,k,l)/LVT_rc%nensem
                      enddo
                   enddo
                enddo
                call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                     LVT_rc%ngrid, model_value_avg, &
                     LVT_rc%ngrid, obs_value_avg)
                   
                deallocate(obs_value_avg)
             else
                   
                call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_avg)            
                
             endif
             deallocate(model_value_avg)
          endif
#endif
       endif
    endif
  end subroutine computeSingleTrend

!BOP
! 
! !ROUTINE: LVT_writeMetric_Trend
! \label(LVT_writeMetric_Trend)
!
! !INTERFACE:
  subroutine LVT_writeMetric_Trend(pass,final,vlevels,stats,obs)
! 
! !USES:   
    use LVT_statsMod, only : LVT_writeSummaryStats
    use LVT_pluginIndices
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   
!   This subroutine writes the trend values (slope of the 
!   trendline) to an external file.
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs

    integer                 :: k,l,m,tind

    integer,    allocatable :: count_model_value_total(:,:)
    real,    allocatable    :: model_value_slope(:,:)
    real,    allocatable    :: model_value_ci(:)

    integer, allocatable    :: count_obs_value_total(:,:)
    real,    allocatable    :: obs_value_slope(:,:)
    real,    allocatable    :: obs_value_ci(:)

    if(pass.eq.LVT_metrics%trend%npass) then
       if(final.eq.1) then
          if(pass.eq.LVT_metrics%trend%npass) then
             if(stats%selectOpt.eq.1) then

                allocate(model_value_slope(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(count_model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(model_value_ci(LVT_rc%nensem))

                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   allocate(obs_value_slope(LVT_rc%ngrid,LVT_rc%nensem))
                   allocate(count_obs_value_total(LVT_rc%ngrid,LVT_rc%nensem))
                   allocate(obs_value_ci(LVT_rc%nensem))
                endif
                
                do k=1,vlevels
                   do m=1,LVT_rc%nensem
                      model_value_slope(:,m) = &
                           stats%trend(m)%model_slope(:,k)
                      count_model_value_total(:,m) = &
                           stats%trend(m)%count_model_value_total(:,k)
                      model_value_ci(m) = stats%trend(m)%model_value_ci(k)
                      
                      if(LVT_rc%obssource(2).ne."none".and.&
                           obs%selectNlevs.ge.1) then 
                         obs_value_slope(:,m) = &
                              stats%trend(m)%obs_slope(:,k)
                         count_obs_value_total(:,m) = &
                              stats%trend(m)%count_obs_value_total(:,k)
                         obs_value_ci(m) = stats%trend(m)%obs_value_ci(k)

                      endif
                      if(LVT_metrics%trend%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%trend%ftn_total, &
                              model_value_slope(:,:),&
                              stats%vid_total(LVT_Trendid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%trend%ftn_total, &
                              real(count_model_value_total(:,:)),&
                              stats%vid_count_total(LVT_Trendid,1),k)

                         if(LVT_rc%obssource(2).ne."none".and.&
                              obs%selectNlevs.ge.1) then 
                            call LVT_writevar_gridded(LVT_metrics%trend%ftn_total, &
                                 obs_value_slope(:,:),&
                                 stats%vid_total(LVT_Trendid,2),k)
                            call LVT_writevar_gridded(LVT_metrics%trend%ftn_total, &
                                 real(count_obs_value_total(:,:)),&
                                 stats%vid_count_total(LVT_Trendid,2),k)
                         endif

                         call LVT_writeSummaryStats(&
                              LVT_metrics%trend%ftn_summ,&
                              l,&
                              LVT_metrics%trend%short_name,&
                              LVT_rc%ngrid,&
                              model_value_slope(:,:), &
                              count_model_value_total(:,:),&
                              stats%standard_name,&
                              model_value_ci(:))
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            call LVT_writeSummaryStats(&
                                 LVT_metrics%trend%ftn_summ,&
                                 l,&
                                 LVT_metrics%trend%short_name,&
                                 LVT_rc%ngrid,&
                                 obs_value_slope(:,:), &
                                 count_obs_value_total(:,:),&
                                 "DS2_"//trim(stats%standard_name),&
                                 obs_value_ci(:))
                         endif

                      endif
                   enddo
                enddo

                deallocate(model_value_slope)
                deallocate(count_model_value_total)
                deallocate(model_value_ci)

                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   deallocate(obs_value_slope)
                   deallocate(count_obs_value_total)
                   deallocate(obs_value_ci)
                endif

             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_Trend

!BOP
! 
! !ROUTINE: LVT_resetMetric_Trend
! \label(LVT_resetMetric_Trend)
!
! !INTERFACE:
  subroutine LVT_resetMetric_Trend(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical                :: alarm
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine resets the variables used in the 
!   trend calculation. By design, this routine is 
!   empty. 
!
! !REVISION HISTORY: 
! 
!EOP

  end subroutine LVT_resetMetric_Trend

!BOP
! 
! !ROUTINE: LVT_writerestart_Trend
! 
! !INTERFACE:
  subroutine LVT_writerestart_Trend(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Trend metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,l,m


    if(LVT_metrics%trend%selectOpt.eq.1) then 
       
       print*, 'WARNING: The writerestart method is not implemented for the Trend metric'

    end if
  end subroutine LVT_writerestart_Trend


!BOP
! 
! !ROUTINE: LVT_readrestart_Trend
! 
! !INTERFACE:
  subroutine LVT_readrestart_Trend(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Trend metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    integer              :: k,l,m,index

    if(LVT_metrics%trend%selectOpt.eq.1) then 
       
       print*, 'WARNING: The readrestart method is not implemented for the Trend metric'

    end if
  end subroutine LVT_readrestart_Trend

!BOP
! 
! !ROUTINE: mann_kendall
! 
! !INTERFACE: 
  subroutine mann_kendall(yvar,nvar,missing,threshold, slope,sig)

! !DESCRIPTION: 
!  This subroutine computes the linear trendline slope and 
!  applies a statistical significance test. 
!  
! The arguments are: 
! yvar = single-point 1-D time series array
! nvar = number of times in the time series
! missing = value for missing data
! slope = output trend/slope from Mann-Kendall
! sig = Mann-Kendall significance test
!           (significant = 1; NOT significant = 0)
!EOP

    use LVT_SortingMod
    
    implicit none
    
    integer, intent(in)  :: nvar
    real, intent(in)     :: yvar(nvar)
    real, intent(in)     :: missing
    real, intent(in)     :: threshold
    real, intent(out)    :: slope, sig
    real, allocatable    :: slopes(:)
    real    :: z_090,z_095,z_0975,var_s,c,zval
    integer :: nummiss,npair,count,sgn,z_mk
    integer :: t,n,i,j,m1,m2
    
    nummiss = 0
    do t = 1,nvar
       if (yvar(t).eq.missing) nummiss = nummiss + 1
    enddo
    n = nvar - nummiss
    npair = n*(n-1)/2
    allocate(slopes(npair))
    
! Mann-Kendall trend calculation
    if(npair.gt.0) then 
       count = 1
       sgn = 0
       do i = 1,n-1
          do j = i+1,n
             if ((yvar(j).ne.missing).and.(yvar(i).ne.missing)) then
                if (yvar(j).gt.yvar(i)) then
                   sgn = sgn + 1
                elseif (yvar(j).lt.yvar(i)) then
                   sgn = sgn - 1
                endif
                slopes(count) = (yvar(j)-yvar(i))/((j-i)*1.0)
                count = count + 1
             endif
          enddo
       enddo
    
! Find the median value
!       call LVT_sort(slopes,count)
       call LVT_quicksort(slopes, count, 1, count)
       if (mod(count,2).eq.0) then
          slope = (slopes(count/2)+slopes((count+2)/2))/2.0
       else
          slope = slopes((count+1)/2)
       endif
       
! Calculate the confidence levelx
!      z_0975 = 1.96
!      z_0900 = 1.64
!      z_095 = 2.00
       if(threshold.eq.0.2) then 
          zval = 1.29
       elseif(threshold.eq.0.1) then 
          zval = 1.65
       elseif(threshold.eq.0.05) then 
          zval = 1.96
       endif

!    z_090  = 1.29 ! (alpha=0.2, 80% confidence interval)
!    z_095  = 1.65 ! (alpha=0.1, 90% confidence interval)
!    z_0975 = 1.96 ! (alpha=0.05,95% confidence interval)
       var_s = (n*(n-1)*(2*n+5))/18.0
       c = zval*sqrt(var_s)

       m1 = int((count-c)/2.0)
       m2 = int((count+c)/2.0)+1
       
    ! Mann-Kendall significance test
       if (sgn.gt.0) then
          z_mk = (sgn-1)/sqrt(var_s)
       else
          z_mk = (sgn+1)/sqrt(var_s)
       endif
       sig = 0
       !      if (abs(z_mk).gt.z_0900) then
       if(threshold.ne.0) then 
          if (abs(z_mk).ge.zval) then
             sig = 1
          endif
       else
          sig = 1
       endif

       if(sig.eq.0) slope = missing
    else
       slope = missing
    endif
    deallocate(slopes)
    
    return
  end subroutine mann_kendall
  


end module LVT_TrendMod
