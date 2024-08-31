!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LVT_VarianceMod
!
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
! !MODULE: LVT_VarianceMod
! 
!  !DESCRIPTION: 
!   This module handles the computations required to compute variance values
!   of desired variables from the two datastreams. 
!
!   The metric calculations are supported across 
!   the entire analysis time period, stratified to a user defined
!   temporal averaging interval, average seasonal and diurnal cycle
!   and stratified to external datasets. 
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
!
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initVariance
  public :: LVT_diagnoseVariance
  public :: LVT_computeVariance
  public :: LVT_writeMetric_Variance
  public :: LVT_resetMetric_Variance
  public :: LVT_writerestart_Variance
  public :: LVT_readrestart_Variance
!EOP
  
  private

contains
  subroutine LVT_initVariance(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!EOP
    integer                 :: m

    allocate(stats%variance(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%variance(m)%model_value_total(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          allocate(stats%variance(m)%model_value_total_sxsx(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%variance(m)%model_value_total = 0.0
          stats%variance(m)%model_value_total_sxsx = 0.0
          allocate(stats%variance(m)%count_model_value_total(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%variance(m)%count_model_value_total = 0
          allocate(stats%variance(m)%model_value_ci(selectNlevs(1),LVT_rc%strat_nlevels))
          stats%variance(m)%model_value_ci = LVT_rc%udef

          if(metric%computeSC.eq.1) then 
             allocate(stats%variance(m)%model_value_asc(LVT_rc%ngrid, selectNlevs(1),&
                  LVT_rc%nasc))
             stats%variance(m)%model_value_asc = 0.0
             allocate(stats%variance(m)%model_value_asc_sxsx(LVT_rc%ngrid, selectNlevs(1),&
                  LVT_rc%nasc))
             stats%variance(m)%model_value_asc_sxsx = 0.0
             allocate(stats%variance(m)%count_model_value_asc(LVT_rc%ngrid, selectNlevs(1),&
                  LVT_rc%nasc))
             stats%variance(m)%count_model_value_asc = 0
          endif
          if(metric%computeADC.eq.1) then 
             allocate(stats%variance(m)%model_value_adc(LVT_rc%ngrid, selectNlevs(2),&
                  LVT_rc%nadc))
             stats%variance(m)%model_value_adc = 0.0
             allocate(stats%variance(m)%model_value_adc_sxsx(LVT_rc%ngrid, selectNlevs(2),&
                  LVT_rc%nadc))
             stats%variance(m)%model_value_adc_sxsx = 0.0
             allocate(stats%variance(m)%count_model_value_adc(LVT_rc%ngrid, selectNlevs(2),&
                  LVT_rc%nadc))
             stats%variance(m)%count_model_value_adc = 0
          endif

          if(LVT_rc%obssource(2).ne."none") then 
             allocate(stats%variance(m)%obs_value_total(LVT_rc%ngrid, selectNlevs(2), &
                  LVT_rc%strat_nlevels))
             stats%variance(m)%obs_value_total = 0.0
             allocate(stats%variance(m)%obs_value_total_sxsx(LVT_rc%ngrid, selectNlevs(2), &
                  LVT_rc%strat_nlevels))
             stats%variance(m)%obs_value_total_sxsx = 0.0
             allocate(stats%variance(m)%count_obs_value_total(LVT_rc%ngrid, selectNlevs(2), &
                  LVT_rc%strat_nlevels))
             stats%variance(m)%count_obs_value_total = 0
             allocate(stats%variance(m)%obs_value_ci(selectNlevs(2),LVT_rc%strat_nlevels))
             stats%variance(m)%obs_value_ci = LVT_rc%udef

             if(metric%computeSC.eq.1) then 
                allocate(stats%variance(m)%obs_value_asc(LVT_rc%ngrid, selectNlevs(2),&
                     LVT_rc%nasc))
                stats%variance(m)%obs_value_asc = 0.0
                allocate(stats%variance(m)%obs_value_asc_sxsx(LVT_rc%ngrid, selectNlevs(2),&
                     LVT_rc%nasc))
                stats%variance(m)%obs_value_asc_sxsx = 0.0
                allocate(stats%variance(m)%count_obs_value_asc(LVT_rc%ngrid, selectNlevs(2),&
                     LVT_rc%nasc))
                stats%variance(m)%count_obs_value_asc = 0
             endif
             if(metric%computeADC.eq.1) then 
                allocate(stats%variance(m)%obs_value_adc(LVT_rc%ngrid, selectNlevs(2),&
                     LVT_rc%nadc))
                stats%variance(m)%obs_value_adc = 0.0
                allocate(stats%variance(m)%obs_value_adc_sxsx(LVT_rc%ngrid, selectNlevs(2),&
                     LVT_rc%nadc))
                stats%variance(m)%obs_value_adc_sxsx = 0.0
                allocate(stats%variance(m)%count_obs_value_adc(LVT_rc%ngrid, selectNlevs(2),&
                     LVT_rc%nadc))
                stats%variance(m)%count_obs_value_adc = 0
             endif
          endif
       endif
       if(metric%timeOpt.eq.1) then 
          allocate(stats%variance(m)%model_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%variance(m)%model_value_ts = 0.0
          allocate(stats%variance(m)%model_value_ts_sxsx(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%variance(m)%model_value_ts_sxsx = 0.0
          allocate(stats%variance(m)%count_model_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%variance(m)%count_model_value_ts = 0

          allocate(stats%variance(m)%tavg_model_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%variance(m)%tavg_model_value_ts = 0.0
          allocate(stats%variance(m)%tavg_count_model_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%variance(m)%tavg_count_model_value_ts = 0

          if(LVT_rc%obssource(2).ne."none") then 
             allocate(stats%variance(m)%obs_value_ts(LVT_rc%ngrid, selectNlevs(1), &
                  LVT_rc%strat_nlevels))
             stats%variance(m)%obs_value_ts = 0.0
             allocate(stats%variance(m)%obs_value_ts_sxsx(LVT_rc%ngrid, selectNlevs(1), &
                  LVT_rc%strat_nlevels))
             stats%variance(m)%obs_value_ts_sxsx = 0.0
             allocate(stats%variance(m)%count_obs_value_ts(LVT_rc%ngrid, selectNlevs(1), &
                  LVT_rc%strat_nlevels))
             stats%variance(m)%count_obs_value_ts = 0

             allocate(stats%variance(m)%tavg_obs_value_ts(LVT_rc%ngrid, selectNlevs(1), &
                  LVT_rc%strat_nlevels))
             stats%variance(m)%tavg_obs_value_ts = 0.0
             allocate(stats%variance(m)%tavg_count_obs_value_ts(LVT_rc%ngrid, selectNlevs(1), &
                  LVT_rc%strat_nlevels))
             stats%variance(m)%tavg_count_obs_value_ts = 0
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
  end subroutine LVT_initVariance

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
! !ROUTINE: LVT_diagnoseVariance
! \label{LVT_diagnoseVariance}
! 
! !INTERFACE: 
  subroutine LVT_diagnoseVariance(pass)
! !USES:     

    implicit none

    integer       :: pass
! 
! !DESCRIPTION: 
!   This subroutine issues the calls to update the computations for 
!   calculating the std of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelVariance](\ref{diagnoseSingleModelVariance})
!     updates the std computation for a single variable 
!   \end{description}
! 
!EOP

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%variance%selectOpt.eq.1.or.&
            LVT_metrics%variance%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)
          
          do while(associated(model))

             call diagnoseSingleVariance(model,obs,stats,&
                  LVT_metrics%variance)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseVariance

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
! !ROUTINE: diagnoseSingleVariance
! \label{diagnoseSingleVariance}
! 
! !INTERFACE: 
  subroutine diagnoseSingleVariance(model, obs, stats,metric)
! !USES: 

    implicit none
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metadataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!   This routine updates the std computation of the 
!   specified variable. 
!
!  The arguments are: 
!
!  \begin{description}
!   \item[model] model variable object
!   \item[stats] object to hold the updated statistics
!  \end{description}
!EOP
    integer    :: t,k,m,m_k,o_k,tind

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             do m=1,LVT_rc%nensem
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(model%count(t,m,m_k).gt.0) then
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         stats%variance(m)%model_value_total(t,k,1) = &
                              stats%variance(m)%model_value_total(t,k,1) + &
                              model%value(t,m,m_k)
                         stats%variance(m)%model_value_total_sxsx(t,k,1) = &
                              stats%variance(m)%model_value_total_sxsx(t,k,1) + &
                              model%value(t,m,m_k)*model%value(t,m,m_k)
                         stats%variance(m)%count_model_value_total(t,k,1) = &
                              stats%variance(m)%count_model_value_total(t,k,1) + 1
                      endif
                   endif
                   
                   if(metric%timeOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         stats%variance(m)%model_value_ts(t,k,1) = &
                              stats%variance(m)%model_value_ts(t,k,1)+&
                              model%value(t,m,m_k)
                         stats%variance(m)%model_value_ts_sxsx(t,k,1) = &
                              stats%variance(m)%model_value_ts_sxsx(t,k,1)+&
                              model%value(t,m,m_k)*model%value(t,m,m_k)
                         stats%variance(m)%count_model_value_ts(t,k,1) = & 
                              stats%variance(m)%count_model_value_ts(t,k,1)+1
                      endif
                   endif
                   if(metric%computeSC.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                              tind)
                         stats%variance(m)%model_value_asc(t,k,tind) = &
                              stats%variance(m)%model_value_asc(t,k,tind)+&
                              model%value(t,m,m_k)
                         stats%variance(m)%model_value_asc_sxsx(t,k,tind) = &
                              stats%variance(m)%model_value_asc_sxsx(t,k,tind)+&
                              model%value(t,m,m_k)*model%value(t,m,m_k)
                         stats%variance(m)%count_model_value_asc(t,k,tind) = &
                              stats%variance(m)%count_model_value_asc(t,k,tind) + 1
                      endif
                   endif
                   if(metric%computeADC.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%variance(m)%model_value_adc(t,k,tind) = &
                              stats%variance(m)%model_value_adc(t,k,tind)+&
                              model%value(t,m,m_k)
                         stats%variance(m)%model_value_adc_sxsx(t,k,tind) = &
                              stats%variance(m)%model_value_adc_sxsx(t,k,tind)+&
                              model%value(t,m,m_k)*model%value(t,m,m_k)
                         stats%variance(m)%count_model_value_adc(t,k,tind) = &
                              stats%variance(m)%count_model_value_adc(t,k,tind) + 1
                      endif
                   endif
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,m,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         if(metric%selectOpt.eq.1) then
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then  
                               stats%variance(m)%model_value_total(t,k,2) = &
                                    stats%variance(m)%model_value_total(t,k,2) + &
                                    model%value(t,m,m_k)
                               stats%variance(m)%model_value_total_sxsx(t,k,2) = &
                                    stats%variance(m)%model_value_total_sxsx(t,k,2) + &
                                    model%value(t,m,m_k)*model%value(t,m,m_k)
                               stats%variance(m)%count_model_value_total(t,k,2) = &
                                    stats%variance(m)%count_model_value_total(t,k,2) + 1
                            endif
                         endif
                         if(metric%timeOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               stats%variance(m)%model_value_ts(t,k,2) = &
                                    stats%variance(m)%model_value_ts(t,k,2)+&
                                    model%value(t,m,m_k)
                               stats%variance(m)%model_value_ts_sxsx(t,k,2) = &
                                    stats%variance(m)%model_value_ts_sxsx(t,k,2)+&
                                    model%value(t,m,m_k)*model%value(t,m,m_k)
                               stats%variance(m)%count_model_value_ts(t,k,2) = & 
                                    stats%variance(m)%count_model_value_ts(t,k,2)+1
                            endif
                         endif
                      elseif(LVT_stats%strat_var(t,m,k).le.&
                           LVT_rc%strat_var_threshold) then
                         if(metric%selectOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               stats%variance(m)%model_value_total(t,k,3) = &
                                    stats%variance(m)%model_value_total(t,k,3) + &
                                    model%value(t,m,m_k)
                               stats%variance(m)%model_value_total_sxsx(t,k,3) = &
                                    stats%variance(m)%model_value_total_sxsx(t,k,3) + &
                                    model%value(t,m,m_k)*model%value(t,m,m_k)
                               stats%variance(m)%count_model_value_total(t,k,3) = &
                                    stats%variance(m)%count_model_value_total(t,k,3) + 1
                            endif
                         endif
                         if(metric%timeOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               stats%variance(m)%model_value_ts(t,k,3) = &
                                    stats%variance(m)%model_value_ts(t,k,3)+&
                                    model%value(t,m,m_k)
                               stats%variance(m)%model_value_ts_sxsx(t,k,3) = &
                                    stats%variance(m)%model_value_ts_sxsx(t,k,3)+&
                                    model%value(t,m,m_k)*model%value(t,m,m_k)
                               stats%variance(m)%count_model_value_ts(t,k,3) = & 
                                    stats%variance(m)%count_model_value_ts(t,k,3)+1
                            endif
                         endif
                      endif
                   endif
                endif
             end do
          enddo
          do k=1,obs%selectNlevs
             do m=1,LVT_rc%nensem
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      if(obs%count(t,m,o_k).gt.0) then 
                         if(metric%selectOpt.eq.1.and.&
                              obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            stats%variance(m)%obs_value_total(t,k,1) = &
                                 stats%variance(m)%obs_value_total(t,k,1) + &
                                 obs%value(t,m,o_k)
                            stats%variance(m)%obs_value_total_sxsx(t,k,1) = &
                                 stats%variance(m)%obs_value_total_sxsx(t,k,1) + &
                                 obs%value(t,m,o_k)*obs%value(t,m,o_k)
                            stats%variance(m)%count_obs_value_total(t,k,1) = &
                                 stats%variance(m)%count_obs_value_total(t,k,1) + 1
                         endif
                         
                         if(metric%timeOpt.eq.1.and.&
                              obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            stats%variance(m)%obs_value_ts(t,k,1) = &
                                 stats%variance(m)%obs_value_ts(t,k,1)+&
                                 obs%value(t,m,o_k)
                            stats%variance(m)%obs_value_ts_sxsx(t,k,1) = &
                                 stats%variance(m)%obs_value_ts_sxsx(t,k,1)+&
                                 obs%value(t,m,o_k)*obs%value(t,m,o_k)
                            stats%variance(m)%count_obs_value_ts(t,k,1) = & 
                                 stats%variance(m)%count_obs_value_ts(t,k,1)+1
                         endif
                         if(metric%computeSC.eq.1.and.&
                              obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,tind)
                            stats%variance(m)%obs_value_asc(t,k,tind) = &
                                 stats%variance(m)%obs_value_asc(t,k,tind)+&
                                 obs%value(t,m,o_k)
                            stats%variance(m)%obs_value_asc_sxsx(t,k,tind) = &
                                 stats%variance(m)%obs_value_asc_sxsx(t,k,tind)+&
                                 obs%value(t,m,o_k)*obs%value(t,m,o_k)
                            stats%variance(m)%count_obs_value_asc(t,k,tind) = &
                                 stats%variance(m)%count_obs_value_asc(t,k,tind) + 1
                         endif
                         if(metric%computeADC.eq.1.and.&
                              obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            call LVT_getADCTimeIndex(tind)
                            stats%variance(m)%obs_value_adc(t,k,tind) = &
                                 stats%variance(m)%obs_value_adc(t,k,tind)+&
                                 obs%value(t,m,o_k)
                            stats%variance(m)%obs_value_adc_sxsx(t,k,tind) = &
                                 stats%variance(m)%obs_value_adc_sxsx(t,k,tind)+&
                                 obs%value(t,m,o_k)*obs%value(t,m,o_k)
                            stats%variance(m)%count_obs_value_adc(t,k,tind) = &
                                 stats%variance(m)%count_obs_value_adc(t,k,tind) + 1
                         endif
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then
                               if(metric%selectOpt.eq.1 &
                                    .and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%variance(m)%obs_value_total(t,k,2) = &
                                       stats%variance(m)%obs_value_total(t,k,2) + &
                                       obs%value(t,m,o_k)
                                  stats%variance(m)%obs_value_total_sxsx(t,k,2) = &
                                       stats%variance(m)%obs_value_total_sxsx(t,k,2) + &
                                       obs%value(t,m,o_k)*obs%value(t,m,o_k)
                                  stats%variance(m)%count_obs_value_total(t,k,2) = &
                                       stats%variance(m)%count_obs_value_total(t,k,2) + 1
                               endif

                               if(metric%timeOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%variance(m)%obs_value_ts(t,k,2) = &
                                       stats%variance(m)%obs_value_ts(t,k,2)+&
                                       obs%value(t,m,o_k)
                                  stats%variance(m)%obs_value_ts_sxsx(t,k,2) = &
                                       stats%variance(m)%obs_value_ts_sxsx(t,k,2)+&
                                       obs%value(t,m,o_k)*obs%value(t,m,o_k)
                                  stats%variance(m)%count_obs_value_ts(t,k,2) = & 
                                       stats%variance(m)%count_obs_value_ts(t,k,2)+1
                               endif
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then
                               if(metric%selectOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%variance(m)%obs_value_total(t,k,3) = &
                                       stats%variance(m)%obs_value_total(t,k,3) + &
                                       obs%value(t,m,o_k)
                                  stats%variance(m)%obs_value_total_sxsx(t,k,3) = &
                                       stats%variance(m)%obs_value_total_sxsx(t,k,3) + &
                                       obs%value(t,m,o_k)*obs%value(t,m,o_k)
                                  stats%variance(m)%count_obs_value_total(t,k,3) = &
                                       stats%variance(m)%count_obs_value_total(t,k,3) + 1
                               endif

                               if(metric%timeOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%variance(m)%obs_value_ts(t,k,3) = &
                                       stats%variance(m)%obs_value_ts(t,k,3)+&
                                       obs%value(t,m,o_k)
                                  stats%variance(m)%obs_value_ts_sxsx(t,k,3) = &
                                       stats%variance(m)%obs_value_ts_sxsx(t,k,3)+&
                                       obs%value(t,m,o_k)*obs%value(t,m,o_k)
                                  stats%variance(m)%count_obs_value_ts(t,k,3) = & 
                                       stats%variance(m)%count_obs_value_ts(t,k,3)+1
                               endif

                            endif
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
  end subroutine diagnoseSingleVariance

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
! !ROUTINE: LVT_computeVariance
! \label{LVT_computeVariance}
! 
! !INTERFACE: 
  subroutine LVT_computeVariance(pass,alarm)
! !USES: 

    implicit none

    integer               :: pass
    logical               :: alarm

! 
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the std values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelVariance](\ref{computeSingleModelVariance})
!     computes the std values for a single variable
!   \end{description}
! 
!EOP
    integer     :: i,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%variance%selectOpt.eq.1.or.&
            LVT_metrics%variance%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%variance%timeOpt.eq.1.and.&
                  LVT_metrics%variance%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%variance%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%variance%ftn_ts_loc(i,1),200,advance='no') &
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
             call computeSingleVariance(alarm,model,obs,stats,&
                  LVT_metrics%variance)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
          
          if(alarm) then 
             if(LVT_metrics%variance%timeOpt.eq.1.and.&
                  LVT_metrics%variance%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%variance%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%variance%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif
             endif
          endif
       endif
    endif
  end subroutine LVT_computeVariance
  

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
! !ROUTINE: computeSingleVariance
! \label{computeSingleVariance}
! 
! !INTERFACE: 
  subroutine computeSingleVariance(alarm,model,obs,stats,metric)
! !USES: 

    implicit none
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
! 
! !DESCRIPTION: 
! 
! !DESCRIPTION: 
!  This routine computes the std values
!  The arguments are: 
!
!  \begin{description}
!    \item[model] model variable object
!    \item[stats] object to hold the updated statistics
!  \end{description}
!EOP

    integer  :: t,l,k,m

    real,    allocatable :: tavg_model_value_ts(:,:,:)
    real,    allocatable :: tavg_obs_value_ts(:,:,:)
    
    real,    allocatable :: model_value_asc(:,:,:)
    integer, allocatable :: count_model_value_asc(:,:,:)
    real,    allocatable :: obs_value_asc(:,:,:)

    real,    allocatable :: model_value_adc(:,:,:)
    real,    allocatable :: obs_value_adc(:,:,:)
    
    real,    allocatable :: model_value_avg(:,:,:)
    real,    allocatable :: obs_value_avg(:,:,:)   

    real     :: term1, term2

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%variance(m)%count_model_value_ts(t,k,l).gt.0) then 
                         term1 = (stats%variance(m)%model_value_ts_sxsx(t,k,l)/&
                              stats%variance(m)%count_model_value_ts(t,k,l))
                         term2 =(stats%variance(m)%model_value_ts(t,k,l)/&
                              stats%variance(m)%count_model_value_ts(t,k,l))**2 
                         
                         if(term1.gt.term2) then 
                            stats%variance(m)%model_value_ts(t,k,l) = &
                                 (term1 -term2)
                            
                            stats%variance(m)%tavg_model_value_ts(t,k,l) = & 
                                 stats%variance(m)%tavg_model_value_ts(t,k,l) + & 
                                 stats%variance(m)%model_value_ts(t,k,l)
                            stats%variance(m)%tavg_count_model_value_ts(t,k,l) = & 
                                 stats%variance(m)%tavg_count_model_value_ts(t,k,l) + 1
                         else
                            stats%variance(m)%model_value_ts(t,k,l) = 0.0
                         endif
                         
                      else
                         stats%variance(m)%model_value_ts(t,k,l) = LVT_rc%udef
                      endif

                      if(LVT_rc%obssource(2).ne."none") then 
                         if(stats%variance(m)%count_obs_value_ts(t,k,l).gt.0) then 
                            term1 = (stats%variance(m)%obs_value_ts_sxsx(t,k,l)/&
                                 stats%variance(m)%count_obs_value_ts(t,k,l))
                            term2 = (stats%variance(m)%obs_value_ts(t,k,l)/&
                                 stats%variance(m)%count_obs_value_ts(t,k,l))**2
                            if(term1.gt.term2) then 
                               stats%variance(m)%obs_value_ts(t,k,l) = & 
                                    (term1 - term2)
                               stats%variance(m)%tavg_obs_value_ts(t,k,l) = & 
                                    stats%variance(m)%tavg_obs_value_ts(t,k,l) + & 
                                    stats%variance(m)%obs_value_ts(t,k,l)
                               stats%variance(m)%tavg_count_obs_value_ts(t,k,l) = & 
                                    stats%variance(m)%tavg_count_obs_value_ts(t,k,l) + 1
                            else
                               stats%variance(m)%obs_value_ts(t,k,l) = 0.0
                            endif
                         else
                            stats%variance(m)%obs_value_ts(t,k,l) = LVT_rc%udef
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo

          if(alarm) then 

             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         if(stats%variance(m)%tavg_count_model_value_ts(t,k,l).gt.0) then 
                            stats%variance(m)%tavg_model_value_ts(t,k,l) = & 
                                 stats%variance(m)%tavg_model_value_ts(t,k,l) / & 
                                 stats%variance(m)%tavg_count_model_value_ts(t,k,l) 
                         else
                            stats%variance(m)%tavg_model_value_ts(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   enddo
                   if(LVT_rc%obssource(2).ne."none") then 
                      do k=1,model%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            if(stats%variance(m)%tavg_count_obs_value_ts(t,k,l).gt.0) then 
                               stats%variance(m)%tavg_obs_value_ts(t,k,l) = & 
                                    stats%variance(m)%tavg_obs_value_ts(t,k,l) / & 
                                    stats%variance(m)%tavg_count_obs_value_ts(t,k,l) 
                            else
                               stats%variance(m)%tavg_obs_value_ts(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      enddo
                      
                   endif
                enddo
             enddo

             if(metric%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then
                         call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                              model,&
                              LVT_rc%ngrid,&
                              stats%variance(m)%tavg_model_value_ts,&
                              stats%variance(m)%tavg_count_model_value_ts,&
                              LVT_rc%ngrid,&
                              stats%variance(m)%tavg_obs_value_ts,&
                              stats%variance(m)%tavg_count_obs_value_ts)
                      else
                         call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                              model,&
                              LVT_rc%ngrid,&
                              stats%variance(m)%tavg_model_value_ts,&
                              stats%variance(m)%tavg_count_model_value_ts)
                      endif
                   enddo
                else
                   allocate(tavg_model_value_ts(LVT_rc%ngrid, &
                        model%vlevels, &
                        LVT_rc%strat_nlevels))
                   tavg_model_value_ts = 0.0

                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 

                      allocate(tavg_obs_value_ts(LVT_rc%ngrid, &
                           obs%vlevels, &
                           LVT_rc%strat_nlevels))
                      tavg_obs_value_ts = 0.0

                         
                      do m=1,LVT_rc%nensem
                         do t=1,LVT_rc%ngrid
                            do k=1,model%selectNlevs
                               do l=1,LVT_rc%strat_nlevels
                                  tavg_model_value_ts(t,k,l) = &
                                       tavg_model_value_ts(t,k,l) + &
                                       stats%variance(m)%tavg_model_value_ts(t,k,l)
                                  tavg_obs_value_ts(t,k,l) = &
                                       tavg_obs_value_ts(t,k,l) + &
                                       stats%variance(m)%tavg_obs_value_ts(t,k,l)
                               enddo
                            enddo
                         enddo
                      enddo
                            
                      do t=1,LVT_rc%ngrid
                         do k=1,model%selectNlevs
                            do l=1,LVT_rc%strat_nlevels
                               tavg_model_value_ts(t,k,l) = &
                                    tavg_model_value_ts(t,k,l)/LVT_rc%nensem
                               tavg_obs_value_ts(t,k,l) = &
                                    tavg_obs_value_ts(t,k,l)/LVT_rc%nensem
                            enddo
                         enddo
                      enddo
                         
                      call LVT_writeTSinfo(metric%ftn_ts_loc(:,1),&
                           model,&
                           LVT_rc%ngrid,&
                           tavg_model_value_ts,&
                           stats%variance(1)%tavg_count_model_value_ts,&
                           LVT_rc%ngrid,&
                           tavg_obs_value_ts,&
                           stats%variance(1)%tavg_count_obs_value_ts)
                   else
                      do m=1,LVT_rc%nensem
                         do t=1,LVT_rc%ngrid
                            do k=1,model%selectNlevs
                               do l=1,LVT_rc%strat_nlevels
                                  tavg_model_value_ts(t,k,l) = &
                                       tavg_model_value_ts(t,k,l) + &
                                       stats%variance(m)%tavg_model_value_ts(t,k,l)
                               enddo
                            enddo
                         enddo
                      enddo
                      
                      do t=1,LVT_rc%ngrid
                         do k=1,model%selectNlevs
                            do l=1,LVT_rc%strat_nlevels
                               tavg_model_value_ts(t,k,l) = &
                                    tavg_model_value_ts(t,k,l)/LVT_rc%nensem
                            enddo
                         enddo
                      enddo
                      
                      call LVT_writeTSinfo(metric%ftn_ts_loc(:,1),&
                           model,&
                           LVT_rc%ngrid,&
                           tavg_model_value_ts,&
                           stats%variance(1)%tavg_count_model_value_ts)
                   endif
                   deallocate(tavg_model_value_ts)
                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                      deallocate(tavg_obs_value_ts)
                   endif
                endif
             endif
          endif
       endif
    endif

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%variance(m)%count_model_value_total(t,k,l).gt.&
                           LVT_rc%obsCountThreshold) then 
                         term1 = (stats%variance(m)%model_value_total_sxsx(t,k,l)/&
                              stats%variance(m)%count_model_value_total(t,k,l)) - & 
                              (stats%variance(m)%model_value_total(t,k,l)/&
                              stats%variance(m)%count_model_value_total(t,k,l))**2
                         if(term1.gt.0) then 
                            stats%variance(m)%model_value_total(t,k,l) = &
                                 (term1)
                         else
                            stats%variance(m)%model_value_total(t,k,l) = LVT_rc%udef
                         endif
                      else
                         stats%variance(m)%model_value_total(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                   if(metric%computeSC.eq.1) then 
                      do l=1,LVT_rc%nasc
                         if(stats%variance(m)%count_model_value_asc(t,k,l).gt.&
                              LVT_rc%SCCountThreshold) then 
                            stats%variance(m)%model_value_asc(t,k,l) = &
                                 ((stats%variance(m)%model_value_asc_sxsx(t,k,l)/&
                                 stats%variance(m)%count_model_value_asc(t,k,l)) - & 
                                 (stats%variance(m)%model_value_asc(t,k,l)/&
                                 stats%variance(m)%count_model_value_asc(t,k,l))**2)
                         else
                            stats%variance(m)%model_value_asc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   if(metric%computeADC.eq.1) then 
                      do l=1,LVT_rc%nadc
                         if(stats%variance(m)%count_model_value_adc(t,k,l).gt.&
                              LVT_rc%ADCCountThreshold) then 
                            stats%variance(m)%model_value_adc(t,k,l) = &
                                 ((stats%variance(m)%model_value_adc_sxsx(t,k,l)/&
                                 stats%variance(m)%count_model_value_adc(t,k,l)) - & 
                                 (stats%variance(m)%model_value_adc(t,k,l)/&
                                 stats%variance(m)%count_model_value_adc(t,k,l))**2)
                         else
                            stats%variance(m)%model_value_adc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   
                   if(LVT_rc%obssource(2).ne."none") then 
                      do l=1,LVT_rc%strat_nlevels
                         if(stats%variance(m)%count_obs_value_total(t,k,l).gt.&
                              LVT_rc%obsCountThreshold) then 
                            term1 = (stats%variance(m)%obs_value_total_sxsx(t,k,l)/&
                                 stats%variance(m)%count_obs_value_total(t,k,l)) - & 
                                 (stats%variance(m)%obs_value_total(t,k,l)/&
                                 stats%variance(m)%count_obs_value_total(t,k,l))**2
                            if(term1.gt.0) then 
                               stats%variance(m)%obs_value_total(t,k,l) = &
                                    (term1)
                            else
                               stats%variance(m)%obs_value_total(t,k,l)= LVT_rc%udef
                            endif
                         else
                            stats%variance(m)%obs_value_total(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                      
                      if(metric%computeSC.eq.1) then
                         do l=1,LVT_rc%nasc
                            if(stats%variance(m)%count_obs_value_asc(t,k,l).gt.&
                                 LVT_rc%SCCountThreshold) then 
                               
                               term1 = stats%variance(m)%obs_value_asc_sxsx(t,k,l)/&
                                    stats%variance(m)%count_obs_value_asc(t,k,l)
                               term2 = (stats%variance(m)%obs_value_asc(t,k,l)/&
                                    stats%variance(m)%count_obs_value_asc(t,k,l))**2
                               if(term1.gt.term2) then 
                                  stats%variance(m)%obs_value_asc(t,k,l) = &
                                       (term1-term2)
                               elseif(abs(term1-term2).lt.0.00001) then 
                                  stats%variance(m)%obs_value_asc(t,k,l) = 0.00
                               else
                                  print*, 'Error in Variance calculation..'
                                  stop
                               endif
                            else
                               stats%variance(m)%obs_value_asc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                      
                      if(metric%computeADC.eq.1) then 
                         do l=1,LVT_rc%nadc
                            if(stats%variance(m)%count_obs_value_adc(t,k,l).gt.&
                                 LVT_rc%ADCCountThreshold) then 
                               stats%variance(m)%obs_value_adc(t,k,l) = &
                                    ((stats%variance(m)%obs_value_adc_sxsx(t,k,l)/&
                                    stats%variance(m)%count_obs_value_adc(t,k,l)) - & 
                                    (stats%variance(m)%obs_value_adc(t,k,l)/&
                                    stats%variance(m)%count_obs_value_adc(t,k,l))**2)
                            else
                               stats%variance(m)%obs_value_adc(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                   endif
                enddo
             enddo
          enddo

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
                              stats%variance(m)%model_value_total(t,k,l)
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
                                 stats%variance(m)%obs_value_total(t,k,l)
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

          if(metric%computeSC.eq.1) then
             allocate(model_value_asc(LVT_rc%ngrid,&
                  model%vlevels,&
                  LVT_rc%nasc))
             allocate(count_model_value_asc(LVT_rc%ngrid,&
                  model%vlevels,&
                  LVT_rc%nasc))
             model_value_asc = 0.0
             count_model_value_asc = 0 

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%nasc
                      do m=1,LVT_rc%nensem        
                         if(stats%variance(m)%model_value_asc(t,k,l).ne.LVT_rc%udef) then 
                            model_value_asc(t,k,l) = &
                                 model_value_asc(t,k,l) + & 
                                 stats%variance(m)%model_value_asc(t,k,l)
                            count_model_value_asc(t,k,l) = & 
                                 count_model_value_asc(t,k,l) + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%nasc
                      if(count_model_value_asc(t,k,l).gt.0) then 
                         model_value_asc(t,k,l) = &
                              model_value_asc(t,k,l)/&
                              count_model_value_asc(t,k,l)
                      else
                         model_value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                allocate(obs_value_asc(LVT_rc%ngrid,&
                     obs%vlevels,&
                     LVT_rc%nasc))
                
                obs_value_asc = 0.0

                do t=1,LVT_rc%ngrid
                   do k=1,obs%selectNlevs                
                      do l=1, LVT_rc%nasc
                         do m=1,LVT_rc%nensem          
                            obs_value_asc(t,k,l) = &
                                 obs_value_asc(t,k,l) + & 
                                 stats%variance(m)%obs_value_asc(t,k,l)
                         enddo
                      enddo
                   enddo
                enddo
                
                do t=1,LVT_rc%ngrid
                   do k=1,obs%selectNlevs                
                      do l=1, LVT_rc%nasc
                         obs_value_asc(t,k,l) = &
                              obs_value_asc(t,k,l)/LVT_rc%nensem
                      enddo
                   enddo
                enddo
             
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_asc,&
                     stats%variance(1)%count_model_value_asc,&
                     LVT_rc%ngrid,obs_value_asc,&
                     stats%variance(1)%count_obs_value_asc)   

                deallocate(obs_value_asc)
       
             else
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_asc,&
                     stats%variance(1)%count_model_value_asc)

                deallocate(model_value_asc)
             endif
          endif

          if(metric%computeADC.eq.1) then
             allocate(model_value_adc(LVT_rc%ngrid,&
                  model%vlevels,&
                  LVT_rc%nadc))
             model_value_adc = 0.0

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%nadc
                      do m=1,LVT_rc%nensem          
                         model_value_adc(t,k,l) = &
                              model_value_adc(t,k,l) + & 
                              stats%variance(m)%model_value_adc(t,k,l)
                      enddo
                   enddo
                enddo
             enddo
             
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%nadc
                      model_value_adc(t,k,l) = &
                           model_value_adc(t,k,l)/LVT_rc%nensem
                   enddo
                enddo
             enddo
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                allocate(obs_value_adc(LVT_rc%ngrid,&
                     obs%vlevels,&
                     LVT_rc%nadc))
                
                obs_value_adc = 0.0

                do t=1,LVT_rc%ngrid
                   do k=1,obs%selectNlevs                
                      do l=1, LVT_rc%nadc
                         do m=1,LVT_rc%nensem          
                            obs_value_adc(t,k,l) = &
                                 obs_value_adc(t,k,l) + & 
                                 stats%variance(m)%obs_value_adc(t,k,l)
                         enddo
                      enddo
                   enddo
                enddo
                
                do t=1,LVT_rc%ngrid
                   do k=1,obs%selectNlevs                
                      do l=1, LVT_rc%nadc
                         obs_value_adc(t,k,l) = &
                              obs_value_adc(t,k,l)/LVT_rc%nensem
                      enddo
                   enddo
                enddo
             
                call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_adc,&
                     stats%variance(1)%count_model_value_adc,&
                     LVT_rc%ngrid,obs_value_adc,&
                     stats%variance(1)%count_obs_value_adc)      
                
                deallocate(obs_value_adc)
             else
                call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_adc,&
                     stats%variance(m)%count_model_value_adc)

                deallocate(model_value_adc)
             endif
          endif
       endif
    endif
  end subroutine computeSingleVariance


  subroutine LVT_writeMetric_Variance(pass,final,vlevels,stats,obs)
! !USES:
    use LVT_statsMod, only : LVT_writeSummaryStats
    use LVT_pluginIndices
! !ARGUMENTS: 

    implicit none

    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)   :: stats
    type(LVT_metaDataEntry) :: obs

    integer                 :: k,l,m,tind

    real,    allocatable    :: model_value_total(:,:)
    integer, allocatable    :: count_model_value_total(:,:)
    real,    allocatable    :: model_value_ci(:)

    real,    allocatable    :: obs_value_total(:,:)
    integer, allocatable    :: count_obs_value_total(:,:)
    real,    allocatable    :: obs_value_ci(:)

    real,    allocatable    :: model_value_ts(:,:)
    integer, allocatable    :: count_model_value_ts(:,:)

    real,    allocatable    :: obs_value_ts(:,:)
    integer, allocatable    :: count_obs_value_ts(:,:)

    real,    allocatable    :: model_value_adc(:,:,:)
    real,    allocatable    :: model_value_asc(:,:,:)

    real,    allocatable    :: obs_value_adc(:,:,:)
    real,    allocatable    :: obs_value_asc(:,:,:)

    if(pass.eq.LVT_metrics%variance%npass) then
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              
              allocate(model_value_ts(LVT_rc%ngrid,LVT_rc%nensem))
              allocate(count_model_value_ts(LVT_rc%ngrid,LVT_rc%nensem))

              if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                 allocate(obs_value_ts(LVT_rc%ngrid,LVT_rc%nensem))
                 allocate(count_obs_value_ts(LVT_rc%ngrid,LVT_rc%nensem))
              endif

              do k=1,vlevels
                 do l=1,LVT_rc%strat_nlevels
                    do m=1,LVT_rc%nensem
                         model_value_ts(:,m) = &
                              stats%variance(m)%tavg_model_value_ts(:,k,l)
                         count_model_value_ts(:,m) = & 
                              stats%variance(m)%tavg_count_model_value_ts(:,k,l)
                       if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                          obs_value_ts(:,m) = &
                               stats%variance(m)%tavg_obs_value_ts(:,k,l)
                          count_obs_value_ts(:,m) = & 
                               stats%variance(m)%tavg_count_obs_value_ts(:,k,l)
                       endif
                    enddo
                    if(LVT_metrics%variance%timeOpt.eq.1) then 
                       
                       call LVT_writevar_gridded(LVT_metrics%variance%ftn_ts, &
                            model_value_ts(:,:),&
                            stats%vid_ts(LVT_VARIANCEid,1),k)
                       
                       call LVT_writevar_gridded(LVT_metrics%variance%ftn_ts, &
                            real(count_model_value_ts(:,:)),&
                            stats%vid_count_ts(LVT_VARIANCEid,1),k)
                    endif
                    
                    if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                       call LVT_writevar_gridded(LVT_metrics%variance%ftn_ts, &
                            obs_value_ts(:,:),&
                            stats%vid_ts(LVT_VARIANCEid,2),k)
                       call LVT_writevar_gridded(LVT_metrics%variance%ftn_ts, &
                            real(count_obs_value_ts(:,:)),&
                            stats%vid_count_ts(LVT_VARIANCEid,2),k)
                    endif
                 enddo
              enddo

              deallocate(model_value_ts)
              deallocate(count_model_value_ts)

              if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                 deallocate(obs_value_ts)
                 deallocate(count_obs_value_ts)
              endif
           endif
       else
          if(pass.eq.LVT_metrics%variance%npass) then
             if(stats%selectOpt.eq.1) then

              allocate(model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
              allocate(count_model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
              allocate(model_value_ci(LVT_rc%nensem))

              if(LVT_metrics%variance%computeSC.eq.1) then 
                 allocate(model_value_asc(LVT_rc%ngrid,&
                      LVT_rc%nensem,&
                      LVT_rc%nasc))
                 if(LVT_rc%obssource(2).ne."none".and.&
                      obs%selectNlevs.ge.1) then 
                    allocate(obs_value_asc(LVT_rc%ngrid,&
                         LVT_rc%nensem,&
                         LVT_rc%nasc))
                 endif
              endif
              if(LVT_metrics%variance%computeADC.eq.1) then 
                 allocate(model_value_adc(LVT_rc%ngrid,&
                      LVT_rc%nensem,&
                      LVT_rc%nadc))        
                 if(LVT_rc%obssource(2).ne."none".and.&
                      obs%selectNlevs.ge.1) then 
                    allocate(obs_value_adc(LVT_rc%ngrid,&
                         LVT_rc%nensem,&
                         LVT_rc%nadc))
                 endif

              endif

              if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                 allocate(obs_value_total(LVT_rc%ngrid,LVT_rc%nensem))
                 allocate(count_obs_value_total(LVT_rc%ngrid,LVT_rc%nensem))
                 allocate(obs_value_ci(LVT_rc%nensem))
              endif                

                do k=1,vlevels
                   do l=1,LVT_rc%strat_nlevels
                      do m=1,LVT_rc%nensem
                         model_value_total(:,m) = &
                              stats%variance(m)%model_value_total(:,k,l)
                         count_model_value_total(:,m) = & 
                              stats%variance(m)%count_model_value_total(:,k,l)
                         model_value_ci(m) = stats%variance(m)%model_value_ci(k,l)
                         
                         if(LVT_metrics%variance%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               model_value_asc(:,m,tind) = & 
                                    stats%variance(m)%model_value_asc(:,k,tind)
                            enddo
                         endif

                         if(LVT_metrics%variance%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               model_value_adc(:,m,tind) = & 
                                    stats%variance(m)%model_value_adc(:,k,tind)
                            enddo
                         endif

                         if(LVT_rc%obssource(2).ne."none".and.&
                              obs%selectNlevs.ge.1) then 
                            obs_value_total(:,m) = &
                                 stats%variance(m)%obs_value_total(:,k,l)
                            count_obs_value_total(:,m) = & 
                                 stats%variance(m)%count_obs_value_total(:,k,l)
                            obs_value_ci(m) = stats%variance(m)%obs_value_ci(k,l)

                            if(LVT_metrics%variance%computeSC.eq.1) then 
                               do tind = 1,LVT_rc%nasc
                                  obs_value_asc(:,m,tind) = & 
                                       stats%variance(m)%obs_value_asc(:,k,tind)
                               enddo
                            endif
                            if(LVT_metrics%variance%computeADC.eq.1) then 
                               do tind = 1,LVT_rc%nadc
                                  obs_value_adc(:,m,tind) = & 
                                       stats%variance(m)%obs_value_adc(:,k,tind)
                               enddo
                            endif
                            
                         endif
                      enddo
                      if(LVT_metrics%variance%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%variance%ftn_total, &
                              model_value_total(:,:),&
                              stats%vid_total(LVT_VARIANCEid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%variance%ftn_total, &
                              real(count_model_value_total(:,:)),&
                              stats%vid_count_total(LVT_VARIANCEid,1),k)

                         if(LVT_rc%obssource(2).ne."none".and.&
                              obs%selectNlevs.ge.1) then 
                            call LVT_writevar_gridded(LVT_metrics%variance%ftn_total, &
                                 obs_value_total(:,:),&
                                 stats%vid_total(LVT_VARIANCEid,2),k)
                            call LVT_writevar_gridded(LVT_metrics%variance%ftn_total, &
                                 real(count_obs_value_total(:,:)),&
                                 stats%vid_count_total(LVT_VARIANCEid,2),k)
                         endif
                         
                         if(LVT_metrics%variance%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%variance%ftn_total,&
                                    model_value_asc(:,:,tind),&
                                    stats%vid_sc_total(tind,LVT_VARIANCEid,1),k)
                            enddo
                            if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                               do tind = 1,LVT_rc%nasc
                                  call LVT_writevar_gridded(&
                                       LVT_metrics%variance%ftn_total,&
                                       obs_value_asc(:,:,tind),&
                                       stats%vid_sc_total(tind,LVT_VARIANCEid,2),k)
                               enddo
                            endif
                         endif
                         if(LVT_metrics%variance%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%variance%ftn_total,&
                                    model_value_adc(:,:,tind),&
                                    stats%vid_adc_total(tind,LVT_VARIANCEid,1),k)
                            enddo
                            if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                               do tind = 1,LVT_rc%nadc
                                  call LVT_writevar_gridded(&
                                       LVT_metrics%variance%ftn_total,&
                                       obs_value_adc(:,:,tind),&
                                       stats%vid_adc_total(tind,LVT_VARIANCEid,2),k)
                               enddo
                            endif
                         endif
                         call LVT_writeSummaryStats(&
                              LVT_metrics%variance%ftn_summ,&
                              l,&
                              LVT_metrics%variance%short_name,&
                              LVT_rc%ngrid,&
                              model_value_total(:,:), &
                              count_model_value_total(:,:),&
                              stats%standard_name,&
                              model_value_ci(:))
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            call LVT_writeSummaryStats(&
                                 LVT_metrics%variance%ftn_summ,&
                                 l,&
                                 LVT_metrics%variance%short_name,&
                                 LVT_rc%ngrid,&
                                 obs_value_total(:,:), &
                                 count_obs_value_total(:,:),&
                                 "DS2_"//trim(stats%standard_name),&
                                 obs_value_ci(:))
                         endif
                      endif
                   enddo
                enddo

                deallocate(model_value_total)
                deallocate(count_model_value_total)
                deallocate(model_value_ci)

                if(LVT_metrics%variance%computeSC.eq.1) then 
                   deallocate(model_value_asc)
                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                      deallocate(obs_value_asc)
                   endif
                endif
                if(LVT_metrics%variance%computeADC.eq.1) then 
                   deallocate(model_value_adc)
                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                      deallocate(obs_value_adc)
                   endif
                endif

                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   deallocate(obs_value_total)
                   deallocate(count_obs_value_total)
                   deallocate(obs_value_ci)
                endif

             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_Variance

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
  subroutine LVT_resetMetric_Variance(alarm)

    logical                 :: alarm
    integer                 :: vlevels
    integer                :: i,k,l,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats


    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                if(LVT_metrics%variance%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%variance(m)%model_value_ts(:,k,l) = 0.0
                      stats%variance(m)%model_value_ts_sxsx(:,k,l) = 0.0
                      stats%variance(m)%count_model_value_ts(:,k,l)=0 
                      if(alarm) then 
                         stats%variance(m)%tavg_model_value_ts(:,k,l) = 0.0
                         stats%variance(m)%tavg_count_model_value_ts(:,k,l)=0 
                      endif
                   enddo
                   if(LVT_rc%obssource(2).ne."none") then 
                      do l=1,LVT_rc%strat_nlevels
                         stats%variance(m)%obs_value_ts(:,k,l) = 0.0
                         stats%variance(m)%obs_value_ts_sxsx(:,k,l) = 0.0
                         stats%variance(m)%count_obs_value_ts(:,k,l)=0 
                         if(alarm) then 
                            stats%variance(m)%tavg_obs_value_ts(:,k,l) = 0.0
                            stats%variance(m)%tavg_count_obs_value_ts(:,k,l)=0 
                         endif
                      enddo
                      
                   endif
                endif
             enddo
          enddo
       endif

       model => model%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_Variance


!BOP
! 
! !ROUTINE: LVT_writerestart_Variance
! 
! !INTERFACE:
  subroutine LVT_writerestart_Variance(ftn, pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Variance metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,l,m
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%variance%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_writevar_restart(ftn,&
                           stats%variance(m)%model_value_total(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%variance(m)%model_value_total_sxsx(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%variance(m)%count_model_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%variance%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_writevar_restart(ftn, &
                              stats%variance(m)%model_value_asc(:,k,l))
                         call LVT_writevar_restart(ftn, &
                              stats%variance(m)%model_value_asc_sxsx(:,k,l))
                         call LVT_writevar_restart(ftn, &
                              stats%variance(m)%count_model_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%variance%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_writevar_restart(ftn,&
                              stats%variance(m)%model_value_adc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%variance(m)%model_value_adc_sxsx(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%variance(m)%count_model_value_adc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels                
                            call LVT_writevar_restart(ftn,&
                                 stats%variance(m)%obs_value_total(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%variance(m)%obs_value_total_sxsx(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%variance(m)%count_obs_value_total(:,k,l))
                         enddo
                      enddo
                      
                      if(LVT_metrics%variance%computeSC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_writevar_restart(ftn,&
                                    stats%variance(m)%obs_value_asc(:,k,l))
                               call LVT_writevar_restart(ftn,&
                                    stats%variance(m)%obs_value_asc_sxsx(:,k,l))
                               call LVT_writevar_restart(ftn,&
                                    stats%variance(m)%count_obs_value_asc(:,k,l))
                            enddo
                         enddo
                      endif
                      if(LVT_metrics%variance%computeADC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nadc
                               call LVT_writevar_restart(ftn,&
                                    stats%variance(m)%obs_value_adc(:,k,l))
                               call LVT_writevar_restart(ftn, &
                                    stats%variance(m)%obs_value_adc_sxsx(:,k,l))
                               call LVT_writevar_restart(ftn,&
                                    stats%variance(m)%count_obs_value_adc(:,k,l))
                            enddo
                         enddo
                      endif
                   endif
                endif
             enddo
          endif
       endif
       
       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo
  end subroutine LVT_writerestart_Variance


!BOP
! 
! !ROUTINE: LVT_readrestart_Variance
! 
! !INTERFACE:
  subroutine LVT_readrestart_Variance(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Variance metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
   
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats
    integer              :: k,l,m


    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%variance%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      call LVT_readvar_restart(ftn,&
                           stats%variance(m)%model_value_total(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%variance(m)%model_value_total_sxsx(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%variance(m)%count_model_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%variance%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_readvar_restart(ftn, &
                              stats%variance(m)%model_value_asc(:,k,l))
                         call LVT_readvar_restart(ftn, &
                              stats%variance(m)%model_value_asc_sxsx(:,k,l))
                         call LVT_readvar_restart(ftn, &
                              stats%variance(m)%count_model_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%variance%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_readvar_restart(ftn,&
                              stats%variance(m)%model_value_adc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%variance(m)%model_value_adc_sxsx(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%variance(m)%count_model_value_adc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels                
                            call LVT_readvar_restart(ftn,&
                                 stats%variance(m)%obs_value_total(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%variance(m)%obs_value_total_sxsx(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%variance(m)%count_obs_value_total(:,k,l))
                         enddo
                      enddo

                      if(LVT_metrics%variance%computeSC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_readvar_restart(ftn,&
                                    stats%variance(m)%obs_value_asc(:,k,l))
                               call LVT_readvar_restart(ftn,&
                                    stats%variance(m)%obs_value_asc_sxsx(:,k,l))
                               call LVT_readvar_restart(ftn,&
                                    stats%variance(m)%count_obs_value_asc(:,k,l))
                            enddo
                         enddo
                      endif
                      if(LVT_metrics%variance%computeADC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nadc
                               call LVT_readvar_restart(ftn,&
                                    stats%variance(m)%obs_value_adc(:,k,l))
                               call LVT_readvar_restart(ftn, &
                                    stats%variance(m)%obs_value_adc_sxsx(:,k,l))
                               call LVT_readvar_restart(ftn,&
                                    stats%variance(m)%count_obs_value_adc(:,k,l))
                            enddo
                         enddo
                      endif
                   endif
                end if
             enddo
          endif
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo
  end subroutine LVT_readrestart_Variance
end module LVT_VarianceMod
