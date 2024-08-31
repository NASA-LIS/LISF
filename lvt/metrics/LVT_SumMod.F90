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
! !MODULE: LVT_SumMod
! \label(LVT_SumMod)
!
! !INTERFACE:
module LVT_SumMod
! 
! !USES:   
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the computations required to compute sum values
!  of desired variables from the two datastreams. 
! 
!   The metric calculations are supported across 
!   the entire analysis time period, stratified to a user defined
!   temporal averaging interval, average seasonal and diurnal cycle
!   and stratified to external datasets. 
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initSum
  public :: LVT_diagnoseSum
  public :: LVT_computeSum
  public :: LVT_writeMetric_Sum
  public :: LVT_resetMetric_Sum
  public :: LVT_writerestart_Sum
  public :: LVT_readrestart_Sum
!EOP
  
  private

contains
  subroutine LVT_initSum(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!EOP
    integer                  :: m
    
    allocate(stats%sum(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%sum(m)%model_value_total(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%sum(m)%model_value_total = 0.0
          allocate(stats%sum(m)%count_model_value_total(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%sum(m)%count_model_value_total = 0
          allocate(stats%sum(m)%model_value_ci(selectNlevs(1),LVT_rc%strat_nlevels))
          stats%sum(m)%model_value_ci = LVT_rc%udef
          
          if(metric%computeSC.eq.1) then 
             allocate(stats%sum(m)%model_value_asc(LVT_rc%ngrid, selectNlevs(1),&
                  LVT_rc%nasc))
             stats%sum(m)%model_value_asc = 0.0
             allocate(stats%sum(m)%count_model_value_asc(LVT_rc%ngrid, selectNlevs(1),&
                  LVT_rc%nasc))
             stats%sum(m)%count_model_value_asc = 0
          endif
          if(metric%computeADC.eq.1) then 
             allocate(stats%sum(m)%model_value_adc(LVT_rc%ngrid, selectNlevs(2),&
                  LVT_rc%nadc))
             stats%sum(m)%model_value_adc = 0.0
             allocate(stats%sum(m)%count_model_value_adc(LVT_rc%ngrid, selectNlevs(2),&
                  LVT_rc%nadc))
             stats%sum(m)%count_model_value_adc = 0
          endif
          
          if(LVT_rc%obssource(2).ne."none") then 
             if(selectNlevs(2).ge.1) then 
                allocate(stats%sum(m)%obs_value_total(LVT_rc%ngrid, selectNlevs(2), &
                     LVT_rc%strat_nlevels))
                stats%sum(m)%obs_value_total = 0.0
                allocate(stats%sum(m)%count_obs_value_total(LVT_rc%ngrid, selectNlevs(2), &
                     LVT_rc%strat_nlevels))
                stats%sum(m)%count_obs_value_total = 0
                allocate(stats%sum(m)%obs_value_ci(selectNlevs(2),LVT_rc%strat_nlevels))
                stats%sum(m)%obs_value_ci = LVT_rc%udef
                
                if(metric%computeSC.eq.1) then 
                   allocate(stats%sum(m)%obs_value_asc(LVT_rc%ngrid, selectNlevs(2),&
                        LVT_rc%nasc))
                   stats%sum(m)%obs_value_asc = 0.0
                   allocate(stats%sum(m)%count_obs_value_asc(LVT_rc%ngrid, selectNlevs(2),&
                        LVT_rc%nasc))
                   stats%sum(m)%count_obs_value_asc = 0
                endif
                if(metric%computeADC.eq.1) then 
                   allocate(stats%sum(m)%obs_value_adc(LVT_rc%ngrid, selectNlevs(2),&
                        LVT_rc%nadc))
                   stats%sum(m)%obs_value_adc = 0.0
                   allocate(stats%sum(m)%count_obs_value_adc(LVT_rc%ngrid, selectNlevs(2),&
                        LVT_rc%nadc))
                   stats%sum(m)%count_obs_value_adc = 0
                endif
             endif
          endif
       endif
       if(metric%timeOpt.eq.1) then 
          allocate(stats%sum(m)%model_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%sum(m)%model_value_ts = 0.0
          allocate(stats%sum(m)%count_model_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%sum(m)%count_model_value_ts = 0
          
          allocate(stats%sum(m)%tavg_model_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%sum(m)%tavg_model_value_ts = 0.0
          allocate(stats%sum(m)%tavg_count_model_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%sum(m)%tavg_count_model_value_ts = 0
          
          if(LVT_rc%obssource(2).ne."none") then 
             if(selectNlevs(2).ge.1) then 
                allocate(stats%sum(m)%obs_value_ts(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%strat_nlevels))
                stats%sum(m)%obs_value_ts = 0.0
                allocate(stats%sum(m)%count_obs_value_ts(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%strat_nlevels))
                stats%sum(m)%count_obs_value_ts = 0
                
                allocate(stats%sum(m)%tavg_obs_value_ts(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%strat_nlevels))
                stats%sum(m)%tavg_obs_value_ts = 0.0
                allocate(stats%sum(m)%tavg_count_obs_value_ts(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%strat_nlevels))
                stats%sum(m)%tavg_count_obs_value_ts = 0
                
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

  end subroutine LVT_initSum

!BOP
! 
! !ROUTINE: LVT_diagnoseSum
! \label{LVT_diagnoseSum}
!
! !INTERFACE: 
  subroutine LVT_diagnoseSum(pass)
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
!   calculating the sum of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelSum](\ref{diagnoseSingleModelSum})
!     updates the sum computation for a single variable 
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

    if(pass.eq.1) then 
       if(LVT_metrics%sum%selectOpt.eq.1.or.&
            LVT_metrics%sum%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleSum(model,obs,stats,&
                  LVT_metrics%sum)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseSum

!BOP
! 
! !ROUTINE: diagnoseSingleSum
! \label{diagnoseSingleSum}
!
! !INTERFACE: 
  subroutine diagnoseSingleSum(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the sum computation of the 
!   specified variable. 
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
                         stats%sum(m)%model_value_total(t,k,1) = &
                              stats%sum(m)%model_value_total(t,k,1) + &
                              model%value(t,m,m_k)*model%count(t,m,m_k)
                         stats%sum(m)%count_model_value_total(t,k,1) = &
                              stats%sum(m)%count_model_value_total(t,k,1) + 1
                      endif
                   endif

                   if(metric%timeOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         stats%sum(m)%model_value_ts(t,k,1) = &
                              stats%sum(m)%model_value_ts(t,k,1)+&
                              model%value(t,m,m_k)*model%count(t,m,m_k)
                         stats%sum(m)%count_model_value_ts(t,k,1) = & 
                              stats%sum(m)%count_model_value_ts(t,k,1)+1
                      endif
                   endif
                   if(metric%computeSC.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                              tind)
                         stats%sum(m)%model_value_asc(t,k,tind) = &
                              stats%sum(m)%model_value_asc(t,k,tind)+&
                              model%value(t,m,m_k)*model%count(t,m,m_k)
                         stats%sum(m)%count_model_value_asc(t,k,tind) = &
                              stats%sum(m)%count_model_value_asc(t,k,tind) + 1
                      endif
                   endif
                   if(metric%computeADC.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%sum(m)%model_value_adc(t,k,tind) = &
                              stats%sum(m)%model_value_adc(t,k,tind)+&
                              model%value(t,m,m_k)*model%count(t,m,m_k)
                         stats%sum(m)%count_model_value_adc(t,k,tind) = &
                              stats%sum(m)%count_model_value_adc(t,k,tind) + 1
                      endif
                   endif
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,m,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         if(metric%selectOpt.eq.1) then
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then  
                               stats%sum(m)%model_value_total(t,k,2) = &
                                    stats%sum(m)%model_value_total(t,k,2) + &
                                    model%value(t,m,m_k)*model%count(t,m,m_k)
                               stats%sum(m)%count_model_value_total(t,k,2) = &
                                    stats%sum(m)%count_model_value_total(t,k,2) + 1
                            endif
                         endif
                         if(metric%timeOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               stats%sum(m)%model_value_ts(t,k,2) = &
                                    stats%sum(m)%model_value_ts(t,k,2)+&
                                    model%value(t,m,m_k)*model%count(t,m,m_k)
                               stats%sum(m)%count_model_value_ts(t,k,2) = & 
                                    stats%sum(m)%count_model_value_ts(t,k,2)+1
                            endif
                         endif
                      elseif(LVT_stats%strat_var(t,m,k).le.&
                           LVT_rc%strat_var_threshold) then
                         if(metric%selectOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               stats%sum(m)%model_value_total(t,k,3) = &
                                    stats%sum(m)%model_value_total(t,k,3) + &
                                    model%value(t,m,m_k)*model%count(t,m,m_k)
                               stats%sum(m)%count_model_value_total(t,k,3) = &
                                    stats%sum(m)%count_model_value_total(t,k,3) + 1
                            endif
                         endif
                         if(metric%timeOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               stats%sum(m)%model_value_ts(t,k,3) = &
                                    stats%sum(m)%model_value_ts(t,k,3)+&
                                    model%value(t,m,m_k)*model%count(t,m,m_k)
                               stats%sum(m)%count_model_value_ts(t,k,3) = & 
                                    stats%sum(m)%count_model_value_ts(t,k,3)+1
                            endif
                         endif
                      endif
                   endif
                endif

                if(LVT_rc%obssource(2).ne."none") then
                   if(obs%selectNlevs.ge.1) then
                      if(obs%count(t,m,o_k).gt.0) then 
                         if(metric%selectOpt.eq.1.and.&
                              obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            stats%sum(m)%obs_value_total(t,k,1) = &
                                 stats%sum(m)%obs_value_total(t,k,1) + &
                                 obs%value(t,m,o_k)*obs%count(t,m,o_k)
                            stats%sum(m)%count_obs_value_total(t,k,1) = &
                                 stats%sum(m)%count_obs_value_total(t,k,1) + 1
                         endif

                         if(metric%timeOpt.eq.1.and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            stats%sum(m)%obs_value_ts(t,k,1) = &
                                 stats%sum(m)%obs_value_ts(t,k,1)+&
                                 obs%value(t,m,o_k)*obs%count(t,m,o_k)
                            stats%sum(m)%count_obs_value_ts(t,k,1) = & 
                                 stats%sum(m)%count_obs_value_ts(t,k,1)+1
                         endif
                         if(metric%computeSC.eq.1.and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,tind)
                            stats%sum(m)%obs_value_asc(t,k,tind) = &
                                 stats%sum(m)%obs_value_asc(t,k,tind)+&
                                 obs%value(t,m,o_k)*obs%count(t,m,o_k)
                            stats%sum(m)%count_obs_value_asc(t,k,tind) = &
                                 stats%sum(m)%count_obs_value_asc(t,k,tind) + 1
                         endif
                         if(metric%computeADC.eq.1.and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            call LVT_getADCTimeIndex(tind)
                            stats%sum(m)%obs_value_adc(t,k,tind) = &
                                 stats%sum(m)%obs_value_adc(t,k,tind)+&
                                 obs%value(t,m,o_k)*obs%count(t,m,o_k)
                            stats%sum(m)%count_obs_value_adc(t,k,tind) = &
                                 stats%sum(m)%count_obs_value_adc(t,k,tind) + 1
                         endif
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then
                               if(metric%selectOpt.eq.1 &
                                    .and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%sum(m)%obs_value_total(t,k,2) = &
                                       stats%sum(m)%obs_value_total(t,k,2) + &
                                       obs%value(t,m,o_k)*obs%count(t,m,o_k)
                                  stats%sum(m)%count_obs_value_total(t,k,2) = &
                                       stats%sum(m)%count_obs_value_total(t,k,2) + 1
                               endif

                               if(metric%timeOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%sum(m)%obs_value_ts(t,k,2) = &
                                       stats%sum(m)%obs_value_ts(t,k,2)+&
                                       obs%value(t,m,o_k)*obs%count(t,m,o_k)
                                  stats%sum(m)%count_obs_value_ts(t,k,2) = & 
                                       stats%sum(m)%count_obs_value_ts(t,k,2)+1
                               endif
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then
                               if(metric%selectOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%sum(m)%obs_value_total(t,k,3) = &
                                       stats%sum(m)%obs_value_total(t,k,3) + &
                                       obs%value(t,m,o_k)*obs%count(t,m,o_k)
                                  stats%sum(m)%count_obs_value_total(t,k,3) = &
                                       stats%sum(m)%count_obs_value_total(t,k,3) + 1
                               endif

                               if(metric%timeOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%sum(m)%obs_value_ts(t,k,3) = &
                                       stats%sum(m)%obs_value_ts(t,k,3)+&
                                       obs%value(t,m,o_k)*obs%count(t,m,o_k)
                                  stats%sum(m)%count_obs_value_ts(t,k,3) = & 
                                       stats%sum(m)%count_obs_value_ts(t,k,3)+1
                               endif

                            endif
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       end do
    end if
  end subroutine diagnoseSingleSum

!BOP
! 
! !ROUTINE: LVT_computeSum
! \label{LVT_computeSum}
!
! !INTERFACE: 
  subroutine LVT_computeSum(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the sum values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelSum](\ref{computeSingleModelSum})
!     computes the sum values for a single variable
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
    
    if(pass.eq.1) then 
       if(LVT_metrics%sum%selectOpt.eq.1.or.&
            LVT_metrics%sum%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%sum%timeOpt.eq.1.and.&
                  LVT_metrics%sum%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sum%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%sum%ftn_ts_loc(i,1),200,advance='no') &
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
             
             call computeSingleSum(alarm,model,obs,stats,&
                  LVT_metrics%sum)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%sum%timeOpt.eq.1.and.&
                  LVT_metrics%sum%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sum%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%sum%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif
             endif

          endif
       endif
    endif
  end subroutine LVT_computeSum
  

!BOP
! 
! !ROUTINE: computeSingleSum
! \label{computeSingleSum}
!
! !INTERFACE: 
  subroutine computeSingleSum(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the sum values
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

    real,    allocatable :: tavg_model_value_ts(:,:,:)
    real,    allocatable :: tavg_obs_value_ts(:,:,:)
    real,    allocatable :: tavg_obs_value_stdev_ts(:,:,:)
    
    real,    allocatable :: model_value_asc(:,:,:)
    integer, allocatable :: count_model_value_asc(:,:,:)
    real,    allocatable :: obs_value_asc(:,:,:)

    real,    allocatable :: model_value_adc(:,:,:)
    real,    allocatable :: obs_value_adc(:,:,:)
    
    real,    allocatable :: model_value_avg(:,:,:)
    real,    allocatable :: obs_value_avg(:,:,:)   


    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%sum(m)%count_model_value_ts(t,k,l).gt.0) then 
                         stats%sum(m)%model_value_ts(t,k,l) = &
                              stats%sum(m)%model_value_ts(t,k,l)
                         
                         stats%sum(m)%tavg_model_value_ts(t,k,l) = & 
                              stats%sum(m)%tavg_model_value_ts(t,k,l) + & 
                              stats%sum(m)%model_value_ts(t,k,l)
                         stats%sum(m)%tavg_count_model_value_ts(t,k,l) = & 
                              stats%sum(m)%tavg_count_model_value_ts(t,k,l) + 1 
                      else
                         stats%sum(m)%model_value_ts(t,k,l) = LVT_rc%udef
                      endif
                      if(LVT_rc%obssource(2).ne."none") then 
                         if(obs%selectNlevs.ge.1) then 
                            if(stats%sum(m)%count_obs_value_ts(t,k,l).gt.0) then 
                               stats%sum(m)%obs_value_ts(t,k,l) = &
                                    stats%sum(m)%obs_value_ts(t,k,l)
                               
                               stats%sum(m)%tavg_obs_value_ts(t,k,l) = & 
                                    stats%sum(m)%tavg_obs_value_ts(t,k,l) + & 
                                    stats%sum(m)%obs_value_ts(t,k,l)
                               stats%sum(m)%tavg_count_obs_value_ts(t,k,l) = & 
                                    stats%sum(m)%tavg_count_obs_value_ts(t,k,l) + 1 
                            else
                               stats%sum(m)%obs_value_ts(t,k,l) = LVT_rc%udef
                            endif
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
                         if(stats%sum(m)%tavg_count_model_value_ts(t,k,l).gt.0) then 
                            stats%sum(m)%tavg_model_value_ts(t,k,l) = & 
                                 stats%sum(m)%tavg_model_value_ts(t,k,l) / & 
                                 stats%sum(m)%tavg_count_model_value_ts(t,k,l) 
                         endif
                      enddo
                   enddo
                enddo
             enddo
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                do t=1,LVT_rc%ngrid
                   do m=1,LVT_rc%nensem                
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            if(stats%sum(m)%tavg_count_obs_value_ts(t,k,l).gt.0) then 
                               stats%sum(m)%tavg_obs_value_ts(t,k,l) = & 
                                    stats%sum(m)%tavg_obs_value_ts(t,k,l) / & 
                                    stats%sum(m)%tavg_count_obs_value_ts(t,k,l) 
                            endif
                         enddo
                      enddo
                   enddo
                enddo
                
             endif

             if(metric%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                         call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                              model,&
                              LVT_rc%ngrid,&
                              stats%sum(m)%tavg_model_value_ts,&
                              stats%sum(m)%tavg_count_model_value_ts,&
                              LVT_rc%ngrid,&
                              stats%sum(m)%tavg_obs_value_ts,&
                              stats%sum(m)%tavg_count_obs_value_ts)
                      else
                         call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                              model,&
                              LVT_rc%ngrid,&
                              stats%sum(m)%tavg_model_value_ts,&
                              stats%sum(m)%tavg_count_model_value_ts)
                      endif
                   enddo
                else
                   allocate(tavg_model_value_ts(LVT_rc%ngrid,&
                        model%vlevels,&
                        LVT_rc%strat_nlevels))
                   tavg_model_value_ts = 0.0
                   
                   if(LVT_rc%obssource(2).ne."none".and.&
                        obs%selectNlevs.ge.1) then 
                      
                      allocate(tavg_obs_value_ts(LVT_rc%ngrid,&
                           obs%vlevels,&
                           LVT_rc%strat_nlevels))
                      tavg_obs_value_ts = 0.0
                      
                      do m=1,LVT_rc%nensem
                         do t=1,LVT_rc%ngrid
                            do k=1,model%selectNlevs
                               do l=1,LVT_rc%strat_nlevels
                                  tavg_model_value_ts(t,k,l) = &
                                       tavg_model_value_ts(t,k,l) + &
                                       stats%sum(m)%tavg_model_value_ts(t,k,l)
                                  tavg_obs_value_ts(t,k,l) = &
                                       tavg_obs_value_ts(t,k,l) + &
                                       stats%sum(m)%tavg_obs_value_ts(t,k,l)
                                  
                                  
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
                           stats%sum(1)%tavg_count_model_value_ts,&
                           LVT_rc%ngrid,&
                           tavg_obs_value_ts,&
                           stats%sum(m)%tavg_count_obs_value_ts)
                   else
                      do m=1,LVT_rc%nensem
                         do t=1,LVT_rc%ngrid
                            do k=1,model%selectNlevs
                               do l=1,LVT_rc%strat_nlevels
                                  tavg_model_value_ts(t,k,l) = &
                                       tavg_model_value_ts(t,k,l) + &
                                       stats%sum(m)%tavg_model_value_ts(t,k,l)
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
                           stats%sum(1)%tavg_count_model_value_ts)
                      
                   endif
                   deallocate(tavg_model_value_ts)

                   if(LVT_rc%obssource(2).ne."none".and.&
                        obs%selectNlevs.ge.1) then 
                      deallocate(tavg_obs_value_ts)
                   endif
                endif
             endif
          end if
       endif
    end if
       
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%sum(m)%count_model_value_total(t,k,l).gt.&
                           LVT_rc%obsCountThreshold) then 
                         stats%sum(m)%model_value_total(t,k,l) = &
                              stats%sum(m)%model_value_total(t,k,l)
                      else
                         stats%sum(m)%model_value_total(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                   if(metric%computeSC.eq.1) then 
                      do l=1, LVT_rc%nasc
                         if(stats%sum(m)%count_model_value_asc(t,k,l).gt.&
                              LVT_rc%SCCountThreshold) then 
                            stats%sum(m)%model_value_asc(t,k,l) = &
                                 stats%sum(m)%model_value_asc(t,k,l)
                         else
                            stats%sum(m)%model_value_asc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   
                   if(metric%computeADC.eq.1) then 
                      do l=1, LVT_rc%nadc
                         if(stats%sum(m)%count_model_value_adc(t,k,l).gt.&
                              LVT_rc%ADCCountThreshold) then 
                            stats%sum(m)%model_value_adc(t,k,l) = &
                                 stats%sum(m)%model_value_adc(t,k,l)
                         else
                            stats%sum(m)%model_value_adc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then
                         do l=1,LVT_rc%strat_nlevels                      
                            if(stats%sum(m)%count_obs_value_total(t,k,l).gt.&
                                 LVT_rc%obsCountThreshold) then 
                               stats%sum(m)%obs_value_total(t,k,l) = &
                                    stats%sum(m)%obs_value_total(t,k,l)
                            else
                               stats%sum(m)%obs_value_total(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                         if(metric%computeSC.eq.1) then
                            do l=1, LVT_rc%nasc 
                               if(stats%sum(m)%count_obs_value_asc(t,k,l).gt.&
                                    LVT_rc%SCCountThreshold) then 
                                  stats%sum(m)%obs_value_asc(t,k,l) = &
                                       stats%sum(m)%obs_value_asc(t,k,l)
                               else
                                  stats%sum(m)%obs_value_asc(t,k,l) = LVT_rc%udef
                               endif
                            enddo
                         endif
                         
                         if(metric%computeADC.eq.1) then
                            do l=1, LVT_rc%nadc  
                               if(stats%sum(m)%count_obs_value_adc(t,k,l).gt.&
                                    LVT_rc%ADCCountThreshold) then 
                                  stats%sum(m)%obs_value_adc(t,k,l) = &
                                       stats%sum(m)%obs_value_adc(t,k,l)
                               else
                                  stats%sum(m)%obs_value_adc(t,k,l) = LVT_rc%udef
                               endif
                            enddo
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%sum(m)%model_value_total(:,k,l),&
                        LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%sum(m)%model_value_ci(k,l))
                enddo
             enddo
             if(LVT_rc%obssource(2).ne."none") then 
                do k=1,obs%selectNlevs
                   do l=1, LVT_rc%strat_nlevels
                      call LVT_computeCI(stats%sum(m)%obs_value_total(:,k,l),&
                           LVT_rc%ngrid,&
                           LVT_rc%pval_CI,stats%sum(m)%obs_value_ci(k,l))
                   enddo
                enddo
             endif
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
                              stats%sum(m)%model_value_total(t,k,l)
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
                                 stats%sum(m)%obs_value_total(t,k,l)
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
                         if(stats%sum(m)%model_value_asc(t,k,l).ne.LVT_rc%udef) then 
                            model_value_asc(t,k,l) = &
                                 model_value_asc(t,k,l) + & 
                                 stats%sum(m)%model_value_asc(t,k,l)
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
                                 stats%sum(m)%obs_value_asc(t,k,l)
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
                     stats%sum(1)%count_model_value_asc,&
                     LVT_rc%ngrid,obs_value_asc,&
                     stats%sum(1)%count_obs_value_asc)   

                deallocate(obs_value_asc)
       
             else
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_asc,&
                     stats%sum(1)%count_model_value_asc)

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
                              stats%sum(m)%model_value_adc(t,k,l)
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
                                 stats%sum(m)%obs_value_adc(t,k,l)
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
                     stats%sum(1)%count_model_value_adc,&
                     LVT_rc%ngrid,obs_value_adc,&
                     stats%sum(1)%count_obs_value_adc)      
                
                deallocate(obs_value_adc)
             else
                call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_adc,&
                     stats%sum(1)%count_model_value_adc)

                deallocate(model_value_adc)
             endif
          endif
       endif
    endif

  end subroutine computeSingleSum

!BOP
! 
! !ROUTINE: LVT_writeMetric_Sum
! \label(LVT_writeMetric_Sum)
!
! !INTERFACE:
  subroutine LVT_writeMetric_Sum(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%sum%npass) then 

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
                              stats%sum(m)%tavg_model_value_ts(:,k,l)
                         count_model_value_ts(:,m) = & 
                              stats%sum(m)%tavg_count_model_value_ts(:,k,l)
                       if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                          obs_value_ts(:,m) = &
                               stats%sum(m)%tavg_obs_value_ts(:,k,l)
                          count_obs_value_ts(:,m) = & 
                               stats%sum(m)%tavg_count_obs_value_ts(:,k,l)
                       endif
                    enddo
                    if(LVT_metrics%sum%timeOpt.eq.1) then 
                       
                       call LVT_writevar_gridded(LVT_metrics%sum%ftn_ts, &
                            model_value_ts(:,:),&
                            stats%vid_ts(LVT_SUMid,1),k)
                       
                       call LVT_writevar_gridded(LVT_metrics%sum%ftn_ts, &
                            real(count_model_value_ts(:,:)),&
                            stats%vid_count_ts(LVT_SUMid,1),k)
                    endif
                    
                    if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                       call LVT_writevar_gridded(LVT_metrics%sum%ftn_ts, &
                            obs_value_ts(:,:),&
                            stats%vid_ts(LVT_SUMid,2),k)
                       call LVT_writevar_gridded(LVT_metrics%sum%ftn_ts, &
                            real(count_obs_value_ts(:,:)),&
                            stats%vid_count_ts(LVT_SUMid,2),k)
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
          if(pass.eq.LVT_metrics%sum%npass) then
             if(stats%selectOpt.eq.1) then

              allocate(model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
              allocate(count_model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
              allocate(model_value_ci(LVT_rc%nensem))

              if(LVT_metrics%sum%computeSC.eq.1) then 
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
              if(LVT_metrics%sum%computeADC.eq.1) then 
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
                              stats%sum(m)%model_value_total(:,k,l)
                         count_model_value_total(:,m) = & 
                              stats%sum(m)%count_model_value_total(:,k,l)
                         model_value_ci(m) = stats%sum(m)%model_value_ci(k,l)
                         
                         if(LVT_metrics%sum%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               model_value_asc(:,m,tind) = & 
                                    stats%sum(m)%model_value_asc(:,k,tind)
                            enddo
                         endif

                         if(LVT_metrics%sum%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               model_value_adc(:,m,tind) = & 
                                    stats%sum(m)%model_value_adc(:,k,tind)
                            enddo
                         endif

                         if(LVT_rc%obssource(2).ne."none".and.&
                              obs%selectNlevs.ge.1) then 
                            obs_value_total(:,m) = &
                                 stats%sum(m)%obs_value_total(:,k,l)
                            count_obs_value_total(:,m) = & 
                                 stats%sum(m)%count_obs_value_total(:,k,l)
                            obs_value_ci(m) = stats%sum(m)%obs_value_ci(k,l)

                            if(LVT_metrics%sum%computeSC.eq.1) then 
                               do tind = 1,LVT_rc%nasc
                                  obs_value_asc(:,m,tind) = & 
                                       stats%sum(m)%obs_value_asc(:,k,tind)
                               enddo
                            endif
                            if(LVT_metrics%sum%computeADC.eq.1) then 
                               do tind = 1,LVT_rc%nadc
                                  obs_value_adc(:,m,tind) = & 
                                       stats%sum(m)%obs_value_adc(:,k,tind)
                               enddo
                            endif
                            
                         endif
                      enddo
                      if(LVT_metrics%sum%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%sum%ftn_total, &
                              model_value_total(:,:),&
                              stats%vid_total(LVT_SUMid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%sum%ftn_total, &
                              real(count_model_value_total(:,:)),&
                              stats%vid_count_total(LVT_SUMid,1),k)

                         if(LVT_rc%obssource(2).ne."none".and.&
                              obs%selectNlevs.ge.1) then 
                            call LVT_writevar_gridded(LVT_metrics%sum%ftn_total, &
                                 obs_value_total(:,:),&
                                 stats%vid_total(LVT_SUMid,2),k)
                            call LVT_writevar_gridded(LVT_metrics%sum%ftn_total, &
                                 real(count_obs_value_total(:,:)),&
                                 stats%vid_count_total(LVT_SUMid,2),k)
                         endif
                         
                         if(LVT_metrics%sum%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%sum%ftn_total,&
                                    model_value_asc(:,:,tind),&
                                    stats%vid_sc_total(tind,LVT_SUMid,1),k)
                            enddo
                            if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                               do tind = 1,LVT_rc%nasc
                                  call LVT_writevar_gridded(&
                                       LVT_metrics%sum%ftn_total,&
                                       obs_value_asc(:,:,tind),&
                                       stats%vid_sc_total(tind,LVT_SUMid,2),k)
                               enddo
                            endif
                         endif
                         if(LVT_metrics%sum%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%sum%ftn_total,&
                                    model_value_adc(:,:,tind),&
                                    stats%vid_adc_total(tind,LVT_SUMid,1),k)
                            enddo
                            if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                               do tind = 1,LVT_rc%nadc
                                  call LVT_writevar_gridded(&
                                       LVT_metrics%sum%ftn_total,&
                                       obs_value_adc(:,:,tind),&
                                       stats%vid_adc_total(tind,LVT_SUMid,2),k)
                               enddo
                            endif
                         endif
                         call LVT_writeSummaryStats(&
                              LVT_metrics%sum%ftn_summ,&
                              l,&
                              LVT_metrics%sum%short_name,&
                              LVT_rc%ngrid,&
                              model_value_total(:,:), &
                              count_model_value_total(:,:),&
                              stats%standard_name,&
                              model_value_ci(:))
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            call LVT_writeSummaryStats(&
                                 LVT_metrics%sum%ftn_summ,&
                                 l,&
                                 LVT_metrics%sum%short_name,&
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

                if(LVT_metrics%sum%computeSC.eq.1) then 
                   deallocate(model_value_asc)
                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                      deallocate(obs_value_asc)
                   endif
                endif
                if(LVT_metrics%sum%computeADC.eq.1) then 
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

  end subroutine LVT_writeMetric_Sum

!BOP
! 
! !ROUTINE: LVT_resetMetric_Sum
! \label(LVT_resetMetric_Sum)
!
! !INTERFACE:
  subroutine LVT_resetMetric_Sum(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical                  :: alarm
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

    integer                :: i,k,l,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats


    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream1Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                if(LVT_metrics%sum%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%sum(m)%model_value_ts(:,k,l) = 0.0
                      stats%sum(m)%count_model_value_ts(:,k,l)=0 
                      if(alarm) then
                         stats%sum(m)%tavg_model_value_ts(:,k,l) = 0.0
                         stats%sum(m)%tavg_count_model_value_ts(:,k,l)=0 
                      endif
                   enddo
                   if(LVT_rc%obssource(2).ne."none".and.&
                        obs%selectNlevs.ge.1) then 
                      do l=1,LVT_rc%strat_nlevels
                         stats%sum(m)%obs_value_ts(:,k,l) = 0.0
                         stats%sum(m)%count_obs_value_ts(:,k,l)=0 
                         if(alarm) then 
                            stats%sum(m)%tavg_obs_value_ts(:,k,l) = 0.0
                            stats%sum(m)%tavg_count_obs_value_ts(:,k,l)=0 
                         endif
                      enddo
                   endif
                endif
             enddo
          enddo
       endif
       
       model => model%next
       obs => obs%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_Sum

!BOP
! 
! !ROUTINE: LVT_writerestart_Sum
! 
! !INTERFACE:
  subroutine LVT_writerestart_Sum(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Sum metric computations
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
       if(LVT_metrics%sum%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_writevar_restart(ftn,&
                           stats%sum(m)%model_value_total(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%sum(m)%count_model_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%sum%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc         
                         call LVT_writevar_restart(ftn,&
                              stats%sum(m)%model_value_asc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%sum(m)%count_model_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%sum%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc         
                         call LVT_writevar_restart(ftn,&
                              stats%sum(m)%model_value_adc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%sum(m)%count_model_value_adc(:,k,l))
                      enddo
                   enddo
                endif
                
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            call LVT_writevar_restart(ftn,&
                                 stats%sum(m)%obs_value_total(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%sum(m)%count_obs_value_total(:,k,l))
                         enddo
                      enddo
                      
                      if(LVT_metrics%sum%computeSC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_writevar_restart(ftn,&
                                    stats%sum(m)%obs_value_asc(:,k,l))
                               call LVT_writevar_restart(ftn,&
                                    stats%sum(m)%count_obs_value_asc(:,k,l))
                            enddo
                         enddo
                      endif
                      if(LVT_metrics%sum%computeADC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_writevar_restart(ftn,&
                                    stats%sum(m)%obs_value_adc(:,k,l))
                               call LVT_writevar_restart(ftn,&
                                    stats%sum(m)%count_obs_value_adc(:,k,l))
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
    end do
  end subroutine LVT_writerestart_Sum


!BOP
! 
! !ROUTINE: LVT_readrestart_Sum
! 
! !INTERFACE:
  subroutine LVT_readrestart_Sum(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Sum metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    integer              :: k,l,m,index
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%sum%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_readvar_restart(ftn,&
                           stats%sum(m)%model_value_total(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%sum(m)%count_model_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%sum%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc         
                         call LVT_readvar_restart(ftn,&
                              stats%sum(m)%model_value_asc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%sum(m)%count_model_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%sum%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc         
                         call LVT_readvar_restart(ftn,&
                              stats%sum(m)%model_value_adc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%sum(m)%count_model_value_adc(:,k,l))
                      enddo
                   enddo
                endif

                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            call LVT_readvar_restart(ftn,&
                                 stats%sum(m)%obs_value_total(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%sum(m)%count_obs_value_total(:,k,l))
                         enddo
                      enddo

                      if(LVT_metrics%sum%computeSC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_readvar_restart(ftn,&
                                    stats%sum(m)%obs_value_asc(:,k,l))
                               call LVT_readvar_restart(ftn,&
                                    stats%sum(m)%count_obs_value_asc(:,k,l))
                            enddo
                         enddo
                      endif
                      if(LVT_metrics%sum%computeADC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_readvar_restart(ftn,&
                                    stats%sum(m)%obs_value_adc(:,k,l))
                               call LVT_readvar_restart(ftn,&
                                    stats%sum(m)%count_obs_value_adc(:,k,l))
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
    end do
  end subroutine LVT_readrestart_Sum


end module LVT_SumMod
