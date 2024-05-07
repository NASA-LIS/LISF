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
! !MODULE: LVT_MaxtimeMod
! \label(LVT_MaxtimeMod)
!
! !INTERFACE:
module LVT_MaxtimeMod
! 
! !USES:  
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_timeMgrMod
  use LVT_CIMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the computations required to compute maximum values
!  (temporally) of desired variables from the LIS output
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
  public :: LVT_initMaxtime
  public :: LVT_diagnoseMaxtime
  public :: LVT_computeMaxtime
  public :: LVT_writeMetric_Maxtime
  public :: LVT_resetMetric_Maxtime
  public :: LVT_writerestart_Maxtime
  public :: LVT_readrestart_Maxtime
!EOP
  
  type, public :: maxtimedec
    character*5          :: time_option
  end type maxtimedec

  type(maxtimedec), save :: LVT_maxtime_struc

  private

contains
  subroutine LVT_initMaxtime(model, obs, stats,metric)
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!  This routine initializes the datastructures required to support
!  the Maxtime computations. 
!
!EOP
    integer      :: rc,m
    if(metric%selectOpt.eq.1) then 

       call ESMF_ConfigGetAttribute(LVT_config,LVT_maxtime_struc%time_option, &
            label="Time specification option for MaxTime metric:",rc=rc)
       if(rc.ne.0) then 
          write(LVT_logunit,*) "[ERR] "
          write(LVT_logunit,*) "[ERR] Time specification option for MaxTime metric: not specified"
          write(LVT_logunit,*) "[ERR] options are 'doy/da/hr/mn/ss' "
          write(LVT_logunit,*) "[ERR] "
          call LVT_endrun()
       endif

       allocate(stats%maxtime(LVT_rc%nensem))

       do m=1,LVT_rc%nensem       
          allocate(stats%maxtime(m)%model_value_total(LVT_rc%ngrid, &
               model%selectNlevs, &
               LVT_rc%strat_nlevels))
          stats%maxtime(m)%model_value_total = -1E10
          allocate(stats%maxtime(m)%model_value_max_total(LVT_rc%ngrid, &
               model%selectNlevs, &
               LVT_rc%strat_nlevels))
          stats%maxtime(m)%model_value_max_total = -1E10
          allocate(stats%maxtime(m)%model_value_ci(model%selectNlevs,&
               LVT_rc%strat_nlevels))
          stats%maxtime(m)%model_value_ci = LVT_rc%udef
          
          if(metric%computeSC.eq.1) then 
             allocate(stats%maxtime(m)%model_value_asc(LVT_rc%ngrid, &
                  model%selectNlevs,&
                  LVT_rc%nasc))
             stats%maxtime(m)%model_value_asc = -1E10
             allocate(stats%maxtime(m)%model_value_max_asc(LVT_rc%ngrid, &
                  model%selectNlevs,&
                  LVT_rc%nasc))
             stats%maxtime(m)%model_value_max_asc = -1E10
          endif
          if(metric%computeADC.eq.1) then 
             allocate(stats%maxtime(m)%model_value_adc(LVT_rc%ngrid, &
                  obs%selectNlevs,&
                  LVT_rc%nadc))
             stats%maxtime(m)%model_value_adc = -1E10
             allocate(stats%maxtime(m)%model_value_adc(LVT_rc%ngrid, &
                  obs%selectNlevs,&
                  LVT_rc%nadc))
             stats%maxtime(m)%model_value_max_adc = -1E10
          endif
          
          if(LVT_rc%obssource(2).ne."none") then 
             if(obs%selectNlevs.ge.1) then 
                allocate(stats%maxtime(m)%obs_value_total(LVT_rc%ngrid, &
                     obs%selectNlevs, &
                     LVT_rc%strat_nlevels))
                stats%maxtime(m)%obs_value_total = -1E10
                
                allocate(stats%maxtime(m)%obs_value_max_total(LVT_rc%ngrid,&
                     obs%selectNlevs, &
                     LVT_rc%strat_nlevels))
                stats%maxtime(m)%obs_value_max_total = -1E10
                allocate(stats%maxtime(m)%obs_value_ci(obs%selectNlevs,&
                     LVT_rc%strat_nlevels))
                stats%maxtime(m)%obs_value_ci = LVT_rc%udef
                
                if(metric%computeSC.eq.1) then 
                   allocate(stats%maxtime(m)%obs_value_asc(LVT_rc%ngrid, &
                        obs%selectNlevs,&
                        LVT_rc%nasc))
                   stats%maxtime(m)%obs_value_asc = -1E10
                   allocate(stats%maxtime(m)%obs_value_max_asc(LVT_rc%ngrid,&
                        obs%selectNlevs,&
                        LVT_rc%nasc))
                   stats%maxtime(m)%obs_value_max_asc = -1E10
                endif
                if(metric%computeADC.eq.1) then 
                   allocate(stats%maxtime(m)%obs_value_adc(LVT_rc%ngrid,&
                        obs%selectNlevs,&
                        LVT_rc%nadc))
                   stats%maxtime(m)%obs_value_adc = -1E10
                   allocate(stats%maxtime(m)%obs_value_max_adc(LVT_rc%ngrid,&
                        obs%selectNlevs,&
                        LVT_rc%nadc))
                   stats%maxtime(m)%obs_value_max_adc = -1E10
                endif
             endif
          endif
          if(metric%timeOpt.eq.1) then 
             allocate(stats%maxtime(m)%model_value_ts(LVT_rc%ngrid, &
                  model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%maxtime(m)%model_value_ts = -1E10
             
             allocate(stats%maxtime(m)%model_value_max_ts(LVT_rc%ngrid,&
                  model%selectNlevs, &
                  LVT_rc%strat_nlevels))
             stats%maxtime(m)%model_value_max_ts = -1E10
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   allocate(stats%maxtime(m)%obs_value_ts(LVT_rc%ngrid, &
                        model%selectNlevs, &
                        LVT_rc%strat_nlevels))
                   stats%maxtime(m)%obs_value_ts = -1E10
                   
                   allocate(stats%maxtime(m)%obs_value_max_ts(LVT_rc%ngrid, &
                        model%selectNlevs, &
                        LVT_rc%strat_nlevels))
                   stats%maxtime(m)%obs_value_max_ts = -1E10
                   
                endif
             endif
             
          endif
       enddo
    end if
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

  end subroutine LVT_initMaxtime

!BOP
! 
! !ROUTINE: LVT_diagnoseMaxtime
! \label{LVT_diagnoseMaxtime}
!
! !INTERFACE: 
  subroutine LVT_diagnoseMaxtime(pass)
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
!   calculating the maxtime of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelMaxtime](\ref{diagnoseSingleModelMaxtime})
!     updates the maxtime computation for a single variable 
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
       if(LVT_metrics%maxtime%selectOpt.eq.1.or.&
            LVT_metrics%maxtime%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMaxtime(model,obs,stats,&
                  LVT_metrics%maxtime)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseMaxtime

!BOP
! 
! !ROUTINE: diagnoseSingleMaxtime
! \label{diagnoseSingleMaxtime}
!
! !INTERFACE: 
  subroutine diagnoseSingleMaxtime(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the maxtime computation of the 
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
    integer    :: t,k,m,m_k,o_k,tind,kk
    integer    :: time

!----------------------------------------------------------------
!It is assumed temporal lag is not chosen in this setup
!----------------------------------------------------------------
    do kk=1,LVT_rc%nDataStreams
       if(LVT_rc%tlag(kk).gt.0) then 
          write(LVT_logunit,*) "[ERR] "
          write(LVT_logunit,*) "[ERR] Non-zero temporal lag specification is not "
          write(LVT_logunit,*) "[ERR] supported for MaxTime metric"
          call LVT_endrun()
       endif
    enddo
    
    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
       if(LVT_maxtime_struc%time_option.eq."doy") then 
          time = LVT_rc%doy
       elseif(LVT_maxtime_struc%time_option.eq."da") then 
          time = LVT_rc%da
       elseif(LVT_maxtime_struc%time_option.eq."hr") then 
          time = LVT_rc%hr
       elseif(LVT_maxtime_struc%time_option.eq."mn") then 
          time = LVT_rc%mn
       elseif(LVT_maxtime_struc%time_option.eq."ss") then 
          time = LVT_rc%ss
       endif
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                
                if(model%count(t,m,m_k).gt.0) then
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         if(model%value(t,m,m_k).gt.stats%maxtime(m)%model_value_max_total(t,k,1)) then 
                            stats%maxtime(m)%model_value_max_total(t,k,1) = &
                                 model%value(t,m,m_k)
                            stats%maxtime(m)%model_value_total(t,k,1) = &
                                 time
                         endif
                      endif
                   endif
                   
                   if(metric%timeOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         if(model%value(t,m,m_k).gt.stats%maxtime(m)%model_value_max_ts(t,k,1)) then 
                            stats%maxtime(m)%model_value_max_ts(t,k,1) = model%value(t,m,m_k)
                            stats%maxtime(m)%model_value_ts(t,k,1) = time
                         endif
                      endif
                   endif
                   if(metric%computeSC.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                              tind)
                         if(model%value(t,m,m_k).gt.&
                              stats%maxtime(m)%model_value_max_asc(t,k,tind)) then
                            stats%maxtime(m)%model_value_max_asc(t,k,tind) = &
                                 model%value(t,m,m_k)
                            stats%maxtime(m)%model_value_asc(t,k,tind) = &
                                 time
                         endif
                      endif
                   endif
                   if(metric%computeADC.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         call LVT_getADCTimeIndex(tind)
                         if(model%value(t,m,m_k).gt.&
                              stats%maxtime(m)%model_value_max_adc(t,k,tind)) then 
                            stats%maxtime(m)%model_value_max_adc(t,k,tind) = & 
                                 model%value(t,m,m_k)
                            stats%maxtime(m)%model_value_adc(t,k,tind) = &
                              time
                         endif
                      endif
                   endif
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,m,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         if(metric%selectOpt.eq.1) then
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then  
                               if(model%value(t,m,m_k).gt.&
                                    stats%maxtime(m)%model_value_max_total(t,k,2)) then
                                  stats%maxtime(m)%model_value_max_total(t,k,2) = & 
                                       model%value(t,m,m_k)
                                  stats%maxtime(m)%model_value_total(t,k,2) = &
                                       time
                               endif
                            endif
                         endif
                         if(metric%timeOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               if(model%value(t,m,m_k).gt.&
                                    stats%maxtime(m)%model_value_max_ts(t,k,2)) then 
                                  stats%maxtime(m)%model_value_max_ts(t,k,2) = & 
                                       model%value(t,m,m_k)
                                  stats%maxtime(m)%model_value_ts(t,k,2) = &
                                       time
                               endif
                            endif
                         endif
                      elseif(LVT_stats%strat_var(t,m,k).le.&
                           LVT_rc%strat_var_threshold) then
                         if(metric%selectOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               if(model%value(t,m,m_k).gt.&
                                    stats%maxtime(m)%model_value_max_total(t,k,3)) then 
                                  stats%maxtime(m)%model_value_max_total(t,k,3) = & 
                                       model%value(t,m,m_k)
                                  stats%maxtime(m)%model_value_total(t,k,3) = &
                                       time
                               endif
                            endif
                         endif
                         if(metric%timeOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               if(model%value(t,m,m_k).gt.&
                                    stats%maxtime(m)%model_value_max_ts(t,k,3)) then 
                                  stats%maxtime(m)%model_value_max_ts(t,k,3) = & 
                                       model%value(t,m,m_k)
                                  stats%maxtime(m)%model_value_ts(t,k,3) = &
                                       time
                               endif
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
                            if(obs%value(t,m,o_k).gt.&
                                 stats%maxtime(m)%obs_value_max_total(t,k,1)) then 
                               stats%maxtime(m)%obs_value_max_total(t,k,1) = & 
                                    obs%value(t,m,o_k)
                               stats%maxtime(m)%obs_value_total(t,k,1) = &
                                    time
                            endif
                         endif
                         
                         if(metric%timeOpt.eq.1.and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            if(obs%value(t,m,o_k).gt.&
                                 stats%maxtime(m)%obs_value_max_ts(t,k,1)) then 
                               stats%maxtime(m)%obs_value_max_ts(t,k,1) = obs%value(t,m,o_k)
                               stats%maxtime(m)%obs_value_ts(t,k,1) = time
                            endif
                         endif
                         if(metric%computeSC.eq.1.and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,tind)
                            if(obs%value(t,m,o_k).gt.stats%maxtime(m)%obs_value_max_asc(t,k,tind)) then 
                               stats%maxtime(m)%obs_value_max_asc(t,k,tind) = &
                                    obs%value(t,m,o_k)
                               stats%maxtime(m)%obs_value_asc(t,k,tind) = &
                                    time
                            endif
                         endif
                         if(metric%computeADC.eq.1.and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            call LVT_getADCTimeIndex(tind)
                            if(obs%value(t,m,o_k).gt.stats%maxtime(m)%obs_value_max_adc(t,k,tind)) then 
                               stats%maxtime(m)%obs_value_max_adc(t,k,tind) = &
                                    obs%value(t,m,o_k)
                               stats%maxtime(m)%obs_value_adc(t,k,tind) = &
                                    time
                            endif
                         endif
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then
                               if(metric%selectOpt.eq.1 &
                                    .and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  if(obs%value(t,m,o_k).gt.stats%maxtime(m)%obs_value_max_total(t,k,2)) then 
                                     stats%maxtime(m)%obs_value_max_total(t,k,2) = &
                                          obs%value(t,m,o_k)
                                     stats%maxtime(m)%obs_value_total(t,k,2) = &
                                          time
                                  endif
                               endif
                               
                               if(metric%timeOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  if(obs%value(t,m,o_k).gt.stats%maxtime(m)%obs_value_max_ts(t,k,2)) then 
                                     stats%maxtime(m)%obs_value_max_ts(t,k,2) = & 
                                          obs%value(t,m,o_k)
                                     stats%maxtime(m)%obs_value_ts(t,k,2) = &
                                          time
                                  endif
                               endif
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then
                               if(metric%selectOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  if(obs%value(t,m,o_k).gt.stats%maxtime(m)%obs_value_max_total(t,k,3)) then 
                                     stats%maxtime(m)%obs_value_max_total(t,k,3) = & 
                                          obs%value(t,m,o_k)
                                     stats%maxtime(m)%obs_value_total(t,k,3) = &
                                          time
                                  endif
                               endif
                               
                               if(metric%timeOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  if(obs%value(t,m,o_k).gt.stats%maxtime(m)%obs_value_max_ts(t,k,3)) then 
                                     stats%maxtime(m)%obs_value_max_ts(t,k,3) = & 
                                          obs%value(t,m,o_k)
                                     stats%maxtime(m)%obs_value_ts(t,k,3) = &
                                          time
                                  endif
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
  end subroutine diagnoseSingleMaxtime

!BOP
! 
! !ROUTINE: LVT_computeMaxtime
! \label{LVT_computeMaxtime}
!
! !INTERFACE: 
  subroutine LVT_computeMaxtime(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the maxtime values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelMaxtime](\ref{computeSingleModelMaxtime})
!     computes the maxtime values for a single variable
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
       if(LVT_metrics%maxtime%selectOpt.eq.1.or.&
            LVT_metrics%maxtime%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%maxtime%timeOpt.eq.1.and.&
                  LVT_metrics%maxtime%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%maxtime%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%maxtime%ftn_ts_loc(i,1),200,advance='no') &
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
             
             call computeSingleMaxtime(alarm,model,obs,stats,&
                  LVT_metrics%maxtime)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%maxtime%timeOpt.eq.1.and.&
                  LVT_metrics%maxtime%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%maxtime%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%maxtime%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif

             endif
          endif
       endif
    endif
  end subroutine LVT_computeMaxtime
  

!BOP
! 
! !ROUTINE: computeSingleMaxtime
! \label{computeSingleMaxtime}
!
! !INTERFACE: 
  subroutine computeSingleMaxtime(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the maxtime values
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
    integer     :: count_model_asc(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nasc)
    integer     :: count_obs_asc(LVT_rc%ngrid,obs%selectNlevs,LVT_rc%nasc)
    integer     :: count_model_adc(LVT_rc%ngrid,model%selectNlevs,LVT_rc%nadc)
    integer     :: count_obs_adc(LVT_rc%ngrid,obs%selectNlevs,LVT_rc%nadc)
    
    integer    :: count_model(LVT_rc%ngrid, model%selectNlevs, LVT_rc%strat_nlevels)
    integer    :: count_obs(LVT_rc%ngrid, model%selectNlevs, LVT_rc%strat_nlevels)

    real,    allocatable :: model_value_ts(:,:,:)
    real,    allocatable :: obs_value_ts(:,:,:)
    
    real,    allocatable :: model_value_asc(:,:,:)
    integer, allocatable :: count_model_value_asc(:,:,:)
    real,    allocatable :: obs_value_asc(:,:,:)

    real,    allocatable :: model_value_adc(:,:,:)
    real,    allocatable :: obs_value_adc(:,:,:)
    
    real,    allocatable :: model_value_avg(:,:,:)
    real,    allocatable :: obs_value_avg(:,:,:)   

    count_model = 1
    count_obs = 1
    count_model_asc = 1
    count_model_adc = 1
    count_obs_adc = 1

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%maxtime(m)%model_value_ts(t,k,l).gt.-1E10) then 
                         
                      else
                         stats%maxtime(m)%model_value_ts(t,k,l) = LVT_rc%udef
                      endif
                      if(LVT_rc%obssource(2).ne."none") then 
                         if(obs%selectNlevs.ge.1) then 
                            if(stats%maxtime(m)%obs_value_ts(t,k,l).gt.-1E10) then 
                               
                            else
                               stats%maxtime(m)%obs_value_ts(t,k,l) = LVT_rc%udef
                            endif
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo
          if(alarm) then 

             if(metric%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then
                         call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                              model,&
                              LVT_rc%ngrid,&
                              stats%maxtime(m)%model_value_ts,&
                              count_model,&
                              LVT_rc%ngrid,&
                              stats%maxtime(m)%obs_value_ts,&
                              count_obs)
                      else
                         call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                              model,&
                              LVT_rc%ngrid,&
                              stats%maxtime(m)%model_value_ts,&
                              count_model)
                      endif
                   enddo
                else
                   allocate(model_value_ts(LVT_rc%ngrid, &
                        model%vlevels, &
                        LVT_rc%strat_nlevels))
                   model_value_ts = 0.0

                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 

                      allocate(obs_value_ts(LVT_rc%ngrid, &
                           obs%vlevels, &
                           LVT_rc%strat_nlevels))
                      obs_value_ts = 0.0

                         
                      do m=1,LVT_rc%nensem
                         do t=1,LVT_rc%ngrid
                            do k=1,model%selectNlevs
                               do l=1,LVT_rc%strat_nlevels
                                  model_value_ts(t,k,l) = &
                                       model_value_ts(t,k,l) + &
                                       stats%maxtime(m)%model_value_ts(t,k,l)
                                  obs_value_ts(t,k,l) = &
                                       obs_value_ts(t,k,l) + &
                                       stats%maxtime(m)%obs_value_ts(t,k,l)
                               enddo
                            enddo
                         enddo
                      enddo
                            
                      do t=1,LVT_rc%ngrid
                         do k=1,model%selectNlevs
                            do l=1,LVT_rc%strat_nlevels
                               model_value_ts(t,k,l) = &
                                    model_value_ts(t,k,l)/LVT_rc%nensem
                               obs_value_ts(t,k,l) = &
                                    obs_value_ts(t,k,l)/LVT_rc%nensem
                            enddo
                         enddo
                      enddo
                         
                      call LVT_writeTSinfo(metric%ftn_ts_loc(:,1),&
                           model,&
                           LVT_rc%ngrid,&
                           model_value_ts,&
                           count_model,&
                           LVT_rc%ngrid,&
                           obs_value_ts,&
                           count_obs)
                   else
                      do m=1,LVT_rc%nensem
                         do t=1,LVT_rc%ngrid
                            do k=1,model%selectNlevs
                               do l=1,LVT_rc%strat_nlevels
                                  model_value_ts(t,k,l) = &
                                       model_value_ts(t,k,l) + &
                                       stats%maxtime(m)%model_value_ts(t,k,l)
                               enddo
                            enddo
                         enddo
                      enddo
                      
                      do t=1,LVT_rc%ngrid
                         do k=1,model%selectNlevs
                            do l=1,LVT_rc%strat_nlevels
                               model_value_ts(t,k,l) = &
                                    model_value_ts(t,k,l)/LVT_rc%nensem
                            enddo
                         enddo
                      enddo
                      
                      call LVT_writeTSinfo(metric%ftn_ts_loc(:,1),&
                           model,&
                           LVT_rc%ngrid,&
                           model_value_ts,&
                           count_model)
                   endif
                   deallocate(model_value_ts)
                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                      deallocate(obs_value_ts)
                   endif
                endif

             end if
             
          endif
       end if
    end if

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%maxtime(m)%model_value_total(t,k,l).gt.-1E10) then 
                         stats%maxtime(m)%model_value_total(t,k,l) = &
                              stats%maxtime(m)%model_value_total(t,k,l)
                      else
                         stats%maxtime(m)%model_value_total(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                   if(metric%computeSC.eq.1) then 
                      do l=1, LVT_rc%nasc
                         if(stats%maxtime(m)%model_value_asc(t,k,l).gt.-1E10) then 
                            stats%maxtime(m)%model_value_asc(t,k,l) = &
                                 count_model_asc(t,1,l)           
                         else
                            stats%maxtime(m)%model_value_asc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   
                   if(metric%computeADC.eq.1) then 
                      do l=1, LVT_rc%nadc
                         if(stats%maxtime(m)%model_value_adc(t,k,l).gt.-1E10) then 
                            stats%maxtime(m)%model_value_adc(t,k,l) = &
                                 stats%maxtime(m)%model_value_adc(t,k,l)
                         else
                            stats%maxtime(m)%model_value_adc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then
                         do l=1,LVT_rc%strat_nlevels                      
                            if(stats%maxtime(m)%obs_value_total(t,k,l).gt.-1E10) then 
                               stats%maxtime(m)%obs_value_total(t,k,l) = &
                                    stats%maxtime(m)%obs_value_total(t,k,l)
                            else
                               stats%maxtime(m)%obs_value_total(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                         if(metric%computeSC.eq.1) then
                            do l=1, LVT_rc%nasc 
                               if(stats%maxtime(m)%obs_value_asc(t,k,l).gt.-1E10) then 
                                  stats%maxtime(m)%obs_value_asc(t,k,l) = &
                                       stats%maxtime(m)%obs_value_asc(t,k,l)
                               else
                                  stats%maxtime(m)%obs_value_asc(t,k,l) = LVT_rc%udef
                               endif
                            enddo
                         endif
                         
                         if(metric%computeADC.eq.1) then
                            do l=1, LVT_rc%nadc  
                               if(stats%maxtime(m)%obs_value_adc(t,k,l).gt.-1E10) then
                                  stats%maxtime(m)%obs_value_adc(t,k,l) = &
                                       stats%maxtime(m)%obs_value_adc(t,k,l)
                               else
                                  stats%maxtime(m)%obs_value_adc(t,k,l) = LVT_rc%udef
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
                   call LVT_computeCI(stats%maxtime(m)%model_value_total(:,k,l),LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%maxtime(m)%model_value_ci(k,l))
                enddo
             enddo

             if(LVT_rc%obssource(2).ne."none") then 
                do k=1,obs%selectNlevs
                   do l=1, LVT_rc%strat_nlevels
                      call LVT_computeCI(stats%maxtime(m)%obs_value_total(:,k,l),LVT_rc%ngrid,&
                           LVT_rc%pval_CI,stats%maxtime(m)%obs_value_ci(k,l))
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
                              stats%maxtime(m)%model_value_total(t,k,l)
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
                                 stats%maxtime(m)%obs_value_total(t,k,l)
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
                         if(stats%maxtime(m)%model_value_asc(t,k,l).ne.LVT_rc%udef) then 
                            model_value_asc(t,k,l) = &
                                 model_value_asc(t,k,l) + & 
                                 stats%maxtime(m)%model_value_asc(t,k,l)
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
                                 stats%maxtime(m)%obs_value_asc(t,k,l)
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
                     count_model_asc,&
                     LVT_rc%ngrid,obs_value_asc,&
                     count_obs_asc)   

                deallocate(obs_value_asc)
       
             else
                call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_asc,&
                     count_model_asc)

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
                              stats%maxtime(m)%model_value_adc(t,k,l)
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
                                 stats%maxtime(m)%obs_value_adc(t,k,l)
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
                     count_model_adc,&
                     LVT_rc%ngrid,obs_value_adc,&
                     count_obs_adc)      
                
                deallocate(obs_value_adc)
             else
                call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_adc,&
                     count_model_adc)

                deallocate(model_value_adc)
             endif
          endif
       endif
    endif
  end subroutine computeSingleMaxtime

!BOP
! 
! !ROUTINE: LVT_writeMetric_Maxtime
! \label(LVT_writeMetric_Maxtime)
!
! !INTERFACE:
  subroutine LVT_writeMetric_Maxtime(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%maxtime%npass) then
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
                              stats%maxtime(m)%model_value_ts(:,k,l)
                         count_model_value_ts(:,m) = 1
                       if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                          obs_value_ts(:,m) = &
                               stats%maxtime(m)%obs_value_ts(:,k,l)
                          count_obs_value_ts(:,m) = 1
                       endif
                    enddo
                    if(LVT_metrics%maxtime%timeOpt.eq.1) then 
                       
                       call LVT_writevar_gridded(LVT_metrics%maxtime%ftn_ts, &
                            model_value_ts(:,:),&
                            stats%vid_ts(LVT_MAXTIMEid,1),k)
                       
                       call LVT_writevar_gridded(LVT_metrics%maxtime%ftn_ts, &
                            real(count_model_value_ts(:,:)),&
                            stats%vid_count_ts(LVT_MAXTIMEid,1),k)
                    endif
                    
                    if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                       call LVT_writevar_gridded(LVT_metrics%maxtime%ftn_ts, &
                            obs_value_ts(:,:),&
                            stats%vid_ts(LVT_MAXTIMEid,2),k)
                       call LVT_writevar_gridded(LVT_metrics%maxtime%ftn_ts, &
                            real(count_obs_value_ts(:,:)),&
                            stats%vid_count_ts(LVT_MAXTIMEid,2),k)
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
          if(pass.eq.LVT_metrics%maxtime%npass) then
             if(stats%selectOpt.eq.1) then

              allocate(model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
              allocate(count_model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
              allocate(model_value_ci(LVT_rc%nensem))

              if(LVT_metrics%maxtime%computeSC.eq.1) then 
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
              if(LVT_metrics%maxtime%computeADC.eq.1) then 
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
                              stats%maxtime(m)%model_value_total(:,k,l)
                         count_model_value_total(:,m) = 1
                         model_value_ci(m) = stats%maxtime(m)%model_value_ci(k,l)
                         
                         if(LVT_metrics%maxtime%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               model_value_asc(:,m,tind) = & 
                                    stats%maxtime(m)%model_value_asc(:,k,tind)
                            enddo
                         endif

                         if(LVT_metrics%maxtime%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               model_value_adc(:,m,tind) = & 
                                    stats%maxtime(m)%model_value_adc(:,k,tind)
                            enddo
                         endif

                         if(LVT_rc%obssource(2).ne."none".and.&
                              obs%selectNlevs.ge.1) then 
                            obs_value_total(:,m) = &
                                 stats%maxtime(m)%obs_value_total(:,k,l)
                            count_obs_value_total(:,m) = 1
                            obs_value_ci(m) = stats%maxtime(m)%obs_value_ci(k,l)

                            if(LVT_metrics%maxtime%computeSC.eq.1) then 
                               do tind = 1,LVT_rc%nasc
                                  obs_value_asc(:,m,tind) = & 
                                       stats%maxtime(m)%obs_value_asc(:,k,tind)
                               enddo
                            endif
                            if(LVT_metrics%maxtime%computeADC.eq.1) then 
                               do tind = 1,LVT_rc%nadc
                                  obs_value_adc(:,m,tind) = & 
                                       stats%maxtime(m)%obs_value_adc(:,k,tind)
                               enddo
                            endif
                            
                         endif
                      enddo
                      if(LVT_metrics%maxtime%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%maxtime%ftn_total, &
                              model_value_total(:,:),&
                              stats%vid_total(LVT_MAXTIMEid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%maxtime%ftn_total, &
                              real(count_model_value_total(:,:)),&
                              stats%vid_count_total(LVT_MAXTIMEid,1),k)

                         if(LVT_rc%obssource(2).ne."none".and.&
                              obs%selectNlevs.ge.1) then 
                            call LVT_writevar_gridded(LVT_metrics%maxtime%ftn_total, &
                                 obs_value_total(:,:),&
                                 stats%vid_total(LVT_MAXTIMEid,2),k)
                            call LVT_writevar_gridded(LVT_metrics%maxtime%ftn_total, &
                                 real(count_obs_value_total(:,:)),&
                                 stats%vid_count_total(LVT_MAXTIMEid,2),k)
                         endif
                         
                         if(LVT_metrics%maxtime%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%maxtime%ftn_total,&
                                    model_value_asc(:,:,tind),&
                                    stats%vid_sc_total(tind,LVT_MAXTIMEid,1),k)
                            enddo
                            if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                               do tind = 1,LVT_rc%nasc
                                  call LVT_writevar_gridded(&
                                       LVT_metrics%maxtime%ftn_total,&
                                       obs_value_asc(:,:,tind),&
                                       stats%vid_sc_total(tind,LVT_MAXTIMEid,2),k)
                               enddo
                            endif
                         endif
                         if(LVT_metrics%maxtime%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%maxtime%ftn_total,&
                                    model_value_adc(:,:,tind),&
                                    stats%vid_adc_total(tind,LVT_MAXTIMEid,1),k)
                            enddo
                            if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                               do tind = 1,LVT_rc%nadc
                                  call LVT_writevar_gridded(&
                                       LVT_metrics%maxtime%ftn_total,&
                                       obs_value_adc(:,:,tind),&
                                       stats%vid_adc_total(tind,LVT_MAXTIMEid,2),k)
                               enddo
                            endif
                         endif
                         call LVT_writeSummaryStats(&
                              LVT_metrics%maxtime%ftn_summ,&
                              l,&
                              LVT_metrics%maxtime%short_name,&
                              LVT_rc%ngrid,&
                              model_value_total(:,:), &
                              count_model_value_total(:,:),&
                              stats%standard_name,&
                              model_value_ci(:))
                         if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                            call LVT_writeSummaryStats(&
                                 LVT_metrics%maxtime%ftn_summ,&
                                 l,&
                                 LVT_metrics%maxtime%short_name,&
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

                if(LVT_metrics%maxtime%computeSC.eq.1) then 
                   deallocate(model_value_asc)
                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                      deallocate(obs_value_asc)
                   endif
                endif
                if(LVT_metrics%maxtime%computeADC.eq.1) then 
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
  end subroutine LVT_writeMetric_Maxtime

!BOP
! 
! !ROUTINE: LVT_resetMetric_Maxtime
! \label(LVT_resetMetric_Maxtime)
!
! !INTERFACE:
  subroutine LVT_resetMetric_Maxtime(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical          :: alarm
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
                if(LVT_metrics%maxtime%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      if(alarm) then 
                         stats%maxtime(m)%model_value_ts(:,k,l) = -1E10
                         stats%maxtime(m)%model_value_max_ts(:,k,l) = -1E10 
                      endif
                   enddo
                   if(LVT_rc%obssource(2).ne."none".and.&
                        obs%selectNlevs.ge.1) then 
                      do l=1,LVT_rc%strat_nlevels
                         if(alarm) then 
                            stats%maxtime(m)%obs_value_ts(:,k,l) = -1E10
                            stats%maxtime(m)%obs_value_max_ts(:,k,l) = -1E10
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

  end subroutine LVT_resetMetric_Maxtime

!BOP
! 
! !ROUTINE: LVT_writerestart_Maxtime
! 
! !INTERFACE:
  subroutine LVT_writerestart_Maxtime(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Maxtime metric computations
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
       if(LVT_metrics%maxtime%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_writevar_restart(ftn,&
                           stats%maxtime(m)%model_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%maxtime%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc         
                         call LVT_writevar_restart(ftn,&
                              stats%maxtime(m)%model_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%maxtime%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc         
                         call LVT_writevar_restart(ftn,&
                              stats%maxtime(m)%model_value_adc(:,k,l))
                      enddo
                   enddo
                endif
                
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            call LVT_writevar_restart(ftn,&
                                 stats%maxtime(m)%obs_value_total(:,k,l))
                         enddo
                      enddo
                      
                      if(LVT_metrics%maxtime%computeSC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_writevar_restart(ftn,&
                                    stats%maxtime(m)%obs_value_asc(:,k,l))
                            enddo
                         enddo
                      endif
                      if(LVT_metrics%maxtime%computeADC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_writevar_restart(ftn,&
                                    stats%maxtime(m)%obs_value_adc(:,k,l))
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
  end subroutine LVT_writerestart_Maxtime


!BOP
! 
! !ROUTINE: LVT_readrestart_Maxtime
! 
! !INTERFACE:
  subroutine LVT_readrestart_Maxtime(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Maxtime metric computations
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
       if(LVT_metrics%maxtime%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_readvar_restart(ftn,&
                           stats%maxtime(m)%model_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%maxtime%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc         
                         call LVT_readvar_restart(ftn,&
                              stats%maxtime(m)%model_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%maxtime%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc         
                         call LVT_readvar_restart(ftn,&
                              stats%maxtime(m)%model_value_adc(:,k,l))
                      enddo
                   enddo
                endif
                
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            call LVT_readvar_restart(ftn,&
                                 stats%maxtime(m)%obs_value_total(:,k,l))
                         enddo
                      enddo
                      
                      if(LVT_metrics%maxtime%computeSC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_readvar_restart(ftn,&
                                    stats%maxtime(m)%obs_value_asc(:,k,l))
                            enddo
                         enddo
                      endif
                      if(LVT_metrics%maxtime%computeADC.eq.1) then 
                         do k=1,obs%selectNlevs
                            do l=1,LVT_rc%nasc
                               call LVT_readvar_restart(ftn,&
                                    stats%maxtime(m)%obs_value_adc(:,k,l))
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
  end subroutine LVT_readrestart_Maxtime


end module LVT_MaxtimeMod
