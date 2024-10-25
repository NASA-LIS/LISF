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
! !MODULE: LVT_HNMod
! \label(LVT_HNMod)
!
! !INTERFACE:
module LVT_HNMod
! 
! !USES:   
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
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
!  This module handles the Hausdorff Norm (HN)
!  computations by comparing 
!  the LIS output to the specified observations. 
! 
! HS = max (min(d(A,B)))
! 
! !FILES USED:
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

  private

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initHN
  public :: LVT_diagnoseHN
  public :: LVT_computeHN 
  public :: LVT_writeMetric_HN
  public :: LVT_resetMetric_HN
  public :: LVT_writerestart_HN
  public :: LVT_readrestart_HN
contains
  
!BOP
! 
! !ROUTINE: LVT_initHN
! \label{LVT_initHN}
!
! !INTERFACE: 
  subroutine LVT_initHN(selectNlevs,stats,metric)
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
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP
    integer                 :: m

    allocate(stats%hn(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       
       if(metric%selectOpt.eq.1.or.metric%timeopt.eq.1) then 
          allocate(stats%hn(m)%value_model_total(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))
          allocate(stats%hn(m)%count_value_model_total(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))   
          stats%hn(m)%value_model_total = 0.0
          stats%hn(m)%count_value_model_total = 0

          allocate(stats%hn(m)%value_obs_total(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))
          allocate(stats%hn(m)%count_value_obs_total(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))   
          stats%hn(m)%value_obs_total = 0.0
          stats%hn(m)%count_value_obs_total = 0

          allocate(stats%hn(m)%value_ts(selectNlevs(1), &
               LVT_rc%strat_nlevels))
          allocate(stats%hn(m)%count_value_ts(selectNlevs(1),&
               LVT_rc%strat_nlevels))       

          stats%hn(m)%value_ts = 0.0
          stats%hn(m)%count_value_ts = 0 

          allocate(stats%hn(m)%tavg_value_ts(selectNlevs(1), &
               LVT_rc%strat_nlevels))
          allocate(stats%hn(m)%tavg_count_value_ts(selectNlevs(1),&
               LVT_rc%strat_nlevels))       

          stats%hn(m)%tavg_value_ts = 0.0
          stats%hn(m)%tavg_count_value_ts = 0 

       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initHN
  
!BOP
! 
! !ROUTINE: LVT_diagnoseHN
! \label{LVT_diagnoseHN}
!
! !INTERFACE: 
  subroutine LVT_diagnoseHN(pass)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the HN calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleHN](\ref{diagnoseSingleHN})
!     updates the HN computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: pass
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%hn%selectOpt.eq.1.or.&
            LVT_metrics%hn%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleHN(obs, model, stats,&
                  LVT_metrics%hn)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseHN

!BOP
! 
! !ROUTINE: diagnoseSingleHN
! \label{diagnoseSingleHN}
!
! !INTERFACE: 
  subroutine diagnoseSingleHN(obs, model, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the HN computation (updates the running 
!  sum calculations of the squared error) 
!  The arguments are: 
!
!  \begin{description}
!   \item[obs] observation object
!   \item[model] model variable object
!   \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,l,m,m_k,o_k,t1
    real       :: hnval(model%selectNlevs,LVT_rc%strat_nlevels)
    real       :: hn_minv 
    real       :: dist
    real       :: minval 

    minval = -100000.0

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       if(metric%selectOpt.eq.1) then          
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                if(trim(obs%units).eq.trim(model%units)) then 
                   do l=1,LVT_rc%strat_nlevels
                      do t=1,LVT_rc%ngrid
                         m_k = k+model%startNlevs -1
                         o_k = k+obs%startNlevs -1
                         if(model%value(t,m,m_k).ne.0) then 
                            !for final comp
                            stats%hn(m)%value_model_total(t,k,l) = &
                                 stats%hn(m)%value_model_total(t,k,l) + model%value(t,m,m_k)
                            stats%hn(m)%count_value_model_total(t,k,l) = & 
                                 stats%hn(m)%count_value_model_total(t,k,l) + 1
                         endif
                         if(obs%value(t,m,o_k).ne.0) then 
                            stats%hn(m)%value_obs_total(t,k,l) = &
                                 stats%hn(m)%value_obs_total(t,k,l) + obs%value(t,m,o_k)
                            stats%hn(m)%count_value_obs_total(t,k,l) = & 
                                 stats%hn(m)%count_value_obs_total(t,k,l) + 1
                         endif
                      enddo
                   enddo
                endif
             enddo
          enddo
       endif
       if(metric%timeOpt.eq.1) then         
          do m=1,LVT_rc%nensem 
             do k=1,model%selectNlevs
                if(trim(obs%units).eq.trim(model%units)) then 
                   hnval(k,1) = minval
                   do t=1,LVT_rc%ngrid
                      if(obs%value(t,m,o_k).ne.0.and. &
                           model%value(t,m,m_k).ne.0) then  
                         hn_minv = 10000.0 
                         do t1=1,LVT_rc%ngrid
                            if(obs%count(t1,m,k).ne.0) then
                               dist = sqrt((model%value(t,m,m_k)-obs%value(t1,m,k))**2)
                               if(dist.lt.hn_minv) then 
                                  hn_minv= dist
                               endif
                            endif
                         enddo
                         if(hn_minv.gt.hnval(k,1)) then 
                            hnval(k,1) = hn_minv 
                         endif
                      endif
                   enddo
                   if(hnval(k,1).eq.minval) then 
                      hnval(k,1) = LVT_rc%udef
                   else
                      stats%hn(m)%count_value_ts(k,1) = stats%hn(m)%count_value_ts(k,1)+1
                   endif
                   
                   stats%hn(m)%value_ts(k,1) = hnval(k,1)
                   
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,m,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         
                         hnval(k,2) = minval
                         do t=1,LVT_rc%ngrid
                            if(obs%value(t,m,o_k).ne.0.and. &
                                 model%value(t,m,m_k).ne.0) then                    
                               hn_minv = 10000.0 
                               do t1=1,LVT_rc%ngrid
                                  if(obs%count(t1,m,k).ne.0) then 
                                     dist = sqrt((model%value(t,m,m_k)-&
                                          obs%value(t1,m,k))**2)
                                     if(dist.lt.hn_minv) then 
                                        hn_minv= dist
                                     endif
                                  endif
                               enddo
                               if(hn_minv.gt.hnval(k,1)) then 
                                  hnval(k,2) = hn_minv                      
                               endif
                            endif
                         enddo
                         if(hnval(k,2).eq.minval) then 
                            hnval(k,2) = LVT_rc%udef
                         else
                            stats%hn(m)%count_value_ts(k,2) = stats%hn(m)%count_value_ts(k,2)+1
                         endif
                         
                         stats%hn(m)%value_ts(k,2) = hnval(k,2)
                         
                      elseif(LVT_stats%strat_var(t,m,k).le.&
                           LVT_rc%strat_var_threshold) then
                         
                         hnval(k,3) = minval
                         do t=1,LVT_rc%ngrid
                            if(obs%value(t,m,o_k).ne.0.and. &
                                 model%value(t,m,m_k).ne.0) then                    
                               hn_minv = 10000.0 
                               do t1=1,LVT_rc%ngrid
                                  if(obs%count(t1,m,k).ne.0) then 
                                     dist = sqrt((model%value(t,m,m_k)-&
                                          obs%value(t1,m,k))**2)
                                     if(dist.lt.hn_minv) then 
                                        hn_minv= dist
                                     endif
                                  endif
                               enddo
                               if(hn_minv.gt.hnval(k,1)) then 
                                  hnval(k,3) = hn_minv                      
                               endif
                            endif
                         enddo
                         if(hnval(k,3).eq.minval) then 
                            hnval(k,3) = LVT_rc%udef
                         else
                            stats%hn(m)%count_value_ts(k,3) = stats%hn(m)%count_value_ts(k,3)+1
                         endif
                         stats%hn(m)%value_ts(k,3) = hnval(k,3)
                      endif
                   endif
                else
                   write(LVT_logunit,*) 'For variable ',trim(model%standard_name)
                   write(LVT_logunit,*) 'observations are in ',trim(obs%units)
                   write(LVT_logunit,*) 'and LIS output is in ',trim(model%units)
                   write(LVT_logunit,*) 'please add the support of ',&
                        trim(model%units), ' in the observation plugin'
                   call LVT_endrun 
                endif
             enddo
          enddo
       endif
    end if
  end subroutine diagnoseSingleHN

!BOP
! 
! !ROUTINE: LVT_computeHN
! \label{LVT_computeHN}
!
! !INTERFACE: 
  subroutine LVT_computeHN(pass,alarm)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute HN values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleHN](\ref{computeSingleHN})
!     updates the HN computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     HN computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: pass
    logical     :: alarm
    real        :: dist
    integer     :: i,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%hn%selectOpt.eq.1.or.LVT_metrics%hn%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%hn%timeOpt.eq.1.and.&
                  LVT_metrics%hn%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%bias%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%bias%ftn_ts_loc(i,1),200,advance='no') &
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

             call computeSingleHN(alarm,obs, model, stats,&
                  LVT_metrics%hn)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       
          if(alarm) then 
             if(LVT_metrics%hn%timeOpt.eq.1.and.&
                  LVT_metrics%hn%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%bias%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%bias%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif

             endif
          endif
       end if
    end if
  end subroutine LVT_ComputeHN

!BOP
! 
! !ROUTINE: computeSingleHN
! \label{computeSingleHN}
!
! !INTERFACE: 
  subroutine computeSingleHN(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the HN values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     HN computation has been reached
!    \item[obs] observation object
!    \item[model] model variable object
!    \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none
    logical                 :: alarm
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer  :: t,l,k,m,tind,t1
    real       :: hnval(model%selectNlevs,LVT_rc%strat_nlevels)
    real       :: minval 
    real       :: dist, hn_minv

    minval = -100000.0

    if((metric%selectOpt.eq.1.or.metric%timeOpt.eq.1)) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%hn(m)%count_value_ts(k,l).ne.0) then 
                      stats%hn(m)%value_ts(k,l) = stats%hn(m)%value_ts(k,l)/&
                           (stats%hn(m)%count_value_ts(k,l))
                      stats%hn(m)%tavg_value_ts(k,l) = stats%hn(m)%tavg_value_ts(k,l)+ & 
                           stats%hn(m)%value_ts(k,l)
                      stats%hn(m)%tavg_count_value_ts(k,l) = &
                           stats%hn(m)%tavg_count_value_ts(k,l)+1
                   else
                      stats%hn(m)%value_ts(k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
             if(alarm) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%hn(m)%tavg_count_value_ts(k,l).ne.0) then 
                         stats%hn(m)%tavg_value_ts(k,l) = stats%hn(m)%tavg_value_ts(k,l)/&
                              (stats%hn(m)%tavg_count_value_ts(k,l))
                      endif
                   enddo
                enddo
                
                if(metric%extractTS.eq.1) then 
                   do k=1,model%selectNlevs
                      write(LVT_metrics%hn%ftn_ts_loc(1,m),201,advance='no') &
                           stats%hn(m)%tavg_value_ts(k,1) 
                   enddo
                endif
             endif
          enddo
       endif
    endif

201 format (E14.6)

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels                      
                   if(trim(obs%units).eq.trim(model%units)) then 
                      hnval(k,l) = minval
                      do t=1,LVT_rc%ngrid
                         if(stats%hn(m)%count_value_obs_total(t,k,l).ne.0.and. &
                              stats%hn(m)%count_value_model_total(t,k,l).ne.0) then 
                            stats%hn(m)%value_obs_total(t,k,l) = &
                                 stats%hn(m)%value_obs_total(t,k,l)/&
                                 stats%hn(m)%count_value_obs_total(t,k,l)
                            stats%hn(m)%value_model_total(t,k,l) = &
                                 stats%hn(m)%value_model_total(t,k,l)/&
                                 stats%hn(m)%count_value_model_total(t,k,l)
                         endif
                         
                         hn_minv = 10000.0 
                         do t1=1,LVT_rc%ngrid
                            if(stats%hn(m)%count_value_obs_total(t1,k,l).ne.0) then 
                               dist = sqrt((stats%hn(m)%value_model_total(t,k,l)-&
                                    stats%hn(m)%value_obs_total(t1,k,l))**2)
                               if(dist.lt.hn_minv) then 
                                  hn_minv= dist
                               endif
                            endif
                         enddo
                         if(hn_minv.gt.hnval(k,l)) then 
                            hnval(k,l) = hn_minv                      
                         endif
                      enddo
                   endif
                   if(hnval(k,l).eq.minval) then 
                      hnval(k,l) = LVT_rc%udef
                   endif
                   
                   stats%hn(m)%value_total(k,l) = hnval(k,l)
                enddo
             enddo

!          do k=1,model%selectNlevs
!             do l=1, LVT_rc%strat_nlevels
!                call LVT_computeCI(stats%hn(m)%value_total(:,k,l),LVT_rc%ngrid,&
!                     LVT_rc%pval_CI,stats%hn(m)%value_ci(k,l))
!             enddo
!          enddo
               
          enddo
       endif
    endif

  end subroutine computeSingleHN

!BOP
! 
! !ROUTINE: LVT_writeMetric_HN
! \label{LVT_writeMetric_HN}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_HN(pass,final,vlevels,stats,obs)
! 
! !USES:   
    use LVT_statsMod, only : LVT_writeSummaryStats
    use LVT_pluginIndices
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
! !ARGUMENTS: 
    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)   :: stats
    type(LVT_metaDataEntry) :: obs
!EOP
    integer                 :: l,m,tind
    integer                 :: k

    if(pass.eq.LVT_metrics%hn%npass) then 
       if(final.ne.1) then 
          if(LVT_metrics%hn%selectOpt.eq.1) then
             if(stats%selectOpt.eq.1) then 
                do m=1,LVT_rc%nensem
                   do k=1,vlevels
                      do l=1,LVT_rc%strat_nlevels
                         write(LVT_metrics%hn%ftn_summ,*) & 
                              k,l,stats%hn(m)%value_total(k,l)
                      enddo
                   enddo
                enddo
             endif
          endif
       end if
    endif
  end subroutine LVT_writeMetric_HN


!BOP
! 
! !ROUTINE: LVT_resetMetric_HN
! \label(LVT_resetMetric_HN)
!
! !INTERFACE:
  subroutine LVT_resetMetric_HN(alarm)
! 
! !INPUT PARAMETERS: 
    logical             :: alarm
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

    integer                 :: vlevels
    integer                 :: i,k,l,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_statsEntry)   , pointer :: stats


    call LVT_getDataStream1Ptr(model)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                if(LVT_metrics%hn%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%hn(m)%value_ts(k,l) = 0.0
                      stats%hn(m)%count_value_ts(k,l)=0 
                      if(alarm) then 
                         stats%hn(m)%tavg_value_ts(k,l) = 0.0
                         stats%hn(m)%tavg_count_value_ts(k,l)=0 
                      endif
                   enddo
                endif
             enddo
          enddo
       endif
       model => model%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_HN

!BOP
! 
! !ROUTINE: LVT_writerestart_HN
! 
! !INTERFACE:
  subroutine LVT_writerestart_HN(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for HN metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%hn%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for HN'
       stop
    end if
    
  end subroutine LVT_writerestart_HN

!BOP
! 
! !ROUTINE: LVT_readrestart_HN
! 
! !INTERFACE:
  subroutine LVT_readrestart_HN(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for HN metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%hn%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for HN'
       stop
    end if
    
  end subroutine LVT_readrestart_HN

end module LVT_HNMod
