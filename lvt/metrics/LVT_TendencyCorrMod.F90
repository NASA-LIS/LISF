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
! !MODULE: LVT_TendencyCorrMod
! \label(LVT_TendencyCorrMod)
!
! !INTERFACE:
module LVT_TendencyCorrMod
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
!  This module handles the computations required to compute the
!  correlation of tendencies from the values of the datastreams. 
!  Tendency is defined as the change in the values between 
!  consecutive averaging intervals. 
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
  public :: LVT_initTendencyCorr
  public :: LVT_diagnoseTendencyCorr
  public :: LVT_computeTendencyCorr
  public :: LVT_writeMetric_TendencyCorr
  public :: LVT_resetMetric_TendencyCorr
  public :: LVT_writerestart_TendencyCorr
  public :: LVT_readrestart_TendencyCorr
!EOP
  
  private

contains
  subroutine LVT_initTendencyCorr(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!EOP
    integer                 :: m 
    
    allocate(stats%tendencycorr(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%tendencycorr(m)%sxy_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr(m)%sxx_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr(m)%syy_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr(m)%sx_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr(m)%sy_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr(m)%rval_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr(m)%count_value(LVT_rc%ngrid, selectNlevs(1), LVT_rc%strat_nlevels))       
          stats%tendencycorr(m)%sxy_r = 0 
          stats%tendencycorr(m)%sxx_r = 0 
          stats%tendencycorr(m)%syy_r = 0 
          stats%tendencycorr(m)%sx_r = 0 
          stats%tendencycorr(m)%sy_r = 0 
          stats%tendencycorr(m)%rval_r = 0
          stats%tendencycorr(m)%count_value = 0 

          allocate(stats%tendencycorr(m)%value_ci(selectNlevs(1),LVT_rc%strat_nlevels))
          stats%tendencycorr(m)%value_ci = LVT_rc%udef

          if(metric%timeOpt.eq.1) then 
             allocate(stats%tendencycorr(m)%sxy_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%tendencycorr(m)%sxx_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%tendencycorr(m)%syy_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%tendencycorr(m)%sx_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%tendencycorr(m)%sy_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%tendencycorr(m)%rval_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%tendencycorr(m)%count_value_ts(LVT_rc%ngrid, selectNlevs(1), LVT_rc%strat_nlevels))       
             stats%tendencycorr(m)%sxy_ts_r = 0 
             stats%tendencycorr(m)%sxx_ts_r = 0 
             stats%tendencycorr(m)%syy_ts_r = 0 
             stats%tendencycorr(m)%sx_ts_r = 0 
             stats%tendencycorr(m)%sy_ts_r = 0 
             stats%tendencycorr(m)%rval_ts_r = 0
             stats%tendencycorr(m)%count_value_ts = 0 

             allocate(stats%tendencycorr(m)%tavg_value_ts(LVT_rc%ngrid, &
                  selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%tendencycorr(m)%tavg_count_value_ts(LVT_rc%ngrid, &
                  selectNlevs(1),LVT_rc%strat_nlevels))
             stats%tendencycorr(m)%tavg_value_ts = 0.0
             stats%tendencycorr(m)%tavg_count_value_ts=0 

          endif

          allocate(stats%tendencycorr(m)%model_value_cval_ts(LVT_rc%ngrid,selectNlevs(1), &
               LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr(m)%model_value_pval_ts(LVT_rc%ngrid,selectNlevs(1), &
               LVT_rc%strat_nlevels))

          stats%tendencycorr(m)%model_value_cval_ts = LVT_rc%udef
          stats%tendencycorr(m)%model_value_pval_ts = LVT_rc%udef

          allocate(stats%tendencycorr(m)%obs_value_cval_ts(LVT_rc%ngrid,selectNlevs(2), &
               LVT_rc%strat_nlevels))
          allocate(stats%tendencycorr(m)%obs_value_pval_ts(LVT_rc%ngrid,selectNlevs(2), &
               LVT_rc%strat_nlevels))

          stats%tendencycorr(m)%obs_value_cval_ts = LVT_rc%udef
          stats%tendencycorr(m)%obs_value_pval_ts = LVT_rc%udef

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

  end subroutine LVT_initTendencyCorr

!BOP
! 
! !ROUTINE: LVT_diagnoseTendencyCorr
! \label{LVT_diagnoseTendencyCorr}
!
! !INTERFACE: 
  subroutine LVT_diagnoseTendencyCorr(pass)
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
!   calculating the tendency of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelTendencyCorr](\ref{diagnoseSingleModelTendencyCorr})
!     updates the tendency computation for a single variable 
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
       if(LVT_metrics%tendencycorr%selectOpt.eq.1.or.&
            LVT_metrics%tendencycorr%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleTendencyCorr(model,obs,stats,&
                  LVT_metrics%tendencycorr)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseTendencyCorr

!BOP
! 
! !ROUTINE: diagnoseSingleTendencyCorr
! \label{diagnoseSingleTendencyCorr}
!
! !INTERFACE: 
  subroutine diagnoseSingleTendencyCorr(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the tendency computation of the 
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
    integer    :: t,k,m,tind,m_k,o_k
    real       :: m_diffval, o_diffval

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1.and.obs%selectNlevs.ge.1) then        
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(metric%selectOpt.eq.1) then
                   stats%tendencycorr(m)%model_value_pval_ts(t,k,1) = &
                        stats%tendencycorr(m)%model_value_cval_ts(t,k,1) 
                   stats%tendencycorr(m)%model_value_cval_ts(t,k,1) = &
                        model%value(t,m,m_k)
                   
                   stats%tendencycorr(m)%obs_value_pval_ts(t,k,1) = &
                        stats%tendencycorr(m)%obs_value_cval_ts(t,k,1) 
                   stats%tendencycorr(m)%obs_value_cval_ts(t,k,1) = &
                        obs%value(t,m,o_k)
                   
                   if(stats%tendencycorr(m)%model_value_cval_ts(t,k,1).ne.LVT_rc%udef.and.&
                        stats%tendencycorr(m)%model_value_pval_ts(t,k,1).ne.LVT_rc%udef.and.&
                        stats%tendencycorr(m)%obs_value_cval_ts(t,k,1).ne.LVT_rc%udef.and.&
                        stats%tendencycorr(m)%obs_value_pval_ts(t,k,1).ne.LVT_rc%udef) then
                      
                      m_diffval = & 
                           (stats%tendencycorr(m)%model_value_cval_ts(t,k,1) - &
                           stats%tendencycorr(m)%model_value_pval_ts(t,k,1))
                      
                      o_diffval = & 
                           (stats%tendencycorr(m)%obs_value_cval_ts(t,k,1) - &
                           stats%tendencycorr(m)%obs_value_pval_ts(t,k,1))
                      
                      stats%tendencycorr(m)%sxy_r(t,k,1) =  stats%tendencycorr(m)%sxy_r(t,k,1) + & 
                           m_diffval*o_diffval
                      
                      stats%tendencycorr(m)%sxx_r(t,k,1) =  stats%tendencycorr(m)%sxx_r(t,k,1) + & 
                           m_diffval*m_diffval
                      
                      stats%tendencycorr(m)%syy_r(t,k,1) =  stats%tendencycorr(m)%syy_r(t,k,1) + & 
                           o_diffval*o_diffval
                      
                      stats%tendencycorr(m)%sx_r(t,k,1) =  stats%tendencycorr(m)%sx_r(t,k,1) + & 
                           m_diffval
                      
                      stats%tendencycorr(m)%sy_r(t,k,1) =  stats%tendencycorr(m)%sy_r(t,k,1) + & 
                           o_diffval
                      
                      stats%tendencycorr(m)%count_value(t,k,1) = & 
                           stats%tendencycorr(m)%count_value(t,k,1)+1
                      
                      if(LVT_rc%strat_nlevels.gt.1) then
                         if(LVT_stats%strat_var(t,m,k).gt.&
                              LVT_rc%strat_var_threshold) then
                            m_diffval = & 
                                 (stats%tendencycorr(m)%model_value_cval_ts(t,k,2) - &
                                 stats%tendencycorr(m)%model_value_pval_ts(t,k,2))
                            
                            o_diffval = & 
                                 (stats%tendencycorr(m)%obs_value_cval_ts(t,k,2) - &
                                 stats%tendencycorr(m)%obs_value_pval_ts(t,k,2))
                            
                            stats%tendencycorr(m)%sxy_r(t,k,2) =  stats%tendencycorr(m)%sxy_r(t,k,2) + & 
                                 m_diffval*o_diffval
                            
                            stats%tendencycorr(m)%sxx_r(t,k,2) =  stats%tendencycorr(m)%sxx_r(t,k,2) + & 
                                 m_diffval*m_diffval
                            
                            stats%tendencycorr(m)%syy_r(t,k,2) =  stats%tendencycorr(m)%syy_r(t,k,2) + & 
                                 o_diffval*o_diffval
                            
                            stats%tendencycorr(m)%sx_r(t,k,2) =  stats%tendencycorr(m)%sx_r(t,k,2) + & 
                                 m_diffval
                            
                            stats%tendencycorr(m)%sy_r(t,k,2) =  stats%tendencycorr(m)%sy_r(t,k,2) + & 
                                 o_diffval
                            
                            stats%tendencycorr(m)%count_value(t,k,2) = & 
                                 stats%tendencycorr(m)%count_value(t,k,2)+1
                            
                         elseif(LVT_stats%strat_var(t,m,k).le.&
                              LVT_rc%strat_var_threshold) then
                            
                            m_diffval = & 
                                 (stats%tendencycorr(m)%model_value_cval_ts(t,k,3) - &
                                 stats%tendencycorr(m)%model_value_pval_ts(t,k,3))
                            
                            o_diffval = & 
                                 (stats%tendencycorr(m)%obs_value_cval_ts(t,k,3) - &
                                 stats%tendencycorr(m)%obs_value_pval_ts(t,k,3))
                            
                            stats%tendencycorr(m)%sxy_r(t,k,3) =  stats%tendencycorr(m)%sxy_r(t,k,3) + & 
                                 m_diffval*o_diffval
                            
                            stats%tendencycorr(m)%sxx_r(t,k,3) =  stats%tendencycorr(m)%sxx_r(t,k,3) + & 
                                 m_diffval*m_diffval
                            
                            stats%tendencycorr(m)%syy_r(t,k,3) =  stats%tendencycorr(m)%syy_r(t,k,3) + & 
                                 o_diffval*o_diffval
                            
                            stats%tendencycorr(m)%sx_r(t,k,3) =  stats%tendencycorr(m)%sx_r(t,k,3) + & 
                                 m_diffval
                            
                            stats%tendencycorr(m)%sy_r(t,k,3) =  stats%tendencycorr(m)%sy_r(t,k,3) + & 
                                 o_diffval
                            
                            stats%tendencycorr(m)%count_value(t,k,3) = & 
                                 stats%tendencycorr(m)%count_value(t,k,3)+1
                         endif
                      end if
                      if(metric%timeOpt.eq.1) then 
                         m_diffval = & 
                              (stats%tendencycorr(m)%model_value_cval_ts(t,k,1) - &
                              stats%tendencycorr(m)%model_value_pval_ts(t,k,1))
                         
                         o_diffval = & 
                              (stats%tendencycorr(m)%obs_value_cval_ts(t,k,1) - &
                              stats%tendencycorr(m)%obs_value_pval_ts(t,k,1))
                         
                         stats%tendencycorr(m)%sxy_ts_r(t,k,1) =  stats%tendencycorr(m)%sxy_ts_r(t,k,1) + & 
                              m_diffval*o_diffval
                         
                         stats%tendencycorr(m)%sxx_ts_r(t,k,1) =  stats%tendencycorr(m)%sxx_ts_r(t,k,1) + & 
                              m_diffval*m_diffval
                         
                         stats%tendencycorr(m)%syy_ts_r(t,k,1) =  stats%tendencycorr(m)%syy_ts_r(t,k,1) + & 
                              o_diffval*o_diffval
                         
                         stats%tendencycorr(m)%sx_ts_r(t,k,1) =  stats%tendencycorr(m)%sx_ts_r(t,k,1) + & 
                              m_diffval
                         
                         stats%tendencycorr(m)%sy_ts_r(t,k,1) =  stats%tendencycorr(m)%sy_ts_r(t,k,1) + & 
                              o_diffval
                         
                         stats%tendencycorr(m)%count_value_ts(t,k,1) = & 
                              stats%tendencycorr(m)%count_value_ts(t,k,1)+1
                         
                         if(LVT_rc%strat_nlevels.gt.1) then
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then
                               m_diffval = & 
                                    (stats%tendencycorr(m)%model_value_cval_ts(t,k,2) - &
                                    stats%tendencycorr(m)%model_value_pval_ts(t,k,2))
                               
                               o_diffval = & 
                                    (stats%tendencycorr(m)%obs_value_cval_ts(t,k,2) - &
                                    stats%tendencycorr(m)%obs_value_pval_ts(t,k,2))
                               
                               stats%tendencycorr(m)%sxy_ts_r(t,k,2) =  stats%tendencycorr(m)%sxy_ts_r(t,k,2) + & 
                                    m_diffval*o_diffval
                               
                               stats%tendencycorr(m)%sxx_ts_r(t,k,2) =  stats%tendencycorr(m)%sxx_ts_r(t,k,2) + & 
                                    m_diffval*m_diffval
                               
                               stats%tendencycorr(m)%syy_ts_r(t,k,2) =  stats%tendencycorr(m)%syy_ts_r(t,k,2) + & 
                                    o_diffval*o_diffval
                               
                               stats%tendencycorr(m)%sx_ts_r(t,k,2) =  stats%tendencycorr(m)%sx_ts_r(t,k,2) + & 
                                    m_diffval
                               
                               stats%tendencycorr(m)%sy_ts_r(t,k,2) =  stats%tendencycorr(m)%sy_ts_r(t,k,2) + & 
                                    o_diffval
                               
                               stats%tendencycorr(m)%count_value_ts(t,k,2) = & 
                                    stats%tendencycorr(m)%count_value_ts(t,k,2)+1
                               
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then
                               
                               m_diffval = & 
                                    (stats%tendencycorr(m)%model_value_cval_ts(t,k,3) - &
                                    stats%tendencycorr(m)%model_value_pval_ts(t,k,3))
                               
                               o_diffval = & 
                                    (stats%tendencycorr(m)%obs_value_cval_ts(t,k,3) - &
                                    stats%tendencycorr(m)%obs_value_pval_ts(t,k,3))
                               
                               stats%tendencycorr(m)%sxy_ts_r(t,k,3) =  stats%tendencycorr(m)%sxy_ts_r(t,k,3) + & 
                                    m_diffval*o_diffval
                               
                               stats%tendencycorr(m)%sxx_ts_r(t,k,3) =  stats%tendencycorr(m)%sxx_ts_r(t,k,3) + & 
                                    m_diffval*m_diffval
                               
                               stats%tendencycorr(m)%syy_ts_r(t,k,3) =  stats%tendencycorr(m)%syy_ts_r(t,k,3) + & 
                                    o_diffval*o_diffval
                               
                               stats%tendencycorr(m)%sx_ts_r(t,k,3) =  stats%tendencycorr(m)%sx_ts_r(t,k,3) + & 
                                    m_diffval
                               
                               stats%tendencycorr(m)%sy_ts_r(t,k,3) =  stats%tendencycorr(m)%sy_ts_r(t,k,3) + & 
                                    o_diffval
                               
                               stats%tendencycorr(m)%count_value_ts(t,k,3) = & 
                                    stats%tendencycorr(m)%count_value_ts(t,k,3)+1
                               
                            endif
                         endif
                      endif
                   end if
                end if
             end do
          enddo
       end do
    end if
  end subroutine diagnoseSingleTendencyCorr

!BOP
! 
! !ROUTINE: LVT_computeTendencyCorr
! \label{LVT_computeTendencyCorr}
!
! !INTERFACE: 
  subroutine LVT_computeTendencyCorr(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the tendency values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelTendencyCorr](\ref{computeSingleModelTendencyCorr})
!     computes the tendency values for a single variable
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
       if(LVT_metrics%tendencycorr%selectOpt.eq.1.or.&
            LVT_metrics%tendencycorr%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%tendencycorr%timeOpt.eq.1.and.&
                  LVT_metrics%tendencycorr%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%tendencycorr%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%tendencycorr%ftn_ts_loc(i,1),200,advance='no') &
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
             
             call computeSingleTendencyCorr(alarm,model,obs,stats,&
                  LVT_metrics%tendencycorr)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%tendencycorr%timeOpt.eq.1.and.&
                  LVT_metrics%tendencycorr%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%tendencycorr%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%tendencycorr%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif

             endif
          endif
       endif
    endif
  end subroutine LVT_computeTendencyCorr
  

!BOP
! 
! !ROUTINE: computeSingleTendencyCorr
! \label{computeSingleTendencyCorr}
!
! !INTERFACE: 
  subroutine computeSingleTendencyCorr(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the tendency values
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
    real     :: numer, denom, denom_sq

    real,    allocatable :: tavg_value_ts(:,:,:)
    real,    allocatable :: value_avg(:,:,:)

    if(metric%timeOpt.eq.1) then 
       if(alarm) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%tendencycorr(m)%count_value_ts(t,k,l).ne.0.and.&
                           stats%tendencycorr(m)%count_value_ts(t,k,l) &
                           .gt.LVT_rc%obsCountThreshold) then 
                      
                         numer = (float(stats%tendencycorr(m)%count_value_ts(t,k,l))* &
                              stats%tendencycorr(m)%sxy_ts_r(t,k,l) - &
                              stats%tendencycorr(m)%sx_ts_r(t,k,l)*stats%tendencycorr(m)%sy_ts_r(t,k,l))
                         denom_sq =  (float(stats%tendencycorr(m)%count_value_ts(t,k,l))* &
                              stats%tendencycorr(m)%sxx_ts_r(t,k,l)-&
                              stats%tendencycorr(m)%sx_ts_r(t,k,l)**2)* &
                              (float(stats%tendencycorr(m)%count_value_ts(t,k,l))*&
                              stats%tendencycorr(m)%syy_ts_r(t,k,l)-&
                              stats%tendencycorr(m)%sy_ts_r(t,k,l)**2)
                         
                         if(denom_sq.gt.0) then 
                            denom = sqrt(denom_sq)
                         else
                            denom = 0.0
                         endif
                         
                         if(denom.ne.0) then 
                            stats%tendencycorr(m)%rval_ts_r(t,k,l) = numer/denom
                            stats%tendencycorr(m)%tavg_value_ts(t,k,l) = & 
                                 stats%tendencycorr(m)%tavg_value_ts(t,k,l) + & 
                                 stats%tendencycorr(m)%rval_ts_r(t,k,l)
                            stats%tendencycorr(m)%tavg_count_value_ts(t,k,l) = & 
                                 stats%tendencycorr(m)%tavg_count_value_ts(t,k,l) + 1
                         else
                            stats%tendencycorr(m)%rval_ts_r(t,k,l) = LVT_rc%udef
                         endif
                         
                      endif
                   enddo
                enddo
             enddo
          enddo
          if(metric%extractTS.eq.1) then 
             if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                do m=1,LVT_rc%nensem
                   call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                        model,&
                        LVT_rc%ngrid,&
                        stats%tendencycorr(m)%tavg_value_ts,&
                        stats%tendencycorr(m)%tavg_count_value_ts)
                enddo
             else
                allocate(tavg_value_ts(LVT_rc%ngrid, &
                     model%vlevels, &
                     LVT_rc%strat_nlevels))
                tavg_value_ts = 0.0
                
                
                do m=1,LVT_rc%nensem
                   do t=1,LVT_rc%ngrid
                      do k=1,model%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            tavg_value_ts(t,k,l) = &
                                 tavg_value_ts(t,k,l) + &
                                 stats%tendencycorr(m)%tavg_value_ts(t,k,l)
                         enddo
                      enddo
                   enddo
                enddo
                
                do t=1,LVT_rc%ngrid
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         tavg_value_ts(t,k,l) = &
                              tavg_value_ts(t,k,l)/LVT_rc%nensem
                      enddo
                   enddo
                enddo
                
                call LVT_writeTSinfo(metric%ftn_ts_loc(:,1),&
                     model,&
                     LVT_rc%ngrid,&
                     tavg_value_ts,&
                     stats%tendencycorr(1)%tavg_count_value_ts)
                deallocate(tavg_value_ts)                   
             endif
          endif
          
       endif
    endif
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   if(stats%tendencycorr(m)%count_value(t,k,l).ne.0.and.&
                        stats%tendencycorr(m)%count_value(t,k,l) &
                        .gt.LVT_rc%obsCountThreshold) then 

                      numer = (float(stats%tendencycorr(m)%count_value(t,k,l))* &
                           stats%tendencycorr(m)%sxy_r(t,k,l) - &
                           stats%tendencycorr(m)%sx_r(t,k,l)*stats%tendencycorr(m)%sy_r(t,k,l))
                      denom_sq =  (float(stats%tendencycorr(m)%count_value(t,k,l))* &
                           stats%tendencycorr(m)%sxx_r(t,k,l)-&
                           stats%tendencycorr(m)%sx_r(t,k,l)**2)* &
                           (float(stats%tendencycorr(m)%count_value(t,k,l))*&
                           stats%tendencycorr(m)%syy_r(t,k,l)-&
                           stats%tendencycorr(m)%sy_r(t,k,l)**2)
                      if(denom_sq.gt.0) then 
                         denom = sqrt(denom_sq)
                      else
                         denom = 0.0
                      endif

                      if(denom.ne.0) then 
                         stats%tendencycorr(m)%rval_r(t,k,l) = numer/denom
                      else
                         stats%tendencycorr(m)%rval_r(t,k,l) = LVT_rc%udef
                         stats%tendencycorr(m)%count_value(t,k,l) = LVT_rc%udef
                      endif
                   else
                      stats%tendencycorr(m)%rval_r(t,k,l) = LVT_rc%udef
                      stats%tendencycorr(m)%count_value(t,k,l) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          do k=1,model%selectNlevs
             do l=1, LVT_rc%strat_nlevels
                call LVT_computeCI(stats%tendencycorr(m)%rval_r(:,k,l),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%tendencycorr(m)%value_ci(k,l))
             enddo
          enddo
       
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid,stats%tendencycorr(m)%rval_r)      
       endif
    endif

  end subroutine computeSingleTendencyCorr

!BOP
! 
! !ROUTINE: LVT_writeMetric_TendencyCorr
! \label(LVT_writeMetric_TendencyCorr)
!
! !INTERFACE:
  subroutine LVT_writeMetric_TendencyCorr(pass,final,vlevels,stats,obs)
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
!EOP
    integer                 :: k,l,m,tind

    real,    allocatable    :: value_total(:,:)
    integer, allocatable    :: count_value_total(:,:)
    real,    allocatable    :: value_ci(:)

    real,    allocatable    :: value_ts(:,:)
    integer, allocatable    :: count_value_ts(:,:)

    if(pass.eq.LVT_metrics%tendencycorr%npass) then
       if(final.ne.1) then
          if(stats%selectOpt.eq.1) then 

             allocate(value_ts(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_ts(LVT_rc%ngrid,LVT_rc%nensem))

             do k=1,vlevels
                do l=1,LVT_rc%strat_nlevels
                   do m=1,LVT_rc%nensem
                      value_ts(:,m) = &
                           stats%tendencycorr(m)%tavg_value_ts(:,k,l)
                      count_value_ts(:,m) = & 
                           stats%tendencycorr(m)%tavg_count_value_ts(:,k,l)
                   enddo
                   if(LVT_metrics%tendencycorr%timeOpt.eq.1) then 

                      call LVT_writevar_gridded(LVT_metrics%tendencycorr%ftn_ts, &
                           value_ts(:,:),&
                           stats%vid_ts(LVT_TENDENCYCORRid,1),k)

                      call LVT_writevar_gridded(LVT_metrics%tendencycorr%ftn_ts, &
                           real(count_value_ts(:,:)),&
                           stats%vid_count_ts(LVT_TENDENCYCORRid,1),k)
                   endif

                enddo
             enddo

             deallocate(value_ts)
             deallocate(count_value_ts)

          endif
       else
          if(pass.eq.LVT_metrics%tendencycorr%npass) then
             if(stats%selectOpt.eq.1) then

                allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(count_value_total(LVT_rc%ngrid,LVT_rc%nensem))

                allocate(value_ci(LVT_rc%nensem))

                do k=1,vlevels
                   do l=1,LVT_rc%strat_nlevels
                      do m=1,LVT_rc%nensem
                         value_total(:,m) = &
                              stats%tendencycorr(m)%rval_r(:,k,l)
                         count_value_total(:,m) = & 
                              stats%tendencycorr(m)%count_value(:,k,l)
                         value_ci(m) = stats%tendencycorr(m)%value_ci(k,l)

                      enddo
                      if(LVT_metrics%tendencycorr%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%tendencycorr%ftn_total, &
                              value_total(:,:),&
                              stats%vid_total(LVT_TENDENCYCORRid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%tendencycorr%ftn_total, &
                              real(count_value_total(:,:)),&
                              stats%vid_count_total(LVT_TENDENCYCORRid,1),k)

                         call LVT_writeSummaryStats(&
                              LVT_metrics%tendencycorr%ftn_summ,&
                              l,&
                              LVT_metrics%tendencycorr%short_name,&
                              LVT_rc%ngrid,&
                              value_total(:,:), &
                              count_value_total(:,:),&
                              stats%standard_name,&
                              value_ci(:))
                      endif
                   enddo
                enddo
                deallocate(value_total)
                deallocate(count_value_total)

             endif
          endif
       endif
    endif

  end subroutine LVT_writeMetric_TendencyCorr

!BOP
! 
! !ROUTINE: LVT_resetMetric_TendencyCorr
! \label(LVT_resetMetric_TendencyCorr)
!
! !INTERFACE:
  subroutine LVT_resetMetric_TendencyCorr(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical                :: alarm
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
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       
       if(stats%selectOpt.eq.1) then 
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                if(LVT_metrics%tendencycorr%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      if(alarm) then
                         stats%tendencycorr(m)%sxy_ts_r(:,k,l) = 0.0
                         stats%tendencycorr(m)%sxx_ts_r(:,k,l) = 0.0
                         stats%tendencycorr(m)%syy_ts_r(:,k,l) = 0.0
                         stats%tendencycorr(m)%sx_ts_r(:,k,l) = 0.0
                         stats%tendencycorr(m)%sy_ts_r(:,k,l) = 0.0
                         stats%tendencycorr(m)%rval_ts_r(:,k,l) = 0.0
                         stats%tendencycorr(m)%count_value_ts(:,k,l)=0 
                         
                         stats%tendencycorr(m)%tavg_value_ts(:,k,l) = 0.0
                         stats%tendencycorr(m)%tavg_count_value_ts(:,k,l)=0 
                      endif
                      
                   enddo
                endif
                
             enddo
          enddo
       endif
       model => model%next
       stats => stats%next

    enddo
  end subroutine LVT_resetMetric_TendencyCorr

!BOP
! 
! !ROUTINE: LVT_writerestart_TendencyCorr
! 
! !INTERFACE:
  subroutine LVT_writerestart_TendencyCorr(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for TendencyCorr metric computations
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

    if(LVT_metrics%tendencycorr%selectOpt.eq.1) then 
       
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)

       do while(associated(model))
          
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels     
                      call LVT_writevar_restart(ftn,stats%tendencycorr(m)%sxy_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%tendencycorr(m)%sxx_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%tendencycorr(m)%syy_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%tendencycorr(m)%sx_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%tendencycorr(m)%sy_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%tendencycorr(m)%rval_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%tendencycorr(m)%count_value(:,k,l))
                   enddo
                enddo
             enddo
          endif

          model => model%next
          obs   => obs%next
          stats => stats%next
          
       enddo
    end if
    
  end subroutine LVT_writerestart_TendencyCorr


!BOP
! 
! !ROUTINE: LVT_readrestart_TendencyCorr
! 
! !INTERFACE:
  subroutine LVT_readrestart_TendencyCorr(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for TendencyCorr metric computations
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

    if(LVT_metrics%tendencycorr%selectOpt.eq.1) then 
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_readvar_restart(ftn,stats%tendencycorr(m)%sxy_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%tendencycorr(m)%sxx_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%tendencycorr(m)%syy_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%tendencycorr(m)%sx_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%tendencycorr(m)%sy_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%tendencycorr(m)%rval_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%tendencycorr(m)%count_value(:,k,l))
                   enddo
                enddo
             enddo
          endif
          
          model => model%next
          obs   => obs%next
          stats => stats%next

       enddo
    end if

  end subroutine LVT_readrestart_TendencyCorr


end module LVT_TendencyCorrMod
