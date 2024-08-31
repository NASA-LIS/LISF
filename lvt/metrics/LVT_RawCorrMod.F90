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
! !MODULE: LVT_RawCorrMod
! \label(LVT_RawCorrMod)
!
! !INTERFACE:
module LVT_RawCorrMod
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

  private
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!  !DESCRIPTION: 
!   This module handles the raw correlation (pearson correlation
!   coefficient)  computations by comparing the values from 
!   the two datastreams
!
!   The metric calculations are supported across 
!   the entire analysis time period, stratified to a user defined
!   temporal averaging interval, average seasonal and diurnal cycle
!   and stratified to external datasets. 
!
! !FILES USED:
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initRawCorr
  public :: LVT_diagnoseRawCorr
  public :: LVT_computeRawCorr 
  public :: LVT_writeMetric_RawCorr
  public :: LVT_resetMetric_RawCorr
  public :: LVT_writerestart_RawCorr
  public :: LVT_readrestart_RawCorr

contains
  
!BOP
! 
! !ROUTINE: LVT_initRawCorr
! \label{LVT_initRawCorr}
!
! !INTERFACE: 
  subroutine LVT_initRawCorr(selectNlevs,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  This routine initializes the datastructures required to support
!  the pearson correlation computations. 
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

    allocate(stats%rcorr(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       
       if(metric%selectOpt.eq.1) then 
          allocate(stats%rcorr(m)%sxy_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%rcorr(m)%sxx_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%rcorr(m)%syy_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%rcorr(m)%sx_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%rcorr(m)%sy_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%rcorr(m)%rval_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
          allocate(stats%rcorr(m)%count_value(LVT_rc%ngrid, selectNlevs(1), LVT_rc%strat_nlevels))       
          stats%rcorr(m)%sxy_r = 0 
          stats%rcorr(m)%sxx_r = 0 
          stats%rcorr(m)%syy_r = 0 
          stats%rcorr(m)%sx_r = 0 
          stats%rcorr(m)%sy_r = 0 
          stats%rcorr(m)%rval_r = 0
          stats%rcorr(m)%count_value = 0 

          allocate(stats%rcorr(m)%value_ci(selectNlevs(1),LVT_rc%strat_nlevels))
          stats%rcorr(m)%value_ci = LVT_rc%udef

          if(metric%timeOpt.eq.1) then 
             allocate(stats%rcorr(m)%sxy_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%rcorr(m)%sxx_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%rcorr(m)%syy_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%rcorr(m)%sx_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%rcorr(m)%sy_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%rcorr(m)%rval_ts_r(LVT_rc%ngrid, selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%rcorr(m)%count_value_ts(LVT_rc%ngrid, selectNlevs(1), LVT_rc%strat_nlevels))       
             stats%rcorr(m)%sxy_ts_r = 0 
             stats%rcorr(m)%sxx_ts_r = 0 
             stats%rcorr(m)%syy_ts_r = 0 
             stats%rcorr(m)%sx_ts_r = 0 
             stats%rcorr(m)%sy_ts_r = 0 
             stats%rcorr(m)%rval_ts_r = 0
             stats%rcorr(m)%count_value_ts = 0 

             allocate(stats%rcorr(m)%tavg_value_ts(LVT_rc%ngrid, &
                  selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%rcorr(m)%tavg_count_value_ts(LVT_rc%ngrid, &
                  selectNlevs(1),LVT_rc%strat_nlevels))
             stats%rcorr(m)%tavg_value_ts = 0.0
             stats%rcorr(m)%tavg_count_value_ts=0 

             if(metric%computeSC.eq.1) then 
                allocate(stats%rcorr(m)%value_asc(LVT_rc%ngrid, selectNlevs(1), LVT_rc%nasc))
                stats%rcorr(m)%value_asc = 0.0
                allocate(stats%rcorr(m)%count_value_asc(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%nasc))
                stats%rcorr(m)%count_value_asc = 0
             endif

             if(metric%computeADC.eq.1) then 
                allocate(stats%rcorr(m)%value_adc(LVT_rc%ngrid, selectNlevs(1), LVT_rc%nadc))
                stats%rcorr(m)%value_adc = 0.0
                allocate(stats%rcorr(m)%count_value_adc(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%nadc))
                stats%rcorr(m)%count_value_adc = 0
             endif

          endif

       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initRawCorr
  
!BOP
! 
! !ROUTINE: LVT_diagnoseRawCorr
! \label{LVT_diagnoseRawCorr}
!
! !INTERFACE:  
  subroutine LVT_diagnoseRawCorr(pass)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the RawCorr calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleRawCorr](\ref{diagnoseSingleRawCorr})
!     updates the RawCorr computation for a single variable 
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
       if(LVT_metrics%rcorr%selectOpt.eq.1.or.&
            LVT_metrics%rcorr%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleRawCorr(obs,model,stats, &
                  LVT_metrics%rcorr)
             
             model => model%next
             obs   => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseRawCorr

!BOP
! 
! !ROUTINE: diagnoseSingleRawCorr
! \label{diagnoseSingleRawCorr}
!
! !INTERFACE: 
  subroutine diagnoseSingleRawCorr(obs, model, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the RawCorr computation (updates the running 
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

    integer    :: t,k,m,m_k,o_k

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                
                if(obs%count(t,m,o_k).ne.0.and. &
                     model%count(t,m,m_k).ne.0) then 
                   if(metric%selectOpt.eq.1) then
                      stats%rcorr(m)%sxy_r(t,k,1) = stats%rcorr(m)%sxy_r(t,k,1) + &
                           (model%value(t,m,m_k))*&
                           (obs%value(t,m,o_k))
                      stats%rcorr(m)%sx_r(t,k,1) = stats%rcorr(m)%sx_r(t,k,1) +&
                           (model%value(t,m,m_k))
                      stats%rcorr(m)%sy_r(t,k,1) = stats%rcorr(m)%sy_r(t,k,1) +&
                           (obs%value(t,m,o_k))
                      stats%rcorr(m)%sxx_r(t,k,1) = stats%rcorr(m)%sxx_r(t,k,1) +&
                           (model%value(t,m,m_k))*&
                           (model%value(t,m,m_k))
                      stats%rcorr(m)%syy_r(t,k,1) = stats%rcorr(m)%syy_r(t,k,1) +&
                           (obs%value(t,m,o_k))*&
                           (obs%value(t,m,o_k))    
                      stats%rcorr(m)%count_value(t,k,1) = &
                           stats%rcorr(m)%count_value(t,k,1) + 1
                      
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,m,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%rcorr(m)%sxy_r(t,k,2) = stats%rcorr(m)%sxy_r(t,k,2) + &
                                 (model%value(t,m,m_k))*&
                                 (obs%value(t,m,o_k))
                            stats%rcorr(m)%sx_r(t,k,2) = stats%rcorr(m)%sx_r(t,k,2) +&
                                 (model%value(t,m,m_k))
                            stats%rcorr(m)%sy_r(t,k,2) = stats%rcorr(m)%sy_r(t,k,2) +&
                                 (obs%value(t,m,o_k))
                            stats%rcorr(m)%sxx_r(t,k,2) = stats%rcorr(m)%sxx_r(t,k,2) +&
                                 (model%value(t,m,m_k))*&
                                 (model%value(t,m,m_k))
                            stats%rcorr(m)%syy_r(t,k,2) = stats%rcorr(m)%syy_r(t,k,2) +&
                                 (obs%value(t,m,o_k))*&
                                 (obs%value(t,m,o_k))               
                            stats%rcorr(m)%count_value(t,k,2) = stats%rcorr(m)%count_value(t,k,2) + 1

                         elseif(LVT_stats%strat_var(t,m,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%rcorr(m)%sxy_r(t,k,3) = stats%rcorr(m)%sxy_r(t,k,3) + &
                                 (model%value(t,m,m_k))*&
                                 (obs%value(t,m,o_k))
                            stats%rcorr(m)%sx_r(t,k,3) = stats%rcorr(m)%sx_r(t,k,3) +&
                                 (model%value(t,m,m_k))
                            stats%rcorr(m)%sy_r(t,k,3) = stats%rcorr(m)%sy_r(t,k,3) +&
                                 (obs%value(t,m,o_k))
                            stats%rcorr(m)%sxx_r(t,k,3) = stats%rcorr(m)%sxx_r(t,k,3) +&
                                 (model%value(t,m,m_k))*&
                                 (model%value(t,m,m_k))
                            stats%rcorr(m)%syy_r(t,k,3) = stats%rcorr(m)%syy_r(t,k,3) +&
                                 (obs%value(t,m,o_k))*&
                                 (obs%value(t,m,o_k))               
                            stats%rcorr(m)%count_value(t,k,3) = stats%rcorr(m)%count_value(t,k,3) + 1

                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         stats%rcorr(m)%sxy_ts_r(t,k,1) = stats%rcorr(m)%sxy_ts_r(t,k,1) + &
                              (model%value(t,m,m_k))*&
                              (obs%value(t,m,o_k))
                         stats%rcorr(m)%sx_ts_r(t,k,1) = stats%rcorr(m)%sx_ts_r(t,k,1) +&
                              (model%value(t,m,m_k))
                         stats%rcorr(m)%sy_ts_r(t,k,1) = stats%rcorr(m)%sy_ts_r(t,k,1) +&
                              (obs%value(t,m,o_k))
                         stats%rcorr(m)%sxx_ts_r(t,k,1) = stats%rcorr(m)%sxx_ts_r(t,k,1) +&
                              (model%value(t,m,m_k))*&
                              (model%value(t,m,m_k))
                         stats%rcorr(m)%syy_ts_r(t,k,1) = stats%rcorr(m)%syy_ts_r(t,k,1) +&
                              (obs%value(t,m,o_k))*&
                              (obs%value(t,m,o_k))     
                         stats%rcorr(m)%count_value_ts(t,k,1) = stats%rcorr(m)%count_value_ts(t,k,1) + 1
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%rcorr(m)%sxy_ts_r(t,k,2) = stats%rcorr(m)%sxy_ts_r(t,k,2) + &
                                    (model%value(t,m,m_k))*&
                                    (obs%value(t,m,o_k))
                               stats%rcorr(m)%sx_ts_r(t,k,2) = stats%rcorr(m)%sx_ts_r(t,k,2) +&
                                    (model%value(t,m,m_k))
                               stats%rcorr(m)%sy_ts_r(t,k,2) = stats%rcorr(m)%sy_ts_r(t,k,2) +&
                                    (obs%value(t,m,o_k))
                               stats%rcorr(m)%sxx_ts_r(t,k,2) = stats%rcorr(m)%sxx_ts_r(t,k,2) +&
                                    (model%value(t,m,m_k))*&
                                    (model%value(t,m,m_k))
                               stats%rcorr(m)%syy_ts_r(t,k,2) = stats%rcorr(m)%syy_ts_r(t,k,2) +&
                                    (obs%value(t,m,o_k))*&
                                    (obs%value(t,m,o_k))               
                               stats%rcorr(m)%count_value_ts(t,k,2) = &
                                    stats%rcorr(m)%count_value_ts(t,k,2) + 1

                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%rcorr(m)%sxy_ts_r(t,k,3) = stats%rcorr(m)%sxy_ts_r(t,k,3) + &
                                    (model%value(t,m,m_k))*&
                                    (obs%value(t,m,o_k))
                               stats%rcorr(m)%sx_ts_r(t,k,3) = stats%rcorr(m)%sx_ts_r(t,k,3) +&
                                    (model%value(t,m,m_k))
                               stats%rcorr(m)%sy_ts_r(t,k,3) = stats%rcorr(m)%sy_ts_r(t,k,3) +&
                                    (obs%value(t,m,o_k))
                               stats%rcorr(m)%sxx_ts_r(t,k,3) = stats%rcorr(m)%sxx_ts_r(t,k,3) +&
                                    (model%value(t,m,m_k))*&
                                    (model%value(t,m,m_k))
                               stats%rcorr(m)%syy_ts_r(t,k,3) = stats%rcorr(m)%syy_ts_r(t,k,3) +&
                                    (obs%value(t,m,o_k))*&
                                    (obs%value(t,m,o_k))               
                               stats%rcorr(m)%count_value_ts(t,k,3) = &
                                    stats%rcorr(m)%count_value_ts(t,k,3) + 1

                            endif
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleRawCorr


!BOP
! 
! !ROUTINE: LVT_computeRawCorr
! \label{LVT_computeRawCorr}
!
! !INTERFACE: 
  subroutine LVT_computeRawCorr(pass,alarm)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute RawCorr values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleRawCorr](\ref{computeSingleRawCorr})
!     updates the RawCorr computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     RawCorr computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: pass
    logical               :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    integer               :: i,m

    if(pass.eq.1) then 

       if(LVT_metrics%rcorr%selectOpt.eq.1.or.&
            LVT_metrics%rcorr%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%rcorr%timeOpt.eq.1.and.&
                  LVT_metrics%rcorr%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%rcorr%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%rcorr%ftn_ts_loc(i,1),200,advance='no') &
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
             call computeSingleRawCorr(alarm,obs,model,stats,&
                  LVT_metrics%rcorr)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo

          if(alarm) then 
             if(LVT_metrics%rcorr%timeOpt.eq.1.and.&
                  LVT_metrics%rcorr%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%rcorr%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%rcorr%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif

             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeRawCorr

!BOP
! 
! !ROUTINE: computeSingleRawCorr
! \label{computeSingleRawCorr}
!
! !INTERFACE: 
  subroutine computeSingleRawCorr(alarm,obs, model,stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the RawCorr values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     RawCorr computation has been reached
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

    logical                 :: alarm
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer  :: t,l,k,m,tind
    real     :: numer
    real     :: denom_sq
    real     :: denom

    real,    allocatable :: tavg_value_ts(:,:,:)
    real,    allocatable :: value_asc(:,:,:)
    integer, allocatable :: count_value_asc(:,:,:)

    real,    allocatable :: value_adc(:,:,:)
    real,    allocatable :: value_avg(:,:,:)


    if(metric%timeOpt.eq.1) then 
       if(alarm) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%rcorr(m)%count_value_ts(t,k,l).ne.0.and.&
                           stats%rcorr(m)%count_value_ts(t,k,l) &
                           .gt.LVT_rc%obsCountThreshold) then 
                         
                         numer = (float(stats%rcorr(m)%count_value_ts(t,k,l))* &
                              stats%rcorr(m)%sxy_ts_r(t,k,l) - &
                              stats%rcorr(m)%sx_ts_r(t,k,l)*stats%rcorr(m)%sy_ts_r(t,k,l))
                         denom_sq =  (float(stats%rcorr(m)%count_value_ts(t,k,l))* &
                              stats%rcorr(m)%sxx_ts_r(t,k,l)-&
                              stats%rcorr(m)%sx_ts_r(t,k,l)**2)* &
                              (float(stats%rcorr(m)%count_value_ts(t,k,l))*&
                              stats%rcorr(m)%syy_ts_r(t,k,l)-&
                              stats%rcorr(m)%sy_ts_r(t,k,l)**2)

                         if(denom_sq.gt.0) then 
                            denom = sqrt(denom_sq)
                         else
                            denom = 0.0
                         endif
                         if(denom.ne.0) then 
                            stats%rcorr(m)%rval_ts_r(t,k,l) = numer/denom
                            stats%rcorr(m)%tavg_value_ts(t,k,l) = & 
                                 stats%rcorr(m)%tavg_value_ts(t,k,l) + & 
                                 stats%rcorr(m)%rval_ts_r(t,k,l)
                            stats%rcorr(m)%tavg_count_value_ts(t,k,l) = & 
                                 stats%rcorr(m)%tavg_count_value_ts(t,k,l) + 1
                         else
                            stats%rcorr(m)%rval_ts_r(t,k,l) = LVT_rc%udef
                         endif
                         
                      endif
                   enddo
                   
                   if(metric%computeSC.eq.1) then 
                      if(stats%rcorr(m)%rval_ts_r(t,k,l).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                              tind)
                         stats%rcorr(m)%value_asc(t,k,tind) = & 
                              stats%rcorr(m)%value_asc(t,k,tind) + & 
                              stats%rcorr(m)%rval_ts_r(t,k,l)
                         stats%rcorr(m)%count_value_asc(t,k,tind) = & 
                              stats%rcorr(m)%count_value_asc(t,k,tind) + 1
                      endif
                   endif
                   
                   if(metric%computeADC.eq.1) then 
                      if(stats%rcorr(m)%rval_ts_r(t,k,l).ne.LVT_rc%udef) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%rcorr(m)%value_adc(t,k,tind) = & 
                              stats%rcorr(m)%value_adc(t,k,tind) + & 
                              stats%rcorr(m)%rval_ts_r(t,k,l)
                         stats%rcorr(m)%count_value_adc(t,k,tind) = & 
                              stats%rcorr(m)%count_value_adc(t,k,tind) + 1
                      endif
                   endif
                   
                enddo
             enddo
          enddo
             
          if(metric%extractTS.eq.1) then 
             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         if(stats%rcorr(m)%tavg_count_value_ts(t,k,l).gt.0) then 
                            stats%rcorr(m)%tavg_value_ts(t,k,l) = &
                                 stats%rcorr(m)%tavg_value_ts(t,k,l)/&
                                 stats%rcorr(m)%tavg_count_value_ts(t,k,l)
                         else
                            stats%rcorr(m)%tavg_value_ts(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   enddo
                enddo
             enddo

             if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                do m=1,LVT_rc%nensem
                   call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                        model,&
                        LVT_rc%ngrid,&
                        stats%rcorr(m)%tavg_value_ts,&
                        stats%rcorr(m)%tavg_count_value_ts)
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
                                 stats%rcorr(m)%tavg_value_ts(t,k,l)
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
                     stats%rcorr(1)%tavg_count_value_ts)
                deallocate(tavg_value_ts)                   
             endif

          endif
          
       endif
    endif
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%rcorr(m)%count_value(t,k,l).ne.0.and.&
                           stats%rcorr(m)%count_value(t,k,l) &
                           .gt.LVT_rc%obsCountThreshold) then 
                         
                         numer = (float(stats%rcorr(m)%count_value(t,k,l))* &
                              stats%rcorr(m)%sxy_r(t,k,l) - &
                              stats%rcorr(m)%sx_r(t,k,l)*stats%rcorr(m)%sy_r(t,k,l))
                         denom_sq =  (float(stats%rcorr(m)%count_value(t,k,l))* &
                              stats%rcorr(m)%sxx_r(t,k,l)-&
                              stats%rcorr(m)%sx_r(t,k,l)**2)* &
                              (float(stats%rcorr(m)%count_value(t,k,l))*&
                              stats%rcorr(m)%syy_r(t,k,l)-&
                              stats%rcorr(m)%sy_r(t,k,l)**2)
                         if(denom_sq.gt.0) then 
                            denom = sqrt(denom_sq)
                         else
                            denom = 0.0
                         endif
                         
                         if(denom.ne.0) then 
                            stats%rcorr(m)%rval_r(t,k,l) = numer/denom
                         else
                            stats%rcorr(m)%rval_r(t,k,l) = LVT_rc%udef
                            stats%rcorr(m)%count_value(t,k,l) = LVT_rc%udef
                         endif
                      else
                         stats%rcorr(m)%rval_r(t,k,l) = LVT_rc%udef
                         stats%rcorr(m)%count_value(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                   if(metric%computeSC.eq.1) then 
                      do l=1,LVT_rc%nasc
                         if(stats%rcorr(m)%count_value_asc(t,k,l).gt.&
                              LVT_rc%SCCountThreshold) then 
                            stats%rcorr(m)%value_asc(t,k,l) = &
                                 stats%rcorr(m)%value_asc(t,k,l)/&
                                 stats%rcorr(m)%count_value_asc(t,k,l) 
                         else
                            stats%rcorr(m)%value_asc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   if(metric%computeADC.eq.1) then 
                      do l=1,LVT_rc%nadc
                         if(stats%rcorr(m)%count_value_adc(t,k,l).gt.&
                              LVT_rc%ADCCountThreshold) then 
                            stats%rcorr(m)%value_adc(t,k,l) = &
                                 stats%rcorr(m)%value_adc(t,k,l)/&
                                 stats%rcorr(m)%count_value_adc(t,k,l) 
                         else
                            stats%rcorr(m)%value_adc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   

                enddo
             enddo
          enddo
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%rcorr(m)%rval_r(:,k,l),LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%rcorr(m)%value_ci(k,l))
                enddo
             enddo
          enddo

!----------------------------------------------------------------------
!  External data based stratification 
!----------------------------------------------------------------------
          if(LVT_rc%data_based_strat.eq.1) then 

             allocate(value_avg(LVT_rc%ngrid,&
                  model%vlevels,&
                  LVT_rc%strat_nlevels))
             value_avg = 0.0

             do m=1,LVT_rc%nensem        
                do t=1,LVT_rc%ngrid                      
                   do k=1,model%selectNlevs                
                      do l=1, LVT_rc%strat_nlevels
                         value_avg(t,k,l) = &
                              value_avg(t,k,l) + & 
                              stats%rcorr(m)%rval_r(t,k,l)
                      enddo
                   enddo
                enddo
             enddo
             
             do t=1,LVT_rc%ngrid                      
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%strat_nlevels
                      value_avg(t,k,l) = &
                           value_avg(t,k,l)/LVT_rc%nensem
                   enddo
                enddo
             enddo
                   
             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                  LVT_rc%ngrid,value_avg)            
             
             deallocate(value_avg)
          endif

          if(metric%computeSC.eq.1) then
             allocate(value_asc(LVT_rc%ngrid,&
                  model%vlevels,&
                  LVT_rc%nasc))
             allocate(count_value_asc(LVT_rc%ngrid,&
                  model%vlevels,&
                  LVT_rc%nasc))
             value_asc = 0.0
             count_value_asc = 0 

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%nasc
                      do m=1,LVT_rc%nensem        
                         if(stats%rcorr(m)%value_asc(t,k,l).ne.LVT_rc%udef) then 
                            value_asc(t,k,l) = &
                                 value_asc(t,k,l) + & 
                                 stats%rcorr(m)%value_asc(t,k,l)
                            count_value_asc(t,k,l) = & 
                                 count_value_asc(t,k,l) + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%nasc
                      if(count_value_asc(t,k,l).gt.0) then 
                         value_asc(t,k,l) = &
                              value_asc(t,k,l)/&
                              count_value_asc(t,k,l)
                      else
                         value_asc(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             call LVT_writeSeasonalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,value_asc,&
                  stats%rcorr(1)%count_value_asc)
             deallocate(value_asc)
             deallocate(count_value_asc)
            
          endif

          if(metric%computeADC.eq.1) then
             allocate(value_adc(LVT_rc%ngrid,&
                  model%vlevels,&
                  LVT_rc%nadc))
             value_adc = 0.0

             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%nadc
                      do m=1,LVT_rc%nensem          
                         value_adc(t,k,l) = &
                              value_adc(t,k,l) + & 
                              stats%rcorr(m)%value_adc(t,k,l)
                      enddo
                   enddo
                enddo
             enddo
             
             do t=1,LVT_rc%ngrid
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%nadc
                      value_adc(t,k,l) = &
                           value_adc(t,k,l)/LVT_rc%nensem
                   enddo
                enddo
             enddo
             call LVT_writeAvgDiurnalCycleInfo(model,obs,stats,metric,&
                  LVT_rc%ngrid,value_adc,&
                  stats%rcorr(m)%count_value_adc)
             
             deallocate(value_adc)
          endif


       endif
    endif

  end subroutine computeSingleRawCorr

!BOP
! 
! !ROUTINE: LVT_writeMetric_RawCorr
! \label{LVT_writeMetric_RawCorr}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_RawCorr(pass,final,vlevels,stats,obs)
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
    integer                 :: k,l,m,tind

    real,    allocatable    :: value_total(:,:)
    integer, allocatable    :: count_value_total(:,:)
    real,    allocatable    :: value_ci(:)

    real,    allocatable    :: value_ts(:,:)
    integer, allocatable    :: count_value_ts(:,:)

    real,    allocatable    :: value_adc(:,:,:)
    real,    allocatable    :: value_asc(:,:,:)

    if(pass.eq.LVT_metrics%rcorr%npass) then
       if(final.ne.1) then
          if(stats%selectOpt.eq.1) then 

             allocate(value_ts(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_ts(LVT_rc%ngrid,LVT_rc%nensem))

             do k=1,vlevels
                do l=1,LVT_rc%strat_nlevels
                   do m=1,LVT_rc%nensem
                      value_ts(:,m) = &
                           stats%rcorr(m)%tavg_value_ts(:,k,l)
                      count_value_ts(:,m) = & 
                           stats%rcorr(m)%tavg_count_value_ts(:,k,l)
                   enddo
                   if(LVT_metrics%rcorr%timeOpt.eq.1) then 

                      call LVT_writevar_gridded(LVT_metrics%rcorr%ftn_ts, &
                           value_ts(:,:),&
                           stats%vid_ts(LVT_RCORRid,1),k)

                      call LVT_writevar_gridded(LVT_metrics%rcorr%ftn_ts, &
                           real(count_value_ts(:,:)),&
                           stats%vid_count_ts(LVT_RCORRid,1),k)
                   endif

                enddo
             enddo

             deallocate(value_ts)
             deallocate(count_value_ts)

          endif
       else
          if(pass.eq.LVT_metrics%rcorr%npass) then
             if(stats%selectOpt.eq.1) then

                allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(count_value_total(LVT_rc%ngrid,LVT_rc%nensem))

                allocate(value_ci(LVT_rc%nensem))

                if(LVT_metrics%rcorr%computeSC.eq.1) then 
                   allocate(value_asc(LVT_rc%ngrid,&
                        LVT_rc%nensem,&
                        LVT_rc%nasc))
                endif
                if(LVT_metrics%rcorr%computeADC.eq.1) then 
                   allocate(value_adc(LVT_rc%ngrid,&
                        LVT_rc%nensem,&
                        LVT_rc%nadc))        

                endif

                do k=1,vlevels
                   do l=1,LVT_rc%strat_nlevels
                      do m=1,LVT_rc%nensem
                         value_total(:,m) = &
                              stats%rcorr(m)%rval_r(:,k,l)
                         count_value_total(:,m) = & 
                              stats%rcorr(m)%count_value(:,k,l)
                         value_ci(m) = stats%rcorr(m)%value_ci(k,l)

                         if(LVT_metrics%rcorr%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               value_asc(:,m,tind) = & 
                                    stats%rcorr(m)%value_asc(:,k,tind)
                            enddo
                         endif

                         if(LVT_metrics%rcorr%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               value_adc(:,m,tind) = & 
                                    stats%rcorr(m)%value_adc(:,k,tind)
                            enddo
                         endif
                      enddo
                      if(LVT_metrics%rcorr%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%rcorr%ftn_total, &
                              value_total(:,:),&
                              stats%vid_total(LVT_RCORRid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%rcorr%ftn_total, &
                              real(count_value_total(:,:)),&
                              stats%vid_count_total(LVT_RCORRid,1),k)

                         if(LVT_metrics%rcorr%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%rcorr%ftn_total,&
                                    value_asc(:,:,tind),&
                                    stats%vid_sc_total(tind,LVT_RCORRid,1),k)
                            enddo
                         endif
                         if(LVT_metrics%rcorr%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%rcorr%ftn_total,&
                                    value_adc(:,:,tind),&
                                    stats%vid_adc_total(tind,LVT_RCORRid,1),k)
                            enddo
                         endif
                         call LVT_writeSummaryStats(&
                              LVT_metrics%rcorr%ftn_summ,&
                              l,&
                              LVT_metrics%rcorr%short_name,&
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

                if(LVT_metrics%rcorr%computeSC.eq.1) then 
                   deallocate(value_asc)
                endif
                if(LVT_metrics%rcorr%computeADC.eq.1) then 
                   deallocate(value_adc)
                endif

             endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_RawCorr


!BOP
! 
! !ROUTINE: LVT_resetMetric_RawCorr
! \label(LVT_resetMetric_RawCorr)
!
! !INTERFACE:
  subroutine LVT_resetMetric_RawCorr(alarm)

! !INPUT PARAMETERS: 
    logical           :: alarm
 

!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

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
                if(LVT_metrics%rcorr%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      
                      if(alarm) then 
                         stats%rcorr(m)%sxy_ts_r(:,k,l) = 0.0
                         stats%rcorr(m)%sxx_ts_r(:,k,l) = 0.0
                         stats%rcorr(m)%syy_ts_r(:,k,l) = 0.0
                         stats%rcorr(m)%sx_ts_r(:,k,l) = 0.0
                         stats%rcorr(m)%sy_ts_r(:,k,l) = 0.0
                         stats%rcorr(m)%rval_ts_r(:,k,l) = 0.0
                         stats%rcorr(m)%count_value_ts(:,k,l)=0 
                         
                         stats%rcorr(m)%tavg_value_ts(:,k,l) = 0.0
                         stats%rcorr(m)%tavg_count_value_ts(:,k,l)=0 
                      endif
                      
                   enddo
                endif
                
             enddo
          enddo
       endif
       model => model%next
       stats => stats%next

    enddo
  end subroutine LVT_resetMetric_RawCorr


!BOP
! 
! !ROUTINE: LVT_writerestart_RawCorr
! 
! !INTERFACE:
  subroutine LVT_writerestart_RawCorr(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for RawCorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,l,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_metrics%rcorr%selectOpt.eq.1) then 
       
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)

       do while(associated(model))
          
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels     
                      call LVT_writevar_restart(ftn,stats%rcorr(m)%sxy_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%rcorr(m)%sxx_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%rcorr(m)%syy_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%rcorr(m)%sx_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%rcorr(m)%sy_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%rcorr(m)%rval_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%rcorr(m)%count_value(:,k,l))
                   enddo
                enddo
                
                if(LVT_metrics%rcorr%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_writevar_restart(ftn,&
                              stats%rcorr(m)%value_asc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%rcorr(m)%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%rcorr%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_writevar_restart(ftn,&
                              stats%rcorr(m)%value_adc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%rcorr(m)%count_value_adc(:,k,l))
                      enddo
                   enddo
                endif
             enddo
          endif

          model => model%next
          obs   => obs%next
          stats => stats%next
          
       enddo
    end if
    
  end subroutine LVT_writerestart_RawCorr

!BOP
! 
! !ROUTINE: LVT_readrestart_RawCorr
! 
! !INTERFACE:
  subroutine LVT_readrestart_RawCorr(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for RawCorr metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    integer              :: k,l,m

    if(LVT_metrics%rcorr%selectOpt.eq.1) then 
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_readvar_restart(ftn,stats%rcorr(m)%sxy_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%rcorr(m)%sxx_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%rcorr(m)%syy_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%rcorr(m)%sx_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%rcorr(m)%sy_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%rcorr(m)%rval_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%rcorr(m)%count_value(:,k,l))
                   enddo
                enddo
                
                if(LVT_metrics%rcorr%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_readvar_restart(ftn,&
                              stats%rcorr(m)%value_asc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%rcorr(m)%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%rcorr%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_readvar_restart(ftn,&
                              stats%rcorr(m)%value_adc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%rcorr(m)%count_value_adc(:,k,l))
                      enddo
                   enddo
                endif
             enddo
          endif
          
          model => model%next
          obs   => obs%next
          stats => stats%next

       enddo
    end if

  end subroutine LVT_readrestart_RawCorr
end module LVT_RawCorrMod
