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
! !MODULE: LVT_RnkCorrMod
! \label(LVT_RnkCorrMod)
!
! !INTERFACE:
module LVT_RnkCorrMod
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod
  use LVT_timeMgrMod

  implicit none

  private
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!  !DESCRIPTION: 
!   This module handles the rank correlation computations by comparing 
!   the values from the two datastreams. 
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
  public :: LVT_initRnkCorr
  public :: LVT_diagnoseRnkCorr
  public :: LVT_computeRnkCorr 
  public :: LVT_writeMetric_RnkCorr
  public :: LVT_resetMetric_RnkCorr
  public :: LVT_writerestart_RnkCorr
  public :: LVT_readrestart_RnkCorr

 type, public :: rnkcorrdec
    integer          :: nsize_total
    integer          :: nsize_season
    integer          :: tscale
 end type rnkcorrdec

  type(rnkcorrdec) :: LVT_rnkcorr_struc
contains
  
!BOP
! 
! !ROUTINE: LVT_initRnkCorr
! \label{LVT_initRnkCorr}
!
! !INTERFACE: 
  subroutine LVT_initRnkCorr(selectNlevs,stats,metric)
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
!  the rank correlation computations. 
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

    allocate(stats%rnkcorr(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          call computeAnalysisLength(LVT_rnkcorr_struc%nsize_total, &
               LVT_rc%tavgInterval)

          allocate(stats%rnkcorr(m)%model_value_final(LVT_rc%ngrid, &
               selectNlevs(1),LVT_rnkcorr_struc%nsize_total,&
               LVT_rc%strat_nlevels))
          allocate(stats%rnkcorr(m)%obs_value_final(LVT_rc%ngrid, &
               selectNlevs(1),LVT_rnkcorr_struc%nsize_total,&
               LVT_rc%strat_nlevels))
          allocate(stats%rnkcorr(m)%count_value(LVT_rc%ngrid, &
               selectNlevs(1),LVT_rc%strat_nlevels))

          allocate(stats%rnkcorr(m)%rval_r(LVT_rc%ngrid, &
               selectNlevs(1),LVT_rc%strat_nlevels))

          stats%rnkcorr(m)%model_value_final = LVT_rc%udef 
          stats%rnkcorr(m)%obs_value_final   = LVT_rc%udef 
          stats%rnkcorr(m)%count_value       = 0 
          stats%rnkcorr(m)%rval_r            = LVT_rc%udef

          allocate(stats%rnkcorr(m)%value_ci(selectNlevs(1),LVT_rc%strat_nlevels))
          stats%rnkcorr(m)%value_ci = LVT_rc%udef

          if(metric%timeOpt.eq.1) then 

             write(LVT_logunit,*) '[ERR] Temporal calculations are not '
             write(LVT_logunit,*) '[ERR] supported for the rank correlation metric.'
             write(LVT_logunit,*) '[ERR] Program stopping ..'
             call LVT_endrun()

             allocate(stats%rnkcorr(m)%model_value_ts(LVT_rc%ngrid, &
                  selectNlevs(1),&
                  LVT_rnkcorr_struc%nsize_season,LVT_rc%strat_nlevels))
             allocate(stats%rnkcorr(m)%obs_value_ts(LVT_rc%ngrid, &
                  selectNlevs(1),&
                  LVT_rnkcorr_struc%nsize_season,LVT_rc%strat_nlevels))
             allocate(stats%rnkcorr(m)%count_value_ts(LVT_rc%ngrid, &
                  selectNlevs(1),LVT_rc%strat_nlevels))

             allocate(stats%rnkcorr(m)%rval_ts_r(LVT_rc%ngrid, &
                  selectNlevs(1),LVT_rc%strat_nlevels))

             stats%rnkcorr(m)%model_value_ts    = LVT_rc%udef 
             stats%rnkcorr(m)%obs_value_ts      = LVT_rc%udef 
             stats%rnkcorr(m)%count_value_ts    = 0 
             stats%rnkcorr(m)%rval_ts_r            = LVT_rc%udef

             allocate(stats%rnkcorr(m)%tavg_value_ts(LVT_rc%ngrid, &
                  selectNlevs(1),LVT_rc%strat_nlevels))
             allocate(stats%rnkcorr(m)%tavg_count_value_ts(LVT_rc%ngrid, &
                  selectNlevs(1),LVT_rc%strat_nlevels))
             stats%rnkcorr(m)%tavg_value_ts = 0.0
             stats%rnkcorr(m)%tavg_count_value_ts=0 

             if(metric%computeSC.eq.1) then
                allocate(stats%rnkcorr(m)%value_asc(LVT_rc%ngrid, &
                     selectNlevs(1), LVT_rc%nasc))
                stats%rnkcorr(m)%value_asc = 0.0
                allocate(stats%rnkcorr(m)%count_value_asc(LVT_rc%ngrid, &
                     selectNlevs(1), &
                     LVT_rc%nasc))
                stats%rnkcorr(m)%count_value_asc = 0
             endif

             if(metric%computeADC.eq.1) then 
                allocate(stats%rnkcorr(m)%value_adc(LVT_rc%ngrid, &
                     selectNlevs(1), LVT_rc%nadc))
                stats%rnkcorr(m)%value_adc = 0.0
                allocate(stats%rnkcorr(m)%count_value_adc(LVT_rc%ngrid, &
                     selectNlevs(1), &
                     LVT_rc%nadc))
                stats%rnkcorr(m)%count_value_adc = 0
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

  end subroutine LVT_initRnkCorr
  
!BOP
! 
! !ROUTINE: LVT_diagnoseRnkCorr
! \label{LVT_diagnoseRnkCorr}
!
! !INTERFACE:  
  subroutine LVT_diagnoseRnkCorr(pass)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the RnkCorr calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleRnkCorr](\ref{diagnoseSingleRnkCorr})
!     updates the RnkCorr computation for a single variable 
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
       if(LVT_metrics%rnkcorr%selectOpt.eq.1.or.&
            LVT_metrics%rnkcorr%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleRnkCorr(obs,model,stats, &
                  LVT_metrics%rnkcorr)
             
             model => model%next
             obs   => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseRnkCorr

!BOP
! 
! !ROUTINE: diagnoseSingleRnkCorr
! \label{diagnoseSingleRnkCorr}
!
! !INTERFACE: 
  subroutine diagnoseSingleRnkCorr(obs, model, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the RnkCorr computation (updates the running 
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

    integer    :: t,k,m
    integer    :: tindex_f,m_k,o_k

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 

       call getDataTimeIndex(tindex_f, LVT_rc%tavgInterval)

       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(metric%selectOpt.eq.1) then
                   if(model%count(t,m,m_k).ne.0.and.&
                        obs%count(t,m,o_k).ne.0) then 
                      stats%rnkcorr(m)%model_value_final(t,k,tindex_f,1) = & 
                           model%value(t,m,m_k)
                      stats%rnkcorr(m)%obs_value_final(t,k,tindex_f,1) = & 
                           obs%value(t,m,o_k)
                      stats%rnkcorr(m)%count_value(t,k,1) = & 
                           stats%rnkcorr(m)%count_value(t,k,1) + 1
                      
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,m,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%rnkcorr(m)%model_value_final(t,k,tindex_f,2) = & 
                                 model%value(t,m,m_k)
                            stats%rnkcorr(m)%obs_value_final(t,k,tindex_f,2) = & 
                                 obs%value(t,m,o_k)
                            stats%rnkcorr(m)%count_value(t,k,2) = & 
                                 stats%rnkcorr(m)%count_value(t,k,2) + 1
                            
                         elseif(LVT_stats%strat_var(t,m,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            
                            stats%rnkcorr(m)%model_value_final(t,k,tindex_f,3) = & 
                                 model%value(t,m,m_k)
                            stats%rnkcorr(m)%obs_value_final(t,k,tindex_f,3) = & 
                                 obs%value(t,m,o_k)
                            stats%rnkcorr(m)%count_value(t,k,3) = & 
                                 stats%rnkcorr(m)%count_value(t,k,2) + 1
                         endif
                      endif
                   endif
                   
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleRnkCorr


!BOP
! 
! !ROUTINE: LVT_computeRnkCorr
! \label{LVT_computeRnkCorr}
!
! !INTERFACE: 
  subroutine LVT_computeRnkCorr(pass,alarm)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute RnkCorr values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleRnkCorr](\ref{computeSingleRnkCorr})
!     updates the RnkCorr computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     RnkCorr computation has been reached
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

       if(LVT_metrics%rnkcorr%selectOpt.eq.1.or.&
            LVT_metrics%rnkcorr%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%rnkcorr%timeOpt.eq.1.and.&
                  LVT_metrics%rnkcorr%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%rnkcorr%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%rnkcorr%ftn_ts_loc(i,1),200,advance='no') &
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
             call computeSingleRnkCorr(alarm,obs,model,stats,&
                  LVT_metrics%rnkcorr)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo

          if(alarm) then 
             if(LVT_metrics%rnkcorr%timeOpt.eq.1.and.&
                  LVT_metrics%rnkcorr%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%rnkcorr%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%rnkcorr%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif

             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeRnkCorr

!BOP
! 
! !ROUTINE: computeSingleRnkCorr
! \label{computeSingleRnkCorr}
!
! !INTERFACE: 
  subroutine computeSingleRnkCorr(alarm,obs, model,stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the RnkCorr values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     RnkCorr computation has been reached
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

    integer  :: t,l,k,m,i,tind
    real     :: numer
    real     :: denom_sq
    real     :: denom

    real, allocatable :: model_rnks(:)
    real, allocatable :: obs_rnks(:)
    real              :: sxy_r, sxx_r, syy_r, sy_r, sx_r
    integer           :: count_r

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
                      if(stats%rnkcorr(m)%count_value_ts(t,k,l).ne.0.and.&
                           stats%rnkcorr(m)%count_value_ts(t,k,l) &
                           .gt.LVT_rc%obsCountThreshold) then 
                         
                         call computeRnkArray(LVT_rnkcorr_struc%nsize_season,&
                              stats%rnkcorr(m)%model_value_ts(t,k,:,l),&
                              model_rnks(:))
                         
                         call computeRnkArray(LVT_rnkcorr_struc%nsize_season,&
                              stats%rnkcorr(m)%obs_value_ts(t,k,:,l),&
                              obs_rnks(:))                   
                         sxy_r = 0 
                         sxx_r = 0 
                         syy_r = 0 
                         sx_r = 0 
                         sy_r = 0 
                         count_r = 0
                         
                         do i = 1, LVT_rnkcorr_struc%nsize_season
                            if(model_rnks(i).ne.LVT_rc%udef.and.&
                                 obs_rnks(i).ne.LVT_rc%udef) then 
                               sxy_r = sxy_r + model_rnks(i)*obs_rnks(i)
                               sxx_r = sxx_r + model_rnks(i)*model_rnks(i)
                               syy_r = syy_r + obs_rnks(i)*obs_rnks(i)
                               sx_r = sx_r + model_rnks(i)
                               sy_r = sy_r + obs_rnks(i)
                               count_r = count_r + 1
                            endif
                         enddo
                         
                         numer = (float(count_r)*sxy_r - sx_r*sy_r)
                         denom_sq =  (float(count_r)* sxx_r-sx_r**2)* &
                              (float(count_r)*syy_r-sy_r**2)
                         if(denom_sq.gt.0) then 
                            denom = sqrt(denom_sq)
                         else
                            denom = 0.0
                         endif
                         
                         if(denom.ne.0) then 
                            stats%rnkcorr(m)%rval_ts_r(t,k,l) = numer/denom
                         else
                            stats%rnkcorr(m)%rval_ts_r(t,k,l) = LVT_rc%udef
                         endif
                      end if
                   enddo
                   if(metric%computeSC.eq.1) then 
                      if(stats%rnkcorr(m)%rval_ts_r(t,k,l).ne.LVT_rc%udef) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                              tind)
                         stats%rnkcorr(m)%value_asc(t,k,tind) = & 
                              stats%rnkcorr(m)%value_asc(t,k,tind) + & 
                              stats%rnkcorr(m)%rval_ts_r(t,k,l)
                         stats%rnkcorr(m)%count_value_asc(t,k,tind) = & 
                              stats%rnkcorr(m)%count_value_asc(t,k,tind) + 1
                      endif
                   endif
                   if(metric%computeADC.eq.1) then 
                      if(stats%rnkcorr(m)%rval_ts_r(t,k,l).ne.LVT_rc%udef) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%rnkcorr(m)%value_adc(t,k,tind) = & 
                              stats%rnkcorr(m)%value_adc(t,k,tind) + & 
                              stats%rnkcorr(m)%rval_ts_r(t,k,l)
                         stats%rnkcorr(m)%count_value_adc(t,k,tind) = & 
                              stats%rnkcorr(m)%count_value_adc(t,k,tind) + 1
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
                         if(stats%rnkcorr(m)%tavg_count_value_ts(t,k,l).gt.0) then 
                            stats%rnkcorr(m)%tavg_value_ts(t,k,l) = &
                                 stats%rnkcorr(m)%tavg_value_ts(t,k,l)/&
                                 stats%rnkcorr(m)%tavg_count_value_ts(t,k,l)
                         else
                            stats%rnkcorr(m)%tavg_value_ts(t,k,l) = LVT_rc%udef
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
                        stats%rnkcorr(m)%tavg_value_ts,&
                        stats%rnkcorr(m)%tavg_count_value_ts)
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
                                 stats%rnkcorr(m)%tavg_value_ts(t,k,l)
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
                     stats%rnkcorr(1)%tavg_count_value_ts)
                deallocate(tavg_value_ts)                   
             endif
          endif
          
       endif
    endif                  
 
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 

          allocate(model_rnks(LVT_rnkcorr_struc%nsize_total))
          allocate(obs_rnks(LVT_rnkcorr_struc%nsize_total))

          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%rnkcorr(m)%count_value(t,k,l).ne.0.and.&
                           stats%rnkcorr(m)%count_value(t,k,l) &
                           .gt.LVT_rc%obsCountThreshold) then 
                         
                         call computeRnkArray(LVT_rnkcorr_struc%nsize_total,&
                              stats%rnkcorr(m)%model_value_final(t,k,:,l),&
                              model_rnks(:))
                         
                         call computeRnkArray(LVT_rnkcorr_struc%nsize_total,&
                              stats%rnkcorr(m)%obs_value_final(t,k,:,l),&
                              obs_rnks(:))
                         
                         sxy_r = 0 
                         sxx_r = 0 
                         syy_r = 0 
                         sx_r = 0 
                         sy_r = 0 
                         count_r = 0
                         
                         do i = 1, LVT_rnkcorr_struc%nsize_total
                            if(model_rnks(i).ne.LVT_rc%udef.and.&
                                 obs_rnks(i).ne.LVT_rc%udef) then 
                               sxy_r = sxy_r + model_rnks(i)*obs_rnks(i)
                               sxx_r = sxx_r + model_rnks(i)*model_rnks(i)
                               syy_r = syy_r + obs_rnks(i)*obs_rnks(i)
                               sx_r = sx_r + model_rnks(i)
                               sy_r = sy_r + obs_rnks(i)
                               count_r = count_r + 1
                            endif
                         enddo
                         
                         numer = (float(count_r)*sxy_r - sx_r*sy_r)
                         denom_sq =  (float(count_r)* sxx_r-sx_r**2)* &
                              (float(count_r)*syy_r-sy_r**2)
                         if(denom_sq.gt.0) then 
                            denom = sqrt(denom_sq)
                         else
                            denom = 0.0
                         endif
                         
                         if(denom.ne.0) then 
                            stats%rnkcorr(m)%rval_r(t,k,l) = numer/denom
                         else
                            stats%rnkcorr(m)%rval_r(t,k,l) = LVT_rc%udef
                         endif
                      endif
                   enddo
                   if(metric%computeSC.eq.1) then 
                      do l=1,LVT_rc%nasc
                         if(stats%rnkcorr(m)%count_value_asc(t,k,l).gt.&
                              LVT_rc%SCCountThreshold) then 
                            stats%rnkcorr(m)%value_asc(t,k,l) = &
                                 stats%rnkcorr(m)%value_asc(t,k,l)/&
                                 stats%rnkcorr(m)%count_value_asc(t,k,l) 
                         else
                            stats%rnkcorr(m)%value_asc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   if(metric%computeADC.eq.1) then 
                      do l=1,LVT_rc%nadc
                         if(stats%rnkcorr(m)%count_value_adc(t,k,l).gt.&
                              LVT_rc%ADCCountThreshold) then 
                            stats%rnkcorr(m)%value_adc(t,k,l) = &
                                 stats%rnkcorr(m)%value_adc(t,k,l)/&
                                 stats%rnkcorr(m)%count_value_adc(t,k,l) 
                         else
                            stats%rnkcorr(m)%value_adc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                
                enddo
             enddo
          enddo
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%rnkcorr(m)%rval_r(:,k,l),LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%rnkcorr(m)%value_ci(k,l))
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
                              stats%rnkcorr(m)%rval_r(t,k,l)
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
                         if(stats%rnkcorr(m)%value_asc(t,k,l).ne.LVT_rc%udef) then 
                            value_asc(t,k,l) = &
                                 value_asc(t,k,l) + & 
                                 stats%rnkcorr(m)%value_asc(t,k,l)
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
                  stats%rnkcorr(m)%count_value_asc)
             
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
                              stats%rnkcorr(m)%value_adc(t,k,l)
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
                  stats%rnkcorr(m)%count_value_adc)
             
             deallocate(value_adc)
          endif


       endif
    endif

  end subroutine computeSingleRnkCorr

!BOP
! 
! !ROUTINE: LVT_writeMetric_RnkCorr
! \label{LVT_writeMetric_RnkCorr}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_RnkCorr(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%rnkcorr%npass) then
       if(final.ne.1) then
          if(stats%selectOpt.eq.1) then 

             allocate(value_ts(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_ts(LVT_rc%ngrid,LVT_rc%nensem))

             do k=1,vlevels
                do l=1,LVT_rc%strat_nlevels
                   do m=1,LVT_rc%nensem
                      value_ts(:,m) = &
                           stats%rnkcorr(m)%tavg_value_ts(:,k,l)
                      count_value_ts(:,m) = & 
                           stats%rnkcorr(m)%tavg_count_value_ts(:,k,l)
                   enddo
                   if(LVT_metrics%rnkcorr%timeOpt.eq.1) then 

                      call LVT_writevar_gridded(LVT_metrics%rnkcorr%ftn_ts, &
                           value_ts(:,:),&
                           stats%vid_ts(LVT_RNKCORRid,1),k)

                      call LVT_writevar_gridded(LVT_metrics%rnkcorr%ftn_ts, &
                           real(count_value_ts(:,:)),&
                           stats%vid_count_ts(LVT_RNKCORRid,1),k)
                   endif

                enddo
             enddo

             deallocate(value_ts)
             deallocate(count_value_ts)

          endif
       else
          if(pass.eq.LVT_metrics%rnkcorr%npass) then
             if(stats%selectOpt.eq.1) then

                allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(count_value_total(LVT_rc%ngrid,LVT_rc%nensem))

                allocate(value_ci(LVT_rc%nensem))

                if(LVT_metrics%rnkcorr%computeSC.eq.1) then 
                   allocate(value_asc(LVT_rc%ngrid,&
                        LVT_rc%nensem,&
                        LVT_rc%nasc))
                endif
                if(LVT_metrics%rnkcorr%computeADC.eq.1) then 
                   allocate(value_adc(LVT_rc%ngrid,&
                        LVT_rc%nensem,&
                        LVT_rc%nadc))        

                endif

                do k=1,vlevels
                   do l=1,LVT_rc%strat_nlevels
                      do m=1,LVT_rc%nensem
                         value_total(:,m) = &
                              stats%rnkcorr(m)%rval_r(:,k,l)
                         count_value_total(:,m) = & 
                              stats%rnkcorr(m)%count_value(:,k,l)
                         value_ci(m) = stats%rnkcorr(m)%value_ci(k,l)

                         if(LVT_metrics%rnkcorr%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               value_asc(:,m,tind) = & 
                                    stats%rnkcorr(m)%value_asc(:,k,tind)
                            enddo
                         endif

                         if(LVT_metrics%rnkcorr%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               value_adc(:,m,tind) = & 
                                    stats%rnkcorr(m)%value_adc(:,k,tind)
                            enddo
                         endif
                      enddo
                      if(LVT_metrics%rnkcorr%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%rnkcorr%ftn_total, &
                              value_total(:,:),&
                              stats%vid_total(LVT_RNKCORRid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%rnkcorr%ftn_total, &
                              real(count_value_total(:,:)),&
                              stats%vid_count_total(LVT_RNKCORRid,1),k)

                         if(LVT_metrics%rnkcorr%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%rnkcorr%ftn_total,&
                                    value_asc(:,:,tind),&
                                    stats%vid_sc_total(tind,LVT_RNKCORRid,1),k)
                            enddo
                         endif
                         if(LVT_metrics%rnkcorr%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%rnkcorr%ftn_total,&
                                    value_adc(:,:,tind),&
                                    stats%vid_adc_total(tind,LVT_RNKCORRid,1),k)
                            enddo
                         endif
                         call LVT_writeSummaryStats(&
                              LVT_metrics%rnkcorr%ftn_summ,&
                              l,&
                              LVT_metrics%rnkcorr%short_name,&
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

                if(LVT_metrics%rnkcorr%computeSC.eq.1) then 
                   deallocate(value_asc)
                endif
                if(LVT_metrics%rnkcorr%computeADC.eq.1) then 
                   deallocate(value_adc)
                endif

             endif
          endif
       endif
    endif


  end subroutine LVT_writeMetric_RnkCorr


!BOP
! 
! !ROUTINE: LVT_resetMetric_RnkCorr
! \label(LVT_resetMetric_RnkCorr)
!
! !INTERFACE:
  subroutine LVT_resetMetric_RnkCorr(alarm)

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
                if(LVT_metrics%rnkcorr%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      
                      if(alarm) then 
                         stats%rnkcorr(m)%rval_ts_r(:,k,l) = 0.0
                         stats%rnkcorr(m)%count_value_ts(:,k,l)=0 
                         
                         stats%rnkcorr(m)%tavg_value_ts(:,k,l) = 0.0
                         stats%rnkcorr(m)%tavg_count_value_ts(:,k,l)=0 
                      endif
                      
                   enddo
                endif
                
             enddo
          enddo
       endif
       model => model%next
       stats => stats%next

    enddo
  end subroutine LVT_resetMetric_RnkCorr


!BOP
! 
! !ROUTINE: LVT_writerestart_RnkCorr
! 
! !INTERFACE:
  subroutine LVT_writerestart_RnkCorr(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for RnkCorr metric computations
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

    if(LVT_metrics%rnkcorr%selectOpt.eq.1) then 
       
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)

       do while(associated(model))
          
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels     
                      call LVT_writevar_restart(ftn,stats%rnkcorr(m)%rval_r(:,k,l))
                      call LVT_writevar_restart(ftn,stats%rnkcorr(m)%count_value(:,k,l))
                   enddo
                enddo
                
                if(LVT_metrics%rnkcorr%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_writevar_restart(ftn,&
                              stats%rnkcorr(m)%value_asc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%rnkcorr(m)%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%rnkcorr%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_writevar_restart(ftn,&
                              stats%rnkcorr(m)%value_adc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%rnkcorr(m)%count_value_adc(:,k,l))
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
    
  end subroutine LVT_writerestart_RnkCorr

!BOP
! 
! !ROUTINE: LVT_readrestart_RnkCorr
! 
! !INTERFACE:
  subroutine LVT_readrestart_RnkCorr(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for RnkCorr metric computations
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

    if(LVT_metrics%rnkcorr%selectOpt.eq.1) then 
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_readvar_restart(ftn,stats%rnkcorr(m)%rval_r(:,k,l))
                      call LVT_readvar_restart(ftn,stats%rnkcorr(m)%count_value(:,k,l))
                   enddo
                enddo
                
                if(LVT_metrics%rnkcorr%computeSC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc
                         call LVT_readvar_restart(ftn,&
                              stats%rnkcorr(m)%value_asc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%rnkcorr(m)%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%rnkcorr%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_readvar_restart(ftn,&
                              stats%rnkcorr(m)%value_adc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%rnkcorr(m)%count_value_adc(:,k,l))
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

  end subroutine LVT_readrestart_RnkCorr

!BOP
! 
! !ROUTINE: computeAnalysisLength
! \label{computeAnalysisLength}
! 
! !INTERFACE: 
  subroutine computeAnalysisLength(nsize, tavgInterval)
! !ARGUMENTS:     
    integer             :: nsize
    integer, intent(in) :: tavgInterval
!
! !DESCRIPTION: 
!  This subroutine computes the length of the analysis time period
! 
!  The arguments are:
!  \begin{description}
!   \item[nsize] 
!     number of months computed by the routine
!   \item[tavgInterval] 
!     temporal averaging interval (expected to be a multiple of months)
!  \end{description}
!EOP    
    type(ESMF_Time)  :: startTime, stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer          :: status

    call ESMF_TimeSet(startTime, yy = LVT_rc%syr, &
         mm = LVT_rc%smo, &
         dd = LVT_rc%sda, &
         h  = LVT_rc%shr, &
         m  = LVT_rc%smn, & 
         s  = LVT_rc%sss, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status, 'ESMF_TimeSet failed')

    call ESMF_TimeSet(stopTime, yy = LVT_rc%eyr, &
         mm = LVT_rc%emo, &
         dd = LVT_rc%eda, &
         h  = LVT_rc%ehr, &
         m  = LVT_rc%emn, & 
         s  = LVT_rc%ess, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status, 'ESMF_TimeSet failed')

    call ESMF_TimeIntervalSet(timeStep, s = tavgInterval, &
         rc=status)
    call LVT_verify(status,'ESMF_TimeIntervalSet failed')

    nsize = nint((stopTime - startTime)/timestep) + 1

  end subroutine computeAnalysisLength

!BOP
! 
! !ROUTINE: getDataTimeIndex
! \label{getDataTimeIndex}
! 
! !INTERFACE: 
  subroutine getDataTimeIndex(nsize, tavgInterval)
! !ARGUMENTS:     
    integer             :: nsize
    integer, intent(in) :: tavgInterval
!
! !DESCRIPTION: 
!  This subroutine computes the number of months from the start time
!  to the current time, based on the time interval
! 
!  The arguments are:
!  \begin{description}
!   \item[nsize] 
!     number of months computed by the routine
!   \item[tavgInterval] 
!     temporal averaging interval (expected to be a multiple of months)
!  \end{description}
!EOP    

    type(ESMF_Time)  :: currTime, startTime
    type(ESMF_TimeInterval) :: timeStep
    integer          :: status

    call ESMF_TimeSet(currTime, yy = LVT_rc%nyr, &
         mm = LVT_rc%nmo, &
         dd = LVT_rc%nda, &
         h  = LVT_rc%nhr, &
         m  = LVT_rc%nmn, & 
         s  = LVT_rc%nss, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status, 'ESMF_TimeSet failed')

    call ESMF_TimeSet(startTime, yy = LVT_rc%syr, &
         mm = LVT_rc%smo, &
         dd = LVT_rc%sda, &
         h  = LVT_rc%shr, &
         m  = LVT_rc%smn, & 
         s  = LVT_rc%sss, &
         calendar = LVT_calendar, & 
         rc = status)
    call LVT_verify(status, 'ESMF_TimeSet failed')

    call ESMF_TimeIntervalSet(timeStep, s = tavgInterval, &
         rc=status)
    call LVT_verify(status, 'ESMF_TimeIntervalSet failed')

    nsize  = nint((currTime - startTime)/timeStep) + 1

  end subroutine getDataTimeIndex

  subroutine computeRnkArray(nsize,data_value, data_rnks)

    integer  :: nsize
    real     :: data_value(nsize)
    real     :: data_rnks(nsize)

    integer  :: i,j
    integer  :: tmprnk
    
    do i=1,nsize
       tmprnk = 0
       if(data_value(i).ne.LVT_rc%udef) then 
          do j=1,nsize
             if(data_value(j).ne.LVT_rc%udef) then 
                if(data_value(j).lt.data_value(i)) then 
                   tmprnk = tmprnk + 1
                endif
             endif
          enddo
          data_rnks(i) = tmprnk + 1
       else
          data_rnks(i)= LVT_rc%udef
       endif
    enddo
  end subroutine computeRnkArray
end module LVT_RnkCorrMod
