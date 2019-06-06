!-----------------------BEGIN NOTICE -- DO NOT EDIT----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT----------------------------
!
!BOP
! 
! !MODULE: LVT_TMBMod
! \label(LVT_TMBMod)
!
! !INTERFACE:
module LVT_TMBMod
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
!   This module handles the computation of the Tian Miss Bias metric. 
!   TMB values are calculated as: 
!
!     TMB = (datastream1_value - datastream2_value)
!   
!   provided that the "observation" data stream has a value greater than a 
!   threshold, and the "forecast" data stream does not (a "miss").
!
!   REFERENCE:
!   Tian Y, and Coauthors, 2009:  Component analysis of errors in
!     satellite-based precipititation estimates.  J Geophys Res, 114,
!     D24101, doi:10.1029/2009JD011949.
!
!   The metric calculations are supported across the entire analysis time 
!   period, stratified to a user defined temporal averaging interval, average 
!   seasonal and diurnal cycle and stratified to external datasets. 
! 
! !FILES USED:
!
!  !REVISION HISTORY: 
!  29 Oct 2018    Eric Kemp   Initial specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initTMB
  public :: LVT_diagnoseTMB
  public :: LVT_computeTMB 
  public :: LVT_writeMetric_TMB
  public :: LVT_resetMetric_TMB
  public :: LVT_writerestart_TMB
  public :: LVT_readrestart_TMB

contains

!BOP
! 
! !ROUTINE: LVT_initTMB
! \label{LVT_initTMB}
!
! !INTERFACE: 
  subroutine LVT_initTMB(selectNlevs,stats,metric)
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
!  the TMB computations. 
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

    allocate(stats%TMB(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       
       if(metric%selectOpt.eq.1) then 
          allocate(stats%TMB(m)%value_total(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))
          stats%TMB(m)%value_total = 0.0

          allocate(stats%TMB(m)%count_value_total(LVT_rc%ngrid, &
               selectNlevs(1), &
               LVT_rc%strat_nlevels))       
          stats%TMB(m)%count_value_total = 0

          allocate(stats%TMB(m)%value_ci(selectNlevs(1),LVT_rc%strat_nlevels))
          stats%TMB(m)%value_ci = LVT_rc%udef
       endif

       if(metric%timeopt.eq.1) then 
          allocate(stats%TMB(m)%value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%TMB(m)%value_ts = 0.0

          allocate(stats%TMB(m)%count_value_ts(LVT_rc%ngrid, &
               selectNlevs(1),&
               LVT_rc%strat_nlevels))       
          stats%TMB(m)%count_value_ts = 0 

          allocate(stats%TMB(m)%tavg_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%TMB(m)%tavg_value_ts = 0.0

          allocate(stats%TMB(m)%tavg_count_value_ts(LVT_rc%ngrid, &
               selectNlevs(1),&
               LVT_rc%strat_nlevels))       
          stats%TMB(m)%tavg_count_value_ts = 0 

          if(metric%computeSC.eq.1) then 
             allocate(stats%TMB(m)%value_asc(LVT_rc%ngrid, selectNlevs(1),&
                  LVT_rc%nasc))
             stats%TMB(m)%value_asc = 0.0

             allocate(stats%TMB(m)%count_value_asc(LVT_rc%ngrid, &
                  selectNlevs(1),&
                  LVT_rc%nasc))
             stats%TMB(m)%count_value_asc = 0
          endif
          if(metric%computeADC.eq.1) then 
             allocate(stats%TMB(m)%value_adc(LVT_rc%ngrid, selectNlevs(1),&
                  LVT_rc%nadc))
             stats%TMB(m)%value_adc = 0.0

             allocate(stats%TMB(m)%count_value_adc(LVT_rc%ngrid, &
                  selectNlevs(1),&
                  LVT_rc%nadc))
             stats%TMB(m)%count_value_adc = 0
          endif
       endif
    enddo

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

 end subroutine LVT_initTMB
 
!BOP
! 
! !ROUTINE: LVT_diagnoseTMB
! \label{LVT_diagnoseTMB}
!
! !INTERFACE: 
  subroutine LVT_diagnoseTMB(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Tian Hib Bias calculation 
!   for desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleTMB](\ref{diagnoseSingleTMB})
!     updates the TMB computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
     implicit none

     integer                 :: pass

     type(LVT_metadataEntry), pointer :: model
     type(LVT_metadataEntry), pointer :: obs
     type(LVT_statsEntry)   , pointer :: stats

     if(pass.eq.1) then 
       if(LVT_metrics%TMB%selectOpt.eq.1.or.&
            LVT_metrics%TMB%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)          
          do while(associated(model))
             call diagnoseSingleTMB(obs, model, stats, &
                  LVT_metrics%TMB)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
 end subroutine LVT_diagnoseTMB

!BOP
! 
! !ROUTINE: diagnoseSingleTMB
! \label{diagnoseSingleTMB}
!
! !INTERFACE: 
  subroutine diagnoseSingleTMB(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the Tian Miss Bias computation by 
!  updating the running sum of the differences in values
!  from datastream1 and datastream2.
!
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
     implicit none

    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k
    logical :: is_miss

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                
                if(trim(obs%units).eq.trim(model%units)) then 
                   if(obs%count(t,m,o_k).ne.0.and. &
                        model%count(t,m,m_k).ne.0) then     

                      ! See if this is a miss
                      is_miss = .false.
                      if(.not.&
                           (model%value(t,m,m_k).gt.metric%minthreshold.and.&
                           model%value(t,m,m_k).lt.metric%maxthreshold).and.&
                           (obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                           obs%value(t,m,o_k).lt.metric%maxthreshold)) then
                         is_miss = .true.
                      end if

                      if(metric%selectOpt.eq.1) then

                         if (is_miss) then
                            stats%TMB(m)%value_total(t,k,1) = &
                                 stats%TMB(m)%value_total(t,k,1) + &
                                 (-obs%value(t,m,o_k) )
                         end if

                         stats%TMB(m)%count_value_total(t,k,1) = &
                              stats%TMB(m)%count_value_total(t,k,1) + 1
                        if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then 
                               if (is_miss) then
                                  stats%TMB(m)%value_total(t,k,2) = &
                                       stats%TMB(m)%value_total(t,k,2) + &
                                       (-obs%value(t,m,o_k) )
                               end if
                               stats%TMB(m)%count_value_total(t,k,2) = &
                                    stats%TMB(m)%count_value_total(t,k,2) + 1
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then 
                               if (is_miss) then
                                  stats%TMB(m)%value_total(t,k,3) = &
                                       stats%TMB(m)%value_total(t,k,3) + &
                                       (-obs%value(t,m,o_k) )
                               end if
                               stats%TMB(m)%count_value_total(t,k,3) = & 
                                    stats%TMB(m)%count_value_total(t,k,3) + 1
                            endif
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if (is_miss) then
                            stats%TMB(m)%value_ts(t,k,1) = &
                                 stats%TMB(m)%value_ts(t,k,1) + &
                                 (-obs%value(t,m,o_k) )
                         end if
                         stats%TMB(m)%count_value_ts(t,k,1) = &
                              stats%TMB(m)%count_value_ts(t,k,1)+1
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then
                               if (is_miss) then
                                  stats%TMB(m)%value_ts(t,k,2) = &
                                       stats%TMB(m)%value_ts(t,k,2) + &
                                       (-obs%value(t,m,o_k) )
                               end if
                               stats%TMB(m)%count_value_ts(t,k,2) = &
                                    stats%TMB(m)%count_value_ts(t,k,2)+1
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then 
                               if (is_miss) then
                                  stats%TMB(m)%value_ts(t,k,3) = &
                                       stats%TMB(m)%value_ts(t,k,3) + &
                                       (-obs%value(t,m,o_k))
                               end if
                               stats%TMB(m)%count_value_ts(t,k,3) = &
                                    stats%TMB(m)%count_value_ts(t,k,3)+1
                            endif
                         endif
                      endif
                   endif
                else
                   write(LVT_logunit,*) 'For variable ',&
                        trim(model%standard_name)
                   write(LVT_logunit,*) 'observations are in ',trim(obs%units)
                   write(LVT_logunit,*) 'and LIS output is in ',&
                        trim(model%units)
                   write(LVT_logunit,*) 'please add the support of ',&
                        trim(model%units), ' in the observation plugin'
                   call LVT_endrun
                endif
             enddo
          enddo
       enddo
    endif
    
 end subroutine diagnoseSingleTMB

    
!BOP
! 
! !ROUTINE: LVT_computeTMB
! \label{LVT_computeTMB}
!
! !INTERFACE: 
  subroutine LVT_computeTMB(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute Tian Miss Bias values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleTMB](\ref{computeSingleTMB})
!     updates the TMB computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     iteration instance through the time series. 
!    \item[alarm]
!     boolean flag indicating if the specified interval for 
!     TMB computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    integer               :: i,m,pass
    logical     :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%TMB%selectOpt.eq.1.or.LVT_metrics%TMB%timeOpt.eq.1) then
 
          if(alarm) then 
             if(LVT_metrics%TMB%timeOpt.eq.1.and.&
                  LVT_metrics%TMB%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%TMB%ftn_ts_loc(i,m),200, &
                              advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%TMB%ftn_ts_loc(i,1),200,advance='no') &
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
             
             call computeSingleTMB(alarm,obs, model, stats, &
                  LVT_metrics%TMB)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%TMB%timeOpt.eq.1.and.&
                  LVT_metrics%TMB%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%TMB%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%TMB%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif

             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeTMB

!BOP
! 
! !ROUTINE: computeSingleTMB
! \label{computeSingleTMB}
!
! !INTERFACE: 
  subroutine computeSingleTMB(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the Tian Miss Bias values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[alarm]
!     boolean flag indicating if the specified interval for 
!     TMB computation has been reached
!    \item[obs] observation object (datastream 2)
!    \item[model] model variable object (datastream 1)
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

    integer  :: t,l,m,k,tind

    real,    allocatable :: tavg_value_ts(:,:,:)
    real,    allocatable :: value_asc(:,:,:)
    integer, allocatable :: count_value_asc(:,:,:)

    real,    allocatable :: value_adc(:,:,:)
    real,    allocatable :: value_avg(:,:,:)

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%TMB(m)%count_value_ts(t,k,l).ne.0) then 
                         stats%TMB(m)%value_ts(t,k,l) = &
                              stats%TMB(m)%value_ts(t,k,l)/&
                              stats%TMB(m)%count_value_ts(t,k,l)
                         stats%TMB(m)%tavg_value_ts(t,k,l) = & 
                              stats%TMB(m)%tavg_value_ts(t,k,l) + &
                              stats%TMB(m)%value_ts(t,k,l)
                         stats%TMB(m)%tavg_count_value_ts(t,k,l) = & 
                              stats%TMB(m)%tavg_count_value_ts(t,k,l) + 1
                      else
                         stats%TMB(m)%value_ts(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                   if(metric%computeSC.eq.1) then 
                      if(stats%TMB(m)%count_value_ts(t,k,1).ne.0) then 
                         call LVT_getSeasonalCycleTimeIndex(LVT_rc%scInterval,&
                              tind)
                         stats%TMB(m)%value_asc(t,k,tind) = &
                              stats%TMB(m)%value_asc(t,k,tind)+ &
                              stats%TMB(m)%value_ts(t,k,1)
                         stats%TMB(m)%count_value_asc(t,k,tind) = &
                              stats%TMB(m)%count_value_asc(t,k,tind)+ 1
                      endif
                   endif
                 if(metric%computeADC.eq.1) then 
                      if(stats%TMB(m)%count_value_ts(t,k,1).ne.0) then 
                         call LVT_getADCTimeIndex(tind)
                         stats%TMB(m)%value_adc(t,k,tind) = &
                              stats%TMB(m)%value_adc(t,k,tind)+ &
                              stats%TMB(m)%value_ts(t,k,1)
                         stats%TMB(m)%count_value_adc(t,k,tind) = &
                              stats%TMB(m)%count_value_adc(t,k,tind)+ 1
                      endif
                   endif
                enddo
             enddo
          enddo

         if(alarm) then 
             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         if(stats%TMB(m)%count_value_ts(t,k,l).ne.0) then 
                            stats%TMB(m)%tavg_value_ts(t,k,l) = & 
                                 stats%TMB(m)%tavg_value_ts(t,k,l) / & 
                                 stats%TMB(m)%tavg_count_value_ts(t,k,l) 
                         else
                            stats%TMB(m)%tavg_value_ts(t,k,l) = LVT_rc%udef
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
                           stats%TMB(m)%tavg_value_ts,&
                           stats%TMB(m)%tavg_count_value_ts)
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
                                    stats%TMB(m)%tavg_value_ts(t,k,l)
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
                        stats%TMB(1)%tavg_count_value_ts)
                   deallocate(tavg_value_ts)                   
                endif
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
                      if(stats%TMB(m)%count_value_total(t,k,l).ne.0.and.&
                           stats%TMB(m)%count_value_total(t,k,l)&
                           .gt.LVT_rc%obsCountThreshold) then 
                         
                         stats%TMB(m)%value_total(t,k,l) = &
                              stats%TMB(m)%value_total(t,k,l)/&
                              stats%TMB(m)%count_value_total(t,k,l) 
                      else
                         stats%TMB(m)%value_total(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                   if(metric%computeSC.eq.1) then 
                      do l=1,LVT_rc%nasc
                         if(stats%TMB(m)%count_value_asc(t,k,l).gt.&
                              LVT_rc%SCCountThreshold) then
                            stats%TMB(m)%value_asc(t,k,l) = &
                                 stats%TMB(m)%value_asc(t,k,l)/&
                                 stats%TMB(m)%count_value_asc(t,k,l)
                         else
                            stats%TMB(m)%value_asc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                   if(metric%computeADC.eq.1) then 
                      do l=1,LVT_rc%nadc
                         if(stats%TMB(m)%count_value_adc(t,k,l).gt.&
                              LVT_rc%ADCCountThreshold) then
                            stats%TMB(m)%value_adc(t,k,l) = &
                                 stats%TMB(m)%value_adc(t,k,l)/&
                                 stats%TMB(m)%count_value_adc(t,k,l)
                         else
                            stats%TMB(m)%value_adc(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   endif
                enddo
             enddo
          enddo

          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%TMB(m)%value_total(:,k,l), &
                        LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%TMB(m)%value_ci(k,l))
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
                              stats%TMB(m)%value_total(t,k,l)
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
                         if(stats%TMB(m)%value_asc(t,k,l).ne.LVT_rc%udef) then 
                            value_asc(t,k,l) = &
                                 value_asc(t,k,l) + & 
                                 stats%TMB(m)%value_asc(t,k,l)
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
                  stats%TMB(1)%count_value_asc)
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
                              stats%TMB(m)%value_adc(t,k,l)
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
                  stats%TMB(m)%count_value_adc)
             
             deallocate(value_adc)
          endif

       endif
    endif

 end subroutine computeSingleTMB

!BOP
! 
! !ROUTINE: LVT_writeMetric_TMB
! \label{LVT_writeMetric_TMB}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_TMB(pass,final,vlevels,stats,obs)
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
!   This routine writes the computed Tian Miss Bias values to 
!   an external file. 
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

    if(pass.eq.LVT_metrics%TMB%npass) then
       if(final.ne.1) then
          if(stats%selectOpt.eq.1) then 

             allocate(value_ts(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_ts(LVT_rc%ngrid,LVT_rc%nensem))

             do k=1,vlevels
                do l=1,LVT_rc%strat_nlevels
                   do m=1,LVT_rc%nensem
                      value_ts(:,m) = &
                           stats%TMB(m)%tavg_value_ts(:,k,l)
                      count_value_ts(:,m) = & 
                           stats%TMB(m)%tavg_count_value_ts(:,k,l)
                   enddo
                   if(LVT_metrics%TMB%timeOpt.eq.1) then 

                      call LVT_writevar_gridded(LVT_metrics%TMB%ftn_ts, &
                           value_ts(:,:),&
                           stats%vid_ts(LVT_TMBid,1),k)
                      call LVT_writevar_gridded(LVT_metrics%TMB%ftn_ts, &
                           real(count_value_ts(:,:)),&
                           stats%vid_count_ts(LVT_TMBid,1),k)
                   endif

                enddo
             enddo

             deallocate(value_ts)
             deallocate(count_value_ts)

          endif
       else
          if(pass.eq.LVT_metrics%TMB%npass) then
             if(stats%selectOpt.eq.1) then

                allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(count_value_total(LVT_rc%ngrid,LVT_rc%nensem))

                allocate(value_ci(LVT_rc%nensem))

                if(LVT_metrics%TMB%computeSC.eq.1) then 
                   allocate(value_asc(LVT_rc%ngrid,&
                        LVT_rc%nensem,&
                        LVT_rc%nasc))
                endif
                if(LVT_metrics%TMB%computeADC.eq.1) then 
                   allocate(value_adc(LVT_rc%ngrid,&
                        LVT_rc%nensem,&
                        LVT_rc%nadc))        

                endif

                do k=1,vlevels
                   do l=1,LVT_rc%strat_nlevels
                      do m=1,LVT_rc%nensem
                         value_total(:,m) = &
                              stats%TMB(m)%value_total(:,k,l)
                         count_value_total(:,m) = & 
                              stats%TMB(m)%count_value_total(:,k,l)
                         value_ci(m) = stats%TMB(m)%value_ci(k,l)

                         if(LVT_metrics%TMB%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               value_asc(:,m,tind) = & 
                                    stats%TMB(m)%value_asc(:,k,tind)
                            enddo
                         endif

                         if(LVT_metrics%TMB%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               value_adc(:,m,tind) = & 
                                    stats%TMB(m)%value_adc(:,k,tind)
                            enddo
                         endif
                      enddo
                     if(LVT_metrics%TMB%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%TMB%ftn_total, &
                              value_total(:,:),&
                              stats%vid_total(LVT_TMBid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%TMB%ftn_total, &
                              real(count_value_total(:,:)),&
                              stats%vid_count_total(LVT_TMBid,1),k)

                         if(LVT_metrics%TMB%computeSC.eq.1) then 
                            do tind = 1,LVT_rc%nasc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%TMB%ftn_total,&
                                    value_asc(:,:,tind),&
                                    stats%vid_sc_total(tind,LVT_TMBid,1),k)
                            enddo
                         endif
                         if(LVT_metrics%TMB%computeADC.eq.1) then 
                            do tind = 1,LVT_rc%nadc
                               call LVT_writevar_gridded(&
                                    LVT_metrics%TMB%ftn_total,&
                                    value_adc(:,:,tind),&
                                    stats%vid_adc_total(tind,LVT_TMBid,1),k)
                            enddo
                         endif
                         call LVT_writeSummaryStats(&
                              LVT_metrics%TMB%ftn_summ,&
                              l,&
                              LVT_metrics%TMB%short_name,&
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

                if(LVT_metrics%TMB%computeSC.eq.1) then 
                   deallocate(value_asc)
                endif
                if(LVT_metrics%TMB%computeADC.eq.1) then 
                   deallocate(value_adc)
                endif

             endif
          endif
       endif
    endif

 end subroutine LVT_writeMetric_TMB

!BOP
! 
! !ROUTINE: LVT_resetMetric_TMB
! \label{LVT_resetMetric_TMB}
!
! !INTERFACE: 
  subroutine LVT_resetMetric_TMB(alarm)
! 
! !INPUT PARAMETERS: 
    logical          :: alarm
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine resets required variables to support the
!  temporal computation of Tian Miss Bias values. 
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
                if(LVT_metrics%TMB%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%TMB(m)%value_ts(:,k,l) = 0.0
                      stats%TMB(m)%count_value_ts(:,k,l)=0 
                      if(alarm) then 
                         stats%TMB(m)%tavg_value_ts(:,k,l) = 0.0
                         stats%TMB(m)%tavg_count_value_ts(:,k,l)=0 
                      endif
                   enddo
                endif
             enddo
          enddo
       endif
       model => model%next
       stats => stats%next
    enddo

 end subroutine LVT_resetMetric_TMB

!BOP
! 
! !ROUTINE: LVT_writerestart_TMB
! 
! !INTERFACE:
  subroutine LVT_writerestart_TMB(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Tian Miss Bias metric computations
! 
!EOP

    integer              :: k,l,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))

       if(LVT_metrics%TMB%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_writevar_restart(ftn,&
                           stats%TMB(m)%value_total(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%TMB(m)%count_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%TMB%computeSC.eq.1) then
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc 
                         call LVT_writevar_restart(ftn,&
                              stats%TMB(m)%value_asc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%TMB(m)%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%TMB%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_writevar_restart(ftn,&
                              stats%TMB(m)%value_adc(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%TMB(m)%count_value_adc(:,k,l))
                      enddo
                   enddo
                endif
             enddo
          endif
       endif
       
       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo
  end subroutine LVT_writerestart_TMB

!BOP
! 
! !ROUTINE: LVT_readrestart_TMB
! 
! !INTERFACE:
  subroutine LVT_readrestart_TMB(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Tian Miss Bias metric computations
! 
!EOP

    integer              :: k,l,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%TMB%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_readvar_restart(ftn,&
                           stats%TMB(m)%value_total(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%TMB(m)%count_value_total(:,k,l))
                   enddo
                enddo
                if(LVT_metrics%TMB%computeSC.eq.1) then
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nasc 
                         call LVT_readvar_restart(ftn,&
                              stats%TMB(m)%value_asc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%TMB(m)%count_value_asc(:,k,l))
                      enddo
                   enddo
                endif
                if(LVT_metrics%TMB%computeADC.eq.1) then 
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%nadc
                         call LVT_readvar_restart(ftn,&
                              stats%TMB(m)%value_adc(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%TMB(m)%count_value_adc(:,k,l))
                      enddo
                   enddo
                endif
             enddo
          endif
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo
 end subroutine LVT_readrestart_TMB
                            
end module LVT_TMBMod
