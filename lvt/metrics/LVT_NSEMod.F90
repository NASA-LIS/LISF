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
! !MODULE: LVT_NSEMod
! \label(LVT_NSEMod)
!
! !INTERFACE:
module LVT_NSEMod
! 
! !USES:   
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod

  private
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This module handles the Nash-SutCliffe-Efficiency (NSE)
!   computations by comparing the LIS output to the specified 
!   observations. 
!   
!   NSE is defined as : 
!    
!    NSE = (1-(sum((obs-model)^2)/sum((obs-mean(obs))^2)))
!
!  NOTES: Stratification by seasons is not supported for NSE computation
!  TODO: add codes for stratification, confidence intervals
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 2 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initNSE
  public :: LVT_diagnoseNSE
  public :: LVT_computeNSE
  public :: LVT_writeMetric_NSE
  public :: LVT_resetMetric_NSE
  public :: LVT_writerestart_NSE
  public :: LVT_readrestart_NSE


contains
  
!BOP
! 
! !ROUTINE: LVT_initNSE
! \label{LVT_initNSE}
!
! !INTERFACE: 
  subroutine LVT_initNSE(selectNlevs,stats,metric)
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

    allocate(stats%nse(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       
       if(metric%selectOpt.eq.1) then 
          allocate(stats%nse(m)%value(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))
          stats%nse(m)%value = 0.0
          allocate(stats%nse(m)%value_numer(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))       
          stats%nse(m)%value_numer = 0
          allocate(stats%nse(m)%value_denom(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))       
          stats%nse(m)%value_denom = 0
          allocate(stats%nse(m)%value_obs_mean(LVT_rc%ngrid,selectNlevs(1),&
               LVT_rc%strat_nlevels))
          allocate(stats%nse(m)%count_value_obs_mean(LVT_rc%ngrid,selectNlevs(1),&
               LVT_rc%strat_nlevels))
          stats%nse(m)%value_obs_mean = 0 
          stats%nse(m)%count_value_obs_mean = 0 

          allocate(stats%nse(m)%value_ci(selectNlevs(1),LVT_rc%strat_nlevels))
          stats%nse(m)%value_ci = LVT_rc%udef
       endif

       if(metric%timeopt.eq.1) then 
          allocate(stats%nse(m)%value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%nse(m)%value_ts = 0.0

          allocate(stats%nse(m)%tavg_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%nse(m)%tavg_value_ts = 0.0
          allocate(stats%nse(m)%tavg_count_value_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%nse(m)%tavg_count_value_ts = 0

          allocate(stats%nse(m)%value_numer_ts(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))       
          stats%nse(m)%value_numer_ts = 0
          allocate(stats%nse(m)%value_denom_ts(LVT_rc%ngrid, selectNlevs(1),&
               LVT_rc%strat_nlevels))       
          stats%nse(m)%value_denom_ts = 0
          allocate(stats%nse(m)%value_obs_mean_ts(LVT_rc%ngrid,selectNlevs(1),&
               LVT_rc%strat_nlevels))
          allocate(stats%nse(m)%count_value_obs_mean_ts(LVT_rc%ngrid,selectNlevs(1),&
               LVT_rc%strat_nlevels))
          stats%nse(m)%value_obs_mean_ts = 0 
          stats%nse(m)%count_value_obs_mean_ts = 0 

       end if
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initNSE

!BOP
! 
! !ROUTINE: LVT_diagnoseNSE
! \label{LVT_diagnoseNSE}
!
! !INTERFACE: 
  subroutine LVT_diagnoseNSE(pass)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the NSE calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleNSE](\ref{diagnoseSingleNSE})
!     updates the anomaly correlation computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer :: pass
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_metrics%nse%selectOpt.eq.1.or.&
         LVT_metrics%nse%timeOpt.eq.1) then 
       if(pass.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleNSEobsMean(&
                  obs,model,stats)
             
             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       elseif(pass.eq.2) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)
          
          do while(associated(model))
             call diagnoseSingleNSE(obs,model,stats,&
                  LVT_metrics%nse)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseNSE

!BOP
! 
! !ROUTINE: diagnoseSingleNSEobsMean
! \label{diagnoseSingleNSEobsMean}
!
! !INTERFACE: 
  subroutine diagnoseSingleNSEobsMean(obs, model, stats)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the computation of the 
!   observation mean for a specified variable. 
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
!BOP
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
!EOP
    integer    :: t,k,m,m_k,o_k
    integer    :: c,r

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                
                if(trim(obs%units).eq.trim(model%units)) then                 
                   if(obs%count(t,m,o_k).ne.0) then  
                      stats%nse(m)%value_obs_mean(t,k,1) = &
                           stats%nse(m)%value_obs_mean(t,k,1) + &
                           obs%value(t,m,o_k)
                      stats%nse(m)%count_value_obs_mean(t,k,1) = &
                           stats%nse(m)%count_value_obs_mean(t,k,1) + 1
                      if(LVT_rc%strat_nlevels.gt.1) then 
                         if(LVT_stats%strat_var(t,m,k).gt.&
                              LVT_rc%strat_var_threshold) then 
                            stats%nse(m)%value_obs_mean(t,k,2) = &
                                 stats%nse(m)%value_obs_mean(t,k,2) + &
                                 obs%value(t,m,o_k)
                            stats%nse(m)%count_value_obs_mean(t,k,2) = &
                                 stats%nse(m)%count_value_obs_mean(t,k,2) + 1
                         elseif(LVT_stats%strat_var(t,m,k).le.&
                              LVT_rc%strat_var_threshold) then 
                            stats%nse(m)%value_obs_mean(t,k,3) = &
                                 stats%nse(m)%value_obs_mean(t,k,3) + &
                                 obs%value(t,m,o_k)
                            stats%nse(m)%count_value_obs_mean(t,k,3) = &
                                 stats%nse(m)%count_value_obs_mean(t,k,3) + 1
                         endif
                      endif
                      if(LVT_metrics%nse%timeOpt.eq.1) then 
                         stats%nse(m)%value_obs_mean_ts(t,k,1) = &
                              stats%nse(m)%value_obs_mean_ts(t,k,1) + &
                              obs%value(t,m,o_k)
                         stats%nse(m)%count_value_obs_mean_ts(t,k,1) = &
                              stats%nse(m)%count_value_obs_mean_ts(t,k,1) + 1
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%nse(m)%value_obs_mean_ts(t,k,2) = &
                                    stats%nse(m)%value_obs_mean_ts(t,k,2) + &
                                    obs%value(t,m,o_k)
                               stats%nse(m)%count_value_obs_mean_ts(t,k,2) = &
                                    stats%nse(m)%count_value_obs_mean_ts(t,k,2) + 1
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%nse(m)%value_obs_mean_ts(t,k,3) = &
                                    stats%nse(m)%value_obs_mean_ts(t,k,3) + &
                                    obs%value(t,m,o_k)
                               stats%nse(m)%count_value_obs_mean_ts(t,k,3) = &
                                    stats%nse(m)%count_value_obs_mean_ts(t,k,3) + 1
                            endif
                         endif
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
       enddo
    endif
    
  end subroutine diagnoseSingleNSEobsMean
!BOP
! 
! !ROUTINE: diagnoseSingleNSE
! \label{diagnoseSingleNSE}
!
! !INTERFACE: 
  subroutine diagnoseSingleNSE(obs, model, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the NSE computation for a single variable. 
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
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k,l
    integer    :: c,r

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1, LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1

                if(trim(obs%units).eq.trim(model%units)) then 
                   if(model%count(t,m,m_k).ne.0.and.&
                        obs%count(t,m,o_k).ne.0) then 
                      if(metric%selectOpt.eq.1) then 
                         if(stats%nse(m)%count_value_obs_mean(t,k,1).ne.0) then 
                            stats%nse(m)%value_numer(t,k,1) = stats%nse(m)%value_numer(t,k,1)+&
                                 (obs%value(t,m,o_k)-model%value(t,m,m_k))**2
                            stats%nse(m)%value_denom(t,k,1) = stats%nse(m)%value_denom(t,k,1)+&
                                 (obs%value(t,m,o_k)-stats%nse(m)%value_obs_mean(t,k,1))**2
                         endif
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%nse(m)%value_numer(t,k,2) = stats%nse(m)%value_numer(t,k,2)+&
                                    (obs%value(t,m,o_k)-model%value(t,m,m_k))**2
                               stats%nse(m)%value_denom(t,k,2) = stats%nse(m)%value_denom(t,k,2)+&
                                    (obs%value(t,m,o_k)-stats%nse(m)%value_obs_mean(t,k,2))**2
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%nse(m)%value_numer(t,k,3) = stats%nse(m)%value_numer(t,k,3)+&
                                    (obs%value(t,m,o_k)-model%value(t,m,m_k))**2
                               stats%nse(m)%value_denom(t,k,3) = stats%nse(m)%value_denom(t,k,3)+&
                                    (obs%value(t,m,o_k)-stats%nse(m)%value_obs_mean(t,k,3))**2
                            endif
                         endif
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(stats%nse(m)%count_value_obs_mean_ts(t,k,1).ne.0) then 
                            stats%nse(m)%value_numer_ts(t,k,1) = stats%nse(m)%value_numer_ts(t,k,1)+&
                                 (obs%value(t,m,o_k)-model%value(t,m,m_k))**2
                            stats%nse(m)%value_denom_ts(t,k,1) = stats%nse(m)%value_denom_ts(t,k,1)+&
                                 (obs%value(t,m,o_k)-stats%nse(m)%value_obs_mean(t,k,1))**2
                         endif
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%nse(m)%value_numer_ts(t,k,2) = stats%nse(m)%value_numer_ts(t,k,2)+&
                                    (obs%value(t,m,o_k)-model%value(t,m,m_k))**2
                               stats%nse(m)%value_denom_ts(t,k,2) = stats%nse(m)%value_denom_ts(t,k,2)+&
                                    (obs%value(t,m,o_k)-stats%nse(m)%value_obs_mean(t,k,2))**2
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then 
                               stats%nse(m)%value_numer_ts(t,k,3) = stats%nse(m)%value_numer_ts(t,k,3)+&
                                    (obs%value(t,m,o_k)-model%value(t,m,m_k))**2
                               stats%nse(m)%value_denom_ts(t,k,3) = stats%nse(m)%value_denom_ts(t,k,3)+&
                                    (obs%value(t,m,o_k)-stats%nse(m)%value_obs_mean(t,k,3))**2
                            endif
                         endif
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
       enddo
    endif
  end subroutine diagnoseSingleNSE


!BOP
! 
! !ROUTINE: LVT_computeNSE
! \label{LVT_computeNSE}
!
! !INTERFACE: 
  subroutine LVT_computeNSE(pass,alarm)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine invokes the method to compute NSE values, 
!  for each specified variable. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    
    integer :: pass
    logical :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    integer                          :: i,m

    if(LVT_metrics%nse%selectOpt.eq.1) then
       if(pass.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call computeSingleNSEobsMean(alarm, obs,model,stats)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       elseif(pass.eq.2) then 
          if(LVT_metrics%nse%selectOpt.eq.1.or.&
               LVT_metrics%nse%timeOpt.eq.1) then 
             if(alarm) then 
                if(LVT_metrics%nse%timeOpt.eq.1.and.&
                     LVT_metrics%nse%extractTS.eq.1) then 

                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%nse%ftn_ts_loc(i,m),200,advance='no') &
                                 LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                                 LVT_rc%hr,'',LVT_rc%mn, '' 
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%nse%ftn_ts_loc(i,1),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   endif
                   
                endif
             endif
200          format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)
             
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model))
                
                call computeSingleNSE(alarm,obs,model,stats,&
                     LVT_metrics%nse)
                
                model => model%next
                obs => obs%next
                stats => stats%next
             enddo
             
             if(alarm) then 
                if(LVT_metrics%nse%timeOpt.eq.1.and.&
                     LVT_metrics%nse%extractTS.eq.1) then 

                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%nse%ftn_ts_loc(i,m),fmt='(a1)') ''
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%nse%ftn_ts_loc(i,1),fmt='(a1)') ''
                      enddo
                   endif
                   
                endif
             end if
          endif
       endif
    endif
  end subroutine LVT_computeNSE


!BOP
! 
! !ROUTINE: computeSingleNSEobsMean
! \label{computeSingleNSEobsMean}
!
! !INTERFACE: 
  subroutine computeSingleNSEobsMean(alarm, obs, model, stats)
! 
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_domain
    use LVT_histDataMod 
    use LVT_statsDataMod
    use LVT_logMod,     only : LVT_logunit

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the mean observation value for a specified variable
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
!BOP
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
!EOP
    integer    :: t,k,l,m

    integer    :: c,r

    if(alarm) then 
       if(LVT_metrics%nse%timeOpt.eq.1.and.&
            stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1, LVT_rc%strat_nlevels
                      if(stats%nse(m)%count_value_obs_mean_ts(t,k,l).ne.0) then 
                         stats%nse(m)%value_obs_mean_ts(t,k,l) = &
                              stats%nse(m)%value_obs_mean_ts(t,k,l) /&
                              stats%nse(m)%count_value_obs_mean_ts(t,k,l)
                      endif
                   enddo
                enddo
             enddo
          enddo
       endif
    endif

    if(LVT_rc%endtime.eq.1) then 
       write(LVT_logunit,*) '[INFO] Computing NSE obs Mean '
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1, LVT_rc%strat_nlevels
                      if(stats%nse(m)%count_value_obs_mean(t,k,l).ne.0) then 
                         stats%nse(m)%value_obs_mean(t,k,l) = &
                              stats%nse(m)%value_obs_mean(t,k,l) /&
                              stats%nse(m)%count_value_obs_mean(t,k,l)
                      endif
                   enddo
                enddo
             enddo
          enddo
       endif
    endif
    
  end subroutine computeSingleNSEobsMean

!BOP
! 
! !ROUTINE: computeSingleNSE
! \label{computeSingleNSE}
!
! !INTERFACE: 
  subroutine computeSingleNSE(alarm,obs, model, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the NSE values for a single variable
!  The arguments are: 
!
!  \begin{description}
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
    real     :: numer,denom

    real,    allocatable :: tavg_value_ts(:,:,:)

    real,    allocatable :: value_avg(:,:,:)

    if(metric%timeOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%nse(m)%count_value_obs_mean_ts(t,k,l).ne.0) then 
                         numer = stats%nse(m)%value_numer_ts(t,k,l)
                         denom = stats%nse(m)%value_denom_ts(t,k,l)
                         if(denom.ne.0) then 
                            stats%nse(m)%value_ts(t,k,l) = 1 - numer/denom
                            
                            stats%nse(m)%tavg_value_ts(t,k,l) = & 
                                 stats%nse(m)%tavg_value_ts(t,k,l) + & 
                                 stats%nse(m)%value_ts(t,k,l)
                            stats%nse(m)%tavg_count_value_ts(t,k,l) = & 
                                 stats%nse(m)%tavg_count_value_ts(t,k,l) + 1
                            
                         else
                            stats%nse(m)%value_ts(t,k,l) = LVT_rc%udef
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
                         if(stats%nse(m)%tavg_count_value_ts(t,k,l).gt.0) then 
                            stats%nse(m)%tavg_value_ts(t,k,l) = & 
                                 stats%nse(m)%tavg_value_ts(t,k,l)/&
                                 stats%nse(m)%tavg_count_value_ts(t,k,l)
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
                           stats%nse(m)%tavg_value_ts,&
                           stats%nse(m)%tavg_count_value_ts)
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
                                    stats%nse(m)%tavg_value_ts(t,k,l)
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
                        stats%nse(1)%tavg_count_value_ts)
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
                      if(stats%nse(m)%count_value_obs_mean(t,k,l).ne.0) then 
                         numer = stats%nse(m)%value_numer(t,k,l)
                         denom = stats%nse(m)%value_denom(t,k,l)
                         if(denom.ne.0) then 
                            stats%nse(m)%value(t,k,l) = 1 - numer/denom
                         else
                            stats%nse(m)%value(t,k,l) = LVT_rc%udef
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo

          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%nse(m)%value(:,k,l),LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%nse(m)%value_ci(k,l))
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
                              stats%nse(m)%value(t,k,l)
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

       endif
    endif
  end subroutine computeSingleNSE

!BOP
! 
! !ROUTINE: LVT_writeMetric_NSE
! \label{LVT_writeMetric_NSE}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_NSE(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%nse%npass) then
       if(final.ne.1) then
          if(stats%selectOpt.eq.1) then 

             allocate(value_ts(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_ts(LVT_rc%ngrid,LVT_rc%nensem))

             do k=1,vlevels
                do l=1,LVT_rc%strat_nlevels
                   do m=1,LVT_rc%nensem
                      value_ts(:,m) = &
                           stats%nse(m)%tavg_value_ts(:,k,l)
                      count_value_ts(:,m) = & 
                           stats%nse(m)%tavg_count_value_ts(:,k,l)
                   enddo
                   if(LVT_metrics%nse%timeOpt.eq.1) then 

                      call LVT_writevar_gridded(LVT_metrics%nse%ftn_ts, &
                           value_ts(:,:),&
                           stats%vid_ts(LVT_NSEid,1),k)

                      call LVT_writevar_gridded(LVT_metrics%nse%ftn_ts, &
                           real(count_value_ts(:,:)),&
                           stats%vid_count_ts(LVT_NSEid,1),k)
                   endif

                enddo
             enddo

             deallocate(value_ts)
             deallocate(count_value_ts)

          endif
       else
          if(pass.eq.LVT_metrics%nse%npass) then
             if(stats%selectOpt.eq.1) then

                allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(count_value_total(LVT_rc%ngrid,LVT_rc%nensem))

                allocate(value_ci(LVT_rc%nensem))

                do k=1,vlevels
                   do l=1,LVT_rc%strat_nlevels
                      do m=1,LVT_rc%nensem
                         value_total(:,m) = &
                              stats%nse(m)%value(:,k,l)
                         count_value_total(:,m) = & 
                              stats%nse(m)%count_value_obs_mean(:,k,l)
                         value_ci(m) = stats%nse(m)%value_ci(k,l)

                      enddo
                      if(LVT_metrics%nse%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%nse%ftn_total, &
                              value_total(:,:),&
                              stats%vid_total(LVT_NSEid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%nse%ftn_total, &
                              real(count_value_total(:,:)),&
                              stats%vid_count_total(LVT_NSEid,1),k)

                         call LVT_writeSummaryStats(&
                              LVT_metrics%nse%ftn_summ,&
                              l,&
                              LVT_metrics%nse%short_name,&
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


  end subroutine LVT_writeMetric_NSE


  subroutine LVT_resetMetric_NSE(alarm)

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
                if(LVT_metrics%nse%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%nse(m)%value_ts(:,k,l) = 0
                      stats%nse(m)%value_numer_ts(:,k,l) = 0
                      stats%nse(m)%value_denom_ts(:,k,l) = 0
                      stats%nse(m)%value_obs_mean_ts(:,k,l) = 0 
                      stats%nse(m)%count_value_obs_mean_ts(:,k,l) = 0 
                      if(alarm) then 
                         stats%nse(m)%tavg_value_ts(:,k,l) = 0 
                         stats%nse(m)%tavg_count_value_ts(:,k,l) = 0 
                      endif
                   enddo
                endif
                
             enddo
          enddo
       endif
       model => model%next
       stats => stats%next

    enddo


  end subroutine LVT_resetMetric_NSE


!BOP
! 
! !ROUTINE: LVT_writerestart_NSE
! 
! !INTERFACE:
  subroutine LVT_writerestart_NSE(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for NSE metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%nse%selectOpt.eq.1) then 
       
       print*, 'WARNING: The writerestart method is not implemented for NSE'

    end if
    
  end subroutine LVT_writerestart_NSE

!BOP
! 
! !ROUTINE: LVT_readrestart_NSE
! 
! !INTERFACE:
  subroutine LVT_readrestart_NSE(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for NSE metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%nse%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for NSE'
       stop
    end if
    
  end subroutine LVT_readrestart_NSE

end module LVT_NSEMod
