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
! !MODULE: LVT_ConditionalEntropyMod
! \label(LVT_ConditionalEntropyMod)
!
! !INTERFACE:
module LVT_ConditionalEntropyMod
! 
! !USES:   
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod
  use LVT_NumericalRecipesMod

  implicit none

  private
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This module handles the computation of the Conditional Entropy metric. 
!   Conditional Entropy values are calculated as: 
!
!      CE(X;Y)=Sum_x Sum_y PXY(x,y)log [PXY(x,y)/(PY(y))]
!   
! !NOTES: 
!   * Temporal/time series calculations are not supported for this metric
!   * This metric should only be used for analyzing one variable at a time
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Apr 2020    Sujay Kumar  Initial Specification
! 
!EOP


!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initConditionalEntropy
  public :: LVT_diagnoseConditionalEntropy
  public :: LVT_computeConditionalEntropy 
  public :: LVT_writeMetric_ConditionalEntropy
  public :: LVT_resetMetric_ConditionalEntropy
  public :: LVT_writerestart_ConditionalEntropy
  public :: LVT_readrestart_ConditionalEntropy

  integer, parameter :: CE_nbins = 10 
contains
  
!BOP
! 
! !ROUTINE: LVT_initConditionalEntropy
! \label{LVT_initConditionalEntropy}
!
! !INTERFACE: 
  subroutine LVT_initConditionalEntropy(selectNlevs,stats,metric)
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
!  the Conditional Entropy metric computations. 
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

    allocate(stats%ce(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%ce(m)%xmaxval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ce(m)%xminval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ce(m)%ymaxval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ce(m)%yminval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ce(m)%xdelta(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ce(m)%ydelta(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ce(m)%pxy(LVT_rc%ngrid, selectNlevs(1),&
               CE_nbins,CE_nbins))

          allocate(stats%ce(m)%count_total(LVT_rc%ngrid, selectNlevs(1)))
          allocate(stats%ce(m)%value_total(LVT_rc%ngrid, selectNlevs(1)))
          allocate(stats%ce(m)%value_ci(selectNlevs(1)))
         
          stats%ce(m)%pxy = 0 
          stats%ce(m)%value_total = 0 
          stats%ce(m)%count_total = 0 
          stats%ce(m)%xmaxval = -1E12
          stats%ce(m)%xminval = 1E12
          stats%ce(m)%ymaxval = -1E12
          stats%ce(m)%yminval = 1E12
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initConditionalEntropy
  
!BOP
! 
! !ROUTINE: LVT_diagnoseConditionalEntropy
! \label{LVT_diagnoseConditionalEntropy}
!
! !INTERFACE: 
  subroutine LVT_diagnoseConditionalEntropy(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Conditional Entropy 
!   calculation for desired variables. During the first pass, the 
!   range of the variables (max/min) is captured. During the second 
!   pass, the marginal and joint probabilities are computed. 
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleConditionalEntropy](\ref{diagnoseSingleConditionalEntropy})
!     updates the ConditionalEntropy computation for a single variable 
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
       if(LVT_metrics%ce%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMaxMin(obs,&
                  model, stats, &
                  LVT_metrics%ce)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%ce%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleConditionalEntropy(obs,&
                  model, stats, &
                  LVT_metrics%ce)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseConditionalEntropy

!BOP
! 
! !ROUTINE: diagnoseSingleMaxMin
! \label{diagnoseSingleMaxMin}
!
! !INTERFACE: 
  subroutine diagnoseSingleMaxMin(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the calculations for estimating 
!  max/min values for each grid point. 
!
!  The arguments are: 
!
!  \begin{description}
!   \item[obs] observation object (Datastream 2)
!   \item[model] model variable object (Datastream 1)
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
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k
    integer    :: xbinval, ybinval

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(obs%count(t,m,o_k).ne.0.and. &
                     model%count(t,m,m_k).ne.0) then      
                   if(metric%selectOpt.eq.1) then
                      
                      if(model%value(t,m,m_k).gt.stats%ce(m)%xmaxval(t,k)) then 
                         stats%ce(m)%xmaxval(t,k) = model%value(t,m,m_k)
                      endif
                      if(obs%value(t,m,o_k).gt.stats%ce(m)%ymaxval(t,k)) then 
                         stats%ce(m)%ymaxval(t,k) = obs%value(t,m,o_k)
                      endif
                      if(model%value(t,m,m_k).lt.stats%ce(m)%xminval(t,k)) then 
                         stats%ce(m)%xminval(t,k) = model%value(t,m,m_k)
                      endif
                      if(obs%value(t,m,o_k).lt.stats%ce(m)%yminval(t,k)) then 
                         stats%ce(m)%yminval(t,k) = obs%value(t,m,o_k)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleMaxMin

!BOP
! 
! !ROUTINE: diagnoseSingleConditionalEntropy
! \label{diagnoseSingleConditionalEntropy}
!
! !INTERFACE: 
  subroutine diagnoseSingleConditionalEntropy(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the calculations for computing the 
!  conditional entropy metric. 
!
!  The arguments are: 
!
!  \begin{description}
!   \item[obs] observation object (Datastream 2)
!   \item[model] model variable object (Datastream 1)
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
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k
    integer    :: xbinval, ybinval

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(obs%count(t,m,o_k).ne.0.and. &
                     model%count(t,m,m_k).ne.0) then      
                   if(metric%selectOpt.eq.1.and.&
                        stats%ce(m)%xdelta(t,k).gt.0.and.&
                        stats%ce(m)%ydelta(t,k).gt.0) then
                      
                      xbinval = nint((model%value(t,m,m_k) - &
                           stats%ce(m)%xminval(t,k))/&
                           stats%ce(m)%xdelta(t,k)) + 1
                      
                      if(xbinval.le.0) xbinval = 1
                      if(xbinval.gt.CE_nbins) xbinval = CE_nbins
                      
                      ybinval = nint((obs%value(t,m,o_k) - &
                           stats%ce(m)%yminval(t,k))/&
                           stats%ce(m)%ydelta(t,k)) + 1
                      
                      if(ybinval.le.0) ybinval = 1
                      if(ybinval.gt.CE_nbins) ybinval = CE_nbins
                      
                      stats%ce(m)%pxy(t,k,xbinval,ybinval) = &
                           stats%ce(m)%pxy(t,k,xbinval,ybinval) + 1
                      
                      stats%ce(m)%count_total(t,k) = &
                           stats%ce(m)%count_total(t,k) + 1
                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleConditionalEntropy


!BOP
! 
! !ROUTINE: LVT_computeConditionalEntropy
! \label{LVT_computeConditionalEntropy}
!
! !INTERFACE: 
  subroutine LVT_computeConditionalEntropy(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute 
!   Conditional Entropy values for the 
!   desired variables. During the first pass, the max min
!   values are computed. During the second pass, the actual
!   conditional entropy metric is calculated.
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleConditionalEntropy](\ref{computeSingleConditionalEntropy})
!     updates the ConditionalEntropy computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     iteration instance through the time series. 
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     ConditionalEntropy computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    integer               :: i, m,pass
    logical               :: alarm

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then
       if(LVT_metrics%ce%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call computeSingleMaxMin(alarm,obs,&
                  model, stats, LVT_metrics%ce)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%ce%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call computeSingleConditionalEntropy(alarm,obs,&
                  model, stats, LVT_metrics%ce)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
       endif
    endif

  end subroutine LVT_ComputeConditionalEntropy

!BOP
! 
! !ROUTINE: computeSingleMaxMin
! \label{computeSingleMaxMin}
!
! !INTERFACE: 
  subroutine computeSingleMaxMin(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes max/min values for each grid point.
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     ConditionalEntropy computation has been reached
!    \item[obs] observation object (datastream 2)
!    \item[model] model variable object(datastream 1)
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
    integer  :: kk,pp


    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(stats%ce(m)%xmaxval(t,k).ne.-1E12.and.&
                        stats%ce(m)%xminval(t,k).ne.1E12) then 

                      stats%ce(m)%xdelta(t,k) = &
                           (stats%ce(m)%xmaxval(t,k) - &
                           stats%ce(m)%xminval(t,k))/CE_nbins
                   else
                      stats%ce(m)%xdelta(t,k) = LVT_rc%udef
                   endif

                   if(stats%ce(m)%ymaxval(t,k).ne.-1E12.and.&
                        stats%ce(m)%yminval(t,k).ne.1E12) then 

                      stats%ce(m)%ydelta(t,k) = &
                           (stats%ce(m)%ymaxval(t,k) - &
                           stats%ce(m)%yminval(t,k))/CE_nbins
                   else
                      stats%ce(m)%ydelta(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
       endif
    endif


  end subroutine computeSingleMaxMin

!BOP
! 
! !ROUTINE: computeSingleConditionalEntropy
! \label{computeSingleConditionalEntropy}
!
! !INTERFACE: 
  subroutine computeSingleConditionalEntropy(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the Conditional Entropy values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     ConditionalEntropy computation has been reached
!    \item[obs] observation object (datastream 2)
!    \item[model] model variable object(datastream 1)
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
    integer  :: kk,pp
    real     :: sumval
    real     :: py(CE_nbins)

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(stats%ce(m)%count_total(t,k)&
                        .gt.LVT_rc%obsCountThreshold) then 

                      sumval = 0.0
                      do kk=1,CE_nbins
                         do pp=1,CE_nbins
                            sumval = sumval + & 
                                 stats%ce(m)%pxy(t,k,kk,pp)
                         enddo
                      enddo

                      do kk=1,CE_nbins
                         do pp=1,CE_nbins
                            if(sumval.gt.0) then 
                               stats%ce(m)%pxy(t,k,kk,pp) = & 
                                    stats%ce(m)%pxy(t,k,kk,pp)/&
                                    sumval
                            else
                               stats%ce(m)%pxy(t,k,kk,pp) =LVT_rc%udef
                            endif
                         enddo
                      enddo
                      
                      py = 0 

                      do pp=1,CE_nbins
                         py(pp) = py(pp) + sum(stats%ce(m)%pxy(t,k,:,pp))
                      enddo


                      stats%ce(m)%value_total(t,k) = 0 
                      do kk=1,CE_nbins
                         do pp=1,CE_nbins
                            if(py(pp).ne.LVT_rc%udef.and.&
                                 stats%ce(m)%pxy(t,k,kk,pp).ne.LVT_rc%udef.and.&
                                 py(pp).ne.0.and.&
                                 stats%ce(m)%pxy(t,k,kk,pp).ne.0) then
                               stats%ce(m)%value_total(t,k)= &
                                    stats%ce(m)%value_total(t,k)- &
                                    stats%ce(m)%pxy(t,k,kk,pp)*&
                                    log2(stats%ce(m)%pxy(t,k,kk,pp)/&
                                    py(pp))
                            endif
                         enddo
                      enddo

                   else
                      stats%ce(m)%value_total(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%ce(m)%value_total(:,k),&
                        LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%ce(m)%value_ci(k))
                enddo
             enddo
          enddo

       endif
    endif


  end subroutine computeSingleConditionalEntropy

!BOP
! 
! !ROUTINE: LVT_writeMetric_ConditionalEntropy
! \label{LVT_writeMetric_ConditionalEntropy}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_ConditionalEntropy(pass,final,vlevels,stats,obs)
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
!   This routine writes the computed MI values to 
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
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs

    integer                 :: k,l,m,tind

    real,    allocatable    :: value_total(:,:)
    integer, allocatable    :: count_value_total(:,:)
    real,    allocatable    :: value_ci(:)

    real,    allocatable    :: value_ts(:,:)
    integer, allocatable    :: count_value_ts(:,:)

    real,    allocatable    :: value_adc(:,:,:)
    real,    allocatable    :: value_asc(:,:,:)

    if(pass.eq.LVT_metrics%ce%npass) then
       if(pass.eq.LVT_metrics%ce%npass) then
          if(stats%selectOpt.eq.1) then
             
             allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_total(LVT_rc%ngrid,LVT_rc%nensem))
             
             allocate(value_ci(LVT_rc%nensem))
             
             do k=1,vlevels
                do m=1,LVT_rc%nensem
                   value_total(:,m) = &
                        stats%ce(m)%value_total(:,k)
                   count_value_total(:,m) = & 
                        stats%ce(m)%count_total(:,k)
                   value_ci(m) = stats%ce(m)%value_ci(k)
                   
                enddo
                if(LVT_metrics%ce%selectOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%ce%ftn_total, &
                        value_total(:,:),&
                        stats%vid_total(LVT_CEid,1),k)
                   call LVT_writevar_gridded(LVT_metrics%ce%ftn_total, &
                        real(count_value_total(:,:)),&
                        stats%vid_count_total(LVT_CEid,1),k)
                   
                   
                   call LVT_writeSummaryStats(&
                        LVT_metrics%ce%ftn_summ,&
                        1,&
                        LVT_metrics%ce%short_name,&
                        LVT_rc%ngrid,&
                        value_total(:,:), &
                        count_value_total(:,:),&
                        stats%standard_name,&
                        value_ci(:))
                endif
             enddo
             deallocate(value_total)
             deallocate(count_value_total)

             deallocate(value_ci)

          endif
       endif
    endif

  end subroutine LVT_writeMetric_ConditionalEntropy

!BOP
! 
! !ROUTINE: LVT_resetMetric_ConditionalEntropy
! \label(LVT_resetMetric_ConditionalEntropy)
!
! !INTERFACE:
  subroutine LVT_resetMetric_ConditionalEntropy(alarm)
! !INPUT PARAMETERS: 
    logical           :: alarm
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  This routine resets required variables to support the
!  temporal computation of ConditionalEntropy values. 
!
! !REVISION HISTORY: 
! 
!EOP


  end subroutine LVT_resetMetric_ConditionalEntropy

!BOP
! 
! !ROUTINE: LVT_writerestart_ConditionalEntropy
! 
! !INTERFACE:
  subroutine LVT_writerestart_ConditionalEntropy(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for ConditionalEntropy metric computations
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
       
       if(LVT_metrics%ce%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   
                   call LVT_writevar_restart(ftn,&
                        stats%ce(m)%value_total(:,k))
                   call LVT_writevar_restart(ftn,&
                        stats%ce(m)%count_total(:,k))
                enddo
             enddo
          end if
       end if
       
       model => model%next
       obs   => obs%next
       stats => stats%next

    end do
             
  end subroutine LVT_writerestart_ConditionalEntropy

!BOP
! 
! !ROUTINE: LVT_readrestart_ConditionalEntropy
! 
! !INTERFACE:
  subroutine LVT_readrestart_ConditionalEntropy(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for ConditionalEntropy metric computations
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
       
       if(LVT_metrics%ce%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   
                   call LVT_readvar_restart(ftn,&
                        stats%ce(m)%value_total(:,k))
                   call LVT_readvar_restart(ftn,&
                        stats%ce(m)%count_total(:,k))
                enddo
             enddo
          end if
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo

             
  end subroutine LVT_readrestart_ConditionalEntropy

end module LVT_ConditionalEntropyMod
