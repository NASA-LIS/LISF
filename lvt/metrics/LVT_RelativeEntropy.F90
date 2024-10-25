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
! !MODULE: LVT_RelativeEntropyMod
! \label(LVT_RelativeEntropyMod)
!
! !INTERFACE:
module LVT_RelativeEntropyMod
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
!   This module handles the computation of the Relative Entropy or 
!   Kullback-Leibler divergence metric. This metric is a measure of 
!   how one probability distribution is different from another probability
!   distribution
! 
!   Relative Entropy values are calculated as: 
!
!      H(X)= Sum_x  P(x)log [P(x)/Q(x)]
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
  public :: LVT_initRelativeEntropy
  public :: LVT_diagnoseRelativeEntropy
  public :: LVT_computeRelativeEntropy 
  public :: LVT_writeMetric_RelativeEntropy
  public :: LVT_resetMetric_RelativeEntropy
  public :: LVT_writerestart_RelativeEntropy
  public :: LVT_readrestart_RelativeEntropy

  integer, parameter :: RE_nbins = 10 
contains
  
!BOP
! 
! !ROUTINE: LVT_initRelativeEntropy
! \label{LVT_initRelativeEntropy}
!
! !INTERFACE: 
  subroutine LVT_initRelativeEntropy(selectNlevs,stats,metric)
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
!  the Information Entropy metric computations. 
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

    allocate(stats%re(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%re(m)%xmaxval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%re(m)%xminval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%re(m)%xdelta(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%re(m)%px(LVT_rc%ngrid, selectNlevs(1),&
               RE_nbins))
          allocate(stats%re(m)%count_total(LVT_rc%ngrid, selectNlevs(1)))
          allocate(stats%re(m)%value_total(LVT_rc%ngrid, selectNlevs(1)))
          allocate(stats%re(m)%value_ci(selectNlevs(1)))

          allocate(stats%re(m)%ymaxval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%re(m)%yminval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%re(m)%ydelta(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%re(m)%py(LVT_rc%ngrid, selectNlevs(1),&
               RE_nbins))

          stats%re(m)%px = 0 
          stats%re(m)%value_total = 0 
          stats%re(m)%count_total = 0 
          stats%re(m)%xmaxval = -1E12
          stats%re(m)%xminval = 1E12

          stats%re(m)%py = 0 
          stats%re(m)%ymaxval = -1E12
          stats%re(m)%yminval = 1E12
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    metric%obsData = .true.
    metric%stdevFlag = .false. 

  end subroutine LVT_initRelativeEntropy
  
!BOP
! 
! !ROUTINE: LVT_diagnoseRelativeEntropy
! \label{LVT_diagnoseRelativeEntropy}
!
! !INTERFACE: 
  subroutine LVT_diagnoseRelativeEntropy(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Information Entropy 
!   calculation for desired variables. During the first pass, the 
!   range of the variables (max/min) is captured. During the second 
!   pass, the probabilities are computed. 
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleRelativeEntropy](\ref{diagnoseSingleRelativeEntropy})
!     updates the Information Entropy computation for a single variable 
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
       if(LVT_metrics%re%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMaxMin(obs,&
                  model, stats, &
                  LVT_metrics%re)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%re%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleRelativeEntropy(obs,&
                  model, stats, &
                  LVT_metrics%re)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseRelativeEntropy

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
                if(model%count(t,m,m_k).ne.0.and.&
                     obs%count(t,m,o_k).ne.0) then      
                   if(metric%selectOpt.eq.1) then
                      
                      if(model%value(t,m,m_k).gt.stats%re(m)%xmaxval(t,k)) then 
                         stats%re(m)%xmaxval(t,k) = model%value(t,m,m_k)
                      endif
                      if(model%value(t,m,m_k).lt.stats%re(m)%xminval(t,k)) then 
                         stats%re(m)%xminval(t,k) = model%value(t,m,m_k)
                      endif

                      if(obs%value(t,m,o_k).gt.stats%re(m)%ymaxval(t,k)) then 
                         stats%re(m)%ymaxval(t,k) = obs%value(t,m,o_k)
                      endif
                      if(obs%value(t,m,o_k).lt.stats%re(m)%yminval(t,k)) then 
                         stats%re(m)%yminval(t,k) = obs%value(t,m,o_k)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    end if
  end subroutine diagnoseSingleMaxMin

!BOP
! 
! !ROUTINE: diagnoseSingleRelativeEntropy
! \label{diagnoseSingleRelativeEntropy}
!
! !INTERFACE: 
  subroutine diagnoseSingleRelativeEntropy(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the calculations for computing the 
!  Information Entropy metric. 
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

                if(model%count(t,m,m_k).ne.0.and.&
                     obs%count(t,m,o_k).ne.0) then      
                   if(metric%selectOpt.eq.1.and.&
                        stats%re(m)%xdelta(t,k).gt.0.and.&
                        stats%re(m)%ydelta(t,k).gt.0) then 
                      xbinval = nint((model%value(t,m,m_k) - &
                           stats%re(m)%xminval(t,k))/&
                           stats%re(m)%xdelta(t,k)) + 1
                      
                      if(xbinval.le.0) xbinval = 1
                      if(xbinval.gt.RE_nbins) xbinval = RE_nbins
                      
                      stats%re(m)%px(t,k,xbinval) = &
                           stats%re(m)%px(t,k,xbinval) + 1

                      ybinval = nint((obs%value(t,m,o_k) - &
                           stats%re(m)%yminval(t,k))/&
                           stats%re(m)%ydelta(t,k)) + 1
                      
                      if(ybinval.le.0) ybinval = 1
                      if(ybinval.gt.RE_nbins) ybinval = RE_nbins
                      
                      stats%re(m)%py(t,k,ybinval) = &
                           stats%re(m)%py(t,k,ybinval) + 1

                      stats%re(m)%count_total(t,k) = & 
                           stats%re(m)%count_total(t,k) + 1

                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleRelativeEntropy


!BOP
! 
! !ROUTINE: LVT_computeRelativeEntropy
! \label{LVT_computeRelativeEntropy}
!
! !INTERFACE: 
  subroutine LVT_computeRelativeEntropy(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute 
!   Information Entropy values for the 
!   desired variables. During the first pass, the max min
!   values are computed. During the second pass, the actual
!   mutual information metric is calculated.
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleRelativeEntropy](\ref{computeSingleRelativeEntropy})
!     updates the Informatio xsnEntropy computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     iteration instance through the time series. 
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     Information Entropy computation has been reached
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
       if(LVT_metrics%re%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call computeSingleMaxMin(alarm,obs,&
                  model, stats, LVT_metrics%re)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%re%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call computeSingleRelativeEntropy(alarm,obs,&
                  model, stats, LVT_metrics%re)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
       endif
    endif

  end subroutine LVT_ComputeRelativeEntropy

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
!     RelativeEntropy computation has been reached
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
                   if(stats%re(m)%xmaxval(t,k).ne.-1E12.and.&
                        stats%re(m)%xminval(t,k).ne.1E12) then 

                      stats%re(m)%xdelta(t,k) = &
                           (stats%re(m)%xmaxval(t,k) - &
                           stats%re(m)%xminval(t,k))/RE_nbins
                   else
                      stats%re(m)%xdelta(t,k) = LVT_rc%udef
                   endif

                   if(stats%re(m)%ymaxval(t,k).ne.-1E12.and.&
                        stats%re(m)%yminval(t,k).ne.1E12) then 
                      
                      stats%re(m)%ydelta(t,k) = &
                           (stats%re(m)%ymaxval(t,k) - &
                           stats%re(m)%yminval(t,k))/RE_nbins
                   else
                      stats%re(m)%ydelta(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
       endif
    endif


  end subroutine computeSingleMaxMin

!BOP
! 
! !ROUTINE: computeSingleRelativeEntropy
! \label{computeSingleRelativeEntropy}
!
! !INTERFACE: 
  subroutine computeSingleRelativeEntropy(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the Information Entropy values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     RelativeEntropy computation has been reached
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
    real     :: px(RE_nbins)
    real     :: py(RE_nbins)

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(stats%re(m)%count_total(t,k)&
                        .gt.LVT_rc%obsCountThreshold) then 

                      sumval = 0.0
                      do kk=1,RE_nbins
                         sumval = sumval + & 
                              stats%re(m)%px(t,k,kk)
                      enddo

                      do kk=1,RE_nbins
                         if(sumval.gt.0) then 
                            stats%re(m)%px(t,k,kk) = & 
                                 stats%re(m)%px(t,k,kk)/&
                                 sumval
                         else
                            stats%re(m)%px(t,k,kk) =LVT_rc%udef
                         endif
                      enddo

                      
                      sumval = 0.0
                      do kk=1,RE_nbins
                         sumval = sumval + & 
                              stats%re(m)%py(t,k,kk)
                      enddo

                      do kk=1,RE_nbins
                         if(sumval.gt.0) then 
                            stats%re(m)%py(t,k,kk) = & 
                                 stats%re(m)%py(t,k,kk)/&
                                 sumval
                         else
                            stats%re(m)%py(t,k,kk) =LVT_rc%udef
                         endif
                      enddo                      
                      stats%re(m)%value_total(t,k) = 0 
                      do kk=1,RE_nbins
                         if(stats%re(m)%px(t,k,kk).ne.LVT_rc%udef.and.&
                              stats%re(m)%px(t,k,kk).ne.0.and.&
                              stats%re(m)%py(t,k,kk).ne.LVT_rc%udef.and.&
                              stats%re(m)%py(t,k,kk).ne.0) then 
                            stats%re(m)%value_total(t,k) = & 
                                 stats%re(m)%value_total(t,k) + &
                                 stats%re(m)%px(t,k,kk)*log2(&
                                 (stats%re(m)%px(t,k,kk))/&
                                 stats%re(m)%py(t,k,kk))                        
                         endif
                      enddo
                   else
                      stats%re(m)%value_total(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo

          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                call LVT_computeCI(stats%re(m)%value_total(:,k),&
                     LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%re(m)%value_ci(k))
             enddo
          enddo
       endif
    endif


  end subroutine computeSingleRelativeEntropy

!BOP
! 
! !ROUTINE: LVT_writeMetric_RelativeEntropy
! \label{LVT_writeMetric_RelativeEntropy}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_RelativeEntropy(pass,final,vlevels,stats,obs)
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
!   This routine writes the computed information entropy values to 
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

    real,    allocatable    :: value_adc(:,:,:)
    real,    allocatable    :: value_asc(:,:,:)

    if(pass.eq.LVT_metrics%re%npass) then
       if(pass.eq.LVT_metrics%re%npass) then
          if(stats%selectOpt.eq.1) then
             
             allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_total(LVT_rc%ngrid,LVT_rc%nensem))
             
             allocate(value_ci(LVT_rc%nensem))
             
             do k=1,vlevels
                do m=1,LVT_rc%nensem
                   value_total(:,m) = &
                        stats%re(m)%value_total(:,k)
                   count_value_total(:,m) = & 
                        stats%re(m)%count_total(:,k)
                   value_ci(m) = stats%re(m)%value_ci(k)
                   
                enddo

                if(LVT_metrics%re%selectOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%re%ftn_total, &
                        value_total(:,:),&
                        stats%vid_total(LVT_REid,1),k)
                   call LVT_writevar_gridded(LVT_metrics%re%ftn_total, &
                        real(count_value_total(:,:)),&
                        stats%vid_count_total(LVT_REid,1),k)
                   
                   call LVT_writeSummaryStats(&
                        LVT_metrics%re%ftn_summ,&
                        1,&
                        LVT_metrics%re%short_name,&
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

  end subroutine LVT_writeMetric_RelativeEntropy

!BOP
! 
! !ROUTINE: LVT_resetMetric_RelativeEntropy
! \label(LVT_resetMetric_RelativeEntropy)
!
! !INTERFACE:
  subroutine LVT_resetMetric_RelativeEntropy(alarm)
! !INPUT PARAMETERS: 
    logical           :: alarm
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  This routine resets required variables to support the
!  temporal computation of Information Entropy values. 
!
! !REVISION HISTORY: 
! 
!EOP


  end subroutine LVT_resetMetric_RelativeEntropy

!BOP
! 
! !ROUTINE: LVT_writerestart_RelativeEntropy
! 
! !INTERFACE:
  subroutine LVT_writerestart_RelativeEntropy(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Information Entropy metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

             
  end subroutine LVT_writerestart_RelativeEntropy

!BOP
! 
! !ROUTINE: LVT_readrestart_RelativeEntropy
! 
! !INTERFACE:
  subroutine LVT_readrestart_RelativeEntropy(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Information Entropy metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

  end subroutine LVT_readrestart_RelativeEntropy


end module LVT_RelativeEntropyMod
