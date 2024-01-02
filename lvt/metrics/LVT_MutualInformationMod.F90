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
! !MODULE: LVT_MutualInformationMod
! \label(LVT_MutualInformationMod)
!
! !INTERFACE:
module LVT_MutualInformationMod
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
!   This module handles the computation of the MutualInformation metric. 
!   MutualInformation values are calculated as: 
!
!      I(X;Y)=Sum_x Sum_y PXY(x,y)log [PXY(x,y)/(PX(x)PY(y))]
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
  public :: LVT_initMutualInformation
  public :: LVT_diagnoseMutualInformation
  public :: LVT_computeMutualInformation 
  public :: LVT_writeMetric_MutualInformation
  public :: LVT_resetMetric_MutualInformation
  public :: LVT_writerestart_MutualInformation
  public :: LVT_readrestart_MutualInformation

  integer, parameter :: MI_nbins = 10 
contains
  
!BOP
! 
! !ROUTINE: LVT_initMutualInformation
! \label{LVT_initMutualInformation}
!
! !INTERFACE: 
  subroutine LVT_initMutualInformation(selectNlevs,stats,metric)
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
!  the Mutual Information metric computations. 
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

    allocate(stats%mi(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%mi(m)%xmaxval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%mi(m)%xminval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%mi(m)%ymaxval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%mi(m)%yminval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%mi(m)%xdelta(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%mi(m)%ydelta(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%mi(m)%pxy(LVT_rc%ngrid, selectNlevs(1),&
               MI_nbins,MI_nbins))

          allocate(stats%mi(m)%count_total(LVT_rc%ngrid, selectNlevs(1)))
          allocate(stats%mi(m)%value_total(LVT_rc%ngrid, selectNlevs(1)))
          allocate(stats%mi(m)%value_ci(selectNlevs(1)))
         
          stats%mi(m)%pxy = 0 
          stats%mi(m)%value_total = 0 
          stats%mi(m)%count_total = 0 
          stats%mi(m)%xmaxval = -1E12
          stats%mi(m)%xminval = 1E12
          stats%mi(m)%ymaxval = -1E12
          stats%mi(m)%yminval = 1E12
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    metric%obsData = .false. 
    metric%stdevFlag = .true. 

  end subroutine LVT_initMutualInformation
  
!BOP
! 
! !ROUTINE: LVT_diagnoseMutualInformation
! \label{LVT_diagnoseMutualInformation}
!
! !INTERFACE: 
  subroutine LVT_diagnoseMutualInformation(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Mutual Information 
!   calculation for desired variables. During the first pass, the 
!   range of the variables (max/min) is captured. During the second 
!   pass, the marginal and joint probabilities are computed. 
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleMutualInformation](\ref{diagnoseSingleMutualInformation})
!     updates the MutualInformation computation for a single variable 
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
       if(LVT_metrics%mi%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMaxMin(obs,&
                  model, stats, &
                  LVT_metrics%mi)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%mi%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMutualInformation(obs,&
                  model, stats, &
                  LVT_metrics%mi)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseMutualInformation

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
                      
                      if(model%value(t,m,m_k).gt.stats%mi(m)%xmaxval(t,k)) then 
                         stats%mi(m)%xmaxval(t,k) = model%value(t,m,m_k)
                      endif
                      if(obs%value(t,m,o_k).gt.stats%mi(m)%ymaxval(t,k)) then 
                         stats%mi(m)%ymaxval(t,k) = obs%value(t,m,o_k)
                      endif
                      if(model%value(t,m,m_k).lt.stats%mi(m)%xminval(t,k)) then 
                         stats%mi(m)%xminval(t,k) = model%value(t,m,m_k)
                      endif
                      if(obs%value(t,m,o_k).lt.stats%mi(m)%yminval(t,k)) then 
                         stats%mi(m)%yminval(t,k) = obs%value(t,m,o_k)
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
! !ROUTINE: diagnoseSingleMutualInformation
! \label{diagnoseSingleMutualInformation}
!
! !INTERFACE: 
  subroutine diagnoseSingleMutualInformation(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the calculations for computing the 
!  mutual information metric. 
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
                        stats%mi(m)%xdelta(t,k).gt.0.and.&
                        stats%mi(m)%ydelta(t,k).gt.0) then
                      
                      xbinval = nint((model%value(t,m,m_k) - &
                           stats%mi(m)%xminval(t,k))/&
                           stats%mi(m)%xdelta(t,k)) + 1
                      
                      if(xbinval.le.0) xbinval = 1
                      if(xbinval.gt.MI_nbins) xbinval = MI_nbins
                      
                      ybinval = nint((obs%value(t,m,o_k) - &
                           stats%mi(m)%yminval(t,k))/&
                           stats%mi(m)%ydelta(t,k)) + 1
                      
                      if(ybinval.le.0) ybinval = 1
                      if(ybinval.gt.MI_nbins) ybinval = MI_nbins
                      
                      stats%mi(m)%pxy(t,k,xbinval,ybinval) = &
                           stats%mi(m)%pxy(t,k,xbinval,ybinval) + 1
                      
                      stats%mi(m)%count_total(t,k) = &
                           stats%mi(m)%count_total(t,k) + 1
                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleMutualInformation


!BOP
! 
! !ROUTINE: LVT_computeMutualInformation
! \label{LVT_computeMutualInformation}
!
! !INTERFACE: 
  subroutine LVT_computeMutualInformation(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute 
!   Mutual Information values for the 
!   desired variables. During the first pass, the max min
!   values are computed. During the second pass, the actual
!   mutual information metric is calculated.
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleMutualInformation](\ref{computeSingleMutualInformation})
!     updates the MutualInformation computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     iteration instance through the time series. 
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     MutualInformation computation has been reached
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
       if(LVT_metrics%mi%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call computeSingleMaxMin(alarm,obs,&
                  model, stats, LVT_metrics%mi)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%mi%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call computeSingleMutualInformation(alarm,obs,&
                  model, stats, LVT_metrics%mi)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
       endif
    endif

  end subroutine LVT_ComputeMutualInformation

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
!     MutualInformation computation has been reached
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
                   if(stats%mi(m)%xmaxval(t,k).ne.-1E12.and.&
                        stats%mi(m)%xminval(t,k).ne.1E12) then 

                      stats%mi(m)%xdelta(t,k) = &
                           (stats%mi(m)%xmaxval(t,k) - &
                           stats%mi(m)%xminval(t,k))/MI_nbins
                   else
                      stats%mi(m)%xdelta(t,k) = LVT_rc%udef
                   endif

                   if(stats%mi(m)%ymaxval(t,k).ne.-1E12.and.&
                        stats%mi(m)%yminval(t,k).ne.1E12) then 

                      stats%mi(m)%ydelta(t,k) = &
                           (stats%mi(m)%ymaxval(t,k) - &
                           stats%mi(m)%yminval(t,k))/MI_nbins
                   else
                      stats%mi(m)%ydelta(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
       endif
    endif


  end subroutine computeSingleMaxMin

!BOP
! 
! !ROUTINE: computeSingleMutualInformation
! \label{computeSingleMutualInformation}
!
! !INTERFACE: 
  subroutine computeSingleMutualInformation(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the Mutual Information values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     MutualInformation computation has been reached
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
    real     :: px(MI_nbins)
    real     :: py(MI_nbins)

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(stats%mi(m)%count_total(t,k)&
                        .gt.LVT_rc%obsCountThreshold) then 

                      sumval = 0.0
                      do kk=1,MI_nbins
                         do pp=1,MI_nbins
                            sumval = sumval + & 
                                 stats%mi(m)%pxy(t,k,kk,pp)
                         enddo
                      enddo

                      do kk=1,MI_nbins
                         do pp=1,MI_nbins
                            if(sumval.gt.0) then 
                               stats%mi(m)%pxy(t,k,kk,pp) = & 
                                    stats%mi(m)%pxy(t,k,kk,pp)/&
                                    sumval
                            else
                               stats%mi(m)%pxy(t,k,kk,pp) =LVT_rc%udef
                            endif
                         enddo
                      enddo
                      
                      px = 0 
                      py = 0 
                      do kk=1,MI_nbins
                         px(kk) = px(kk) + sum(stats%mi(m)%pxy(t,k,kk,:))
                      enddo

                      do pp=1,MI_nbins
                         py(pp) = py(pp) + sum(stats%mi(m)%pxy(t,k,:,pp))
                      enddo


                      stats%mi(m)%value_total(t,k) = 0 
                      do kk=1,MI_nbins
                         do pp=1,MI_nbins
                            if(px(kk).ne.LVT_rc%udef.and.&
                                 py(pp).ne.LVT_rc%udef.and.&
                                 stats%mi(m)%pxy(t,k,kk,pp).ne.LVT_rc%udef.and.&
                                 py(pp).ne.0.and.&
                                 stats%mi(m)%pxy(t,k,kk,pp).ne.0) then
                               stats%mi(m)%value_total(t,k)= &
                                    stats%mi(m)%value_total(t,k)+ &
                                    stats%mi(m)%pxy(t,k,kk,pp)*&
                                    log2(stats%mi(m)%pxy(t,k,kk,pp)/&
                                    (px(kk)*py(pp)))
                            endif
                         enddo
                      enddo

                   else
                      stats%mi(m)%value_total(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%mi(m)%value_total(:,k),&
                        LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%mi(m)%value_ci(k))
                enddo
             enddo
          enddo

       endif
    endif


  end subroutine computeSingleMutualInformation

!BOP
! 
! !ROUTINE: LVT_writeMetric_MutualInformation
! \label{LVT_writeMetric_MutualInformation}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_MutualInformation(pass,final,vlevels,stats,obs)
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

    if(pass.eq.LVT_metrics%mi%npass) then
       if(pass.eq.LVT_metrics%mi%npass) then
          if(stats%selectOpt.eq.1) then
             
             allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_total(LVT_rc%ngrid,LVT_rc%nensem))
             
             allocate(value_ci(LVT_rc%nensem))
             
             do k=1,vlevels
                do m=1,LVT_rc%nensem
                   value_total(:,m) = &
                        stats%mi(m)%value_total(:,k)
                   count_value_total(:,m) = & 
                        stats%mi(m)%count_total(:,k)
                   value_ci(m) = stats%mi(m)%value_ci(k)
                   
                enddo
                if(LVT_metrics%mi%selectOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%mi%ftn_total, &
                        value_total(:,:),&
                        stats%vid_total(LVT_MIid,1),k)
                   call LVT_writevar_gridded(LVT_metrics%mi%ftn_total, &
                        real(count_value_total(:,:)),&
                        stats%vid_count_total(LVT_MIid,1),k)
                   
                   
                   call LVT_writeSummaryStats(&
                        LVT_metrics%mi%ftn_summ,&
                        1,&
                        LVT_metrics%mi%short_name,&
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

  end subroutine LVT_writeMetric_MutualInformation

!BOP
! 
! !ROUTINE: LVT_resetMetric_MutualInformation
! \label(LVT_resetMetric_MutualInformation)
!
! !INTERFACE:
  subroutine LVT_resetMetric_MutualInformation(alarm)
! !INPUT PARAMETERS: 
    logical           :: alarm
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  This routine resets required variables to support the
!  temporal computation of MutualInformation values. 
!
! !REVISION HISTORY: 
! 
!EOP


  end subroutine LVT_resetMetric_MutualInformation

!BOP
! 
! !ROUTINE: LVT_writerestart_MutualInformation
! 
! !INTERFACE:
  subroutine LVT_writerestart_MutualInformation(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for MutualInformation metric computations
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
       
       if(LVT_metrics%mi%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   
                   call LVT_writevar_restart(ftn,&
                        stats%mi(m)%value_total(:,k))
                   call LVT_writevar_restart(ftn,&
                        stats%mi(m)%count_total(:,k))
                enddo
             enddo
          end if
       end if
       
       model => model%next
       obs   => obs%next
       stats => stats%next

    end do
             
  end subroutine LVT_writerestart_MutualInformation

!BOP
! 
! !ROUTINE: LVT_readrestart_MutualInformation
! 
! !INTERFACE:
  subroutine LVT_readrestart_MutualInformation(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for MutualInformation metric computations
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
       
       if(LVT_metrics%mi%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   
                   call LVT_readvar_restart(ftn,&
                        stats%mi(m)%value_total(:,k))
                   call LVT_readvar_restart(ftn,&
                        stats%mi(m)%count_total(:,k))
                enddo
             enddo
          end if
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo

             
  end subroutine LVT_readrestart_MutualInformation

end module LVT_MutualInformationMod
