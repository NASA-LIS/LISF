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
! !MODULE: LVT_InformationEntropyMod
! \label(LVT_InformationEntropyMod)
!
! !INTERFACE:
module LVT_InformationEntropyMod
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
!   This module handles the computation of the Information Entropy metric. 
!   Information Entropy values are calculated as: 
!
!      H(X)= -Sum_x  P(x)log [P(x)]
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
  public :: LVT_initInformationEntropy
  public :: LVT_diagnoseInformationEntropy
  public :: LVT_computeInformationEntropy 
  public :: LVT_writeMetric_InformationEntropy
  public :: LVT_resetMetric_InformationEntropy
  public :: LVT_writerestart_InformationEntropy
  public :: LVT_readrestart_InformationEntropy

  integer, parameter :: IE_nbins = 10 
contains
  
!BOP
! 
! !ROUTINE: LVT_initInformationEntropy
! \label{LVT_initInformationEntropy}
!
! !INTERFACE: 
  subroutine LVT_initInformationEntropy(selectNlevs,stats,metric)
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

    allocate(stats%ie(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%ie(m)%xmaxval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ie(m)%xminval(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ie(m)%xdelta(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ie(m)%px(LVT_rc%ngrid, selectNlevs(1),&
               IE_nbins))
          allocate(stats%ie(m)%model_count_total(LVT_rc%ngrid, selectNlevs(1)))
          allocate(stats%ie(m)%model_value_total(LVT_rc%ngrid, selectNlevs(1)))
          allocate(stats%ie(m)%model_value_ci(selectNlevs(1)))

          stats%ie(m)%px = 0 
          stats%ie(m)%model_value_total = 0 
          stats%ie(m)%model_count_total = 0 
          stats%ie(m)%xmaxval = -1E12
          stats%ie(m)%xminval = 1E12

          if(LVT_rc%obssource(2).ne."none") then
             if(selectNlevs(2).ge.1) then
                allocate(stats%ie(m)%ymaxval(LVT_rc%ngrid,selectNlevs(2)))
                allocate(stats%ie(m)%yminval(LVT_rc%ngrid,selectNlevs(2)))
                allocate(stats%ie(m)%ydelta(LVT_rc%ngrid,selectNlevs(2)))
                allocate(stats%ie(m)%py(LVT_rc%ngrid, selectNlevs(2),&
                     IE_nbins))
                allocate(stats%ie(m)%obs_count_total(LVT_rc%ngrid, &
                     selectNlevs(2)))
                allocate(stats%ie(m)%obs_value_total(LVT_rc%ngrid, &
                     selectNlevs(2)))
                allocate(stats%ie(m)%obs_value_ci(selectNlevs(2)))
                stats%ie(m)%py = 0 
                
                stats%ie(m)%ymaxval = -1E12
                stats%ie(m)%yminval = 1E12
             endif
          endif
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    if(LVT_rc%obssource(2).ne."none") then
       metric%obsData = .true.
    else
       metric%obsData = .false.
    endif
    
    metric%stdevFlag = .false. 

  end subroutine LVT_initInformationEntropy
  
!BOP
! 
! !ROUTINE: LVT_diagnoseInformationEntropy
! \label{LVT_diagnoseInformationEntropy}
!
! !INTERFACE: 
  subroutine LVT_diagnoseInformationEntropy(pass)
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
!    \item[diagnoseSingleInformationEntropy](\ref{diagnoseSingleInformationEntropy})
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
       if(LVT_metrics%ie%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMaxMin(obs,&
                  model, stats, &
                  LVT_metrics%ie)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%ie%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleInformationEntropy(obs,&
                  model, stats, &
                  LVT_metrics%ie)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseInformationEntropy

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

    if(stats%selectOpt.eq.1.and.model%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                if(model%count(t,m,m_k).ne.0) then      
                   if(metric%selectOpt.eq.1) then
                      
                      if(model%value(t,m,m_k).gt.stats%ie(m)%xmaxval(t,k)) then 
                         stats%ie(m)%xmaxval(t,k) = model%value(t,m,m_k)
                      endif
                      if(model%value(t,m,m_k).lt.stats%ie(m)%xminval(t,k)) then 
                         stats%ie(m)%xminval(t,k) = model%value(t,m,m_k)
                      endif
                   endif
                endif                
             enddo
          enddo
          do k=1,obs%selectNlevs
             do m=1,LVT_rc%nensem
                o_k = k+obs%startNlevs -1
                if(LVT_rc%obssource(2).ne."none") then
                   if(obs%selectNlevs.ge.1) then
                      
                      if(obs%value(t,m,o_k).gt.stats%ie(m)%ymaxval(t,k)) then 
                         stats%ie(m)%ymaxval(t,k) = obs%value(t,m,o_k)
                      endif
                      if(obs%value(t,m,o_k).lt.stats%ie(m)%yminval(t,k)) then 
                         stats%ie(m)%yminval(t,k) = obs%value(t,m,o_k)
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
! !ROUTINE: diagnoseSingleInformationEntropy
! \label{diagnoseSingleInformationEntropy}
!
! !INTERFACE: 
  subroutine diagnoseSingleInformationEntropy(obs, model, stats,metric)
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
                if(model%count(t,m,m_k).ne.0) then      
                   if(metric%selectOpt.eq.1.and.&
                        stats%ie(m)%xdelta(t,k).gt.0) then 
                      xbinval = nint((model%value(t,m,m_k) - &
                           stats%ie(m)%xminval(t,k))/&
                           stats%ie(m)%xdelta(t,k)) + 1
                      
                      if(xbinval.le.0) xbinval = 1
                      if(xbinval.gt.IE_nbins) xbinval = IE_nbins
                      
                      stats%ie(m)%px(t,k,xbinval) = &
                           stats%ie(m)%px(t,k,xbinval) + 1
                      stats%ie(m)%model_count_total(t,k) = & 
                           stats%ie(m)%model_count_total(t,k) + 1

                   endif
                endif
             enddo
             do k=1,obs%selectNlevs
                o_k = k+obs%startNlevs -1
                if(obs%count(t,m,o_k).ne.0) then      
                   if(metric%selectOpt.eq.1.and.&
                        stats%ie(m)%ydelta(t,k).gt.0) then 
                      ybinval = nint((obs%value(t,m,o_k) - &
                           stats%ie(m)%yminval(t,k))/&
                           stats%ie(m)%ydelta(t,k)) + 1
                      
                      if(ybinval.le.0) ybinval = 1
                      if(ybinval.gt.IE_nbins) ybinval = IE_nbins
                      
                      stats%ie(m)%py(t,k,ybinval) = &
                           stats%ie(m)%py(t,k,ybinval) + 1
                      stats%ie(m)%obs_count_total(t,k) = & 
                           stats%ie(m)%obs_count_total(t,k) + 1

                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleInformationEntropy


!BOP
! 
! !ROUTINE: LVT_computeInformationEntropy
! \label{LVT_computeInformationEntropy}
!
! !INTERFACE: 
  subroutine LVT_computeInformationEntropy(pass,alarm)
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
!    \item[computeSingleInformationEntropy](\ref{computeSingleInformationEntropy})
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
       if(LVT_metrics%ie%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call computeSingleMaxMin(alarm,obs,&
                  model, stats, LVT_metrics%ie)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
       endif
    elseif(pass.eq.2) then 
       if(LVT_metrics%ie%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call computeSingleInformationEntropy(alarm,obs,&
                  model, stats, LVT_metrics%ie)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
       endif
    endif

  end subroutine LVT_ComputeInformationEntropy

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
!     InformationEntropy computation has been reached
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
       if(stats%selectOpt.eq.1.and.model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(stats%ie(m)%xmaxval(t,k).ne.-1E12.and.&
                        stats%ie(m)%xminval(t,k).ne.1E12) then 

                      stats%ie(m)%xdelta(t,k) = &
                           (stats%ie(m)%xmaxval(t,k) - &
                           stats%ie(m)%xminval(t,k))/IE_nbins
                   else
                      stats%ie(m)%xdelta(t,k) = LVT_rc%udef
                   endif
                end do
             enddo
          enddo
       endif
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1.and.&
            LVT_rc%obssource(2).ne."none") then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,obs%selectNlevs
                   if(stats%ie(m)%ymaxval(t,k).ne.-1E12.and.&
                        stats%ie(m)%yminval(t,k).ne.1E12) then 
                      
                      stats%ie(m)%ydelta(t,k) = &
                           (stats%ie(m)%ymaxval(t,k) - &
                           stats%ie(m)%yminval(t,k))/IE_nbins
                   else
                      stats%ie(m)%ydelta(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
       endif
    endif


  end subroutine computeSingleMaxMin

!BOP
! 
! !ROUTINE: computeSingleInformationEntropy
! \label{computeSingleInformationEntropy}
!
! !INTERFACE: 
  subroutine computeSingleInformationEntropy(alarm,obs, model,stats,metric)
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
!     InformationEntropy computation has been reached
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
    real     :: px(IE_nbins)
    real     :: py(IE_nbins)

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(stats%ie(m)%model_count_total(t,k)&
                        .gt.LVT_rc%obsCountThreshold) then 

                      sumval = 0.0
                      do kk=1,IE_nbins
                         sumval = sumval + & 
                              stats%ie(m)%px(t,k,kk)
                      enddo

                      do kk=1,IE_nbins
                         if(sumval.gt.0) then 
                            stats%ie(m)%px(t,k,kk) = & 
                                 stats%ie(m)%px(t,k,kk)/&
                                 sumval
                         else
                            stats%ie(m)%px(t,k,kk) =LVT_rc%udef
                         endif
                      enddo
                      
                      stats%ie(m)%model_value_total(t,k) = 0 
                      do kk=1,IE_nbins
                         if(stats%ie(m)%px(t,k,kk).ne.LVT_rc%udef.and.&
                              stats%ie(m)%px(t,k,kk).ne.0) then 
                            stats%ie(m)%model_value_total(t,k) = & 
                                 stats%ie(m)%model_value_total(t,k) - &
                                 stats%ie(m)%px(t,k,kk)*log2(&
                                 stats%ie(m)%px(t,k,kk))
                         endif
                      enddo
                   else
                      stats%ie(m)%model_value_total(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo

          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                call LVT_computeCI(stats%ie(m)%model_value_total(:,k),&
                     LVT_rc%ngrid,&
                     LVT_rc%pval_CI,stats%ie(m)%model_value_ci(k))
             enddo
          enddo
       endif

       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1.and.&
            LVT_rc%obssource(2).ne."none") then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,obs%selectNlevs
                   if(stats%ie(m)%obs_count_total(t,k)&
                        .gt.LVT_rc%obsCountThreshold) then 

                      sumval = 0.0
                      do kk=1,IE_nbins
                         sumval = sumval + & 
                              stats%ie(m)%py(t,k,kk)
                      enddo

                      do kk=1,IE_nbins
                         if(sumval.gt.0) then 
                            stats%ie(m)%py(t,k,kk) = & 
                                 stats%ie(m)%py(t,k,kk)/&
                                 sumval
                         else
                            stats%ie(m)%py(t,k,kk) =LVT_rc%udef
                         endif
                      enddo
                      
                      stats%ie(m)%obs_value_total(t,k) = 0 
                      do kk=1,IE_nbins
                         if(stats%ie(m)%py(t,k,kk).ne.LVT_rc%udef.and.&
                              stats%ie(m)%py(t,k,kk).ne.0) then 
                            stats%ie(m)%obs_value_total(t,k) = & 
                                 stats%ie(m)%obs_value_total(t,k) - &
                                 stats%ie(m)%py(t,k,kk)*log2(&
                                 stats%ie(m)%py(t,k,kk))
                         endif
                      enddo
                   else
                      stats%ie(m)%obs_value_total(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          do m=1,LVT_rc%nensem
            do k=1,obs%selectNlevs
               call LVT_computeCI(stats%ie(m)%obs_value_total(:,k),&
                    LVT_rc%ngrid,&
                    LVT_rc%pval_CI,stats%ie(m)%obs_value_ci(k))
            enddo
         enddo

       endif
    endif


  end subroutine computeSingleInformationEntropy

!BOP
! 
! !ROUTINE: LVT_writeMetric_InformationEntropy
! \label{LVT_writeMetric_InformationEntropy}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_InformationEntropy(pass,final,vlevels,stats,obs)
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

    real,    allocatable    :: model_value_total(:,:)
    integer, allocatable    :: count_model_value_total(:,:)
    real,    allocatable    :: model_value_ci(:)

    real,    allocatable    :: obs_value_total(:,:)
    integer, allocatable    :: count_obs_value_total(:,:)
    real,    allocatable    :: obs_value_ci(:)


    real,    allocatable    :: value_adc(:,:,:)
    real,    allocatable    :: value_asc(:,:,:)

    if(pass.eq.LVT_metrics%ie%npass) then
       if(pass.eq.LVT_metrics%ie%npass) then
          if(stats%selectOpt.eq.1) then
             
             allocate(model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
             
             allocate(model_value_ci(LVT_rc%nensem))
             
             if(LVT_rc%obssource(2).ne."none".and.&
                  obs%selectNlevs.ge.1) then 
                
                allocate(obs_value_total(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(count_obs_value_total(LVT_rc%ngrid,LVT_rc%nensem))
                
                allocate(obs_value_ci(LVT_rc%nensem))
             endif

             do k=1,vlevels
                do m=1,LVT_rc%nensem
                   model_value_total(:,m) = &
                        stats%ie(m)%model_value_total(:,k)
                   count_model_value_total(:,m) = & 
                        stats%ie(m)%model_count_total(:,k)
                   model_value_ci(m) = stats%ie(m)%model_value_ci(k)
                   
                enddo

                if(LVT_rc%obssource(2).ne."none".and.&
                     obs%selectNlevs.ge.1) then 
                   do m=1,LVT_rc%nensem
                      obs_value_total(:,m) = &
                           stats%ie(m)%obs_value_total(:,k)
                      count_obs_value_total(:,m) = & 
                           stats%ie(m)%obs_count_total(:,k)
                      obs_value_ci(m) = stats%ie(m)%obs_value_ci(k)
                      
                   enddo
                endif

                if(LVT_metrics%ie%selectOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%ie%ftn_total, &
                        model_value_total(:,:),&
                        stats%vid_total(LVT_IEid,1),k)
                   call LVT_writevar_gridded(LVT_metrics%ie%ftn_total, &
                        real(count_model_value_total(:,:)),&
                        stats%vid_count_total(LVT_IEid,1),k)
                   
                   if(LVT_rc%obssource(2).ne."none".and.&
                        obs%selectNlevs.ge.1) then 


                      call LVT_writevar_gridded(LVT_metrics%ie%ftn_total, &
                        obs_value_total(:,:),&
                        stats%vid_total(LVT_IEid,2),k)
                      call LVT_writevar_gridded(LVT_metrics%ie%ftn_total, &
                           real(count_obs_value_total(:,:)),&
                           stats%vid_count_total(LVT_IEid,2),k)
                   endif

                   call LVT_writeSummaryStats(&
                        LVT_metrics%ie%ftn_summ,&
                        1,&
                        LVT_metrics%ie%short_name,&
                        LVT_rc%ngrid,&
                        model_value_total(:,:), &
                        count_model_value_total(:,:),&
                        stats%standard_name,&
                        model_value_ci(:))
                   
                   if(LVT_rc%obssource(2).ne."none".and.&
                        obs%selectNlevs.ge.1) then 
                      call LVT_writeSummaryStats(&
                           LVT_metrics%ie%ftn_summ,&
                           1,&
                           LVT_metrics%ie%short_name,&
                           LVT_rc%ngrid,&
                           obs_value_total(:,:), &
                           count_obs_value_total(:,:),&
                           "DS2_"//trim(stats%standard_name),&
                           obs_value_ci(:))
                      
                   endif
                endif
             enddo
             deallocate(model_value_total)
             deallocate(count_model_value_total)
             deallocate(model_value_ci)

             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then
                deallocate(obs_value_total)
                deallocate(count_obs_value_total)
                deallocate(obs_value_ci)
             endif

          endif
       endif
    endif

  end subroutine LVT_writeMetric_InformationEntropy

!BOP
! 
! !ROUTINE: LVT_resetMetric_InformationEntropy
! \label(LVT_resetMetric_InformationEntropy)
!
! !INTERFACE:
  subroutine LVT_resetMetric_InformationEntropy(alarm)
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


  end subroutine LVT_resetMetric_InformationEntropy

!BOP
! 
! !ROUTINE: LVT_writerestart_InformationEntropy
! 
! !INTERFACE:
  subroutine LVT_writerestart_InformationEntropy(ftn,pass)
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

             
  end subroutine LVT_writerestart_InformationEntropy

!BOP
! 
! !ROUTINE: LVT_readrestart_InformationEntropy
! 
! !INTERFACE:
  subroutine LVT_readrestart_InformationEntropy(ftn)
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

  end subroutine LVT_readrestart_InformationEntropy


end module LVT_InformationEntropyMod
