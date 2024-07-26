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
! !MODULE: LVT_JointEntropyMod
! \label(LVT_JointEntropyMod)
!
! !INTERFACE:
module LVT_JointEntropyMod
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
!   This module handles the computation of the Joint Entropy metric. 
!   Joint Entropy values are calculated as: 
!
!      JE(X,Y)=-Sum_x Sum_y P(x,y)log [P(x,y,z)]
!      JE(X,Y,Z)= -Sum_x Sum_y Sum_z P(x,y,z)log [P(x,y,z)]
!   
! !NOTES: 
!   * Temporal/time series calculations are not supported for this metric
!   * This metric should only be used for analyzing one variable at a time
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  27 May 2020    Sujay Kumar  Initial Specification
! 
!EOP


!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initJointEntropy
  public :: LVT_diagnoseJointEntropy
  public :: LVT_computeJointEntropy 
  public :: LVT_writeMetric_JointEntropy
  public :: LVT_resetMetric_JointEntropy
  public :: LVT_writerestart_JointEntropy
  public :: LVT_readrestart_JointEntropy

  integer, parameter :: JE_nbins = 10 
contains
  
!BOP
! 
! !ROUTINE: LVT_initJointEntropy
! \label{LVT_initJointEntropy}
!
! !INTERFACE: 
  subroutine LVT_initJointEntropy(selectNlevs,stats,metric)
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
!  the Joint Entropy metric computations. 
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

    allocate(stats%je(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          if(LVT_rc%nDataStreams.eq.2) then 
             allocate(stats%je(m)%xmaxval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%xminval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%ymaxval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%yminval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%xdelta(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%ydelta(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%pxy(LVT_rc%ngrid, selectNlevs(1),&
                  JE_nbins,JE_nbins))
             
             allocate(stats%je(m)%count_total(LVT_rc%ngrid, selectNlevs(1)))
             allocate(stats%je(m)%value_total(LVT_rc%ngrid, selectNlevs(1)))
             allocate(stats%je(m)%value_ci(selectNlevs(1)))
             
             stats%je(m)%pxy = 0 
             stats%je(m)%value_total = 0 
             stats%je(m)%count_total = 0 
             stats%je(m)%xmaxval = -1E12
             stats%je(m)%xminval = 1E12
             stats%je(m)%ymaxval = -1E12
             stats%je(m)%yminval = 1E12
          elseif(LVT_rc%nDataStreams.eq.3) then 
             allocate(stats%je(m)%xmaxval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%xminval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%ymaxval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%yminval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%zmaxval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%zminval(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%xdelta(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%ydelta(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%zdelta(LVT_rc%ngrid,selectNlevs(1)))
             allocate(stats%je(m)%pxyz(LVT_rc%ngrid, selectNlevs(1),&
                  JE_nbins,JE_nbins,JE_nbins))
             
             allocate(stats%je(m)%count_total(LVT_rc%ngrid, selectNlevs(1)))
             allocate(stats%je(m)%value_total(LVT_rc%ngrid, selectNlevs(1)))
             allocate(stats%je(m)%value_ci(selectNlevs(1)))
             
             stats%je(m)%pxyz = 0 
             stats%je(m)%value_total = 0 
             stats%je(m)%count_total = 0 
             stats%je(m)%xmaxval = -1E12
             stats%je(m)%xminval = 1E12
             stats%je(m)%ymaxval = -1E12
             stats%je(m)%yminval = 1E12
             stats%je(m)%zmaxval = -1E12
             stats%je(m)%zminval = 1E12
          endif
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 2
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initJointEntropy
  
!BOP
! 
! !ROUTINE: LVT_diagnoseJointEntropy
! \label{LVT_diagnoseJointEntropy}
!
! !INTERFACE: 
  subroutine LVT_diagnoseJointEntropy(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Joint Entropy 
!   calculation for desired variables. During the first pass, the 
!   range of the variables (max/min) is captured. During the second 
!   pass, the marginal and joint probabilities are computed. 
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleJointEntropy](\ref{diagnoseSingleJointEntropy})
!     updates the JointEntropy computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none
    integer                 :: pass 
    type(LVT_metadataEntry), pointer :: ds1
    type(LVT_metadataEntry), pointer :: ds2
    type(LVT_metadataEntry), pointer :: ds3
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_rc%nDataStreams.eq.2) then 
       if(pass.eq.1) then
          if(LVT_metrics%je%selectOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(ds1)
             call LVT_getDataStream2Ptr(ds2)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(ds1))
                call diagnoseSingleMaxMin2ds(ds1,&
                     ds2, stats, &
                     LVT_metrics%je)
                
                ds1 => ds1%next
                ds2 => ds2%next
                stats => stats%next
             enddo
          endif
       elseif(pass.eq.2) then 
          if(LVT_metrics%je%selectOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(ds1)
             call LVT_getDataStream2Ptr(ds2)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(ds1))
                call diagnoseSingleJointEntropy2ds(ds1,&
                     ds2, stats, &
                     LVT_metrics%je)
                
                ds1 => ds1%next
                ds2 => ds2%next
                stats => stats%next
             enddo
          endif
       endif
    elseif(LVT_rc%nDataStreams.eq.3) then 
       if(pass.eq.1) then
          if(LVT_metrics%je%selectOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(ds1)
             call LVT_getDataStream2Ptr(ds2)
             call LVT_getDataStream3Ptr(ds3)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(ds1))
                call diagnoseSingleMaxMin3ds(ds1,&
                     ds2, ds3, stats, &
                     LVT_metrics%je)
                
                ds1 => ds1%next
                ds2 => ds2%next
                ds3 => ds3%next
                stats => stats%next
             enddo
          endif
       elseif(pass.eq.2) then 
          if(LVT_metrics%je%selectOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(ds1)
             call LVT_getDataStream2Ptr(ds2)
             call LVT_getDataStream3Ptr(ds3)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(ds1))
                call diagnoseSingleJointEntropy3ds(ds1,&
                     ds2, ds3, stats, &
                     LVT_metrics%je)
                
                ds1 => ds1%next
                ds2 => ds2%next
                ds3 => ds3%next
                stats => stats%next
             enddo
          endif
       endif
    endif

  end subroutine LVT_diagnoseJointEntropy

!BOP
! 
! !ROUTINE: diagnoseSingleMaxMin2ds
! \label{diagnoseSingleMaxMin2ds}
!
! !INTERFACE: 
  subroutine diagnoseSingleMaxMin2ds(ds1,ds2, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the calculations for estimating 
!  max/min values for each grid point when two datastreams
!  are enabled
!
!  The arguments are: 
!
!  \begin{description}
!   \item[ds1] datastream object1
!   \item[ds2] datastream object2
!   \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    type(LVT_metaDataEntry) :: ds1
    type(LVT_metaDataEntry) :: ds2
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k
    integer    :: xbinval, ybinval

    if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,ds1%selectNlevs
                m_k = k+ds1%startNlevs -1
                o_k = k+ds2%startNlevs -1
                if(ds2%count(t,m,o_k).ne.0.and. &
                     ds1%count(t,m,m_k).ne.0) then      
                   if(metric%selectOpt.eq.1) then
                      
                      if(ds1%value(t,m,m_k).gt.stats%je(m)%xmaxval(t,k)) then 
                         stats%je(m)%xmaxval(t,k) = ds1%value(t,m,m_k)
                      endif
                      if(ds2%value(t,m,o_k).gt.stats%je(m)%ymaxval(t,k)) then 
                         stats%je(m)%ymaxval(t,k) = ds2%value(t,m,o_k)
                      endif
                      if(ds1%value(t,m,m_k).lt.stats%je(m)%xminval(t,k)) then 
                         stats%je(m)%xminval(t,k) = ds1%value(t,m,m_k)
                      endif
                      if(ds2%value(t,m,o_k).lt.stats%je(m)%yminval(t,k)) then 
                         stats%je(m)%yminval(t,k) = ds2%value(t,m,o_k)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleMaxMin2ds

!BOP
! 
! !ROUTINE: diagnoseSingleMaxMin3ds
! \label{diagnoseSingleMaxMin3ds}
!
! !INTERFACE: 
  subroutine diagnoseSingleMaxMin3ds(ds1,ds2,ds3,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the calculations for estimating 
!  max/min values for each grid point when three datastreams
!  are enabled
!
!  The arguments are: 
!
!  \begin{description}
!   \item[ds1] datastream object1
!   \item[ds2] datastream object2
!   \item[ds3] datastream object3
!   \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    type(LVT_metaDataEntry) :: ds1
    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds3
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k,t_k
    integer    :: xbinval, ybinval

    if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,ds1%selectNlevs
                m_k = k+ds1%startNlevs -1
                o_k = k+ds2%startNlevs -1
                t_k = k+ds3%startNlevs -1
                if(ds3%count(t,m,t_k).ne.0.and.&
                     ds2%count(t,m,o_k).ne.0.and. &
                     ds1%count(t,m,m_k).ne.0) then      
                   if(metric%selectOpt.eq.1) then
                      
                      if(ds1%value(t,m,m_k).gt.stats%je(m)%xmaxval(t,k)) then 
                         stats%je(m)%xmaxval(t,k) = ds1%value(t,m,m_k)
                      endif
                      if(ds2%value(t,m,o_k).gt.stats%je(m)%ymaxval(t,k)) then 
                         stats%je(m)%ymaxval(t,k) = ds2%value(t,m,o_k)
                      endif
                      if(ds3%value(t,m,t_k).gt.stats%je(m)%zmaxval(t,k)) then 
                         stats%je(m)%zmaxval(t,k) = ds3%value(t,m,t_k)
                      endif
                      if(ds1%value(t,m,m_k).lt.stats%je(m)%xminval(t,k)) then 
                         stats%je(m)%xminval(t,k) = ds1%value(t,m,m_k)
                      endif
                      if(ds2%value(t,m,o_k).lt.stats%je(m)%yminval(t,k)) then 
                         stats%je(m)%yminval(t,k) = ds2%value(t,m,o_k)
                      endif
                      if(ds3%value(t,m,t_k).lt.stats%je(m)%zminval(t,k)) then 
                         stats%je(m)%zminval(t,k) = ds3%value(t,m,t_k)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleMaxMin3ds
!BOP
! 
! !ROUTINE: diagnoseSingleJointEntropy2ds
! \label{diagnoseSingleJointEntropy}
!
! !INTERFACE: 
  subroutine diagnoseSingleJointEntropy2ds(ds1, ds2, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the calculations for computing the 
!  conditional entropy metric for two data streams
!
!  The arguments are: 
!
!  \begin{description}
!   \item[ds1] datastream object 1
!   \item[ds2] datastream object 2
!   \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds1
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k
    integer    :: xbinval, ybinval

    if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,ds1%selectNlevs
                m_k = k+ds1%startNlevs -1
                o_k = k+ds2%startNlevs -1
                if(ds2%count(t,m,o_k).ne.0.and. &
                     ds1%count(t,m,m_k).ne.0) then      
                   if(metric%selectOpt.eq.1.and.&
                        stats%je(m)%xdelta(t,k).gt.0.and.&
                        stats%je(m)%ydelta(t,k).gt.0) then
                      
                      xbinval = nint((ds1%value(t,m,m_k) - &
                           stats%je(m)%xminval(t,k))/&
                           stats%je(m)%xdelta(t,k)) + 1
                      
                      if(xbinval.le.0) xbinval = 1
                      if(xbinval.gt.JE_nbins) xbinval = JE_nbins
                      
                      ybinval = nint((ds2%value(t,m,o_k) - &
                           stats%je(m)%yminval(t,k))/&
                           stats%je(m)%ydelta(t,k)) + 1
                      
                      if(ybinval.le.0) ybinval = 1
                      if(ybinval.gt.JE_nbins) ybinval = JE_nbins
                      
                      stats%je(m)%pxy(t,k,xbinval,ybinval) = &
                           stats%je(m)%pxy(t,k,xbinval,ybinval) + 1
                      
                      stats%je(m)%count_total(t,k) = &
                           stats%je(m)%count_total(t,k) + 1
                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleJointEntropy2ds

!BOP
! 
! !ROUTINE: diagnoseSingleJointEntropy3ds
! \label{diagnoseSingleJointEntropy}
!
! !INTERFACE: 
  subroutine diagnoseSingleJointEntropy3ds(ds1, ds2, ds3, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the calculations for computing the 
!  conditional entropy metric for two data streams
!
!  The arguments are: 
!
!  \begin{description}
!   \item[ds1] datastream object 1
!   \item[ds2] datastream object 2
!   \item[ds3] datastream object 3
!   \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    type(LVT_metaDataEntry) :: ds1
    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds3

    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k,t_k
    integer    :: xbinval, ybinval,zbinval

    if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,ds1%selectNlevs
                m_k = k+ds1%startNlevs -1
                o_k = k+ds2%startNlevs -1
                t_k = k+ds3%startNlevs -1
                if(ds1%count(t,m,m_k).ne.0.and. &
                     ds2%count(t,m,o_k).ne.0.and.&
                     ds3%count(t,m,t_k).ne.0) then      
                   if(metric%selectOpt.eq.1.and.&
                        stats%je(m)%xdelta(t,k).gt.0.and.&
                        stats%je(m)%ydelta(t,k).gt.0.and.&
                        stats%je(m)%zdelta(t,k).gt.0) then
                      
                      xbinval = nint((ds1%value(t,m,m_k) - &
                           stats%je(m)%xminval(t,k))/&
                           stats%je(m)%xdelta(t,k)) + 1
                      
                      if(xbinval.le.0) xbinval = 1
                      if(xbinval.gt.JE_nbins) xbinval = JE_nbins
                      
                      ybinval = nint((ds2%value(t,m,o_k) - &
                           stats%je(m)%yminval(t,k))/&
                           stats%je(m)%ydelta(t,k)) + 1
                      
                      if(ybinval.le.0) ybinval = 1
                      if(ybinval.gt.JE_nbins) ybinval = JE_nbins

                      zbinval = nint((ds3%value(t,m,t_k) - &
                           stats%je(m)%zminval(t,k))/&
                           stats%je(m)%zdelta(t,k)) + 1
                      
                      if(zbinval.le.0) zbinval = 1
                      if(zbinval.gt.JE_nbins) zbinval = JE_nbins
                      
                      stats%je(m)%pxyz(t,k,xbinval,ybinval,zbinval) = &
                           stats%je(m)%pxyz(t,k,xbinval,ybinval,zbinval) + 1
                      
                      stats%je(m)%count_total(t,k) = &
                           stats%je(m)%count_total(t,k) + 1
                   endif
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleJointEntropy3ds

!BOP
! 
! !ROUTINE: LVT_computeJointEntropy
! \label{LVT_computeJointEntropy}
!
! !INTERFACE: 
  subroutine LVT_computeJointEntropy(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute 
!   Joint Entropy values for the 
!   desired variables. During the first pass, the max min
!   values are computed. During the second pass, the actual
!   conditional entropy metric is calculated.
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleJointEntropy](\ref{computeSingleJointEntropy})
!     updates the JointEntropy computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     iteration instance through the time series. 
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     JointEntropy computation has been reached
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

    type(LVT_metadataEntry), pointer :: ds1
    type(LVT_metadataEntry), pointer :: ds2
    type(LVT_metadataEntry), pointer :: ds3
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_rc%nDataStreams.eq.2) then 
       if(pass.eq.1) then
          if(LVT_metrics%je%selectOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(ds1)
             call LVT_getDataStream2Ptr(ds2)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(ds1))
                
                call computeSingleMaxMin2ds(alarm,ds1,&
                     ds2, stats, LVT_metrics%je)
                
                ds1 => ds1%next
                ds2 => ds2%next
                stats => stats%next
             enddo
             
          endif
       elseif(pass.eq.2) then 
          if(LVT_metrics%je%selectOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(ds1)
             call LVT_getDataStream2Ptr(ds2)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(ds1))
                
                call computeSingleJointEntropy2ds(alarm,ds1,&
                     ds2, stats, LVT_metrics%je)
                
                ds1 => ds1%next
                ds2 => ds2%next
                stats => stats%next
             enddo
             
          endif
       endif
    elseif(LVT_rc%nDataStreams.eq.3) then 
       if(pass.eq.1) then
          if(LVT_metrics%je%selectOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(ds1)
             call LVT_getDataStream2Ptr(ds2)
             call LVT_getDataStream3Ptr(ds3)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(ds1))
                
                call computeSingleMaxMin3ds(alarm,ds1,&
                     ds2, ds3, stats, LVT_metrics%je)
                
                ds1 => ds1%next
                ds2 => ds2%next
                ds3 => ds3%next
                stats => stats%next
             enddo
             
          endif
       elseif(pass.eq.2) then 
          if(LVT_metrics%je%selectOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(ds1)
             call LVT_getDataStream2Ptr(ds2)
             call LVT_getDataStream3Ptr(ds3)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(ds1))
                
                call computeSingleJointEntropy3ds(alarm,ds1,&
                     ds2, ds3, stats, LVT_metrics%je)
                
                ds1 => ds1%next
                ds2 => ds2%next
                ds3 => ds3%next
                stats => stats%next
             enddo
             
          endif
       endif

    endif
  end subroutine LVT_ComputeJointEntropy

!BOP
! 
! !ROUTINE: computeSingleMaxMin2ds
! \label{computeSingleMaxMin2ds}
!
! !INTERFACE: 
  subroutine computeSingleMaxMin2ds(alarm,ds1, ds2,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes max/min values for each grid point when
!  two data streams are present
! 
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     JointEntropy computation has been reached
!    \item[ds2] datastream 1 object
!    \item[ds1] datastream 2 object 
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
    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds1
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer  :: t,l,m,k,tind
    integer  :: kk,pp


    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,ds1%selectNlevs
                   if(stats%je(m)%xmaxval(t,k).ne.-1E12.and.&
                        stats%je(m)%xminval(t,k).ne.1E12) then 

                      stats%je(m)%xdelta(t,k) = &
                           (stats%je(m)%xmaxval(t,k) - &
                           stats%je(m)%xminval(t,k))/JE_nbins
                   else
                      stats%je(m)%xdelta(t,k) = LVT_rc%udef
                   endif

                   if(stats%je(m)%ymaxval(t,k).ne.-1E12.and.&
                        stats%je(m)%yminval(t,k).ne.1E12) then 

                      stats%je(m)%ydelta(t,k) = &
                           (stats%je(m)%ymaxval(t,k) - &
                           stats%je(m)%yminval(t,k))/JE_nbins
                   else
                      stats%je(m)%ydelta(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
       endif
    endif


  end subroutine computeSingleMaxMin2ds

!BOP
! 
! !ROUTINE: computeSingleMaxMin3ds
! \label{computeSingleMaxMin3ds}
!
! !INTERFACE: 
  subroutine computeSingleMaxMin3ds(alarm,ds1, ds2,ds3,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes max/min values for each grid point when
!  three data streams are present
! 
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     JointEntropy computation has been reached
!    \item[ds2] datastream 1 object
!    \item[ds1] datastream 2 object 
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
    type(LVT_metaDataEntry) :: ds1
    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds3
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer  :: t,l,m,k,tind
    integer  :: kk,pp


    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,ds1%selectNlevs
                   if(stats%je(m)%xmaxval(t,k).ne.-1E12.and.&
                        stats%je(m)%xminval(t,k).ne.1E12) then 

                      stats%je(m)%xdelta(t,k) = &
                           (stats%je(m)%xmaxval(t,k) - &
                           stats%je(m)%xminval(t,k))/JE_nbins
                   else
                      stats%je(m)%xdelta(t,k) = LVT_rc%udef
                   endif

                   if(stats%je(m)%ymaxval(t,k).ne.-1E12.and.&
                        stats%je(m)%yminval(t,k).ne.1E12) then 

                      stats%je(m)%ydelta(t,k) = &
                           (stats%je(m)%ymaxval(t,k) - &
                           stats%je(m)%yminval(t,k))/JE_nbins
                   else
                      stats%je(m)%ydelta(t,k) = LVT_rc%udef
                   endif

                   if(stats%je(m)%zmaxval(t,k).ne.-1E12.and.&
                        stats%je(m)%zminval(t,k).ne.1E12) then 

                      stats%je(m)%zdelta(t,k) = &
                           (stats%je(m)%zmaxval(t,k) - &
                           stats%je(m)%zminval(t,k))/JE_nbins
                   else
                      stats%je(m)%zdelta(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
       endif
    endif


  end subroutine computeSingleMaxMin3ds
  
!BOP
! 
! !ROUTINE: computeSingleJointEntropy2ds
! \label{computeSingleJointEntropy2ds}
!
! !INTERFACE: 
  subroutine computeSingleJointEntropy2ds(alarm,ds1, ds2,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the Joint Entropy values for a single variable
!  when two datastreams are present
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     JointEntropy computation has been reached
!    \item[ds1] datastream 1 object
!    \item[ds2] datastream 2 object
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
    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds1
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer  :: t,l,m,k,tind
    integer  :: kk,pp
    real     :: sumval

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,ds1%selectNlevs
                   if(stats%je(m)%count_total(t,k)&
                        .gt.LVT_rc%obsCountThreshold) then 

                      sumval = 0.0
                      do kk=1,JE_nbins
                         do pp=1,JE_nbins
                            sumval = sumval + & 
                                 stats%je(m)%pxy(t,k,kk,pp)
                         enddo
                      enddo

                      do kk=1,JE_nbins
                         do pp=1,JE_nbins
                            if(sumval.gt.0) then 
                               stats%je(m)%pxy(t,k,kk,pp) = & 
                                    stats%je(m)%pxy(t,k,kk,pp)/&
                                    sumval
                            else
                               stats%je(m)%pxy(t,k,kk,pp) =LVT_rc%udef
                            endif
                         enddo
                      enddo

                      stats%je(m)%value_total(t,k) = 0 
                      do kk=1,JE_nbins
                         do pp=1,JE_nbins
                            if(stats%je(m)%pxy(t,k,kk,pp).ne.LVT_rc%udef.and.&
                                 stats%je(m)%pxy(t,k,kk,pp).ne.0) then
                               stats%je(m)%value_total(t,k)= &
                                    stats%je(m)%value_total(t,k)- &
                                    stats%je(m)%pxy(t,k,kk,pp)*&
                                    log2(stats%je(m)%pxy(t,k,kk,pp))
                            endif
                         enddo
                      enddo

                   else
                      stats%je(m)%value_total(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          
          do m=1,LVT_rc%nensem
             do k=1,ds1%selectNlevs
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%je(m)%value_total(:,k),&
                        LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%je(m)%value_ci(k))
                enddo
             enddo
          enddo

       endif
    endif


  end subroutine computeSingleJointEntropy2ds

  
!BOP
! 
! !ROUTINE: computeSingleJointEntropy3ds
! \label{computeSingleJointEntropy3ds}
!
! !INTERFACE: 
  subroutine computeSingleJointEntropy3ds(alarm,ds1,ds2,ds3,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the Joint Entropy values for a single variable
!  when three datastreams are present
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     JointEntropy computation has been reached
!    \item[ds1] datastream 1 object
!    \item[ds2] datastream 2 object
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
    type(LVT_metaDataEntry) :: ds1
    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds3
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer  :: t,l,m,k,tind
    integer  :: kk,pp,tt
    real     :: sumval

    if(LVT_rc%nDataStreams.eq.2) then 
       if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem
                   do k=1,ds1%selectNlevs
                      if(stats%je(m)%count_total(t,k)&
                           .gt.LVT_rc%obsCountThreshold) then 
                         
                         sumval = 0.0
                         do kk=1,JE_nbins
                            do pp=1,JE_nbins
                               sumval = sumval + & 
                                    stats%je(m)%pxy(t,k,kk,pp)
                            enddo
                         enddo
                         
                         do kk=1,JE_nbins
                            do pp=1,JE_nbins
                               if(sumval.gt.0) then 
                                  stats%je(m)%pxy(t,k,kk,pp) = & 
                                       stats%je(m)%pxy(t,k,kk,pp)/&
                                       sumval
                               else
                                  stats%je(m)%pxy(t,k,kk,pp) =LVT_rc%udef
                               endif
                            enddo
                         enddo

                         stats%je(m)%value_total(t,k) = 0 
                         do kk=1,JE_nbins
                            do pp=1,JE_nbins
                               if(stats%je(m)%pxy(t,k,kk,pp).ne.LVT_rc%udef.and.&
                                    stats%je(m)%pxy(t,k,kk,pp).ne.0) then
                                  stats%je(m)%value_total(t,k)= &
                                       stats%je(m)%value_total(t,k)- &
                                       stats%je(m)%pxy(t,k,kk,pp)*&
                                       log2(stats%je(m)%pxy(t,k,kk,pp))
                               endif
                            enddo
                         enddo

                      else
                         stats%je(m)%value_total(t,k) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             do m=1,LVT_rc%nensem
                do k=1,ds1%selectNlevs
                   do l=1, LVT_rc%strat_nlevels
                      call LVT_computeCI(stats%je(m)%value_total(:,k),&
                           LVT_rc%ngrid,&
                           LVT_rc%pval_CI,stats%je(m)%value_ci(k))
                   enddo
                enddo
             enddo

          endif
       endif
    elseif(LVT_rc%nDataStreams.eq.3) then 
       if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem
                   do k=1,ds1%selectNlevs
                      if(stats%je(m)%count_total(t,k)&
                           .gt.LVT_rc%obsCountThreshold) then 
                         
                         sumval = 0.0
                         do kk=1,JE_nbins
                            do pp=1,JE_nbins
                               do tt=1,JE_nbins
                                  sumval = sumval + & 
                                       stats%je(m)%pxyz(t,k,kk,pp,tt)
                               enddo
                            enddo
                         enddo
                         do kk=1,JE_nbins
                            do pp=1,JE_nbins
                               do tt=1,JE_nbins
                                  if(sumval.gt.0) then 
                                     stats%je(m)%pxyz(t,k,kk,pp,tt) = & 
                                          stats%je(m)%pxyz(t,k,kk,pp,tt)/&
                                          sumval
                                  else
                                     stats%je(m)%pxyz(t,k,kk,pp,tt) =LVT_rc%udef
                                  endif
                               enddo
                            enddo
                         enddo

                         stats%je(m)%value_total(t,k) = 0 
                         do kk=1,JE_nbins
                            do pp=1,JE_nbins
                               do tt=1,JE_nbins
                                  if(stats%je(m)%pxyz(t,k,kk,pp,tt).ne.LVT_rc%udef.and.&
                                       stats%je(m)%pxyz(t,k,kk,pp,tt).ne.0) then
                                     stats%je(m)%value_total(t,k)= &
                                          stats%je(m)%value_total(t,k)- &
                                          stats%je(m)%pxyz(t,k,kk,pp,tt)*&
                                          log2(stats%je(m)%pxyz(t,k,kk,pp,tt))
                                  endif
                               enddo
                            enddo
                         enddo
                      else
                         stats%je(m)%value_total(t,k) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo

             do m=1,LVT_rc%nensem
                do k=1,ds1%selectNlevs
                   do l=1, LVT_rc%strat_nlevels
                      call LVT_computeCI(stats%je(m)%value_total(:,k),&
                           LVT_rc%ngrid,&
                           LVT_rc%pval_CI,stats%je(m)%value_ci(k))
                   enddo
                enddo
             enddo

          endif
       endif
    endif

  end subroutine computeSingleJointEntropy3ds

!BOP
! 
! !ROUTINE: LVT_writeMetric_JointEntropy
! \label{LVT_writeMetric_JointEntropy}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_JointEntropy(pass,final,vlevels,stats,ds2)
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
    type(LVT_metaDataEntry) :: ds2

    integer                 :: k,l,m,tind

    real,    allocatable    :: value_total(:,:)
    integer, allocatable    :: count_value_total(:,:)
    real,    allocatable    :: value_ci(:)

    real,    allocatable    :: value_ts(:,:)
    integer, allocatable    :: count_value_ts(:,:)

    real,    allocatable    :: value_adc(:,:,:)
    real,    allocatable    :: value_asc(:,:,:)

    if(pass.eq.LVT_metrics%je%npass) then
       if(pass.eq.LVT_metrics%je%npass) then
          if(stats%selectOpt.eq.1) then
             
             allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_total(LVT_rc%ngrid,LVT_rc%nensem))
             
             allocate(value_ci(LVT_rc%nensem))
             
             do k=1,vlevels
                do m=1,LVT_rc%nensem
                   value_total(:,m) = &
                        stats%je(m)%value_total(:,k)
                   count_value_total(:,m) = & 
                        stats%je(m)%count_total(:,k)
                   value_ci(m) = stats%je(m)%value_ci(k)
                   
                enddo
                if(LVT_metrics%je%selectOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%je%ftn_total, &
                        value_total(:,:),&
                        stats%vid_total(LVT_JEid,1),k)
                   call LVT_writevar_gridded(LVT_metrics%je%ftn_total, &
                        real(count_value_total(:,:)),&
                        stats%vid_count_total(LVT_JEid,1),k)
                   
                   
                   call LVT_writeSummaryStats(&
                        LVT_metrics%je%ftn_summ,&
                        1,&
                        LVT_metrics%je%short_name,&
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

  end subroutine LVT_writeMetric_JointEntropy

!BOP
! 
! !ROUTINE: LVT_resetMetric_JointEntropy
! \label(LVT_resetMetric_JointEntropy)
!
! !INTERFACE:
  subroutine LVT_resetMetric_JointEntropy(alarm)
! !INPUT PARAMETERS: 
    logical           :: alarm
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  This routine resets required variables to support the
!  temporal computation of JointEntropy values. 
!
! !REVISION HISTORY: 
! 
!EOP


  end subroutine LVT_resetMetric_JointEntropy

!BOP
! 
! !ROUTINE: LVT_writerestart_JointEntropy
! 
! !INTERFACE:
  subroutine LVT_writerestart_JointEntropy(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for JointEntropy metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    integer              :: k,l,m
    type(LVT_metaDataEntry), pointer :: ds1
    type(LVT_metaDataEntry), pointer :: ds2
    type(LVT_statsEntry),    pointer :: stats

    call LVT_getDataStream1Ptr(ds1)
    call LVT_getDataStream2Ptr(ds2)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(ds1))
       
       if(LVT_metrics%je%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,ds1%selectNlevs
                   
                   call LVT_writevar_restart(ftn,&
                        stats%je(m)%value_total(:,k))
                   call LVT_writevar_restart(ftn,&
                        stats%je(m)%count_total(:,k))
                enddo
             enddo
          end if
       end if
       
       ds1 => ds1%next
       ds2   => ds2%next
       stats => stats%next

    end do
             
  end subroutine LVT_writerestart_JointEntropy

!BOP
! 
! !ROUTINE: LVT_readrestart_JointEntropy
! 
! !INTERFACE:
  subroutine LVT_readrestart_JointEntropy(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for JointEntropy metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    type(LVT_metaDataEntry), pointer :: ds1
    type(LVT_metaDataEntry), pointer :: ds2
    type(LVT_statsEntry),    pointer :: stats
    integer              :: k,l,m


    call LVT_getDataStream1Ptr(ds1)
    call LVT_getDataStream2Ptr(ds2)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(ds1))
       
       if(LVT_metrics%je%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,ds1%selectNlevs
                   
                   call LVT_readvar_restart(ftn,&
                        stats%je(m)%value_total(:,k))
                   call LVT_readvar_restart(ftn,&
                        stats%je(m)%count_total(:,k))
                enddo
             enddo
          end if
       endif

       ds1 => ds1%next
       ds2   => ds2%next
       stats => stats%next

    enddo

             
  end subroutine LVT_readrestart_JointEntropy

end module LVT_JointEntropyMod
