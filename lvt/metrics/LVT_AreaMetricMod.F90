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
! !MODULE: LVT_AreaMetricMod
! !\label(LVT_AreaMetricMod)
!
! !INTERFACE:
module LVT_AreaMetricMod
! 
! !USES:   
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
!  !DESCRIPTION: 
!   This module computes the area (domain-wide sum) of 
!   the model variables and the observation variables. 
!   By definition, this is a time-series computation. 
! 
!  !NOTES: Average diurnal and seasonal cycle computations
!  are not currently supported. 
! 
! !FILES USED:
!
!  !REVISION HISTORY: 
!  27 Apr 2011    Sujay Kumar  Initial Specification
! 
!EOP

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initAreaMetric
  public :: LVT_diagnoseAreaMetric
  public :: LVT_computeAreaMetric
  public :: LVT_writeMetric_AreaMetric
  public :: LVT_resetMetric_AreaMetric
  public :: LVT_writerestart_AreaMetric
  public :: LVT_readrestart_AreaMetric
  private
  
contains

!BOP
! 
! !ROUTINE: LVT_initAreaMetric
! \label{LVT_initAreaMetric}
!
! !INTERFACE: 
  subroutine LVT_initAreaMetric(selectNlevs,stats,metric)
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
!  the area computations. 
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
    
    allocate(stats%area(LVT_rc%nensem))

    if(metric%selectOpt.eq.1.or.&
         metric%timeOpt.eq.1) then 
       allocate(stats%area(m)%model_area(selectNlevs(1)))
       allocate(stats%area(m)%obs_area(selectNlevs(1)))
       
       stats%area(m)%model_area = 0.0
       stats%area(m)%obs_area = 0.0 
    endif
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1
    metric%obsData = .false. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initAreaMetric

!BOP
! 
! !ROUTINE: LVT_diagnoseAreaMetric
! \label{LVT_diagnoseAreaMetric}
!
! !INTERFACE: 
  subroutine LVT_diagnoseAreaMetric(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Area calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleAreaMetric](\ref{diagnoseSingleAreaMetric})
!     updates the AreaMetric computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
    implicit none
    integer        :: pass
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%area%selectOpt.eq.1.or.&
            LVT_metrics%area%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleAreaMetric(obs,model,stats,&
                  LVT_metrics%area)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif

  end subroutine LVT_diagnoseAreaMetric

!BOP
! 
! !ROUTINE: diagnoseSingleAreaMetric
! \label{diagnoseSingleAreaMetric}
!
! !INTERFACE: 
  subroutine diagnoseSingleAreaMetric(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the AreaMetric computation (updates the running 
!  sum calculations )
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
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                
                if(trim(obs%units).eq.trim(model%units)) then 
                   if(metric%selectOpt.eq.1.or.metric%timeOpt.eq.1) then
                      if(obs%count(t,m,o_k).ne.0) then 
                         stats%area(m)%obs_area(k) = stats%area(m)%obs_area(k) + & 
                              obs%value(t,m,o_k)
                      endif
                      if(model%value(t,m,m_k).ne.0) then 
                         stats%area(m)%model_area(k) = stats%area(m)%model_area(k) + &
                              model%value(t,m,m_k)
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
    
  end subroutine diagnoseSingleAreaMetric


!BOP
! 
! !ROUTINE: LVT_computeAreaMetric
! \label{LVT_computeAreaMetric}
!
! !INTERFACE: 
  subroutine LVT_computeAreaMetric(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute Area values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleAreaMetric](\ref{computeSingleAreaMetric})
!     updates the AreaMetric computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    integer     :: i,m,pass
    logical     :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%area%selectOpt.eq.1.or.LVT_metrics%area%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%area%timeOpt.eq.1.and.&
                  LVT_metrics%area%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%area%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%area%ftn_ts_loc(i,1),200,advance='no') &
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
             call computeSingleAreaMetric(alarm,&
                  obs, model,stats,&
                  LVT_metrics%area)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
          
          if(alarm) then 
             if(LVT_metrics%area%timeOpt.eq.1.and.&
                  LVT_metrics%area%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%area%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%area%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif
             endif
          endif
       endif
    endif
  end subroutine LVT_ComputeAreaMetric

!BOP
! 
! !ROUTINE: computeSingleAreaMetric
! \label{computeSingleAreaMetric}
!
! !INTERFACE: 
  subroutine computeSingleAreaMetric(alarm, obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the Area Metric values for a single variable
!  The arguments are: 
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     Area computation has been reached
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

    implicit none
    logical                 :: alarm
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer  :: t,l,k,m,tind

    if(metric%timeOpt.eq.1.and.alarm) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                stats%area(m)%obs_area(k) = stats%area(m)%obs_area(k)
                stats%area(m)%model_area(k) = stats%area(m)%model_area(k)
             enddo
          enddo
       endif
    endif
  end subroutine computeSingleAreaMetric

!BOP
! 
! !ROUTINE: LVT_writeMetric_AreaMetric
! \label{LVT_writeMetric_AreaMetric}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_AreaMetric(pass,final,vlevels,stats,obs)
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
!   This routine writes the computed area metric values to 
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
  
  end subroutine LVT_writeMetric_AreaMetric

!BOP
! 
! !ROUTINE: LVT_resetMetric_AreaMetric
! \label{LVT_resetMetric_AreaMetric}
!
! !INTERFACE: 
  subroutine LVT_resetMetric_AreaMetric()
! 
! !USES:
    use LVT_coreMod,   only : LVT_rc
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  This routine resets required variables to support the
!  temporal computation of the area metric values. 
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
                if(LVT_metrics%area%timeOpt.eq.1) then 
                   stats%area(m)%model_area(k) = 0.0
                   stats%area(m)%obs_area(k)=0 
                endif
             enddo
          enddo
       endif
       
       model => model%next
       stats = stats%next
    enddo
  end subroutine LVT_resetMetric_AreaMetric

!BOP
! 
! !ROUTINE: LVT_writerestart_AreaMetric
! 
! !INTERFACE:
  subroutine LVT_writerestart_AreaMetric(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Area Metric metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%Area%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for AreaMetric'
       stop
    end if
    
  end subroutine LVT_writerestart_AreaMetric

!BOP
! 
! !ROUTINE: LVT_readrestart_AreaMetric
! 
! !INTERFACE:
  subroutine LVT_readrestart_AreaMetric(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Area Metric metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%Area%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for AreaMetric'
       stop
    end if
    
  end subroutine LVT_readrestart_AreaMetric

  
end module LVT_AreaMetricMod
