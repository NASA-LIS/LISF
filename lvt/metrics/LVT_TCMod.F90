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
! !MODULE: LVT_TCMod
! \label(LVT_TCMod)
!
! !INTERFACE:
module LVT_TCMod
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
!   This module handles the triple collocation computations by comparing
!   the values from the three datastreams. The primary output of the 
!   calculations is the variance of the noise error. 
!   In addition, this module also includes the extended triple
!   collocation (ETC; McColl et al. GRL, 2014) to calculate the 
!   correlation coefficients. 
!
!  McColl, K.A., J. Vogelzang, A.G. Konings, D. Entekhabi, M. Piles, A.
!  Stoffelen (2014). Extended Triple Collocation: Estimating errors and 
!  correlation coefficients with respect to an unknown target. Geophysical 
!  Research Letters 41:6229-6236
!
!  The code uses the multiplicative error model in logarithmic scale based
!  on Alemohammad et al. (2015): Characterization of precipitation product
!  errors across the United States using multiplicative triple collocation, 
!  HESS, 19, 3489-3503.
!
! !NOTES: Currently calculations stratified temporally are not supported. 
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
  public :: LVT_initTC
  public :: LVT_diagnoseTC
  public :: LVT_computeTC 
  public :: LVT_writeMetric_TC
  public :: LVT_resetMetric_TC
  public :: LVT_writerestart_TC
  public :: LVT_readrestart_TC

 type, public :: tcdec
    integer          :: nsize_total
    integer          :: nsize_season
    integer          :: tscale
 end type tcdec

  type(tcdec) :: LVT_tc_struc
contains
  
!BOP
! 
! !ROUTINE: LVT_initTC
! \label{LVT_initTC}
!
! !INTERFACE: 
  subroutine LVT_initTC(selectNlevs,stats,metric)
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
!  the triple collocation computations. 
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

    if(LVT_rc%nDataStreams.ne.3) then 
       write(LVT_logunit,*) '[ERR] Triple collocation requires the use of three'
       write(LVT_logunit,*) '[ERR] data streams. Program stopping ...'
       call LVT_endrun()
    endif

    allocate(stats%tc(LVT_rc%nensem))
    
    do m=1,LVT_rc%nensem

       if(metric%selectOpt.eq.1) then           
          call computeAnalysisLength(LVT_tc_struc%nsize_total, &
               LVT_rc%tavgInterval)

          allocate(stats%tc(m)%value_final(LVT_rc%ngrid, &
               selectNlevs(1),LVT_tc_struc%nsize_total,&
               3))
          allocate(stats%tc(m)%count_value(LVT_rc%ngrid, &
               selectNlevs(1)))

          allocate(stats%tc(m)%errVar(LVT_rc%ngrid, &
               selectNlevs(1),3))

          allocate(stats%tc(m)%rval(LVT_rc%ngrid, &
               selectNlevs(1),3))

          allocate(stats%tc(m)%sx_ds12(LVT_rc%ngrid,&
               selectNlevs(1)))
          allocate(stats%tc(m)%sx_ds13(LVT_rc%ngrid,&
               selectNlevs(1)))
          allocate(stats%tc(m)%sx_ds23(LVT_rc%ngrid,&
               selectNlevs(1)))


          stats%tc(m)%sx_ds12       = 0 
          stats%tc(m)%sx_ds13       = 0 
          stats%tc(m)%sx_ds23       = 0 

          stats%tc(m)%value_final   = LVT_rc%udef 
          stats%tc(m)%count_value   = 0 
          stats%tc(m)%errVar        = LVT_rc%udef
          stats%tc(m)%rval         = LVT_rc%udef

          allocate(stats%tc(m)%value_ci(selectNlevs(1)))
          stats%tc(m)%value_ci = LVT_rc%udef

       endif
    enddo
    
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass       = 2
    metric%nLevs       = 3
    metric%obsData     = .false. 
    metric%stdevFlag   = .false. 

    metric%customNames = .true. 
    metric%nFields     = 2  ! errVar and rho

    allocate(metric%mName(metric%nFields))
    metric%mName(1) = "err_var"
    metric%mName(2) = "rho"

  end subroutine LVT_initTC
  
!BOP
! 
! !ROUTINE: LVT_diagnoseTC
! \label{LVT_diagnoseTC}
!
! !INTERFACE:  
  subroutine LVT_diagnoseTC(pass)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the TC calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleTC](\ref{diagnoseSingleTC})
!     updates the TC computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: pass
    type(LVT_metadataEntry), pointer :: ds1
    type(LVT_metadataEntry), pointer :: ds2
    type(LVT_metadataEntry), pointer :: ds3
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%tc%selectOpt.eq.1.or.&
            LVT_metrics%tc%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(ds1)
          call LVT_getDataStream2Ptr(ds2)
          call LVT_getDataStream3Ptr(ds3)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(ds1))
             call diagnoseSingleTCsum(ds1,ds2,ds3,stats, &
                  LVT_metrics%tc)
             
             ds1 => ds1%next
             ds2 => ds2%next
             ds3 => ds3%next
             stats => stats%next

          enddo
       endif

    elseif(pass.eq.2) then 
       if(LVT_metrics%tc%selectOpt.eq.1.or.&
            LVT_metrics%tc%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(ds1)
          call LVT_getDataStream2Ptr(ds2)
          call LVT_getDataStream3Ptr(ds3)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(ds1))
             call diagnoseSingleTC(ds1,ds2,ds3,stats, &
                  LVT_metrics%tc)
             
             ds1 => ds1%next
             ds2 => ds2%next
             ds3 => ds3%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseTC

!BOP
! 
! !ROUTINE: diagnoseSingleTCsum
! \label{diagnoseSingleTCsum}
!
! !INTERFACE: 
  subroutine diagnoseSingleTCsum(ds1,ds2,ds3, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the TC computation (updates the running 
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

    type(LVT_metaDataEntry) :: ds1
    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds3
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m
    integer    :: tindex_f,m_k,o_k,t_k

    if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1.and.&
         ds3%selectNlevs.ge.1) then 

       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,ds1%selectNlevs

                m_k = k+ds1%startNlevs -1
                o_k = k+ds2%startNlevs -1
                t_k = k+ds3%startNlevs -1

                if(metric%selectOpt.eq.1) then
                   if(ds1%count(t,m,m_k).ne.0.and.&
                        ds2%count(t,m,o_k).ne.0.and.&
                        ds3%count(t,m,t_k).ne.0) then 
                      stats%tc(m)%sx_ds12(t,k) = & 
                           stats%tc(m)%sx_ds12(t,k) + & 
                           ds1%value(t,m,m_k)*ds2%value(t,m,o_k)
                      stats%tc(m)%sx_ds13(t,k) = & 
                           stats%tc(m)%sx_ds13(t,k) + & 
                           ds1%value(t,m,m_k)*ds3%value(t,m,t_k)
                      stats%tc(m)%sx_ds23(t,k) = & 
                           stats%tc(m)%sx_ds23(t,k) + & 
                           ds2%value(t,m,o_k)*ds3%value(t,m,t_k)
                      
                   endif                   
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleTCsum


!BOP
! 
! !ROUTINE: diagnoseSingleTC
! \label{diagnoseSingleTC}
!
! !INTERFACE: 
  subroutine diagnoseSingleTC(ds1,ds2,ds3, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the TC computation (updates the running 
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

    type(LVT_metaDataEntry) :: ds1
    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds3
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m
    integer    :: tindex_f,m_k,o_k,t_k

    if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1.and.&
         ds3%selectNlevs.ge.1) then 

       call getDataTimeIndex(tindex_f, LVT_rc%tavgInterval)

       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,ds1%selectNlevs

                m_k = k+ds1%startNlevs -1
                o_k = k+ds2%startNlevs -1
                t_k = k+ds3%startNlevs -1

                if(metric%selectOpt.eq.1) then
                   if(ds1%count(t,m,m_k).ne.0.and.&
                        ds2%count(t,m,o_k).ne.0.and.&
                        ds3%count(t,m,t_k).ne.0.and.&
                        stats%tc(m)%sx_ds23(t,k).gt.1E-10) then 
                      stats%tc(m)%value_final(t,k,tindex_f,1) = & 
                           ds1%value(t,m,m_k)
                      stats%tc(m)%value_final(t,k,tindex_f,2) = & 
                           ds2%value(t,m,o_k)* & 
                           stats%tc(m)%sx_ds13(t,k)/&
                           stats%tc(m)%sx_ds23(t,k)
                      stats%tc(m)%value_final(t,k,tindex_f,3) = & 
                           ds3%value(t,m,t_k)* & 
                           stats%tc(m)%sx_ds12(t,k)/&
                           stats%tc(m)%sx_ds23(t,k)
                      stats%tc(m)%count_value(t,k) = & 
                           stats%tc(m)%count_value(t,k) + 1
                   endif                   
                endif
             enddo
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleTC


!BOP
! 
! !ROUTINE: LVT_computeTC
! \label{LVT_computeTC}
!
! !INTERFACE: 
  subroutine LVT_computeTC(pass,alarm)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute TC values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleTC](\ref{computeSingleTC})
!     updates the TC computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     TC computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: pass
    logical               :: alarm
    type(LVT_metadataEntry), pointer :: ds1
    type(LVT_metadataEntry), pointer :: ds2
    type(LVT_metadataEntry), pointer :: ds3
    type(LVT_statsEntry)   , pointer :: stats

    integer               :: i,m

    if(pass.eq.2) then 

       if(LVT_metrics%tc%selectOpt.eq.1.or.&
            LVT_metrics%tc%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%tc%timeOpt.eq.1.and.&
                  LVT_metrics%tc%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%tc%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%tc%ftn_ts_loc(i,1),200,advance='no') &
                           LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                           LVT_rc%hr,'',LVT_rc%mn, '' 
                   enddo
                endif


             endif
          endif
200       format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)


          call LVT_getDataStream1Ptr(ds1)
          call LVT_getDataStream2Ptr(ds2)
          call LVT_getDataStream3Ptr(ds3)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(ds1))
             call computeSingleTC(alarm,ds1,ds2,ds3,stats,&
                  LVT_metrics%tc)

             ds1 => ds1%next
             ds2 => ds2%next
             ds3 => ds3%next
             stats => stats%next
          enddo

          if(alarm) then 
             if(LVT_metrics%tc%timeOpt.eq.1.and.&
                  LVT_metrics%tc%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%tc%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%tc%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif

             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeTC

!BOP
! 
! !ROUTINE: computeSingleTC
! \label{computeSingleTC}
!
! !INTERFACE: 
  subroutine computeSingleTC(alarm,ds1,ds2,ds3,stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the TC values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     TC computation has been reached
!    \item[obs] observation object
!    \item[ds1] ds1 variable object
!    \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    logical                 :: alarm
    type(LVT_metaDataEntry) :: ds1
    type(LVT_metaDataEntry) :: ds2
    type(LVT_metaDataEntry) :: ds3
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer                 :: t,l,k,m,i,tind
    integer                 :: trmLength
    real,  allocatable      :: value_final_trm(:,:)

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.ds2%selectNlevs.ge.1) then 

          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,ds1%selectNlevs
                   if(stats%tc(m)%count_value(t,k).ne.0.and.&
                        stats%tc(m)%count_value(t,k) &
                        .gt.LVT_rc%obsCountThreshold) then 
                         
                      call findTrimmedTimeSeriesLength(&
                           LVT_tc_struc%nsize_total,&
                           stats%tc(m)%value_final(t,k,:,:),&
                           trmLength)

                      allocate(value_final_trm(trmLength,3))
                      value_final_trm = LVT_rc%udef
                      call trimTimeSeries(LVT_tc_struc%nsize_total,&
                           stats%tc(m)%value_final(t,k,:,:),&
                           trmLength,&
                           value_final_trm)
                           
                      call computeTCmetrics(&
                           trmLength,&
                           value_final_trm,&
                           stats%tc(m)%errVar(t,k,:),&
                           stats%tc(m)%rval(t,k,:))

                      deallocate(value_final_trm)
                   endif
                
                enddo
             enddo
          enddo
       endif
    endif

  end subroutine computeSingleTC


  subroutine findTrimmedTimeSeriesLength(inputLen,value_final,outputLen)

    integer            :: inputLen, outputLen
    real               :: value_final(inputLen,3)
    
    integer            :: k
    
    outputLen = 0 
    do k=1,inputLen
       if(value_final(k,1).ne.LVT_rc%udef) then 
          outputLen = outputLen + 1          
       endif
    enddo

  end subroutine findTrimmedTimeSeriesLength

  subroutine trimTimeSeries(inputLen,value_final,&
       outputLen, value_final_trm)

    integer            :: inputLen, outputLen
    real               :: value_final(inputLen,3)
    real               :: value_final_trm(outputLen,3)
    
    integer            :: j,k
    
    j = 0 
    do k=1,inputLen
       if(value_final(k,1).ne.LVT_rc%udef) then 
          j = j + 1
          value_final_trm(j,:) = value_final(k,:)
       endif
    enddo

  end subroutine trimTimeSeries
  
!BOP
! 
! !ROUTINE: LVT_writeMetric_TC
! \label{LVT_writeMetric_TC}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_TC(pass,final,vlevels,stats,obs)
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
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs
!EOP

    integer                 :: k,l,m,tind

    real,    allocatable    :: errVar(:,:,:)
    real,    allocatable    :: rho(:,:,:)
    integer, allocatable    :: count_errVar(:,:,:)
    real,    allocatable    :: value_ci(:)

    real,    allocatable    :: value_ts(:,:)
    integer, allocatable    :: count_value_ts(:,:)

    real,    allocatable    :: value_adc(:,:,:)
    real,    allocatable    :: value_asc(:,:,:)

    if(pass.eq.LVT_metrics%tc%npass) then
       if(final.eq.1) then
          if(pass.eq.LVT_metrics%tc%npass) then
             if(stats%selectOpt.eq.1) then

                allocate(errVar(LVT_rc%ngrid,3, LVT_rc%nensem))
                allocate(rho(LVT_rc%ngrid,3, LVT_rc%nensem))
                allocate(count_errVar(LVT_rc%ngrid,3, LVT_rc%nensem))

                allocate(value_ci(LVT_rc%nensem))

                do k=1,vlevels
                   do m=1,LVT_rc%nensem
                      do l=1,3
                         errVar(:,l,m) = &
                              stats%tc(m)%errVar(:,k,l)
                         rho(:,l,m) = &
                              stats%tc(m)%rval(:,k,l)
                         count_errVar(:,l,m) = & 
                              stats%tc(m)%count_value(:,k)
                      enddo
                      value_ci(m) = stats%tc(m)%value_ci(k)

                      if(LVT_metrics%tc%selectOpt.eq.1) then 
                         call LVT_writevar_gridded(LVT_metrics%tc%ftn_total, &
                              3, &
                              errVar,&
                              stats%vid_total(LVT_TCid,1),k)
                         call LVT_writevar_gridded(LVT_metrics%tc%ftn_total, &
                              3, &
                              real(count_errVar),&
                              stats%vid_count_total(LVT_TCid,1),k)

                         call LVT_writevar_gridded(LVT_metrics%tc%ftn_total, &
                              3, &
                              rho,&
                              stats%vid_total(LVT_TCid,3),k)
                         call LVT_writevar_gridded(LVT_metrics%tc%ftn_total, &
                              3, &
                              real(count_errVar),&
                              stats%vid_count_total(LVT_TCid,3),k)

                         call LVT_writeSummaryStats(&
                              LVT_metrics%tc%ftn_summ,&
                              l,&
                              LVT_metrics%tc%short_name,&
                              LVT_rc%ngrid,&
                              errVar(:,1,:), &
                              count_errVar(:,1,:),&
                              stats%standard_name,&
                              value_ci(:))
                      endif
                   enddo
                enddo
                deallocate(errVar)
                deallocate(count_errVar)
                deallocate(rho)
                deallocate(value_ci)

             endif
          endif
       endif
    endif


  end subroutine LVT_writeMetric_TC


!BOP
! 
! !ROUTINE: LVT_resetMetric_TC
! \label(LVT_resetMetric_TC)
!
! !INTERFACE:
  subroutine LVT_resetMetric_TC(alarm)

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

  end subroutine LVT_resetMetric_TC


!BOP
! 
! !ROUTINE: LVT_writerestart_TC
! 
! !INTERFACE:
  subroutine LVT_writerestart_TC(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for TC metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    
  end subroutine LVT_writerestart_TC

!BOP
! 
! !ROUTINE: LVT_readrestart_TC
! 
! !INTERFACE:
  subroutine LVT_readrestart_TC(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for TC metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

  end subroutine LVT_readrestart_TC

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

  subroutine computeTCmetrics (N, A, errVar, rho)
    
    integer                 :: N
    real                    :: A (N,3)
    real                    :: errVar(3)
    real                    :: rho(3)

    real                    :: B (3,3)

    real                    :: Abar (1,3)
    real                    :: Adev (N,3)
    real                    :: AdevHat(3,N)
    integer                 :: i,j,k

    Abar = 0 
    do i=1,N
       do j=1,3
          Abar(1,j) = Abar(1,j) + A(i,j)
       enddo
    enddo

    do j=1,3
       Abar(1,j) = Abar(1,j)/N
    enddo

    do i=1,N
       do j=1,3
          Adev(i,j) = A(i,j)-Abar(1,j)
       enddo
    enddo

!    call transpose(N,3,Adev, AdevHat)
    B = matmul(transpose(Adev), Adev)/N

    if(B(2,3).ne.0) then 
       errVar(1) = (B(1,1) - (B(1,2)*B(1,3)/B(2,3)))
    else
       errVar(1) = LVT_rc%udef
    endif
    if(B(1,3).ne.0) then 
       errVar(2) = (B(2,2) - (B(1,2)*B(2,3)/B(1,3)))
    else
       errVar(2) = LVT_rc%udef
    endif
    if(B(1,2).ne.0) then 
       errVar(3) = (B(3,3) - (B(1,3)*B(2,3)/B(1,2))) 
    else
       errVar(3) = LVT_rc%udef
    endif    
  
    do k=1,3
       if(errVar(k).lt.0) then 
          write(LVT_logunit,*) '[WARN] calculated error variance from TC is negative',errVar(k)
          errVar(k) = LVT_rc%udef
       endif       
    enddo

    if(B(1,1)*B(2,3).ne.0) then 
       rho(1) = (B(1,2)*B(1,3)/(B(1,1)*B(2,3)))
    else
       rho(1) = LVT_rc%udef
    endif
    if(B(2,2)*B(1,3).ne.0) then 
       rho(2) = (B(1,2)*B(2,3)/(B(2,2)*B(1,3)))
    else
       rho(2) = LVT_rc%udef
    endif

    if(B(3,3)*B(1,2).ne.0) then 
       rho(3) = (B(1,3)*B(2,3)/(B(3,3)*B(1,2)))
    else
       rho(3) = LVT_rc%udef
    endif

  end subroutine computeTCmetrics

end module LVT_TCMod
