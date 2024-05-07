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
! !MODULE: LVT_waveletStatMod
! \label(LVT_waveletStatMod)
!
! !INTERFACE:
module LVT_waveletStatMod
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

  public :: LVT_initwaveletStat
  public :: LVT_diagnoseWaveletStat
  public :: LVT_computeWaveletStat
  public :: LVT_writeMetric_waveletStat
  public :: LVT_resetMetric_waveletStat
  public :: LVT_writerestart_waveletStat
  public :: LVT_readrestart_waveletStat

!
! !DESCRIPTION: 
!   This module implements the intensity-scale approach of 
!   Casati et al. (2004) for spatial verification. 
!   The evaluation is done using the mean squared error (MSE) of
!   binary images, obtained by thresholding the model output and
!   the observation using different threshold values. A two-dimensional
!   discrete Haar wavelet decomposition of binary error images is
!   employed to evaluate errors at different spatial scales. 
!
!   Reference: Casati, B., Ross, G. and Stephenson, D.B. "A new 
!   intensity-scale approach for the verification of spatial 
!   precipitation forecasts", Meteorol. Appl. 11, 141--154, 2004. 
!EOP

contains

  subroutine LVT_initwaveletStat(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!EOP
    integer                 :: m

    allocate(stats%wvt(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       
       if(metric%selectOpt.eq.1) then 
          allocate(stats%wvt(m)%value_model_mean_total(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%wvt(m)%value_model_mean_total = 0.0
          allocate(stats%wvt(m)%count_value_model_mean_total(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%wvt(m)%count_value_model_mean_total = 0

          if(LVT_rc%obssource(2).ne."none") then 
             if(selectNlevs(2).ge.1) then 
                allocate(stats%wvt(m)%value_obs_mean_total(LVT_rc%ngrid, selectNlevs(2), &
                     LVT_rc%strat_nlevels))
                stats%wvt(m)%value_obs_mean_total = 0.0
                allocate(stats%wvt(m)%count_value_obs_mean_total(LVT_rc%ngrid, &
                     selectNlevs(2), &
                     LVT_rc%strat_nlevels))
                stats%wvt(m)%count_value_obs_mean_total = 0
             endif
          endif
          allocate(stats%wvt(m)%value_mse_total(LVT_rc%nscales))
          allocate(stats%wvt(m)%value_mse_pct_total(LVT_rc%nscales))
       endif
       if(metric%timeOpt.eq.1) then 
          allocate(stats%wvt(m)%value_model_mean_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%wvt(m)%value_model_mean_ts = 0.0
          allocate(stats%wvt(m)%count_value_model_mean_ts(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%wvt(m)%count_value_model_mean_ts = 0

          if(LVT_rc%obssource(2).ne."none") then 
             if(selectNlevs(2).ge.1) then 
                allocate(stats%wvt(m)%value_obs_mean_ts(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%strat_nlevels))
                stats%wvt(m)%value_obs_mean_ts = 0.0
                allocate(stats%wvt(m)%count_value_obs_mean_ts(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%strat_nlevels))
                stats%wvt(m)%count_value_obs_mean_ts = 0
             endif
          endif
          allocate(stats%wvt(m)%value_mse_ts(LVT_rc%nscales))
          allocate(stats%wvt(m)%value_mse_pct_ts(LVT_rc%nscales))

       endif
    enddo

!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    metric%obsData = .true. 

  end subroutine LVT_initwaveletStat

!BOP
! 
! !ROUTINE: LVT_diagnoseWaveletstat
! \label{LVT_diagnoseWaveletstat}
!
! !INTERFACE: 
  subroutine LVT_diagnoseWaveletstat(pass)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Waveletstat calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleWaveletstat](\ref{diagnoseSingleWaveletstat})
!     updates the Waveletstat computation for a single variable 
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
       if(LVT_metrics%waveletStat%selectOpt.eq.1.or.&
            LVT_metrics%waveletStat%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleWaveletstat(&
                  obs,model,stats,&
                  LVT_metrics%waveletStat)

             model => model%next
             obs   => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseWaveletstat

!BOP
! 
! !ROUTINE: diagnoseSinglewaveletStat
! \label{diagnoseSinglewaveletStat}
!
! !INTERFACE: 
  subroutine diagnoseSinglewaveletStat(obs, model, stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine updates the waveletStat computation (updates the running 
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

    implicit none

    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry) :: stats
    type(LVT_metricEntry)   :: metric

    integer    :: t,k,m,m_k,o_k

    if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(trim(obs%units).eq.trim(model%units)) then 
                   if(metric%selectOpt.eq.1) then
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         stats%wvt(m)%value_model_mean_total(t,k,1) = &
                              stats%wvt(m)%value_model_mean_total(t,k,1) + &
                              model%value(t,m,m_k)
                         stats%wvt(m)%count_value_model_mean_total(t,k,1) = &
                              stats%wvt(m)%count_value_model_mean_total(t,k,1) + 1
                      endif
                      if(obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                         stats%wvt(m)%value_obs_mean_total(t,k,1) = &
                              stats%wvt(m)%value_obs_mean_total(t,k,1) + &
                              obs%value(t,m,o_k)
                         stats%wvt(m)%count_value_obs_mean_total(t,k,1) = &
                              stats%wvt(m)%count_value_obs_mean_total(t,k,1) + 1
                      endif
                      if(metric%timeOpt.eq.1) then 
                         if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                            stats%wvt(m)%value_model_mean_ts(t,k,1) = &
                                 stats%wvt(m)%value_model_mean_ts(t,k,1)+&
                                 model%value(t,m,m_k)
                            stats%wvt(m)%count_value_model_mean_ts(t,k,1) = & 
                                 stats%wvt(m)%count_value_model_mean_ts(t,k,1)+1
                         endif
                         if(obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            stats%wvt(m)%value_obs_mean_ts(t,k,1) = &
                                 stats%wvt(m)%value_obs_mean_ts(t,k,1)+&
                                 obs%value(t,m,o_k)
                            stats%wvt(m)%count_value_obs_mean_ts(t,k,1) = & 
                                 stats%wvt(m)%count_value_obs_mean_ts(t,k,1)+1
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
  end subroutine diagnoseSinglewaveletStat

!BOP
! 
! !ROUTINE: LVT_computeWaveletStat
! \label{LVT_computeWaveletStat}
!
! !INTERFACE: 
  subroutine LVT_computeWaveletStat(pass,alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute WaveletStat values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleWaveletStat](\ref{computeSingleWaveletStat})
!     updates the WaveletStat computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     WaveletStat computation has been reached
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none

    integer               :: pass
    logical     :: alarm
    integer     :: i,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%waveletStat%selectOpt.eq.1.or.&
            LVT_metrics%waveletStat%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%waveletStat%timeOpt.eq.1.and.&
                  LVT_metrics%waveletStat%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%bias%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%bias%ftn_ts_loc(i,1),200,advance='no') &
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
             call computeSingleWaveletStat(alarm,&
                  obs, model,stats,&
                  LVT_metrics%waveletStat)

             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%waveletStat%timeOpt.eq.1.and.&
                  LVT_metrics%waveletStat%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%bias%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%bias%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif


             endif
          endif
       endif
    endif

  end subroutine LVT_ComputeWaveletStat

!BOP
! 
! !ROUTINE: computeSingleWaveletStat
! \label{computeSingleWaveletStat}
!
! !INTERFACE: 
  subroutine computeSingleWaveletStat(alarm,obs, model,stats,metric)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the WaveletStat values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[check]
!     boolean flag indicating if the specified interval for 
!     WaveletStat computation has been reached
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

    integer  :: t,i,l,k,m,tind
    real     :: Ix(LVT_rc%ngrid), Iy(LVT_rc%ngrid),Zval(LVT_rc%ngrid)


    if(metric%timeOpt.eq.1.and.alarm) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%wvt(m)%count_value_model_mean_ts(t,k,l).ne.0) then 
                         stats%wvt(m)%value_model_mean_ts(t,k,l) = &
                              stats%wvt(m)%value_model_mean_ts(t,k,l)/&
                              stats%wvt(m)%count_value_model_mean_ts(t,k,l) 
                      else
                         stats%wvt(m)%value_model_mean_ts(t,k,l) = LVT_rc%udef
                      endif
                      
                      if(stats%wvt(m)%count_value_obs_mean_ts(t,k,l).ne.0) then 
                         stats%wvt(m)%value_obs_mean_ts(t,k,l) = &
                              stats%wvt(m)%value_obs_mean_ts(t,k,l)/&
                              stats%wvt(m)%count_value_obs_mean_ts(t,k,l) 
                      else
                         stats%wvt(m)%value_obs_mean_ts(t,k,l) = LVT_rc%udef
                      endif
                      
                   enddo
                enddo
             enddo
          enddo
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   do t=1,LVT_rc%ngrid
                      if(stats%wvt(m)%value_model_mean_ts(t,k,l).gt.metric%minthreshold) then 
                         Ix(t) = 1
                      else
                         Ix(t) = 0
                      endif
                      if(stats%wvt(m)%value_obs_mean_ts(t,k,l).gt.metric%minthreshold) then 
                         Iy(t) = 1
                      else
                         Iy(t) = 0
                      endif
                      
                      Zval(t) = Iy(t)-Ix(t)
                   enddo
                   call waveletScaleDecomp(Zval, stats%wvt(m)%value_mse_ts, &
                        stats%wvt(m)%value_mse_pct_ts)
                enddo
             enddo
          enddo
          if(metric%extractTS.eq.1) then 
             do m=1,LVT_rc%nensem
                do i=1, LVT_rc%ntslocs
                   do k=1,LVT_rc%nscales
                      write(metric%ftn_ts_loc(i,m),203,advance='no') &
                           stats%wvt(m)%value_mse_ts(k), stats%wvt(m)%value_mse_pct_ts(k)
                   enddo
                enddo
             end do
          endif
       endif
    endif
203 format(2E14.6)

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.obs%selectNlevs.ge.1) then
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%wvt(m)%count_value_model_mean_total(t,k,l).ne.0) then 
                         stats%wvt(m)%value_model_mean_total(t,k,l) = &
                              stats%wvt(m)%value_model_mean_total(t,k,l)/&
                              stats%wvt(m)%count_value_model_mean_total(t,k,l) 
                      else
                         stats%wvt(m)%value_model_mean_total(t,k,l) = LVT_rc%udef
                      endif
                      
                      if(stats%wvt(m)%count_value_obs_mean_total(t,k,l).ne.0) then 
                         stats%wvt(m)%value_obs_mean_total(t,k,l) = &
                              stats%wvt(m)%value_obs_mean_total(t,k,l)/&
                              stats%wvt(m)%count_value_obs_mean_total(t,k,l) 
                      else
                         stats%wvt(m)%value_obs_mean_total(t,k,l) = LVT_rc%udef
                      endif
                      
                   enddo
                enddo
             enddo
          enddo
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels
                   do t=1,LVT_rc%ngrid
                      if(stats%wvt(m)%value_model_mean_total(t,k,l).gt.metric%minthreshold) then 
                         Ix(t) = 1
                      else
                         Ix(t) = 0
                      endif
                      
                      if(stats%wvt(m)%value_obs_mean_total(t,k,l).gt.metric%minthreshold) then 
                         Iy(t) = 1
                      else
                         Iy(t) = 0
                      endif
                      Zval(t) = Iy(t)-Ix(t)
                   enddo
                   
                   call waveletScaleDecomp(Zval, stats%wvt(m)%value_mse_total,&
                        stats%wvt(m)%value_mse_pct_total)
                enddo
             enddo
          enddo
       endif
    endif
  end subroutine computeSingleWaveletStat


!BOP
! 
! !ROUTINE: LVT_writeMetric_waveletstat
! \label(LVT_writeMetric_waveletstat)
!
! !INTERFACE:
  subroutine LVT_writeMetric_waveletstat(pass,final,vlevels,stats,obs)
! 
! !USES:   
    use LVT_logMod, only : LVT_getNextUnitNumber, LVT_releaseUnitNumber
    use LVT_statsMod, only : LVT_writeSummaryStats
    use LVT_pluginIndices
    implicit none
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

    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs

    integer                 :: i ,k,m
    integer                 :: ftn
    character*100           :: filename
    

    if(LVT_metrics%waveletstat%selectOpt.eq.1) then
       filename = trim(LVT_rc%statsodir)//'/LVT_Waveletstat_'//&
            trim(stats%standard_name)//'.dat'
       
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=trim(filename),form='formatted')
       

       if(pass.eq.LVT_metrics%waveletstat%npass) then 
          if(final.eq.1) then 
             if(pass.eq.1) then 
                if(stats%selectOpt.eq.1) then 
                   do m=1,LVT_rc%nensem
                      do i=1, LVT_rc%ntslocs
                         do k=1,LVT_rc%nscales
                            write(ftn,204) &
                                 k, stats%wvt(m)%value_mse_total(k), &
                                 stats%wvt(m)%value_mse_pct_total(k)
                         enddo
                      enddo
                   enddo
                endif
             endif
          endif
       end if
       call LVT_releaseUnitNumber(ftn)

    endif

204 format(I2.2, E14.6, E14.2)

  end subroutine LVT_writeMetric_waveletstat

!BOP
! 
! !ROUTINE: LVT_resetMetric_waveletstat
! \label(LVT_resetMetric_waveletstat)
!
! !INTERFACE:
  subroutine LVT_resetMetric_waveletstat
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

    integer                :: i,k,l,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats


    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream1Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                if(LVT_metrics%waveletstat%timeOpt.eq.1) then 
                   do l=1,LVT_rc%strat_nlevels
                      stats%wvt(m)%value_model_mean_ts(:,k,l) = 0.0
                      stats%wvt(m)%count_value_model_mean_ts(:,k,l)=0 
                   enddo
                   if(LVT_rc%obssource(2).ne."none".and.&
                        obs%selectNlevs.ge.1) then 
                      do l=1,LVT_rc%strat_nlevels
                         stats%wvt(m)%value_obs_mean_ts(:,k,l) = 0.0
                         stats%wvt(m)%count_value_obs_mean_ts(:,k,l)=0 
                      enddo
                   endif
                endif
             enddo
          enddo
       endif
       
       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo

  end subroutine LVT_resetMetric_waveletstat


!BOP
! 
! !ROUTINE: LVT_writerestart_waveletStat
! 
! !INTERFACE:
  subroutine LVT_writerestart_waveletStat(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for waveletStat metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%waveletstat%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for waveletStat'
       stop
    end if
    
  end subroutine LVT_writerestart_waveletStat

!BOP
! 
! !ROUTINE: LVT_readrestart_waveletStat
! 
! !INTERFACE:
  subroutine LVT_readrestart_waveletStat(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for waveletStat metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%waveletstat%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for waveletStat'
       stop
    end if
    
  end subroutine LVT_readrestart_waveletStat


  subroutine waveletScaleDecomp(Zval,MSE, MSE_PCT)
    
    use LVT_coreMod, only : LVT_rc, LVT_domain

    real, intent(in) :: Zval(LVT_rc%ngrid)
    real             :: MSE(LVT_rc%nscales)
    real             :: MSE_PCT(LVT_rc%nscales)

    real             :: Zval_2d(LVT_rc%lnc, LVT_rc%lnr)
    real             :: Zval_pad_2d(LVT_rc%lnc_sc(1),&
         LVT_rc%lnr_sc(1))
    real             :: fthr_wvt(LVT_rc%lnc_sc(1),&
         LVT_rc%lnr_sc(1),LVT_rc%nscales)
    real             :: mthr_wvt(LVT_rc%lnc_sc(1),&
         LVT_rc%lnr_sc(1),LVT_rc%nscales)
    integer          :: c,r,t,m
    integer          :: offset1_1,offset1_2,offset2_1,offset2_2
    integer          :: gid
    real             :: MSE0
    real             :: thresh


!convert Zval to 2d 
    zval_2d = 0
    fthr_wvt = 0 
    mthr_wvt = 0 
    do t=1,LVT_rc%ngrid

       c = LVT_domain%grid(t)%col
       r = LVT_domain%grid(t)%row

       gid = LVT_domain%gindex(c,r)

       if(gid.ne.-1) then 
          Zval_2d(c,r) = Zval(gid)          
       endif

    end do

!    open(100,file='var_orig.bin',form='unformatted')
!    write(100) zval_2d
!    close(100)

!pad to the highest scale.., the original data should be in the
!center of the domain.---redo
    offset1_1 = nint((LVT_rc%lnr_sc(1)-LVT_rc%lnr)/2.0)
    offset1_2 = LVT_rc%lnr_sc(1)-offset1_1-LVT_rc%gnr
    offset2_1 = nint((LVT_rc%lnc_sc(1)-LVT_rc%lnc)/2.0)
    offset2_2 = LVT_rc%lnc_sc(1)-offset2_1-LVT_rc%gnc
    
!    print*, LVT_rc%lnc_sc(1), LVT_rc%lnr_sc(1)
    
    Zval_pad_2d = 0.0
    Zval_pad_2d(offset2_1:offset2_1+LVT_rc%gnc-1,&
         offset1_1:offset1_1+LVT_rc%gnr-1) = Zval_2d(:,:)

!    open(100,file='var_raw.bin',form='unformatted')
!    write(100) Zval_pad_2d
!    close(100)
!    stop

    call wvt_decomp1(zval_pad_2d,LVT_rc%lnc_sc(1),LVT_rc%lnr_sc(1), &
         fthr_wvt(:,:,1), mthr_wvt(:,:,1))

    MSE0 = 0 
    do r=1,LVT_rc%lnr_sc(1)
       do c=1,LVT_rc%lnc_sc(1)
          MSE0 = MSE0 + zval_pad_2d(c,r)**2
       enddo
    enddo
    MSE0 = MSE0/float(LVT_rc%lnc_sc(1)*LVT_rc%lnr_sc(1))

    if(MSE0.ne.0) then 
       MSE = 0 
       MSE_PCT = 0 
       MSE(1) = 0 
       do r=1,LVT_rc%lnr_sc(1)
          do c=1,LVT_rc%lnc_sc(1)
             MSE(1) = MSE(1) + mthr_wvt(c,r,1)**2
          enddo
       enddo
       MSE(1) = MSE(1)/float(LVT_rc%lnc_sc(1)*LVT_rc%lnr_sc(1))
       mse_pct(1) = MSE(1)*100.0/MSE0
       
       do m=2,LVT_rc%nscales
          call wvt_decomp(zval_pad_2d,LVT_rc%lnc_sc(1),LVT_rc%lnr_sc(1), &
               m, fthr_wvt,mthr_wvt)
          
          do r=1,LVT_rc%lnr_sc(1)
             do c=1,LVT_rc%lnc_sc(1)
                MSE(m) = MSE(m) + mthr_wvt(c,r,m)**2
             enddo
          enddo
          MSE(m) = MSE(m)/float(LVT_rc%lnc_sc(1)*LVT_rc%lnr_sc(1))
          mse_pct(m) = MSE(m)*100.0/MSE0
       enddo
    else
       MSE = LVT_rc%udef
       MSE_PCT = LVT_rc%udef
    endif
  end subroutine WaveletScaleDecomp

  subroutine wvt_decomp1(data_in, nc,nr,fthr_wvt,mthr_wvt)

    integer    :: nc,nr
    real       :: data_in(nc,nr)

    real       :: fthr_wvt(nc,nr)
    real       :: mthr_wvt(nc,nr)

    integer       :: count
    real          :: tempval
    real, allocatable :: x(:,:)

    integer    :: c,r,c1,r1,dim,dx,i,j
    
    dx = 2
    dim = 2**(LVT_rc%nscales)/2
    allocate(x(dim,dim))

    c1 = 1
    r1 = 1
    do i=1,dim
       do j=1,dim
          tempval = 0 
          count = 0 
          do c=c1,c1+dx-1
             do r=r1,r1+dx-1
                tempval = tempval + data_in(c,r)
                count = count + 1
             enddo
          enddo
          x(i,j) = tempval/count
          c1 = c1+dx
       enddo
       c1 = 1
       r1 = r1+dx
    enddo
    
    c1 = 1
    r1 = 1

    do i=1,dim
       do j=1,dim
          fthr_wvt(c1:c1+dx-1,r1:r1+dx-1) = x(i,j)
          c1 = c1+dx
       enddo
       c1 = 1
       r1 = r1+dx
    enddo

    do r=1,nr
       do c=1,nc
          mthr_wvt(c,r) = data_in(c,r) - fthr_wvt(c,r)
       enddo
    enddo

    deallocate(x)    
    
    
  end subroutine wvt_decomp1

  subroutine wvt_decomp(data_in, nc,nr,scale, fthr_wvt,mthr_wvt)

    integer    :: nc,nr
    real       :: data_in(nc,nr)
    integer    :: scale

    real       :: fthr_wvt(nc,nr,LVT_rc%nscales)    
    real       :: mthr_wvt(nc,nr,LVT_rc%nscales)

    integer       :: count
    real          :: tempval
    real, allocatable :: x(:,:)

    integer    :: c,r,c1,r1,dim,dx,i,j
    
    dx = 2**(scale)
    dim = 2**(LVT_rc%nscales)/2**(scale)
    allocate(x(dim,dim))

    c1 = 1
    r1 = 1
    do i=1,dim
       do j=1,dim
          tempval = 0 
          count = 0 
          do c=c1,c1+dx-1
             do r=r1,r1+dx-1
                tempval = tempval + data_in(c,r)
                count = count + 1
             enddo
          enddo
          x(i,j) = tempval/count
          c1 = c1+dx
       enddo
       c1 = 1
       r1 = r1+dx
    enddo
    
    c1 = 1
    r1 = 1

    do i=1,dim
       do j=1,dim
          fthr_wvt(c1:c1+dx-1,r1:r1+dx-1,scale) = x(i,j)
          c1 = c1+dx
       enddo
       c1 = 1
       r1 = r1+dx
    enddo

    do r=1,nr
       do c=1,nc
          mthr_wvt(c,r,scale) = fthr_wvt(c,r,scale)-fthr_wvt(c,r,scale-1)
       enddo
    enddo

    deallocate(x)    

  end subroutine wvt_decomp

#if 0 
  subroutine fwt_2d(Zval, nc,nr)
    
    implicit none
    integer :: nc, nr
    real    :: Zval(nc,nr)
    
    integer :: c,r
    
    do r=1,nr
       call fwt_1d(Zval(:,r), nc, 1,nc/2) !forward
    enddo
    do c=1,nc
       call fwt_1d(Zval(c,:),nr,1,nr/2)
    enddo
  end subroutine fwt_2d
  
  subroutine fwti_2d(Zval, nc,nr,scale)
    
    implicit none
    
    integer :: nc, nr
    real    :: Zval(nc,nr)
    integer :: scale

    integer :: c,r
    
    do r=1,nr
       call fwt_1d(Zval(:,r), nc, -1,scale) !inverse
    enddo
    do c=1,nc
       call fwt_1d(Zval(c,:),nr,-1,scale)
    enddo
    
  end subroutine fwti_2d

  subroutine fwt_1d(data_1d, N, direction, scale)
    
    implicit none
    
    integer :: N
    real    :: data_1d(N)
    integer :: direction
    integer :: scale

    integer :: thisn

    
    if(N>2) then 
       if(direction.eq.1) then  !forward
          thisn = N
          do while(thisn>=2) 
             call wavelet(data_1d,thisn,direction)
             thisn = thisn/2
          enddo
       elseif(direction.eq.-1) then !inverse
          thisn = 2
          do while(thisn .le.scale) 
             call wavelet(data_1d, thisn, direction)
             thisn = thisn*2
          enddo
       endif
    endif
  end subroutine fwt_1d
  
  subroutine wavelet(data_1d, N, direction)
    
    implicit none
    
    integer :: N 
    integer :: direction
    real    :: data_1d(N)
    
    call HaarWT(data_1d, N, direction)
    
  end subroutine wavelet
  
  subroutine HaarWT(data_1d, N, direction)
  
    implicit none
    
    integer :: N 
    integer :: direction
    real    :: data_1d(N)
    
    real    :: h0
    real    :: h1
    integer :: nover2
    integer :: i,j
    real    :: tmp(N)
    
    h0 = 0.5 
    h1 = 0.5 
    nover2 = N/2
    
    if(direction.eq.1) then 
       i = 1
       do j=1,N,2
          tmp(i)        = h0*data_1d(j) + h1*data_1d(j+1)
          tmp(i+nover2) = h0*data_1d(j) - h1*data_1d(j+1)
          i = i+1
       enddo
    elseif(direction.eq.-1) then 
       h1 = 1.0
       h0 = 1.0 
       i = 1
       j = 1
       do i=1,nover2
          tmp(j)   = h0*data_1d(i) + h1*data_1d(i+nover2)
          tmp(j+1) = h0*data_1d(i) - h1*data_1d(i+nover2)
          j = j +2
       enddo
    endif
    do i=1,N
       data_1d(i) = tmp(i)
    enddo
  end subroutine HaarWT
#endif
  
end module LVT_waveletStatMod
