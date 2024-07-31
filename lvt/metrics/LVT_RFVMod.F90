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
! !MODULE: LVT_RFVMod
! \label(LVT_RFVMod)
!
! !INTERFACE:
module LVT_RFVMod
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
! !DESCRIPTION: 
!  This module handles the computations required to compute the
!  river flow variate (RFV) metric. This metric computes the day of the year
!  on which 50% of the total water flow for the year has been reached. 
!  The metric is computed based on the streamflow variable. 
!
!  Ref: Barnett et al., Science (2008)
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  25 Jun 2013   Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initRFV
  public :: LVT_diagnoseRFV
  public :: LVT_computeRFV
  public :: LVT_writeMetric_RFV
  public :: LVT_resetMetric_RFV
  public :: LVT_writerestart_RFV
  public :: LVT_readrestart_RFV
!EOP

  type, public :: rfvdec
     real,  allocatable :: model_cdf(:,:)
     real,  allocatable :: obs_cdf(:,:)
     real,  allocatable :: model_rfv_val(:,:,:)
     real,  allocatable :: obs_rfv_val(:,:,:)
     integer        :: prev_yr_tavg
     logical        :: change_of_year
     real           :: model_rfv_ci(1,1)
     real           :: obs_rfv_ci(1,1)
  end type rfvdec
  type(rfvdec), allocatable  :: LVT_rfv_struc(:)

  private

contains
  subroutine LVT_initRFV(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!EOP
    integer                 :: m
    
    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          if(LVT_rc%tavgInterval.ne.86400) then 
             write(LVT_logunit,*) 'The time averaging interval should be a day'
             write(LVT_logunit,*) 'when using the river flow variate metric '
             call LVT_endrun()
          endif
          if(LVT_rc%statswriteint.ne.31536000) then 
             write(LVT_logunit,*) 'The stats output intervalinterval should be a year'
             write(LVT_logunit,*) 'when using the river flow variate metric '
             call LVT_endrun()
          endif

!          if(model%short_name.eq."Streamflow") then 
             allocate(LVT_rfv_struc(m)%model_cdf(LVT_rc%ngrid,366))
             allocate(LVT_rfv_struc(m)%obs_cdf(LVT_rc%ngrid,366))

             LVT_rfv_struc(m)%model_cdf = 0 
             LVT_rfv_struc(m)%obs_cdf = 0 

             allocate(LVT_rfv_struc(m)%model_rfv_val(LVT_rc%ngrid,1,1))
             allocate(LVT_rfv_struc(m)%obs_rfv_val(LVT_rc%ngrid,1,1))

             LVT_rfv_struc(m)%model_rfv_val = LVT_rc%udef
             LVT_rfv_struc(m)%obs_rfv_val = LVT_rc%udef
!          endif
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    LVT_rfv_struc(:)%prev_yr_tavg = LVT_rc%syr
    LVT_rfv_struc(:)%change_of_year = .false. 

    if(LVT_rc%obssource(2).ne."none") then 
       metric%obsData = .true. 
    else
       metric%obsData = .false. 
    endif
    metric%stdevFlag = .false. 

  end subroutine LVT_initRFV

!BOP
! 
! !ROUTINE: LVT_diagnoseRFV
! \label{LVT_diagnoseRFV}
!
! !INTERFACE: 
  subroutine LVT_diagnoseRFV(pass)
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the computations for 
!   calculating the rfv of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelRFV](\ref{diagnoseSingleModelRFV})
!     updates the rfv computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer       :: pass
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%rfv%selectOpt.eq.1.or.&
            LVT_metrics%rfv%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleRFV(model,obs,stats,&
                  LVT_metrics%rfv)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseRFV

!BOP
! 
! !ROUTINE: diagnoseSingleRFV
! \label{diagnoseSingleRFV}
!
! !INTERFACE: 
  subroutine diagnoseSingleRFV(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the rfv computation of the 
!   specified variable. 
!
!  The arguments are: 
!
!  \begin{description}
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
    type(LVT_metaDataEntry) :: model
    type(LVT_metadataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP
    integer    :: t,k,m,tind,m_k,o_k

    if(stats%selectOpt.eq.1.and.&
         model%short_name.eq."Streamflow".and.&
         model%selectNlevs.ge.1) then        
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(model%count(t,m,m_k).gt.0) then
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         if(LVT_rc%doy.gt.1) then 
                            LVT_rfv_struc(m)%model_cdf(t,LVT_rc%doy) = &
                                 LVT_rfv_struc(m)%model_cdf(t,LVT_rc%doy-1) + & 
                                 model%value(t,m,m_k)
                         endif
                      endif
                   endif
                endif
                
                if(LVT_rc%obssource(2).ne."none") then
                   if(obs%selectNlevs.ge.1.and.&
                        model%short_name.eq."Streamflow") then
                      if(obs%count(t,m,o_k).gt.0) then 
                         if(metric%selectOpt.eq.1.and.&
                              obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            if(LVT_rc%doy.gt.1) then 
                               LVT_rfv_struc(m)%obs_cdf(t,LVT_rc%doy) = & 
                                    LVT_rfv_struc(m)%obs_cdf(t,LVT_rc%doy-1) + & 
                                    obs%value(t,m,o_k)
                            endif
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
  end subroutine diagnoseSingleRFV

!BOP
! 
! !ROUTINE: LVT_computeRFV
! \label{LVT_computeRFV}
!
! !INTERFACE: 
  subroutine LVT_computeRFV(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the rfv values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelRFV](\ref{computeSingleModelRFV})
!     computes the rfv values for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer               :: pass
    logical               :: alarm
!EOP
    integer     :: i,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    
    if(pass.eq.1) then 
       if(LVT_metrics%rfv%selectOpt.eq.1.or.&
            LVT_metrics%rfv%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%rfv%timeOpt.eq.1.and.&
                  LVT_metrics%rfv%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%rfv%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%rfv%ftn_ts_loc(i,1),200,advance='no') &
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
             
             call computeSingleRFV(alarm,model,obs,stats,&
                  LVT_metrics%rfv)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%rfv%timeOpt.eq.1.and.&
                  LVT_metrics%rfv%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%rfv%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%rfv%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif

             endif
          endif
       endif
    endif
  end subroutine LVT_computeRFV
  

!BOP
! 
! !ROUTINE: computeSingleRFV
! \label{computeSingleRFV}
!
! !INTERFACE: 
  subroutine computeSingleRFV(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the rfv values
!  The arguments are: 
!
!  \begin{description}
!    \item[model] model variable object
!    \item[stats] object to hold the updated statistics
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
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP

    integer  :: t,l,k,m
    real     :: maxv, thresh_val
    integer  :: count_model_rfv(LVT_rc%ngrid,1,1)
    integer  :: count_obs_rfv(LVT_rc%ngrid,1,1)

    LVT_rfv_struc(:)%change_of_year = .false. 
       
    count_model_rfv = 1
    count_obs_rfv   = 1
    
    do m=1,LVT_rc%nensem
       if(LVT_rc%timeAvgOpt.eq.0) then 
          if(LVT_rc%nyr.ne.LVT_rfv_struc(m)%prev_yr_tavg) then 
             LVT_rfv_struc(m)%prev_yr_tavg = LVT_rc%nyr
             LVT_rfv_struc(m)%change_of_year = .true. 
          endif
       else
          if(LVT_rc%yr.ne.LVT_rfv_struc(m)%prev_yr_tavg) then 
             LVT_rfv_struc(m)%prev_yr_tavg = LVT_rc%yr
             LVT_rfv_struc(m)%change_of_year = .true. 
          endif
       endif
    enddo

    if(LVT_rfv_struc(1)%change_of_year) then 
       if(metric%timeOpt.eq.1.and.alarm) then 
          if(stats%selectOpt.eq.1.and.&
               model%short_name.eq."Streamflow".and.& 
               model%selectNlevs.ge.1) then 
             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem
                   maxv = -1E10
                   do l=1,366
                      if(LVT_rfv_struc(m)%model_cdf(t,l).gt.maxv) then 
                         maxv = LVT_rfv_struc(m)%model_cdf(t,l)
                      endif
                   enddo
                   thresh_val = maxv*0.50
                   LVT_rfv_struc(m)%model_rfv_val(t,1,1) = LVT_rc%udef
                   if(thresh_val.gt.0) then 
                      do l=1,366
                         if(LVT_rfv_struc(m)%model_cdf(t,l).gt.thresh_val) then 
                            LVT_rfv_struc(m)%model_rfv_val(t,1,1) = l
                            exit
                         endif
                      enddo
                   endif
                enddo
             enddo
             
             if(LVT_rc%obssource(2).ne."none") then 
                do t=1,LVT_rc%ngrid
                   do m=1,LVT_rc%nensem
                      maxv = -1E10
                      do l=1,366
                         if(LVT_rfv_struc(m)%obs_cdf(t,l).gt.maxv) then 
                            maxv = LVT_rfv_struc(m)%obs_cdf(t,l)
                         endif
                      enddo
                      
                      thresh_val = maxv*0.50
                      LVT_rfv_struc(m)%obs_rfv_val(t,1,1) = LVT_rc%udef
                      if(thresh_val.gt.0) then 
                         do l=1,366
                            if(LVT_rfv_struc(m)%obs_cdf(t,l).gt.thresh_val) then 
                               LVT_rfv_struc(m)%obs_rfv_val(t,1,1) = l
                               exit
                            endif
                         enddo
                      endif
                   enddo
                enddo
             endif
                
                
             if(metric%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                         call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                              model,&
                              LVT_rc%ngrid,&
                              LVT_rfv_struc(m)%model_rfv_val,&
                              count_model_rfv,&
                              LVT_rc%ngrid,&
                              LVT_rfv_struc(m)%obs_rfv_val,&
                              count_obs_rfv)
                      else
                         call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                              model,&
                              LVT_rc%ngrid,&
                              LVT_rfv_struc(m)%model_rfv_val,&
                              count_model_rfv)
                      endif
                   enddo
                   
                else
!TBD
                endif
             endif
          endif
          do m=1,LVT_rc%nensem
             call LVT_computeCI(LVT_rfv_struc(m)%model_rfv_val(:,1,1),LVT_rc%ngrid,&
                  LVT_rc%pval_CI,LVT_rfv_struc(m)%model_rfv_ci(1,1))
             
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                call LVT_computeCI(LVT_rfv_struc(m)%obs_rfv_val(:,1,1),LVT_rc%ngrid,&
                     LVT_rc%pval_CI,LVT_rfv_struc(m)%obs_rfv_ci(1,1))
             endif
          enddo

!          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
!               LVT_rc%ngrid,LVT_rfv_struc(m)%model_rfv_val)
!          
!          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
!             call LVT_writeDataBasedStrat(model,obs,stats,metric,&
!                  LVT_rc%ngrid, LVT_rfv_struc(m)%model_rfv_val, LVT_rc%ngrid,&
!                  LVT_rfv_struc(m)%obs_rfv_val)
!          endif
       endif
    endif
  end subroutine computeSingleRFV

!BOP
! 
! !ROUTINE: LVT_writeMetric_RFV
! \label(LVT_writeMetric_RFV)
!
! !INTERFACE:
  subroutine LVT_writeMetric_RFV(pass,final,vlevels,stats,obs)
! 
! !USES:   
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

    integer                 :: k,l,m,tind

    if(pass.eq.LVT_metrics%rfv%npass) then
       if(stats%selectOpt.eq.1) then 
          print*, 'Writemetric needs to be implemented for RFV'
          stop
       endif
    endif
#if 0  

    if(LVT_rfv_struc(m)%change_of_year) then 
       if(pass.eq.LVT_metrics%rfv%npass) then 
          if(final.ne.1) then
             if(stats%selectOpt.eq.1) then 
                if(LVT_metrics%rfv%timeOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%rfv%ftn_ts, &
                        LVT_rfv_struc(m)%model_rfv_val(:,1,1),&
                        stats%vid_ts(LVT_RFVid,1),1)
                   if(LVT_rc%obssource(2).ne."none") then 
                      call LVT_writevar_gridded(LVT_metrics%rfv%ftn_ts, &
                           LVT_rfv_struc(m)%obs_rfv_val(:,1,1),&
                           stats%vid_ts(LVT_RFVid,2),1)
                   endif
                endif
             endif
          endif
       endif
       LVT_rfv_struc(m)%model_rfv_val = LVT_rc%udef
       LVT_rfv_struc(m)%obs_rfv_val = LVT_rc%udef
    endif
#endif
   end subroutine LVT_writeMetric_RFV

!BOP
! 
! !ROUTINE: LVT_resetMetric_RFV
! \label(LVT_resetMetric_RFV)
!
! !INTERFACE:
  subroutine LVT_resetMetric_RFV
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


    LVT_rfv_struc(m)%prev_yr_tavg = LVT_rc%syr

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream1Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(stats%selectOpt.eq.1) then 
          do m=1,LVT_rc%nensem
             if(LVT_metrics%rfv%timeOpt.eq.1) then 
                LVT_rfv_struc(m)%model_rfv_val(:,1,1) = 0.0
                if(LVT_rc%obssource(2).ne."none".and.&
                     obs%selectNlevs.ge.1) then 
                   LVT_rfv_struc(m)%obs_rfv_val(:,1,1) = 0.0
                endif
             endif
          enddo
       endif
       
       model => model%next
       obs => obs%next
       stats => stats%next
    enddo

  end subroutine LVT_resetMetric_RFV

!BOP
! 
! !ROUTINE: LVT_writerestart_RFV
! 
! !INTERFACE:
  subroutine LVT_writerestart_RFV(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for RFV metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,l
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

#if 0 
    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%rfv%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_writevar_restart(ftn,&
                        stats%model_rfv_total(:,k,l))
                   call LVT_writevar_restart(ftn,&
                        stats%count_model_rfv_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%rfv%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_writevar_restart(ftn,&
                           stats%model_rfv_asc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%count_model_rfv_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%rfv%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_writevar_restart(ftn,&
                           stats%model_rfv_adc(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%count_model_rfv_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%vlevels
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_writevar_restart(ftn,&
                              stats%obs_rfv_total(:,k,l))
                         call LVT_writevar_restart(ftn,&
                              stats%count_obs_rfv_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%rfv%computeSC.eq.1) then 
                      do k=1,obs%vlevels
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%obs_rfv_asc(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%count_obs_rfv_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%rfv%computeADC.eq.1) then 
                      do k=1,obs%vlevels
                         do l=1,LVT_rc%nasc
                            call LVT_writevar_restart(ftn,&
                                 stats%obs_rfv_adc(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%count_obs_rfv_adc(:,k,l))
                         enddo
                      enddo
                   endif
                endif
             endif
          endif
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
    end do
#endif
  end subroutine LVT_writerestart_RFV


!BOP
! 
! !ROUTINE: LVT_readrestart_RFV
! 
! !INTERFACE:
  subroutine LVT_readrestart_RFV(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for RFV metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    integer              :: k,l,index
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

#if 0 
    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%rfv%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 

             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels         
                   call LVT_readvar_restart(ftn,&
                        stats%model_rfv_total(:,k,l))
                   call LVT_readvar_restart(ftn,&
                        stats%count_model_rfv_total(:,k,l))
                enddo
             enddo
             if(LVT_metrics%rfv%computeSC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nasc         
                      call LVT_readvar_restart(ftn,&
                           stats%model_rfv_asc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%count_model_rfv_asc(:,k,l))
                   enddo
                enddo
             endif
             if(LVT_metrics%rfv%computeADC.eq.1) then 
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%nadc         
                      call LVT_readvar_restart(ftn,&
                           stats%model_rfv_adc(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%count_model_rfv_adc(:,k,l))
                   enddo
                enddo
             endif
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   do k=1,obs%vlevels
                      do l=1,LVT_rc%strat_nlevels
                         call LVT_readvar_restart(ftn,&
                              stats%obs_rfv_total(:,k,l))
                         call LVT_readvar_restart(ftn,&
                              stats%count_obs_rfv_total(:,k,l))
                      enddo
                   enddo
                   
                   if(LVT_metrics%rfv%computeSC.eq.1) then 
                      do k=1,obs%vlevels
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%obs_rfv_asc(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%count_obs_rfv_asc(:,k,l))
                         enddo
                      enddo
                   endif
                   if(LVT_metrics%rfv%computeADC.eq.1) then 
                      do k=1,obs%vlevels
                         do l=1,LVT_rc%nasc
                            call LVT_readvar_restart(ftn,&
                                 stats%obs_rfv_adc(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%count_obs_rfv_adc(:,k,l))
                         enddo
                      enddo
                   endif
                endif
             endif
          endif
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
    end do
#endif
  end subroutine LVT_readrestart_RFV


end module LVT_RFVMod
