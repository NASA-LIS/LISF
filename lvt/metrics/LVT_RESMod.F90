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
! !MODULE: LVT_RESMod
! \label(LVT_RESMod)
!
! !INTERFACE:
module LVT_RESMod
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

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the computations required to compute resilience based
!  on a specific variable. Resilience is defined to be an indicator of how 
!  quickly a system returns to a satisfactory condition after an unsatisfactory
!  condition. 
! 
!      RES = number of times a satisfactory codnition follows an unsatisfactory one/
!             total number of unsatisfactory conditions
!
!  A satisfactory condition is defined as a case when the variable is 
!  greater than the min threshold and less than the max threshold. 
! 
!  Reference on the metric: 
!  Thomas et al. 2017, Global assessment of groundwater sustainability based
!  on storage anomalies, GRL. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  26 Dec 17    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initRES
  public :: LVT_diagnoseRES
  public :: LVT_computeRES
  public :: LVT_writeMetric_RES
  public :: LVT_resetMetric_RES
  public :: LVT_writerestart_RES
  public :: LVT_readrestart_RES
!EOP
  
  private

contains
  subroutine LVT_initRES(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
! 
!EOP
    integer                 :: m

    allocate(stats%res(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%res(m)%model_value_total(&
               LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%res(m)%model_value_total = 0.0
          allocate(stats%res(m)%count_model_value_total(&
               LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%res(m)%count_model_value_total = 0

          allocate(stats%res(m)%model_value_prev(&
               LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%res(m)%model_value_prev = 0.0

          allocate(stats%res(m)%model_value_ci(&
               selectNlevs(1),LVT_rc%strat_nlevels))
          stats%res(m)%model_value_ci = LVT_rc%udef
          
          if(LVT_rc%obssource(2).ne."none") then 
             if(selectNlevs(2).ge.1) then 
                allocate(stats%res(m)%obs_value_total(&
                     LVT_rc%ngrid, selectNlevs(2), &
                     LVT_rc%strat_nlevels))
                stats%res(m)%obs_value_total = 0.0
                allocate(stats%res(m)%count_obs_value_total(&
                     LVT_rc%ngrid, selectNlevs(2), &
                     LVT_rc%strat_nlevels))
                stats%res(m)%count_obs_value_total = 0

                allocate(stats%res(m)%obs_value_prev(&
                     LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%strat_nlevels))
                stats%res(m)%obs_value_prev = 0.0

                allocate(stats%res(m)%obs_value_ci(&
                     selectNlevs(2),LVT_rc%strat_nlevels))
                stats%res(m)%obs_value_ci = LVT_rc%udef

             endif
          endif
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    if(LVT_rc%obssource(2).ne."none") then 
       metric%obsData = .true. 
    else
       metric%obsData = .false. 
    endif
    metric%stdevFlag = .false. 

  end subroutine LVT_initRES

!BOP
! 
! !ROUTINE: LVT_diagnoseRES
! \label{LVT_diagnoseRES}
!
! !INTERFACE: 
  subroutine LVT_diagnoseRES(pass)
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
!   calculating the res of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelRES](\ref{diagnoseSingleModelRES})
!     updates the res computation for a single variable 
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

    if(pass.eq.LVT_metrics%res%npass) then
       if(LVT_metrics%res%selectOpt.eq.1.or.&
            LVT_metrics%res%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleRES(model,obs,stats,&
                  LVT_metrics%res)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseRES

!BOP
! 
! !ROUTINE: diagnoseSingleRES
! \label{diagnoseSingleRES}
!
! !INTERFACE: 
  subroutine diagnoseSingleRES(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the res computation of the 
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
    integer    :: t,k,n,m,m_k,o_k,tind

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             do m=1,LVT_rc%nensem
                m_k = k+model%startNlevs -1
                if(model%count(t,m,m_k).gt.0) then
                   if(model%value(t,m,m_k).ne.LVT_rc%udef) then
                      
                      if((model%value(t,m,m_k).gt.metric%minthreshold.and.&
                           model%value(t,m,m_k).lt.metric%maxthreshold)) then
                         !unsatisfactory condition                         
                         stats%res(m)%count_model_value_total(t,k,1) = &
                              stats%res(m)%count_model_value_total(t,k,1) + 1
                      endif
                      !satisfactory condition following an unsatisfactory one 
                      if((stats%res(m)%model_value_prev(t,k,1).gt.0).and.&
                           (.not.(model%value(t,m,m_k).gt.metric%minthreshold.and.&
                           model%value(t,m,m_k).lt.metric%maxthreshold))) then 
                      
                         stats%res(m)%model_value_total(t,k,1) = &
                              stats%res(m)%model_value_total(t,k,1) + 1
                      endif

                      if((model%value(t,m,m_k).gt.metric%minthreshold.and.&
                           model%value(t,m,m_k).lt.metric%maxthreshold)) then
                         stats%res(m)%model_value_prev(t,k,1) = 1.0
                      else
                         stats%res(m)%model_value_prev(t,k,1) = 0.0
                      endif
                   endif
                
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,m,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         if(model%value(t,m,m_k).ne.LVT_rc%udef) then  
                            if((model%value(t,m,m_k).gt.metric%minthreshold.and.&
                                 model%value(t,m,m_k).lt.metric%maxthreshold)) then
                               !unsatisfactory condition                         
                               stats%res(m)%count_model_value_total(t,k,2) = &
                                    stats%res(m)%count_model_value_total(t,k,2) + 1
                            endif
                            !satisfactory condition following an unsatisfactory one 
                            if((stats%res(m)%model_value_prev(t,k,2).gt.0).and.&
                                 (.not.(model%value(t,m,m_k).gt.metric%minthreshold.and.&
                                 model%value(t,m,m_k).lt.metric%maxthreshold))) then 
                               
                               stats%res(m)%model_value_total(t,k,2) = &
                                    stats%res(m)%model_value_total(t,k,2) + 1
                            endif
                            
                            if((model%value(t,m,m_k).gt.metric%minthreshold.and.&
                                 model%value(t,m,m_k).lt.metric%maxthreshold)) then
                               stats%res(m)%model_value_prev(t,k,2) = 1.0
                            else
                               stats%res(m)%model_value_prev(t,k,2) = 0.0
                            endif

                         endif
                      elseif(LVT_stats%strat_var(t,m,k).le.&
                           LVT_rc%strat_var_threshold) then
                         if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                            if((model%value(t,m,m_k).gt.metric%minthreshold.and.&
                                 model%value(t,m,m_k).lt.metric%maxthreshold)) then
                               !unsatisfactory condition                         
                               stats%res(m)%count_model_value_total(t,k,3) = &
                                    stats%res(m)%count_model_value_total(t,k,3) + 1
                            endif
                            !satisfactory condition following an unsatisfactory one 
                            if((stats%res(m)%model_value_prev(t,k,3).gt.0).and.&
                                 (.not.(model%value(t,m,m_k).gt.metric%minthreshold.and.&
                                 model%value(t,m,m_k).lt.metric%maxthreshold))) then 
                               
                               stats%res(m)%model_value_total(t,k,3) = &
                                    stats%res(m)%model_value_total(t,k,3) + 1
                            endif
                            
                            if((model%value(t,m,m_k).gt.metric%minthreshold.and.&
                                 model%value(t,m,m_k).lt.metric%maxthreshold)) then
                               stats%res(m)%model_value_prev(t,k,3) = 1.0
                            else
                               stats%res(m)%model_value_prev(t,k,3) = 0.0
                            endif
                         endif
                      endif
                   endif
                endif
             end do
          enddo
          do k=1,obs%selectNlevs
             do m=1,LVT_rc%nensem
                o_k = k+obs%startNlevs -1
                if(LVT_rc%obssource(2).ne."none") then
                   if(obs%selectNlevs.ge.1) then
                      if(obs%count(t,m,o_k).gt.0) then 

                         if((obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                              obs%value(t,m,o_k).lt.metric%maxthreshold)) then
                            !unsatisfactory condition                         
                            stats%res(m)%count_obs_value_total(t,k,3) = &
                                 stats%res(m)%count_obs_value_total(t,k,3) + 1
                         endif
                         !satisfactory condition following an unsatisfactory one 
                         if((stats%res(m)%obs_value_prev(t,k,1).gt.0).and.&
                              (.not.(obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                              obs%value(t,m,o_k).lt.metric%maxthreshold))) then 
                            
                            stats%res(m)%obs_value_total(t,k,1) = &
                                 stats%res(m)%obs_value_total(t,k,1) + 1
                         endif
                         
                         if((obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                              obs%value(t,m,o_k).lt.metric%maxthreshold)) then
                            stats%res(m)%obs_value_prev(t,k,1) = 1.0
                         else
                            stats%res(m)%obs_value_prev(t,k,1) = 0.0
                         endif
                         
                         
                         if(LVT_rc%strat_nlevels.gt.1) then 

                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then
                               
                               if((obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                                    obs%value(t,m,o_k).lt.metric%maxthreshold)) then
                                  !unsatisfactory condition                         
                                  stats%res(m)%count_obs_value_total(t,k,2) = &
                                       stats%res(m)%count_obs_value_total(t,k,2) + 1
                               endif
                               !satisfactory condition following an unsatisfactory one 
                               if((stats%res(m)%obs_value_prev(t,k,2).gt.0).and.&
                                    (.not.(obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                                    obs%value(t,m,o_k).lt.metric%maxthreshold))) then 
                                  
                                  stats%res(m)%obs_value_total(t,k,2) = &
                                       stats%res(m)%obs_value_total(t,k,2) + 1
                               endif
                               
                               if((obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                                    obs%value(t,m,o_k).lt.metric%maxthreshold)) then
                                  stats%res(m)%obs_value_prev(t,k,2) = 1.0
                               else
                                  stats%res(m)%obs_value_prev(t,k,2) = 0.0
                               endif
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then
                               if(obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  

                                  if((obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                                       obs%value(t,m,o_k).lt.metric%maxthreshold)) then
                                     !unsatisfactory condition                         
                                     stats%res(m)%count_obs_value_total(t,k,3) = &
                                          stats%res(m)%count_obs_value_total(t,k,3) + 1
                                  endif
                                  !satisfactory condition following an unsatisfactory one 
                                  if((stats%res(m)%obs_value_prev(t,k,3).gt.0).and.&
                                       (.not.(obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                                       obs%value(t,m,o_k).lt.metric%maxthreshold))) then 
                            
                                     stats%res(m)%obs_value_total(t,k,3) = &
                                          stats%res(m)%obs_value_total(t,k,3) + 1
                                  endif
                                  
                                  if((obs%value(t,m,o_k).gt.metric%minthreshold.and.&
                                       obs%value(t,m,o_k).lt.metric%maxthreshold)) then
                                     stats%res(m)%obs_value_prev(t,k,3) = 1.0
                                  else
                                     stats%res(m)%obs_value_prev(t,k,3) = 0.0
                                  endif
                               endif
                            endif
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine diagnoseSingleRES

!BOP
! 
! !ROUTINE: LVT_computeRES
! \label{LVT_computeRES}
!
! !INTERFACE: 
  subroutine LVT_computeRES(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the res values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelRES](\ref{computeSingleModelRES})
!     computes the res values for a single variable
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
    
    if(pass.eq.LVT_metrics%res%npass) then
       if(LVT_metrics%res%selectOpt.eq.1.or.&
            LVT_metrics%res%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%res%timeOpt.eq.1.and.&
                  LVT_metrics%res%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%res%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%res%ftn_ts_loc(i,1),200,advance='no') &
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
             
             call computeSingleRES(alarm,model,obs,stats,&
                  LVT_metrics%res)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%res%timeOpt.eq.1.and.&
                  LVT_metrics%res%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%res%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%res%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif
             endif
          endif
       endif
    endif
  end subroutine LVT_computeRES
  

!BOP
! 
! !ROUTINE: computeSingleRES
! \label{computeSingleRES}
!
! !INTERFACE: 
  subroutine computeSingleRES(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the res values
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
    real,    allocatable :: model_value_avg(:,:,:)
    real,    allocatable :: obs_value_avg(:,:,:)   
       
    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%res(m)%count_model_value_total(t,k,l).gt.&
                           LVT_rc%obsCountThreshold) then 
                         
                         stats%res(m)%model_value_total(t,k,l) = &
                              stats%res(m)%model_value_total(t,k,l)/&
                              stats%res(m)%count_model_value_total(t,k,l)           
                      else
                         stats%res(m)%model_value_total(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
                do k=1,obs%selectNlevs
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then
                         do l=1,LVT_rc%strat_nlevels                      
                            if(stats%res(m)%count_obs_value_total(t,k,l).gt.&
                                 LVT_rc%obsCountThreshold) then 
                               stats%res(m)%obs_value_total(t,k,l) = &
                                    stats%res(m)%obs_value_total(t,k,l)/&
                                    stats%res(m)%count_obs_value_total(t,k,l)       
                            else
                               stats%res(m)%obs_value_total(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                   endif
                enddo
             enddo
          enddo

          do m=1,LVT_rc%nensem          
             do k=1,model%selectNlevs                
                do l=1, LVT_rc%strat_nlevels
                   call LVT_computeCI(stats%res(m)%model_value_total(:,k,l),&
                        LVT_rc%ngrid,&
                        LVT_rc%pval_CI,stats%res(m)%model_value_ci(k,l))
                enddo
             enddo
             if(LVT_rc%obssource(2).ne."none") then 
                do k=1,obs%selectNlevs
                   do l=1, LVT_rc%strat_nlevels
                      call LVT_computeCI(stats%res(m)%obs_value_total(:,k,l),&
                           LVT_rc%ngrid,&
                           LVT_rc%pval_CI,stats%res(m)%obs_value_ci(k,l))
                   enddo
                enddo
             endif
          enddo

!----------------------------------------------------------------------
!  External data based stratification 
!----------------------------------------------------------------------
          if(LVT_rc%data_based_strat.eq.1) then 

             allocate(model_value_avg(LVT_rc%ngrid,&
                  model%vlevels,&
                  LVT_rc%strat_nlevels))
             model_value_avg = 0.0

             do m=1,LVT_rc%nensem        
                do t=1,LVT_rc%ngrid                      
                   do k=1,model%selectNlevs                
                      do l=1, LVT_rc%strat_nlevels
                         model_value_avg(t,k,l) = &
                              model_value_avg(t,k,l) + & 
                              stats%res(m)%model_value_total(t,k,l)
                      enddo
                   enddo
                enddo
             enddo
             
             do t=1,LVT_rc%ngrid                      
                do k=1,model%selectNlevs                
                   do l=1, LVT_rc%strat_nlevels
                      model_value_avg(t,k,l) = &
                           model_value_avg(t,k,l)/LVT_rc%nensem
                   enddo
                enddo
             enddo
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   
                allocate(obs_value_avg(LVT_rc%ngrid,&
                     obs%vlevels,&
                     LVT_rc%strat_nlevels))
                obs_value_avg = 0.0

                do m=1,LVT_rc%nensem        
                   do t=1,LVT_rc%ngrid
                      do k=1,obs%selectNlevs                
                         do l=1, LVT_rc%strat_nlevels
                            obs_value_avg(t,k,l) = &
                                 obs_value_avg(t,k,l) + & 
                                 stats%res(m)%obs_value_total(t,k,l)
                         enddo
                      enddo
                   enddo
                enddo
                do t=1,LVT_rc%ngrid
                   do k=1,obs%selectNlevs                
                      do l=1, LVT_rc%strat_nlevels                
                         obs_value_avg(t,k,l) = &
                              obs_value_avg(t,k,l)/LVT_rc%nensem
                      enddo
                   enddo
                enddo
                call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                     LVT_rc%ngrid, model_value_avg, &
                     LVT_rc%ngrid, obs_value_avg)
                   
                deallocate(obs_value_avg)
             else
                   
                call LVT_writeDataBasedStrat(model,obs,stats,metric,&
                     LVT_rc%ngrid,model_value_avg)            
                
             endif
             deallocate(model_value_avg)
          endif

       endif
    endif
  end subroutine computeSingleRES

!BOP
! 
! !ROUTINE: LVT_writeMetric_RES
! \label(LVT_writeMetric_RES)
!
! !INTERFACE:
  subroutine LVT_writeMetric_RES(pass,final,vlevels,stats,obs)
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

    real,    allocatable    :: model_value_total(:,:)
    integer, allocatable    :: count_model_value_total(:,:)
    real,    allocatable    :: model_value_ci(:)

    real,    allocatable    :: obs_value_total(:,:)
    integer, allocatable    :: count_obs_value_total(:,:)
    real,    allocatable    :: obs_value_ci(:)

    real,    allocatable    :: model_value_ts(:,:)
    integer, allocatable    :: count_model_value_ts(:,:)

    real,    allocatable    :: obs_value_ts(:,:)
    integer, allocatable    :: count_obs_value_ts(:,:)

    real,    allocatable    :: model_value_adc(:,:,:)
    real,    allocatable    :: model_value_asc(:,:,:)

    real,    allocatable    :: obs_value_adc(:,:,:)
    real,    allocatable    :: obs_value_asc(:,:,:)

    if(pass.eq.LVT_metrics%res%npass) then
       if(stats%selectOpt.eq.1) then
          
          allocate(model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
          allocate(count_model_value_total(LVT_rc%ngrid,LVT_rc%nensem))
          allocate(model_value_ci(LVT_rc%nensem))
          
          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
             allocate(obs_value_total(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_obs_value_total(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(obs_value_ci(LVT_rc%nensem))
          endif
          
          do k=1,vlevels
             do l=1,LVT_rc%strat_nlevels
                do m=1,LVT_rc%nensem
                   model_value_total(:,m) = &
                        stats%res(m)%model_value_total(:,k,l)
                   count_model_value_total(:,m) = & 
                        stats%res(m)%count_model_value_total(:,k,l)
                   model_value_ci(m) = stats%res(m)%model_value_ci(k,l)
                   
                enddo
                if(LVT_metrics%res%selectOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%res%ftn_total, &
                        model_value_total(:,:),&
                        stats%vid_total(LVT_RESid,1),k)
                   call LVT_writevar_gridded(LVT_metrics%res%ftn_total, &
                        real(count_model_value_total(:,:)),&
                        stats%vid_count_total(LVT_RESid,1),k)
                   
                   if(LVT_rc%obssource(2).ne."none".and.&
                        obs%selectNlevs.ge.1) then 
                      call LVT_writevar_gridded(LVT_metrics%res%ftn_total, &
                           obs_value_total(:,:),&
                           stats%vid_total(LVT_RESid,2),k)
                      call LVT_writevar_gridded(LVT_metrics%res%ftn_total, &
                           real(count_obs_value_total(:,:)),&
                           stats%vid_count_total(LVT_RESid,2),k)
                   endif
                   
                   call LVT_writeSummaryStats(&
                        LVT_metrics%res%ftn_summ,&
                        l,&
                        LVT_metrics%res%short_name,&
                        LVT_rc%ngrid,&
                        model_value_total(:,:), &
                        count_model_value_total(:,:),&
                        stats%standard_name,&
                        model_value_ci(:))
                   if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                      call LVT_writeSummaryStats(&
                           LVT_metrics%res%ftn_summ,&
                           l,&
                           LVT_metrics%res%short_name,&
                           LVT_rc%ngrid,&
                           obs_value_total(:,:), &
                           count_obs_value_total(:,:),&
                           "DS2_"//trim(stats%standard_name),&
                           obs_value_ci(:))
                   endif
                endif
             enddo
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

  end subroutine LVT_writeMetric_RES

!BOP
! 
! !ROUTINE: LVT_resetMetric_RES
! \label(LVT_resetMetric_RES)
!
! !INTERFACE:
  subroutine LVT_resetMetric_RES(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical                :: alarm
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

  end subroutine LVT_resetMetric_RES

!BOP
! 
! !ROUTINE: LVT_writerestart_RES
! 
! !INTERFACE:
  subroutine LVT_writerestart_RES(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for RES metric computations
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
       if(LVT_metrics%res%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_writevar_restart(ftn,&
                           stats%res(m)%model_value_total(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%res(m)%count_model_value_total(:,k,l))
                   enddo
                enddo

                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            call LVT_writevar_restart(ftn,&
                                 stats%res(m)%obs_value_total(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%res(m)%count_obs_value_total(:,k,l))
                         enddo
                      enddo
                   
                   endif
                endif
             enddo
          endif
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
    end do
  end subroutine LVT_writerestart_RES


!BOP
! 
! !ROUTINE: LVT_readrestart_RES
! 
! !INTERFACE:
  subroutine LVT_readrestart_RES(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for RES metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    integer              :: k,l,m,index
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%res%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_readvar_restart(ftn,&
                           stats%res(m)%model_value_total(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%res(m)%count_model_value_total(:,k,l))
                   enddo
                enddo

                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            call LVT_readvar_restart(ftn,&
                                 stats%res(m)%obs_value_total(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%res(m)%count_obs_value_total(:,k,l))
                         enddo
                      enddo

                   endif
                endif
             end do
          endif
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
    end do
  end subroutine LVT_readrestart_RES


end module LVT_RESMod
