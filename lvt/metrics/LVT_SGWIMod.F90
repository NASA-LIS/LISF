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
! !MODULE: LVT_SGWIMod
! \label(LVT_SGWIMod)
!
! !INTERFACE:
module LVT_SGWIMod
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_timeMgrMod
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
!  !DESCRIPTION: 
!   This module handles the Standardized soil water index (SGWI)
!   metric based on the root zone soil moisture outputs in LIS. 
!   
!   Reference: McKee, T. B., N. J. Doesken, and J. Kleist, 1993:
!   The re- lationship of drought frequency and duration to time scales. 
!   Preprints, Eighth Conf. on Applied Climatol- ogy, Anaheim, CA, 
!   Amer. Meteor. Soc., 179–184.
! 
! !NOTES: 
!   The standard practice is to compute SGWI on monthly scale. Since the 
!   code needs to store the precip values in time to compute SGWI, 
!   currently LVT only supports SGWI computations when the time 
!   averaging interval is set to a multiple of a month (2592000 seconds). 
!
!   Only monthly SGWI computations have been tested so far
!
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  22 Mar 2012    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initSGWI
  public :: LVT_diagnoseSGWI
  public :: LVT_computeSGWI
  public :: LVT_writeMetric_SGWI
  public :: LVT_resetMetric_SGWI
  public :: LVT_writerestart_SGWI
  public :: LVT_readrestart_SGWI

!EOP

  type, public :: sgwidec

     real,    allocatable :: model_rzsm_final(:,:)
     real,    allocatable :: obs_rzsm_final(:,:)

     real,    allocatable :: model_mu(:,:)
     real,    allocatable :: model_sigma(:,:)
     real,    allocatable :: model_sgwi(:)

     real,    allocatable :: obs_mu(:,:)
     real,    allocatable :: obs_sigma(:,:)
     real,    allocatable :: obs_sgwi(:)

     real             :: model_sgwi_ci(1,1)
     real             :: obs_sgwi_ci(1,1)

     integer          :: nsize_total
     integer          :: nsize_season
     integer          :: nasc
  end type sgwidec

  type(sgwidec), allocatable :: LVT_sgwi_struc(:)
  private
  integer, parameter :: maxiter = 100
  real,    parameter :: epsilon = 3.0e-7


contains
!BOP
!
! !ROUTINE: LVT_initSGWI
! \label{LVT_initSGWI}
! 
! !INTERFACE: 
  subroutine LVT_initSGWI(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This routine initializes the data structures required
!  for the SGWI computations. 
! 
!EOP

    type(ESMF_Time)         :: startTime
    type(ESMF_Time)         :: stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: rc

    integer                 :: m 

    allocate(LVT_sgwi_struc(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(LVT_metrics%sgwi%selectOpt.eq.1.or.&
            LVT_metrics%sgwi%timeOpt.eq.1) then 
          
          if(mod(LVT_rc%tavgInterval,2592000).ne.0) then 
             write(LVT_logunit,*) 'SGWI values are computed at a monthly timescale'
             write(LVT_logunit,*) 'Please set the time averaging interval to a '
             write(LVT_logunit,*) 'multiple of a month (2592000 seconds)'
             write(LVT_logunit,*) 'program stopping...'
             call LVT_endrun()
          endif
          
          if(LVT_rc%tavgInterval/2592000.eq.1) then 
             LVT_sgwi_struc(m)%nasc = 12
          elseif(LVT_rc%tavgInterval/2592000.eq.3) then 
             LVT_sgwi_struc(m)%nasc = 4
          elseif(LVT_rc%tavgInterval/2592000.eq.6) then 
             LVT_sgwi_struc(m)%nasc = 2
          elseif(LVT_rc%tavgInterval/2592000.eq.12) then 
             LVT_sgwi_struc(m)%nasc = 1
          else
             write(LVT_logunit,*) 'The interval ',LVT_rc%tavgInterval,&
                  ' is not supported for SGWI calculations '
             write(LVT_logunit,*) 'Program stopping...'
             call LVT_endrun()
          endif
          
          call computeNmonths(LVT_sgwi_struc(m)%nsize_total, LVT_rc%tavgInterval)
          LVT_sgwi_struc(m)%nsize_season = &
               LVT_sgwi_struc(m)%nsize_total/LVT_sgwi_struc(m)%nasc
          LVT_sgwi_struc(m)%nsize_total = &
               LVT_sgwi_struc(m)%nasc*LVT_sgwi_struc(m)%nsize_season
          
          if(metric%selectOpt.eq.1) then 
             allocate(LVT_sgwi_struc(m)%model_rzsm_final(LVT_LIS_rc(1)%ntiles,&
                  LVT_sgwi_struc(m)%nsize_total))
             LVT_sgwi_struc(m)%model_rzsm_final = 0.0
             
             allocate(LVT_sgwi_struc(m)%model_mu(LVT_LIS_rc(1)%ntiles,&
                  LVT_sgwi_struc(m)%nasc))
             allocate(LVT_sgwi_struc(m)%model_sigma(LVT_LIS_rc(1)%ntiles,&
                  LVT_sgwi_struc(m)%nasc))
             allocate(LVT_sgwi_struc(m)%model_sgwi(LVT_LIS_rc(1)%ntiles))
             LVT_sgwi_struc(m)%model_sgwi = LVT_rc%udef
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(selectNlevs(2).ge.1) then 
                   allocate(LVT_sgwi_struc(m)%obs_rzsm_final(LVT_rc%ngrid,&
                        LVT_sgwi_struc(m)%nsize_total))
                   LVT_sgwi_struc(m)%obs_rzsm_final = 0.0
                   
                   allocate(LVT_sgwi_struc(m)%obs_mu(LVT_rc%ngrid,&
                        LVT_sgwi_struc(m)%nasc))
                   allocate(LVT_sgwi_struc(m)%obs_sigma(LVT_rc%ngrid,&
                        LVT_sgwi_struc(m)%nasc))
                   allocate(LVT_sgwi_struc(m)%obs_sgwi(LVT_rc%ngrid))
                   LVT_sgwi_struc(m)%obs_sgwi = LVT_rc%udef
                   
                endif
             endif
          endif
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------
    if(LVT_rc%startmode.eq."coldstart") then           
       metric%npass = 2   
    else
       metric%npass = 1
    endif
    
    metric%obsData = .true. 
    metric%stdevFlag = .false. 

  end subroutine LVT_initSGWI

!BOP
! 
! !ROUTINE: LVT_diagnoseSGWI
! \label{LVT_diagnoseSGWI}
!
! !INTERFACE: 
  subroutine LVT_diagnoseSGWI(pass)
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
!   calculating the sgwi of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelSGWI](\ref{diagnoseSingleModelSGWI})
!     updates the sgwi computation for a single variable. This routine
!     stores the precip values to be used later for computing SGWI, 
!     during the first pass through the data. 
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
       if(LVT_metrics%sgwi%selectOpt.eq.1.or.&
            LVT_metrics%sgwi%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))   
             
             call diagnoseSingleSGWI(model, obs, stats,&
                  LVT_metrics%sgwi)
             
             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseSGWI


!BOP
! 
! !ROUTINE: diagnoseSingleSGWI
! \label{diagnoseSingleSGWI}
!
! !INTERFACE: 
  subroutine diagnoseSingleSGWI(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the sgwi computation of the 
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
    integer                 :: t,k,tind
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_TimeInterval) :: ts
    integer                 :: tindex_f,tindex_t
    integer                 :: rc,m,m_k,o_k

    if(stats%selectOpt.eq.1.and.&
         (model%short_name.eq."TWS").and.&
         model%selectNlevs.ge.1) then        

       call getMonthlyTimeIndex(tindex_f, LVT_rc%tavgInterval)
       
       do t=1,LVT_LIS_rc(1)%ntiles
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(model%count(t,m,m_k).gt.0) then
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         LVT_sgwi_struc(m)%model_rzsm_final(t,tindex_f) = &
                              model%value(t,m,m_k)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo

       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                if(LVT_rc%obssource(2).ne."none") then
                   if(obs%selectNlevs.ge.1) then
                      if(obs%count(t,m,o_k).gt.0) then 
                         if(metric%selectOpt.eq.1.and.&
                              obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            LVT_sgwi_struc(m)%obs_rzsm_final(t,tindex_f) = &
                                 obs%value(t,m,o_k)
                         endif
                         
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
  end subroutine diagnoseSingleSGWI

!BOP
! 
! !ROUTINE: LVT_computeSGWI
! \label{LVT_computeSGWI}
!
! !INTERFACE: 
  subroutine LVT_computeSGWI(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the SGWI values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelSGWI](\ref{computeSingleModelSGWI})
!     computes the SGWI values for a single variable
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
    integer         :: i,m
    integer         :: index
    type(ESMF_Time) :: currTime
    integer         :: rc
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(LVT_rc%startmode.eq."coldstart") then 
       if(pass.eq.1) then 
          if(LVT_metrics%sgwi%selectOpt.eq.1.or.LVT_metrics%sgwi%timeOpt.eq.1) then 
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
          
             do while(associated(model))             
                call computeSingleSGWIparams(alarm,&
                     model,obs,stats,LVT_metrics%sgwi)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
          endif
       elseif(pass.eq.2) then        
          if(LVT_metrics%sgwi%selectOpt.eq.1.or.LVT_metrics%sgwi%timeOpt.eq.1) then 
             if(alarm) then 
                if(LVT_metrics%sgwi%timeOpt.eq.1.and.&
                     LVT_metrics%sgwi%extractTS.eq.1) then 
                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%sgwi%ftn_ts_loc(i,m),200,advance='no') &
                                 LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                                 LVT_rc%hr,'',LVT_rc%mn, '' 
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sgwi%ftn_ts_loc(i,1),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   endif
                   
                endif
             endif
200          format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)
             
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model)) 
                call computeSingleSGWI(alarm,&
                     model,obs,stats,LVT_metrics%sgwi)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
             
             if(alarm) then 
                if(LVT_metrics%sgwi%timeOpt.eq.1.and.&
                     LVT_metrics%sgwi%extractTS.eq.1) then 

                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%sgwi%ftn_ts_loc(i,m),fmt='(a1)') ''
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sgwi%ftn_ts_loc(i,1),fmt='(a1)') ''
                      enddo
                   endif
                   
                endif
             endif
          endif
       endif
    elseif(LVT_rc%startmode.eq."restart") then 
       if(LVT_metrics%sgwi%selectOpt.eq.1.or.LVT_metrics%sgwi%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%sgwi%timeOpt.eq.1.and.&
                  LVT_metrics%sgwi%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sgwi%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%sgwi%ftn_ts_loc(i,1),200,advance='no') &
                           LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                           LVT_rc%hr,'',LVT_rc%mn, '' 
                   enddo
                endif
             endif
          endif

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)
          
          do while(associated(model)) 
             call computeSingleSGWI(alarm,&
                  model,obs,stats,LVT_metrics%sgwi)
             
             model => model%next
             obs => obs%next
             stats => stats%next
             
          enddo
          
          if(alarm) then 
             if(LVT_metrics%sgwi%timeOpt.eq.1.and.&
                  LVT_metrics%sgwi%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sgwi%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%sgwi%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif
             endif
          endif
       endif
    endif

  end subroutine LVT_computeSGWI
  

!BOP
! 
! !ROUTINE: computeSingleSGWIparams
! \label{computeSingleSGWIparams}
!
! !INTERFACE: 
  subroutine computeSingleSGWIparams(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the SGWI values for each grid cell. 
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

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            (model%short_name.eq."TWS").and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_LIS_rc(1)%ntiles
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_sgwi_struc(m)%nasc
                      call sgwi_normal(t,l,m,LVT_sgwi_struc(m)%nsize_total,&
                           LVT_sgwi_struc(m)%nsize_season,&
                           LVT_sgwi_struc(m)%model_rzsm_final(t,:),&
                           LVT_sgwi_struc(m)%model_mu(t,l),&
                           LVT_sgwi_struc(m)%model_sigma(t,l))
                   enddo
                enddo
             enddo
          enddo
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then 
                         do l=1,LVT_sgwi_struc(m)%nasc
                            call sgwi_normal(t,l,m,LVT_sgwi_struc(m)%nsize_total,&
                                 LVT_sgwi_struc(m)%nsize_season,&
                                 LVT_sgwi_struc(m)%obs_rzsm_final(t,:),&
                                 LVT_sgwi_struc(m)%obs_mu(t,l),&
                                 LVT_sgwi_struc(m)%obs_sigma(t,l))
                         enddo
                      endif
                   endif
                enddo
             enddo
          enddo
       endif
    endif
  end subroutine computeSingleSGWIparams

!BOP
! !ROUTINE:sgwi_normal
! \label{sgwi_normal}
! 
! !INTERFACE: 
  subroutine sgwi_normal(t,l,m,nsize,nsize_season,rzsm, mu, sigma)

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: t
    integer, intent(in) :: l
    integer, intent(in) :: m
    integer, intent(in) :: nsize
    integer, intent(in) :: nsize_season
    real,    intent(in) :: rzsm(nsize)
    real                :: mu
    real                :: sigma
!
! !DESCRIPTION: 
!   This subroutine fits a gamma distribution to the given precipitation
!   data and returns the parameters of the distribution
! 
!   \begin{description}
!    \item[t]
!      index of the grid/tile
!    \item[l]
!      index of the season/month
!    \item[nsize]
!      size of the precipitation input array
!    \item[nsize\_season]
!      size of the seasonally stratified precipitation input
!    \item[rzsm]
!      array of precipitation values
!    \item[mu]
!      mu parameter of the fitted distribution
!    \item[sigma]
!      sigma parameter of the fitted distribution
!   \end{description}  
!EOP  
    integer             :: i,k
    real                :: rzsm_season(nsize_season)
    
    rzsm_season = 0.0
    k = 0 
    do i=l,nsize,LVT_sgwi_struc(m)%nasc
       k = k+1
       rzsm_season(k) = rzsm(i)
    enddo

    call normal_fit(nsize_season, rzsm_season, &
         mu, sigma)

  end subroutine sgwi_normal

!BOP
! 
! !ROUTINE: computeSingleSGWI
! \label{computeSingleSGWI}
!
! !INTERFACE: 
  subroutine computeSingleSGWI(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the SGWI values for each grid cell. 
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

    integer  :: t,l,k,i,m
    integer  :: tindex_f, tind
    real     :: model_sgwi(LVT_LIS_rc(1)%ntiles,1,1)
    real     :: obs_sgwi(LVT_rc%ngrid,1,1)
    integer  :: count_model_sgwi(LVT_LIS_rc(1)%ntiles,1,1)
    integer  :: count_obs_sgwi(LVT_rc%ngrid,1,1)

    if(stats%selectOpt.eq.1.and.&
         (model%short_name.eq."TWS").and.&
         model%selectNlevs.ge.1) then 

       call getMonthlyTimeIndex(tindex_f, LVT_rc%tavgInterval)
       call getSeasonalTimeIndex(tind,LVT_rc%tavgInterval)

       do t=1,LVT_LIS_rc(1)%ntiles
          do m=1,LVT_rc%nensem
             call sgwi_zval(t,&
                  LVT_sgwi_struc(m)%model_rzsm_final(t,tindex_f),&
                  LVT_sgwi_struc(m)%model_mu(t,tind), &
                  LVT_sgwi_struc(m)%model_sigma(t,tind), &
                  LVT_sgwi_struc(m)%model_sgwi(t))
          enddo
       enddo

       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             if(LVT_rc%obssource(2).ne."none") then 
                if(obs%selectNlevs.ge.1) then 
                   call sgwi_zval(t,&
                        LVT_sgwi_struc(m)%obs_rzsm_final(t,tindex_f),&
                        LVT_sgwi_struc(m)%obs_mu(t,tind), &
                        LVT_sgwi_struc(m)%obs_sigma(t,tind), &
                        LVT_sgwi_struc(m)%obs_sgwi(t))
                   
                endif
             endif
          enddo
       enddo

       count_model_sgwi = 1
       count_obs_sgwi = 1       
       if(metric%extractTS.eq.1) then 
          if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
             do m=1,LVT_rc%nensem
                model_sgwi(:,1,1) = LVT_sgwi_struc(m)%model_sgwi(:)
                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   obs_sgwi(:,1,1) = LVT_sgwi_struc(m)%model_sgwi(:)
                   call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                        model,&
                        LVT_LIS_rc(1)%ntiles,&
                        model_sgwi,&
                        count_model_sgwi,&
                        LVT_rc%ngrid,&
                        obs_sgwi,&
                        count_obs_sgwi)
                   
                else
                   call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                        model,&
                        LVT_LIS_rc(1)%ntiles,&
                        model_sgwi,&
                        count_model_sgwi)
                endif
             enddo
          else
!TBD
          endif
       endif


       do m=1,LVT_rc%nensem
          call LVT_computeCI(LVT_sgwi_struc(m)%model_sgwi(:),LVT_LIS_rc(1)%ntiles,&
               LVT_rc%pval_CI, LVT_sgwi_struc(m)%model_sgwi_ci(1,1))
          
          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
             call LVT_computeCI(LVT_sgwi_struc(m)%obs_sgwi(:),LVT_rc%ngrid,&
                  LVT_rc%pval_CI,LVT_sgwi_struc(m)%obs_sgwi_ci(1,1))
          endif
       enddo

#if 0 
       call LVT_writeDataBasedStrat(model,obs,stats,metric,&
            LVT_LIS_rc(1)%ntiles, model_sgwi)

       if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_LIS_rc(1)%ntiles, model_sgwi, &
               LVT_rc%ngrid, obs_sgwi)
       endif
#endif          
    endif
  end subroutine computeSingleSGWI

!BOP
! 
! !ROUTINE: LVT_writeMetric_SGWI
! \label(LVT_writeMetric_SGWI)
!
! !INTERFACE:
  subroutine LVT_writeMetric_SGWI(pass,final,vlevels,stats,obs)
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
!    real                    :: model_sgwi(LVT_LIS_rc(1)%ntiles,1,1)
!    real                    :: obs_sgwi(LVT_rc%ngrid,1,1)
!    integer                 :: count_model_sgwi(LVT_LIS_rc(1)%ntiles,1,1)
!    integer                 :: count_obs_sgwi(LVT_rc%ngrid,1,1)

!    count_model_sgwi = 1
!    count_obs_sgwi = 1

    if(pass.eq.LVT_metrics%sgwi%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              do m=1,LVT_rc%nensem
                 do k=1,vlevels
                    if(LVT_metrics%sgwi%timeOpt.eq.1) then 
                       call LVT_writevar_gridded(LVT_metrics%sgwi%ftn_ts, &
                            LVT_sgwi_struc(m)%model_sgwi(:),&
                            stats%vid_ts(LVT_SGWIid,1),1)
                       if(LVT_rc%obssource(2).ne."none") then 
                          call LVT_writevar_gridded(LVT_metrics%sgwi%ftn_ts, &
                               LVT_sgwi_struc(m)%obs_sgwi(:),&
                               stats%vid_ts(LVT_SGWIid,2),1)
                       endif
                    endif
                 enddo
              enddo
#if 0 
             model_sgwi(:,1,1) = LVT_sgwi_struc(m)%model_sgwi(:)
             call LVT_writeSummaryStats(&
                  LVT_metrics%sgwi%ftn_summ,&
                  LVT_metrics%sgwi%short_name,&
                  LVT_LIS_rc(1)%ntiles,&
                  model_sgwi, &
                  count_model_sgwi,&
                  stats%short_name,&
                  LVT_sgwi_struc(m)%model_sgwi_ci)
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                obs_sgwi(:,1,1) = LVT_sgwi_struc(m)%obs_sgwi(:)
                call LVT_writeSummaryStats(&
                     LVT_metrics%sgwi%ftn_summ,&
                     LVT_metrics%sgwi%short_name,&
                     LVT_rc%ngrid,&
                     obs_sgwi, &
                     count_obs_sgwi,&
                     "OBS_"//trim(stats%short_name),&
                     LVT_sgwi_struc(m)%obs_sgwi_ci)
             endif
#endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_SGWI

!BOP
! 
! !ROUTINE: LVT_resetMetric_SGWI
! \label(LVT_resetMetric_SGWI)
!
! !INTERFACE:
  subroutine LVT_resetMetric_SGWI
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine resets the arrays that stores the precipitation
!   values. The reset is done at every stats output writing interval
!   to get the arrays reinitialized for the next set of time series
!   computations.
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
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(stats%selectOpt.eq.1.and.&
            model%short_name.eq.&
            "TWS") then 
          do m=1,LVT_rc%nensem
             if(LVT_metrics%sgwi%timeOpt.eq.1) then 
                LVT_sgwi_struc(m)%model_sgwi  = LVT_rc%udef 
             endif
             if(LVT_rc%obssource(2).ne."none".and.&
                  obs%selectNlevs.ge.1) then 
                LVT_sgwi_struc(m)%obs_sgwi  = LVT_rc%udef 
             endif
          enddo
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo

  end subroutine LVT_resetMetric_SGWI

!BOP
! 
! !ROUTINE: LVT_writerestart_SGWI
! 
! !INTERFACE:
  subroutine LVT_writerestart_SGWI(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for SGWI metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer                          :: k,l,m
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

!    if(LVT_rc%endtime.eq.1) then 
       call LVT_getDataStream1Ptr(model)
       call LVT_getDataStream2Ptr(obs)
       call LVT_getstatsEntryPtr(stats)
       
       do while(associated(model))
          if(model%short_name.eq."TWS") then 
             if(LVT_metrics%sgwi%selectOpt.eq.1) then 
                do m=1,LVT_rc%nensem
                   do l=1,LVT_sgwi_struc(m)%nsize_total
                      call LVT_writevar_restart(ftn,&
                           LVT_sgwi_struc(m)%model_rzsm_final(:,l),tileflag=1)
                   enddo
                   do l=1,LVT_sgwi_struc(m)%nasc
                      call LVT_writevar_restart(ftn,&
                           LVT_sgwi_struc(m)%model_mu(:,l),tileflag=1)
                      call LVT_writevar_restart(ftn,&
                           LVT_sgwi_struc(m)%model_sigma(:,l),tileflag=1)
                   enddo
                   if(LVT_rc%obssource(2).ne."none") then 
                      do l=1,LVT_sgwi_struc(m)%nsize_total
                         call LVT_writevar_restart(ftn,&
                              LVT_sgwi_struc(m)%obs_rzsm_final(:,l))
                      enddo
                      do l=1,LVT_sgwi_struc(m)%nasc
                         call LVT_writevar_restart(ftn,&
                              LVT_sgwi_struc(m)%obs_mu(:,l))
                         call LVT_writevar_restart(ftn,&
                              LVT_sgwi_struc(m)%obs_sigma(:,l))
                      enddo
                   endif
                enddo
             end if
          endif
          model => model%next
          obs   => obs%next
          stats => stats%next

       enddo
!    endif

  end subroutine LVT_writerestart_SGWI

!BOP
! 
! !ROUTINE: LVT_readrestart_SGWI
! 
! !INTERFACE:
  subroutine LVT_readrestart_SGWI(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for SGWI metric computations
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
       if(model%short_name.eq."TWS") then 
          if(LVT_metrics%sgwi%selectOpt.eq.1) then 
             do m=1,LVT_rc%nensem
                do l=1,LVT_sgwi_struc(m)%nsize_total
                   call LVT_readvar_restart(ftn,&
                        LVT_sgwi_struc(m)%model_rzsm_final(:,l),tileflag=1)
                enddo
                do l=1,LVT_sgwi_struc(m)%nasc
                   call LVT_readvar_restart(ftn,&
                        LVT_sgwi_struc(m)%model_mu(:,l),tileflag=1)
                   call LVT_readvar_restart(ftn,&
                        LVT_sgwi_struc(m)%model_sigma(:,l),tileflag=1)
                enddo
                if(LVT_rc%obssource(2).ne."none") then 
                   do l=1,LVT_sgwi_struc(m)%nsize_total
                      call LVT_readvar_restart(ftn,&
                           LVT_sgwi_struc(m)%obs_rzsm_final(:,l))
                   enddo
                   do l=1,LVT_sgwi_struc(m)%nasc
                      call LVT_readvar_restart(ftn,&
                           LVT_sgwi_struc(m)%obs_mu(:,l))
                      call LVT_readvar_restart(ftn,&
                           LVT_sgwi_struc(m)%obs_sigma(:,l))
                   enddo
                endif
             enddo
          end if
       end if
       
       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo

  end subroutine LVT_readrestart_SGWI

!BOP
! 
! !ROUTINE: sgwi_zval
!  \label{sgwi_zval}
! 
! !INTERFACE: 
  subroutine sgwi_zval(t,rzsm, mu, sigma,zval)
! !ARGUMENTS:
    integer, intent(in) :: t
    real,    intent(in) :: rzsm
    real                :: mu, sigma
    real                :: rval
    real                :: zval
! 
! !DESCRPTION:
!   This subroutine computes the SGWI values. The cumulative 
!   probability of a given event is computed from the given 
!   normal distribution. The cumulative probability is then 
!   transformed to the standard normal random variable (SGWI) 
!   with mean zero and variance of one, which the value of 
!   SGWI. 
!
!   \begin{description}
!    \item[t]
!      index of the grid/tile
!    \item[accum\_rzsm]
!      accumulated precipitation value
!    \item[mu]
!      mu parameter of the fitted distribution
!    \item[sigma]
!      sigma parameter of the fitted distribution
!    \item[pzero]
!      probability of zero precipitation
!    \item[zval]
!     transformed standard normal random variable
!   \end{description}  
!EOP

    zval = (rzsm - mu)/sigma

  end subroutine sgwi_zval
 
!BOP
! 
! !ROUTINE: computeNmonths
! \label{computeNmonths}
! 
! !INTERFACE: 
  subroutine computeNmonths(nsize, tavgInterval)
! !ARGUMENTS:     
    integer             :: nsize
    integer, intent(in) :: tavgInterval
!
! !DESCRIPTION: 
!  This subroutine computes the number of months within the 
!  start and stop times based on the time interval
! 
!  The arguments are:
!  \begin{description}
!   \item[nsize] 
!     number of months computed by the routine
!   \item[tavgInterval] 
!     temporal averaging interval (expected to be a multiple of months)
!  \end{description}
!EOP    
    integer             :: yr,mo
    logical             :: togo
    integer             :: mfactor

    mfactor = tavgInterval/2592000

    nsize = 1
    yr = LVT_rc%syr
    mo = LVT_rc%smo
    togo = .true. 

    do while(togo) 
       mo = mo + mfactor

       if(yr.ge.LVT_rc%eyr.and.mo.ge.LVT_rc%emo) then 
          togo = .false. 
       endif
       
       if(mo.gt.12) then 
          mo = mo-12
          yr = yr + 1
       endif

       nsize = nsize + 1
    enddo

  end subroutine computeNmonths
!BOP
! 
! !ROUTINE: getMonthlyTimeIndex
! \label{getMonthlyTimeIndex}
! 
! !INTERFACE: 
  subroutine getMonthlyTimeIndex(nsize, tavgInterval)
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
    integer             :: yr,mo
    logical             :: togo
    integer             :: mfactor

    mfactor = tavgInterval/2592000

    nsize = 1
    yr = LVT_rc%syr
    mo = LVT_rc%smo
    togo = .true. 

    do while(togo) 
       mo = mo + mfactor
       
       if(mo.gt.12) then 
          mo = mo-12
          yr = yr + 1
       endif

       if(yr.ge.LVT_rc%yr.and.mo.ge.LVT_rc%mo) then 
          togo = .false.           
          exit
       endif
       nsize = nsize + 1
    enddo

  end subroutine getMonthlyTimeIndex

!BOP
! 
! !ROUTINE: getSeasonalTimeIndex
!  \label{getSeasonalTimeIndex}
! 
! !INTERFACE: 
  subroutine getSeasonalTimeIndex(tind,tavgInterval)
! !ARGUMENTS:
    integer             :: tind
    integer, intent(in) :: tavgInterval
!
! !DESCRIPTION:
!   
!  This subroutine computes the seasonal time index of the current time,
!  based on the time interval
! 
!  The arguments are:
!  \begin{description}
!   \item[tind] 
!     seasonal time index
!   \item[tavgInterval] 
!     temporal averaging interval (expected to be a multiple of months)
!  \end{description}
!EOP
    if(tavgInterval/2592000.eq.1) then 
       tind = LVT_rc%mo -1
       if(tind.eq.0) then 
          tind = 12
       endif
    elseif(tavgInterval/2592000.eq.3) then
       if(LVT_rc%mo.eq.1.or.LVT_rc%mo.eq.2.or.LVT_rc%mo.eq.3) then
           !DJF
           tind = 1
        elseif(LVT_rc%mo.eq.4.or.LVT_rc%mo.eq.5.or.LVT_rc%mo.eq.6) then
           !MAM
           tind = 2
        elseif(LVT_rc%mo.eq.7.or.LVT_rc%mo.eq.8.or.LVT_rc%mo.eq.9) then
           !JJA
           tind = 3
        elseif(LVT_rc%mo.eq.10.or.LVT_rc%mo.eq.11.or.LVT_rc%mo.eq.12) then
           !SON
           tind = 4
        endif
    elseif(tavgInterval/2592000.eq.6) then
       if(LVT_rc%mo.ge.1.or.LVT_rc%mo.le.6) then 
           tind = 1
        else
           tind = 2
        endif
     else 
        write(LVT_logunit,*) 'getSeasonalTimeIndex needs to be implemented'
        call LVT_endrun()
    endif
  end subroutine getSeasonalTimeIndex

!BOP
! 
! !ROUTINE: normal_fit
! \label{normal_fit}
! 
! !INTERFACE: 
  subroutine normal_fit(nsize, rzsm, mu, sigma)

    implicit none
! !ARGUMENTS:     
    integer,         intent(in) :: nsize
    real,            intent(in) :: rzsm(nsize)
    real                        :: mu
    real                        :: sigma
! 
! !DESCRIPTION: 
!   This subroutine estimates incomplete gamma distribution 
!   parameters. 
! 
!EOP    
    real                        :: sumv
    integer                     :: nsumv
    integer                     :: i

    sumv = 0.0
    nsumv = 0

    do i=1,nsize
       if(rzsm(i).gt.0) then 
          sumv = sumv + rzsm(i)
          nsumv = nsumv + 1
       endif
    enddo

    if(nsumv.ne.0) then 
       mu = sumv/nsumv
    endif  

    sumv = 0 
    nsumv = 0 

    do i=1,nsize
       if(rzsm(i).gt.0) then 
          sumv = sumv + (rzsm(i)-mu)**2
          nsumv = nsumv + 1
       endif
    enddo
    
    if(nsumv.ne.0) then 
       sigma = sqrt(sumv/nsumv)
    else
       mu = LVT_rc%udef
       sigma = LVT_rc%udef
    endif

  end subroutine normal_fit

end module LVT_SGWIMod
