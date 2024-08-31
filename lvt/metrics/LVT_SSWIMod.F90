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
! !MODULE: LVT_SSWIMod
! \label(LVT_SSWIMod)
!
! !INTERFACE:
module LVT_SSWIMod
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
!   This module handles the Standardized soil water index (SSWI)
!   metric based on the root zone soil moisture outputs in LIS. 
!   
!   Reference: McKee, T. B., N. J. Doesken, and J. Kleist, 1993:
!   The relationship of drought frequency and duration to time scales. 
!   Preprints, Eighth Conf. on Applied Climatology, Anaheim, CA, 
!   Amer. Meteor. Soc., 179–184.
! 
! !NOTES: 
!   The standard practice is to compute SSWI on monthly scale. Since the 
!   code needs to store the soil moisture values in time to compute SSWI, 
!   currently LVT only supports SSWI computations when the time 
!   averaging interval is set to a multiple of a month (2592000 seconds). 
!
!   Only monthly SSWI computations have been tested so far
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
  public :: LVT_initSSWI
  public :: LVT_diagnoseSSWI
  public :: LVT_computeSSWI
  public :: LVT_writeMetric_SSWI
  public :: LVT_resetMetric_SSWI
  public :: LVT_writerestart_SSWI
  public :: LVT_readrestart_SSWI

!EOP

  type, public :: sswidec
     real,    allocatable :: model_rzsm_final(:,:)
     real,    allocatable :: obs_rzsm_final(:,:)

     real,    allocatable :: model_mu(:,:)
     real,    allocatable :: model_sigma(:,:)
     real,    allocatable :: model_sswi(:)

     real,    allocatable :: obs_mu(:,:)
     real,    allocatable :: obs_sigma(:,:)
     real,    allocatable :: obs_sswi(:)
     real,    allocatable :: model_current(:,:)
     real,    allocatable :: obs_current(:,:)

     real             :: model_sswi_ci(1,1)
     real             :: obs_sswi_ci(1,1)

     integer          :: nsize_total
     integer          :: nsize_season
     integer          :: nasc
     integer          :: tscale
  end type sswidec

  type(sswidec), allocatable :: LVT_sswi_struc(:)
  private
  integer, parameter :: maxiter = 100
  real,    parameter :: epsilon = 3.0e-7


contains
!BOP
!
! !ROUTINE: LVT_initSSWI
! \label{LVT_initSSWI}
! 
! !INTERFACE: 
  subroutine LVT_initSSWI(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This routine initializes the data structures required
!  for the SSWI computations. 
! 
!EOP

    type(ESMF_Time)         :: startTime
    type(ESMF_Time)         :: stopTime
    type(ESMF_TimeInterval) :: timeStep
    character*10            :: time
    integer                 :: rc
    integer                 :: m

    allocate(LVT_sswi_struc(LVT_rc%nensem))

    if(LVT_metrics%sswi%selectOpt.eq.1.or.&
         LVT_metrics%sswi%timeOpt.eq.1) then 
       call ESMF_ConfigGetAttribute(LVT_config,time,&
            label="SSWI timescale of computation:",rc=rc)
       if(rc.ne.0) then 
          write(LVT_logunit,*) "[ERR] 'SSWI timescale of computation:' not defined"
          write(LVT_logunit,*) "[ERR] Options are 1mo, 3mo, 6mo, 9mo, 12mo or 24mo"
          write(LVT_logunit,*) "[ERR] Corresponding to 1-month SSWI, 3-month SSWI and so on"
          call LVT_endrun()
       endif
    endif

    if(LVT_metrics%sswi%selectOpt.eq.1.or.&
         LVT_metrics%sswi%timeOpt.eq.1) then 
       
       do m=1,LVT_rc%nensem
          call LVT_parseTimeString(time,LVT_sswi_struc(m)%tscale)
          LVT_sswi_struc(m)%tscale = LVT_sswi_struc(m)%tscale/2592000          
          
          if(LVT_rc%tavgInterval.ne.2592000) then 
             write(LVT_logunit,*)'[ERR] SSWI values are computed at a monthly timescale'
             write(LVT_logunit,*)'[ERR] Please set the time averaging interval to a '
             write(LVT_logunit,*)'  multiple of a month (2592000 seconds).'
             write(LVT_logunit,*)'Program stopping ...'
             call LVT_endrun()
          endif
          
          LVT_sswi_struc(m)%nasc = 12
          
          call computeNmonths(LVT_sswi_struc(m)%nsize_total, LVT_rc%tavgInterval)
          LVT_sswi_struc(m)%nsize_season = &
               LVT_sswi_struc(m)%nsize_total/LVT_sswi_struc(m)%nasc
          
          if(metric%selectOpt.eq.1) then 
             allocate(LVT_sswi_struc(m)%model_rzsm_final(LVT_rc%ngrid,&
                  LVT_sswi_struc(m)%nsize_total))
             LVT_sswi_struc(m)%model_rzsm_final = 0.0
             
             allocate(LVT_sswi_struc(m)%model_mu(LVT_rc%ngrid,&
                  LVT_sswi_struc(m)%nasc))
             allocate(LVT_sswi_struc(m)%model_sigma(LVT_rc%ngrid,&
                  LVT_sswi_struc(m)%nasc))
             allocate(LVT_sswi_struc(m)%model_sswi(LVT_rc%ngrid))
             LVT_sswi_struc(m)%model_sswi = LVT_rc%udef

             allocate(LVT_sswi_struc(m)%model_current(LVT_rc%ngrid,&
                  LVT_sswi_struc(m)%tscale))
             LVT_sswi_struc(m)%model_current = LVT_rc%udef
             
             if(LVT_rc%obssource(2).ne."none") then 
                if(selectNlevs(2).ge.1) then 
                   allocate(LVT_sswi_struc(m)%obs_rzsm_final(LVT_rc%ngrid,&
                        LVT_sswi_struc(m)%nsize_total))
                   LVT_sswi_struc(m)%obs_rzsm_final = 0.0
                   
                   allocate(LVT_sswi_struc(m)%obs_mu(LVT_rc%ngrid,&
                        LVT_sswi_struc(m)%nasc))
                   allocate(LVT_sswi_struc(m)%obs_sigma(LVT_rc%ngrid,&
                        LVT_sswi_struc(m)%nasc))
                   allocate(LVT_sswi_struc(m)%obs_sswi(LVT_rc%ngrid))
                   LVT_sswi_struc(m)%obs_sswi = LVT_rc%udef
                  
                   allocate(LVT_sswi_struc(m)%obs_current(LVT_rc%ngrid,1))
                   LVT_sswi_struc(m)%obs_current = LVT_rc%udef
                   
                endif
             endif
          endif
       enddo
    endif
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

  end subroutine LVT_initSSWI

!BOP
! 
! !ROUTINE: LVT_diagnoseSSWI
! \label{LVT_diagnoseSSWI}
!
! !INTERFACE: 
  subroutine LVT_diagnoseSSWI(pass)
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
!   calculating the sswi of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelSSWI](\ref{diagnoseSingleModelSSWI})
!     updates the sswi computation for a single variable. This routine
!     stores the precip values to be used later for computing SSWI, 
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
       if(LVT_metrics%sswi%selectOpt.eq.1.or.&
            LVT_metrics%sswi%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))   
             
             call diagnoseSingleSSWI(model, obs, stats,&
                  LVT_metrics%sswi)
             
             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseSSWI


!BOP
! 
! !ROUTINE: diagnoseSingleSSWI
! \label{diagnoseSingleSSWI}
!
! !INTERFACE: 
  subroutine diagnoseSingleSSWI(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the sswi computation of the 
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
    integer                 :: t,k,tind,m
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_TimeInterval) :: ts
    integer                 :: tindex_f,tindex_t
    integer                 :: rc,m_k,o_k

    if(stats%selectOpt.eq.1.and.&
         ((model%short_name.eq."RootMoist").or.&
         (model%short_name.eq."TWS").or.&
         (model%short_name.eq."SWE")).and.&
         model%selectNlevs.ge.1) then        

       call getMonthlyTimeIndex(tindex_f, LVT_rc%tavgInterval)
       
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(model%count(t,m,m_k).gt.0) then
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         LVT_sswi_struc(m)%model_rzsm_final(t,tindex_f) = &
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
                            LVT_sswi_struc(m)%obs_rzsm_final(t,tindex_f) = &
                                 obs%value(t,m,o_k)
                         endif
                         
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
  end subroutine diagnoseSingleSSWI

!BOP
! 
! !ROUTINE: LVT_computeSSWI
! \label{LVT_computeSSWI}
!
! !INTERFACE: 
  subroutine LVT_computeSSWI(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the SSWI values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelSSWI](\ref{computeSingleModelSSWI})
!     computes the SSWI values for a single variable
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
          if(LVT_metrics%sswi%selectOpt.eq.1.or.LVT_metrics%sswi%timeOpt.eq.1) then 
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
          
             do while(associated(model))             
                call computeSingleSSWIparams(alarm,&
                     model,obs,stats,LVT_metrics%sswi)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
          endif
       elseif(pass.eq.2) then        
          if(LVT_metrics%sswi%selectOpt.eq.1.or.LVT_metrics%sswi%timeOpt.eq.1) then 
             if(alarm) then 
                if(LVT_metrics%sswi%timeOpt.eq.1.and.&
                     LVT_metrics%sswi%extractTS.eq.1) then 

                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%sswi%ftn_ts_loc(i,m),200,advance='no') &
                                 LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                                 LVT_rc%hr,'',LVT_rc%mn, '' 
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sswi%ftn_ts_loc(i,1),200,advance='no') &
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
                call computeSingleSSWI(alarm,&
                     model,obs,stats,LVT_metrics%sswi)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
             
             if(alarm) then 
                if(LVT_metrics%sswi%timeOpt.eq.1.and.&
                     LVT_metrics%sswi%extractTS.eq.1) then 
                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%sswi%ftn_ts_loc(i,m),fmt='(a1)') ''
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sswi%ftn_ts_loc(i,1),fmt='(a1)') ''
                      enddo
                   endif

                endif
             endif
          endif
       endif
    elseif(LVT_rc%startmode.eq."restart") then 
       if(LVT_metrics%sswi%selectOpt.eq.1.or.LVT_metrics%sswi%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%sswi%timeOpt.eq.1.and.&
                  LVT_metrics%sswi%extractTS.eq.1) then 

                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%sswi%ftn_ts_loc(i,m),200,advance='no') &
                                 LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                                 LVT_rc%hr,'',LVT_rc%mn, '' 
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sswi%ftn_ts_loc(i,1),200,advance='no') &
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
             call computeSingleSSWI(alarm,&
                  model,obs,stats,LVT_metrics%sswi)
             
             model => model%next
             obs => obs%next
             stats => stats%next
             
          enddo
          
          if(alarm) then 
             if(LVT_metrics%sswi%timeOpt.eq.1.and.&
                  LVT_metrics%sswi%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%sswi%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%sswi%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif
                
             endif
          endif
       endif
    endif

  end subroutine LVT_computeSSWI
  

!BOP
! 
! !ROUTINE: computeSingleSSWIparams
! \label{computeSingleSSWIparams}
!
! !INTERFACE: 
  subroutine computeSingleSSWIparams(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the SSWI values for each grid cell. 
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
            ((model%short_name.eq."RootMoist").or.&
            (model%short_name.eq."TWS").or.&
            (model%short_name.eq."SWE")).and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_sswi_struc(m)%nasc
                      call sswi_normal(t,l,m,LVT_sswi_struc(m)%nsize_total,&
                           LVT_sswi_struc(m)%nsize_season,&
                           LVT_sswi_struc(m)%model_rzsm_final(t,:),&
                           LVT_sswi_struc(m)%model_mu(t,l),&
                           LVT_sswi_struc(m)%model_sigma(t,l))
                   enddo
                enddo
             enddo
          enddo
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then 
                         do l=1,LVT_sswi_struc(m)%nasc
                            call sswi_normal(t,l,m,LVT_sswi_struc(m)%nsize_total,&
                                 LVT_sswi_struc(m)%nsize_season,&
                                 LVT_sswi_struc(m)%obs_rzsm_final(t,:),&
                                 LVT_sswi_struc(m)%obs_mu(t,l),&
                                 LVT_sswi_struc(m)%obs_sigma(t,l))
                         enddo
                      endif
                   endif
                enddo
             enddo
          enddo
       endif
    endif
  end subroutine computeSingleSSWIparams

!BOP
! !ROUTINE:sswi_normal
! \label{sswi_normal}
! 
! !INTERFACE: 
  subroutine sswi_normal(t,l,kk,nsize,nsize_season,rzsm, mu, sigma)

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: t
    integer, intent(in) :: l
    integer, intent(in) :: kk
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
    integer             :: i,k,m,t1_ind
    real                :: rzsm_season(nsize_season)
    real                :: tmprzsm
    integer             :: rzsm_count

    rzsm_season = 0.0
    k = 0 
    do i=l,nsize,LVT_sswi_struc(kk)%nasc
       k = k+1

       tmprzsm = 0.0
       rzsm_count = 0 
       do m=1,LVT_sswi_struc(kk)%tscale
          t1_ind = (i-LVT_sswi_struc(kk)%tscale+m)
          if(t1_ind.gt.0) then 
             tmprzsm = tmprzsm + rzsm(t1_ind)
             rzsm_count = rzsm_count + 1
          endif
       enddo
       rzsm_season(k) = tmprzsm/rzsm_count
    enddo

    call normal_fit(nsize_season, rzsm_season, &
         mu, sigma)

  end subroutine sswi_normal

!BOP
! 
! !ROUTINE: computeSingleSSWI
! \label{computeSingleSSWI}
!
! !INTERFACE: 
  subroutine computeSingleSSWI(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the SSWI values for each grid cell. 
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

    integer  :: t,l,k,i,m,ii
    integer  :: tindex_f, tind
    real     :: model_sswi(LVT_rc%ngrid,1,1)
    real     :: obs_sswi(LVT_rc%ngrid,1,1)
    integer  :: count_model_sswi(LVT_rc%ngrid,1,1)
    integer  :: count_obs_sswi(LVT_rc%ngrid,1,1)
    real                :: tmpval
    integer             :: val_count
    real,    allocatable :: tavg_model_value_ts(:,:,:)
    real,    allocatable :: tavg_obs_value_ts(:,:,:)
    integer, allocatable :: count_tavg_model_value_ts(:,:,:)

    if(stats%selectOpt.eq.1.and.&
         ((model%short_name.eq."RootMoist").or.&
         (model%short_name.eq."TWS").or.&
         (model%short_name.eq."SWE")).and.&
         model%selectNlevs.ge.1) then 

       call getMonthlyTimeIndex(tindex_f, LVT_rc%tavgInterval)
       call getSeasonalTimeIndex(tind,LVT_rc%tavgInterval)

       do ii=1,LVT_rc%nensem
          do m=LVT_sswi_struc(ii)%tscale,1,-1
             if(m.ne.1) then 
                LVT_sswi_struc(ii)%model_current(:,m) = &
                     LVT_sswi_struc(ii)%model_current(:,m-1)
             else
                LVT_sswi_struc(ii)%model_current(:,m) = & 
                     model%value(:,m,1)
             endif
          enddo
       enddo

       do t=1,LVT_rc%ngrid
          do ii=1,LVT_rc%nensem
             tmpval = 0.0
             val_count = 0 
             
             do m=1,LVT_sswi_struc(ii)%tscale
                if(LVT_sswi_struc(ii)%model_current(t,m).ne.LVT_rc%udef) then 
                   tmpval = tmpval + LVT_sswi_struc(ii)%model_current(t,m)
                   val_count = val_count + 1
                endif
             enddo
             tmpval = tmpval/val_count


             call sswi_zval(t,&
                  tmpval,&
                  LVT_sswi_struc(ii)%model_mu(t,tind), &
                  LVT_sswi_struc(ii)%model_sigma(t,tind), &
                  LVT_sswi_struc(ii)%model_sswi(t))
             
          enddo
       enddo

       if(LVT_rc%obssource(2).ne."none") then  
          if(obs%selectNlevs.ge.1) then 
             do ii=1,LVT_rc%nensem
                do m=LVT_sswi_struc(ii)%tscale,1,-1
                   if(m.ne.1) then 
                      LVT_sswi_struc(ii)%obs_current(:,m) = &
                           LVT_sswi_struc(ii)%obs_current(:,m-1)
                   else
                      LVT_sswi_struc(ii)%obs_current(:,m) = & 
                           obs%value(:,m,1)
                   endif
                enddo
             enddo

             do t=1,LVT_rc%ngrid
                do ii=1,LVT_rc%nensem
                   tmpval = 0.0
                   val_count = 0 
                   
                   do m=1,LVT_sswi_struc(ii)%tscale
                      if(LVT_sswi_struc(ii)%obs_current(t,m).ne.LVT_rc%udef) then 
                         tmpval = tmpval + LVT_sswi_struc(ii)%obs_current(t,m)
                         val_count = val_count + 1
                      endif
                   enddo
                   if(val_count.gt.0) then 
                      tmpval = tmpval/val_count
                      
                      call sswi_zval(t,&
                           tmpval,&
                           LVT_sswi_struc(ii)%obs_mu(t,tind), &
                           LVT_sswi_struc(ii)%obs_sigma(t,tind), &
                           LVT_sswi_struc(ii)%obs_sswi(t))
                   else
                      LVT_sswi_struc(ii)%obs_sswi(t) = LVT_rc%udef
                   endif
                enddo
             enddo
          endif
       endif

       count_model_sswi = 1
       count_obs_sswi = 1       
       if(metric%extractTS.eq.1) then
          if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
             do ii=1,LVT_rc%nensem
                model_sswi(:,1,1) = LVT_sswi_struc(ii)%model_sswi(:)
                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   obs_sswi(:,1,1) = LVT_sswi_struc(ii)%model_sswi(:)
                   call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                        model,&
                        LVT_rc%ngrid,&
                        model_sswi,&
                        count_model_sswi,&
                        LVT_rc%ngrid,&
                        obs_sswi,&
                        count_obs_sswi)
                   
                else
                   call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                        model,&
                        LVT_rc%ngrid,&
                        model_sswi,&
                        count_model_sswi)
                endif
             enddo
          else

             allocate(tavg_model_value_ts(LVT_rc%ngrid,1,1))
             allocate(count_tavg_model_value_ts(LVT_rc%ngrid,1,1))
             tavg_model_value_ts = 0.0
             count_tavg_model_value_ts = 0.0
             
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                
                allocate(tavg_obs_value_ts(LVT_rc%ngrid,1,1))
                tavg_obs_value_ts = 0.0
                
                do ii=1,LVT_rc%nensem
                   do t=1,LVT_rc%ngrid
                      if(LVT_sswi_struc(ii)%model_sswi(t).ne.LVT_rc%udef) then 
                         tavg_model_value_ts(t,1,1) = &
                              tavg_model_value_ts(t,1,1) + &
                              LVT_sswi_struc(ii)%model_sswi(t)
                         tavg_obs_value_ts(t,1,1) = &
                              tavg_obs_value_ts(t,1,1) + &
                              LVT_sswi_struc(ii)%obs_sswi(t)
                         count_tavg_model_value_ts(t,1,1) = & 
                              count_tavg_model_value_ts(t,1,1) + 1
                              
                      endif
                   enddo
                enddo
                
                do t=1,LVT_rc%ngrid
                   if(count_tavg_model_value_ts(t,1,1).gt.0) then 
                      tavg_model_value_ts(t,1,1) = &
                           tavg_model_value_ts(t,1,1)/&
                           count_tavg_model_value_ts(t,1,1)
                      tavg_obs_value_ts(t,1,1) = &
                           tavg_obs_value_ts(t,1,1)/&
                           count_tavg_model_value_ts(t,1,1)
                   else
                      tavg_model_value_ts(t,1,1) = LVT_rc%udef
                      tavg_obs_value_ts(t,1,1) = LVT_rc%udef

                   endif
                enddo
                         
                call LVT_writeTSinfo(metric%ftn_ts_loc(:,1),&
                     model,&
                     LVT_rc%ngrid,&
                     tavg_model_value_ts,&
                     count_model_sswi,&
                     LVT_rc%ngrid,&
                     tavg_obs_value_ts,&
                     count_obs_sswi)
             else
                do ii=1,LVT_rc%nensem
                   do t=1,LVT_rc%ngrid
                      tavg_model_value_ts(t,1,1) = &
                           tavg_model_value_ts(t,1,1) + &
                           LVT_sswi_struc(ii)%model_sswi(t)
                      count_tavg_model_value_ts(t,1,1) = & 
                           count_tavg_model_value_ts(t,1,1) + 1
                   enddo
                enddo
                do t=1,LVT_rc%ngrid
                   if(count_tavg_model_value_ts(t,1,1).gt.0) then 
                      tavg_model_value_ts(t,1,1) = &
                           tavg_model_value_ts(t,1,1)/&
                           count_tavg_model_value_ts(t,1,1)
                   else
                      tavg_model_value_ts(t,1,1) = LVT_rc%udef
                   endif
                enddo

                call LVT_writeTSinfo(metric%ftn_ts_loc(:,1),&
                     model,&
                     LVT_rc%ngrid,&
                     tavg_model_value_ts,&
                     count_model_sswi)
             endif
             deallocate(tavg_model_value_ts)
             deallocate(count_tavg_model_value_ts)
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                deallocate(tavg_obs_value_ts)
             endif


          end if
       end if

       do ii=1,LVT_rc%nensem
          
          call LVT_computeCI(LVT_sswi_struc(ii)%model_sswi(:),LVT_rc%ngrid,&
               LVT_rc%pval_CI, LVT_sswi_struc(ii)%model_sswi_ci(1,1))
          
          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
             call LVT_computeCI(LVT_sswi_struc(ii)%obs_sswi(:),LVT_rc%ngrid,&
                  LVT_rc%pval_CI,LVT_sswi_struc(ii)%obs_sswi_ci(1,1))
          endif
       enddo

!       call LVT_writeDataBasedStrat(model,obs,stats,metric,&
!            LVT_rc%ngrid, model_sswi)
!
!       if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
!          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
!               LVT_rc%ngrid, model_sswi, &
!               LVT_rc%ngrid, obs_sswi)
          
    endif

  end subroutine computeSingleSSWI

!BOP
! 
! !ROUTINE: LVT_writeMetric_SSWI
! \label(LVT_writeMetric_SSWI)
!
! !INTERFACE:
  subroutine LVT_writeMetric_SSWI(pass,final,vlevels,stats,obs)
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
!    real                    :: model_sswi(LVT_rc%ngrid,1,1)
!    real                    :: obs_sswi(LVT_rc%ngrid,1,1)
!    integer                 :: count_model_sswi(LVT_rc%ngrid,1,1)
!    integer                 :: count_obs_sswi(LVT_rc%ngrid,1,1)

!    count_model_sswi = 1
!    count_obs_sswi = 1

    if(pass.eq.LVT_metrics%sswi%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              do m=1,LVT_rc%nensem
                 do k=1,vlevels
                    if(LVT_metrics%sswi%timeOpt.eq.1) then 
                       call LVT_writevar_gridded(LVT_metrics%sswi%ftn_ts, &
                            LVT_sswi_struc(m)%model_sswi(:),&
                            stats%vid_ts(LVT_SSWIid,1),1)
                       if(LVT_rc%obssource(2).ne."none") then 
                          call LVT_writevar_gridded(LVT_metrics%sswi%ftn_ts, &
                               LVT_sswi_struc(m)%obs_sswi(:),&
                               stats%vid_ts(LVT_SSWIid,2),1)
                       endif
                    endif
                 enddo
              enddo
#if 0 
             model_sswi(:,1,1) = LVT_sswi_struc(m)%model_sswi(:)
             call LVT_writeSummaryStats(&
                  LVT_metrics%sswi%ftn_summ,&
                  LVT_metrics%sswi%short_name,&
                  LVT_rc%ngrid,&
                  model_sswi, &
                  count_model_sswi,&
                  stats%short_name,&
                  LVT_sswi_struc(m)%model_sswi_ci)
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                obs_sswi(:,1,1) = LVT_sswi_struc(m)%obs_sswi(:)
                call LVT_writeSummaryStats(&
                     LVT_metrics%sswi%ftn_summ,&
                     LVT_metrics%sswi%short_name,&
                     LVT_rc%ngrid,&
                     obs_sswi, &
                     count_obs_sswi,&
                     "OBS_"//trim(stats%short_name),&
                     LVT_sswi_struc(m)%obs_sswi_ci)
             endif
#endif
          endif
       endif
    endif
  end subroutine LVT_writeMetric_SSWI

!BOP
! 
! !ROUTINE: LVT_resetMetric_SSWI
! \label(LVT_resetMetric_SSWI)
!
! !INTERFACE:
  subroutine LVT_resetMetric_SSWI
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
            ((model%short_name.eq."RootMoist").or.&
            (model%short_name.eq."TWS").or.&
            (model%short_name.eq."SWE"))) then 
          do m=1,LVT_rc%nensem
             if(LVT_metrics%sswi%timeOpt.eq.1) then 
                LVT_sswi_struc(m)%model_sswi  = LVT_rc%udef 
             endif
             if(LVT_rc%obssource(2).ne."none".and.&
                  obs%selectNlevs.ge.1) then 
                LVT_sswi_struc(m)%obs_sswi  = LVT_rc%udef 
             endif
          enddo
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo

  end subroutine LVT_resetMetric_SSWI

!BOP
! 
! !ROUTINE: LVT_writerestart_SSWI
! 
! !INTERFACE:
  subroutine LVT_writerestart_SSWI(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for SSWI metric computations
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

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if((model%short_name.eq."RootMoist").or.&
            (model%short_name.eq."TWS").or.&
            (model%short_name.eq."SWE")) then 
          if(LVT_metrics%sswi%selectOpt.eq.1) then 
             do m=1,LVT_rc%nensem
                do l=1,LVT_sswi_struc(m)%nasc
                   call LVT_writevar_restart(ftn,&
                        LVT_sswi_struc(m)%model_mu(:,l))
                   call LVT_writevar_restart(ftn,&
                        LVT_sswi_struc(m)%model_sigma(:,l))
                enddo
                if(LVT_rc%obssource(2).ne."none") then 
                   do l=1,LVT_sswi_struc(m)%nasc
                      call LVT_writevar_restart(ftn,&
                           LVT_sswi_struc(m)%obs_mu(:,l))
                      call LVT_writevar_restart(ftn,&
                           LVT_sswi_struc(m)%obs_sigma(:,l))
                   enddo
                endif
             enddo
          end if
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
       
    enddo

  end subroutine LVT_writerestart_SSWI

!BOP
! 
! !ROUTINE: LVT_readrestart_SSWI
! 
! !INTERFACE:
  subroutine LVT_readrestart_SSWI(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for SSWI metric computations
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
       if((model%short_name.eq."RootMoist").or.&
            (model%short_name.eq."TWS").or.&
            (model%short_name.eq."SWE")) then 
          if(LVT_metrics%sswi%selectOpt.eq.1) then 
             do m=1,LVT_rc%nensem
                do l=1,LVT_sswi_struc(m)%nasc
                   call LVT_readvar_restart(ftn,&
                        LVT_sswi_struc(m)%model_mu(:,l))
                   call LVT_readvar_restart(ftn,&
                        LVT_sswi_struc(m)%model_sigma(:,l))
                enddo
                if(LVT_rc%obssource(2).ne."none") then 
                   do l=1,LVT_sswi_struc(m)%nasc
                      call LVT_readvar_restart(ftn,&
                           LVT_sswi_struc(m)%obs_mu(:,l))
                      call LVT_readvar_restart(ftn,&
                           LVT_sswi_struc(m)%obs_sigma(:,l))
                   enddo
                endif
             enddo
          end if
       end if
       
       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo

  end subroutine LVT_readrestart_SSWI

!BOP
! 
! !ROUTINE: sswi_zval
!  \label{sswi_zval}
! 
! !INTERFACE: 
  subroutine sswi_zval(t,rzsm, mu, sigma,zval)
! !ARGUMENTS:
    integer, intent(in) :: t
    real,    intent(in) :: rzsm
    real                :: mu, sigma
    real                :: rval
    real                :: zval
! 
! !DESCRPTION:
!   This subroutine computes the SSWI values. The cumulative 
!   probability of a given event is computed from the given 
!   normal distribution. The cumulative probability is then 
!   transformed to the standard normal random variable (SSWI) 
!   with mean zero and variance of one, which the value of 
!   SSWI. 
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
    
    if(sigma.ne.0) then 
       zval = (rzsm - mu)/sigma
    endif

  end subroutine sswi_zval
 
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
       tind = LVT_rc%mo - 1
       if(tind.lt.1) then 
          tind = 12
       endif
    else
       write(LVT_logunit,*) '[ERR] getSeasonalTimeIndex needs to be implemented'
       write(LVT_logunit,*) 'Program stopping ...'
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

end module LVT_SSWIMod
