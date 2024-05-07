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
! !MODULE: LVT_SPIMod
! \label(LVT_SPIMod)
!
! !INTERFACE:
module LVT_SPIMod
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
!   This module handles the Standard Precipitation Index (SPI)
!   metric based on the precipitation outputs in LIS. 
!   
!   Reference: McKee, T. B., N. J. Doesken, and J. Kleist, 1993:
!   The relationship of drought frequency and duration to time scales. 
!   Preprints, Eighth Conf. on Applied Climatology, Anaheim, CA, 
!   Amer. Meteor. Soc., 179–184.
! 
! !NOTES: 
!   The standard practice is to compute SPI on monthly scale. Since the 
!   code needs to store the precip values in time to compute SPI, 
!   currently LVT only supports SPI computations when the time 
!   averaging interval is set to a multiple of a month (2592000 seconds). 
!
!   Interpretation of SPI at different timescales: the 6-month SPI represents
!   how the most recent 6 months of precipitation
!   compares to that same 6-month period historically.
!   This allows for a cumulative depiction of precipitation
!   deficits or excesses. The SPIs at longer timescales provide estimates
!   of precip patterns at longer timescales. 
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
  public :: LVT_initSPI
  public :: LVT_diagnoseSPI
  public :: LVT_computeSPI
  public :: LVT_writeMetric_SPI
  public :: LVT_resetMetric_SPI
  public :: LVT_writerestart_SPI
  public :: LVT_readrestart_SPI

!EOP

  type, public :: spidec
     real,    allocatable :: model_pcp_final(:,:)
     real,    allocatable :: obs_pcp_final(:,:)

     real,    allocatable :: model_beta(:,:)
     real,    allocatable :: model_gamma(:,:)
     real,    allocatable :: model_pzero(:,:)
     real,    allocatable :: model_spi(:)

     real,    allocatable :: obs_beta(:,:)
     real,    allocatable :: obs_gamma(:,:)
     real,    allocatable :: obs_pzero(:,:)
     real,    allocatable :: obs_spi(:)
     real,    allocatable :: model_current(:,:)
     real,    allocatable :: obs_current(:,:)

     real             :: model_spi_ci(1,1)
     real             :: obs_spi_ci(1,1)

     integer          :: nsize_total
     integer          :: nsize_season
     integer          :: nasc
     integer          :: tscale
  end type spidec

  type(spidec), allocatable :: LVT_spi_struc(:)
  private
  integer, parameter :: maxiter = 100
  real,    parameter :: epsilon = 3.0e-7


contains
!BOP
!
! !ROUTINE: LVT_initSPI
! \label{LVT_initSPI}
! 
! !INTERFACE: 
  subroutine LVT_initSPI(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This routine initializes the data structures required
!  for the SPI computations. 
! 
!EOP

    type(ESMF_Time)         :: startTime
    type(ESMF_Time)         :: stopTime
    type(ESMF_TimeInterval) :: timeStep
    character*10            :: time
    integer                 :: rc

    integer                 :: m 

    allocate(LVT_spi_struc(LVT_rc%nensem))
    
    if(LVT_metrics%spi%selectOpt.eq.1.or.&
         LVT_metrics%spi%timeOpt.eq.1) then 
       call ESMF_ConfigGetAttribute(LVT_config,time,&
            label="SPI timescale of computation:",rc=rc)
       if(rc.ne.0) then 
          write(LVT_logunit,*) "[ERR] 'SPI timescale of computation:' not defined"
          write(LVT_logunit,*) "[ERR] Options are 1mo, 3mo, 6mo, 9mo, 12mo or 24mo"
          write(LVT_logunit,*) "[ERR] Corresponding to 1-month SPI, 3-month SPI and so on"
          write(LVT_logunit,*) ' Programg stopping ... '
          call LVT_endrun()
       endif
    endif
    if(LVT_metrics%spi%selectOpt.eq.1.or.&
         LVT_metrics%spi%timeOpt.eq.1) then 

       do m=1,LVT_rc%nensem
          
          call LVT_parseTimeString(time,LVT_spi_struc(m)%tscale)
          LVT_spi_struc(m)%tscale = LVT_spi_struc(m)%tscale/2592000

          if(LVT_rc%tavgInterval.ne.2592000) then 
             write(LVT_logunit,*) '[ERR] '
             write(LVT_logunit,*) '[ERR] SPI values are computed at a monthly timescale'
             write(LVT_logunit,*) '[ERR] Please set the time averaging interval to a month '
             write(LVT_logunit,*) '[ERR] program stopping...'
             write(LVT_logunit,*) '[ERR] '
             call LVT_endrun()
          endif
          LVT_spi_struc(m)%nasc = 12
          
          call computeNmonths(LVT_spi_struc(m)%nsize_total, LVT_rc%tavgInterval)
          LVT_spi_struc(m)%nsize_season = &
               LVT_spi_struc(m)%nsize_total/LVT_spi_struc(m)%nasc
          
          if(metric%selectOpt.eq.1) then 
             allocate(LVT_spi_struc(m)%model_pcp_final(LVT_rc%ngrid,&
                  LVT_spi_struc(m)%nsize_total))
             LVT_spi_struc(m)%model_pcp_final = 0.0
             
             allocate(LVT_spi_struc(m)%model_beta(LVT_rc%ngrid,&
                  LVT_spi_struc(m)%nasc))
             allocate(LVT_spi_struc(m)%model_gamma(LVT_rc%ngrid,&
                  LVT_spi_struc(m)%nasc))
             allocate(LVT_spi_struc(m)%model_pzero(LVT_rc%ngrid,&
                  LVT_spi_struc(m)%nasc))
             allocate(LVT_spi_struc(m)%model_spi(LVT_rc%ngrid))
             LVT_spi_struc(m)%model_spi = LVT_rc%udef
             
             allocate(LVT_spi_struc(m)%model_current(LVT_rc%ngrid,&
                  LVT_spi_struc(m)%tscale))
             LVT_spi_struc(m)%model_current = LVT_rc%udef

             if(LVT_rc%obssource(2).ne."none") then 
                if(selectNlevs(2).ge.1) then 
                   allocate(LVT_spi_struc(m)%obs_pcp_final(LVT_rc%ngrid,&
                        LVT_spi_struc(m)%nsize_total))
                   LVT_spi_struc(m)%obs_pcp_final = 0.0
                   
                   allocate(LVT_spi_struc(m)%obs_beta(LVT_rc%ngrid,&
                        LVT_spi_struc(m)%nasc))
                   allocate(LVT_spi_struc(m)%obs_gamma(LVT_rc%ngrid,&
                        LVT_spi_struc(m)%nasc))
                   allocate(LVT_spi_struc(m)%obs_pzero(LVT_rc%ngrid,&
                        LVT_spi_struc(m)%nasc))
                   allocate(LVT_spi_struc(m)%obs_spi(LVT_rc%ngrid))
                   LVT_spi_struc(m)%obs_spi = LVT_rc%udef
                   
                   allocate(LVT_spi_struc(m)%obs_current(LVT_rc%ngrid,1))
                   LVT_spi_struc(m)%obs_current = LVT_rc%udef

                endif
             endif
             
          endif
       end do
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------
       if(LVT_rc%startmode.eq."coldstart") then           
          metric%npass = 2   
       else
          metric%npass = 1
       endif
       
       if(LVT_rc%obssource(2).ne."none") then 
          metric%obsData = .true. 
       else
          metric%obsData = .false. 
       endif

       metric%stdevFlag = .false. 

    endif

  end subroutine LVT_initSPI

!BOP
! 
! !ROUTINE: LVT_diagnoseSPI
! \label{LVT_diagnoseSPI}
!
! !INTERFACE: 
  subroutine LVT_diagnoseSPI(pass)
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
!   calculating the spi of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelSPI](\ref{diagnoseSingleModelSPI})
!     updates the spi computation for a single variable. This routine
!     stores the precip values to be used later for computing SPI, 
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
       if(LVT_metrics%spi%selectOpt.eq.1.or.&
            LVT_metrics%spi%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))   
             
             call diagnoseSingleSPI(model, obs, stats,&
                  LVT_metrics%spi)
             
             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif
  end subroutine LVT_diagnoseSPI


!BOP
! 
! !ROUTINE: diagnoseSingleSPI
! \label{diagnoseSingleSPI}
!
! !INTERFACE: 
  subroutine diagnoseSingleSPI(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the spi computation of the 
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
    integer                 :: t,k,m,m_k,o_k,tind
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_TimeInterval) :: ts
    integer                 :: tindex_f,tindex_t
    integer                 :: rc

    if(stats%selectOpt.eq.1.and.&
         ((trim(model%short_name).eq."TotalPrecip").or.&
         (trim(model%short_name).eq."Rainf")).and.&
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
                         
                         LVT_spi_struc(m)%model_pcp_final(t,tindex_f) = &
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
                            LVT_spi_struc(m)%obs_pcp_final(t,tindex_f) = &
                                 obs%value(t,m,o_k)
                         endif
                         
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
  end subroutine diagnoseSingleSPI

!BOP
! 
! !ROUTINE: LVT_computeSPI
! \label{LVT_computeSPI}
!
! !INTERFACE: 
  subroutine LVT_computeSPI(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the SPI values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelSPI](\ref{computeSingleModelSPI})
!     computes the SPI values for a single variable
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
          if(LVT_metrics%spi%selectOpt.eq.1.or.LVT_metrics%spi%timeOpt.eq.1) then 
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model))             
                call computeSingleSPIparams(alarm,&
                     model,obs,stats,LVT_metrics%spi)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
          endif
       elseif(pass.eq.2) then        
          if(LVT_metrics%spi%selectOpt.eq.1.or.LVT_metrics%spi%timeOpt.eq.1) then 
             if(alarm) then 
                if(LVT_metrics%spi%timeOpt.eq.1.and.&
                     LVT_metrics%spi%extractTS.eq.1) then 

                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%spi%ftn_ts_loc(i,m),200,advance='no') &
                                 LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                                 LVT_rc%hr,'',LVT_rc%mn, '' 
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%spi%ftn_ts_loc(i,1),200,advance='no') &
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
                call computeSingleSPI(alarm,&
                     model,obs,stats,LVT_metrics%spi)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
             
             if(alarm) then 
                if(LVT_metrics%spi%timeOpt.eq.1.and.&
                     LVT_metrics%spi%extractTS.eq.1) then 

                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%spi%ftn_ts_loc(i,m),fmt='(a1)') ''
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%spi%ftn_ts_loc(i,1),fmt='(a1)') ''
                      enddo
                   endif

                endif
             endif
          endif

       endif
    elseif(LVT_rc%startmode.eq."restart") then 
       if(LVT_metrics%spi%selectOpt.eq.1.or.LVT_metrics%spi%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%spi%timeOpt.eq.1.and.&
                  LVT_metrics%spi%extractTS.eq.1) then 

                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%spi%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%spi%ftn_ts_loc(i,1),200,advance='no') &
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
             call computeSingleSPI(alarm,&
                  model,obs,stats,LVT_metrics%spi)
             
             model => model%next
             obs => obs%next
             stats => stats%next
             
          enddo
          
          if(alarm) then 
             if(LVT_metrics%spi%timeOpt.eq.1.and.&
                  LVT_metrics%spi%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%spi%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%spi%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif
             endif
          endif
       endif
    endif


  end subroutine LVT_computeSPI
  

!BOP
! 
! !ROUTINE: computeSingleSPIparams
! \label{computeSingleSPIparams}
!
! !INTERFACE: 
  subroutine computeSingleSPIparams(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the SPI values for each grid cell. 
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
            ((trim(model%short_name).eq."TotalPrecip".or.&
            trim(model%short_name).eq."Rainf")).and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_spi_struc(m)%nasc
                      call spi_gamma(t,l,m,LVT_spi_struc(m)%nsize_total,&
                           LVT_spi_struc(m)%nsize_season,&
                           LVT_spi_struc(m)%model_pcp_final(t,:),&
                           LVT_spi_struc(m)%model_beta(t,l),&
                           LVT_spi_struc(m)%model_gamma(t,l),&
                           LVT_spi_struc(m)%model_pzero(t,l))
                   enddo
                enddo
             enddo
          enddo
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then 
                         do l=1,LVT_spi_struc(m)%nasc
                            call spi_gamma(t,l,m,LVT_spi_struc(m)%nsize_total,&
                                 LVT_spi_struc(m)%nsize_season,&
                                 LVT_spi_struc(m)%obs_pcp_final(t,:),&
                                 LVT_spi_struc(m)%obs_beta(t,l),&
                                 LVT_spi_struc(m)%obs_gamma(t,l),&
                                 LVT_spi_struc(m)%obs_pzero(t,l))
                         enddo
                      endif
                   endif
                enddo
             enddo
          enddo
       endif
    endif
  end subroutine computeSingleSPIparams

!BOP
! !ROUTINE:spi_gamma
! \label{spi_gamma}
! 
! !INTERFACE: 
  subroutine spi_gamma(t,l,kk,nsize,nsize_season,pcp,beta,gamma, pzero)

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: t
    integer, intent(in) :: l
    integer, intent(in) :: kk
    integer, intent(in) :: nsize
    integer, intent(in) :: nsize_season
    real,    intent(in) :: pcp(nsize)
    real                :: beta
    real                :: gamma
    real                :: pzero
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
!    \item[pcp]
!      array of precipitation values
!    \item[beta]
!      beta parameter of the fitted distribution
!    \item[gamma]
!      gamma parameter of the fitted distribution
!    \item[pzero]
!      probability of zero precipitation
!   \end{description}  
!EOP  
    integer             :: i,k,m,t1_ind
    real                :: pcp_season(nsize_season)
    real                :: tmppcp
    integer             :: pcp_count
    real                :: alpha 
    
!    real test1(10), rval

    pcp_season = 0.0
    k = 0 
    do i=l,nsize,LVT_spi_struc(kk)%nasc
       k = k+1

       tmppcp = 0.0
       pcp_count = 0 
       do m=1,LVT_spi_struc(kk)%tscale
          t1_ind = (i-LVT_spi_struc(kk)%tscale+m)
          if(t1_ind.gt.0) then 
             tmppcp = tmppcp + pcp(t1_ind)
             pcp_count = pcp_count + 1
          endif
       enddo
       pcp_season(k) = tmppcp/pcp_count
    enddo
    call gamma_fit(nsize_season, pcp_season, &
         alpha, beta,gamma, pzero)

#if 0
    print*, 'here'
    test1 = (/0.8147, 0.9058, 0.1270, 0.9134,0.6324,0.0975,0.2785,0.5469, 0.9575,0.9649/)
!    
    call gamma_fit(10, test1,&
         alpha, beta,gamma, pzero)

    print*, alpha, beta, gamma, pzero
!    beta = 0.2885
!    gamma =2.1627 
    call gamma_cdf(beta, gamma, pzero,test1(6),rval)
    print*, inv_normal(rval)
    stop
#endif
  end subroutine spi_gamma


!BOP
! 
! !ROUTINE: computeSingleSPI
! \label{computeSingleSPI}
!
! !INTERFACE: 
  subroutine computeSingleSPI(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the SPI values for each grid cell. 
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

    integer  :: t,l,k,i
    integer  :: tindex_f, tind
    real     :: model_spi(LVT_rc%ngrid,1,1)
    real     :: obs_spi(LVT_rc%ngrid,1,1)
    integer  :: count_model_spi(LVT_rc%ngrid,1,1)
    integer  :: count_obs_spi(LVT_rc%ngrid,1,1)

    real,    allocatable :: tavg_model_value_ts(:,:,:)
    real,    allocatable :: tavg_obs_value_ts(:,:,:)
    integer, allocatable :: count_tavg_model_value_ts(:,:,:)

    real     :: tmppcp
    integer  :: pcp_count
    integer  :: t1_ind,m

    if(stats%selectOpt.eq.1.and.&
         ((trim(model%short_name).eq."TotalPrecip".or.&
            trim(model%short_name).eq."Rainf")).and.&
         model%selectNlevs.ge.1) then 

       call getMonthlyTimeIndex(tindex_f, LVT_rc%tavgInterval)
       call getSeasonalTimeIndex(tind,LVT_rc%tavgInterval)

       do m=1,LVT_rc%nensem
          do i=LVT_spi_struc(m)%tscale,1,-1
             if(i.ne.1) then 
                LVT_spi_struc(m)%model_current(:,i) = &
                     LVT_spi_struc(m)%model_current(:,i-1)
             else
                LVT_spi_struc(m)%model_current(:,i) = & 
                     model%value(:,m,1)
             endif
          enddo
       enddo
       
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem         
             tmppcp = 0.0
             pcp_count = 0 

             do i=1,LVT_spi_struc(m)%tscale
                if(LVT_spi_struc(m)%model_current(t,i).ne.LVT_rc%udef) then 
                   tmppcp = tmppcp + LVT_spi_struc(m)%model_current(t,i)
                   pcp_count = pcp_count + 1
                endif
             enddo
             tmppcp = tmppcp/pcp_count
          
             call spi_zval(t,&
                  tmppcp, &
                  LVT_spi_struc(m)%model_beta(t,tind), &
                  LVT_spi_struc(m)%model_gamma(t,tind), &
                  LVT_spi_struc(m)%model_pzero(t,tind), &
                  LVT_spi_struc(m)%model_spi(t))
          enddo
       enddo
       if(LVT_rc%obssource(2).ne."none") then 
          if(obs%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem         
                do i=LVT_spi_struc(m)%tscale,1,-1
                   if(i.ne.1) then 
                      LVT_spi_struc(m)%obs_current(:,i) = &
                           LVT_spi_struc(m)%obs_current(:,i-1)
                   else
                      LVT_spi_struc(m)%obs_current(:,i) = & 
                           obs%value(:,m,1)
                   endif
                enddo
             enddo

             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem                         
                   tmppcp = 0.0
                   pcp_count = 0 
                   
                   do i=1,LVT_spi_struc(m)%tscale
                      if(LVT_spi_struc(m)%obs_current(t,i).ne.LVT_rc%udef) then 
                         tmppcp = tmppcp + LVT_spi_struc(m)%obs_current(t,i)
                         pcp_count = pcp_count + 1
                      endif
                   enddo
                   tmppcp = tmppcp/pcp_count
                   
                   call spi_zval(t,&
                        tmppcp,&
                        LVT_spi_struc(m)%obs_beta(t,tind), &
                        LVT_spi_struc(m)%obs_gamma(t,tind), &
                        LVT_spi_struc(m)%obs_pzero(t,tind), &
                        LVT_spi_struc(m)%obs_spi(t))
                   
                enddo
             enddo
          endif          
       endif

       count_model_spi = 1
       count_obs_spi = 1       
       if(metric%extractTS.eq.1) then 
          if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
             do m=1,LVT_rc%nensem
                model_spi(:,1,1) = LVT_spi_struc(m)%model_spi(:)
                if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                   obs_spi(:,1,1) = LVT_spi_struc(m)%obs_spi(:)
                   call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                        model,&
                        LVT_rc%ngrid,&
                        model_spi,&
                        count_model_spi,&
                        LVT_rc%ngrid,&
                        obs_spi,&
                        count_obs_spi)
                   
                else
                   call LVT_writeTSinfo(metric%ftn_ts_loc(:,m),&
                        model,&
                        LVT_rc%ngrid,&
                        model_spi,&
                        count_model_spi)
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
                
                do m=1,LVT_rc%nensem
                   do t=1,LVT_rc%ngrid
                      if(LVT_spi_struc(m)%model_spi(t).ne.LVT_rc%udef) then 
                         tavg_model_value_ts(t,1,1) = &
                              tavg_model_value_ts(t,1,1) + &
                              LVT_spi_struc(m)%model_spi(t)
                         tavg_obs_value_ts(t,1,1) = &
                              tavg_obs_value_ts(t,1,1) + &
                              LVT_spi_struc(m)%obs_spi(t)
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
                     count_model_spi,&
                     LVT_rc%ngrid,&
                     tavg_obs_value_ts,&
                     count_obs_spi)
             else
                do m=1,LVT_rc%nensem
                   do t=1,LVT_rc%ngrid
                      tavg_model_value_ts(t,1,1) = &
                           tavg_model_value_ts(t,1,1) + &
                           LVT_spi_struc(m)%model_spi(t)
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
                     count_model_spi)
             endif
             deallocate(tavg_model_value_ts)
             deallocate(count_tavg_model_value_ts)
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                deallocate(tavg_obs_value_ts)
             endif

          endif
       endif
         
       do m=1,LVT_rc%nensem
          call LVT_computeCI(LVT_spi_struc(m)%model_spi(:),LVT_rc%ngrid,&
               LVT_rc%pval_CI, LVT_spi_struc(m)%model_spi_ci(1,1))
          
          if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
             call LVT_computeCI(LVT_spi_struc(m)%obs_spi(:),LVT_rc%ngrid,&
                  LVT_rc%pval_CI,LVT_spi_struc(m)%obs_spi_ci(1,1))
          endif
       enddo

#if 0 
       call LVT_writeDataBasedStrat(model,obs,stats,metric,&
            LVT_rc%ngrid, model_spi)

       if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
               LVT_rc%ngrid, model_spi, &
               LVT_rc%ngrid, obs_spi)
          
       endif
#endif
    endif

  end subroutine computeSingleSPI

!BOP
! 
! !ROUTINE: LVT_writeMetric_SPI
! \label(LVT_writeMetric_SPI)
!
! !INTERFACE:
  subroutine LVT_writeMetric_SPI(pass,final,vlevels,stats,obs)
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
!    real                    :: model_spi(LVT_rc%ngrid,1,1)
!    real                    :: obs_spi(LVT_rc%ngrid,1,1)

    real, allocatable       :: value_ts(:,:)
    integer, allocatable    :: count_value_ts(:,:)

    integer                 :: count_model_spi(LVT_rc%ngrid)
    integer                 :: count_obs_spi(LVT_rc%ngrid)

    count_model_spi = LVT_rc%udef
    count_obs_spi = LVT_rc%udef

    if(pass.eq.LVT_metrics%spi%npass) then 
       if(final.ne.1) then
          if(stats%selectOpt.eq.1) then 
             
             allocate(value_ts(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(count_value_ts(LVT_rc%ngrid,LVT_rc%nensem))
             
             if(LVT_metrics%spi%timeOpt.eq.1) then 
                do k=1,vlevels
                   do m=1,LVT_rc%nensem
                      value_ts(:,m) = & 
                           LVT_spi_struc(m)%model_spi(:)
                      count_value_ts(:,m) = & 
                           count_model_spi(:)
                   enddo
                   call LVT_writevar_gridded(LVT_metrics%spi%ftn_ts, &
                        value_ts(:,:),&
                        stats%vid_ts(LVT_SPIid,1),k)
                   call LVT_writevar_gridded(LVT_metrics%spi%ftn_ts, &
                        real(count_value_ts(:,:)),&
                        stats%vid_count_ts(LVT_SPIid,1),k)
                enddo

                if(LVT_rc%obssource(2).ne."none") then 
                   value_ts = LVT_rc%udef
                   count_value_ts = 0
             
                   do k=1,vlevels
                      do m=1,LVT_rc%nensem
                         value_ts(:,m) = & 
                              LVT_spi_struc(m)%obs_spi(:)
                         count_value_ts(:,m) = & 
                              count_obs_spi(:)
                      enddo
                      call LVT_writevar_gridded(LVT_metrics%spi%ftn_ts, &
                           value_ts(:,:),&
                           stats%vid_ts(LVT_SPIid,2),k)
                      call LVT_writevar_gridded(LVT_metrics%spi%ftn_ts, &
                           real(count_value_ts(:,:)),&
                           stats%vid_count_ts(LVT_SPIid,2),k)
                   enddo
                endif

                deallocate(value_ts)
                deallocate(count_value_ts)

             endif
             
          endif
       endif
    end if
#if 0 
             model_spi(:,1,1) = LVT_spi_struc(m)%model_spi(:)
             call LVT_writeSummaryStats(&
                  LVT_metrics%spi%ftn_summ,&
                  LVT_metrics%spi%short_name,&
                  LVT_rc%ngrid,&
                  model_spi, &
                  count_model_spi,&
                  stats%short_name,&
                  LVT_spi_struc(m)%model_spi_ci)
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                obs_spi(:,1,1) = LVT_spi_struc(m)%obs_spi(:)
                call LVT_writeSummaryStats(&
                     LVT_metrics%spi%ftn_summ,&
                     LVT_metrics%spi%short_name,&
                     LVT_rc%ngrid,&
                     obs_spi, &
                     count_obs_spi,&
                     "OBS_"//trim(stats%short_name),&
                     LVT_spi_struc(m)%obs_spi_ci)
             endif
#endif
  end subroutine LVT_writeMetric_SPI

!BOP
! 
! !ROUTINE: LVT_resetMetric_SPI
! \label(LVT_resetMetric_SPI)
!
! !INTERFACE:
  subroutine LVT_resetMetric_SPI
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
            (model%short_name.eq.&
            "TotalPrecip".or.&
            trim(model%short_name).eq."Rainf")) then 
          do m=1,LVT_rc%nensem
             if(LVT_metrics%spi%timeOpt.eq.1) then 
                LVT_spi_struc(m)%model_spi  = LVT_rc%udef 
             endif
             if(LVT_rc%obssource(2).ne."none".and.&
                  obs%selectNlevs.ge.1) then 
                LVT_spi_struc(m)%obs_spi  = LVT_rc%udef 
             endif
          enddo
       endif

       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo

  end subroutine LVT_resetMetric_SPI

!BOP
! 
! !ROUTINE: LVT_writerestart_SPI
! 
! !INTERFACE:
  subroutine LVT_writerestart_SPI(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for SPI metric computations
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
          if((trim(model%short_name).eq."TotalPrecip".or.&
            trim(model%short_name).eq."Rainf")) then 
             do m=1,LVT_rc%nensem
                if(LVT_metrics%spi%selectOpt.eq.1) then 
                   do l=1,LVT_spi_struc(m)%nasc
                      call LVT_writevar_restart(ftn,&
                           LVT_spi_struc(m)%model_beta(:,l))
                      call LVT_writevar_restart(ftn,&
                           LVT_spi_struc(m)%model_gamma(:,l))
                      call LVT_writevar_restart(ftn,&
                           LVT_spi_struc(m)%model_pzero(:,l))
                   enddo
                   if(LVT_rc%obssource(2).ne."none") then 
                      do l=1,LVT_spi_struc(m)%nasc
                         call LVT_writevar_restart(ftn,&
                              LVT_spi_struc(m)%obs_beta(:,l))
                         call LVT_writevar_restart(ftn,&
                              LVT_spi_struc(m)%obs_gamma(:,l))
                         call LVT_writevar_restart(ftn,&
                              LVT_spi_struc(m)%obs_pzero(:,l))
                      enddo
                   endif                
                end if
             enddo
          endif
          model => model%next
          obs   => obs%next
          stats => stats%next
       enddo
!    endif

  end subroutine LVT_writerestart_SPI

!BOP
! 
! !ROUTINE: LVT_readrestart_SPI
! 
! !INTERFACE:
  subroutine LVT_readrestart_SPI(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for SPI metric computations
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
       if((trim(model%short_name).eq."TotalPrecip".or.&
            trim(model%short_name).eq."Rainf")) then 
          if(LVT_metrics%spi%selectOpt.eq.1) then 
!             do l=1,LVT_spi_struc(m)%nsize_total
!                call LVT_readvar_restart(ftn,&
!                     LVT_spi_struc(m)%model_pcp_final(:,l))
!             enddo
             do m=1,LVT_rc%nensem
                do l=1,LVT_spi_struc(m)%nasc
                   call LVT_readvar_restart(ftn,&
                        LVT_spi_struc(m)%model_beta(:,l))
                   call LVT_readvar_restart(ftn,&
                        LVT_spi_struc(m)%model_gamma(:,l))
                   call LVT_readvar_restart(ftn,&
                        LVT_spi_struc(m)%model_pzero(:,l))
                enddo
                if(LVT_rc%obssource(2).ne."none") then 
!                do l=1,LVT_spi_struc(m)%nsize_total
!                   call LVT_readvar_restart(ftn,&
!                        LVT_spi_struc(m)%obs_pcp_final(:,l))
!                enddo
                   do l=1,LVT_spi_struc(m)%nasc
                      call LVT_readvar_restart(ftn,&
                           LVT_spi_struc(m)%obs_beta(:,l))
                      call LVT_readvar_restart(ftn,&
                           LVT_spi_struc(m)%obs_gamma(:,l))
                      call LVT_readvar_restart(ftn,&
                           LVT_spi_struc(m)%obs_pzero(:,l))
                   enddo
                endif
             enddo
          end if
       end if
          
       model => model%next
       obs   => obs%next
       stats => stats%next

    enddo

  end subroutine LVT_readrestart_SPI

!BOP
! 
! !ROUTINE: spi_zval
!  \label{spi_zval}
! 
! !INTERFACE: 
  subroutine spi_zval(t,pcp, beta, gamma, pzero, zval)
! !ARGUMENTS:
    integer, intent(in) :: t
    real,    intent(in) :: pcp
    real                :: beta, gamma, pzero
    real                :: rval
    real                :: zval
! 
! !DESCRPTION:
!   This subroutine computes the SPI values. The cumulative 
!   probability of a given event is computed from the given 
!   gamma distribution. The cumulative probability is then 
!   transformed to the standard normal random variable (SPI) 
!   with mean zero and variance of one, which the value of 
!   SPI. 
!
!   \begin{description}
!    \item[t]
!      index of the grid/tile
!    \item[accum\_pcp]
!      accumulated precipitation value
!    \item[beta]
!      beta parameter of the fitted distribution
!    \item[gamma]
!      gamma parameter of the fitted distribution
!    \item[pzero]
!      probability of zero precipitation
!    \item[zval]
!     transformed standard normal random variable
!   \end{description}  
!EOP
    if(beta.eq.0) then 
       zval = LVT_rc%udef
    else
       call gamma_cdf(beta, gamma, pzero, pcp,rval)
       zval = inv_normal(rval)
    endif
  end subroutine spi_zval
 
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
       write(LVT_logunit,*) ' Programg stopping ... '
       call LVT_endrun()
    endif
  end subroutine getSeasonalTimeIndex

!BOP
! 
! !ROUTINE: gamma_fit
! \label{gamma_fit}
! 
! !INTERFACE: 
  subroutine gamma_fit(nsize, pcp, alpha, beta, gamma, pzero)

    implicit none
! !ARGUMENTS:     
    integer,         intent(in) :: nsize
    real,            intent(in) :: pcp(nsize)
    real                        :: alpha
    real                        :: beta
    real                        :: gamma
    real                        :: pzero
! 
! !DESCRIPTION: 
!   This subroutine estimates incomplete gamma distribution 
!   parameters. 
! 
!EOP    
    integer                     :: i 
    real                        :: sumlog,sumv
    integer                     :: nsumlog
    real                        :: mn

    sumlog = 0.0
    sumv = 0.0
    nsumlog = 0 
    pzero = 0 

    do i=1,nsize
       if(pcp(i).gt.0) then 
          sumlog = sumlog+log(pcp(i))
          sumv = sumv + pcp(i)
          nsumlog = nsumlog+1
       else
          pzero = pzero + 1
       endif
    enddo

    pzero = pzero/nsize

    if(nsumlog.ne.0) then 
       mn = sumv/nsumlog
    endif

    if(nsumlog.eq.1) then !bogus data, do something reasonable. 
       alpha = 0.0
       beta = mn
       gamma = 1.0
    elseif(nsumlog.eq.0) then !data is all zeroes. 
       alpha = 0.0
       beta = 0.0   
       gamma = 1.0
    else !maximum likelihood estimate
       alpha = log(mn)-sumlog/nsumlog
       if(alpha.ne.0) then 
          gamma = (1.0 + sqrt(1.0+4.0*alpha/3.0)) /(4.0*alpha)
          beta = mn/gamma
       else
          beta = 0.0
          gamma = 1.0
       endif
    endif    

    if(beta.lt.0) then 
       beta = 0.0
       gamma = 1.0
       alpha = 0.0
    endif


  end subroutine gamma_fit

!BOP
! !ROUTINE: gamma_cdf
! \label{gamma_cdf}
! 
! !INTERFACE: 
  subroutine gamma_cdf(beta, gamma, pzero, x,rval)

    implicit none
! !ARGUMENTS: 
    real               :: beta
    real               :: gamma
    real               :: pzero
    real               :: x
    real               :: rval
! 
! !DESCRIPTION: 
!
!EOP
    real               :: tval

    if(x<= 0.0) then
       rval = pzero
    else
       call gamma_P(gamma,x/beta,tval)
       rval = pzero+(1-pzero)*tval
    endif

  end subroutine gamma_cdf
!BOP
! 
! !ROUTINE: gamma_P
! \label{gamma_P}
! 
! !INTERFACE: 
  subroutine gamma_P(a,x,rval)
! 
! !DESCRIPTION: 
!  Evaluate the incomplete gamma function P(a,x), choosing the 
!  most appropriate representation
! 
!EOP   
    implicit none

    real           :: a
    real           :: x
    real           :: rval

    if(x .lt. (a+1.0)) then 
       call gammser(a,x,rval)
    else
       call gammcf(a,x,rval)
       rval = 1.0-rval
    endif

  end subroutine gamma_P

!BOP
! 
! !ROUTINE: gamma_Q
! \label{gamma_Q}
! 
! !INTERFACE: 
  subroutine gamma_Q(a,x,rval)
! 
! !DESCRIPTION: 
!  Evaluate the incomplete gamma function Q(a,x), choosing the 
!  most appropriate representation
! 
!EOP   
    implicit none

    real           :: a
    real           :: x
    real           :: rval

    if(x .lt. (a+1.0)) then 
       call gammser(a,x,rval)
       rval = 1.0 -rval
    else
       call gammcf(a,x,rval)
    endif

  end subroutine gamma_Q

!BOP
! 
! !ROUTINE: gammser
! \label{gammser}
!
! !INTERFACE: 
  subroutine gammser(a,x,rval)

    implicit none
! !ARGUMENTS: 
    
    real             :: a
    real             :: x
    real             :: rval

! 
! !DESCRIPTION: 
!  Evaluate P(a,x) by its series representation
!
!EOP
    integer      :: n
    real         :: gln
    real         :: ap,sumv,del
    integer      :: iwarn
    
    iwarn = 0 
!    gln = log(gamma(a))
    gln = gammln(a)

    if(x.eq.0.0) then 
       rval = 0.0
    else
       ap = a
       sumv = 1.0/a
       del =sumv
       
       do n=1,maxiter
          ap = ap + 1
          del = del*(x/ap)
          sumv = sumv + del
          if(abs(del) .lt. epsilon*abs(sumv)) then
             exit;
          endif          
       enddo
       rval = sumv*exp(-x+a*log(x)-gln)
!       if(abs(a*log(x)-gln-x).gt.10) then 
!          rval = 0.0
!       else
!          print*,sumv, x, a, gln, a*log(x),(-x+a*log(x)-gln)
!       endif

    endif
  end subroutine gammser
  

!BOP
! 
! !ROUTINE: gammcf
! \label{gammcf}
!
! !INTERFACE: 
  subroutine gammcf(a,x,rval)

    implicit none
! !ARGUMENTS: 
    
    real             :: a
    real             :: x
    real             :: rval

! 
! !DESCRIPTION: 
!  Evaluate P(a,x) in its continued fraction representation. 
!
!EOP
    integer           :: n
    real              :: gln,fac,gold
    real              :: a0,a1,an,anf,ana,g,b0,b1

!    gln = log(gamma(a))
    gln = gammln(a)

    g = 0.0
    gold = 0.0
    a0 = 1.0
    a1 = x
    b0 = 0.0
    b1 = 1.0
    fac = 1.0

    do n=1,maxiter
       an = n
       ana = an-a
       a0 = (a1+a0*ana)*fac
       b0 = (b1+b0*ana)*fac
       anf = an*fac
       a1 = x*a0 + anf*a1
       b1 = x*b0 + anf*b1
       
       if(a1 .ne.0.0) then 
          fac = 1.0/a1
          g = b1*fac
          if(abs((g-gold)/g).lt.epsilon) then 
             exit;
          else
             gold = g
          endif
       endif
    enddo
    rval = g*exp(-x+a*log(x)-gln)
!    if(abs(a*log(x)-gln-x).gt.10) then 
!       rval = 0.0
!    else
!       print*, g, x, a, gln, a*log(x),(-x+a*log(x)-gln)
!       rval = g*exp(-x+a*log(x)-gln)
!    endif

  end subroutine gammcf


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!   input prob; return z.
!
!   See Abromowitz and Stegun _Handbook of Mathematical Functions_, p. 933
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  real function inv_normal (prob)
    real         :: prob, prob1
    real         :: c0, c1, c2, d1, d2, d3
    real         :: t, sign

    data c0, c1, c2 /2.515517, 0.802853, 0.010328/
    data d1, d2, d3 /1.432788, 0.189269, 0.001308/
    
    prob1 = prob
    if (prob .gt. 0.5) then
       sign = 1.0
       prob = 1.0 - prob
    else
       sign = -1.0
    endif
    
    if (prob .lt. 0.0) then
       write(0, *) 'Error in inv_normal(). Prob. not in [0,1.0]',prob
       inv_normal = 0.0
       return
    endif

!    if (prob .eq. 0.0) then
!       inv_normal = 1.0e37 * sign
    if (prob .lt.1E-5) then
       inv_normal = 0.0
       return
!       inv_normal = 0.0
!       return
    endif
    
    t = sqrt(alog (1.0 / (prob * prob)))
    inv_normal = (sign * (t - ((((c2 * t) + c1) * t) + c0) / &
         ((((((d3 * t) + d2) * t) + d1) * t) + 1.0)))
    return
  end function inv_normal

!
!     For those who don't have a ln(gamma) function.
!
  function gammln(xx)

    implicit none

    real cof(6)
    data cof /76.18009173, -86.50532033, 24.01409822, -1.231739516, &
         0.120858003e-2, -0.536382e-5/          
    real xx
    real x 
    real tmp
    real ser
    integer j 
    real gammln

    x = xx - 1.0
    tmp = x + 5.5
    tmp = tmp - (x+0.5) * log (tmp)
    ser = 1.0
    do j = 1, 5
       x = x + 1.0
       ser = ser + cof(j) / x
    enddo
    gammln = -tmp + log (2.50662827465 * ser)
    return
  end function gammln
    
  REAL FUNCTION GAMMA(X)
!    DOUBLE PRECISION FUNCTION DGAMMA(X)
!----------------------------------------------------------------------
!
! This routine calculates the GAMMA function for a real argument X.
!   Computation is based on an algorithm outlined in reference 1.
!   The program uses rational functions that approximate the GAMMA
!   function to at least 20 significant decimal digits.  Coefficients
!   for the approximation over the interval (1,2) are unpublished.
!   Those for the approximation for X .GE. 12 are from reference 2.
!   The accuracy achieved depends on the arithmetic system, the
!   compiler, the intrinsic functions, and proper selection of the
!   machine-dependent constants.
!
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
! XBIG   - the largest argument for which GAMMA(X) is representable
!          in the machine, i.e., the solution to the equation
!                  GAMMA(XBIG) = beta**maxexp
! XINF   - the largest machine representable floating-point number;
!          approximately beta**maxexp
! EPS    - the smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0
! XMININ - the smallest positive floating-point number such that
!          1/XMININ is machine representable
!
!     Approximate values for some important machines are:
!
!                            beta       maxexp        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! Cyber 180/855
!   under NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-Format   (D.P.)        2          127        34.844
! VAX G-Format   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! Cyber 180/855
!   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns the value XINF for singularities or
!     when overflow would occur.  The computation is believed
!     to be free of underflow and overflow.
!
!
!  Intrinsic functions required are:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! References: "An Overview of Software Development for Special
!              Functions", W. J. Cody, Lecture Notes in Mathematics,
!              506, Numerical Analysis Dundee, 1975, G. A. Watson
!              (ed.), Springer Verlag, Berlin, 1976.
!
!              Computer Approximations, Hart, Et. Al., Wiley and
!              sons, New York, 1968.
!
!  Latest modification: October 12, 1989
!
!  Authors: W. J. Cody and L. Stoltz
!           Applied Mathematics Division
!           Argonne National Laboratory
!           Argonne, IL 60439
!
!----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
      REAL C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,&
           TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,&
           SQRTPI/0.9189385332046727417803297E0/,&
           PI/3.1415926535897932384626434E0/
!----------------------------------------------------------------------
!  Machine dependent parameters
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,&
           XINF/3.4E38/

!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     approximation over (1,2).
!----------------------------------------------------------------------
    DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,&
         -3.79804256470945635097577E+2,6.29331155312818442661052E+2,&
         8.66966202790413211295064E+2,-3.14512729688483675254357E+4,&
         -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,&
         -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,&
         .25381184209801510330112D+4,4.75584627752788110767815D+3,&
         -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
!----------------------------------------------------------------------
!  Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------
    DATA C/-1.910444077728E-03,8.4171387781295E-04,&
         -5.952379913043012E-04,7.93650793500350248E-04,&
         -2.777777777777681622553E-03,8.333333333333333331554247E-02,&
         .7083835261E-03/

!----------------------------------------------------------------------
!  Statement functions for conversion between integer and float
!----------------------------------------------------------------------
    CONV(I) = REAL(I)
    PARITY = .FALSE.
    FACT = ONE
    N = 0
    Y = X
    IF (Y .LE. ZERO) THEN
!----------------------------------------------------------------------
!  Argument is negative
!----------------------------------------------------------------------
       Y = -X
       Y1 = AINT(Y)
       RES = Y - Y1
       IF (RES .NE. ZERO) THEN
          IF (Y1 .NE. AINT(Y1*HALF)*TWO) PARITY = .TRUE.
          FACT = -PI / SIN(PI*RES)
          Y = Y + ONE
       ELSE
          RES = XINF
          GO TO 900
       END IF
    END IF
!----------------------------------------------------------------------
!  Argument is positive
!----------------------------------------------------------------------
    IF (Y .LT. EPS) THEN
!----------------------------------------------------------------------
!  Argument .LT. EPS
!----------------------------------------------------------------------
       IF (Y .GE. XMININ) THEN
          RES = ONE / Y
       ELSE
          RES = XINF
          GO TO 900
       END IF
    ELSE IF (Y .LT. TWELVE) THEN
       Y1 = Y
       IF (Y .LT. ONE) THEN
!----------------------------------------------------------------------
!  0.0 .LT. argument .LT. 1.0
!----------------------------------------------------------------------
          Z = Y
          Y = Y + ONE
       ELSE
!----------------------------------------------------------------------
!  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
!----------------------------------------------------------------------
          N = INT(Y) - 1
          Y = Y - CONV(N)
          Z = Y - ONE
       END IF
!----------------------------------------------------------------------
!  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
!----------------------------------------------------------------------
       XNUM = ZERO
       XDEN = ONE
       DO I = 1, 8
          XNUM = (XNUM + P(I)) * Z
          XDEN = XDEN * Z + Q(I)
       ENDDO
       RES = XNUM / XDEN + ONE
       IF (Y1 .LT. Y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  0.0 .LT. argument .LT. 1.0
!----------------------------------------------------------------------
          RES = RES / Y1
       ELSE IF (Y1 .GT. Y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  2.0 .LT. argument .LT. 12.0
!----------------------------------------------------------------------
          DO  I = 1, N
             RES = RES * Y
             Y = Y + ONE
          ENDDO
       END IF
    ELSE
!----------------------------------------------------------------------
!  Evaluate for argument .GE. 12.0,
!----------------------------------------------------------------------
       IF (Y .LE. XBIG) THEN
          YSQ = Y * Y
          SUM = C(7)
          DO I = 1, 6
             SUM = SUM / YSQ + C(I)
          ENDDO
          SUM = SUM/Y - Y + SQRTPI
          SUM = SUM + (Y-HALF)*LOG(Y)
          RES = EXP(SUM)
       ELSE
          RES = XINF
          GO TO 900
       END IF
    END IF
!----------------------------------------------------------------------
!  Final adjustments and return
!----------------------------------------------------------------------
    IF (PARITY) RES = -RES
    IF (FACT .NE. ONE) RES = FACT / RES
900 GAMMA = RES
    RETURN
! ---------- Last line of GAMMA ----------
  END FUNCTION GAMMA

end module LVT_SPIMod
