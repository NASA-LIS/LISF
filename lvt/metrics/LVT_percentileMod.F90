!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! 
! !MODULE: LVT_percentileMod
! \label(LVT_percentileMod)
!
! !INTERFACE:
module LVT_percentileMod
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
!  This module handles the computation of the percentiles for a given variable
!  against a corresponding climatology. 
!  
!  The code makes two passes through the data. During the first pass, it
!  generates the percentile climatology and during the second pass, it computes
!  daily percentiles using the percentile climatology. Note that code 
!  requires the use of a time averaging interval of a day. 
! 
!  The general methodology for assembling percentile climatology is as follows:
!  For each grid point, the variable values for the calendar day 
!  across all years is assembled. For example, assume that a 10 year simulation
!  is used for the analysis (say 2000-2010) and we are computing percentiles 
!  for soil moisture. For day 1, 10 soil moisture values are assembled by 
!  using the Jan 1 values for year 2000, 2001, 2002,... to 2010. Similarly 
!  10 soil moisture values are assembled for 365 days. 
!
!  Once the percentile climatology is assembled, this array is then sorted 
!  in the ascending order. During the second pass through the data, each 
!  day's soil moisture value is then ranked against the percentile climatology. 
!  For example, soil moisture value from Jan 1, 2000 is ranked against
!  the 10 sorted percentile climatology values and the percentile is determined.
! 
!  The particular implementation used here follows an extension of the 
!  above-mentioned strategy, where a moving window of 5 days is employed 
!  to improve the sampling density while calculating percentiles. While
!  assembling the percentile climatology, instead of using a single day 
!  across all years, we use 5 days (2 previous days, current day and 2 
!  next days) from each year. For the 10 year example, we will have 
!  5 days for Jan 3 2000, using Jan 1 to Jan 5 values. As a result 
!  across the 10 years, each day will have 5x10 = 50 values. This sample
!  is then used to compute the percentiles as before. 
!
!
!  NOTE: Because of the long temporal dimension, the memory becomes limiting
!  to do these calculations all in memory. LVT gets around this issue by 
!  performing intermediate I/O.   
!
! !FILES USED:
!
! !REVISION HISTORY: 
!  22 Mar 2012    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
  implicit none

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initpercentile
  public :: LVT_diagnosepercentile
  public :: LVT_computepercentile
  public :: LVT_writeMetric_percentile
  public :: LVT_resetMetric_percentile
  public :: LVT_writerestart_percentile
  public :: LVT_readrestart_percentile

!EOP

  type, public :: pctiledec
     integer          :: climo_calc_only
     integer          :: nyears
     integer          :: nwindow
     integer          :: nsize_total
     real             :: d_classes(5)
     integer, allocatable :: ftn_ts_darea(:)
     integer          :: usescaling
     character*100    :: mean_filename1
     character*100    :: mean_filename2
     character*100    :: std_filename1
     character*100    :: std_filename2
     type(ESMF_Time)  :: stime, etime
     real,    allocatable :: mean1(:)
     real,    allocatable :: mean2(:)
     real,    allocatable :: std1(:)
     real,    allocatable :: std2(:)
  end type pctiledec

  type(pctiledec), allocatable, save :: LVT_pctile_struc(:)
  private

contains
!BOP
!
! !ROUTINE: LVT_initpercentile
! \label{LVT_initpercentile}
! 
! !INTERFACE: 
  subroutine LVT_initpercentile(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This routine initializes the data structures required
!  for the percentile computations. 
! 
!EOP

    type(ESMF_Time)         :: startTime
    type(ESMF_Time)         :: stopTime
    type(ESMF_TimeInterval) :: timeStep
    integer                 :: i,k,t,c,r
    character*100           :: filename
    character*3             :: ftime
    character*4             :: fens
    integer                 :: rc
    integer                 :: m
    integer                 :: ftn
    integer                 :: nsize, nyears
    integer                 :: stime(6),etime(6)
    real                    :: var(LVT_rc%lnc,LVT_rc%lnr)
    real, allocatable       :: pctile_model(:,:,:)

    integer          :: usescaling
    character*100    :: mean_filename1
    character*100    :: mean_filename2
    character*100    :: std_filename1
    character*100    :: std_filename2

    allocate(LVT_pctile_struc(LVT_rc%nensem))
    allocate(stats%pctile(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       
       if(LVT_metrics%percentile%selectOpt.eq.1.or.&
            LVT_metrics%percentile%timeOpt.eq.1) then 
          
          metric%npass = 2   

          LVT_pctile_struc(m)%d_classes(1) = 0.30
          LVT_pctile_struc(m)%d_classes(2) = 0.20
          LVT_pctile_struc(m)%d_classes(3) = 0.10
          LVT_pctile_struc(m)%d_classes(4) = 0.05
          LVT_pctile_struc(m)%d_classes(5) = 0.02

          allocate(stats%pctile(m)%value_model_nsize(LVT_rc%ngrid, &
               365))
          stats%pctile(m)%value_model_nsize = 0 

          allocate(stats%pctile(m)%value_model(LVT_rc%ngrid,selectNlevs(1),1))
          stats%pctile(m)%value_model  = 0
          allocate(stats%pctile(m)%value_count_model(LVT_rc%ngrid,selectNlevs(1),1))
          stats%pctile(m)%value_count_model  = 0

          allocate(stats%pctile(m)%tavg_value_model(LVT_rc%ngrid,selectNlevs(1),1))
          stats%pctile(m)%tavg_value_model  = 0
          allocate(stats%pctile(m)%tavg_count_value_model(LVT_rc%ngrid,selectNlevs(1),1))
          stats%pctile(m)%tavg_count_value_model  = 0

          allocate(stats%pctile(m)%value_model_ci(selectNlevs(1), 1))

          if(LVT_rc%tavgInterval.ne.86400) then
             write(LVT_logunit,*) '[ERR] Percentile computations are only supported'
             write(LVT_logunit,*) '[ERR] with a daily average interval...'
             call LVT_endrun()
          endif
          LVT_pctile_struc(m)%nwindow = 5


          call system("mkdir -p "//trim(LVT_rc%statsodir))       
          call system("mkdir -p "//trim(LVT_rc%statsodir)//'/RST')       
       endif

    enddo

    if(LVT_metrics%percentile%selectOpt.eq.1.or.&
         LVT_metrics%percentile%timeOpt.eq.1) then 
       if(LVT_rc%startmode.eq."coldstart") then 
          
          call LVT_computeTimeSpanInYears(nyears)
          LVT_pctile_struc(:)%nyears = nyears
             
          do m=1,LVT_rc%nensem
             LVT_pctile_struc(m)%nsize_total = LVT_pctile_struc(m)%nwindow * & 
               LVT_pctile_struc(m)%nyears
          enddo

          allocate(pctile_model(LVT_rc%ngrid,&
               LVT_rc%nensem,&
               LVT_pctile_struc(1)%nsize_total))
          pctile_model = LVT_rc%udef
          
          do k=1,365
             write(unit=ftime,fmt='(I3.3)') k
             filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//&
                  '_pctile_climo.bin'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=trim(filename),form='unformatted')
             write(ftn) LVT_pctile_struc(1)%nsize_total
             write(ftn) pctile_model
             call LVT_releaseUnitNumber(ftn)
          enddo
       else     
          
          metric%npass = 1
          
          do k=1,365
             write(unit=ftime,fmt='(I3.3)') k
             filename = trim(LVT_rc%statsodir)//'/RST/'//&
                  trim(ftime)//'_pctile_climo.bin'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=trim(filename),form='unformatted')
             read(ftn) nsize
             LVT_pctile_struc(:)%nsize_total = nsize

             if(k.eq.1) then 
                allocate(pctile_model(LVT_rc%ngrid,&
                     LVT_rc%nensem,&
                     LVT_pctile_struc(1)%nsize_total))
                pctile_model = LVT_rc%udef
             endif
             
             read(ftn) pctile_model
             call LVT_releaseUnitNumber(ftn)             
             
             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem
                   call compute_sorted_sizes(&
                        LVT_pctile_struc(m)%nsize_total,&
                        pctile_model(t,m,:),&
                        stats%pctile(m)%value_model_nsize(t,k))
                enddo
             enddo
          enddo
       endif

       deallocate(pctile_model)
    endif
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------
    
    metric%obsData = .true. 

    if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then     
       do m=1,LVT_rc%nensem
          
          allocate(LVT_pctile_struc(m)%ftn_ts_darea(LVT_rc%ntslocs))
          write(fens,fmt='(i4.4)') m
          
          do i=1,LVT_rc%ntslocs
             
             filename = trim(LVT_rc%statsodir)//'/'//&
                  'PercentArea_percentile_'//&
                  trim(LVT_TSobj(i)%tslocname)//'_'//&
                  trim(fens)//&
                  '.dat'
             
             LVT_pctile_struc(m)%ftn_ts_darea(i) = LVT_getNextUnitNumber()
             
             open(LVT_pctile_struc(m)%ftn_ts_darea(i),file=(filename),&
                  form='formatted')
             
          enddo
       enddo
    else
       allocate(LVT_pctile_struc(1)%ftn_ts_darea(LVT_rc%ntslocs))
       do i=1,LVT_rc%ntslocs
          
          filename = trim(LVT_rc%statsodir)//'/'//&
               'PercentArea_percentile_'//&
               trim(LVT_TSobj(i)%tslocname)//&
               '.dat'
          
          LVT_pctile_struc(1)%ftn_ts_darea(i) = LVT_getNextUnitNumber()
          
          open(LVT_pctile_struc(1)%ftn_ts_darea(i),file=(filename),&
               form='formatted')
          
       enddo
    endif
!-------------------------------------------------------------------------
! LVT can be used to compute/establish the climatology alone. If
! LVT is used to compute percentiles based on an established climo, 
! then the 'restart' mode must be used. 
!-------------------------------------------------------------------------

    call ESMF_ConfigGetAttribute(LVT_config,&
         usescaling,&
         label="Scale model data prior to computing percentiles:",rc=rc)
    call LVT_verify(rc,"Scale model data prior to computing percentiles: option not specified in the config file")
    
    LVT_pctile_struc(:)%usescaling = usescaling

    if(LVT_pctile_struc(1)%usescaling.eq.1) then 
       call ESMF_ConfigGetAttribute(LVT_config,&
            mean_filename1,&
            label="Percentile scaling mean (input data) filename:",rc=rc)
       call LVT_verify(rc,"Percentile scaling mean (input data) filename: not specified")
       LVT_pctile_struc(:)%mean_filename1 = mean_filename1

       call ESMF_ConfigGetAttribute(LVT_config,&
            std_filename1,&
            label="Percentile scaling standard deviation (input data) filename:",rc=rc)
       call LVT_verify(rc,"Percentile scaling standard deviation (input data) filename: not specified")
       
       LVT_pctile_struc(:)%std_filename1 = std_filename1

       call ESMF_ConfigGetAttribute(LVT_config,&
            mean_filename2,&
            label="Percentile scaling mean (scaled data) filename:",rc=rc)
       call LVT_verify(rc,"Percentile scaling mean (scaled data) filename: not specified")
       LVT_pctile_struc(:)%mean_filename2 = mean_filename2

       call ESMF_ConfigGetAttribute(LVT_config,&
            std_filename2,&
            label="Percentile scaling standard deviation (scaled data) filename:",rc=rc)
       call LVT_verify(rc,"Percentile scaling standard deviation (scaled data) filename: not specified")
       
       LVT_pctile_struc(:)%std_filename2 = std_filename2

       call ESMF_ConfigFindLabel(LVT_config,&
            "Percentile scaling start time for scaling:",rc=rc)
       do i=1,6
          call ESMF_ConfigGetAttribute(LVT_config,&
               stime(i),rc=rc)
          call LVT_verify(rc,"Percentile scaling start time for scaling: not defined")
       enddo

       do m=1,LVT_rc%nensem
          call ESMF_TimeSet(LVT_pctile_struc(m)%stime, &
               yy=stime(1), mm=stime(2), dd=stime(3),&
               h = stime(4), m = stime(5), s = stime(6), &
               calendar=LVT_calendar, rc=rc)
          call LVT_verify(rc, 'error in ESMF_TimeSet: start time')
       enddo

       call ESMF_ConfigFindLabel(LVT_config,&
            "Percentile scaling end time for scaling:",rc=rc)
       do i=1,6
          call ESMF_ConfigGetAttribute(LVT_config,&
               etime(i),rc=rc)
          call LVT_verify(rc,"Percentile scaling end time for scaling: not defined")
       enddo
       do m=1,LVT_rc%nensem
          call ESMF_TimeSet(LVT_pctile_struc(m)%etime, &
               yy=etime(1), mm=etime(2), dd=etime(3),&
               h = etime(4), m = etime(5), s = etime(6), &
               calendar=LVT_calendar, rc=rc)
          call LVT_verify(rc, 'error in ESMF_TimeSet: end time')
       enddo

       do m=1,LVT_rc%nensem
          allocate(LVT_pctile_struc(m)%mean1(LVT_rc%ngrid))
          allocate(LVT_pctile_struc(m)%mean2(LVT_rc%ngrid))
          allocate(LVT_pctile_struc(m)%std1(LVT_rc%ngrid))
          allocate(LVT_pctile_struc(m)%std2(LVT_rc%ngrid))

          write(LVT_logunit,*) '[INFO] Reading mean file ',&
               trim(LVT_pctile_struc(m)%mean_filename1)
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=LVT_pctile_struc(m)%mean_filename1, form='unformatted')
          read(ftn) var
          call LVT_releaseUnitNumber(ftn)
          
          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(LVT_domain%gindex(c,r).ne.-1) then 
                   LVT_pctile_struc(m)%mean1(LVT_domain%gindex(c,r)) = var(c,r) 
                endif
             enddo
          enddo

          write(LVT_logunit,*) '[INFO] Reading mean file ',&
               trim(LVT_pctile_struc(m)%mean_filename2)
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=LVT_pctile_struc(m)%mean_filename2, form='unformatted')
          read(ftn) var
          call LVT_releaseUnitNumber(ftn)

          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(LVT_domain%gindex(c,r).ne.-1) then 
                   LVT_pctile_struc(m)%mean2(LVT_domain%gindex(c,r)) = var(c,r) 
                endif
             enddo
          enddo

          write(LVT_logunit,*) '[INFO] Reading standard deviation file ',&
               trim(LVT_pctile_struc(m)%std_filename1)
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=LVT_pctile_struc(m)%std_filename1, form='unformatted')
          read(ftn) var
          call LVT_releaseUnitNumber(ftn)

          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(LVT_domain%gindex(c,r).ne.-1) then 
                   LVT_pctile_struc(m)%std1(LVT_domain%gindex(c,r)) = var(c,r) 
                endif
             enddo
          enddo

          write(LVT_logunit,*) '[INFO] Reading standard deviation file ',&
               trim(LVT_pctile_struc(m)%std_filename2)
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=LVT_pctile_struc(m)%std_filename2, form='unformatted')
          read(ftn) var
          call LVT_releaseUnitNumber(ftn)

          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(LVT_domain%gindex(c,r).ne.-1) then 
                   LVT_pctile_struc(m)%std2(LVT_domain%gindex(c,r)) = var(c,r) 
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine LVT_initpercentile

!BOP
! 
! !ROUTINE: LVT_diagnosepercentile
! \label{LVT_diagnosepercentile}
!
! !INTERFACE: 
  subroutine LVT_diagnosepercentile(pass)
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
!   calculating the percentile of desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelpercentile](\ref{diagnoseSingleModelpercentile})
!     updates the percentile computation for a single variable. This routine
!     stores the precip values to be used later for computing percentile, 
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

    if(LVT_rc%startmode.eq."coldstart") then 
       if(pass.eq.1) then 
          if(LVT_metrics%percentile%selectOpt.eq.1.or.&
               LVT_metrics%percentile%timeOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model))   
                call diagnoseSinglepercentileParams(model, obs, stats,&
                     LVT_metrics%percentile)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
          endif
       elseif(pass.eq.2) then 
          if(LVT_metrics%percentile%selectOpt.eq.1.or.&
               LVT_metrics%percentile%timeOpt.eq.1) then 
             
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model))   
                
                call diagnoseSinglepercentile(model, obs, stats,&
                     LVT_metrics%percentile)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
          endif
          
       endif
    elseif(LVT_rc%startmode.eq."restart") then 
       if(LVT_metrics%percentile%selectOpt.eq.1.or.&
            LVT_metrics%percentile%timeOpt.eq.1) then 
          
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)
          
          do while(associated(model))   
             
             call diagnoseSinglepercentile(model, obs, stats,&
                  LVT_metrics%percentile)
             
             model => model%next
             obs => obs%next
             stats => stats%next
             
          enddo
       endif
    endif
  end subroutine LVT_diagnosepercentile


!BOP
! 
! !ROUTINE: diagnoseSinglepercentileParams
! \label{diagnoseSinglepercentileParams}
!
! !INTERFACE: 
  subroutine diagnoseSinglepercentileParams(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the percentile computation of the 
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
    integer                 :: t,k,m,tind
    integer                 :: nsize
    type(ESMF_Time)         :: startTime, currTime
    type(ESMF_TimeInterval) :: ts
    character*3             :: ftime
    character*100           :: filename
    integer                 :: tindex, day1,day2, day3, day4, day5
    integer                 :: yr_index
    integer                 :: p1,p2,p3,p4,p5
    logical                 :: skip_flag, leap_year
    integer                 :: ftn
    integer                 :: rc    
    real                    :: min_value
    real                    :: m_value
    real, allocatable       :: pctile_model_day1(:,:,:)
    real, allocatable       :: pctile_model_day2(:,:,:)
    real, allocatable       :: pctile_model_day3(:,:,:)
    real, allocatable       :: pctile_model_day4(:,:,:)
    real, allocatable       :: pctile_model_day5(:,:,:)

    allocate(pctile_model_day1(LVT_rc%ngrid,&
         LVT_rc%nensem,LVT_pctile_struc(1)%nsize_total))
    allocate(pctile_model_day2(LVT_rc%ngrid,&
         LVT_rc%nensem,LVT_pctile_struc(1)%nsize_total))
    allocate(pctile_model_day3(LVT_rc%ngrid,&
         LVT_rc%nensem,LVT_pctile_struc(1)%nsize_total))
    allocate(pctile_model_day4(LVT_rc%ngrid,&
         LVT_rc%nensem,LVT_pctile_struc(1)%nsize_total))
    allocate(pctile_model_day5(LVT_rc%ngrid,&
         LVT_rc%nensem,LVT_pctile_struc(1)%nsize_total))

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then        

       call LVT_getYearIndex(LVT_rc, yr_index, leap_year)

       skip_flag = .false.
       if(.not.leap_year) then 
          tindex = LVT_rc%doy
       else
          if(LVT_rc%doy.gt.60) then 
             tindex = LVT_rc%doy - 1
          elseif(LVT_rc%doy.eq.60) then  !29th
             skip_flag = .true. 
          else
             tindex = LVT_rc%doy
          endif
       endif
       
       if(.not.skip_flag) then 
          if(tindex-2.lt.1) then 
             day1 = 365+tindex-2
          else                      
             day1 = tindex - 2
          endif
          if(tindex-1.lt.1) then 
             day2 = 365+tindex-1
          else                      
             day2 = tindex - 1
          endif
          
          day3 = tindex
          
          if(tindex+1.gt.365) then 
             day4 = tindex+1-365
          else
             day4 = tindex+1
          endif
          
          if(tindex+2.gt.365) then 
             day5 = tindex+2-365
          else
             day5 = tindex+2
          endif
                      
          p1  = (yr_index -1)*LVT_pctile_struc(1)%nwindow +1
          p2  = (yr_index -1)*LVT_pctile_struc(1)%nwindow +2
          p3  = (yr_index -1)*LVT_pctile_struc(1)%nwindow +3
          p4  = (yr_index -1)*LVT_pctile_struc(1)%nwindow +4
          p5  = (yr_index -1)*LVT_pctile_struc(1)%nwindow +5
!open 5 files corresponding to the 5 days. 

          write(unit=ftime,fmt='(I3.3)') day1
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) nsize
          read(ftn) pctile_model_day1
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day2
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) nsize
          read(ftn) pctile_model_day2
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day3
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) nsize
          read(ftn) pctile_model_day3
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day4
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) nsize
          read(ftn) pctile_model_day4
          call LVT_releaseUnitNumber(ftn)
         
          write(unit=ftime,fmt='(I3.3)') day5
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) nsize
          read(ftn) pctile_model_day5
          call LVT_releaseUnitNumber(ftn)

          if(trim(model%short_name).eq."SWE") then 
             min_value = 1 !1mm
          elseif(trim(model%short_name).eq."Qle".or.&
               trim(model%short_name).eq."Evap") then 
             min_value = -50
          elseif(trim(model%short_name).eq."TWS") then 
             min_value = -10000
          else
             min_value = -100.0 
          endif

          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                if(model%count(t,m,1).gt.0) then
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,m,1).ne.LVT_rc%udef.and.&
                           model%value(t,m,1).gt.min_value) then                    
!ignore if the data goes outside the allotted time window
                         call ESMF_TimeSet(currTime, &
                              yy=LVT_rc%yr, mm=LVT_rc%mo, dd=LVT_rc%da,&
                              h =LVT_rc%hr, m =LVT_rc%mn, s =LVT_rc%ss, &
                              calendar=LVT_calendar, rc=rc)
                         call LVT_verify(rc, 'error in ESMF_TimeSet: curr time')
                         
                         if(LVT_pctile_struc(m)%usescaling.eq.1.and.&
                              currTime.ge.LVT_pctile_struc(m)%stime.and.&
                              currTime.le.LVT_pctile_struc(m)%etime) then 
                            m_value = LVT_pctile_struc(m)%mean2(t) + &
                                 (model%value(t,m,1) - LVT_pctile_struc(m)%mean1(t))* &
                                 (LVT_pctile_struc(m)%std2(t)/&
                                 LVT_pctile_struc(m)%std1(t))
                         else
                            m_value = model%value(t,m,1)
                         endif
                         
                         if(p1.le.LVT_pctile_struc(m)%nsize_total) then 
                            pctile_model_day1(t,m,p1)  = & 
                                 m_value
                         endif
                         
                         
                         if(p2.le.LVT_pctile_struc(m)%nsize_total) then 
                            pctile_model_day2(t,m,p2)   = & 
                                 m_value
                         endif
                         
                         if(p3.le.LVT_pctile_struc(m)%nsize_total) then 
                            pctile_model_day3(t,m,p3)   = & 
                                 m_value
                         endif
                         
                         
                         if(p4.le.LVT_pctile_struc(m)%nsize_total) then 
                            pctile_model_day4(t,m,p4) = &
                                 m_value
                         endif
                         
                         if(p5.le.LVT_pctile_struc(m)%nsize_total) then 
                            pctile_model_day5(t,m,p5) = &
                                 m_value
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
!          open(100,file='test.bin',form='unformatted')
!          write(100) tmp1
!          write(100) tmp2
!          close(100)
!          stop

          write(unit=ftime,fmt='(I3.3)') day1
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) LVT_pctile_struc(1)%nsize_total
          write(ftn) pctile_model_day1
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day2
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) LVT_pctile_struc(1)%nsize_total
          write(ftn) pctile_model_day2
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day3
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) LVT_pctile_struc(1)%nsize_total
          write(ftn) pctile_model_day3
          call LVT_releaseUnitNumber(ftn)

          write(unit=ftime,fmt='(I3.3)') day4
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) LVT_pctile_struc(1)%nsize_total
          write(ftn) pctile_model_day4
          call LVT_releaseUnitNumber(ftn)
         
          write(unit=ftime,fmt='(I3.3)') day5
          filename = trim(LVT_rc%statsodir)//'/RST/'//trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          write(ftn) LVT_pctile_struc(1)%nsize_total
          write(ftn) pctile_model_day5
          call LVT_releaseUnitNumber(ftn)
       endif
       
       deallocate(pctile_model_day1)
       deallocate(pctile_model_day2)
       deallocate(pctile_model_day3)
       deallocate(pctile_model_day4)
       deallocate(pctile_model_day5)

    endif
  end subroutine diagnoseSinglepercentileParams

!BOP
! 
! !ROUTINE: diagnoseSinglepercentile
! \label{diagnoseSinglepercentile}
!
! !INTERFACE: 
  subroutine diagnoseSinglepercentile(model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine diagnoses the percentile values for each grid cell. 
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

    integer                 :: t,l,k,m
    integer                 :: nsize
    type(ESMF_Time)         :: currTime
    real                    :: min_value
    character*3             :: ftime
    character*100           :: filename
    integer  :: yr_index, tindex
    logical  :: leap_year, skip_flag
    real     :: pvalue
    real                    :: m_value
    integer                 :: ftn
    integer                 :: rc
    real, allocatable       :: pctile_model(:,:,:)

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then 
       
       allocate(pctile_model(LVT_rc%ngrid,&
            LVT_rc%nensem,&
            LVT_pctile_struc(1)%nsize_total))
       

       call LVT_getYearIndex(LVT_rc, yr_index, leap_year)

       skip_flag = .false.
       if(.not.leap_year) then 
          tindex = LVT_rc%doy
       else
          if(LVT_rc%doy.gt.60) then 
             tindex = LVT_rc%doy - 1
          elseif(LVT_rc%doy.eq.60) then  !29th
             skip_flag = .true. 
          else
             tindex = LVT_rc%doy
          endif
       endif

       if(.not.skip_flag) then 
!read the day's file to compute percentile values. 
          write(unit=ftime,fmt='(I3.3)') tindex
          filename = trim(LVT_rc%statsodir)//'/RST/'//&
               trim(ftime)//'_pctile_climo.bin'
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=trim(filename),form='unformatted')
          read(ftn) nsize
          read(ftn) pctile_model
          call LVT_releaseUnitNumber(ftn)             
          
          if(trim(model%short_name).eq."SWE") then 
             min_value = 1 !1mm
          elseif(trim(model%short_name).eq."Qle".or.&
               trim(model%short_name).eq."Evap") then 
             min_value = -50
          elseif(trim(model%short_name).eq."TWS") then
             min_value = -10000
          else
             min_value = -100.0
          endif
          
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                if(model%value(t,m,1).gt.min_value) then 
                   call ESMF_TimeSet(currTime, &
                        yy=LVT_rc%yr, mm=LVT_rc%mo, dd=LVT_rc%da,&
                        h =LVT_rc%hr, m =LVT_rc%mn, s =LVT_rc%ss, &
                        calendar=LVT_calendar, rc=rc)
                   call LVT_verify(rc, 'error in ESMF_TimeSet: curr time')
                   
                   if(LVT_pctile_struc(m)%usescaling.eq.1.and.&
                        currTime.ge.LVT_pctile_struc(m)%stime.and.&
                        currTime.le.LVT_pctile_struc(m)%etime) then 
                      m_value = LVT_pctile_struc(m)%mean2(t) + &
                           (model%value(t,m,1) - LVT_pctile_struc(m)%mean1(t))* &
                           (LVT_pctile_struc(m)%std2(t)/&
                           LVT_pctile_struc(m)%std1(t))
                   else
                      m_value = model%value(t,m,1)
                   endif
                   call percentile_value(&
                        m_value,&
                        stats%pctile(m)%value_model_nsize(t,tindex),&                
                        pctile_model(t,m,:),&
                        pvalue)
                   stats%pctile(m)%value_model(t,1,1) = & 
                        stats%pctile(m)%value_model(t,1,1) + pvalue
                   stats%pctile(m)%value_count_model(t,1,1) = &
                        stats%pctile(m)%value_count_model(t,1,1) +1
                endif
             enddo
          enddo
       endif
          
       deallocate(pctile_model)
    endif
  end subroutine diagnoseSinglepercentile

!BOP
! 
! !ROUTINE: LVT_computepercentile
! \label{LVT_computepercentile}
!
! !INTERFACE: 
  subroutine LVT_computepercentile(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to compute the percentile values for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelpercentile](\ref{computeSingleModelpercentile})
!     computes the percentile values for a single variable
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
          if(LVT_metrics%percentile%selectOpt.eq.1.or.&
               LVT_metrics%percentile%timeOpt.eq.1) then 
             call LVT_getDataStream1Ptr(model)
             call LVT_getDataStream2Ptr(obs)
             call LVT_getstatsEntryPtr(stats)
             
             do while(associated(model))             
                call computeSinglepercentileparams(alarm,&
                     model,obs,stats,LVT_metrics%percentile)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
          endif
       elseif(pass.eq.2) then        
          
          if(LVT_metrics%percentile%selectOpt.eq.1.or.&
               LVT_metrics%percentile%timeOpt.eq.1) then 
             if(alarm) then 
                if(LVT_metrics%percentile%timeOpt.eq.1.and.&
                     LVT_metrics%percentile%extractTS.eq.1) then 

                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%percentile%ftn_ts_loc(i,m),200,advance='no') &
                                 LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                                 LVT_rc%hr,'',LVT_rc%mn, '' 

                            write(LVT_pctile_struc(m)%ftn_ts_darea(i),&
                                 200,advance='no') &
                                 LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                                 LVT_rc%hr,'',LVT_rc%mn, '' 
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%percentile%ftn_ts_loc(i,1),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, ''
                         write(LVT_pctile_struc(1)%ftn_ts_darea(i),&
                              200,advance='no') &
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
                call computeSinglepercentile(alarm,&
                     model,obs,stats,LVT_metrics%percentile)
                
                model => model%next
                obs => obs%next
                stats => stats%next
                
             enddo
             
             if(alarm) then 
                if(LVT_metrics%percentile%timeOpt.eq.1.and.&
                     LVT_metrics%percentile%extractTS.eq.1) then

                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%percentile%ftn_ts_loc(i,m),fmt='(a1)') ''
                            write(LVT_pctile_struc(m)%ftn_ts_darea(i),fmt='(a1)') ''
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%percentile%ftn_ts_loc(i,1),fmt='(a1)') ''
                         write(LVT_pctile_struc(1)%ftn_ts_darea(i),fmt='(a1)') ''
                      enddo
                   endif
 
                endif
             endif
          endif
       endif
    elseif(LVT_rc%startmode.eq."restart") then 

       if(LVT_metrics%percentile%selectOpt.eq.1.or.&
            LVT_metrics%percentile%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%percentile%timeOpt.eq.1.and.&
                  LVT_metrics%percentile%extractTS.eq.1) then 
                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%percentile%ftn_ts_loc(i,m),200,advance='no') &
                                 LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                                 LVT_rc%hr,'',LVT_rc%mn, '' 

                            write(LVT_pctile_struc(m)%ftn_ts_darea(i),&
                                 200,advance='no') &
                                 LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                                 LVT_rc%hr,'',LVT_rc%mn, '' 
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%percentile%ftn_ts_loc(i,1),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, ''
                         write(LVT_pctile_struc(1)%ftn_ts_darea(i),&
                              200,advance='no') &
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
             call computeSinglepercentile(alarm,&
                  model,obs,stats,LVT_metrics%percentile)
             
             model => model%next
             obs => obs%next
             stats => stats%next
             
          enddo
          
          if(alarm) then 
             if(LVT_metrics%percentile%timeOpt.eq.1.and.&
                  LVT_metrics%percentile%extractTS.eq.1) then 
                   if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                      do m=1,LVT_rc%nensem
                         do i=1,LVT_rc%ntslocs
                            write(LVT_metrics%percentile%ftn_ts_loc(i,m),fmt='(a1)') ''
                            write(LVT_pctile_struc(m)%ftn_ts_darea(i),fmt='(a1)') ''
                         enddo
                      enddo
                   else
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%percentile%ftn_ts_loc(i,1),fmt='(a1)') ''
                         write(LVT_pctile_struc(1)%ftn_ts_darea(i),fmt='(a1)') ''
                      enddo
                   endif

             endif
          endif
       endif
    endif
  end subroutine LVT_computepercentile
  

!BOP
! 
! !ROUTINE: computeSinglepercentileparams
! \label{computeSinglepercentileparams}
!
! !INTERFACE: 
  subroutine computeSinglepercentileparams(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the percentile values for each grid cell. 
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

    character*100           :: filename
    character*3             :: ftime
    integer                 :: ftn
    integer                 :: nsize
    real, allocatable       :: pctile_model(:,:,:)

    integer  :: t,l,k,m

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       allocate(pctile_model(LVT_rc%ngrid,&
            LVT_rc%nensem,&
            LVT_pctile_struc(1)%nsize_total))
    
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do k=1,365
!read each day's file, sort it and write it back. 
             write(unit=ftime,fmt='(I3.3)') k
             filename = trim(LVT_rc%statsodir)//'/RST/'//&
                  trim(ftime)//'_pctile_climo.bin'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=trim(filename),form='unformatted')
             read(ftn) nsize
             read(ftn) pctile_model
             call LVT_releaseUnitNumber(ftn)             

             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem
                   call compress_and_sort(&
                        LVT_pctile_struc(m)%nsize_total,&
                        pctile_model(t,m,:),&
                        stats%pctile(m)%value_model_nsize(t,k))
                enddo
             enddo

             write(unit=ftime,fmt='(I3.3)') k
             filename = trim(LVT_rc%statsodir)//'/RST/'//&
                  trim(ftime)//'_pctile_climo.bin'
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=trim(filename),form='unformatted')
             write(ftn) LVT_pctile_struc(1)%nsize_total
             write(ftn) pctile_model
             call LVT_releaseUnitNumber(ftn)                          
          enddo
          
       endif
       deallocate(pctile_model)
    endif

  end subroutine computeSinglepercentileparams

!BOP
!
! !ROUTINE: compress_and_sort
! \label{compress_and_sort}
!
! !INTERFACE: 
  subroutine compress_and_sort(size_in, value_in, size_out)
! !USES: 
    use LVT_SortingMod, only : LVT_sort
! !ARGUMENTS: 
    integer             :: size_in
    real                :: value_in(size_in)
    integer             :: size_out
! 
! !DESCRIPTION: 
!   This subroutine compresses the input array after removing
!   undefined values and sorts it in the ascending order. 
! 
!EOP

    real                :: value_out(size_in)
    integer             :: i 

!remove undefs, compute size_out

    value_out = LVT_rc%udef

    size_out = 0 
    do i=1,size_in
       if(value_in(i).ne.LVT_rc%udef) then 
          size_out = size_out + 1
          value_out(size_out) = value_in(i)
       endif
    enddo
!sort value_out and set it back to value_in
    call LVT_sort(value_out, size_out)
    
    value_in = value_out
    
  end subroutine compress_and_sort

!BOP
!
! !ROUTINE: compute_sorted_sizes
! \label{compute_sorted_sizes}
!
! !INTERFACE: 
  subroutine compute_sorted_sizes(size_in, value_in, size_out)
! !USES: 
    use LVT_SortingMod, only : LVT_sort
! !ARGUMENTS: 
    integer             :: size_in
    real                :: value_in(size_in)
    integer             :: size_out
! 
! !DESCRIPTION: 
!   This subroutine compresses the input array after removing
!   undefined values and sorts it in the ascending order. 
! 
!EOP

    integer             :: i 

!remove undefs, compute size_out

    size_out = 0 
    do i=1,size_in
       if(value_in(i).ne.LVT_rc%udef) then 
          size_out = size_out + 1
       endif
    enddo
  end subroutine compute_sorted_sizes

!BOP
! 
! !ROUTINE: computeSinglepercentile
! \label{computeSinglepercentile}
!
! !INTERFACE: 
  subroutine computeSinglepercentile(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine computes the percentile values for each grid cell. 
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

    real     :: pctile_model(LVT_rc%ngrid,&
         LVT_pctile_struc(1)%nsize_total)
    integer  :: t,l,k,i,kk,m,tid,stid,ii
    integer  :: sumv(5, LVT_rc%ntslocs)
    integer  :: sumv_ens(5, LVT_rc%ntslocs,LVT_rc%nensem)
    real     :: percentile_darea(5, LVT_rc%ntslocs)
    real     :: percentile_darea_ens(5, LVT_rc%ntslocs,LVT_rc%nensem)
    real,    allocatable :: tavg_value_ts(:,:,:)

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then 

       do t=1,LVT_rc%ngrid
          do ii=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                if(stats%pctile(ii)%value_count_model(t,k,1).gt.0) then 
                   stats%pctile(ii)%value_model(t,k,1) = &
                        stats%pctile(ii)%value_model(t,k,1)/&
                        stats%pctile(ii)%value_count_model(t,k,1)
                   
                   stats%pctile(ii)%tavg_value_model(t,k,1) = & 
                        stats%pctile(ii)%tavg_value_model(t,k,1) + & 
                        stats%pctile(ii)%value_model(t,k,1) 
                   stats%pctile(ii)%tavg_count_value_model(t,k,1) = &
                        stats%pctile(ii)%tavg_count_value_model(t,k,1) + 1
                else
                   stats%pctile(ii)%value_model(t,k,1) = LVT_rc%udef
                endif
             enddo
          enddo
       enddo

       if(alarm) then 
          
          do t=1,LVT_rc%ngrid
             do ii=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   if(stats%pctile(ii)%tavg_count_value_model(t,k,1).gt.0) then 
                      stats%pctile(ii)%tavg_value_model(t,k,1) = & 
                           stats%pctile(ii)%tavg_value_model(t,k,1) /&
                           stats%pctile(ii)%tavg_count_value_model(t,k,1)  
                   else
                      stats%pctile(ii)%tavg_value_model(t,k,1) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo

          if(metric%extractTS.eq.1) then 

             if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                do ii=1,LVT_rc%nensem
                   call LVT_writeTSinfo(metric%ftn_ts_loc(:,ii),&
                        model,&
                        LVT_rc%ngrid,&
                        stats%pctile(ii)%tavg_value_model,&
                        stats%pctile(ii)%tavg_count_value_model)
                enddo
             else
                allocate(tavg_value_ts(LVT_rc%ngrid, &
                        model%vlevels, &
                        LVT_rc%strat_nlevels))
                tavg_value_ts = 0.0
                do ii=1,LVT_rc%nensem
                   do t=1,LVT_rc%ngrid
                      do k=1,model%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            tavg_value_ts(t,k,l) = &
                                 tavg_value_ts(t,k,l) + &
                                 stats%pctile(ii)%tavg_value_model(t,k,l)
                         enddo
                      enddo
                   enddo
                enddo
                
                do t=1,LVT_rc%ngrid
                   do k=1,model%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         tavg_value_ts(t,k,l) = &
                              tavg_value_ts(t,k,l)/LVT_rc%nensem
                      enddo
                   enddo
                enddo
                      
                call LVT_writeTSinfo(metric%ftn_ts_loc(:,1),&
                     model,&
                     LVT_rc%ngrid,&
                     tavg_value_ts,&
                     stats%pctile(1)%tavg_count_value_model)
                deallocate(tavg_value_ts)    

             endif
          endif
       endif

       do ii=1,LVT_rc%nensem
          do k=1,model%selectNlevs
             call LVT_computeCI(stats%pctile(ii)%tavg_value_model(:,k,1),&
                  LVT_rc%ngrid,&
                  LVT_rc%pval_CI, stats%pctile(ii)%value_model_ci(k,1))
          enddo
       enddo
          
!       do k=1,model%selectNlevs
!          call LVT_writeDataBasedStrat(model,obs,stats,metric,&
!               LVT_rc%ngrid, stats%pctile(ii)%tavg_value_model)
!       enddo

       if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
          do ii=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do i=1,LVT_rc%ntslocs
                   sumv_ens(:,i,ii) = 0
                   do kk=1,LVT_TSobj(i)%npts
                      tid = LVT_TSobj(i)%ts_tindex(kk)
                      if(tid.ne.-1) then 
                         t = (tid-1)*1+1
                         if(stats%pctile(ii)%tavg_value_model(t,k,1).ne.LVT_rc%udef) then 
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(1)) then 
                               sumv_ens(1,i,ii) = sumv_ens(1,i,ii) + 1
                            endif
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(2)) then 
                               sumv_ens(2,i,ii) = sumv_ens(2,i,ii) + 1
                            endif
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(3)) then 
                               sumv_ens(3,i,ii) = sumv_ens(3,i,ii) + 1
                            endif
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(4)) then 
                               sumv_ens(4,i,ii) = sumv_ens(4,i,ii) + 1
                            endif
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(5)) then 
                               sumv_ens(5,i,ii) = sumv_ens(5,i,ii) + 1
                            endif
                         endif
                      endif
                   enddo
                   if(LVT_TSobj(i)%npts.gt.0) then
                      percentile_darea_ens(1,i,ii) = real(sumv_ens(1,i,ii))*100.0/&
                           real(LVT_TSobj(i)%npts)
                      percentile_darea_ens(2,i,ii) = real(sumv_ens(2,i,ii))*100.0/&
                           real(LVT_TSobj(i)%npts)
                      percentile_darea_ens(3,i,ii) = real(sumv_ens(3,i,ii))*100.0/&
                           real(LVT_TSobj(i)%npts)
                      percentile_darea_ens(4,i,ii) = real(sumv_ens(4,i,ii))*100.0/&
                           real(LVT_TSobj(i)%npts)
                      percentile_darea_ens(5,i,ii) = real(sumv_ens(5,i,ii))*100.0/&
                           real(LVT_TSobj(i)%npts)
                   else
                      percentile_darea_ens(:,i,ii) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo

       else
          do ii=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do i=1,LVT_rc%ntslocs
                   sumv(:,i) = 0
                   do kk=1,LVT_TSobj(i)%npts
                      tid = LVT_TSobj(i)%ts_tindex(kk)
                      if(tid.ne.-1) then 
                         t = tid
                         if(stats%pctile(ii)%tavg_value_model(t,k,1).ne.LVT_rc%udef) then 
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(1)) then 
                               sumv(1,i) = sumv(1,i) + 1
                            endif
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(2)) then 
                               sumv(2,i) = sumv(2,i) + 1
                            endif
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(3)) then 
                               sumv(3,i) = sumv(3,i) + 1
                            endif
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(4)) then 
                               sumv(4,i) = sumv(4,i) + 1
                            endif
                            if(stats%pctile(ii)%tavg_value_model(t,k,1).lt.&
                                 LVT_pctile_struc(ii)%d_classes(5)) then 
                               sumv(5,i) = sumv(5,i) + 1
                            endif
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo
          do i=1,LVT_rc%ntslocs
             if(LVT_TSobj(i)%npts.gt.0) then
                percentile_darea(1,i) = real(sumv(1,i))*100.0/&
                     real(LVT_TSobj(i)%npts*LVT_rc%nensem)
                percentile_darea(2,i) = real(sumv(2,i))*100.0/&
                     real(LVT_TSobj(i)%npts*LVT_rc%nensem)
                percentile_darea(3,i) = real(sumv(3,i))*100.0/&
                     real(LVT_TSobj(i)%npts*LVT_rc%nensem)
                percentile_darea(4,i) = real(sumv(4,i))*100.0/&
                     real(LVT_TSobj(i)%npts*LVT_rc%nensem)
                percentile_darea(5,i) = real(sumv(5,i))*100.0/&
                     real(LVT_TSobj(i)%npts*LVT_rc%nensem)
             else
                percentile_darea(:,i) = LVT_rc%udef
             endif
          enddo
       endif
       if(alarm) then 
          if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
             do ii=1,LVT_rc%nensem
                do i=1,LVT_rc%ntslocs
                   write(LVT_pctile_struc(ii)%ftn_ts_darea(i),203,advance='no') &
                        (percentile_darea_ens(k,i,ii),k=1,5)
                enddo
             enddo
          else
             do i=1,LVT_rc%ntslocs
                write(LVT_pctile_struc(1)%ftn_ts_darea(i),203,advance='no') &
                     (percentile_darea(k,i),k=1,5)
             enddo
          endif
       endif
    endif

203 format(5E14.6)

  end subroutine computeSinglepercentile
  
  subroutine percentile_value(& 
       value_in, &
       nsize, &
       value_sorted, &
       percentile)

    real    :: value_in
    integer :: nsize
    real    :: value_sorted(nsize)
    real    :: percentile

    integer :: k 
    integer :: nx

    nx = 1
    k = 1
!    do while(k.le.nsize) 
!       if(value_sorted(k).lt.value_in) then 
!          k = k + 1
!       else
!          exit
!       endif
!    enddo

    do while(k.lt.nsize.and.&
         value_sorted(k).lt.value_in)
       k = k + 1
    enddo

    nx = k
    percentile = real(nx)/real(nsize+1)
  end subroutine percentile_value
!BOP
! 
! !ROUTINE: LVT_writeMetric_percentile
! \label(LVT_writeMetric_percentile)
!
! !INTERFACE:
  subroutine LVT_writeMetric_percentile(pass,final,vlevels,stats,obs)
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
!    real                    :: percentile_model(LVT_rc%ngrid,1,1)
!    real                    :: percentile_obs(LVT_rc%ngrid,1,1)
!    integer                 :: count_percentile_model(LVT_rc%ngrid,1,1)
!    integer                 :: count_percentile_obs(LVT_rc%ngrid,1,1)

    real, allocatable       :: value_ts(:,:)
    integer, allocatable    :: count_value_ts(:,:)

!    count_percentile_model = 1
!    count_percentile_obs = 1

    if(pass.eq.LVT_metrics%percentile%npass) then 
       if(final.ne.1) then
           if(stats%selectOpt.eq.1) then 
              allocate(value_ts(LVT_rc%ngrid,LVT_rc%nensem))
              allocate(count_value_ts(LVT_rc%ngrid, LVT_rc%nensem))
              count_value_ts = 1

              if(LVT_metrics%percentile%timeOpt.eq.1) then 
                 do k=1,vlevels
                    do m=1,LVT_rc%nensem
                       value_ts(:,m) = & 
                            stats%pctile(m)%tavg_value_model(:,k,1)
                    enddo
                    call LVT_writevar_gridded(LVT_metrics%percentile%ftn_ts, &
                         value_ts(:,:),&
                         stats%vid_ts(LVT_percentileid,1),k)
                    call LVT_writevar_gridded(LVT_metrics%percentile%ftn_ts, &
                         real(count_value_ts(:,:)),&
                         stats%vid_count_ts(LVT_Percentileid,1),k)
                 enddo
                 deallocate(value_ts)
                 deallocate(count_value_ts)
              endif

           end if
        end if
     end if
#if 0 
              do m=1,LVT_rc%nensem
                 do k=1,vlevels
                    if(LVT_metrics%percentile%timeOpt.eq.1) then 
                       call LVT_writevar_gridded(LVT_metrics%percentile%ftn_ts, &
                            stats%pctile(m)%tavg_value_model(:,k,1),&
                            stats%vid_ts(LVT_percentileid,1),dim1=k)
                       !                        stats%vid_ts(LVT_percentileid,1),1.0,dim1=k)
                    endif
                 enddo
              enddo
          endif
       endif
    endif
#endif

  end subroutine LVT_writeMetric_percentile

!BOP
! 
! !ROUTINE: LVT_resetMetric_percentile
! \label(LVT_resetMetric_percentile)
!
! !INTERFACE:
  subroutine LVT_resetMetric_percentile(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical                   :: alarm
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
       if(stats%selectOpt.eq.1) then 
          do m=1,LVT_rc%nensem
             if(LVT_metrics%percentile%timeOpt.eq.1) then 
                stats%pctile(m)%value_model  = 0
                stats%pctile(m)%value_count_model = 0 
                if(alarm) then
                   stats%pctile(m)%tavg_value_model  = 0
                   stats%pctile(m)%tavg_count_value_model = 0 
                endif
             endif
          enddo
       endif
       
       model => model%next
       obs   => obs%next
       stats => stats%next
       
    enddo

  end subroutine LVT_resetMetric_percentile

!BOP
! 
! !ROUTINE: LVT_writerestart_percentile
! 
! !INTERFACE:
  subroutine LVT_writerestart_percentile(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for percentile metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer                          :: k,l
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

#if 0 
    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)

    do while(associated(model))
       if(LVT_metrics%percentile%selectOpt.eq.1) then 
          do l=1,LVT_pctile_struc(m)%nsize_total
             call LVT_writevar_restart(ftn,&
                  LVT_pctile_struc(m)%model_final(:,l),tileflag=1)
          enddo
          do l=1,LVT_pctile_struc(m)%nasc
             call LVT_writevar_restart(ftn,&
                  LVT_pctile_struc(m)%model_mu(:,l),tileflag=1)
             call LVT_writevar_restart(ftn,&
                  LVT_pctile_struc(m)%model_sigma(:,l),tileflag=1)
          enddo
       end if

       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo
#endif
  end subroutine LVT_writerestart_percentile

!BOP
! 
! !ROUTINE: LVT_readrestart_percentile
! 
! !INTERFACE:
  subroutine LVT_readrestart_percentile(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for percentile metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats
    integer              :: k,l
#if 0 

    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(model%short_name.eq."RootMoist") then 
          if(LVT_metrics%percentile%selectOpt.eq.1) then 
             do l=1,LVT_pctile_struc(m)%nsize_total
                call LVT_readvar_restart(ftn,&
                     LVT_pctile_struc(m)%model_final(:,l),tileflag=1)
             enddo
             do l=1,LVT_pctile_struc(m)%nasc
                call LVT_readvar_restart(ftn,&
                     LVT_pctile_struc(m)%model_mu(:,l),tileflag=1)
                call LVT_readvar_restart(ftn,&
                     LVT_pctile_struc(m)%model_sigma(:,l),tileflag=1)
             enddo
          end if
       end if
       
       model => model%next
       obs   => obs%next
       stats => stats%next
    enddo
#endif
  end subroutine LVT_readrestart_percentile

end module LVT_percentileMod
