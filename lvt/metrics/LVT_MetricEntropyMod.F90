!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
module LVT_MetricEntropyMod
!
!BOP
! !MODULE: LVT_MetricEntropyMod
! 
!  !DESCRIPTION: 
!   This module handles the computation of the information-theory
!   metric called 'metric entropy'. The program works by converting
!   time series of variables to symbolic string sequences (binary
!   sequences). 
!  
!   References: 
!   Shannon, A mathematical theory of communication, 1984
!   Pachepsky et al., Geoderma, 134, 253-266, 2006. 
!
!  !NOTES: 
!  The computation of the metric in time, averge seasonal cycles, 
!  average diurnal cycles are not available for this metric. 
!
!  !REVISION HISTORY: 
!  1 Aug 2011    Sujay Kumar  Initial Specification
!
!EOP
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod

  implicit none

  PRIVATE

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initMetricEntropy
  public :: LVT_diagnoseMetricEntropy
  public :: LVT_computeMetricEntropy 
  public :: LVT_writeMetric_MetricEntropy
  public :: LVT_resetMetric_MetricEntropy
  public :: LVT_writerestart_MetricEntropy
  public :: LVT_readrestart_MetricEntropy

  type, public :: mentropydec
     integer                 :: ngap_threshold
     integer                 :: nts_threshold
     integer, allocatable    :: nts(:,:)
     real, allocatable       :: value_model_ts(:)
     character*100           :: rstFile
  end type mentropydec

  type(mentropydec), allocatable, save :: LVT_mentropy_struc(:)

contains
  
!BOP
! !ROUTINE: LVT_initMetricEntropy 
! \label{LVT_initMetricEntropy}
!
! !INTERFACE: 
  subroutine LVT_initMetricEntropy(selectNlevs,stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This is the initialization routine for the metric entropy calculations. 
!  This routine allocates memory for the required data structures
!  and sets the number of required passes through the data (=1).
!  
!  The subroutine arguments are: 
!  \begin{description}
!    \item[model]
!      object to hold model variable information    
!    \item[obs]
!      object to hold observation variable information
!    \item[stats]
!     object to hold statistics computation related information
!    \item[metric]
!     object to hold metric-specific information
!   \end{description}
!EOP
    integer             :: m

    allocate(LVT_mentropy_struc(LVT_rc%nensem))
    allocate(stats%mentropy(LVT_rc%nensem))
    
    do m=1,LVT_rc%nensem
       
       if(metric%selectOpt.eq.1) then 
          allocate(stats%mentropy(m)%value_model_ts(LVT_rc%ngrid,selectNlevs,&
               LVT_rc%nts))
          allocate(stats%mentropy(m)%value_final(LVT_rc%ngrid,selectNlevs))
          allocate(stats%mentropy(m)%value_ci(selectNlevs))
          
          stats%mentropy(m)%value_model_ts = LVT_rc%udef
          stats%mentropy(m)%value_final = LVT_rc%udef
          stats%mentropy(m)%value_ci = LVT_rc%udef
       endif
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------
    metric%npass = 1    
    metric%obsData = .true. 
    metric%stdevFlag = .false. 

!-------------------------------------------------------------------------
! These options are not supported
!-------------------------------------------------------------------------
    metric%timeOpt = 0 
    metric%extractTS = 0 
    metric%computeSC = 0 
    metric%computeADC = 0 
    LVT_mentropy_struc(:)%rstFile = trim(LVT_rc%statsodir)//'/ITmetric_RST.bin'

    LVT_mentropy_struc(:)%ngap_threshold= 20 !5 days of gap
    LVT_mentropy_struc(:)%nts_threshold= 10 !minimum number of values in the time series
    do m=1,LVT_rc%nensem
       allocate(LVT_mentropy_struc(m)%nts(LVT_rc%ngrid, selectNlevs))
    enddo

  end subroutine LVT_initMetricEntropy

!BOP
! !ROUTINE: LVT_diagnoseMetricEntropy
! \label{LVT_diagnoseMetricEntropy}
! 
! !INTERFACE: 
  subroutine LVT_diagnoseMetricEntropy(pass)

    implicit none
! !ARGUMENTS:
    integer                 :: pass
! !DESCRIPTION:
!   This routine invokes the call to update the metric entropy calculation
!   for each variable. 
! 
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     current pass index number over the data points
!   \end{description}
!EOP

    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%mentropy%selectOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleMetricEntropy(&
                  obs,model,stats,&
                  LVT_metrics%mentropy)

             model => model%next
             obs => obs%next
             stats => stats%next

          enddo
       endif
    endif

  end subroutine LVT_diagnoseMetricEntropy

!BOP
! !ROUTINE: diagnoseSingleMetricEntropy
! \label{diagnoseSingleMetricEntropy}
! 
! !INTERFACE: 
  subroutine diagnoseSingleMetricEntropy(obs,model,stats,metric)
! !USES: 
    use LVT_coreMod,    only : LVT_rc, LVT_domain
    use LVT_timeMgrMod, only : LVT_clock, LVT_calendar
    use LVT_histDataMod 
    use LVT_statsDataMod
    use LVT_logMod,     only : LVT_logunit, LVT_endrun

    implicit none
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: obs
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
! !DESCRIPTION:  
!  This subroutine gathers the relevant information for the computation
!  of metric entropy during the pass through the data. The routine
!  simply stores the variable values into an array for computing
!  the metric at the end of the analysis time. 
!
!  The subroutine arguments are: 
!  \begin{description}
!    \item[model]
!      object to hold model variable information    
!    \item[obs]
!      object to hold observation variable information
!    \item[stats]
!     object to hold statistics computation related information
!    \item[metric]
!     object to hold metric-specific information
!   \end{description}
!EOP 
    type(ESMF_Time)         :: currTime,startTime
    type(ESMF_TimeInterval) :: timestep
    integer                 :: tindex
    integer                 :: status
    integer                 :: t,k,m, m_k,o_k
    
    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                m_k = k+model%startNlevs -1
                o_k = k+obs%startNlevs -1
                if(model%count(t,m,m_k).gt.0) then
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then
                         
                         call ESMF_TimeSet(currTime, yy = LVT_rc%yr, &
                              mm = LVT_rc%mo, &
                              dd = LVT_rc%da, &
                              h  = LVT_rc%hr, &
                              m  = LVT_rc%mn, & 
                              s  = LVT_rc%ss, &
                              calendar = LVT_calendar, & 
                              rc = status)
                         call LVT_verify(status,&
                              'Error in ESMF_TimeSet diagnoseEntropy')
                         call ESMF_ClockGet(LVT_clock,&
                              startTime = starttime,&
                              rc=status)
                         call LVT_verify(status,&
                              'Error in ESMF_TimeGet diagnoseEntropy')
                         call ESMF_TimeIntervalSet(timestep,&
                              s=LVT_rc%tavgInterval,rc=status)
                         call LVT_verify(status,&
                              'Error in ESMF_TimeIntervalSet diagnoseEntropy')
                         tindex = (nint((currTime - starttime)/timestep )+1)-1
                         
                         stats%mentropy(m)%value_model_ts(t,k,tindex) = model%value(t,m,m_k)

                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine diagnoseSingleMetricEntropy

!BOP
! !ROUTINE: LVT_computeMetricEntropy
! \label{LVT_computeMetricEntropy}
! 
! !INTERFACE: 
  subroutine LVT_computeMetricEntropy(pass,alarm)

    implicit none
! !ARGUMENTS:
    integer     :: pass
    logical     :: alarm
!
! !DESCRIPTION: 
! This subroutine invokes the call to compute the metric entropy 
!  values for each variable. 
!
!   The arguments are: 
!   \begin{description}
!    \item[pass]
!     current pass index number over the data points
!    \item[alarm]
!     flag indicating if the temporal output interval has been reached
!   \end{description}
! 
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%mentropy%selectOpt.eq.1.or.&
            LVT_metrics%mentropy%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call computeSingleMetricEntropy(&
                  alarm,model,obs,stats,&
                  LVT_metrics%mentropy)

             model => model%next
             obs => obs%next
             stats => stats%next             

          enddo
       endif
    endif
  end subroutine LVT_computeMetricEntropy

!BOP
! !ROUTINE: computeSingleMetricEntropy
! \label{computeSingleMetricEntropy}
! 
! !INTERFACE: 
  subroutine computeSingleMetricEntropy(alarm,model,obs,stats,metric)
! !USES: 
    use LVT_informationContentMod

    implicit none
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This subroutine performs the computation of the metric entropy values 
!  for each variable. The time series values are converted to a binary 
!  string based on the computed median (All values above the median are 
!  given a value of 1 and all values below the median are given a value
!  of 0). From this binary string, words are defined basec on a given 
!  word length. The metric entropy is calculcated based on 
!  Pachepsky et al. (2006)
! 
!  The subroutine arguments are: 
!  \begin{description}
!    \item[alarm]
!     flag indicating if the temporal output interval has been reached
!    \item[model]
!      object to hold model variable information    
!    \item[obs]
!      object to hold observation variable information
!    \item[stats]
!     object to hold statistics computation related information
!    \item[metric]
!     object to hold metric-specific information
!   \end{description}
!EOP
    integer :: t,k,kk,i,l,iprev,m
    integer :: binval
    real    :: med_val(LVT_rc%ngrid,model%selectNlevs)
    character*1, allocatable :: bitstr(:)
    integer     :: ts_class(LVT_rc%ngrid, model%selectNlevs)
!number of words in the string : N-L+1
    character*1, allocatable :: words(:,:)
!    character*1 :: words(LVT_rc%nts-LVT_ICwordlength+1,LVT_ICwordlength)
    integer     :: nbins(2**LVT_ICwordlength)
    real        :: pL(2**LVT_ICwordlength)
    real        :: mentropy
    character*4 :: fdim2
    character*7 :: fdim1
    integer     :: ftn
    character*100 :: matlab_command

    ts_class = 1

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then

#if (defined USE_MATLAB_SUPPORT)          
! in the cold start mode, write the time series. Else, simply read the gap filled data. 
          if(LVT_rc%startmode.eq."coldstart") then 
          print*, LVT_rc%ngrid, LVT_rc%nts
             ftn = LVT_getNextUnitNumber()
             open(ftn,file=LVT_mentropy_struc(m)%rstFile,form='unformatted',&
                  access='direct',recl=LVT_rc%ngrid*LVT_rc%nts*4)
             write(ftn,rec=1) stats%mentropy(m)%value_model_ts(:,1,:)
             call LVT_releaseUnitNumber(ftn)
             
             write(unit=fdim1, fmt='(i7.7)') LVT_rc%ngrid
             write(unit=fdim2, fmt='(i4.4)') LVT_rc%nts

             matlab_command='matlab -nosplash -nodisplay -nojvm -r "filldata(' & 
                  //trim(fdim1)//','//trim(fdim2)//','&
                  //"'"//trim(LVT_mentropy_struc(m)%rstFile)//"'"//') "'

! Run the Matlab program
             call system(trim(matlab_command))
          endif

          ftn = LVT_getNextUnitNumber()
          open(ftn,file=LVT_mentropy_struc(m)%rstFile,form='unformatted',&
               access='direct',recl=LVT_rc%ngrid*LVT_rc%nts*4)
          read(ftn,rec=1) stats%mentropy(m)%value_model_ts(:,1,:)
          call LVT_releaseUnitNumber(ftn)          
          
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   call trimTimeSeries (t,k,m,stats%mentropy(m)%value_model_ts(t,k,:))
                   
                   if(LVT_mentropy_struc(m)%nts(t,k).le.LVT_mentropy_struc(m)%nts_threshold) then 
                      ts_class(t,k) = 0
                   else
                      do kk=1,LVT_mentropy_struc(m)%nts(t,k)
                         if(LVT_mentropy_struc(m)%value_model_ts(kk).eq.LVT_rc%udef) then 
                            ts_class(t,k) = 0
                         endif
                      enddo
                   endif
                   
                   if(ts_class(t,k).eq.1) then 
                      med_val(t,k) = median(LVT_mentropy_struc(m)%value_model_ts(:),&
                           LVT_mentropy_struc(m)%nts(t,k))
                   else
                      med_val(t,k) = LVT_rc%udef
                   endif
                   
                   if(LVT_mentropy_struc(m)%nts(t,k).gt.LVT_mentropy_struc(m)%nts_threshold) then 
                      allocate(bitstr(LVT_mentropy_struc(m)%nts(t,k)))
                      allocate(words(LVT_mentropy_struc(m)%nts(t,k)-&
                           LVT_ICwordlength+1,LVT_ICwordlength))
                      
                      if(ts_class(t,k).eq.1) then 
                         do l=1,LVT_mentropy_struc(m)%nts(t,k)
                            if(LVT_mentropy_struc(m)%value_model_ts(l).ne.LVT_rc%udef) then 
                               if(LVT_mentropy_struc(m)%value_model_ts(l).gt.med_val(t,k)) then 
                                  bitstr(l) = '1'
                               else
                                  bitstr(l) = '0'
                               endif
                            endif
                         enddo
                      endif
                      
                      if(ts_class(t,k).eq.1) then 
                         nbins = 0 
                         do kk=1,LVT_mentropy_struc(m)%nts(t,k)-LVT_ICwordlength+1
                            do i=1,LVT_ICwordLength    
                               if(kk.gt.1) then 
                                  l = (i+iprev)
                                  if(i.eq.3) iprev = iprev+1
                               else
                                  l= i
                                  iprev=1
                               endif
                               words(kk,i) = bitstr(l)
                            enddo
                            call LVT_findWordIndex(words(kk,:),binval)
                            if(binval.gt.0) &
                                 nbins(binval) = nbins(binval)+1
                         enddo
                         do l=1,2**LVT_ICwordlength
                            pL(l) = (float(nbins(l)))/float(sum(nbins))
                         enddo
                         !metric entropy 
                         mentropy = 0 
                         do l=1,2**LVT_ICwordlength
                            if(pL(l).ne.0) then 
                               mentropy = mentropy -pL(l)*(log(pL(l))/log(2.0))
                            end if
                         enddo
                         
                         stats%mentropy(m)%value_final(t,k) = mentropy/float(LVT_ICwordLength)
                      else
                         stats%mentropy(m)%value_final(t,k) = LVT_rc%udef
                      endif
                      deallocate(bitstr)
                      deallocate(words)
                   else
                      stats%mentropy(m)%value_final(t,k) = LVT_rc%udef
                   endif
                   
                
                enddo
             enddo
          enddo
#else
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   call trimTimeSeries (t,k,m,stats%mentropy(m)%value_model_ts(t,k,:))
                   if(LVT_mentropy_struc(m)%nts(t,k).gt.LVT_mentropy_struc(m)%nts_threshold) then 
                      call fillTimeSeries (t,k,m,LVT_mentropy_struc(m)%value_model_ts)
                      print*, 'completed ',t*100.0/LVT_rc%ngrid, '% of grids'
                   endif
                   
                   if(LVT_mentropy_struc(m)%nts(t,k).le.LVT_mentropy_struc(m)%nts_threshold) then 
                      ts_class(t,k) = 0
                   else
                      do kk=1,LVT_mentropy_struc(m)%nts(t,k)
                         if(LVT_mentropy_struc(m)%value_model_ts(kk).eq.LVT_rc%udef) then 
                            ts_class(t,k) = 0
                         endif
                      enddo
                   endif
                   
                   if(ts_class(t,k).eq.1) then 
                      med_val(t,k) = median(LVT_mentropy_struc(m)%value_model_ts(:),&
                           LVT_mentropy_struc(m)%nts(t,k))
                   else
                      med_val(t,k) = LVT_rc%udef
                   endif
                   
                   if(LVT_mentropy_struc(m)%nts(t,k).gt.LVT_mentropy_struc(m)%nts_threshold) then 
                      allocate(bitstr(LVT_mentropy_struc(m)%nts(t,k)))
                      allocate(words(LVT_mentropy_struc(m)%nts(t,k)-&
                           LVT_ICwordlength+1,LVT_ICwordlength))
                      
                      if(ts_class(t,k).eq.1) then 
                         do l=1,LVT_mentropy_struc(m)%nts(t,k)
                            if(LVT_mentropy_struc(m)%value_model_ts(l).ne.LVT_rc%udef) then 
                               if(LVT_mentropy_struc(m)%value_model_ts(l).gt.med_val(t,k)) then 
                                  bitstr(l) = '1'
                               else
                                  bitstr(l) = '0'
                               endif
                            endif
                         enddo
                      endif
                      
                      if(ts_class(t,k).eq.1) then 
                         nbins = 0 
                         do kk=1,LVT_mentropy_struc(m)%nts(t,k)-LVT_ICwordlength+1
                            do i=1,LVT_ICwordLength    
                               if(kk.gt.1) then 
                                  l = (i+iprev)
                                  if(i.eq.3) iprev = iprev+1
                               else
                                  l= i
                                  iprev=1
                               endif
                               words(kk,i) = bitstr(l)
                            enddo
                            call LVT_findWordIndex(words(kk,:),binval)
                            if(binval.gt.0) &
                                 nbins(binval) = nbins(binval)+1
                         enddo
                         do l=1,2**LVT_ICwordlength
                            pL(l) = (float(nbins(l)))/float(sum(nbins))
                         enddo
                         !metric entropy 
                         mentropy = 0 
                         do l=1,2**LVT_ICwordlength
                            if(pL(l).ne.0) then 
                               mentropy = mentropy -pL(l)*(log(pL(l))/log(2.0))
                            end if
                         enddo
                         
                         stats%mentropy(m)%value_final(t,k) = mentropy/float(LVT_ICwordLength)
                      else
                         stats%mentropy(m)%value_final(t,k) = LVT_rc%udef
                      endif
                      deallocate(bitstr)
                      deallocate(words)
                   else
                      stats%mentropy(m)%value_final(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                call LVT_computeCI(stats%mentropy(m)%value_final(:,k),&
                     LVT_rc%ngrid, LVT_rc%pval_CI,&
                     stats%mentropy(m)%value_ci(k))
             enddo
          enddo
#endif
       endif
    endif
                   
  end subroutine computeSingleMetricEntropy

  subroutine trimTimeSeries( t, l,m, ts_value)

    integer         :: t,l,m
    real            :: ts_value(LVT_rc%nts)
    integer         :: st_index, en_index, k

    st_index = -1
    en_index = -1
    LVT_mentropy_struc(m)%nts(t,l) = -1
    do k=1,LVT_rc%nts
       if(ts_value(k).ne.LVT_rc%udef) then 
          st_index = k
          exit
       endif
    enddo

    do k=LVT_rc%nts, 1, -1
       if(ts_value(k).ne.LVT_rc%udef) then 
          en_index = k
          exit
       endif
    enddo
    if(st_index.gt.0.and.en_index.gt.0) then 
       LVT_mentropy_struc(m)%nts(t,l) = en_index - st_index + 1
    endif

    if(LVT_mentropy_struc(m)%nts(t,l).gt.LVT_mentropy_struc(m)%nts_threshold) then 
       if(allocated(LVT_mentropy_struc(m)%value_model_ts)) then 
          deallocate(LVT_mentropy_struc(m)%value_model_ts)
       endif
       allocate(LVT_mentropy_struc(m)%value_model_ts(LVT_mentropy_struc(m)%nts(t,l)))
       
       do k=st_index, en_index
          LVT_mentropy_struc(m)%value_model_ts(k-st_index+1) =&
               ts_value(k)
       enddo
    endif
  end subroutine trimTimeSeries

  subroutine fillTimeSeries( t, l, m, ts_value)

    integer         :: t,l,m
    real            :: ts_value(LVT_mentropy_struc(m)%nts(t,l))
    integer         :: k, k_new, k_new1, k_new2
    character*100   :: filloption
    integer         :: ftn
    logical         :: fillcheck
    real, allocatable :: temp_data(:)

    filloption = "Linear" !discrete cosine transform
    
!    if(filloption.eq."Linear") then 
       do k=1,LVT_mentropy_struc(m)%nts(t,l)
          if(ts_value(k).eq.LVT_rc%udef) then 
             call findNearestIndex2d(m, t, l, k,ts_value,k_new1, k_new2)
             if(k_new1.gt.0.and.k_new2.gt.0) then 
                if(((k-k_new1).le.LVT_mentropy_struc(m)%ngap_threshold)&
                     .and.((k_new2-k).le.LVT_mentropy_struc(m)%ngap_threshold)) then 
                   ts_value(k) = (ts_value(k_new1)*(k_new2-k) + &
                        ts_value(k_new2)*(k-k_new1))/&
                     (k_new2-k_new1)
                endif
             endif
          endif
       enddo
#if 0 
    elseif(filloption.eq."DCT") then
       fillcheck = .true. 
! first check to see if the gaps are within predefined thresholds. 
       do k=1,LVT_mentropy_struc(m)%nts(t,l)
          if(ts_value(k).eq.LVT_rc%udef) then 
             call findNearestIndex2d(m,t, l, k,ts_value,k_new1, k_new2)
             if(((k-k_new1).gt.LVT_mentropy_struc(m)%ngap_threshold)&
                  .or.((k_new2-k).gt.LVT_mentropy_struc(m)%ngap_threshold)) then 
                fillcheck = .false.
                exit
             endif
          endif
       enddo
       if(fillcheck) then 
          ftn = LVT_getNextUnitNumber()
          open(ftn,file='matlab_input.txt',form='formatted')
          do k=1,LVT_mentropy_struc(m)%nts(t,l)
             if(ts_value(k).eq.LVT_rc%udef) then 
                write(ftn,*) 'NaN'
             else
                write(ftn,fmt='(E14.6)') ts_value(k)
             endif
          enddo
          call LVT_releaseUnitNumber(ftn)

          call system("matlab -nosplash -nodisplay -nojvm < fill1ddata.m");
          ftn = LVT_getNextUnitNumber()
          open(ftn,file='matlab_output.txt',form='formatted')
          allocate(temp_data(LVT_mentropy_struc(m)%nts(t,l)))

          do k=1,LVT_mentropy_struc(m)%nts(t,l)
             read(ftn,*) temp_data(k)
             if(ts_value(k).eq.LVT_rc%udef) then 
                ts_value(k) = temp_data(k)
             endif
          enddo
          deallocate(temp_data)
          call LVT_releaseUnitNumber(ftn)
       endif
#endif
!    endif
  end subroutine fillTimeSeries


  subroutine findNearestIndex2d(m, tile, level, k, ts_value, k_new1, k_new2)
    
    integer         :: m, tile, level
    real            :: ts_value(LVT_mentropy_struc(m)%nts(tile, level))
    integer         :: k, k_new, k_new1, k_new2
    integer         :: t

    k_new1 = -1
    k_new2 = -1
    do t=k,1, -1
       if(ts_value(t).ne.LVT_rc%udef) then 
          k_new1 = t
          exit
       endif
    enddo

    do t=k,LVT_mentropy_struc(m)%nts(tile,level)
       if(ts_value(t).ne.LVT_rc%udef) then 
          k_new2 = t
          exit
       endif
    enddo

  end subroutine findNearestIndex2d


  real function median(x,n)

    use LVT_SortingMod, only : LVT_sort

    IMPLICIT  NONE

    REAL,    DIMENSION(1:), INTENT(IN) :: X
    INTEGER, INTENT(IN)                :: N
    REAL,    DIMENSION(1:N)            :: Temp
    INTEGER                            :: i
    
    DO i = 1, N                       ! make a copy
       Temp(i) = X(i)
    END DO
    CALL  LVT_Sort(Temp, N)               ! sort the copy
    IF (MOD(N,2) == 0) THEN           ! compute the median
       Median = (Temp(N/2) + Temp(N/2+1)) / 2.0
    ELSE
       Median = Temp(N/2+1)
    END IF
  end function median

!BOP
! !ROUTINE: LVT_writeMetric_MetricEntropy
! \label{LVT_writeMetric_MetricEntropy}
!
! !INTERFACE: 
  subroutine LVT_writeMetric_MetricEntropy(pass,final,vlevels,stats,obs)
! !USES:
    use LVT_statsMod, only : LVT_writeSummaryStats2
    use LVT_pluginIndices
! !ARGUMENTS: 
    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs
!
! !DESCRIPTION:
!   This subroutine writes the computed metric entropy values to an 
!   external file
!
!  The subroutine arguments are: 
!  \begin{description}
!    \item[pass]
!     current pass index number over the data points
!    \item[final]
!     integer flag indicating if the end of the analysis period is reached
!    \item[vlevels]
!     number of vertical levels in the current variable
!    \item[obs]
!      object to hold observation variable information
!    \item[stats]
!     object to hold statistics computation related information
!   \end{description}
!EOP
    integer                 :: k,l,m,tind

    integer                 :: count_mentropy_final(LVT_rc%ngrid,1)
    real,    allocatable    :: value_total(:,:)
    real,    allocatable    :: value_ci(:)

    count_mentropy_final = 1

    if(pass.eq.LVT_metrics%mentropy%npass) then
       if(stats%selectOpt.eq.1) then 

          allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
          allocate(value_ci(LVT_rc%nensem))

          do k=1,vlevels
             do m=1,LVT_rc%nensem
                value_total(:,m) = &
                     stats%mentropy(m)%value_final(:,k)
                value_ci(m) = stats%mentropy(m)%value_ci(k)

                if(LVT_metrics%mentropy%selectOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%mentropy%ftn_total, &
                        value_total(:,:),&
                        stats%vid_total(LVT_MENTROPYid,1),k)
                   call LVT_writevar_gridded(LVT_metrics%mentropy%ftn_total, &
                        real(count_mentropy_final(:,:)),&
                        stats%vid_count_total(LVT_MENTROPYid,1),k)
                   
                   call LVT_writeSummaryStats2(&
                        LVT_metrics%mentropy%ftn_summ,&
                        LVT_metrics%mentropy%short_name,&
                        LVT_rc%ngrid,&
                        value_total(:,:), &
                        count_mentropy_final(:,:),&
                        stats%standard_name,&
                        value_ci(:))
                endif
                
             enddo
          enddo
          
          deallocate(value_total)
          deallocate(value_ci)
          
       endif
    endif

  end subroutine LVT_writeMetric_MetricEntropy

!BOP
! 
! !ROUTINE: LVT_resetMetric_MetricEntropy
! \label{LVT_resetMetric_MetricEntropy}
!
! !INTERFACE: 
  subroutine LVT_resetMetric_MetricEntropy
!
! !DESCRIPTION: 
!  This routine resets the relevant variables between each temporal averaging
!  interval. 
! 
!EOP
  end subroutine LVT_resetMetric_MetricEntropy

!BOP
! 
! !ROUTINE: LVT_writerestart_MetricEntropy
! 
! !INTERFACE:
  subroutine LVT_writerestart_MetricEntropy(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for MetricEntropy metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%mentropy%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for MetricEntropy'
!       stop
    end if
    
  end subroutine LVT_writerestart_MetricEntropy

!BOP
! 
! !ROUTINE: LVT_readrestart_MetricEntropy
! 
! !INTERFACE:
  subroutine LVT_readrestart_MetricEntropy(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for MetricEntropy metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%mentropy%selectOpt.eq.1) then 
       
!       print*, 'The readrestart method is not implemented for MetricEntropy'
!       stop
    end if
    
  end subroutine LVT_readrestart_MetricEntropy

end module LVT_MetricEntropyMod
