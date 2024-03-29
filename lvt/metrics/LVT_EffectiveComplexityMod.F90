!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LVT_EffectiveComplexityMod
!
!BOP
! !MODULE: LVT_EffectiveComplexityMod
! 
!  !DESCRIPTION: 
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
  public :: LVT_initEffectiveComplexity
  public :: LVT_diagnoseEffectiveComplexity
  public :: LVT_computeEffectiveComplexity 
  public :: LVT_writeMetric_EffectiveComplexity
  public :: LVT_resetMetric_EffectiveComplexity
  public :: LVT_writerestart_EffectiveComplexity
  public :: LVT_readrestart_EffectiveComplexity

  type, public :: ecomplexitydec
     integer                 :: ngap_threshold
     integer                 :: nts_threshold
     integer, allocatable    :: nts(:,:)
     real, allocatable       :: value_model_ts(:)
  end type ecomplexitydec

  type(ecomplexitydec), allocatable, save :: LVT_ecomplexity_struc(:)

contains
  
  subroutine LVT_initEffectiveComplexity(selectNlevs,stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
    
    integer                 :: m
    
    allocate(LVT_ecomplexity_struc(LVT_rc%nensem))
    allocate(stats%ecomplexity(LVT_rc%nensem))

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%ecomplexity(m)%value_model_ts(LVT_rc%ngrid,selectNlevs(1),&
            LVT_rc%nts))
          allocate(stats%ecomplexity(m)%value_final(LVT_rc%ngrid,selectNlevs(1)))
          allocate(stats%ecomplexity(m)%value_ci(selectNlevs(1)))
          
          stats%ecomplexity(m)%value_model_ts = LVT_rc%udef
          stats%ecomplexity(m)%value_final = LVT_rc%udef
          stats%ecomplexity(m)%value_ci = LVT_rc%udef
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

    LVT_ecomplexity_struc(:)%ngap_threshold= 20 !5 days of gap
    LVT_ecomplexity_struc(:)%nts_threshold= 10 !minimum number of values in the time series
    do m=1,LVT_rc%nensem
       allocate(LVT_ecomplexity_struc(m)%nts(LVT_rc%ngrid, selectNlevs(1)))
    enddo

  end subroutine LVT_initEffectiveComplexity

  subroutine LVT_diagnoseEffectiveComplexity(pass)

    implicit none
    integer                 :: pass
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.1) then 
       if(LVT_metrics%ecomplexity%selectOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))

             call diagnoseSingleEffectiveComplexity(&
                  obs, model, stats, &
                  LVT_metrics%ecomplexity)
             
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif

  end subroutine LVT_diagnoseEffectiveComplexity

  subroutine diagnoseSingleEffectiveComplexity(obs,model,stats,metric)
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
! 
!  while diagnosing, simply gather the soil moisture values into an array
! 
    type(ESMF_Time)         :: currTime,startTime
    type(ESMF_TimeInterval) :: timestep
    integer                 :: tindex
    integer                 :: status
    integer                 :: t,k,m,m_k
    
    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then
       do t=1,LVT_rc%ngrid
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                
                m_k = k+model%startNlevs -1
                
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
                         
                         stats%ecomplexity(m)%value_model_ts(t,k,tindex) = &
                              model%value(t,m,m_k)
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine diagnoseSingleEffectiveComplexity

  subroutine LVT_computeEffectiveComplexity(pass,alarm)
! !ARGUMENTS:
    implicit none

    integer     :: pass
    logical     :: alarm
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats


    if(pass.eq.1) then 
       if(LVT_metrics%ecomplexity%selectOpt.eq.1.or.&
            LVT_metrics%ecomplexity%timeOpt.eq.1) then 
          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call computeSingleEffectiveComplexity(&
                  alarm,model, obs, stats, &
                  LVT_metrics%ecomplexity)

             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_computeEffectiveComplexity

!BOP
! !ROUTINE: computeSingleEffectiveComplexity
! \label{computeSingleEffectiveComplexity}
! 
! !INTERFACE: 
  subroutine computeSingleEffectiveComplexity(alarm,model,obs,stats,metric)
! !USES: 
    use LVT_informationContentMod

    implicit none
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric

    integer :: t,k,kk,i,j,l,iprev,m
    integer :: col,row
    real    :: med_val(LVT_rc%ngrid,model%selectNlevs)
    character*1, allocatable :: bitstr(:)
    integer     :: ts_class(LVT_rc%ngrid, model%selectNlevs)
!number of words in the string : N-L+1
    character*1, allocatable :: words(:,:)
    integer     :: nbinsij(2**LVT_ICwordlength,2**LVT_ICwordlength)
    real        :: pLij(2**LVT_ICwordlength,2**LVT_ICwordlength)
    real        :: pLcij(2**LVT_ICwordlength,2**LVT_ICwordlength)    
    integer     :: nbins(2**LVT_ICwordlength)
    real        :: pL(2**LVT_ICwordlength)
    integer     :: binval
    real        :: ecomplexity

    ts_class = 1

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   call trimTimeSeries (t,k,m,&
                        stats%ecomplexity(m)%value_model_ts(t,k,:))
                   if(LVT_ecomplexity_struc(m)%nts(t,k).gt.LVT_ecomplexity_struc(m)%nts_threshold) then 
                      call fillTimeSeries (t,k,m,&
                           LVT_ecomplexity_struc(m)%value_model_ts)
                   endif
                   
                   if(LVT_ecomplexity_struc(m)%nts(t,k).le.LVT_ecomplexity_struc(m)%nts_threshold) then 
                      ts_class(t,k) = 0
                   else
                      do kk=1,LVT_ecomplexity_struc(m)%nts(t,k)
                         if(LVT_ecomplexity_struc(m)%value_model_ts(kk).eq.LVT_rc%udef) then 
                            ts_class(t,k) = 0
                         endif
                      enddo
                   endif
                   
                   if(ts_class(t,k).eq.1) then 
                      med_val(t,k) = median(LVT_ecomplexity_struc(m)%value_model_ts(:),&
                           LVT_ecomplexity_struc(m)%nts(t,k))
                   else
                      med_val(t,k) = LVT_rc%udef
                   endif
                   
                   if(LVT_ecomplexity_struc(m)%nts(t,k).gt.LVT_ecomplexity_struc(m)%nts_threshold) then 
                      allocate(bitstr(LVT_ecomplexity_struc(m)%nts(t,k)))
                      allocate(words(LVT_ecomplexity_struc(m)%nts(t,k)-&
                           LVT_ICwordlength+1,LVT_ICwordlength))
                      
                      if(ts_class(t,k).eq.1) then 
                         do l=1,LVT_ecomplexity_struc(m)%nts(t,k)
                            if(LVT_ecomplexity_struc(m)%value_model_ts(l).ne.LVT_rc%udef) then 
                               if(LVT_ecomplexity_struc(m)%value_model_ts(l).gt.med_val(t,k)) then 
                                  bitstr(l) = '1'
                               else
                                  bitstr(l) = '0'
                               endif
                            endif
                         enddo
                      endif
                      
                      if(ts_class(t,k).eq.1) then 
                         nbinsij = 0 
                         nbins = 0 
                         do kk=1,LVT_ecomplexity_struc(m)%nts(t,k)-LVT_ICwordlength+1
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
                         enddo
                         do kk=1,LVT_ecomplexity_struc(m)%nts(t,k)-LVT_ICwordlength+1
                            call LVT_findWordIndex(words(kk,:),binval)
                            if(binval.gt.0) &
                                 nbins(binval) = nbins(binval)+1
                         enddo
                         do l=1,2**LVT_ICwordlength
                            pL(l) = (float(nbins(l)))/float(sum(nbins))
                         enddo
                         do kk=1,LVT_ecomplexity_struc(m)%nts(t,k)-LVT_ICwordlength+1
                            if(kk.gt.1) then
                               call LVT_findWordTransitionIndex(&
                                    words(kk,:),words(kk-1,:), col,row)
                               nbinsij(col,row) = nbinsij(col,row) + 1
                            endif
                         enddo
!probability of transition from the ith to the jth word
                         do i=1,(2**LVT_ICwordlength)
                            do j=1,(2**LVT_ICwordlength)
                               if(sum(nbinsij).ne.0) then 
                                  pLij(i,j) = (float(nbinsij(i,j)))/float(sum(nbinsij))
                               else
                                  pLij(i,j) = 0.0
                               endif
                            enddo
                         enddo
                         !conditional probability from the ith to the jth word, given ith word
                         do i=1,(2**LVT_ICwordlength)
                            do j=1,(2**LVT_ICwordlength)
                               if(sum(nbinsij(i,:)).ne.0) then 
                                  pLcij(i,j) = (float(nbinsij(i,j)))/float(sum(nbinsij(i,:)))
                               else
                                  pLcij(i,j) = 0.0
                               endif
                            enddo
                         enddo
                         
                         ecomplexity = 0 
                         do i=1,2**LVT_ICwordlength
                            do j=1,2**LVT_ICwordlength
                               if(pLcij(i,j).ne.0.and.pL(i).ne.0) then 
                                  if(log(pL(i)).ne.0) then 
                                     ecomplexity = ecomplexity +&
                                          pLij(i,j)*(log(pLcij(i,j))/log(pL(i)))
                                  endif
                               end if
                            enddo
                         enddo
                         
                         stats%ecomplexity(m)%value_final(t,k) = ecomplexity
                      else
                         stats%ecomplexity(m)%value_final(t,k) = LVT_rc%udef
                      endif
                      deallocate(bitstr)
                      deallocate(words)
                   else
                      stats%ecomplexity(m)%value_final(t,k) = LVT_rc%udef
                   endif
                enddo
             enddo
          enddo
          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                call LVT_computeCI(stats%ecomplexity(m)%value_final(:,k),&
                     LVT_rc%ngrid, LVT_rc%pval_CI,&
                     stats%ecomplexity(m)%value_ci(k))
             enddo
          enddo
       endif
    endif
                   
  end subroutine computeSingleEffectiveComplexity

  subroutine trimTimeSeries( t, l,m,ts_value)

    integer         :: t,l,m
    real            :: ts_value(LVT_rc%nts)
    integer         :: st_index, en_index, k

    st_index = -1
    en_index = -1
    LVT_ecomplexity_struc(m)%nts(t,l) = -1
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
       LVT_ecomplexity_struc(m)%nts(t,l) = en_index - st_index + 1
    endif

    if(LVT_ecomplexity_struc(m)%nts(t,l).gt.LVT_ecomplexity_struc(m)%nts_threshold) then 
       if(allocated(LVT_ecomplexity_struc(m)%value_model_ts)) then 
          deallocate(LVT_ecomplexity_struc(m)%value_model_ts)
       endif
       allocate(LVT_ecomplexity_struc(m)%value_model_ts(LVT_ecomplexity_struc(m)%nts(t,l)))
       
       do k=st_index, en_index
          LVT_ecomplexity_struc(m)%value_model_ts(k-st_index+1) =&
               ts_value(k)
       enddo
    endif
  end subroutine trimTimeSeries

  subroutine fillTimeSeries( t, l, m,ts_value)

    integer         :: t,l,m
    real            :: ts_value(LVT_ecomplexity_struc(m)%nts(t,l))
    integer         :: k, k_new, k_new1, k_new2

    do k=1,LVT_ecomplexity_struc(m)%nts(t,l)
       if(ts_value(k).eq.LVT_rc%udef) then 
          call findNearestIndex2d(m,t, l, k,ts_value,k_new1, k_new2)
          if(k_new1.gt.0.and.k_new2.gt.0) then 
             if(((k-k_new1).le.LVT_ecomplexity_struc(m)%ngap_threshold)&
                  .and.((k_new2-k).le.LVT_ecomplexity_struc(m)%ngap_threshold)) then 
                ts_value(k) = (ts_value(k_new1)*(k_new2-k) + &
                     ts_value(k_new2)*(k-k_new1))/&
                     (k_new2-k_new1)
             endif
          endif
       endif
    enddo
    
  end subroutine fillTimeSeries


  subroutine findNearestIndex2d(m,tile, level, k, ts_value, k_new1, k_new2)
    
    integer         :: m
    integer         :: tile, level
    real            :: ts_value(LVT_ecomplexity_struc(m)%nts(tile, level))
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

    do t=k,LVT_ecomplexity_struc(m)%nts(tile,level)
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


  subroutine LVT_writeMetric_EffectiveComplexity(pass,final,vlevels,stats,obs)
! !USES:
    use LVT_statsMod, only : LVT_writeSummaryStats2
    use LVT_pluginIndices
! !ARGUMENTS: 
    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_statsEntry)   :: stats
    type(LVT_metaDataEntry) :: obs
!EOP
    integer                 :: k,l,m,tind

    integer                 :: count_ecomplexity_final(LVT_rc%ngrid,1)
    real,    allocatable    :: value_total(:,:)
    real,    allocatable    :: value_ci(:)

    count_ecomplexity_final = 1

    if(pass.eq.LVT_metrics%ecomplexity%npass) then
       if(stats%selectOpt.eq.1) then 

          allocate(value_total(LVT_rc%ngrid,LVT_rc%nensem))
          allocate(value_ci(LVT_rc%nensem))

          do k=1,vlevels
             do m=1,LVT_rc%nensem
                value_total(:,m) = &
                     stats%ecomplexity(m)%value_final(:,k)
                value_ci(m) = stats%ecomplexity(m)%value_ci(k)

                if(LVT_metrics%ecomplexity%selectOpt.eq.1) then 
                   call LVT_writevar_gridded(LVT_metrics%ecomplexity%ftn_total, &
                        value_total(:,:),&
                        stats%vid_total(LVT_ECOMPLEXITYid,1),k)
                   call LVT_writevar_gridded(LVT_metrics%ecomplexity%ftn_total, &
                        real(count_ecomplexity_final(:,:)),&
                        stats%vid_count_total(LVT_ECOMPLEXITYid,1),k)
                   
                   call LVT_writeSummaryStats2(&
                        LVT_metrics%ecomplexity%ftn_summ,&
                        LVT_metrics%ecomplexity%short_name,&
                        LVT_rc%ngrid,&
                        value_total(:,:), &
                        count_ecomplexity_final(:,:),&
                        stats%standard_name,&
                        value_ci(:))
                endif
                
             enddo
          enddo
          
          deallocate(value_total)
          deallocate(value_ci)
          
       endif
    endif

  end subroutine LVT_writeMetric_EffectiveComplexity

  subroutine LVT_resetMetric_EffectiveComplexity

  end subroutine LVT_resetMetric_EffectiveComplexity


!BOP
! 
! !ROUTINE: LVT_writerestart_EffectiveComplexity
! 
! !INTERFACE:
  subroutine LVT_writerestart_EffectiveComplexity(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for EffectiveComplexity metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ecomplexity%selectOpt.eq.1) then 
       
       print*, 'The writerestart method is not implemented for EffectiveComplexity'

    end if
    
  end subroutine LVT_writerestart_EffectiveComplexity

!BOP
! 
! !ROUTINE: LVT_readrestart_EffectiveComplexity
! 
! !INTERFACE:
  subroutine LVT_readrestart_EffectiveComplexity(ftn)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for EffectiveComplexity metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    if(LVT_metrics%ecomplexity%selectOpt.eq.1) then 
       
       print*, 'The readrestart method is not implemented for EffectiveComplexity'
       stop
    end if
    
  end subroutine LVT_readrestart_EffectiveComplexity

end module LVT_EffectiveComplexityMod
