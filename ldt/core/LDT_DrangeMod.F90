!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_DrangeMod
!
!BOP
! !MODULE: LDT_DrangeMod
! 
!  !DESCRIPTION: 
!   This module handles the dynamic range computations 
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
!  16 Feb 2022   Mahdi Navari; modified to save stratify CDF
!
!EOP
  use LDT_DAobsDataMod

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_diagnoseDrange
  public :: LDT_computeDrange

  private
  
contains

!BOP
! !ROUTINE: LDT_diagnoseDrange
! \label{LDT_diagnoseDrange}
!
! !INTERFACE: 
  subroutine LDT_diagnoseDrange(n)
! !USES:     
    use LDT_coreMod,             only : LDT_rc
    use LDT_DAmetricsDataMod,        only : LDT_DAmetricsPtr
! 
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Drange calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleDrange](\ref{diagnoseSingleDrange})
!     updates the Drange computation for a single variable 
!   \end{description}
! 
!EOP
    implicit none
    integer, intent(in) :: n
    integer :: i, index

    do index=1,LDT_DA_MOC_COUNT
       call diagnoseSingleDrange(n,LDT_DAobsDataPtr(n,index)%dataEntryPtr,&
            LDT_DAmetricsPtr(index)%dataEntryPtr)
    enddo

  end subroutine LDT_diagnoseDrange

!BOP
! !ROUTINE: diagnoseSingleDrange
! \label{diagnoseSingleDrange}
! 
! !INTERFACE: 
  subroutine diagnoseSingleDrange(n,obs, metrics)
! !USES: 
    use LDT_coreMod
    use LDT_DAmetricsDataMod
! 
! !DESCRIPTION: 
!  This routine updates the Drange computation (updates the running 
!  sum calculations of the squared error) 
!  The arguments are: 
!
!  \begin{description}
!   \item[obs] observation object
!   \item[model] model variable object
!   \item[metrics] object to hold the updated statistics
!  \end{description}
!EOP

    implicit none
    integer, intent(in)     :: n
    type(LDT_DAmetaDataEntry) :: obs
    type(DAmetricsEntry) :: metrics

    integer    :: t,j,k, c,r,c1,r1,t1
    integer    :: r_min, r_max
    integer    :: c_min, c_max

    if(LDT_rc%cdf_ntimes.eq.12) then 
       j = LDT_rc%mo
    elseif(LDT_rc%cdf_ntimes.eq.1) then 
       j = 1
    endif

    if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
       do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
             
             t = LDT_domain(n)%gindex(c,r)
             if(t.ge.-1) then 
                r_min = max(r-LDT_rc%sp_sample_cdf_rad,1)
                c_min = max(c-LDT_rc%sp_sample_cdf_rad,1)
                r_max = min(r+LDT_rc%sp_sample_cdf_rad,LDT_rc%lnr(n))
                c_max = min(c+LDT_rc%sp_sample_cdf_rad,LDT_rc%lnc(n))
                
                do r1=r_min, r_max
                   do c1=c_min, c_max
                      
                      t1 = LDT_domain(n)%gindex(c1,r1)
                      if(t1.ne.-1) then 
                         do k=1,obs%vlevels
                            
                            if(obs%count(t1,k).ne.0) then 
                               if(obs%value(t1,k).gt.metrics%maxval(t,j,k)) then 
                                  metrics%maxval(t,j,k) = obs%value(t1,k)
                               endif
                               if(obs%value(t1,k).lt.metrics%minval(t,j,k)) then 
                                  metrics%minval(t,j,k) = obs%value(t1,k)
                               endif
                               metrics%count_drange_total(t,j,k) = &
                                    metrics%count_drange_total(t,j,k) + 1
                               
                            endif
                         enddo
                         
                      endif
                   enddo
                enddo
             endif
          enddo
       enddo
             
    endif
    
  end subroutine diagnoseSingleDrange


!BOP
! 
! !ROUTINE: LDT_computeDrange
! \label{LDT_computeDrange}
! 
! !INTERFACE: 
  subroutine LDT_computeDrange(n)
! !USES: 
    use LDT_coreMod,  only : LDT_rc
    use LDT_DAobsDataMod, only : LDT_DAobsDataPtr
    use LDT_DAmetricsDataMod,        only : LDT_DAmetricsPtr

! 
! !DESCRIPTION: 
!   This subroutine issues the calls to compute Drange values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleDrange](\ref{computeSingleDrange})
!     updates the Drange computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!   \end{description}
!EOP
    implicit none

    integer, intent(in)     :: n
    integer     :: i, index

    do index=1,LDT_DA_MOC_COUNT
       call computeSingleDrange(n,LDT_DAobsDataPtr(n,index)%dataEntryPtr,&
            LDT_DAmetricsPtr(index)%dataEntryPtr)
    enddo

  end subroutine LDT_ComputeDrange

!BOP
! 
! !ROUTINE: computeSingleDrange
! \label{computeSingleDrange}
! 
! !INTERFACE: 
  subroutine computeSingleDrange(n,obs, metrics)
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_DAmetricsDataMod
! 
! !DESCRIPTION: 
!  This routine computes the Drange values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[obs] observation object
!    \item[metrics] object to hold the updated statistics
!  \end{description}
!EOP

    implicit none
    integer, intent(in)     :: n 
    type(LDT_DAmetaDataEntry) :: obs
    type(DAmetricsEntry) :: metrics

    integer  :: t,i,j,k,c,r
    integer  :: sindex,sindex0,sindex1
    real, allocatable :: strat_xrange(:,:,:,:)
    real, allocatable :: strat_delta(:,:,:)
    real, allocatable :: strat_mask(:,:,:)
    real, allocatable :: strat_minval(:,:,:)
    real, allocatable :: strat_maxval(:,:,:)
    real, allocatable :: strat_drange_total(:,:,:)


    if(LDT_rc%endtime.eq.1) then 
       
       if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
          if(LDT_rc%group_cdfs.eq.0 .and. LDT_rc%strat_cdfs.eq.0) then 
             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(metrics%count_drange_total(t,j,k).le.&
                           LDT_rc%obsCountThreshold) then 
                         metrics%maxval(t,j,k) = LDT_rc%udef
                         metrics%minval(t,j,k) = LDT_rc%udef
                         metrics%mask(t,j,k) = LDT_rc%udef
                      else
                         metrics%mask(t,j,k) = metrics%count_drange_total(t,j,k)
                         metrics%delta(t,j,k) = &
                              (metrics%maxval(t,j,k) - metrics%minval(t,j,k))/&
                              (LDT_rc%cdf_nbins-1)
                         metrics%xrange(t,j,k,1) = metrics%minval(t,j,k)

                         do i=2, LDT_rc%cdf_nbins
                            metrics%xrange(t,j,k,i) = &
                                 metrics%xrange(t,j,k,i-1) + & 
                                 metrics%delta(t,j,k)
                         enddo
                      endif
                   enddo
                enddo
             enddo
          elseif(LDT_rc%group_cdfs.eq.1 .and. LDT_rc%strat_cdfs.eq.0) then 
             allocate(strat_drange_total(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_maxval(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_minval(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_mask(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_delta(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_xrange(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels, &
                  LDT_rc%cdf_nbins))

             strat_drange_total = 0 
             strat_maxval = -1000000.0
             strat_minval = 1000000.0
             strat_mask = 0 
             strat_delta = 0 

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex = LDT_rc%cdf_strat_data(t)
                      strat_drange_total(sindex,j,k) = &
                           strat_drange_total(sindex,j,k) + & 
                           metrics%count_drange_total(t,j,k)
                      if(metrics%maxval(t,j,k).gt.&
                           strat_maxval(sindex,j,k)) then 
                         strat_maxval(sindex,j,k) = metrics%maxval(t,j,k) 
                      endif
                      if(metrics%minval(t,j,k).lt.&
                           strat_minval(sindex,j,k)) then 
                         strat_minval(sindex,j,k) = metrics%minval(t,j,k) 
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%group_cdfs_nbins
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(strat_drange_total(t,j,k).le.&
                           LDT_rc%obsCountThreshold) then 
                         strat_maxval(t,j,k) = LDT_rc%udef
                         strat_minval(t,j,k) = LDT_rc%udef
                         strat_mask(t,j,k) = LDT_rc%udef
                      else
                         strat_mask(t,j,k) = strat_drange_total(t,j,k)
                         strat_delta(t,j,k) = &
                              (strat_maxval(t,j,k) - &
                              strat_minval(t,j,k))/&
                              (LDT_rc%cdf_nbins-1)
                         strat_xrange(t,j,k,1) = strat_minval(t,j,k)
                         do i=2, LDT_rc%cdf_nbins                            
                            strat_xrange(t,j,k,i) = &
                                 strat_xrange(t,j,k,i-1) + & 
                                 strat_delta(t,j,k)
                                 
                         enddo
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex = LDT_rc%cdf_strat_data(t)
                      if(strat_mask(sindex,j,k).eq.LDT_rc%udef) then 
                         metrics%maxval(t,j,k) = LDT_rc%udef
                         metrics%minval(t,j,k) = LDT_rc%udef
                         metrics%mask(t,j,k) = LDT_rc%udef
                      else
                         metrics%mask(t,j,k) = strat_mask(sindex,j,k)
                         metrics%delta(t,j,k) = strat_delta(sindex,j,k)
                         metrics%xrange(t,j,k,:) =strat_xrange(sindex,j,k,:)
                         
                         metrics%maxval(t,j,k) = strat_maxval(sindex,j,k)
                         metrics%minval(t,j,k) = strat_minval(sindex,j,k)
                      endif
                   enddo
                enddo
             enddo
             metrics%strat_xrange = strat_xrange
             deallocate(strat_drange_total)
             deallocate(strat_maxval)
             deallocate(strat_minval)
             deallocate(strat_mask)
             deallocate(strat_delta)

!MN: Startification based on monthly precipitation climatology 
!    monthly total precipitation climatology for each pixel are stored in LDT_rc%stratification_data

          elseif(LDT_rc%group_cdfs.eq.0 .and. LDT_rc%strat_cdfs.eq.1) then
             allocate(strat_drange_total(LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_maxval(LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_minval(LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_mask(LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_delta(LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_xrange(LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels, &
                  LDT_rc%cdf_nbins))

             strat_drange_total = 0
             strat_maxval = -1000000.0
             strat_minval = 1000000.0
             strat_mask = 0
             strat_delta = 0
             strat_xrange = 0 

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex = LDT_rc%stratification_data(t,j)
                      strat_drange_total(sindex,j,k) = &
                           strat_drange_total(sindex,j,k) + &
                           metrics%count_drange_total(t,j,k)
                      if(metrics%maxval(t,j,k).gt.&
                           strat_maxval(sindex,j,k)) then
                           strat_maxval(sindex,j,k) = metrics%maxval(t,j,k)
                      endif
                      if(metrics%minval(t,j,k).lt.&
                           strat_minval(sindex,j,k)) then
                           strat_minval(sindex,j,k) = metrics%minval(t,j,k)
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%strat_cdfs_nbins
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(strat_drange_total(t,j,k).le.&
                           LDT_rc%obsCountThreshold) then
                         strat_maxval(t,j,k) = LDT_rc%udef
                         strat_minval(t,j,k) = LDT_rc%udef
                         strat_mask(t,j,k) = LDT_rc%udef
                      else
                         strat_mask(t,j,k) = strat_drange_total(t,j,k)
                         strat_delta(t,j,k) = &
                              (strat_maxval(t,j,k) - &
                              strat_minval(t,j,k))/&
                              (LDT_rc%cdf_nbins-1)
                         strat_xrange(t,j,k,1) = strat_minval(t,j,k)
                         do i=2, LDT_rc%cdf_nbins
                            strat_xrange(t,j,k,i) = &
                                 strat_xrange(t,j,k,i-1) + &
                                 strat_delta(t,j,k)

                         enddo
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex = LDT_rc%stratification_data(t,j)
                      if(strat_mask(sindex,j,k).eq.LDT_rc%udef) then
                         metrics%maxval(t,j,k) = LDT_rc%udef
                         metrics%minval(t,j,k) = LDT_rc%udef
                         metrics%mask(t,j,k) = LDT_rc%udef
                      else
                         metrics%mask(t,j,k) = strat_mask(sindex,j,k)
                         metrics%delta(t,j,k) = strat_delta(sindex,j,k)
                         metrics%xrange(t,j,k,:) =strat_xrange(sindex,j,k,:)

                         metrics%maxval(t,j,k) = strat_maxval(sindex,j,k)
                         metrics%minval(t,j,k) = strat_minval(sindex,j,k)
                      endif
                   enddo
                enddo
             enddo
             metrics%strat_xrange = strat_xrange
             deallocate(strat_drange_total)
             deallocate(strat_maxval)
             deallocate(strat_minval)
             deallocate(strat_mask)
             deallocate(strat_delta)
             deallocate(strat_xrange)

          elseif(LDT_rc%group_cdfs.eq.1 .and. LDT_rc%strat_cdfs.eq.1) then
             allocate(strat_drange_total(LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_maxval(LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_minval(LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_mask(LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_delta(LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_xrange(LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels, &
                  LDT_rc%cdf_nbins))

             strat_drange_total = 0
             strat_maxval = -1000000.0
             strat_minval = 1000000.0
             strat_mask = 0
             strat_delta = 0
             strat_xrange = 0

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex0  = LDT_rc%stratification_data(t,j)
                      sindex1 = LDT_rc%cdf_strat_data(t)
                      sindex = sindex0 + (sindex1 - 1)*LDT_rc%strat_cdfs_nbins
                      strat_drange_total(sindex,j,k) = &
                           strat_drange_total(sindex,j,k) + &
                           metrics%count_drange_total(t,j,k)
                      if(metrics%maxval(t,j,k).gt.&
                           strat_maxval(sindex,j,k)) then
                         strat_maxval(sindex,j,k) = metrics%maxval(t,j,k)
                      endif
                      if(metrics%minval(t,j,k).lt.&
                           strat_minval(sindex,j,k)) then
                         strat_minval(sindex,j,k) = metrics%minval(t,j,k)
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(strat_drange_total(t,j,k).le.&
                           LDT_rc%obsCountThreshold) then
                         strat_maxval(t,j,k) = LDT_rc%udef
                         strat_minval(t,j,k) = LDT_rc%udef
                         strat_mask(t,j,k) = LDT_rc%udef
                      else
                         strat_mask(t,j,k) = strat_drange_total(t,j,k)
                         strat_delta(t,j,k) = &
                              (strat_maxval(t,j,k) - &
                              strat_minval(t,j,k))/&
                              (LDT_rc%cdf_nbins-1)
                         strat_xrange(t,j,k,1) = strat_minval(t,j,k)
                         do i=2, LDT_rc%cdf_nbins
                            strat_xrange(t,j,k,i) = &
                                 strat_xrange(t,j,k,i-1) + &
                                 strat_delta(t,j,k)

                         enddo
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex0  = LDT_rc%stratification_data(t,j)
                      sindex1 = LDT_rc%cdf_strat_data(t)
                      sindex = sindex0 + (sindex1 - 1)*LDT_rc%strat_cdfs_nbins
                      if(strat_mask(sindex,j,k).eq.LDT_rc%udef) then
                         metrics%maxval(t,j,k) = LDT_rc%udef
                         metrics%minval(t,j,k) = LDT_rc%udef
                         metrics%mask(t,j,k) = LDT_rc%udef
                      else
                         metrics%mask(t,j,k) = strat_mask(sindex,j,k)
                         metrics%delta(t,j,k) = strat_delta(sindex,j,k)
                         metrics%xrange(t,j,k,:) =strat_xrange(sindex,j,k,:)

                         metrics%maxval(t,j,k) = strat_maxval(sindex,j,k)
                         metrics%minval(t,j,k) = strat_minval(sindex,j,k)
                      endif
                   enddo
                enddo
             enddo
             metrics%strat_xrange = strat_xrange
             deallocate(strat_drange_total)
             deallocate(strat_maxval)
             deallocate(strat_minval)
             deallocate(strat_mask)
             deallocate(strat_delta)
             deallocate(strat_xrange)

          endif
       endif
    endif

  end subroutine computeSingleDrange
end module LDT_DrangeMod
