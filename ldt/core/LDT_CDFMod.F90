!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_CDFMod
!
!BOP
! !MODULE: LDT_CDFMod
! 
!  !DESCRIPTION: 
!   This module handles the dynamic range computations 
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
!  2 Dec 2021:   Mahdi Navari; modified to stratify CDF based on precipitation
!                              and save the stratified CDF  
!
!EOP
  use LDT_DAobsDataMod

!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_diagnoseCDF
  public :: LDT_computeCDF

  private
  
contains

!BOP
! !ROUTINE: LDT_diagnoseCDF
! \label{LDT_diagnoseCDF}
!
! !INTERFACE: 
  subroutine LDT_diagnoseCDF(n)
! !USES:     
    use LDT_coreMod
    use LDT_DAmetricsDataMod,        only : LDT_DAmetricsPtr
! 
! !DESCRIPTION: 
!   This subroutine issues the calls to update the CDF calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleCDF](\ref{diagnoseSingleCDF})
!     updates the CDF computation for a single variable 
!   \end{description}
! 
!EOP
    implicit none
    integer, intent(in) :: n 
    integer :: i, index 

    do index=1,LDT_DA_MOC_COUNT
       call diagnoseSingleCDF(n,LDT_DAobsDataPtr(n,index)%dataEntryPtr,&
            LDT_DAmetricsPtr(index)%dataEntryPtr)
    enddo

  end subroutine LDT_diagnoseCDF

!BOP
! !ROUTINE: diagnoseSingleCDF
! \label{diagnoseSingleCDF}
! 
! !INTERFACE: 
  subroutine diagnoseSingleCDF(n, obs, metrics)
! !USES: 
    use LDT_coreMod
    use LDT_DAmetricsDataMod
! 
! !DESCRIPTION: 
!  This routine updates the CDF computation (updates the running 
!  sum calculations of the squared error). If spatial sampling is 
!  enabled, then all grid points within a specified radius is 
!  used in the CDF calculations. 
! 
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

    integer    :: t,j,k,c,r,c1,r1,t1
    integer    :: binval
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
             
             if(t.ne.-1) then 
                r_min = max(r-LDT_rc%sp_sample_cdf_rad,1)
                c_min = max(c-LDT_rc%sp_sample_cdf_rad,1)
                r_max = min(r+LDT_rc%sp_sample_cdf_rad,LDT_rc%lnr(n))
                c_max = min(c+LDT_rc%sp_sample_cdf_rad,LDT_rc%lnc(n))

                do r1=r_min, r_max
                   do c1=c_min, c_max             
                   
                      t1 = LDT_domain(n)%gindex(c1,r1)

                      if(t1.ne.-1) then 
                         
                         do k=1,obs%vlevels

                            if(obs%count(t1,k).ne.0.and.&
                                 (.not.(metrics%maxval(t1,j,k).eq.&
                                 metrics%minval(t1,j,k)))) then 

                               if(metrics%delta(t,j,k).gt.0) then 
                                  binval = nint((obs%value(t1,k) - metrics%minval(t,j,k))/&
                                       metrics%delta(t,j,k)) + 1                
                                  metrics%cdf_bincounts(t,j,k,binval) = &
                                       metrics%cdf_bincounts(t,j,k,binval) + 1
                               endif
                            endif
                         enddo
                      endif
                   enddo
                enddo
             endif
          enddo
       enddo
    endif
    
  end subroutine diagnoseSingleCDF


!BOP
! 
! !ROUTINE: LDT_computeCDF
! \label{LDT_computeCDF}
! 
! !INTERFACE: 
  subroutine LDT_computeCDF(n)
! !USES: 
    use LDT_coreMod,  only : LDT_rc
    use LDT_DAmetricsDataMod,        only : LDT_DAmetricsPtr
    use LDT_logMod,                only : LDT_logunit

! 
! !DESCRIPTION: 
!   This subroutine issues the calls to compute CDF values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleCDF](\ref{computeSingleCDF})
!     updates the CDF computation for a single variable 
!   \end{description}
! 
!   The arguments are: 
!   \begin{description}
!   \end{description}
!EOP
    implicit none
    integer, intent(in) :: n 
    integer     :: i, index

    do index=1,LDT_DA_MOC_COUNT
       call computeSingleCDF(n, LDT_DAobsDataPtr(n,index)%dataEntryPtr,&
            LDT_DAmetricsPtr(index)%dataEntryPtr)
    enddo

  end subroutine LDT_ComputeCDF

!BOP
! 
! !ROUTINE: computeSingleCDF
! \label{computeSingleCDF}
! 
! !INTERFACE: 
  subroutine computeSingleCDF(n,obs, metrics)
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_DAmetricsDataMod
! 
! !DESCRIPTION: 
!  This routine computes the CDF values for a single variable
!  The arguments are: 
!
!  \begin{description}
!    \item[obs] observation object
!    \item[metrics] object to hold the updated statistics
!  \end{description}
!EOP

    implicit none
    integer, intent(in) :: n 
    type(LDT_DAmetaDataEntry) :: obs
    type(DAmetricsEntry) :: metrics

    integer              :: t,i,j,k
    integer              :: c,r
    integer              :: sindex,sindex0,sindex1
    integer, allocatable :: strat_bincounts(:,:,:,:)
    real,    allocatable :: strat_cdf(:,:,:,:)
!    real                 :: datamask(LDT_rc%ngrid(n))

    if(LDT_rc%endtime.eq.1) then 
       if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
          if(LDT_rc%group_cdfs.eq.0 .and. LDT_rc%strat_cdfs.eq.0) then 
             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(sum(metrics%cdf_bincounts(t,j,k,:)).ne.0) then 
                         metrics%cdf(t,j,k,1) = metrics%cdf_bincounts(t,j,k,1)/&
                              float(sum(metrics%cdf_bincounts(t,j,k,:)))
                         do i=2,LDT_rc%cdf_nbins
                            metrics%cdf(t,j,k,i) = metrics%cdf(t,j,k,i-1)+&
                                 metrics%cdf_bincounts(t,j,k,i)/&
                                 float(sum(metrics%cdf_bincounts(t,j,k,:)))
                         enddo
                      endif
                   enddo
                enddo
             enddo
          elseif(LDT_rc%group_cdfs.eq.1 .and. LDT_rc%strat_cdfs.eq.0) then 

             allocate(strat_bincounts(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels,&
                  LDT_rc%cdf_nbins))
             allocate(strat_cdf(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels,&
                  LDT_rc%cdf_nbins))

             strat_bincounts = 0
             strat_cdf = 0 

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      do i=1,LDT_rc%cdf_nbins
                         sindex = LDT_rc%cdf_strat_data(t)
                         strat_bincounts(sindex,j,k,i) = &
                              strat_bincounts(sindex,j,k,i) + & 
                              metrics%cdf_bincounts(t,j,k,i)

                      enddo
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%group_cdfs_nbins
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(sum(strat_bincounts(t,j,k,:)).ne.0) then
                         strat_cdf(t,j,k,1) = &
                              strat_cdf(t,j,k,1)/&
                              float(sum(strat_bincounts(t,j,k,:)))
                         do i=2,LDT_rc%cdf_nbins
                            strat_cdf(t,j,k,i) = strat_cdf(t,j,k,i-1) + & 
                                 strat_bincounts(t,j,k,i)/&
                                 float(sum(strat_bincounts(t,j,k,:)))
                         enddo
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%ngrid(n)
!                if(datamask(t).gt.0) then 
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      do i=1,LDT_rc%cdf_nbins
                         sindex = LDT_rc%cdf_strat_data(t)
                         metrics%cdf(t,j,k,i) = strat_cdf(sindex,j,k,i)
                      enddo
                   enddo
                enddo
!                endif
             enddo
             ! MN: save stratified CDF 
             metrics%strat_cdf = strat_cdf
             deallocate(strat_bincounts)
             deallocate(strat_cdf)

!MN: Startification based on monthly precipitation climatology 
!    monthly total precipitation climatology for each pixel are stored in LDT_rc%stratification_data

          elseif(LDT_rc%group_cdfs.eq.0 .and. LDT_rc%strat_cdfs.eq.1 ) then

             allocate(strat_bincounts(LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels,&
                  LDT_rc%cdf_nbins))
             allocate(strat_cdf(LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels,&
                  LDT_rc%cdf_nbins))

             strat_bincounts = 0
             strat_cdf = 0

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      do i=1,LDT_rc%cdf_nbins
                         sindex = LDT_rc%stratification_data(t,j)
                         strat_bincounts(sindex,j,k,i) = &
                              strat_bincounts(sindex,j,k,i) + &
                              metrics%cdf_bincounts(t,j,k,i)

                      enddo
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%strat_cdfs_nbins
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(sum(strat_bincounts(t,j,k,:)).ne.0) then
                         strat_cdf(t,j,k,1) = &
                              strat_bincounts(t,j,k,1)/&  ! strat_cdf(t,j,k,1)/&
                              float(sum(strat_bincounts(t,j,k,:)))
                         do i=2,LDT_rc%cdf_nbins
                            strat_cdf(t,j,k,i) = strat_cdf(t,j,k,i-1) + &
                                 strat_bincounts(t,j,k,i)/&
                                 float(sum(strat_bincounts(t,j,k,:)))
                         enddo
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      do i=1,LDT_rc%cdf_nbins
                         sindex = LDT_rc%stratification_data(t,j)
                         metrics%cdf(t,j,k,i) = strat_cdf(sindex,j,k,i)
                      enddo
                   enddo
                enddo
             enddo
             ! MN: save stratified CDF 
             metrics%strat_cdf = strat_cdf
             deallocate(strat_bincounts)
             deallocate(strat_cdf)

          elseif(LDT_rc%group_cdfs.eq.1 .and. LDT_rc%strat_cdfs.eq.1 ) then

             allocate(strat_bincounts(LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels,&
                  LDT_rc%cdf_nbins))
             allocate(strat_cdf(LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels,&
                  LDT_rc%cdf_nbins))

             strat_bincounts = 0
             strat_cdf = 0

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      do i=1,LDT_rc%cdf_nbins
                         sindex0  = LDT_rc%stratification_data(t,j)
                         sindex1 = LDT_rc%cdf_strat_data(t)
                         sindex = sindex0 + (sindex1 - 1)*LDT_rc%strat_cdfs_nbins
                         strat_bincounts(sindex,j,k,i) = &
                              strat_bincounts(sindex,j,k,i) + &
                              metrics%cdf_bincounts(t,j,k,i)

                      enddo
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%group_cdfs_nbins*LDT_rc%strat_cdfs_nbins
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(sum(strat_bincounts(t,j,k,:)).ne.0) then
                         strat_cdf(t,j,k,1) = &
                              strat_bincounts(t,j,k,1)/&  ! strat_cdf(t,j,k,1)/&
                              float(sum(strat_bincounts(t,j,k,:)))
                         do i=2,LDT_rc%cdf_nbins
                            strat_cdf(t,j,k,i) = strat_cdf(t,j,k,i-1) + &
                                 strat_bincounts(t,j,k,i)/&
                                 float(sum(strat_bincounts(t,j,k,:)))
                         enddo
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      do i=1,LDT_rc%cdf_nbins
                         sindex0  = LDT_rc%stratification_data(t,j)
                         sindex1 = LDT_rc%cdf_strat_data(t)
                         sindex = sindex0 + (sindex1 - 1)*LDT_rc%strat_cdfs_nbins
                         metrics%cdf(t,j,k,i) = strat_cdf(sindex,j,k,i)
                      enddo
                   enddo
                enddo
             enddo

             ! MN: save stratified CDF
             metrics%strat_cdf = strat_cdf
             deallocate(strat_bincounts)
             deallocate(strat_cdf)

          endif
       endif
    endif

  end subroutine computeSingleCDF


end module LDT_CDFMod
