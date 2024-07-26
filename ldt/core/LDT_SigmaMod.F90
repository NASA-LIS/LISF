!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_SigmaMod
!
!BOP
! !MODULE: LDT_SigmaMod
! 
!  !DESCRIPTION: 
!   This module handles the standard deviation computations 
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
!
!EOP
  use LDT_DAobsDataMod
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_diagnoseSigma
  public :: LDT_computeSigma

  private
  
contains

!BOP
! !ROUTINE: LDT_diagnoseSigma
! \label{LDT_diagnoseSigma}
!
! !INTERFACE: 
  subroutine LDT_diagnoseSigma(n)
! !USES:     
    use LDT_coreMod,             only : LDT_rc
    use LDT_DAmetricsDataMod,        only : LDT_DAmetricsPtr
! 
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Sigma calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleSigma](\ref{diagnoseSingleSigma})
!     updates the Sigma computation for a single variable 
!   \end{description}
! 
!EOP
    implicit none
    integer, intent(in) :: n
    integer :: i, index

    do index=1,LDT_DA_MOC_COUNT
       call diagnoseSingleSigma(n,LDT_DAobsDataPtr(n,index)%dataEntryPtr,&
            LDT_DAmetricsPtr(index)%dataEntryPtr)
    enddo

  end subroutine LDT_diagnoseSigma

!BOP
! !ROUTINE: diagnoseSingleSigma
! \label{diagnoseSingleSigma}
! 
! !INTERFACE: 
  subroutine diagnoseSingleSigma(n,obs, metrics)
! !USES: 
    use LDT_coreMod
    use LDT_DAmetricsDataMod
    use LDT_timeMgrMod
! 
! !DESCRIPTION: 
!  This routine updates the Sigma computation (updates the running 
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

    integer    :: t,k,j,c,r,c1,r1,t1
    integer    :: r_min, r_max
    integer    :: c_min, c_max
    real       :: lon                 !Y.Kwon
    real       :: gmt                 !Y.Kwon
    real       :: lhour               !Y.Kwon
    integer    :: zone                !Y.Kwon 
 
    if(LDT_rc%cdf_ntimes.eq.12) then 
       j = LDT_rc%mo
    elseif(LDT_rc%cdf_ntimes.eq.1) then 
       j = 1
    elseif(LDT_rc%cdf_ntimes.eq.24.or.&
           LDT_rc%cdf_ntimes.eq.365) then   !half-monthly (for 9-km operational SMAP; Y.Kwon) 
       if(LDT_rc%mo.eq.2) then
          if(LDT_rc%da.le.14) then
             j = LDT_rc%mo * 2 - 1
          else
             j = LDT_rc%mo * 2
          endif
       else
          if(LDT_rc%da.le.15) then
             j = LDT_rc%mo * 2 - 1
          else
             j = LDT_rc%mo * 2
          endif
       endif
    endif
    
    if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 

       if(LDT_rc%daily_interp_switch.eq.0) then     !Y.Kwon
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
                               if(obs%count(t,k).ne.0) then 
                                  metrics%sxx_sigma(t,j,k) = metrics%sxx_sigma(t,j,k) + &
                                       obs%value(t1,k)*obs%value(t1,k)
                                  metrics%sx_sigma(t,j,k) = metrics%sx_sigma(t,j,k) + &
                                       obs%value(t1,k)
                                  metrics%count_sigma(t,j,k) = &
                                       metrics%count_sigma(t,j,k) + 1
                               endif
                            enddo
                         endif
                      enddo
                   enddo
                endif
             enddo
          enddo

       !Y.Kwon
       elseif(LDT_rc%daily_interp_switch.eq.1) then
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
                                  lon = LDT_domain(n)%lon(c1+(r1-1)*LDT_rc%lnc(n))
                                  gmt = LDT_rc%hr
                                  call LDT_gmt2localtime(gmt, lon, lhour, zone)

                                  !if(lhour.eq.6) then ! Orig
                                  if(lhour > 4 .and. lhour < 8) then ! EMK for 3-hrly data
                                     metrics%sxx_sigma_6am(t,j,k) = metrics%sxx_sigma_6am(t,j,k) + &
                                          obs%value(t1,k)*obs%value(t1,k)
                                     metrics%sx_sigma_6am(t,j,k) = metrics%sx_sigma_6am(t,j,k) + &
                                          obs%value(t1,k)
                                     metrics%count_sigma_6am(t,j,k) = &
                                          metrics%count_sigma_6am(t,j,k) + 1
                                   !elseif(lhour.eq.18) then ! Orig
                                  elseif (lhour > 16 .and. lhour < 20) then ! EMK for 3-hrly data
                                     metrics%sxx_sigma_6pm(t,j,k) = metrics%sxx_sigma_6pm(t,j,k) + &
                                          obs%value(t1,k)*obs%value(t1,k)
                                     metrics%sx_sigma_6pm(t,j,k) = metrics%sx_sigma_6pm(t,j,k) + &
                                          obs%value(t1,k)
                                     metrics%count_sigma_6pm(t,j,k) = &
                                          metrics%count_sigma_6pm(t,j,k) + 1
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
    endif
    
  end subroutine diagnoseSingleSigma


!BOP
! 
! !ROUTINE: LDT_computeSigma
! \label{LDT_computeSigma}
! 
! !INTERFACE: 
  subroutine LDT_computeSigma(n)
! !USES: 
    use LDT_coreMod,  only : LDT_rc
    use LDT_DAmetricsDataMod,        only : LDT_DAmetricsPtr

! 
! !DESCRIPTION: 
!   This subroutine issues the calls to compute Sigma values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleSigma](\ref{computeSingleSigma})
!     updates the Sigma computation for a single variable 
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
       call computeSingleSigma(n,LDT_DAobsDataPtr(n,index)%dataEntryPtr,&
            LDT_DAmetricsPtr(index)%dataEntryPtr)
    enddo

  end subroutine LDT_ComputeSigma

!BOP
! 
! !ROUTINE: computeSingleSigma
! \label{computeSingleSigma}
! 
! !INTERFACE: 
  subroutine computeSingleSigma(n,obs, metrics)
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_DAmetricsDataMod
    use LDT_logMod
! 
! !DESCRIPTION: 
!  This routine computes the Sigma values for a single variable
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

    integer           :: sindex
    real, allocatable :: strat_sxx_sigma(:,:,:)
    real, allocatable :: strat_sx_sigma(:,:,:)
    real, allocatable :: strat_sigma(:,:,:)
    integer, allocatable :: strat_count(:,:,:)
    integer  :: t,i,k,j
    real     :: term1, term2
    real     :: sigma_temp_6am(LDT_rc%ngrid(n),24,obs%vlevels)   !Y.Kwon 
    real     :: sigma_temp_6pm(LDT_rc%ngrid(n),24,obs%vlevels)   !Y.Kwon
    real     :: d_sigma_temp_6am                                 !Y.Kwon
    real     :: d_sigma_temp_6pm                                 !Y.kwon
    integer  :: d_days                                           !Y.Kwon

    if(LDT_rc%endtime.eq.1) then 
       if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
          if(LDT_rc%group_cdfs.eq.0) then 
             do t=1,LDT_rc%ngrid(n)

                !for 9-km operational SMAP (Y.Kwon)
                if(LDT_rc%daily_interp_switch.eq.1) then
                   do j=1,24
                      do k=1,obs%vlevels
                         if(metrics%count_sigma_6am(t,j,k).le.LDT_rc%obsCountThreshold) then
                            sigma_temp_6am(t,j,k) = LDT_rc%udef
                            metrics%mask_6am(t,j,k) = LDT_rc%udef
                         else
                            term1 = metrics%sxx_sigma_6am(t,j,k)/&
                                 metrics%count_sigma_6am(t,j,k)
                            term2 = (metrics%sx_sigma_6am(t,j,k)/metrics%count_sigma_6am(t,j,k))**2
                            if(term1.gt.term2) then
                               sigma_temp_6am(t,j,k) = sqrt(term1-term2)
                               metrics%mask_6am(t,j,k) = metrics%count_sigma_6am(t,j,k)
                            else
                               metrics%mask_6am(t,j,k) = LDT_rc%udef
                               sigma_temp_6am(t,j,k) = 0
                               write(LDT_logunit,*) '[WARN] setting sigma_6am to zero ',t,j,k,&
                                    metrics%count_sigma_6am(t,j,k)
                            endif
                         endif

                         if(metrics%count_sigma_6pm(t,j,k).le.LDT_rc%obsCountThreshold) then
                            sigma_temp_6pm(t,j,k) = LDT_rc%udef
                            metrics%mask_6pm(t,j,k) = LDT_rc%udef
                         else
                            term1 = metrics%sxx_sigma_6pm(t,j,k)/&
                                 metrics%count_sigma_6pm(t,j,k)
                            term2 = (metrics%sx_sigma_6pm(t,j,k)/metrics%count_sigma_6pm(t,j,k))**2
                            if(term1.gt.term2) then
                               sigma_temp_6pm(t,j,k) = sqrt(term1-term2)
                               metrics%mask_6pm(t,j,k) = metrics%count_sigma_6pm(t,j,k)
                            else
                               metrics%mask_6pm(t,j,k) = LDT_rc%udef
                               sigma_temp_6pm(t,j,k) = 0
                               write(LDT_logunit,*) '[WARN] setting sigma_6pm to zero ',t,j,k,&
                                    metrics%count_sigma_6pm(t,j,k)
                            endif
                         endif
                      enddo
                   enddo

                   !daily interpolation using half-monthly values
                   do j=1,LDT_rc%cdf_ntimes    !365
                      do k=1,obs%vlevels
                         if(j<8) then
                            if(sigma_temp_6am(t,1,k).gt.0.and.sigma_temp_6am(t,24,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,1,k) - sigma_temp_6am(t,24,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,24,k) + (d_sigma_temp_6am/d_days) * (j+7)
                            else
                               if(sigma_temp_6am(t,1,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,1,k)
                               elseif(sigma_temp_6am(t,24,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,24,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,1,k).gt.0.and.sigma_temp_6pm(t,24,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,1,k) - sigma_temp_6pm(t,24,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,24,k) + (d_sigma_temp_6pm/d_days) * (j+7)
                            else
                               if(sigma_temp_6pm(t,1,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,1,k)
                               elseif(sigma_temp_6pm(t,24,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,24,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==8) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,1,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,1,k)
                         elseif(j<24) then
                            if(sigma_temp_6am(t,1,k).gt.0.and.sigma_temp_6am(t,2,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,2,k) - sigma_temp_6am(t,1,k)
                               d_days = 16
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,1,k) + (d_sigma_temp_6am/d_days) * (j-8)
                            else
                               if(sigma_temp_6am(t,1,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,1,k)
                               elseif(sigma_temp_6am(t,2,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,2,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,1,k).gt.0.and.sigma_temp_6pm(t,2,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,2,k) - sigma_temp_6pm(t,1,k)
                               d_days = 16
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,1,k) + (d_sigma_temp_6pm/d_days) * (j-8)
                            else
                               if(sigma_temp_6pm(t,1,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,1,k)
                               elseif(sigma_temp_6pm(t,2,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,2,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==24) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,2,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,2,k)
                         elseif(j<39) then
                            if(sigma_temp_6am(t,2,k).gt.0.and.sigma_temp_6am(t,3,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,3,k) - sigma_temp_6am(t,2,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,2,k) + (d_sigma_temp_6am/d_days) * (j-24)
                            else
                               if(sigma_temp_6am(t,2,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,2,k)
                               elseif(sigma_temp_6am(t,3,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,3,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,2,k).gt.0.and.sigma_temp_6pm(t,3,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,3,k) - sigma_temp_6pm(t,2,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,2,k) + (d_sigma_temp_6pm/d_days) * (j-24)
                            else
                               if(sigma_temp_6pm(t,2,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,2,k)
                               elseif(sigma_temp_6pm(t,3,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,3,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==39) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,3,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,3,k)
                         elseif(j<53) then
                            if(sigma_temp_6am(t,3,k).gt.0.and.sigma_temp_6am(t,4,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,4,k) - sigma_temp_6am(t,3,k)
                               d_days = 14
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,3,k) + (d_sigma_temp_6am/d_days) * (j-39)
                            else
                               if(sigma_temp_6am(t,3,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,3,k)
                               elseif(sigma_temp_6am(t,4,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,4,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,3,k).gt.0.and.sigma_temp_6pm(t,4,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,4,k) - sigma_temp_6pm(t,3,k)
                               d_days = 14
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,3,k) + (d_sigma_temp_6pm/d_days) * (j-39)
                            else
                               if(sigma_temp_6pm(t,3,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,3,k)
                               elseif(sigma_temp_6pm(t,4,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,4,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==53) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,4,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,4,k)
                         elseif(j<67) then
                            if(sigma_temp_6am(t,4,k).gt.0.and.sigma_temp_6am(t,5,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,5,k) - sigma_temp_6am(t,4,k)
                               d_days = 14
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,4,k) + (d_sigma_temp_6am/d_days) * (j-53)
                            else
                               if(sigma_temp_6am(t,4,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,4,k)
                               elseif(sigma_temp_6am(t,5,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,5,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,4,k).gt.0.and.sigma_temp_6pm(t,5,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,5,k) - sigma_temp_6pm(t,4,k)
                               d_days = 14
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,4,k) + (d_sigma_temp_6pm/d_days) * (j-53)
                            else
                               if(sigma_temp_6pm(t,4,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,4,k)
                               elseif(sigma_temp_6pm(t,5,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,5,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==67) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,5,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,5,k)
                         elseif(j<83) then
                            if(sigma_temp_6am(t,5,k).gt.0.and.sigma_temp_6am(t,6,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,6,k) - sigma_temp_6am(t,5,k)
                               d_days = 16
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,5,k) + (d_sigma_temp_6am/d_days) * (j-67)
                            else
                               if(sigma_temp_6am(t,5,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,5,k)
                               elseif(sigma_temp_6am(t,6,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,6,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,5,k).gt.0.and.sigma_temp_6pm(t,6,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,6,k) - sigma_temp_6pm(t,5,k)
                               d_days = 16
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,5,k) + (d_sigma_temp_6pm/d_days) * (j-67)
                            else
                               if(sigma_temp_6pm(t,5,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,5,k)
                               elseif(sigma_temp_6pm(t,6,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,6,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==83) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,6,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,6,k)
                         elseif(j<98) then
                            if(sigma_temp_6am(t,6,k).gt.0.and.sigma_temp_6am(t,7,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,7,k) - sigma_temp_6am(t,6,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,6,k) + (d_sigma_temp_6am/d_days) * (j-83)
                            else
                               if(sigma_temp_6am(t,6,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,6,k)
                               elseif(sigma_temp_6am(t,7,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,7,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,6,k).gt.0.and.sigma_temp_6pm(t,7,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,7,k) - sigma_temp_6pm(t,6,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,6,k) + (d_sigma_temp_6pm/d_days) * (j-83)
                            else
                               if(sigma_temp_6pm(t,6,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,6,k)
                               elseif(sigma_temp_6pm(t,7,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,7,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==98) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,7,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,7,k)
                         elseif(j<113) then
                            if(sigma_temp_6am(t,7,k).gt.0.and.sigma_temp_6am(t,8,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,8,k) - sigma_temp_6am(t,7,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,7,k) + (d_sigma_temp_6am/d_days) * (j-98)
                            else
                               if(sigma_temp_6am(t,7,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,7,k)
                               elseif(sigma_temp_6am(t,8,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,8,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,7,k).gt.0.and.sigma_temp_6pm(t,8,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,8,k) - sigma_temp_6pm(t,7,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,7,k) + (d_sigma_temp_6pm/d_days) * (j-98)
                            else
                               if(sigma_temp_6pm(t,7,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,7,k)
                               elseif(sigma_temp_6pm(t,8,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,8,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==113) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,8,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,8,k)
                         elseif(j<128) then
                            if(sigma_temp_6am(t,8,k).gt.0.and.sigma_temp_6am(t,9,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,9,k) - sigma_temp_6am(t,8,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,8,k) + (d_sigma_temp_6am/d_days) * (j-113)
                            else
                               if(sigma_temp_6am(t,8,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,8,k)
                               elseif(sigma_temp_6am(t,9,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,9,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,8,k).gt.0.and.sigma_temp_6pm(t,9,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,9,k) - sigma_temp_6pm(t,8,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,8,k) + (d_sigma_temp_6pm/d_days) * (j-113)
                            else
                               if(sigma_temp_6pm(t,8,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,8,k)
                               elseif(sigma_temp_6pm(t,9,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,9,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==128) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,9,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,9,k)
                         elseif(j<144) then
                            if(sigma_temp_6am(t,9,k).gt.0.and.sigma_temp_6am(t,10,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,10,k) - sigma_temp_6am(t,9,k)
                               d_days = 16
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,9,k) + (d_sigma_temp_6am/d_days) * (j-128)
                            else
                               if(sigma_temp_6am(t,9,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,9,k)
                               elseif(sigma_temp_6am(t,10,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,10,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,9,k).gt.0.and.sigma_temp_6pm(t,10,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,10,k) - sigma_temp_6pm(t,9,k)
                               d_days = 16
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,9,k) + (d_sigma_temp_6pm/d_days) * (j-128)
                            else
                               if(sigma_temp_6pm(t,9,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,9,k)
                               elseif(sigma_temp_6pm(t,10,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,10,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==144) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,10,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,10,k)
                         elseif(j<159) then
                            if(sigma_temp_6am(t,10,k).gt.0.and.sigma_temp_6am(t,11,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,11,k) - sigma_temp_6am(t,10,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,10,k) + (d_sigma_temp_6am/d_days) * (j-144)
                            else
                               if(sigma_temp_6am(t,10,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,10,k)
                               elseif(sigma_temp_6am(t,11,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,11,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,10,k).gt.0.and.sigma_temp_6pm(t,11,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,11,k) - sigma_temp_6pm(t,10,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,10,k) + (d_sigma_temp_6pm/d_days) * (j-144)
                            else
                               if(sigma_temp_6pm(t,10,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,10,k)
                               elseif(sigma_temp_6pm(t,11,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,11,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==159) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,11,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,11,k)
                         elseif(j<174) then
                            if(sigma_temp_6am(t,11,k).gt.0.and.sigma_temp_6am(t,12,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,12,k) - sigma_temp_6am(t,11,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,11,k) + (d_sigma_temp_6am/d_days) * (j-159)
                            else
                               if(sigma_temp_6am(t,11,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,11,k)
                               elseif(sigma_temp_6am(t,12,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,12,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,11,k).gt.0.and.sigma_temp_6pm(t,12,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,12,k) - sigma_temp_6pm(t,11,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,11,k) + (d_sigma_temp_6pm/d_days) * (j-159)
                            else
                               if(sigma_temp_6pm(t,11,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,11,k)
                               elseif(sigma_temp_6pm(t,12,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,12,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==174) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,12,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,12,k)
                         elseif(j<189) then
                            if(sigma_temp_6am(t,12,k).gt.0.and.sigma_temp_6am(t,13,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,13,k) - sigma_temp_6am(t,12,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,12,k) + (d_sigma_temp_6am/d_days) * (j-174)
                            else
                               if(sigma_temp_6am(t,12,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,12,k)
                               elseif(sigma_temp_6am(t,13,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,13,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,12,k).gt.0.and.sigma_temp_6pm(t,13,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,13,k) - sigma_temp_6pm(t,12,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,12,k) + (d_sigma_temp_6pm/d_days) * (j-174)
                            else
                               if(sigma_temp_6pm(t,12,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,12,k)
                               elseif(sigma_temp_6pm(t,13,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,13,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==189) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,13,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,13,k)
                         elseif(j<205) then
                            if(sigma_temp_6am(t,13,k).gt.0.and.sigma_temp_6am(t,14,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,14,k) - sigma_temp_6am(t,13,k)
                               d_days = 16
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,13,k) + (d_sigma_temp_6am/d_days) * (j-189)
                            else
                               if(sigma_temp_6am(t,13,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,13,k)
                               elseif(sigma_temp_6am(t,14,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,14,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,13,k).gt.0.and.sigma_temp_6pm(t,14,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,14,k) - sigma_temp_6pm(t,13,k)
                               d_days = 16
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,13,k) + (d_sigma_temp_6pm/d_days) * (j-189)
                            else
                               if(sigma_temp_6pm(t,13,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,13,k)
                               elseif(sigma_temp_6pm(t,14,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,14,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==205) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,14,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,14,k)
                         elseif(j<220) then
                            if(sigma_temp_6am(t,14,k).gt.0.and.sigma_temp_6am(t,15,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,15,k) - sigma_temp_6am(t,14,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,14,k) + (d_sigma_temp_6am/d_days) * (j-205)
                            else
                               if(sigma_temp_6am(t,14,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,14,k)
                               elseif(sigma_temp_6am(t,15,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,15,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,14,k).gt.0.and.sigma_temp_6pm(t,15,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,15,k) - sigma_temp_6pm(t,14,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,14,k) + (d_sigma_temp_6pm/d_days) * (j-205)
                            else
                               if(sigma_temp_6pm(t,14,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,14,k)
                               elseif(sigma_temp_6pm(t,15,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,15,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==220) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,15,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,15,k)
                         elseif(j<236) then
                            if(sigma_temp_6am(t,15,k).gt.0.and.sigma_temp_6am(t,16,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,16,k) - sigma_temp_6am(t,15,k)
                               d_days = 16
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,15,k) + (d_sigma_temp_6am/d_days) * (j-220)
                            else
                               if(sigma_temp_6am(t,15,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,15,k)
                               elseif(sigma_temp_6am(t,16,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,16,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,15,k).gt.0.and.sigma_temp_6pm(t,16,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,16,k) - sigma_temp_6pm(t,15,k)
                               d_days = 16
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,15,k) + (d_sigma_temp_6pm/d_days) * (j-220)
                            else
                               if(sigma_temp_6pm(t,15,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,15,k)
                               elseif(sigma_temp_6pm(t,16,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,16,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==236) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,16,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,16,k)
                         elseif(j<251) then
                            if(sigma_temp_6am(t,16,k).gt.0.and.sigma_temp_6am(t,17,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,17,k) - sigma_temp_6am(t,16,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,16,k) + (d_sigma_temp_6am/d_days) * (j-236)
                            else
                               if(sigma_temp_6am(t,16,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,16,k)
                               elseif(sigma_temp_6am(t,17,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,17,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,16,k).gt.0.and.sigma_temp_6pm(t,17,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,17,k) - sigma_temp_6pm(t,16,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,16,k) + (d_sigma_temp_6pm/d_days) * (j-236)
                            else
                               if(sigma_temp_6pm(t,16,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,16,k)
                               elseif(sigma_temp_6pm(t,17,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,17,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==251) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,17,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,17,k)
                         elseif(j<266) then
                            if(sigma_temp_6am(t,17,k).gt.0.and.sigma_temp_6am(t,18,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,18,k) - sigma_temp_6am(t,17,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,17,k) + (d_sigma_temp_6am/d_days) * (j-251)
                            else
                               if(sigma_temp_6am(t,17,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,17,k)
                               elseif(sigma_temp_6am(t,18,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,18,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,17,k).gt.0.and.sigma_temp_6pm(t,18,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,18,k) - sigma_temp_6pm(t,17,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,17,k) + (d_sigma_temp_6pm/d_days) * (j-251)
                            else
                               if(sigma_temp_6pm(t,17,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,17,k)
                               elseif(sigma_temp_6pm(t,18,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,18,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==266) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,18,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,18,k)
                         elseif(j<281) then
                            if(sigma_temp_6am(t,18,k).gt.0.and.sigma_temp_6am(t,19,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,19,k) - sigma_temp_6am(t,18,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,18,k) + (d_sigma_temp_6am/d_days) * (j-266)
                            else
                               if(sigma_temp_6am(t,18,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,18,k)
                               elseif(sigma_temp_6am(t,19,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,19,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,18,k).gt.0.and.sigma_temp_6pm(t,19,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,19,k) - sigma_temp_6pm(t,18,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,18,k) + (d_sigma_temp_6pm/d_days) * (j-266)
                            else
                               if(sigma_temp_6pm(t,18,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,18,k)
                               elseif(sigma_temp_6pm(t,19,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,19,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==281) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,19,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,19,k)
                         elseif(j<297) then
                            if(sigma_temp_6am(t,19,k).gt.0.and.sigma_temp_6am(t,20,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,20,k) - sigma_temp_6am(t,19,k)
                               d_days = 16
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,19,k) + (d_sigma_temp_6am/d_days) * (j-281)
                            else
                               if(sigma_temp_6am(t,19,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,19,k)
                               elseif(sigma_temp_6am(t,20,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,20,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,19,k).gt.0.and.sigma_temp_6pm(t,20,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,20,k) - sigma_temp_6pm(t,19,k)
                               d_days = 16
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,19,k) + (d_sigma_temp_6pm/d_days) * (j-281)
                            else
                               if(sigma_temp_6pm(t,19,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,19,k)
                               elseif(sigma_temp_6pm(t,20,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,20,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==297) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,20,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,20,k)
                         elseif(j<312) then
                            if(sigma_temp_6am(t,20,k).gt.0.and.sigma_temp_6am(t,21,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,21,k) - sigma_temp_6am(t,20,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,20,k) + (d_sigma_temp_6am/d_days) * (j-297)
                            else
                               if(sigma_temp_6am(t,20,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,20,k)
                               elseif(sigma_temp_6am(t,21,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,21,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,20,k).gt.0.and.sigma_temp_6pm(t,21,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,21,k) - sigma_temp_6pm(t,20,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,20,k) + (d_sigma_temp_6pm/d_days) * (j-297)
                            else
                               if(sigma_temp_6pm(t,20,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,20,k)
                               elseif(sigma_temp_6pm(t,21,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,21,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==312) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,21,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,21,k)
                         elseif(j<327) then
                            if(sigma_temp_6am(t,21,k).gt.0.and.sigma_temp_6am(t,22,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,22,k) - sigma_temp_6am(t,21,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,21,k) + (d_sigma_temp_6am/d_days) * (j-312)
                            else
                               if(sigma_temp_6am(t,21,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,21,k)
                               elseif(sigma_temp_6am(t,22,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,22,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,21,k).gt.0.and.sigma_temp_6pm(t,22,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,22,k) - sigma_temp_6pm(t,21,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,21,k) + (d_sigma_temp_6pm/d_days) * (j-312)
                            else
                               if(sigma_temp_6pm(t,21,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,21,k)
                               elseif(sigma_temp_6pm(t,22,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,22,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==327) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,22,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,22,k)
                         elseif(j<342) then
                            if(sigma_temp_6am(t,22,k).gt.0.and.sigma_temp_6am(t,23,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,23,k) - sigma_temp_6am(t,22,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,22,k) + (d_sigma_temp_6am/d_days) * (j-327)
                            else
                               if(sigma_temp_6am(t,22,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,22,k)
                               elseif(sigma_temp_6am(t,23,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,23,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,22,k).gt.0.and.sigma_temp_6pm(t,23,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,23,k) - sigma_temp_6pm(t,22,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,22,k) + (d_sigma_temp_6pm/d_days) * (j-327)
                            else
                               if(sigma_temp_6pm(t,22,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,22,k)
                               elseif(sigma_temp_6pm(t,23,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,23,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==342) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,23,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,23,k)
                         elseif(j<358) then
                            if(sigma_temp_6am(t,23,k).gt.0.and.sigma_temp_6am(t,24,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,24,k) - sigma_temp_6am(t,23,k)
                               d_days = 16
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,23,k) + (d_sigma_temp_6am/d_days) * (j-342)
                            else
                               if(sigma_temp_6am(t,23,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,23,k)
                               elseif(sigma_temp_6am(t,24,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,24,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,23,k).gt.0.and.sigma_temp_6pm(t,24,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,24,k) - sigma_temp_6pm(t,23,k)
                               d_days = 16
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,23,k) + (d_sigma_temp_6pm/d_days) * (j-342)
                            else
                               if(sigma_temp_6pm(t,23,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,23,k)
                               elseif(sigma_temp_6pm(t,24,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,24,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==358) then
                            metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,24,k)
                            metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,24,k)
                         elseif(j<=365) then
                            if(sigma_temp_6am(t,24,k).gt.0.and.sigma_temp_6am(t,1,k).gt.0) then
                               d_sigma_temp_6am = sigma_temp_6am(t,1,k) - sigma_temp_6am(t,24,k)
                               d_days = 15
                               metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,24,k) + (d_sigma_temp_6am/d_days) * (j-358)
                            else
                               if(sigma_temp_6am(t,24,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,24,k)
                               elseif(sigma_temp_6am(t,1,k).gt.0) then
                                  metrics%sigma_6am(t,j,k) = sigma_temp_6am(t,1,k)
                               else
                                  metrics%sigma_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(sigma_temp_6pm(t,24,k).gt.0.and.sigma_temp_6pm(t,1,k).gt.0) then
                               d_sigma_temp_6pm = sigma_temp_6pm(t,1,k) - sigma_temp_6pm(t,24,k)
                               d_days = 15
                               metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,24,k) + (d_sigma_temp_6pm/d_days) * (j-358)
                            else
                               if(sigma_temp_6pm(t,24,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,24,k)
                               elseif(sigma_temp_6pm(t,1,k).gt.0) then
                                  metrics%sigma_6pm(t,j,k) = sigma_temp_6pm(t,1,k)
                               else
                                  metrics%sigma_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         endif
                      enddo
                   enddo

                else
                   do j=1,LDT_rc%cdf_ntimes
                      do k=1,obs%vlevels
                
                         if(metrics%count_sigma(t,j,k).le.LDT_rc%obsCountThreshold) then 
                            metrics%sigma(t,j,k) = LDT_rc%udef
                         else
                            term1 = metrics%sxx_sigma(t,j,k)/&
                                 metrics%count_sigma(t,j,k)
                            term2 = (metrics%sx_sigma(t,j,k)/metrics%count_sigma(t,j,k))**2 
                            if(term1.gt.term2) then 
                               metrics%sigma(t,j,k) = sqrt(term1-term2)
                            else
                               metrics%mask(t,j,k) = LDT_rc%udef
                               metrics%sigma(t,j,k) = 0
                               write(LDT_logunit,*) '[WARN] setting sigma to zero ',t,j,k,&
                                    metrics%count_sigma(t,j,k)
                            endif
                         endif
                      enddo
                   enddo
                endif
             enddo
          else
             allocate(strat_sxx_sigma(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_sx_sigma(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_sigma(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_count(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))

             strat_sxx_sigma = 0 
             strat_sx_sigma = 0 
             strat_sigma = 0 
             strat_count = 0 

              do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex = LDT_rc%cdf_strat_data(t)
                      strat_sxx_sigma(sindex,j,k) = &
                           strat_sxx_sigma(sindex,j,k) + &
                           metrics%sxx_sigma(t,j,k)
                      strat_sx_sigma(sindex,j,k) = &
                           strat_sx_sigma(sindex,j,k) + &
                           metrics%sx_sigma(t,j,k)
                      strat_count(sindex,j,k) = & 
                           strat_count(sindex,j,k) + &
                           metrics%count_sigma(t,j,k)
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%group_cdfs_nbins
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(strat_count(t,j,k).le.&
                           LDT_rc%obsCountThreshold) then 
                         strat_sigma(t,j,k) = LDT_rc%udef
                      else
                         term1 = strat_sxx_sigma(t,j,k)/&
                              strat_count(t,j,k)
                         term2 = (strat_sx_sigma(t,j,k)/strat_count(t,j,k))**2 
                         if(term1.gt.term2) then 
                            strat_sigma(t,j,k) = sqrt(term1-term2)
                         else
                            strat_sigma(t,j,k) = 0
                            print*, 'warning: setting sigma to zero ',t,j,k,&
                                 strat_sigma(t,j,k)
                         endif
                         
                      endif

                   enddo
                enddo
             enddo

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex = LDT_rc%cdf_strat_data(t)
                      if(strat_sigma(sindex,j,k).eq.LDT_rc%udef) then 
                         metrics%sigma(t,j,k) = LDT_rc%udef
                      else
                         metrics%sigma(t,j,k) = strat_sigma(sindex,j,k)
                      endif
                   enddo
                enddo
             enddo

             deallocate(strat_sigma)
             deallocate(strat_sxx_sigma)
             deallocate(strat_sx_sigma)
             
          endif
       endif
    endif

  end subroutine computeSingleSigma
end module LDT_SigmaMod
