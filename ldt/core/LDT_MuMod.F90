!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_MuMod
!
!BOP
! !MODULE: LDT_MuMod
! 
!  !DESCRIPTION: 
!   This module handles the mean computations 
!
!  !REVISION HISTORY: 
!  2 Oct 2008    Sujay Kumar  Initial Specification
!
!EOP
  use LDT_DAobsDataMod
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_diagnoseMu
  public :: LDT_computeMu

  private
  
contains

!BOP
! !ROUTINE: LDT_diagnoseMu
! \label{LDT_diagnoseMu}
!
! !INTERFACE: 
  subroutine LDT_diagnoseMu(n)
! !USES:     
    use LDT_coreMod,             only : LDT_rc
    use LDT_DAmetricsDataMod,        only : LDT_DAmetricsPtr
! 
! !DESCRIPTION: 
!   This subroutine issues the calls to update the Mu calculation for 
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleMu](\ref{diagnoseSingleMu})
!     updates the Mu computation for a single variable 
!   \end{description}
! 
!EOP
    implicit none
    integer, intent(in) :: n
    integer :: i, index

    do index=1,LDT_DA_MOC_COUNT
       call diagnoseSingleMu(n,LDT_DAobsDataPtr(n,index)%dataEntryPtr,&
            LDT_DAmetricsPtr(index)%dataEntryPtr)
    enddo

  end subroutine LDT_diagnoseMu

!BOP
! !ROUTINE: diagnoseSingleMu
! \label{diagnoseSingleMu}
! 
! !INTERFACE: 
  subroutine diagnoseSingleMu(n,obs, metrics)
! !USES: 
    use LDT_coreMod
    use LDT_DAmetricsDataMod
    use LDT_timeMgrMod          !Y.Kwon
! 
! !DESCRIPTION: 
!  This routine updates the Mu computation (updates the running 
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
                            
                               if(obs%count(t1,k).ne.0) then 
                                  metrics%sx_mu(t,j,k) = metrics%sx_mu(t,j,k) + &
                                       obs%value(t1,k)
                                  metrics%count_mu(t,j,k) = &
                                       metrics%count_mu(t,j,k) + 1                         
                               endif
                            enddo
                         endif
                      enddo
                   enddo
                endif
             enddo
          enddo

       !------------------------------------------------Y.Kwon
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

                                  if(lhour.eq.6) then
                                     metrics%sx_mu_6am(t,j,k) = metrics%sx_mu_6am(t,j,k) + &
                                          obs%value(t1,k)
                                     metrics%count_mu_6am(t,j,k) = &
                                          metrics%count_mu_6am(t,j,k) + 1
                                  elseif(lhour.eq.18) then
                                     metrics%sx_mu_6pm(t,j,k) = metrics%sx_mu_6pm(t,j,k) + &
                                          obs%value(t1,k)
                                     metrics%count_mu_6pm(t,j,k) = &
                                          metrics%count_mu_6pm(t,j,k) + 1
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
    
  end subroutine diagnoseSingleMu


!BOP
! 
! !ROUTINE: LDT_computeMu
! \label{LDT_computeMu}
! 
! !INTERFACE: 
  subroutine LDT_computeMu(n)
! !USES: 
    use LDT_coreMod,  only : LDT_rc
    use LDT_DAmetricsDataMod,        only : LDT_DAmetricsPtr

! 
! !DESCRIPTION: 
!   This subroutine issues the calls to compute Mu values for the 
!   desired variables
! 
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleMu](\ref{computeSingleMu})
!     updates the Mu computation for a single variable 
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
       call computeSingleMu(n,LDT_DAobsDataPtr(n,index)%dataEntryPtr,&
            LDT_DAmetricsPtr(index)%dataEntryPtr)
    enddo

  end subroutine LDT_ComputeMu

!BOP
! 
! !ROUTINE: computeSingleMu
! \label{computeSingleMu}
! 
! !INTERFACE: 
  subroutine computeSingleMu(n,obs, metrics)
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_domain
    use LDT_DAmetricsDataMod
! 
! !DESCRIPTION: 
!  This routine computes the Mu values for a single variable
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
    real, allocatable :: strat_mu(:,:,:)
    integer, allocatable :: strat_count(:,:,:)
    integer  :: t,i,j,k
    real     :: mu_temp_6am(LDT_rc%ngrid(n),24,obs%vlevels)   !Y.Kwon 
    real     :: mu_temp_6pm(LDT_rc%ngrid(n),24,obs%vlevels)   !Y.Kwon
    real     :: d_mu_temp_6am                                 !Y.Kwon
    real     :: d_mu_temp_6pm                                 !Y.Kwon
    integer  :: d_days                                        !Y.Kwon

    if(LDT_rc%endtime.eq.1) then 
       if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
          if(LDT_rc%group_cdfs.eq.0) then 
             do t=1,LDT_rc%ngrid(n)

                !for 9-km operational SMAP (Y.Kwon)
                if(LDT_rc%daily_interp_switch.eq.1) then
                   do j=1,24    
                      do k=1,obs%vlevels
                         if(metrics%count_mu_6am(t,j,k).le.LDT_rc%obsCountThreshold) then
                            mu_temp_6am(t,j,k) = LDT_rc%udef
                         else
                            mu_temp_6am(t,j,k) = &
                                 (metrics%sx_mu_6am(t,j,k)/metrics%count_mu_6am(t,j,k))
                         endif
                         if(metrics%count_mu_6pm(t,j,k).le.LDT_rc%obsCountThreshold) then
                            mu_temp_6pm(t,j,k) = LDT_rc%udef
                         else
                            mu_temp_6pm(t,j,k) = &
                                 (metrics%sx_mu_6pm(t,j,k)/metrics%count_mu_6pm(t,j,k))
                         endif
                      enddo
                   enddo

                   !daily interpolation using half-monthly values
                   do j=1,LDT_rc%cdf_ntimes    !365
                      do k=1,obs%vlevels
                         if(j<8) then
                            if(mu_temp_6am(t,1,k).gt.0.and.mu_temp_6am(t,24,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,1,k) - mu_temp_6am(t,24,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,24,k) + (d_mu_temp_6am/d_days) * (j+7)
                            else
                               if(mu_temp_6am(t,1,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,1,k)
                               elseif(mu_temp_6am(t,24,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,24,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,1,k).gt.0.and.mu_temp_6pm(t,24,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,1,k) - mu_temp_6pm(t,24,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,24,k) + (d_mu_temp_6pm/d_days) * (j+7)
                            else
                               if(mu_temp_6pm(t,1,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,1,k)
                               elseif(mu_temp_6pm(t,24,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,24,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==8) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,1,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,1,k)
                         elseif(j<24) then
                            if(mu_temp_6am(t,1,k).gt.0.and.mu_temp_6am(t,2,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,2,k) - mu_temp_6am(t,1,k)
                               d_days = 16
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,1,k) + (d_mu_temp_6am/d_days) * (j-8)
                            else
                               if(mu_temp_6am(t,1,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,1,k)
                               elseif(mu_temp_6am(t,2,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,2,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif 
                            endif
                            if(mu_temp_6pm(t,1,k).gt.0.and.mu_temp_6pm(t,2,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,2,k) - mu_temp_6pm(t,1,k)
                               d_days = 16
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,1,k) + (d_mu_temp_6pm/d_days) * (j-8)
                            else
                               if(mu_temp_6pm(t,1,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,1,k)
                               elseif(mu_temp_6pm(t,2,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,2,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==24) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,2,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,2,k)
                         elseif(j<39) then
                            if(mu_temp_6am(t,2,k).gt.0.and.mu_temp_6am(t,3,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,3,k) - mu_temp_6am(t,2,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,2,k) + (d_mu_temp_6am/d_days) * (j-24)
                            else
                               if(mu_temp_6am(t,2,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,2,k)
                               elseif(mu_temp_6am(t,3,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,3,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,2,k).gt.0.and.mu_temp_6pm(t,3,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,3,k) - mu_temp_6pm(t,2,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,2,k) + (d_mu_temp_6pm/d_days) * (j-24)
                            else
                               if(mu_temp_6pm(t,2,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,2,k)
                               elseif(mu_temp_6pm(t,3,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,3,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==39) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,3,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,3,k)
                         elseif(j<53) then
                            if(mu_temp_6am(t,3,k).gt.0.and.mu_temp_6am(t,4,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,4,k) - mu_temp_6am(t,3,k)
                               d_days = 14
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,3,k) + (d_mu_temp_6am/d_days) * (j-39)
                            else
                               if(mu_temp_6am(t,3,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,3,k)
                               elseif(mu_temp_6am(t,4,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,4,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,3,k).gt.0.and.mu_temp_6pm(t,4,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,4,k) - mu_temp_6pm(t,3,k)
                               d_days = 14
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,3,k) + (d_mu_temp_6pm/d_days) * (j-39)
                            else
                               if(mu_temp_6pm(t,3,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,3,k)
                               elseif(mu_temp_6pm(t,4,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,4,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==53) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,4,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,4,k)
                         elseif(j<67) then 
                            if(mu_temp_6am(t,4,k).gt.0.and.mu_temp_6am(t,5,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,5,k) - mu_temp_6am(t,4,k)
                               d_days = 14
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,4,k) + (d_mu_temp_6am/d_days) * (j-53)
                            else
                               if(mu_temp_6am(t,4,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,4,k)
                               elseif(mu_temp_6am(t,5,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,5,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,4,k).gt.0.and.mu_temp_6pm(t,5,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,5,k) - mu_temp_6pm(t,4,k)
                               d_days = 14
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,4,k) + (d_mu_temp_6pm/d_days) * (j-53)
                            else
                               if(mu_temp_6pm(t,4,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,4,k)
                               elseif(mu_temp_6pm(t,5,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,5,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==67) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,5,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,5,k)
                         elseif(j<83) then
                            if(mu_temp_6am(t,5,k).gt.0.and.mu_temp_6am(t,6,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,6,k) - mu_temp_6am(t,5,k)
                               d_days = 16
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,5,k) + (d_mu_temp_6am/d_days) * (j-67)
                            else
                               if(mu_temp_6am(t,5,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,5,k)
                               elseif(mu_temp_6am(t,6,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,6,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,5,k).gt.0.and.mu_temp_6pm(t,6,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,6,k) - mu_temp_6pm(t,5,k)
                               d_days = 16
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,5,k) + (d_mu_temp_6pm/d_days) * (j-67)
                            else
                               if(mu_temp_6pm(t,5,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,5,k)
                               elseif(mu_temp_6pm(t,6,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,6,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==83) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,6,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,6,k)
                         elseif(j<98) then
                            if(mu_temp_6am(t,6,k).gt.0.and.mu_temp_6am(t,7,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,7,k) - mu_temp_6am(t,6,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,6,k) + (d_mu_temp_6am/d_days) * (j-83)
                            else
                               if(mu_temp_6am(t,6,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,6,k)
                               elseif(mu_temp_6am(t,7,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,7,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,6,k).gt.0.and.mu_temp_6pm(t,7,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,7,k) - mu_temp_6pm(t,6,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,6,k) + (d_mu_temp_6pm/d_days) * (j-83)
                            else
                               if(mu_temp_6pm(t,6,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,6,k)
                               elseif(mu_temp_6pm(t,7,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,7,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==98) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,7,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,7,k)
                         elseif(j<113) then
                            if(mu_temp_6am(t,7,k).gt.0.and.mu_temp_6am(t,8,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,8,k) - mu_temp_6am(t,7,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,7,k) + (d_mu_temp_6am/d_days) * (j-98)
                            else
                               if(mu_temp_6am(t,7,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,7,k)
                               elseif(mu_temp_6am(t,8,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,8,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,7,k).gt.0.and.mu_temp_6pm(t,8,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,8,k) - mu_temp_6pm(t,7,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,7,k) + (d_mu_temp_6pm/d_days) * (j-98)
                            else
                               if(mu_temp_6pm(t,7,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,7,k)
                               elseif(mu_temp_6pm(t,8,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,8,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==113) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,8,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,8,k)
                         elseif(j<128) then
                            if(mu_temp_6am(t,8,k).gt.0.and.mu_temp_6am(t,9,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,9,k) - mu_temp_6am(t,8,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,8,k) + (d_mu_temp_6am/d_days) * (j-113)
                            else
                               if(mu_temp_6am(t,8,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,8,k)
                               elseif(mu_temp_6am(t,9,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,9,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,8,k).gt.0.and.mu_temp_6pm(t,9,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,9,k) - mu_temp_6pm(t,8,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,8,k) + (d_mu_temp_6pm/d_days) * (j-113)
                            else
                               if(mu_temp_6pm(t,8,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,8,k)
                               elseif(mu_temp_6pm(t,9,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,9,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==128) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,9,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,9,k)
                         elseif(j<144) then
                            if(mu_temp_6am(t,9,k).gt.0.and.mu_temp_6am(t,10,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,10,k) - mu_temp_6am(t,9,k)
                               d_days = 16
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,9,k) + (d_mu_temp_6am/d_days) * (j-128)
                            else
                               if(mu_temp_6am(t,9,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,9,k)
                               elseif(mu_temp_6am(t,10,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,10,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,9,k).gt.0.and.mu_temp_6pm(t,10,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,10,k) - mu_temp_6pm(t,9,k)
                               d_days = 16
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,9,k) + (d_mu_temp_6pm/d_days) * (j-128)
                            else
                               if(mu_temp_6pm(t,9,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,9,k)
                               elseif(mu_temp_6pm(t,10,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,10,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==144) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,10,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,10,k)
                         elseif(j<159) then
                            if(mu_temp_6am(t,10,k).gt.0.and.mu_temp_6am(t,11,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,11,k) - mu_temp_6am(t,10,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,10,k) + (d_mu_temp_6am/d_days) * (j-144)
                            else
                               if(mu_temp_6am(t,10,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,10,k)
                               elseif(mu_temp_6am(t,11,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,11,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,10,k).gt.0.and.mu_temp_6pm(t,11,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,11,k) - mu_temp_6pm(t,10,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,10,k) + (d_mu_temp_6pm/d_days) * (j-144)
                            else
                               if(mu_temp_6pm(t,10,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,10,k)
                               elseif(mu_temp_6pm(t,11,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,11,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==159) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,11,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,11,k)
                         elseif(j<174) then
                            if(mu_temp_6am(t,11,k).gt.0.and.mu_temp_6am(t,12,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,12,k) - mu_temp_6am(t,11,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,11,k) + (d_mu_temp_6am/d_days) * (j-159)
                            else
                               if(mu_temp_6am(t,11,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,11,k)
                               elseif(mu_temp_6am(t,12,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,12,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,11,k).gt.0.and.mu_temp_6pm(t,12,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,12,k) - mu_temp_6pm(t,11,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,11,k) + (d_mu_temp_6pm/d_days) * (j-159)
                            else
                               if(mu_temp_6pm(t,11,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,11,k)
                               elseif(mu_temp_6pm(t,12,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,12,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==174) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,12,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,12,k)
                         elseif(j<189) then
                            if(mu_temp_6am(t,12,k).gt.0.and.mu_temp_6am(t,13,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,13,k) - mu_temp_6am(t,12,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,12,k) + (d_mu_temp_6am/d_days) * (j-174)
                            else
                               if(mu_temp_6am(t,12,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,12,k)
                               elseif(mu_temp_6am(t,13,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,13,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,12,k).gt.0.and.mu_temp_6pm(t,13,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,13,k) - mu_temp_6pm(t,12,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,12,k) + (d_mu_temp_6pm/d_days) * (j-174)
                            else
                               if(mu_temp_6pm(t,12,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,12,k)
                               elseif(mu_temp_6pm(t,13,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,13,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==189) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,13,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,13,k)
                         elseif(j<205) then
                            if(mu_temp_6am(t,13,k).gt.0.and.mu_temp_6am(t,14,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,14,k) - mu_temp_6am(t,13,k)
                               d_days = 16
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,13,k) + (d_mu_temp_6am/d_days) * (j-189)
                            else
                               if(mu_temp_6am(t,13,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,13,k)
                               elseif(mu_temp_6am(t,14,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,14,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,13,k).gt.0.and.mu_temp_6pm(t,14,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,14,k) - mu_temp_6pm(t,13,k)
                               d_days = 16
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,13,k) + (d_mu_temp_6pm/d_days) * (j-189)
                            else
                               if(mu_temp_6pm(t,13,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,13,k)
                               elseif(mu_temp_6pm(t,14,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,14,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==205) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,14,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,14,k)
                         elseif(j<220) then
                            if(mu_temp_6am(t,14,k).gt.0.and.mu_temp_6am(t,15,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,15,k) - mu_temp_6am(t,14,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,14,k) + (d_mu_temp_6am/d_days) * (j-205)
                            else
                               if(mu_temp_6am(t,14,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,14,k)
                               elseif(mu_temp_6am(t,15,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,15,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,14,k).gt.0.and.mu_temp_6pm(t,15,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,15,k) - mu_temp_6pm(t,14,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,14,k) + (d_mu_temp_6pm/d_days) * (j-205)
                            else
                               if(mu_temp_6pm(t,14,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,14,k)
                               elseif(mu_temp_6pm(t,15,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,15,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==220) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,15,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,15,k)
                         elseif(j<236) then
                            if(mu_temp_6am(t,15,k).gt.0.and.mu_temp_6am(t,16,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,16,k) - mu_temp_6am(t,15,k)
                               d_days = 16
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,15,k) + (d_mu_temp_6am/d_days) * (j-220)
                            else
                               if(mu_temp_6am(t,15,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,15,k)
                               elseif(mu_temp_6am(t,16,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,16,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,15,k).gt.0.and.mu_temp_6pm(t,16,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,16,k) - mu_temp_6pm(t,15,k)
                               d_days = 16
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,15,k) + (d_mu_temp_6pm/d_days) * (j-220)
                            else
                               if(mu_temp_6pm(t,15,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,15,k)
                               elseif(mu_temp_6pm(t,16,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,16,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==236) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,16,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,16,k)
                         elseif(j<251) then
                            if(mu_temp_6am(t,16,k).gt.0.and.mu_temp_6am(t,17,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,17,k) - mu_temp_6am(t,16,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,16,k) + (d_mu_temp_6am/d_days) * (j-236)
                            else
                               if(mu_temp_6am(t,16,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,16,k)
                               elseif(mu_temp_6am(t,17,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,17,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,16,k).gt.0.and.mu_temp_6pm(t,17,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,17,k) - mu_temp_6pm(t,16,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,16,k) + (d_mu_temp_6pm/d_days) * (j-236)
                            else
                               if(mu_temp_6pm(t,16,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,16,k)
                               elseif(mu_temp_6pm(t,17,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,17,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==251) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,17,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,17,k)
                         elseif(j<266) then
                            if(mu_temp_6am(t,17,k).gt.0.and.mu_temp_6am(t,18,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,18,k) - mu_temp_6am(t,17,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,17,k) + (d_mu_temp_6am/d_days) * (j-251)
                            else
                               if(mu_temp_6am(t,17,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,17,k)
                               elseif(mu_temp_6am(t,18,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,18,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,17,k).gt.0.and.mu_temp_6pm(t,18,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,18,k) - mu_temp_6pm(t,17,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,17,k) + (d_mu_temp_6pm/d_days) * (j-251)
                            else
                               if(mu_temp_6pm(t,17,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,17,k)
                               elseif(mu_temp_6pm(t,18,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,18,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==266) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,18,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,18,k)
                         elseif(j<281) then
                            if(mu_temp_6am(t,18,k).gt.0.and.mu_temp_6am(t,19,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,19,k) - mu_temp_6am(t,18,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,18,k) + (d_mu_temp_6am/d_days) * (j-266)
                            else
                               if(mu_temp_6am(t,18,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,18,k)
                               elseif(mu_temp_6am(t,19,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,19,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,18,k).gt.0.and.mu_temp_6pm(t,19,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,19,k) - mu_temp_6pm(t,18,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,18,k) + (d_mu_temp_6pm/d_days) * (j-266)
                            else
                               if(mu_temp_6pm(t,18,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,18,k)
                               elseif(mu_temp_6pm(t,19,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,19,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==281) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,19,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,19,k)
                         elseif(j<297) then
                            if(mu_temp_6am(t,19,k).gt.0.and.mu_temp_6am(t,20,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,20,k) - mu_temp_6am(t,19,k)
                               d_days = 16
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,19,k) + (d_mu_temp_6am/d_days) * (j-281)
                            else
                               if(mu_temp_6am(t,19,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,19,k)
                               elseif(mu_temp_6am(t,20,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,20,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,19,k).gt.0.and.mu_temp_6pm(t,20,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,20,k) - mu_temp_6pm(t,19,k)
                               d_days = 16
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,19,k) + (d_mu_temp_6pm/d_days) * (j-281)
                            else
                               if(mu_temp_6pm(t,19,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,19,k)
                               elseif(mu_temp_6pm(t,20,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,20,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif 
                         elseif(j==297) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,20,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,20,k)
                         elseif(j<312) then
                            if(mu_temp_6am(t,20,k).gt.0.and.mu_temp_6am(t,21,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,21,k) - mu_temp_6am(t,20,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,20,k) + (d_mu_temp_6am/d_days) * (j-297)
                            else
                               if(mu_temp_6am(t,20,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,20,k)
                               elseif(mu_temp_6am(t,21,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,21,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,20,k).gt.0.and.mu_temp_6pm(t,21,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,21,k) - mu_temp_6pm(t,20,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,20,k) + (d_mu_temp_6pm/d_days) * (j-297)
                            else
                               if(mu_temp_6pm(t,20,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,20,k)
                               elseif(mu_temp_6pm(t,21,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,21,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==312) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,21,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,21,k)
                         elseif(j<327) then
                            if(mu_temp_6am(t,21,k).gt.0.and.mu_temp_6am(t,22,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,22,k) - mu_temp_6am(t,21,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,21,k) + (d_mu_temp_6am/d_days) * (j-312)
                            else
                               if(mu_temp_6am(t,21,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,21,k)
                               elseif(mu_temp_6am(t,22,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,22,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,21,k).gt.0.and.mu_temp_6pm(t,22,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,22,k) - mu_temp_6pm(t,21,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,21,k) + (d_mu_temp_6pm/d_days) * (j-312)
                            else
                               if(mu_temp_6pm(t,21,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,21,k)
                               elseif(mu_temp_6pm(t,22,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,22,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==327) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,22,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,22,k)
                         elseif(j<342) then
                            if(mu_temp_6am(t,22,k).gt.0.and.mu_temp_6am(t,23,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,23,k) - mu_temp_6am(t,22,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,22,k) + (d_mu_temp_6am/d_days) * (j-327)
                            else
                               if(mu_temp_6am(t,22,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,22,k)
                               elseif(mu_temp_6am(t,23,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,23,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,22,k).gt.0.and.mu_temp_6pm(t,23,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,23,k) - mu_temp_6pm(t,22,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,22,k) + (d_mu_temp_6pm/d_days) * (j-327)
                            else
                               if(mu_temp_6pm(t,22,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,22,k)
                               elseif(mu_temp_6pm(t,23,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,23,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==342) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,23,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,23,k)
                         elseif(j<358) then
                            if(mu_temp_6am(t,23,k).gt.0.and.mu_temp_6am(t,24,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,24,k) - mu_temp_6am(t,23,k)
                               d_days = 16
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,23,k) + (d_mu_temp_6am/d_days) * (j-342)
                            else
                               if(mu_temp_6am(t,23,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,23,k)
                               elseif(mu_temp_6am(t,24,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,24,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,23,k).gt.0.and.mu_temp_6pm(t,24,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,24,k) - mu_temp_6pm(t,23,k)
                               d_days = 16
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,23,k) + (d_mu_temp_6pm/d_days) * (j-342)
                            else
                               if(mu_temp_6pm(t,23,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,23,k)
                               elseif(mu_temp_6pm(t,24,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,24,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         elseif(j==358) then
                            metrics%mu_6am(t,j,k) = mu_temp_6am(t,24,k)
                            metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,24,k)
                         elseif(j<=365) then
                            if(mu_temp_6am(t,24,k).gt.0.and.mu_temp_6am(t,1,k).gt.0) then
                               d_mu_temp_6am = mu_temp_6am(t,1,k) - mu_temp_6am(t,24,k)
                               d_days = 15
                               metrics%mu_6am(t,j,k) = mu_temp_6am(t,24,k) + (d_mu_temp_6am/d_days) * (j-358)
                            else
                               if(mu_temp_6am(t,24,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,24,k)
                               elseif(mu_temp_6am(t,1,k).gt.0) then
                                  metrics%mu_6am(t,j,k) = mu_temp_6am(t,1,k)
                               else
                                  metrics%mu_6am(t,j,k) = LDT_rc%udef
                               endif
                            endif
                            if(mu_temp_6pm(t,24,k).gt.0.and.mu_temp_6pm(t,1,k).gt.0) then
                               d_mu_temp_6pm = mu_temp_6pm(t,1,k) - mu_temp_6pm(t,24,k)
                               d_days = 15
                               metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,24,k) + (d_mu_temp_6pm/d_days) * (j-358)
                            else
                               if(mu_temp_6pm(t,24,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,24,k)
                               elseif(mu_temp_6pm(t,1,k).gt.0) then
                                  metrics%mu_6pm(t,j,k) = mu_temp_6pm(t,1,k)
                               else
                                  metrics%mu_6pm(t,j,k) = LDT_rc%udef
                               endif
                            endif
                         endif
                      enddo
                   enddo

                else
                   do j=1,LDT_rc%cdf_ntimes
                      do k=1,obs%vlevels
                         if(metrics%count_mu(t,j,k).le.LDT_rc%obsCountThreshold) then
                            metrics%mu(t,j,k) = LDT_rc%udef
                         else
                            metrics%mu(t,j,k) = &
                                 (metrics%sx_mu(t,j,k)/metrics%count_mu(t,j,k))
                         endif
                      enddo
                   enddo
                endif
             enddo
          else
             allocate(strat_mu(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))
             allocate(strat_count(LDT_rc%group_cdfs_nbins, &
                  LDT_rc%cdf_ntimes, &
                  obs%vlevels))

             strat_mu    = 0 
             strat_count = 0 

              do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex = LDT_rc%cdf_strat_data(t)
                      strat_mu(sindex,j,k) = &
                           strat_mu(sindex,j,k) + metrics%sx_mu(t,j,k)
                      strat_count(sindex,j,k) = & 
                           strat_count(sindex,j,k) + metrics%count_mu(t,j,k)
                   enddo
                enddo
             enddo
             
             do t=1,LDT_rc%group_cdfs_nbins
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      if(strat_count(t,j,k).le.&
                           LDT_rc%obsCountThreshold) then 
                         strat_mu(t,j,k) = LDT_rc%udef
                      else
                         strat_mu(t,j,k) = strat_mu(t,j,k)/&
                              strat_count(t,j,k)
                      endif
                   enddo
                enddo
             enddo

             do t=1,LDT_rc%ngrid(n)
                do j=1,LDT_rc%cdf_ntimes
                   do k=1,obs%vlevels
                      sindex = LDT_rc%cdf_strat_data(t)
                      if(strat_mu(sindex,j,k).eq.LDT_rc%udef) then 
                         metrics%mu(t,j,k) = LDT_rc%udef
                      else
                         metrics%mu(t,j,k) = strat_mu(sindex,j,k)
                      endif
                   enddo
                enddo
             enddo

             deallocate(strat_mu)
             deallocate(strat_count)

          endif
       endif
    endif

  end subroutine computeSingleMu
end module LDT_MuMod
