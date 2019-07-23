!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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

    if(LDT_rc%endtime.eq.1) then 
       if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
          if(LDT_rc%group_cdfs.eq.0) then 
             do t=1,LDT_rc%ngrid(n)
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
