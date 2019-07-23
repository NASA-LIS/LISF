!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
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

    if(LDT_rc%endtime.eq.1) then 
       if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
          if(LDT_rc%group_cdfs.eq.0) then 
             do t=1,LDT_rc%ngrid(n)
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
