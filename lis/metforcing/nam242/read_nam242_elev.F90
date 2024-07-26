!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_nam242_elev
! \label{read_nam242_elev}
!
! !REVISION HISTORY:
!     Sep 2012: NOHRSC/NOAA: Initial specification
!
! !INTERFACE:
subroutine read_nam242_elev(n, findex, change)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_metforcingMod, only : LIS_forc
  use LIS_fileIOMod,     only : LIS_read_param
  use LIS_logMod,        only : LIS_logunit, LIS_getNextUnitNumber, &
                                LIS_releaseUnitNumber, LIS_endrun
  use nam242_forcingMod, only : nam242_struc
!  use gaussian_mod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: change
! !DESCRIPTION:
!
!  Opens, reads, and interpolates NAM model elevation to the LIS
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_readData](\ref{LIS_readData}) \newline
!    Abstract method to read the elevation of the forcing
!    data in the same map projection used in LIS.
!  \end{description}
!EOP
  integer :: c,r
  real :: go(LIS_rc%lnc(n),LIS_rc%lnr(n))

  if ( LIS_rc%met_ecor(findex).eq."lapse-rate") then 
     write(LIS_logunit,*) 'Reading the NAM elevation'

     call LIS_read_param(n,"ELEV_NAM242",go)
     
     if ( change == 0 ) then ! period 1980--1991
        ! Note that for this time period we have a difference file
        ! instead of a proper elevation file.  Correct this.
        !
        ! For now set LIS_domain%grid%elev = 0, and set the modelelev to the
        ! additive inverse of the difference data.
        ! We want to compute x - y, but we are reading (x-y) instead of y.
        go = -go
     endif

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              LIS_forc(n,findex)%modelelev(LIS_domain(n)%gindex(c,r)) = go(c,r)
           endif
        enddo
     enddo

  endif

end subroutine read_nam242_elev
