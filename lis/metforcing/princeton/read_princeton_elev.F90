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
! !ROUTINE: read_princeton_elev
! \label{read_princeton_elev}
!
! !REVISION HISTORY:
!
!  1 Feb 2007; Hiroko Kato; Initial Specificaton adopted from 
!                           read_ecmwf_elev.F90
!
! !INTERFACE:
subroutine read_princeton_elev(n,findex)
! !USES:
  use LIS_coreMod,         only : LIS_rc, LIS_domain
  use LIS_metforcingMod,   only : LIS_forc
  use LIS_fileIOMod,       only : LIS_read_param
  use LIS_logMod,          only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber
  use princeton_forcingMod,only : princeton_struc
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
! !DESCRIPTION:
!
!  Opens and reads PRINCETON model elevation to the LIS
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!  The elevation file needs to be preprocessed to fit the running 
!  resolution and domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!  The routines invoked are: 
!   \begin{description}
!   \item[LIS\_readData](\ref{LIS_readData}) \newline
!     call the abstract method to get the elevation data in the
!     LIS projection. 
!   \end{description}
!EOP

  real              :: go(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer           :: c,r
  integer           :: ftn

  if ( trim(LIS_rc%met_ecor(findex)) .ne. "none") then 

     write(LIS_logunit,*) 'Reading the PRINCETON elevation map ...'
     
!     ftn = LIS_getNextUnitNumber()
!     open(ftn,file=princeton_struc(n)%elevfile,form='unformatted', &
!          access='direct',recl=4)
!  
!     call LIS_readData(n,ftn,LIS_rc%topo_gridDesc(n,:),go)
!
!     call LIS_releaseUnitNumber(ftn)

     call LIS_read_param(n,"ELEV_PRINCETON",go)

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              LIS_forc(n,findex)%modelelev(LIS_domain(n)%gindex(c,r)) = go(c,r)
           endif
        enddo
     enddo

  endif

end subroutine read_princeton_elev
