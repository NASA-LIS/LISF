!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_ecmwfreanal_elev
!  \label{read_ecmwfreanal_elev}
!
! !REVISION HISTORY:
!  17Dec2004; Sujay Kumar; Initial Specificaton
!  3  Dec 2007; Sujay Kumar; Added the abstract method to read the 
!               EDAS data.  
!
! !INTERFACE:
subroutine read_ecmwfreanal_elev(n,findex)
! !USES:
  use LIS_coreMod,            only : LIS_rc, LIS_domain
  use LIS_metforcingMod,     only : LIS_forc
  use LIS_logMod,             only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber
  use LIS_fileIOMod,          only : LIS_read_param
  use ecmwfreanal_forcingMod, only : ecmwfreanal_struc

!EOP
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates ECMWF reanalysis model elevation. 
!  The data will be used to perform any topographical adjustments 
!  to the forcing. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!  The routines invoked are: 
!   \begin{description}
!    \item[readforcingelev](\ref{readforcingelev}) \newline
!     abstract method to read EDAS elevation data in the map 
!     projection used in the LIS model grid. 
!   \end{description}
!EOP
  integer :: c,r
  integer :: ftn
  real :: go(LIS_rc%lnc(n),LIS_rc%lnr(n))

  allocate(LIS_forc(n,findex)%modelelev(LIS_rc%ngrid(n)))
  
  write(LIS_logunit,*) 'Reading the ECMWF Reanalysis elevation ',&
       ecmwfreanal_struc(n)%elevfile
!  ftn = LIS_getNextUnitNumber()
!  open(ftn,file=ecmwfreanal_struc(n)%elevfile,form='unformatted',&
!       access='direct',recl=4)

!  call LIS_readData(n,ftn,LIS_rc%topo_gridDesc(n,:),go)
  
!  call LIS_releaseUnitNumber(ftn)
  call LIS_read_param(n,"ELEV_ECMWF reanalysis",go)

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           LIS_forc(n,findex)%modelelev(LIS_domain(n)%gindex(c,r)) = go(c,r)
        endif
     enddo
  enddo

end subroutine read_ecmwfreanal_elev
