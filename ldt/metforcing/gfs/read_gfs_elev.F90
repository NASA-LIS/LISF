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
! !ROUTINE: read_gfs_elev
! \label{read_gfs_elev}
!
! !REVISION HISTORY:
!  16 Mar 2008: Sujay Kumar; Initial specification
!
! !INTERFACE:
subroutine read_gfs_elev(n, findex, change)

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_metforcingMod, only : LDT_forc
!  use LDT_fileIOMod,     only : LDT_readData
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use gfs_forcingMod,    only : gfs_struc

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: change

! !DESCRIPTION:
!
!  Opens, reads, and interpolates GFS model elevation to the LDT
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
!  \item[LDT\_readData](\ref{LDT_readData}) \newline
!    Abstract method to read the elevation of the forcing
!    data in the same map projection used in LDT.
!  \end{description}
!EOP
  integer :: c,r
  real    :: go(LDT_rc%lnc(n),LDT_rc%lnr(n))
  integer :: ftn

#if 0 
  if ( LDT_rc%met_ecor(findex) .ge. 1 ) then 

     if ( change == 0 ) then ! period 1980--1991
        ! Note that for this time period we have a difference file
        ! instead of a proper elevation file.  Correct this.
        gfs_struc(n)%elevfile = trim(gfs_struc(n)%t126elevfile)
     elseif ( change == 1 ) then ! period 1991--2000
        gfs_struc(n)%elevfile = trim(gfs_struc(n)%t126elevfile)
     elseif ( change == 2 ) then ! period 2000--2002
        gfs_struc(n)%elevfile = trim(gfs_struc(n)%t170elevfile)
     elseif ( change == 3 ) then ! period 2002--2005
        gfs_struc(n)%elevfile = trim(gfs_struc(n)%t254elevfile)
     elseif ( change == 4 ) then ! period 2005--2010
        gfs_struc(n)%elevfile = trim(gfs_struc(n)%t382elevfile)
     elseif ( change == 5 ) then ! period 2010--2015
        gfs_struc(n)%elevfile = trim(gfs_struc(n)%t574elevfile)
     elseif ( change == 6 ) then ! period 2015--
        gfs_struc(n)%elevfile = trim(gfs_struc(n)%t1534elevfile)
     else
        write(LDT_logunit,*) 'ERR: invalid update request. ', change, &
                         'is not in {1,2,3,4}'
        call LDT_endrun
     endif

     write(LDT_logunit,*) 'Reading the GFS elevation ',trim(gfs_struc(n)%elevfile)
     
     ftn = LDT_getNextUnitNumber()
     open(ftn,file=trim(gfs_struc(n)%elevfile),&
          form='unformatted',access='direct',recl=4)
     
!     call LDT_readData(n,ftn,LDT_rc%topo_gridDesc(n,:),go)

     call LDT_releaseUnitNumber(ftn)

     if ( change == 0 ) then ! period 1980--1991
        ! Note that for this time period we have a difference file
        ! instead of a proper elevation file.  Correct this.
        !
        ! For now set LDT_domain%grid%elev = 0, and set the modelelev to the
        ! additive inverse of the difference data.
        ! We want to compute x - y, but we are reading (x-y) instead of y.
        go = -go
     endif

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(LDT_domain(n)%gindex(c,r).ne.-1) then 
              LDT_forc(n,findex)%modelelev(LDT_domain(n)%gindex(c,r)) = go(c,r)
           endif
        enddo
     enddo

  endif
#endif

end subroutine read_gfs_elev
