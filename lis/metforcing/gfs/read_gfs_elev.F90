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
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_metforcingMod, only : LIS_forc
  use LIS_fileIOMod,      only : LIS_readData
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber, LIS_endrun
  use gfs_forcingMod,    only : gfs_struc
!  use gaussian_mod

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: change

! !DESCRIPTION:
!
!  Opens, reads, and interpolates GFS model elevation to the LIS
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
  integer :: ftn

#if 0 
  if ( LIS_rc%met_ecor(findex) .ge. 1 ) then 

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
        write(LIS_logunit,*) 'ERR: invalid update request. ', change, &
                         'is not in {1,2,3,4,5,6}'
        call LIS_endrun
     endif

     write(LIS_logunit,*) 'Reading the GFS elevation ',trim(gfs_struc(n)%elevfile)
     
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=trim(gfs_struc(n)%elevfile),&
          form='unformatted',access='direct',recl=4)
     
     call LIS_readData(n,ftn,LIS_rc%topo_gridDesc(n,:),go)

     call LIS_releaseUnitNumber(ftn)

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
#endif
end subroutine read_gfs_elev

