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
! !ROUTINE: read_ecmwf_elev
! \label{read_ecmwf_elev}
!
! !REVISION HISTORY:
!
!  17Dec2004; Sujay Kumar; Initial Specificaton
!
! !INTERFACE:
subroutine read_ecmwf_elev(n, findex, change)
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain
  use LIS_metforcingMod, only : LIS_forc
  use LIS_fileIOMod,     only : LIS_read_param
  use LIS_logMod,         only : LIS_logunit, LIS_endrun
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use ecmwf_forcingMod,   only : ecmwf_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: change

! !DESCRIPTION:
!
!  Opens, reads, and interpolates ECMWF model elevation to the LIS
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!EOP
  integer :: c,r,line1,line2,nc_dom,line
  integer :: glnc, glnr
  real :: go(LIS_rc%lnc(n),LIS_rc%lnr(n))
  character(len=LIS_CONST_PATH_LEN) :: filename

  if ( trim(LIS_rc%met_ecor(findex)) .ne."none" ) then 
     write(LIS_logunit,*) 'Reading the ECMWF elevation ',trim(filename)
     if ( change == 0 ) then ! period 2001-2002
        call LIS_read_param(n,"ELEV_ECMWF_S23R4",go)
     elseif ( change == 1 ) then ! period 04/2002-02/2006
        call LIS_read_param(n,"ELEV_ECMWF_S25R1",go)
     elseif ( change == 2 ) then ! period 02/2006-06/2008
        call LIS_read_param(n,"ELEV_ECMWF_S30R1",go)
     elseif ( change == 3 ) then ! period 06/2008-03/2009
        call LIS_read_param(n,"ELEV_ECMWF_S33R1",go)
     elseif ( change == 4 ) then ! period 03/2009-10/2009
        call LIS_read_param(n,"ELEV_ECMWF_S35R2",go)
     elseif ( change == 5 ) then ! period 10/2009-01/2010
        call LIS_read_param(n,"ELEV_ECMWF_S35R3",go)
     elseif ( change == 6 ) then ! period 01/2010-05/2011
        call LIS_read_param(n,"ELEV_ECMWF_S36R1",go)
     elseif ( change == 7 ) then ! period 05/2011 onwards
        call LIS_read_param(n,"ELEV_ECMWF_S37R2",go)
     else
        write(LIS_logunit,*) 'ERR: invalid update request. ', change, &
                         'is not in {0,1,2}'
        call LIS_endrun
     endif


     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              LIS_forc(n,findex)%modelelev(LIS_domain(n)%gindex(c,r)) = go(c,r)
           endif
        enddo
     enddo

  endif

end subroutine read_ecmwf_elev
