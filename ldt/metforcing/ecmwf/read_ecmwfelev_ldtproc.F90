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
! !ROUTINE: read_ecmwfelev_ldtproc
! \label{read_ecmwfelev_ldtproc}
!
! !REVISION HISTORY:
!
!  17Dec2004; Sujay Kumar; Initial Specificaton
!
! !INTERFACE:
subroutine read_ecmwfelev_ldtproc(n, findex, change)
! !USES:
  use LDT_coreMod,        only : LDT_rc, LDT_domain
  use LDT_metforcingMod,  only : LDT_forc
  use LDT_fileIOMod,      only : LDT_read_param
  use LDT_logMod,         only : LDT_logunit, LDT_endrun
  use ecmwf_forcingMod,   only : ecmwf_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: change

! !DESCRIPTION:
!
!  Opens, reads, and interpolates ECMWF model elevation to the LDT
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
  real    :: go(LDT_rc%lnc(n),LDT_rc%lnr(n))
  character(80) :: filename

  if ( trim(LDT_rc%met_ecor(findex)) .ne."none" ) then 
     write(LDT_logunit,*) " Reading the ECMWF elevation: ",trim(filename)
     if ( change == 0 ) then ! period 2001-2002
        call LDT_read_param(n,"ECMWF S23R4 ELEVATION",go)
     elseif ( change == 1 ) then ! period 04/2002-02/2006
        call LDT_read_param(n,"ECMWF S25R1 ELEVATION",go)
     elseif ( change == 2 ) then ! period 02/2006-06/2008
        call LDT_read_param(n,"ECMWF S30R1 ELEVATION",go)
     elseif ( change == 3 ) then ! period 06/2008-03/2009
        call LDT_read_param(n,"ECMWF S33R1 ELEVATION",go)
     elseif ( change == 4 ) then ! period 03/2009-10/2009
        call LDT_read_param(n,"ECMWF S35R2 ELEVATION",go)
     elseif ( change == 5 ) then ! period 10/2009-01/2010
        call LDT_read_param(n,"ECMWF S35R3 ELEVATION",go)
     elseif ( change == 6 ) then ! period 01/2010-05/2011
        call LDT_read_param(n,"ECMWF S36R1 ELEVATION",go)
     elseif ( change == 7 ) then ! period 05/2011 onwards
        call LDT_read_param(n,"ECMWF S37R2 ELEVATION",go)
     else
        write(LDT_logunit,*) 'ERR: invalid update request. ', change, &
                         'is not in {0,1,2}'
        call LDT_endrun
     endif

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(LDT_domain(n)%gindex(c,r).ne.-1) then 
              LDT_forc(n,findex)%modelelev(LDT_domain(n)%gindex(c,r)) = go(c,r)
           endif
        enddo
     enddo

  endif

end subroutine read_ecmwfelev_ldtproc
