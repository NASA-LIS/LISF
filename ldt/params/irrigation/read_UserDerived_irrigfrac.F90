!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_UserDerived_irrigfrac
!  \label{read_UserDerived_irrigfrac}
!
! !REVISION HISTORY:
!  03 Feb 2020: K. Arsenault; Added Irrigation UserDerived map reader
!
! !INTERFACE:
 subroutine read_UserDerived_irrigfrac(n, fgrd) 

! !USES:
  use LDT_coreMod, only : LDT_rc, LDT_domain
  use LDT_logMod,  only : LDT_logunit, LDT_verify, &
                          LDT_getNextUnitNumber, &
                          LDT_releaseUnitNumber, LDT_endrun
  use LDT_gridmappingMod
  use LDT_fileIOMod,  only : readLISdata
  use LDT_irrigationMod

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION:
!
!  User Derived Irrigation fraction map
!  -- Irrigation fraction map derived by the user
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[fgrd]
!   irrigation fraction (user-derived and readin)
!  \end{description}
!
!EOP      

  integer   :: nc, nr, i, c, r, j, t
  integer   :: ftn, line
  integer   :: isum
  logical   :: file_exists
  integer   :: mi                     ! Total number of input param grid array points
  integer   :: mo                     ! Total number of output LIS grid array points
  integer   :: glpnc, glpnr           ! Parameter (global) total columns and rows
  real      :: param_gridDesc(20)     ! Input parameter grid desc fgrd
  integer   :: noncrop 
  real      :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! ______________________________________________________________

   noncrop = 255

!- Check if file is present:
   inquire(file=trim(LDT_irrig_struc(n)%irrigfracfile), exist=file_exists)
   if(.not.file_exists) then
      write(LDT_logunit,*) "[ERR] Irrigation fraction map, ",&
                          trim(LDT_irrig_struc(n)%irrigfracfile),", not found."
      call LDT_endrun
   endif
   write(LDT_logunit,*) "[INFO] Reading UserDerived irrigation area map"

!- Open the file:
   ftn = LDT_getNextUnitNumber()
   open( ftn, file = LDT_irrig_struc(n)%irrigfracfile,form='unformatted', &
!          access='direct', convert='little_endian', status='old', &
          access='direct', status='old', &
          recl=4 )

  call readLISdata(n, ftn, LDT_irrig_struc(n)%irrig_proj, &
       LDT_irrig_struc(n)%irrigfrac_gridtransform, &
       LDT_irrig_struc(n)%irrig_gridDesc(:), 1, array)  ! 1 indicates 2D layer

  do r = 1,LDT_rc%lnr(n)
     do c = 1,LDT_rc%lnc(n)
        if( array(c,r,1).lt.0 ) then
           fgrd(c,r) = 0.
        else
           fgrd(c,r) = array(c,r,1)
        endif
     enddo
  enddo
!  fgrd(:,:) = array(:,:,1)

  call LDT_releaseUnitNumber(ftn)

  write(LDT_logunit,*) "[INFO] Done reading User derived irrigation gridcell fractions"

end subroutine read_UserDerived_irrigfrac

