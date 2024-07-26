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
! !ROUTINE: read_GTOPO30_curv
! \label{read_GTOPO30_curv}
!
! !REVISION HISTORY:
!  22 Dec 2004: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_GTOPO30_curv(n,curv)
! !USES:
  use LDT_coreMod,     only : LDT_rc
  use LDT_logMod, only      : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,   only : readLISdata 

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  real                :: curv(LDT_rc%lnc(n), LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves static curvature data from the GTOPO30 source
!  and reprojects it to the latlon projection. 
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[curv]
!   output field with the retrieved curvature 
!  \end{description}
!EOP

  integer :: ftn
  logical :: file_exists

  write(LDT_logunit,*) "[INFO] Reading GTOPO30-LIS curvature file ...",&
                       trim(LDT_rc%curvfile(n))
  
  inquire(file=trim(LDT_rc%curvfile(n)), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "Curvature map ",trim(LDT_rc%curvfile(n))," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  
  ftn = LDT_getNextUnitNumber()
  open(ftn,file=LDT_rc%curvfile(n),form='unformatted', &
       access='direct',convert='big_endian',recl=4,status='old')

  call readLISdata(n, ftn, LDT_rc%topo_proj, LDT_rc%topo_gridtransform(n), &
                   LDT_rc%topo_gridDesc(n,:), 1, curv )

  call LDT_releaseUnitNumber(ftn)
  write(LDT_logunit, *) "[INFO] read_GTOPO30_curv -- Done reading curvature file"

end subroutine read_GTOPO30_curv

