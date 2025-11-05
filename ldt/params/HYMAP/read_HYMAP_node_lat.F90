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
!
! !ROUTINE: read_HYMAP_node_lat
! \label{read_HYMAP_node_lat}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!
! !INTERFACE:
subroutine read_HYMAP_node_lat(n, array)

! !USES:
  use ESMF
  use HYMAP_parmsMod
  use LDT_coreMod,      only : LDT_rc
  use LDT_fileIOMod,    only : readLISdata
  use LDT_logMod,       only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun

  implicit none

! !ARGUMENTS:

  integer,          intent(in) :: n
  real,          intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

  integer :: ftn
  logical :: file_exists

  ftn = LDT_getNextUnitNumber()

  inquire(file=trim(HYMAP_struc(n)%nodelatfile), exist=file_exists)
  if(.not.file_exists) then
     write(LDT_logunit,*) '[ERR] HYMAP node latitude map, ',&
           trim(HYMAP_struc(n)%nodelatfile),', not found.'
     write(LDT_logunit,*) 'Program stopping ...'
     call LDT_endrun
  endif

  open(ftn, file=trim(HYMAP_struc(n)%nodelatfile), access='direct',&
       status='old', form="unformatted", convert="big_endian", recl=4)

  call readLISdata(n, ftn, HYMAP_struc(n)%hymap_proj, &
       HYMAP_struc(n)%hymap_gridtransform, &
       HYMAP_struc(n)%hymapparms_gridDesc(:), 1, array)  ! 1 = 2D layer

  call LDT_releaseUnitNumber(ftn)

end subroutine read_HYMAP_node_lat
