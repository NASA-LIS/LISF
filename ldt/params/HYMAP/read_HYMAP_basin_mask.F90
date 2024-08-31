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
! !ROUTINE: read_HYMAP_basin_mask
! \label{read_HYMAP_basin_mask}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!
! !INTERFACE:
subroutine read_HYMAP_basin_mask(n, array)
! !USES:
  use ESMF
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,     only : readLISdata
  use HYMAP_parmsMod

  implicit none
! !ARGUMENTS: 
  integer,     intent(in) :: n
  real,     intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves the bottom temperature climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[mo]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved greenness fraction
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c,r
  logical :: file_exists
  integer :: iarray(LDT_rc%lnc(n),LDT_rc%lnr(n),1) ! precision of binary file is int32

  inquire(file=trim(HYMAP_struc(n)%basinmaskfile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) '[ERR] HYMAP basin mask map, ',&
           trim(HYMAP_struc(n)%basinmaskfile),', not found.'
     write(LDT_logunit,*) 'Program stopping ...'
     call LDT_endrun
  endif

  ftn = LDT_getNextUnitNumber()
  open(ftn, file=trim(HYMAP_struc(n)%basinmaskfile), access='direct',&
       status='old', form="unformatted", convert="big_endian", recl=4)
  
  call readLISdata(n, ftn, HYMAP_struc(n)%hymap_proj, &
       HYMAP_struc(n)%hymap_gridtransform, &
       HYMAP_struc(n)%hymapparms_gridDesc(:), 1, iarray)
  
  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(iarray(c,r,1).lt.0) then
           array(c,r,1) = LDT_rc%udef
        else
           array(c,r,1) = iarray(c,r,1)
        endif
     enddo
  enddo
  call LDT_releaseUnitNumber(ftn)

end subroutine read_HYMAP_basin_mask
