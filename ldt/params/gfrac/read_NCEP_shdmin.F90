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
! !ROUTINE: read_NCEP_shdmin
! \label{read_NCEP_shdmin}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!
! !INTERFACE:
subroutine read_NCEP_shdmin(n, array)

! !USES:
  use ESMF
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_logMod,        only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,     only : readLISdata 
  use LDT_gfracMod

  implicit none
! !ARGUMENTS: 
  integer,   intent(in) :: n
  real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves the greenness fraction climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved greenness fraction
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c,r
  logical :: file_exists
! ____________________________

  inquire(file=trim(LDT_gfrac_struc(n)%shdminfile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) 'SHDMIN map ',trim(LDT_gfrac_struc(n)%shdminfile),' not found'
     write(LDT_logunit,*) 'Program stopping ...'
     call LDT_endrun
  endif
  write(LDT_logunit,*) "[INFO] Reading NCEP-LIS greenness min file: ",&
        trim(LDT_gfrac_struc(n)%shdminfile)

  ftn = LDT_getNextUnitNumber()
  open(ftn, file=trim(LDT_gfrac_struc(n)%shdminfile), access='direct',status='old', &
       form="unformatted", recl=4)
  
  call readLISdata(n, ftn, LDT_gfrac_struc(n)%gfrac_proj, &
       LDT_gfrac_struc(n)%gfrac_gridtransform, &
       LDT_gfrac_struc(n)%gfrac_gridDesc(:), 1, array )  ! 1 indicates 2D layer

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(array(c,r,1).lt.0) then
           array(c,r,1) = LDT_rc%udef
        endif
     enddo
  enddo
  call LDT_releaseUnitNumber(ftn)

end subroutine read_NCEP_shdmin
