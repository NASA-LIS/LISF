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
! !ROUTINE: read_AVHRR_lai
! \label{read_AVHRR_lai}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!
! !INTERFACE:
subroutine read_AVHRR_lai(n, array)

! !USES:
  use ESMF
  use LDT_coreMod,     only : LDT_rc, LDT_domain
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
       LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,   only : readLISdata 
  use LDT_laisaiMod

  implicit none
! !ARGUMENTS: 

  integer,          intent(in) :: n
  real,          intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves the leaf area index (LAI) climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[mo]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved LAI 
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c,r
  logical :: file_exists
! __________________________

  inquire(file=trim(LDT_laisai_struc(n)%laiFile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "LAI map ",trim(LDT_laisai_struc(n)%laiFile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

  ftn = LDT_getNextUnitNumber()
  open(ftn, file=trim(LDT_laisai_struc(n)%laiFile), access='direct',status='old', &
       form="unformatted", recl=4)
  
  call readLISdata(n, ftn, LDT_laisai_struc(n)%laisai_proj, &
       LDT_laisai_struc(n)%laisai_gridtransform, &
       LDT_laisai_struc(n)%laisai_gridDesc(:), 1, array)  ! 1 indicates 2D layer

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(array(c,r,1).lt.0) then
           array(c,r,1) = LDT_rc%udef
        endif
     enddo
  enddo
  call LDT_releaseUnitNumber(ftn)

end subroutine read_AVHRR_lai
