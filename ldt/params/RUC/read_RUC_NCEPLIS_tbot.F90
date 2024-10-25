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
! !ROUTINE: read_RUC_NCEPLIS_tbot
! \label{read_RUC_NCEPLIS_tbot}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_RUC_NCEPLIS_tbot(n, array)

! !USES:
  use ESMF
  use LDT_coreMod,     only : LDT_rc, LDT_domain
  use LDT_logMod,      only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use LDT_fileIOMod,   only : readLISdata 
  use RUC_parmsMod

  implicit none

! !ARGUMENTS: 
  integer,   intent(in) :: n
  real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine retrieves the bottom temperature climatology for the 
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

  inquire(file=trim(RUC_struc(n)%tbotFile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "RUC_NCEPLIS TBOT map ",trim(RUC_struc(n)%tbotFile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif

  write(LDT_logunit,*) "MSG: Reading LIS-based RUC_NCEPLIS Bottom Temperature file: ",&
        trim(RUC_struc(n)%tbotfile)

  ftn = LDT_getNextUnitNumber()
  open(ftn, file=RUC_struc(n)%tbotFile, access='direct',status='old', &
       form="unformatted", recl=4)
  
  call readLISdata(n, ftn, RUC_struc(n)%tbot_proj,  &
                   RUC_struc(n)%tbot_gridtransform, &
                   RUC_struc(n)%tbot_gridDesc(:), 1, array)

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(array(c,r,1).lt.0) then
           array(c,r,1) = LDT_rc%udef
        endif
     enddo
  enddo
  call LDT_releaseUnitNumber(ftn)

end subroutine read_RUC_NCEPLIS_tbot
