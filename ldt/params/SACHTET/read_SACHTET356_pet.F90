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
! !ROUTINE: read_SACHTET356_pet
! \label{read_SACHTET356_pet}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Nov 2013: KR Arsenault; Added PET data support for SACHTET
!
! !INTERFACE:
subroutine read_SACHTET356_pet(n, array)

! !USES:
  use ESMF
  use LDT_coreMod,  only : LDT_rc, LDT_domain
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun
  use LDT_xmrg_reader
  use SACHTET_parmsMod

  implicit none

! !ARGUMENTS: 
  integer,   intent(in) :: n
  real,   intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine retrieves the leaf area index (PET) climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[mo]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved PET 
!  \end{description}
!
!EOP      
  integer :: ftn
  integer :: c,r
  logical :: file_exists

  array = LDT_rc%udef

  inquire(file=trim(SACHTET_struc(n)%petFile), exist=file_exists)
  if(.not.file_exists) then 
     write(LDT_logunit,*) "PET map ",trim(SACHTET_struc(n)%petFile)," not found."
     write(LDT_logunit,*) "Program stopping ..."
     call LDT_endrun
  endif
  select case ( SACHTET_struc(n)%pet_gridtransform )
    case( "none" )
      write(LDT_logunit,*) "[INFO] Reading SAC-HTET climatological PET monthly files "
    case default
      write(LDT_logunit,*) "[ERR] The spatial transform option selected for SAC-HTET "
      write(LDT_logunit,*) "     clim PPT monthly file is not recognized nor recommended."
      write(LDT_logunit,*) "     Please select for now:  none "
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
  end select

  call LDT_transform_xmrgparam( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
           LDT_rc%gridDesc(n,:), SACHTET_struc(n)%petFile, LDT_rc%udef, array )


end subroutine read_SACHTET356_pet
