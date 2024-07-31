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
! !ROUTINE: read_CONSTANT_elev
! \label{read_CONSTANT_elev}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  17 Jun 2014: KR Arsenault, added reader for constant elev 
!
! !INTERFACE:
subroutine read_CONSTANT_elev( n, num_bins, fgrd, elevave )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit,LDT_endrun

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: num_bins
  real,    intent(out):: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out):: elevave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine sets a constant elevation for basic tests
!   or hypothetical test cases.  
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[fgrd]
!   output field with the retrieved elevation data
!  \item[elevave]
!   output field with the retrieved elevation data
!  \end{description}
!EOP
! ____________________________

  fgrd    = LDT_rc%udef    ! Set as elevation water type
  elevave = LDT_rc%udef    ! Set as elevation water type

  write(LDT_logunit, *) "[INFO] Setting a CONSTANT elevation value ..."
  

end subroutine read_CONSTANT_elev
