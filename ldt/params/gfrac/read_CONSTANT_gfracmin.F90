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
! !ROUTINE: read_CONSTANT_gfracmin
! \label{read_CONSTANT_gfracmin}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  17 Jun 2014: KR Arsenault, added reader for constant greenness 
!
! !INTERFACE:
subroutine read_CONSTANT_gfracmin(n, array)

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit,LDT_endrun

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  real, intent(inout) :: array(LDT_rc%lnc(n),LDT_rc%lnr(n))

! !DESCRIPTION:
!  This subroutine sets a constant greenness fraction minimum for basic tests
!   or hypothetical test cases.  
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved greenness fraction minimum data
!  \end{description}
!EOP
! ____________________________

  array = LDT_rc%udef    ! Set as gfracmin water index value 

  write(LDT_logunit, *) "[INFO] Setting a CONSTANT greenness fraction min. value ..."
  

end subroutine read_CONSTANT_gfracmin
