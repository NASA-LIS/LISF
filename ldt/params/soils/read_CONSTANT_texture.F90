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
! !ROUTINE: read_CONSTANT_texture
! \label{read_CONSTANT_texture}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  17 Jun 2014: KR Arsenault, added reader for constant texture 
!
! !INTERFACE:
subroutine read_CONSTANT_texture(n, num_bins, fgrd, texture_layers )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit,LDT_endrun

  implicit none
! !ARGUMENTS: 
  integer, intent(in)   :: n
  integer, intent(in)   :: num_bins   ! Number of soil types
  real,   intent(inout) :: fgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,   intent(inout) :: texture_layers(LDT_rc%lnc(n),LDT_rc%lnr(n),1)

! !DESCRIPTION:
!  This subroutine sets a constant soil texture for basic tests
!   or hypothetical test cases.  
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[fgrd]
!   output field with the retrieved soil texture data
!  \end{description}
!EOP
! ____________________________

  fgrd = 0.
  texture_layers = 0.

  write(LDT_logunit, *) "[INFO] Setting a CONSTANT soil texture value (based on STATSGO range) ..."

  fgrd(:,:,14) = 1.    ! Set as soil texture water index value = 14 for STATSGO
  

end subroutine read_CONSTANT_texture
