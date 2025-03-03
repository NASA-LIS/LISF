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
! !ROUTINE: read_CONSTANT_soilfractions
! \label{read_CONSTANT_soilfractions}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  17 Jun 2014: KR Arsenault, added reader for constant soil fractions 
!
! !INTERFACE:
subroutine read_CONSTANT_soilfractions(n, num_bins, soilsfgrd, &
                                       sandave, clayave, siltave )

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit,LDT_endrun

  implicit none
! !ARGUMENTS: 
  integer, intent(in)   :: n          ! nest index
  integer, intent(in)   :: num_bins   ! number of bins for tiling
  real,    intent(out)  :: soilsfgrd(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: sandave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: clayave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)
  real,    intent(out)  :: siltave(LDT_rc%lnc(n),LDT_rc%lnr(n),num_bins)

! !DESCRIPTION:
!  This subroutine sets a constant soil fractions for basic tests
!   or hypothetical test cases.  
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[array]
!   output field with the retrieved soil fractions data
!  \end{description}
!EOP
! ____________________________

  soilsfgrd = 0.        
  sandave   = LDT_rc%udef   ! Set as soil fractions water index value
  clayave   = LDT_rc%udef
  siltave   = LDT_rc%udef

  write(LDT_logunit, *) "[INFO] Setting a CONSTANT soil fractions value ..."
  

end subroutine read_CONSTANT_soilfractions
