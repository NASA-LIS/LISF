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
! !ROUTINE: cable_setup
! \label{cable_setup}
!
! !REVISION HISTORY:
!  21 Jul 2004: Sujay Kumar, Initial Specification
!  23 Feb 2007: Kristi Arsenault, Updated for LISv5.0
!  25 Jul 2011: David Mocko, CABLE LSM implementation in LISv6.+
!  23 May 2013: David Mocko, latest CABLE v1.4b version for LIS6.2
!
! !INTERFACE:
subroutine cable_setup
! !USES:
  use LIS_coreMod,        only : LIS_rc
  use cable_arrays
  use cable_lsmMod
!
! !DESCRIPTION:
!  This subroutine is the entry point to set up the parameters
!  required for the CABLE LSM.  Right now, it includes just the
!  vegetation and soil parameters.
!
!  The routines invoked are:
!  \begin{description}
!   \item[cable\_setsoilparms](\ref{cable_setsoilparms}) \newline
!    Initializes the soil parameters in CABLE
!   \item[cable\_setvegparms](\ref{cable_setvegparms}) \newline
!    Initializes the vegetation parameters in CABLE
!  \end{description}
!EOP
  implicit none
  
  integer :: n
  
! Must call setsoilparms first because a vegetation parameter
! in setvegparms is a function of the soil layer thicknesses
  call cable_setsoilparms
  call cable_setvegparms

  call alloc_cbm_var(model_structure,1)
  
  do n = 1,LIS_rc%nnest
!ccc         cable_struc(n)%ktau = 1
     model_structure%canopy = cable_struc(n)%canopyflag
     model_structure%photosynthesis = cable_struc(n)%photosynflag
     model_structure%soil = cable_struc(n)%soilflag
     if (any(model_structure%soil.eq.'sli')) then
        model_structure%sli_litter = cable_struc(n)%slilitterflag
        model_structure%sli_isotope = cable_struc(n)%sliisotopeflag
        model_structure%sli_coupled = cable_struc(n)%slicoupledflag
     endif
  enddo

! Allocate CABLE's main variables:
  call allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,       &
       soil,ssoil,sum_flux,veg,model_structure,1)
  
end subroutine cable_setup
