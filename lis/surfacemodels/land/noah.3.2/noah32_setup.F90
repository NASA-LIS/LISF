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
! !ROUTINE: noah32_setup
! \label{noah32_setup}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial version
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
! 
! !INTERFACE:
subroutine noah32_setup()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use noah32_lsmMod
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Noah3.2 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in Noah3.2.
!  
! The routines invoked are: 
! \begin{description}
! \item[noah32\_setvegparms](\ref{noah32_setvegparms}) \newline
!   initializes the vegetation-related parameters in Noah3.2
! \item[noah32\_settbot](\ref{noah32_settbot}) \newline
!   initializes the bottom temperature fields
! \item[noah32\_setsoils](\ref{noah32_setsoils}) \newline
!   initializes the soil parameters
! \end{description}
!EOP

  implicit none

  call noah32_setvegparms(LIS_rc%lsm_index)
  call noah32_settbot(LIS_rc%lsm_index)
  call noah32_setsoils(LIS_rc%lsm_index)

end subroutine noah32_setup
