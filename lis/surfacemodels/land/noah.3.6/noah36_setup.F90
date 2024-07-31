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
! !ROUTINE: noah36_setup
! \label{noah36_setup}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial version
!   8 May 2009: Sujay Kumar; additions for Noah3.1
!  27 Oct 2010: David Mocko, changes for Noah3.1 in LIS6.1
!   7 Nov 2010: David Mocko, changes for Noah3.2 in LIS6.1
!   9 Sep 2011: David Mocko, changes for Noah3.3 in LIS6.1
!  30 Oct 2014: David Mocko, added Noah-3.6 into LIS-7
! 
! !INTERFACE:
subroutine noah36_setup()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use noah36_lsmMod
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Noah-3.6 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in Noah-3.6.
!  
! The routines invoked are: 
! \begin{description}
! \item[noah36\_setvegparms](\ref{noah36_setvegparms}) \newline
!   initializes the vegetation-related parameters in Noah-3.6
! \item[noah36\_settbot](\ref{noah36_settbot}) \newline
!   initializes the bottom temperature fields
! \item[noah36\_setsoils](\ref{noah36_setsoils}) \newline
!   initializes the soil parameters
! \end{description}
!EOP

  implicit none

  call noah36_setvegparms(LIS_rc%lsm_index)
  call noah36_settbot(LIS_rc%lsm_index)
  call noah36_setsoils(LIS_rc%lsm_index)

end subroutine noah36_setup
