!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: noah39_setup
! \label{noah39_setup}
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
subroutine noah39_setup()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use noah39_lsmMod
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Noah-3.9 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in Noah-3.9.
!  
! The routines invoked are: 
! \begin{description}
! \item[noah39\_setvegparms](\ref{noah39_setvegparms}) \newline
!   initializes the vegetation-related parameters in Noah-3.9
! \item[noah39\_settbot](\ref{noah39_settbot}) \newline
!   initializes the bottom temperature fields
! \item[noah39\_setsoils](\ref{noah39_setsoils}) \newline
!   initializes the soil parameters
! \end{description}
!EOP

  implicit none

  call noah39_setvegparms(LIS_rc%lsm_index)
  call noah39_settbot(LIS_rc%lsm_index)
  call noah39_setsoils(LIS_rc%lsm_index)

end subroutine noah39_setup
