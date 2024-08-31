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
! !ROUTINE: noah271_setup
! \label{noah271_setup}
!
! !REVISION HISTORY:
!  28 Apr 2002: Sujay Kumar; Initial version
! 
! !INTERFACE:
subroutine noah271_setup()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use noah271_lsmMod
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Noah2.7.1 LSM.  These include the soils, greenness,
!  albedo, bottom temperature and the initialization of state
!  variables in Noah2.7.1.
!  
! The routines invoked are: 
! \begin{description}
! \item[noah271\_setvegparms](\ref{noah271_setvegparms}) \newline
!   initializes the vegetation-related parameters in Noah2.7.1
! \item[noah271\_settbot](\ref{noah271_settbot}) \newline
!   initializes the bottom temperature fields
! \item[noah271\_setsoils](\ref{noah271_setsoils}) \newline
!   initializes the soil parameters
! \end{description}
!EOP

  implicit none

! to generate the restart file from AGRMET output. 
! call noah271_read2drst()
! open(40,file='noah271_dom2_15km.rst',status='unknown',form='unformatted')
! call noah271_dump_restart(40)
! close(40)
! stop

  call noah271_setvegparms(LIS_rc%lsm_index)
  call noah271_settbot(LIS_rc%lsm_index)
  call noah271_setsoils(LIS_rc%lsm_index)

end subroutine noah271_setup
