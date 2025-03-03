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
! !ROUTINE: hyssib_dynsetup
! \label{hyssib_dynsetup}
!
! !REVISION HISTORY:
! 15 Apr 2002: Sujay Kumar, Initial Specification
! 23 Apr 2004: David Mocko, Conversion from NOAH to HY-SSiB
! 25 Aug 2007: Chuck Alonge, Updates for LIS 5.0 compliance
!
! !INTERFACE:
subroutine hyssib_dynsetup(n)
! !DESCRIPTION:
!  This routine sets up the time-dependent variables in Hyssib. 
!  Currently this routine is a placeholder
! 
!EOP     
  implicit none
  integer, intent(in) :: n

  !update vegetation data from monthly tables
  call hyssib_gfrac()

  return
end subroutine hyssib_dynsetup

