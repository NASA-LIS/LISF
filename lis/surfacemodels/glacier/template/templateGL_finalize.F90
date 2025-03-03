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
! !ROUTINE: templateGL_finalize
! \label{templateGL_finalize}
!
! !REVISION HISTORY:
!
!   06 Apr 2018: Sujay Kumar, Initial imlementation
!
!
! !INTERFACE:
subroutine templateGL_finalize(n)
! !USES:

!
! !DESCRIPTION:
!
!  This routine cleans up the allocated memory structures in templateGL
!
!EOP
    implicit none   
    
    integer, intent(in) :: n

  end subroutine templateGL_finalize

