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
! !ROUTINE: clsmf25_finalize
! \label{clsmf25_finalize}
!
! !DESCRIPTION:
!  
! Complete the finalize routines for catchment
!
! !REVISION HISTORY:
!
! 11 Feb 2006: Sujay Kumar; Initial Specification
! 
! !INTERFACE:
subroutine clsmf25_finalize()
! !USES:
  use LIS_coreMod, only : LIS_rc
  use clsmf25_lsmMod
!EOP
  implicit none
  integer :: n 
  
  do n=1,LIS_rc%nnest
     deallocate(clsmf25_struc(n)%cat_param)
     deallocate(clsmf25_struc(n)%met_force)
     deallocate(clsmf25_struc(n)%cat_progn)
     deallocate(clsmf25_struc(n)%cat_diagn)
     deallocate(clsmf25_struc(n)%cat_output)
     deallocate(clsmf25_struc(n)%modis_alb_param)
     deallocate(clsmf25_struc(n)%good_forcing_mask)
  enddo
  
  deallocate(clsmf25_struc)

end subroutine clsmf25_finalize
