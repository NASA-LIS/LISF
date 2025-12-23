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
! !ROUTINE: AC72_finalize
! \label{AC72_finalize}
!
! !REVISION HISTORY:
!  04 NOV 2024, Louise Busschaert; initial implementation for AC72
!
! !INTERFACE:
subroutine AC72_finalize(n)

  ! !USES:
  use AC72_lsmMod
  use LIS_coreMod, only : LIS_rc

  !
  ! !DESCRIPTION:
  !
  !  This routine cleans up the allocated memory structures in AC72
  !
  !EOP

  implicit none

  integer, intent(inout) :: n

  do n=1, LIS_rc%nnest
     !deallocate all vars
     ! free memory for ac72, the data at tile level
     deallocate(AC72_struc(n)%ac72)

     ! free memory for constant parameter
     deallocate(AC72_struc(n)%Thickness)

     ! free memory for initial state variable
     deallocate(AC72_struc(n)%init_smc)
  end do ! nest loop

  deallocate(AC72_struc)

end subroutine AC72_finalize

