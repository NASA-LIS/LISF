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
! !ROUTINE: l_to_rc
! \label{l_to_rc}
!
! !REVISION HISTORY:
! ! We treat the 2D state variable as a matrix of 
! ! nrow x ncol, for example ds(nrow, ncol).
! ! The matrix is reorganized into a
! ! 1D array in a column by column manner.
! ! We can determine the column index (c_idx)
! ! and column index (r_idx) according to the index (l) 
! ! of the 1D array. 
!   5/9/16: Shugong Wang; initial implementation for LIS 7 and jules 4.3 
!   3/6/18: Shugong Wang; updated for JULES 5.0
subroutine l_to_rc(lev, nrow, r_idx, c_idx)
    implicit none 
    integer, intent(in)  :: lev     ! index of 1D array 
    integer, intent(in)  :: nrow    ! number of rows in the 2D state variable 
    integer, intent(out) :: r_idx   ! index of row
    integer, intent(out) :: c_idx   ! index of column
    
    c_idx=(lev-1)/nrow+1
    r_idx =lev - (int((lev-1)/nrow) * nrow)
endsubroutine l_to_rc
