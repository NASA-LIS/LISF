!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! 
!BOP
! 
! !MODULE: my_matrix_functions
! 
! this file contains a collection of matrix operation subroutines 
!
! !ROUTINE: 
! reichle,  1 May 01
! reichle, 18 Apr 06 - renamed subroutine sort()
!EOP
module my_matrix_functions
  
  implicit none
  
contains
  
!BOP
! 
! !ROUTINE: row_variance
! \label{row_variance_enkf}
! 
! !INTERFACE:
  subroutine row_variance( M, N, A, var )

    implicit none
! !ARGUMENTS:    
    integer, intent(in) :: M, N
    
    real, intent(in), dimension(M,N) :: A
    
    real, intent(out), dimension(M)  :: var    

! !DESCRIPTION:    
! compute variance of each row of an M-by-N matrix A
! 
!EOP
    
    
    ! locals
    
    integer :: i, j
    
    real :: x2, N_real, N_real_minus_one
    
    ! -------------------------------------------------------------
    
    N_real = real(N)
    
    N_real_minus_one = real(N-1)
    
    do i=1,M
       
       x2 = 0.0
       do j=1,N
          x2 = x2 + A(i,j)*A(i,j)
       end do
       
       var(i) = ( x2 - (sum(A(i,:))**2)/N_real )/N_real_minus_one

    end do
    
    ! deal with possible round-off errors
    ! reichle, 24 Sep 2004

    var = max(var,0.)
    
  end subroutine row_variance
  
  
  ! **********************************************************************
!BOP
! 
! !ROUTINE: row_third_moment
! \label{row_third_moment_enkf}
! 
! !INTERFACE:   
  subroutine row_third_moment( M, N, A, third_moment )

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: M, N
    real, intent(in), dimension(M,N) :: A
    real, intent(out), dimension(M)  :: third_moment    

! !DESCRIPTION:    
! compute third moment of each row of an M-by-N matrix A
!
! \begin{equation} 
! third\_moment = 1/N * sum_{i=1}^N (x_i - mean(x))**3 
! \end{equation}
! 
!EOP
    
    ! locals
    
    integer :: i, j
    
    real :: x3, mx, N_real
    
    ! -------------------------------------------------------------
    
    N_real = real(N)
    
    do i=1,M
       
       mx = sum(A(i,:))/N_real
       
       x3 = 0.0
       do j=1,N
          x3 = x3 + (A(i,j)-mx)**3
       end do
       
       third_moment(i) = x3/N_real
       
    end do
    
  end subroutine row_third_moment

  
  ! **********************************************************************
!BOP
! 
! !ROUTINE: five_number_summary
! \label{five_number_summary_enkf}
!  
! !INTERFACE: 
  subroutine five_number_summary( M, N, A, five_numbers )

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: M, N
    real, intent(in), dimension(M,N)  :: A
    real, intent(out), dimension(M,5)  :: five_numbers
    
!
! !DESCRIPTION:    
! get five number summary (median, lower and upper quartiles, min and max)
! of each row of an M-by-N data matrix A
!
! inputs:
!   M  : number of rows of data = number of different data types
!   N  : number of columns of data = number of ensemble members
!   A  : M-by-N matrix, left unchanged by this function
!
! outputs:
!   five\_numbers : M-by-5 matrix containing statistical summary:
!                  column 1: min
!                  column 2: lower quartile
!                  column 3: median
!                  column 4: upper quartile
!                  column 5: max
!
! Type:   f90
! Author: Rolf Reichle
! Date:   2 May 2001
!EOP
    
    integer i,d
    
    real, dimension(N) :: tmpvec
    
    do i=1,M
       
       ! put i-th row of data into tmpvec
       
       tmpvec = A(i,:)
       
       ! sort tmpvec in ascending order
       
       call nr_sort( N, tmpvec )
     
       ! get min, max, median and quartiles
       
       ! min and max
       
       five_numbers(i,1) = tmpvec(1)       ! min
       
       five_numbers(i,5) = tmpvec(N)   ! max
       
       ! median
       
       if (mod(N,2) == 0) then
          five_numbers(i,3) = .5*(tmpvec(N/2)+tmpvec(N/2+1))
       else
          five_numbers(i,3) = tmpvec(N/2+1)
       end if
       
       ! quartiles 
       ! (follows Robert Johnson, "Elementary Statistics", PWS-Kent, p69, 1988)
       
       if (mod(N,4) == 0) then
          
          d = N/4
          
          five_numbers(i,2) = .5*(tmpvec(d)+tmpvec(d+1))             ! lower
          
          five_numbers(i,4) = .5*(tmpvec(N-d)+tmpvec(N-d+1)) ! upper
          
       else
          
          d = N/4+1         
          
          five_numbers(i,2) = tmpvec(d)           ! lower
          
          five_numbers(i,4) = tmpvec(N-d)     ! upper
          
       end if
       
    end do
    
  end subroutine five_number_summary

  ! ------------------------------------------------------------------
!BOP
! 
! !ROUTINE: adjust_mean
! \label{gmaoenkf_adjust_mean}
! 
! !INTERFACE:   
  subroutine adjust_mean( N_row, N_col, A, M )

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: N_row, N_col
    real, intent(inout), dimension(N_row,N_col)  :: A
    real, intent(in), optional, dimension(N_row) :: M
! 
! !DESCRIPTION:    
! adjust N\_row by N\_col matrix A such that 
! mean over columns for each row is given by the
! corresponding element in vector M of length N\_row
! 
! vector of mean values M is optional input, if not present 
! zero mean is assumed
!EOP    
    
    ! ----------------------------
    
    ! locals
    
    integer i
    
    real, dimension(N_row) :: correction
    
    ! ------------------------------------------------------------
    
    if (present(M)) then
       correction = M - sum(A,2)/real(N_col) 
    else
       correction = - sum(A,2)/real(N_col) 
    end if
    
    do i=1,N_col
       A(:,i) = A(:,i) + correction
    end do
    
  end subroutine adjust_mean
  
  ! ------------------------------------------------------------------
!BOP
! !ROUTINE: adjust_std
!  \label{adjust_std_enkf}
! 
! !INTERFACE:   
  subroutine adjust_std( N_row, N_col, A, std )

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: N_row, N_col
    real, intent(inout), dimension(N_row,N_col)  :: A
    real, intent(in), optional :: std

!
! !DESCRIPTION: 
! adjust N\_row by N\_col matrix A such that (sample) standard deviation
!  of all elements is exactly equal to std
! 
! std is optional input, if not present std=1 is assumed
!
!EOP    
    
    ! ----------------------------
    
    ! locals
    
    integer :: i, j
    
    real :: correction, sample_std
    
    ! ------------------------------------------------------------
    
    ! compute sample std

    call matrix_std( N_row, N_col, A, sample_std )
    
    if (present(std)) then
       correction = std/sample_std
    else
       correction = 1./sample_std
    end if
    
    do i=1,N_row
       do j=1,N_col
          A(i,j) = correction*A(i,j)
       end do
    end do
    
  end subroutine adjust_std
  
  ! ------------------------------------------------------------------
!BOP
! 
! !ROUTINE: matrix_std
! \label{matrix_std_enkf}
!   
! !INTERFACE:
  subroutine matrix_std( N_row, N_col, A, std )

    implicit none
! !ARGUMENTS:    
    integer, intent(in) :: N_row, N_col
    real, intent(inout), dimension(N_row,N_col)  :: A
    real, intent(out) :: std
! 
! !DESCRIPTION:
! compute std of all elements of N\_row by N\_col matrix A
! 
!EOP
    
    ! ----------------------------
    
    ! locals
    
    integer :: i, j
    
    real :: x2, m, N_real, N_real_minus_one
    
    ! ------------------------------------------------------------
    
    N_real = real(N_row)*real(N_col)
    
    N_real_minus_one = N_real - 1.
    
    ! compute sample std
    
    x2 = 0.0
    m  = 0.0
    
    do i=1,N_row
       do j=1,N_col
          m  = m  + A(i,j)
          x2 = x2 + A(i,j)*A(i,j)
       end do
    end do
    
    std = sqrt( ( x2 - m**2/N_real )/N_real_minus_one )
        
  end subroutine matrix_std
    
end module my_matrix_functions


! ***** EOF **************************************************************
