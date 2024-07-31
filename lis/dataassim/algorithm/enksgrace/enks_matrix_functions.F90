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
! !MODULE: enks_matrix_functions
! 
! this file contains a collection of matrix operation subroutines 
!
! !ROUTINE: 
! reichle,  1 May 01
! reichle, 18 Apr 06 - renamed subroutine sort()
!EOP
module enks_matrix_functions
  
  implicit none
  
contains
  
!BOP
! 
! !ROUTINE: row_variance
! \label{row_variance_enksgrace}
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
! \label{row_third_moment_enksgrace}
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
! \label{five_number_summary_enksgrace}
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
! \label{enksgrace_adjust_mean}
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
!  \label{adjust_std_enksgrace}
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
! \label{matrix_std_enksgrace}
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
    
!BOP
! !ROUTINE: nr_sort
!  \label{nr_sort_enksgrace}
!  
!  QUICKSORT algorithm after Numerical Recipes
! 
!  Sorts an array arr[1..n] into ascending numerical order using the
!  Quicksort algorithm. n is input; arr is replaced on output by its
!  sorted rearrangement.
! 
! !INTERFACE:       
SUBROUTINE nr_sort(n,arr)
!EOP
  implicit none

  INTEGER n,M,NSTACK
  REAL arr(n)
  PARAMETER (M=7,NSTACK=50)
  INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
  REAL a,temp
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M)then
     do 12 j=l+1,ir
        a=arr(j)
        do 11 i=j-1,l,-1
           if(arr(i).le.a)goto 2
           arr(i+1)=arr(i)
11         continue
           i=l-1
2          arr(i+1)=a
12         continue
           if(jstack.eq.0)return
           ir=istack(jstack)
           l=istack(jstack-1)
           jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)print*,'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
   end subroutine

end module enks_matrix_functions


! ***** EOF **************************************************************
