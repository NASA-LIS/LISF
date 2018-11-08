!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: pf_general
! 
! this file contains a collection of general particle filter
! subroutines and compact support subroutines
!
! !REVISION HISTORY: 
!
!   10 Aug 2017: Sujay Kumar; Initial Specification
!
!EOP
module pf_general
  
  implicit none
  
  private
  
  public :: pf_analysis
  
contains
  
!BOP
! 
! !ROUTINE: pf_analysis
! \label{pf_analysis}
! 
! !INTERFACE:  
  subroutine pf_analysis( gid, &
       N_state, N_obs, N_ens, &
       Observations, Obs_pred, Obs_err, Obs_cov, &
       State_incr, &
       State_lon, State_lat, xcompact, ycompact )
! !USES:    
    use pf_types
    use my_matrix_functions
    use LIS_logMod, only : LIS_logunit, LIS_endrun

    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: gid
    integer, intent(in) :: N_state, N_obs, N_ens    
    type(pf_obs_type), intent(in), dimension(N_obs) :: Observations     
    real, intent(in), dimension(N_obs,N_ens) :: Obs_pred
    real, intent(in), dimension(N_obs,N_ens) :: Obs_err
    real, intent(in), dimension(N_obs,N_obs) :: Obs_cov        
    real, intent(inout), dimension(N_state,N_ens) :: State_incr
    
    ! optional inputs
    real, dimension(N_state), intent(in), optional :: State_lon, State_lat     
    real, intent(in), optional :: xcompact       ! [deg] longitude
    real, intent(in), optional :: ycompact       ! [deg] latitude

!
! !DESCRIPTION:
!   
! perform Pf update
!
! IMPORTANT:
! on input, State\_incr must contain State\_minus(1:N\_state,1:N\_ens)
! on output, State\_incr contains the increments
!
! if optional inputs State\_lon, State\_lat, xcompact, and ycompact
! are present, Hadamard product is applied to HPHt and PHt
!EOP
    

    ! -----------------------------
    
    ! locals
    
    integer          :: n_e, i, ii, jj, kk

    real :: PHt_ij, dx, dy
    
    real, dimension(N_state,N_ens) :: State_prime
    real, dimension(N_state)       :: State_bar
    real, dimension(N_state)       :: State_incr_tmp
    
    real, dimension(N_obs,N_ens)   :: Obs_pred_prime
    real, dimension(N_obs)         :: Obs_pred_bar
    real, dimension(1,N_obs)       :: rhs
    real, dimension(1,N_obs)       :: rhs_temp
    
    real, dimension(N_ens)         :: weights    
    
    real, dimension(1,1)   :: Repr_matrix
    
    integer,          dimension(N_obs)         :: indx

    integer           :: IPIV(N_obs)   
    real*8            :: Obs_cov_dp(N_obs, N_obs)
    real*8            :: Obs_cov_dp_inv(N_obs, N_obs)
    real*8            :: Work(N_obs)
    integer           :: info

    do n_e=1,N_ens
       
       do i=1,N_obs
          if(Observations(i)%pert_type.eq.0) then 
             rhs(1,i) = Observations(i)%value + Obs_err(i,n_e) - Obs_pred(i,n_e)
          elseif(Observations(i)%pert_type.eq.1) then 
             rhs(1,i) = Observations(i)%value * Obs_err(i,n_e) - Obs_pred(i,n_e)
          endif
       end do

       Obs_cov_dp = Obs_cov
#if (defined LAPACK)
       call ZGETRF(N_obs,N_obs,Obs_cov_dp,N_obs,IPIV,info)       
       if(info.eq.0) then 
          call ZGETRI(N_obs,Obs_cov_dp,N_obs,IPIV,Work,N_obs,info)
       endif
       Obs_cov_dp_inv = Obs_cov_dp
#else
       call Minverse(Obs_cov_dp, Obs_cov_dp_inv,N_obs)
#endif
                  
       rhs_temp = matmul( rhs, Obs_cov_dp_inv)
       
       Repr_matrix = matmul(rhs_temp, transpose(rhs))
       weights(n_e) = exp(-0.5*Repr_matrix(1,1))
    enddo
   
  end subroutine pf_analysis

  subroutine Minverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
    implicit none 
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k
    
    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0
    
    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do
    
    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
       L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
       do i=1,j
          U(i,j) = a(i,j)
       end do
    end do
    
    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
       b(k)=1.0
       d(1) = b(1)
       ! Step 3a: Solve Ld=b using the forward substitution
       do i=2,n
          d(i)=b(i)
          do j=1,i-1
             d(i) = d(i) - L(i,j)*d(j)
          end do
       end do
       ! Step 3b: Solve Ux=d using the back substitution
       x(n)=d(n)/U(n,n)
       do i = n-1,1,-1
          x(i) = d(i)
          do j=n,i+1,-1
             x(i)=x(i)-U(i,j)*x(j)
          end do
          x(i) = x(i)/u(i,i)
       end do
       ! Step 3c: fill the solutions x(n) into column k of C
       do i=1,n
          c(i,k) = x(i)
       end do
       b(k)=0.0
    end do
  end subroutine Minverse

end module pf_general
