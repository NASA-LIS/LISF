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
!  !MODULE: random_fields.F90
! 
!  !DESCRIPTION: 
!  This module contains a random field generator in 2d to 
!  generate a pair of random fields in 2d with zero mean
!
! subroutines rfg2d\_fft() and sqrt\_gauss\_spectrum are translated from 
!  C++ code rfg2d.C written for MIT EnKF work by reichle
!  (see janus:~reichle/nasa/EnKF)
!
! covariance is specified through its spectrum, so far only Gaussian
!
! IMPORTANT: read comments for function rfg2d\_fft()
!
! written for NSIPP - EnKF
! Type:   f90
! Author: Rolf Reichle
! Date:   2 Nov 2001

!   
!  !REVISION HISTORY: 
!  18Feb05 Rolf Reichle  Initial Specification
!                        updated for use with module nr\_ran2\_gasdev
!                        deleted use of module select\_kinds
!  07Jul05 Sujay Kumar   Specification in LIS
! 
!EOP
module random_fields
    
  use LIS_ran2_gasdev
  
  implicit none

  private
  
  public :: rfg2d_fft
  public :: generate_white_field
  public :: get_fft_grid

  real, parameter :: MY_PI = 3.14159265
  
contains
  
 
!BOP
! 
! !ROUTINE: sqrt_gauss_spectrum_2d
! \label{sqrt_gauss_spectrum_2d}
!
! !INTERFACE:   
  subroutine sqrt_gauss_spectrum_2d( &
       variance, N_x, N_y, dkx, dky, lambda_x, lambda_y, &
       sqrt_spectrum_2d )
     
    implicit none
! !ARGUMENTS:     
    real, intent(in) :: variance, dkx, dky, lambda_x, lambda_y
    integer, intent(in) :: N_x, N_y
    real, intent(out), dimension(N_x,N_y) :: sqrt_spectrum_2d
! 
! !DESCRIPTION:
!
! get SQUARE ROOT of 2d Gaussian spectrum (incl volume element)
!
! 2d Gaussian spectrum:
!   
! \begin{equation*}
! S(kx,ky) = variance
!            *
!            lambda_x*lambda_y/(2*pi) 
!            * 
!            exp( -(kx^2*lambda_x^2 + ky^2*lambda_y^2)/2 )
! \end{equation*}
! 
! return: sqrt( S*dkx*dky )
!
! that is return the SQUARE ROOT of the spectrum multiplied with the
!  square root of the volume element d2k=dkx*dky of the ifft integral
!
! spectrum is returned in "wrap-around" order compatible with CXML and 
!  matlab FFT algorithms
! 
! inputs:
!  variance : variance desired for complex field, if pair of real fields 
!             is used each field must eventually be multiplied with sqrt(2)
!  N\_x      : number of nodes in x direction
!  N\_y      : number of nodes in y direction
!  dkx      : wave number spacing in x direction
!  dky      : wave number spacing in y direction
!  lambda\_x : decorrelation length in x direction 
!  lambda\_y : decorrelation length in y direction
! 
!EOP    
    ! -------------------------------------
    
    ! local variables

    integer :: i, j

    real :: fac, lamx2dkx2, lamy2dky2 

    real, dimension(N_x) :: lx2kx2
    real, dimension(N_y) :: ly2ky2
    
    ! ------------------------------------------------------------------

    ! factor includes sqrt of volume element of ifft integral
    
    fac       = sqrt( variance*lambda_x*lambda_y/(2.*MY_PI)*dkx*dky )
    
    lamx2dkx2 = lambda_x*lambda_x*dkx*dkx
    lamy2dky2 = lambda_y*lambda_y*dky*dky
    
  
    ! precompute (lambda_x*k_x)^2 in "wrap-around" order suitable for CXML fft
    
    do i=1,(N_x/2)
       lx2kx2(i) = lamx2dkx2*(i-1)*(i-1)
    end do
    
    do i=(N_x/2+1),N_x
       lx2kx2(i) = lamx2dkx2*(N_x-i+1)*(N_x-i+1) ! minus drops out when squared
    end do
    
    ! precompute (lambda_y*k_y)^2 in "wrap-around" order suitable for CXML fft
    
    do j=1,(N_y/2)
       ly2ky2(j) = lamy2dky2*(j-1)*(j-1)
    end do
    
    do j=(N_y/2+1),N_y
       ly2ky2(j) = lamy2dky2*(N_y-j+1)*(N_y-j+1) ! minus drops out when squared
    end do
    
    ! assemble spectrum in "wrap-around" order
    
    do i=1,N_x
       do j=1,N_y
          sqrt_spectrum_2d(i,j) = fac * exp(-.25*(lx2kx2(i) + ly2ky2(j)) )  
       end do
    end do
    
    return
    
  end subroutine sqrt_gauss_spectrum_2d

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  ! ----------------------------------------------------------------------
  !
  ! subroutine rfg2d_fft()
  !  
  ! generate a pair of 2d zero-mean random fields using FFT method
  ! (so far only Gaussian covariance implemented)
  !
  ! NOTE: implemented with index counters of type int, must have
  !
  !          N_x*N_y < maximum integer on given machine
  !
  !       on yama/alborz (at MIT) int varies from -2147483648...2147483647
  !       -> can handle up to N_x*N_y = (32768)^2
  !       if larger fields are needed, rewrite with type long int etc
  !
  ! NOTE: The fft method for generating random fields produces 
  !       fields that are periodic with the size of the domain,
  !       that is the field is (almost) the same on each boundary.
  !       This introduces unwanted correlations at lags shorter than
  !       the domain size. Therefore, only a part of the generated
  !       field is usable. As a rule of thumb, the fields should be
  !       generated on a grid that is two correlation lengths bigger
  !       than the field on which the grid is desired. Then cut out
  !       fields of the necessary size.
  !       This procedure is included in rfg2d_fft().
  !
  ! NOTE: The variance specified as input is the theoretical variance of 
  !       the complex field that is obtained from the inverse fft of the 
  !       realization dZ.
  !       The sample variance of this *complex* field is a FIXED (non-random) 
  !       number which depends on the size of the domain, the grid spacing,
  !       and the correlation length (but not on the random seed!!).
  !       (This number is non-random because the variance is the integral
  !        of the absolute value of the spectrum. In this integral the 
  !        randomness disappears because we only choose a random phase angle.)
  !       For vanishing discretization and spectral truncation error,
  !       this number converges to the theoretical (input) value.
  !
  !       This function is set up to use the pair of two real fields, where
  !
  !                field1 = sqrt(2)*real(ifft(dZ))
  !                field2 = sqrt(2)*imag(ifft(dZ)).
  !     
  !       The factor sqrt(2) re-scales the variance of field1 and field2
  !       such that for vanishing discretization error and spectral
  !       truncation error each field converges to the specified theoretical
  !       (input) variance. 
  !       NOTE: The sum of the sample variances of the two real fields 
  !       is equal to the (FIXED) sample variance of the complex field
  !       (before re-scaling with sqrt(2)). 
  !       The individual sample variances within each pair vary from 
  !       realization to realization.
  !
  ! -------------------------------------------------------------------
  
  subroutine rfg2d_fft( &
       rseed, variance, N_x, N_y, N_x_fft, N_y_fft, dx, dy, lambda_x, lambda_y, &
       field1, field2 )

#ifndef CXML_FFT_AVAILABLE
    use nr_fft
#endif    

    implicit none
    
    ! use Compaq Extended Math Library (CXML)
    ! (free with Compaq Fortran Compiler)
    
#ifdef CXML_FFT_AVAILABLE
    include 'CXMLDEF.FOR'
#endif
    integer, intent(in) :: N_x_fft, N_y_fft
    integer, dimension(NRANDSEED, N_x_fft, N_y_fft), intent(inout) :: rseed
    
    real, intent(in) :: variance, dx, dy, lambda_x, lambda_y
    
    integer, intent(in) :: N_x, N_y
    
    real, intent(out), dimension(N_x,N_y) :: field1, field2
    
    ! -------------------------------------
    
    ! local variables
    
    integer :: i, j, m!, N_x_fft, N_y_fft
    
    real :: dkx, dky, theta, ran_num
    
    real, dimension(:,:), allocatable :: field1_fft, field2_fft

#ifdef CXML_FFT_AVAILABLE

    ! variables for CXML library routines
    
    integer :: cxml_status
    
#else

    ! using nr_fft
    
    real, dimension(:), allocatable :: tmpdata
    
    integer :: N_xy_fft, k
    
    integer, dimension(2) :: nn

#endif

    ! ----------------------------------------------------------
    !
    ! get fft grid size of suitable power of two and at least 
    ! two correlation lengths larger than N_x*N_y
!SVK DEBUG   
!    call get_fft_grid( N_x, N_y, dx, dy, lambda_x, lambda_y, N_x_fft, N_y_fft )
    
    allocate(field1_fft(N_x_fft,N_y_fft))
    allocate(field2_fft(N_x_fft,N_y_fft))
    
    ! ----------------------------------
    
    dkx = (2.*MY_PI)/(float(N_x_fft)*dx)
    dky = (2.*MY_PI)/(float(N_y_fft)*dy)
    
    ! follow Ruan & McLaughlin, 1998:
    !
    ! compute dZ = H * exp(i*theta) * sqrt(d2k)     
    
    ! start with square root of spectrum (factor H*sqrt(d2k)), put into field1
    
    call sqrt_gauss_spectrum_2d( &
         variance, N_x_fft, N_y_fft, dkx, dky, lambda_x, lambda_y, &
         field1_fft)
    
    ! multiply with random phase angle 

    do i=1,N_x_fft
       do j=1,N_y_fft
          
          ! theta = (2.*MY_PI)*((double) rand()/RAND_MAX)
          
          ! theta = (2.*MY_PI)*ran2(RSEED2)  ! random phase angle
          
          call nr_ran2(rseed(:,i,j), ran_num)          
          
          theta = (2.*MY_PI)*ran_num          ! random phase angle
          
          field2_fft(i,j) = sin( theta )*field1_fft(i,j)
          field1_fft(i,j) = cos( theta )*field1_fft(i,j)
          
       end do
    end do
    
    
    ! force dZ(1,1) to zero (-> zero mean random field)
    
    field1_fft(1,1) = 0.
    field2_fft(1,1) = 0.
  
    ! -----------------------------------------------------------
    !
    ! apply 2d fft
    
#ifdef CXML_FFT_AVAILABLE
    
    ! invoke appropriate CXML 2d fft - SINGLE precision!
    !
    ! order of arguments is wrong in web-documentation!!!
    ! (compare 1d and 3d ffts)
    !
    ! Complex transform in real data format: 
    ! status = {C,Z}FFT_2D (input_format, output_format, direction, 
    !                       in_real, in_imag, out_real, out_imag, ni, nj,
    !                       lda, ni_stride, nj_stride)

    cxml_status = cfft_2d('R', 'R', 'B',                                    & 
         field1_fft, field2_fft, field1_fft, field2_fft, N_x_fft, N_y_fft,  &
         N_x_fft, 1, 1)
    
    if (cxml_status .ne. 0) then
       write (*,*) 'error when using CXML cfft_2d in rfg2d_fft()'
       stop
    end if
    
#else
    
    ! use nr_fft

    N_xy_fft = N_x_fft*N_y_fft
    
    allocate(tmpdata(2*N_xy_fft))
    
    nn(1) = N_x_fft
    nn(2) = N_y_fft
    
    ! fill tmpdata according to Figures 12.2.2 and 12.4.1 of f77 NR book
    
    k=0
    do j=1,N_y_fft
       do i=1,N_x_fft
          k=k+1
          tmpdata(k) = field1_fft(i,j)
          k=k+1
          tmpdata(k) = field2_fft(i,j)
       end do
    end do
    
    ! apply nr_fft
    
    call fourn(2*N_xy_fft,tmpdata,nn,2,1)
    
    ! extract random fields from tmpdata
        
    k=0
    do j=1,N_y_fft
       do i=1,N_x_fft      
          k=k+1
          field1_fft(i,j) = tmpdata(k)/real(N_xy_fft)
          k=k+1
          field2_fft(i,j) = tmpdata(k)/real(N_xy_fft)
       end do
    end do

    deallocate(tmpdata)
    
#endif    
    
    ! multiply with factor sqrt(2) to get correct variance
    ! (see above and p. 388 Ruan and McLaughlin, 1998),
    ! also multiply with N_x_fft*N_y_fft to get correct scaling,
    ! also retain only useable part of field?_fft
    
    field1 = sqrt(2.)*float(N_x_fft)*float(N_y_fft)*field1_fft(1:N_x,1:N_y)
    field2 = sqrt(2.)*float(N_x_fft)*float(N_y_fft)*field2_fft(1:N_x,1:N_y)
    
    deallocate(field1_fft)
    deallocate(field2_fft)
    
  end subroutine rfg2d_fft

  ! ***********************************************************
  
  subroutine get_fft_grid( N_x, N_y, dx, dy, lambda_x, lambda_y, &
       N_x_fft, N_y_fft )
    
    ! get a grid that can be used for 2d fft random number generator
    ! (this fft grid must extend beyond the desired random field
    !  by about two correlations lengths, it must also have N_x_fft and 
    !  N_y_fft that are powers of two)
    
    implicit none
    
    ! ----------------------------------------------------------------
    
    integer, intent(in) :: N_x, N_y
    
    real, intent(in) :: dx, dy, lambda_x, lambda_y
    
    integer, intent(out) :: N_x_fft, N_y_fft
    
    ! ------------------------------
    
    ! local variables
    
    ! specify by how many correlation lengths the fft grid must be 
    ! be larger than the grid2cat grid
    
    real, parameter :: mult_of_xcorr = 2.
    real, parameter :: mult_of_ycorr = 2.
    
    ! ----------------------------------------------------------------
    
    ! add minimum required correlation lengths 

    N_x_fft = N_x + ceiling(mult_of_xcorr*lambda_x/dx)
    N_y_fft = N_y + ceiling(mult_of_ycorr*lambda_y/dy)
    
    ! make sure N_x_fft, N_y_fft are power of two
    
    N_x_fft = 2**ceiling( log(real(N_x_fft))/log(2.) )
    N_y_fft = 2**ceiling( log(real(N_y_fft))/log(2.) )
    
    ! echo findings

#if 0
    write (*,*)
    write (*,*) 'desired random field:'
    write (*,*) 'N_x     = ', N_x,      ' N_y      = ', N_y
    write (*,*) 'dx      = ', dx,       ' dy       = ', dy
    write (*,*) 'xcorr   = ', lambda_x, ' lambda_y = ', lambda_y
    write (*,*)
    write (*,*) 'grid used for fft: '
    write (*,*) 'N_x_fft = ', N_x_fft,  ' N_y_fft = ', N_y_fft
    write (*,*)
#endif    

  end subroutine get_fft_grid

  ! *******************************************************************

  subroutine generate_white_field(N_x, N_y, rseed, rfield )
    
    ! generate standard-normal random field that is white in space
        
    ! note that nr_gasdev always produces a pair of random numbers
    !
    ! do not store random numbers between subsequent calls to 
    ! the random field generator - works best if fields are large
    ! (ie. avoid using this subroutine with N_x=N_y=1)
    
    implicit none
    
    integer, intent(in) :: N_x, N_y
    
    integer, dimension(NRANDSEED,N_x,N_y), intent(inout) :: rseed
    
    real, dimension(N_x,N_y), intent(out) :: rfield
    
    ! local variables
    
    logical :: gauss_stored
    
    integer :: i, j
    
    real, dimension(2) :: tmp_real
    
    ! -----------------------------------------------------
    
    gauss_stored = .false.

    do i=1,N_x
       do j=1,N_y
          
          if (.not. gauss_stored) then
                
             call nr_gasdev(rseed(:,i,j), tmp_real)
             
             rfield(i,j) = tmp_real(1)
             
             gauss_stored = .true.
             
          else
             
             rfield(i,j) = tmp_real(2)
             
             gauss_stored = .false.
             
          end if
          
       end do
    end do
    
  end subroutine generate_white_field
  
  ! *******************************************************************

end module random_fields

! ==============================================================

! driver routine for testing

#if 0

program test_rfg2d
  
!  use random_fields
!  use nr_ran2_gasdev
  
  implicit none
  
  integer :: N_x, N_y, i, j, n_e

  real :: dx, dy, lambda_x, lambda_y, variance
  
  real, allocatable, dimension(:,:) :: field1, field2
  
  character(300) :: file_name
  character(10)  :: n_e_string
  character(100) :: output_format
  character(10)  :: tmp_string

  integer :: RSEEDCONST
  
  integer, dimension(NRANDSEED) :: rseed
  
  character(5) :: fft_tag
    
  ! --------------------------------
  
  RSEEDCONST = -777
  
  rseed(1) = RSEEDCONST
  
  write (*,*) RSEEDCONST
  
  call init_randseed(rseed)
  
  ! ------------------------------------------

  dx = 5000.
  dy = 5000.
  
  N_x      = 144
  N_y      =  91
  
  lambda_x = 45000.
  lambda_y = 45000.

  variance = 1.
  
  allocate(field1(N_x,N_y))
  allocate(field2(N_x,N_y))
  
#ifdef CXML_FFT_AVAILABLE
  fft_tag = 'cxml.'
#else
  fft_tag = 'nrxx.'
#endif

  ! -----------------------------------------------------------
  
  ! get N_e fields

  N_e = 10
  
  do n_e=1,N_e,2
     
     call rfg2d_fft( &
          rseed, &
          variance, &
          N_x, &
          N_y, &
          dx, &
          dy, &
          lambda_x, &
          lambda_y, &
          field1, &
          field2 )
      
     ! write to file
     
     ! field1
     
     write(n_e_string,  '(i3.3)') n_e
     
     file_name = 'rf.'//fft_tag// n_e_string(1:len_trim(n_e_string)) // '.dat'
     
     write(tmp_string, '(i3.3)') N_y
     
     output_format = '(' // tmp_string(1:len_trim(tmp_string)) // '(1x,e13.5))'
     
     open (10,file=file_name,status='unknown')
          
     do i=1,N_x

        write (10,output_format(1:len_trim(output_format)))  &
             (field1(i,j), j=1,N_y)
        
     end do
     close (10,status='keep')


     ! field2
     
     write(n_e_string,  '(i3.3)') (n_e+1)
     
     file_name = 'rf.' //fft_tag// n_e_string(1:len_trim(n_e_string)) // '.dat'
     
     write(tmp_string, '(i3.3)') N_y
     
     output_format = '(' // tmp_string(1:len_trim(tmp_string)) // '(1x,e13.5))'

     open (10,file=file_name,status='unknown')
          
     do i=1,N_x

        write (10,output_format(1:len_trim(output_format)))  &
             (field2(i,j), j=1,N_y)
        
     end do
     close (10,status='keep')
     
  end do

end program test_rfg2d


#endif

