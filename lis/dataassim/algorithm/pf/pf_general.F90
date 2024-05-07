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
! !MODULE: pf_general
! 
! this file contains a collection of general Particle Filter
! subroutines and compact support subroutines
!
! !REVISION HISTORY: 
! Abolafia-Rosenzweig, August 2019 - Implement Particle filter using sequential importance resampling and output state increment for DA

!EOP
module pf_general
  
  use LIS_logMod

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
       State_lon, State_lat, xcompact, ycompact)

! !USES:    
    use pf_types
    use my_matrix_functions
    use LIS_logMod, only : LIS_logunit, LIS_endrun
    use LIS_numerRecipesMod, only : LIS_rand_func
    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: gid
    integer, intent(in) :: N_state, N_obs, N_ens    
    type(obs_type), intent(in), dimension(N_obs) :: Observations     
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
    
    integer          :: n_e, i, ii, jj, kk, ij, jk 

    real :: PHt_ij, dx, dy
    
    real, dimension(N_state,N_ens) :: State_prime
    real, dimension(N_state)       :: State_bar
    real, dimension(N_state)       :: State_incr_tmp
    
    real, dimension(N_obs,N_ens)   :: Obs_pred_prime
    real, dimension(N_obs)         :: Obs_pred_bar
    real, dimension(N_obs)         :: rhs
    
    real, dimension(N_ens)         :: weights    
    real, dimension(N_obs,N_obs)   :: Repr_matrix
    
    integer,          dimension(N_obs)         :: indx
    
    real                                 :: dweight
    logical                              :: apply_hadamard
    real                                 :: det_Obs_cov
    real, dimension(N_obs,N_ens)         :: innovation
    real, dimension(N_obs,N_ens)         :: Pw_raw
    real, dimension(1,N_ens)            :: Pw_norm
    real, dimension(N_ens)                 :: Pw_norm_vec
    real, dimension(1,N_ens)               :: Pw_combined
    real, dimension(N_ens)               :: Pw_cumsum
    real, dimension(N_ens)               :: new_particle_idx
    !real, dimension(N_obs,N_obs)         :: Obs_cov_inv
    real                                 :: ran_face
    real, dimension(N_state,N_ens)       :: updated_state
    real                                 :: current_p
    real                                 :: IDX
    real                                 :: iii
    real                                 :: Pw_sum
    
    ! ------------------------------------------------------------------
    ! find out whether Hadamard product should be applied

    apply_hadamard = (               &
         present(State_lon)    .and. &
         present(State_lat)    .and. &
         present(xcompact)     .and. &
         present(ycompact)              )
    
    ! ----------------------
    ! IMPORTANT: on input, State_incr contains State_minus(1:N_state,1:N_ens)
    ! ----------------------

    !calculate determinate of observation error covariance matrix
    det_Obs_cov=determinant(Obs_cov)
    
    !calcualte innovations, or difference between obs prediction made by the model and observations
    do i=1,N_obs
       do jj=1,N_ens
          innovation(i,jj)= Observations(i)%value - Obs_pred(i,jj)
       enddo
    enddo   
    
    !calculate the likelihood function for each enesmble member
    do i=1,N_obs
       do kk=1,N_ens
          Pw_raw(i,kk)= 1/( (2*3.1415927)**(N_obs/2) *&
               det_Obs_cov**(1/2)) * &
               exp(- ( (innovation(i,kk))**2 / (2* Obs_cov(i,i) ) ) )
       end do
    end do

    !normalize PF weights assigned from each observation between 0-1 
    do i=1,N_obs
       Pw_sum=0
       do kk=1,N_ens
          Pw_sum=Pw_sum+Pw_raw(i,kk)
       end do
       if(Pw_sum.ne.0) then 
          Pw_raw(i,:)=Pw_raw(i,:)/Pw_sum
       else
          Pw_raw(i,:) = 0.0
       endif
    end do
    
    !account for all observations in this time step, similar to 
    !smoothing technique by taking the product of the weights for
    !each ensemble member derived from all assimilated observations 
    do ij=1,N_obs
       if (ij.eq.1) then
          Pw_combined(1,:)=Pw_raw(ij,:)
       else if (ij.gt.1) then
          Pw_combined(1,:)=Pw_combined(1,:)*Pw_raw(ij,:)
       endif
    end do
    

    !normalize into probability "weights" that account for all observations
    if(sum(Pw_combined(1,:)).ne.0) then 
       Pw_norm(1,:)=Pw_combined(1,:)/sum(Pw_combined(1,:))

       !vectorize the normalized weights
       do n_e=1,N_ens
          Pw_norm_vec(n_e)=Pw_norm(1,n_e)
       end do
       
       !resample using Sequential importance resampling (SIR) described 
       !in Weerts and Serafy, 2006 - Section 2.3.1
       
       !calculate the cumulative sum of Pw_combined for random dice generation:
       do n_e=1,N_ens
          Pw_cumsum(n_e)=sum(Pw_norm_vec(1:n_e))
       end do
       

       do n_e=1,N_ens
          !generate random number between 0-1 (random dice roll)
          call LIS_rand_func(1,ran_face)
          !select first insance of the random number 
          !being exceeded by the posterior cumulative
          !sumation of the posterior PDF
          
          !vector of 0's and 1's. 1's begin in the first 
          !instance of ran_face exceeding the Pw_cumsum vector
          do jk=1,N_ens
             if (Pw_cumsum(jk).le.ran_face) then
                new_particle_idx(jk)=0
             else if (Pw_cumsum(jk).gt.ran_face) then
                new_particle_idx(jk)=1
             end if
          end do
          
          !select the first instance of the logical statement 
          !occuring (or the first non-zero instance) and resample
          !from this instance
          iii=0
          do i=1,N_ens
             current_p=new_particle_idx(i)
             if (current_p.gt.0 .AND. iii.eq.0) then
                do jj=1,N_state
                   updated_state(jj,n_e)=State_incr(jj,i)
                end do
                iii=iii+1
             end if
          end do
          
       end do
    
    !the increment is the analysis minus the model's first guess
       State_incr = (updated_state-State_incr)    
    else
       State_incr = 0 
    endif
    
  end subroutine pf_analysis

!BOP
! 
! !ROUTINE: hadamard_for_repr_matrix
! \label{hadamard_for_repr_matrix_pf}
! 
! !INTERFACE:  
  subroutine hadamard_for_repr_matrix( N_obs, Observations,             & 
       xcompact, ycompact,                                              &
       Repr_Matrix )
! !USES:    
    use pf_types                     ! for obs_type
    
    implicit none
! !ARGUMENTS:     
    integer, intent(in) :: N_obs
    
    type(obs_type), intent(in), dimension(N_obs) :: Observations       
    
    real, intent(in) :: xcompact       ! [deg] longitude
    real, intent(in) :: ycompact       ! [deg] latitude
    
    real, dimension(N_obs,N_obs), intent(inout) :: Repr_matrix
!EOP    
    ! locals
    
    integer :: i, j
    
    real :: tmpreal, dx, dy
    
    ! ----------------------------------
    
    do i=1,N_obs
       do j=i+1,N_obs
          
          dx = Observations(i)%lon - Observations(j)%lon
          dy = Observations(i)%lat - Observations(j)%lat
          
          tmpreal = get_gaspari_cohn( dx, dy, xcompact, ycompact ) 
          
          Repr_matrix(i,j) = tmpreal * Repr_matrix(i,j)
          Repr_matrix(j,i) = tmpreal * Repr_matrix(j,i)
          
       end do
    end do

  end subroutine hadamard_for_repr_matrix
  

!BOP
! !ROUTINE: gaspari_cohn
! \label{gaspari_cohn_pf}
!   
! !INTERFACE:  
  function gaspari_cohn( d )
! 
! !DESCRIPTION: 
! evaluate 5th-order polynomial from Gaspari and Cohn, 1999, Eq (4.10)
!
! On input, d = separation distance relative to the distance at which 
! all correlations vanish. In the isotropic case, Gaspari and Cohn, 1999,
! Eq. (4.10)
!
! \begin{verbatim}
!    d = sqrt(dx**2 + dy**2) / (2*c) = |z| / (2*c)    
! \end{verbatim}
!
! or in the anisotropic case
!
! \begin{verbatim}
!    d = sqrt( (dx/xcompact)**2 + (dy/ycompact)**2 )
! \end{verbatim}
!
! \begin{verbatim}
! *** Use |z|/c = 2*d. All correlations vanish for d > 1. ***
! \end{verbatim}
!EOP
    
    implicit none
    
    real :: gaspari_cohn, d, y

    real, parameter :: tol = 1e-3
    
    ! Get rid of possibly negative distances.
    ! Multiply with 2. to return to the Gaspari and Cohn, 1999, notation.
    
    d = 2.*abs(d)
    
    if (d >= 2.) then
       
       y = 0.
       
    else if (d <= tol) then              
       
       y = 1.
       
    else if (d <= 1.) then
       
       ! y = -.25*d**5 + .5*d**4 + .625*d**3 - 5./3.*d**2 + 1.
       
       y = d**2 *( d*( d*( -.25*d + .5) + .625) -5./3.) + 1.
       
    else
       
       ! y = d**5/12. - .5*d**4 + .625*d**3 + 5./3.*d**2 - 5.*d + 4. - 2./3./d
       
       y = d*( d*( d*( d*( d/12. - .5) + .625) + 5./3.) -5.) + 4. - 2./3./d
       
    end if
    
    gaspari_cohn = y
    
    return
    
  end function gaspari_cohn
  
!BOP
! 
! !ROUTINE: get_gaspari_cohn
!  \label{get_gaspari_cohn_pf}
! 
! !INTERFACE:   
  function get_gaspari_cohn( dx, dy, xcompact, ycompact )

! !DESCRIPTION:    
! for a given lat/lon distance, compute the
! anisotropic compact support (Gaspari and Cohn) weights
!
! input coordinates must be in degrees latitude/longitude
!  
!  dx = longitude separation of two points \newline
!  dy = latitude separation of two points
!
!  xcompact = longitude scale of compact support \newline
!  ycompact = latitude scale of compact support
!
! for the isotropic Gaspari and Cohn function, the relative
! distance is in Gaspari and Cohn, 1999, notation (Eq (4.10))
!
!      d = sqrt(dx**2 + dy**2) / 2*c
!
! get\_gaspari\_cohn() uses a generalized anisotropic Gaspari and Cohn
! approach (essentially coordinate stretching, see handwritten
! notes for details)
!
!  d = sqrt((dx/xcompact)**2 + (dy/ycompact)**2  )
!
! All correlations vanish outside of an ellipse with semi-axes 
! xcompact and ycompact, ie Gaspari and Cohn weights vanish 
! for d $>$ 1 (note the factor 2!)
!
! When the anisotropic case is reduced back to the isotropic case,
! (ie if xcompact==ycompact) then c = xcompact/2 = ycompact/2.
!
!EOP
            
    implicit none
    
    real :: get_gaspari_cohn, dx, dy, xcompact, ycompact, d

    ! compute (anisotropic) distance relative to compact support
    
    d = sqrt( (dx/xcompact)**2 + (dy/ycompact)**2 )    
    get_gaspari_cohn = gaspari_cohn(d)

  end function get_gaspari_cohn
  
  ! ************************************************************


!BOP
! 
! !ROUTINE: giinvert
!
! !INTERFACE:
!
  subroutine gjinvert(n, A_orig,A_inv)
! !DESCRIPTION:  
! Gauss-Jordan Method for Matrix Inversion
!
! Inversion of n-square matrix
!
! Loosely based on original code by David Goodrich, USDA ARS, for use with 
!   MQ-B program by Haider (1994)
!
! References
!   Haider, S.K., 1994:  Spatial storm characteristics and basin response.  
!      M.S. Thesis, Dept. of Hydrology and Water Resources, University of 
!      Arizona.  261 pp.
!
! 
!EOP

    integer, intent(IN) :: n
    real, intent(IN) :: A_orig(n,n)
    real, intent(OUT) :: A_inv(n,n)

! local variables
    real :: big,dum,pivinv
    real :: a(n,n),b(n,1)
    integer :: irow,icol
    integer :: i,j,k,l,ll
    integer :: m = 0
    integer :: ipiv(n),indxr(n),indxc(n)

    a = A_orig
    ipiv = 0
    do i = 1,n
       big = 0
       do j = 1,n
          if (ipiv(j).ne.1) then
             do k = 1,n
                if (ipiv(k).eq.0) then
                   if (abs(a(j,k)).ge.big) then
                      big = abs(a(j,k))
                      irow = j
                      icol = k
                   endif
                else if (ipiv(k).gt.1) then
                   write(LIS_logunit,*) '[ERR] mqb_module -- error in matrix inversion'
                   call LIS_endrun()
                endif
             end do
          endif
       end do
       !
       ipiv(icol) = ipiv(icol)+1
       if (irow.ne.icol) then
          do l = 1,n
             dum = a(irow,l)
             a(irow,l) = a(icol,l)
             a(icol,l) = dum
          end do
          do l = 1,m
             dum = b(irow,l)
             b(irow,l) = b(icol,l)
             b(icol,l) = dum
          end do
       endif
       indxr(i)=irow
       indxc(i)=icol
       if (a(icol,icol).eq.0) then
          write(LIS_logunit,*) '[ERR] mqb_module -- singular matrix encountered'
          call LIS_endrun()
       endif
       !
       pivinv = 1 / a(icol,icol)
       a(icol,icol) = 1
       do l = 1,n
          a(icol,l) = a(icol,l) * pivinv
       end do
       do l = 1,m
          b(icol,l) = b(icol,l) * pivinv
       end do
       !
       do ll = 1,n
          if (ll.ne.icol) then
             dum = a(ll,icol)
             a(ll,icol) = 0
             do l = 1,n
                a(ll,l) = a(ll,l) - a(icol,l) * dum
             end do
             do l = 1,m
                b(ll,l) = b(ll,l) - b(icol,l) * dum
             end do
          endif
       end do
    end do
    !
    do j = 1,n
       l = n + 1 - j
       if (indxr(l).ne.indxc(l)) then
          do k = 1,n
             dum = a(k,indxr(l))
             a(k,indxr(l)) = a(k,indxc(l))
             a(k,indxc(l)) = dum
          end do
       endif
    end do
    !
    A_inv = a
    return
  end subroutine gjinvert
  
  recursive function determinant(matrix) result(laplace_det)
    real, dimension(:,:) :: matrix
    integer :: msize(2), i, n
    real :: laplace_det, det
    real, dimension(:,:), allocatable :: cf
    
    msize = shape(matrix) 
    n = msize(1)          
    
    if (n .eq. 1) then  
       det = matrix(1,1)
    else
       det = 0    
       do i=1, n  
          allocate(cf(n-1, n-1))     
          cf = cofactor(matrix, i, 1)
          det = det + ((-1)**(i+1))* matrix(i,1) * determinant(cf)
          deallocate(cf)
       end do
    end if
    laplace_det = det
  end function determinant
    
  function cofactor(matrix, mI, mJ)
    real, dimension(:,:) :: matrix
    integer :: mI, mJ
    integer :: msize(2), i, j, k, l, n
    real, dimension(:,:), allocatable :: cofactor
    msize = shape(matrix)
    n = msize(1)
    
    allocate(cofactor(n-1, n-1))
    l=0
    k = 1
    do i=1, n
       if (i .ne. mI) then
          l = 1
          do j=1, n
             if (j .ne. mJ) then
                cofactor(k,l) = matrix(i,j)
                l = l+ 1
             end if
          end do
          k = k+ 1
       end if
    end do
    return
  end function cofactor

end module pf_general
