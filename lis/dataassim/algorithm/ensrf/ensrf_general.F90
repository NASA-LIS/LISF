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
! !MODULE: ensrf_general
! 
! this file contains a collection of general Ensemble Kalman filter
! subroutines and compact support subroutines
!
! !REVISION HISTORY: 
! reichle, 20 Apr 2001
! reichle, 18 Mar 2004 - optional arguments
! reichle, 27 Jan 2005 - eliminated use of module select_kinds
! reichle, 19 Jul 2005 - merged compact_support.f90 and ensrf_general.f90
! reichle,  1 Aug 2005 - eliminated tile_coord
! reichle, 18 Oct 2005 - return increments instead of updated State
!EOP
module ensrf_general
  
  implicit none
  
  private
  
  public :: ensrf_analysis
  
contains
  
!BOP
! 
! !ROUTINE: ensrf_analysis
! \label{ensrf_analysis}
! 
! !INTERFACE:  
  subroutine ensrf_analysis( gid, &
       N_state, N_obs, N_ens, &
       Observations, Obs_pred, Obs_err, Obs_cov, &
       State_incr, &
       State_lon, State_lat, xcompact, ycompact )
! !USES:    
    use ensrf_types
    use my_matrix_functions
    use LIS_logMod, only : LIS_logunit, LIS_endrun

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
! perform Ensrf update
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
    real, dimension(N_state,N_ens) :: State_vec
    real, dimension(N_state,N_ens) :: State_post
    real, dimension(N_state,N_ens) :: State_prime
    real, dimension(N_state)       :: State_bar
    real, dimension(N_state)       :: State_incr_tmp
    
    real, dimension(N_obs,N_ens)   :: Obs_pred_prime
    real, dimension(N_obs)         :: Obs_pred_bar
    real, dimension(N_obs,N_ens)         :: rhs
    
    real, dimension(N_ens)         :: weights    
    
    real, dimension(N_obs,N_obs)   :: Repr_matrix
    real, dimension(N_obs,N_obs)   :: Repr_matrix_inv
    real, dimension(N_obs,N_obs)   :: Repr_matrix_mult
    real, dimension(N_state, N_obs) :: PHt
    real, dimension(N_state, N_obs) :: Kgain

    integer,          dimension(N_obs)         :: indx
    real, dimension(N_state) :: ens_mean
    real    :: dweight, alpha,rho
    logical :: apply_hadamard
    
    ! ------------------------------------------------------------------

!    if(present(xcompact).and.present(ycompact)) then 
!       write(LIS_logunit,*) 'Currently only 1-d filter is supported...'
!       write(LIS_logunit,*) 'Program stopping...'
!       call LIS_endrun()
!    endif

    ! find out whether Hadamard product should be applied

    apply_hadamard = (               &
         present(State_lon)    .and. &
         present(State_lat)    .and. &
         present(xcompact)     .and. &
         present(ycompact)              )
    
    ! ----------------------
    
    rho = 1.10
    State_vec = State_incr
    print*, State_incr(2,:)
    ! IMPORTANT: on input, State_incr contains State_minus(1:N_state,1:N_ens)
    
    ! compute ensemble mean Ybar at current update time
    
    State_bar = sum( State_incr, 2) / real(N_ens)

    ! finalize matrix Y_prime = Y - Ybar
    
    do n_e=1,N_ens
       
       State_prime(:,n_e) = rho*(State_incr(:,n_e) - State_bar)
       
    end do
        
    ! --------------------
    
    ! compute ensemble mean H*Ybar
    
    Obs_pred_bar = sum( Obs_pred, 2) / real(N_ens)

    ! finalize matrix Q_prime = H*Y - H*ybar
    
    do n_e=1,N_ens
       
       Obs_pred_prime(:,n_e) = Obs_pred(:,n_e) - Obs_pred_bar
       
    end do

    ! --------------------

    ! form repr matrix HPHt = Q_prime*(Q_prime)t/(N_e-1)
    
    Repr_matrix = &
         (matmul(Obs_pred_prime,transpose(Obs_pred_prime))) &
         /real(N_ens-1)

    ! reichle, 18 Mar 2004:
    ! maybe Hadamard product should be applied *after* adding Obs_cov
    ! to representer matrix? only matters if Obs_cov is not diagonal...
    
    if (apply_hadamard) &
         call hadamard_for_repr_matrix( N_obs, Observations, &
         xcompact, ycompact, Repr_Matrix )
    
    ! form matrix W = HPHt+ R
    
    Repr_matrix = Repr_matrix + Obs_cov

    call gjinvert(N_obs, Repr_matrix, Repr_matrix_inv)
    
    Repr_matrix_mult = matmul(Repr_matrix,Repr_matrix_inv)
    ! update each ensemble member

    !assuming scalar for now
    alpha = 1/(1+sqrt(Repr_matrix_mult(1,1)))
    print*, 'alpha',alpha
    do n_e=1,N_ens
       
       ! use random measurement error field, get Zpert = Z + v,
       ! compute right hand side rhs = Zpert - H*y^f for system equation,
       
       do i=1,N_obs
          rhs(i,n_e) = Observations(i)%value  - Obs_pred(i,n_e)
       end do
    enddo

    State_incr = 0.           ! State_incr = analysis - forecast
       
    do ii=1,N_state
       
       do jj=1,N_obs
          
          ! compute [PHt]_ij  (normalize later)
          
          PHt_ij = 0.
          
          do kk=1,N_ens
             
             PHt_ij = PHt_ij + State_prime(ii,kk)*Obs_pred_prime(jj,kk)
             
          end do
          
          ! apply Hadamard factor
          
          dx=State_lon(ii)-Observations(jj)%lon
          dy=State_lat(ii)-Observations(jj)%lat
          
          ! multiply [PHt]_ij with Hadamard factor
          
          dweight = get_gaspari_cohn( dx, dy, xcompact, ycompact )
          
          if(ii.eq.2) then 
             print*, sqrt(dx**2+dy**2), dweight
          endif

          PHt(ii,jj)= PHt_ij * dweight /real(N_ens-1)   
          
       end do
       
    end do
       
    Kgain = matmul(PHt,Repr_matrix_inv)
    
    !bar(x^+) = bar(x^-) + K(Y-H bar(x^-))
    !The N_ens member is assumed to be unperturbed, representing the
    !mean
    ens_mean(:) = matmul(Kgain,rhs(:,N_ens))
    
    !update perturbations
    do n_e = 1,N_ens
       State_post(:,n_e) = State_prime(:,n_e) - alpha*matmul(Kgain,rhs(:,n_e))
    enddo

    do n_e = 1,N_ens
       State_incr(:,n_e) = ens_mean(:) + State_post(:,n_e) - state_vec(:,n_e)
    enddo

  end subroutine ensrf_analysis
#if 0 
!BOP
! 
! !ROUTINE: ensrf_analysis
! \label{ensrf_analysis}
! 
! !INTERFACE:  
  subroutine ensrf_analysis( gid, &
       N_state, N_obs, N_ens, &
       Observations, Obs_pred, Obs_err, Obs_cov, &
       State_incr, &
       State_lon, State_lat, xcompact, ycompact )
! !USES:    
    use ensrf_types
    use my_matrix_functions
    use LIS_logMod, only : LIS_logunit, LIS_endrun

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
! perform Ensrf update
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
    real, dimension(N_obs)         :: rhs
    
    real, dimension(N_ens)         :: weights    
    
    real, dimension(N_obs,N_obs)   :: Repr_matrix
    
    integer,          dimension(N_obs)         :: indx
    
    real    :: dweight
    logical :: apply_hadamard
    
    ! ------------------------------------------------------------------

!    if(present(xcompact).and.present(ycompact)) then 
!       write(LIS_logunit,*) 'Currently only 1-d filter is supported...'
!       write(LIS_logunit,*) 'Program stopping...'
!       call LIS_endrun()
!    endif

    ! find out whether Hadamard product should be applied

    apply_hadamard = (               &
         present(State_lon)    .and. &
         present(State_lat)    .and. &
         present(xcompact)     .and. &
         present(ycompact)              )
    
    ! ----------------------

    print*, State_incr(2,:)
    ! IMPORTANT: on input, State_incr contains State_minus(1:N_state,1:N_ens)
    
    ! compute ensemble mean Ybar at current update time
    
    State_bar = sum( State_incr, 2) / real(N_ens)

    ! finalize matrix Y_prime = Y - Ybar
    
    do n_e=1,N_ens
       
       State_prime(:,n_e) = State_incr(:,n_e) - State_bar
       
    end do
        
    ! --------------------
    
    ! compute ensemble mean H*Ybar
    
    Obs_pred_bar = sum( Obs_pred, 2) / real(N_ens)

    ! finalize matrix Q_prime = H*Y - H*ybar
    
    do n_e=1,N_ens
       
       Obs_pred_prime(:,n_e) = Obs_pred(:,n_e) - Obs_pred_bar
       
    end do

    ! --------------------

    ! form repr matrix HPHt = Q_prime*(Q_prime)t/(N_e-1)
    
    Repr_matrix = &
         (matmul(Obs_pred_prime,transpose(Obs_pred_prime))) &
         /real(N_ens-1)

    ! reichle, 18 Mar 2004:
    ! maybe Hadamard product should be applied *after* adding Obs_cov
    ! to representer matrix? only matters if Obs_cov is not diagonal...
    
    if (apply_hadamard) &
         call hadamard_for_repr_matrix( N_obs, Observations, &
         xcompact, ycompact, Repr_Matrix )
    
    ! form matrix W = HPHt+ R
    
    Repr_matrix = Repr_matrix + Obs_cov

    ! maybe later: save representer matrix into file
    
    ! decompose W once (look into LAPACKs sgelss/dgelss, ask Christian)
    
    call ludcmp( Repr_Matrix, N_obs, indx ) 
   
    ! --------------------------------------
           
    ! update each ensemble member

    do n_e=1,N_ens
       
       ! use random measurement error field, get Zpert = Z + v,
       ! compute right hand side rhs = Zpert - H*y^f for system equation,
       
       do i=1,N_obs
          if(Observations(i)%pert_type.eq.0) then 
             rhs(i) = Observations(i)%value + Obs_err(i,n_e) - Obs_pred(i,n_e)
          elseif(Observations(i)%pert_type.eq.1) then 
             rhs(i) = Observations(i)%value * (1+Obs_err(i,n_e)) - Obs_pred(i,n_e)
          endif
       end do

       ! solve W*b = Zpert - H*y^f
       
       call lubksb( Repr_matrix, N_obs, indx, rhs )
       
       if (.not. apply_hadamard) then
          
          ! update *without* Hadamard product
          
          ! compute w = (Q_prime)t b

          weights = matmul( transpose(Obs_pred_prime), rhs)

          ! compute new initial cond y^a = y^f + (Y-ybar)*w/(N_e-1)
          
          ! start with (Y-ybar)*w, write into Y_Vector
          
          State_incr_tmp    = matmul( State_prime, weights)
          
          State_incr(:,n_e) = State_incr_tmp
!          if(gid.eq.1) then 
!             write(LIS_logunit,*)'State prime ',State_prime
!             write(LIS_logunit,*)'rhs ',rhs
!             write(LIS_logunit,*)'weights ',weights
!             write(LIS_logunit,*)'incr ',State_incr_tmp
!          endif
       else
          
          ! update *with* Hadamard product
          
          State_incr(:,n_e) = 0.           ! State_incr = analysis - forecast
          
          do ii=1,N_state
             
             do jj=1,N_obs
                
                ! compute [PHt]_ij  (normalize later)
                
                PHt_ij = 0.
                
                do kk=1,N_ens
                   
                   PHt_ij = PHt_ij + State_prime(ii,kk)*Obs_pred_prime(jj,kk)

                end do
                
                ! apply Hadamard factor
                   
                dx=State_lon(ii)-Observations(jj)%lon
                dy=State_lat(ii)-Observations(jj)%lat
                
                ! multiply [PHt]_ij with Hadamard factor
                
                dweight = get_gaspari_cohn( dx, dy, xcompact, ycompact )

                PHt_ij = PHt_ij * dweight      
                
                if(ii.eq.2) then 
                   print*, 'st ',n_e,dweight, sqrt(dx**2+dy**2)
                endif
                State_incr(ii,n_e) = State_incr(ii,n_e) + PHt_ij*rhs(jj)

                
             end do
             
          end do
          
       end if
       
       
       ! finish computation of increment for ensemble member n_e
       ! (normalization is NOT ensemble average, see Eq above)
       
       State_incr(:,n_e) = State_incr(:,n_e)/real(N_ens-1)
       
    end do
    
  end subroutine ensrf_analysis
#endif  

!BOP
! 
! !ROUTINE: hadamard_for_repr_matrix
! \label{hadamard_for_repr_matrix_ensrf}
! 
! !INTERFACE:  
  subroutine hadamard_for_repr_matrix( N_obs, Observations,             & 
       xcompact, ycompact,                                              &
       Repr_Matrix )
! !USES:    
    use ensrf_types                     ! for obs_type
    
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
! \label{gaspari_cohn_ensrf}
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
!  \label{get_gaspari_cohn_ensrf}
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

!    real :: dist,maxdist
    ! compute (anisotropic) distance relative to compact support
    
    d = sqrt( (dx/xcompact)**2 + (dy/ycompact)**2 )    
    get_gaspari_cohn = gaspari_cohn(d)

!    dist = dx**2 + dy**2
!    maxdist = xcompact**2 + ycompact**2
!    if(dist.le.maxdist) then 
!       get_gaspari_cohn = (maxdist -dist)/(maxdist+dist)
!    else
!       get_gaspari_cohn = 0.0
!    endif

  end function get_gaspari_cohn
  
  ! ************************************************************

!
!
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
! Version History
!   05 May 2005  Matthew Garcia  Initial specification
!
!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
subroutine gjinvert(n, A_orig,A_inv)
!
! argument variables
  integer, intent(IN) :: n
  real, intent(IN) :: A_orig(n,n)
  real, intent(OUT) :: A_inv(n,n)
!
! local variables
  real :: big,dum,pivinv
  real :: a(n,n),b(n,1)
  integer :: irow,icol
  integer :: i,j,k,l,ll
  integer :: m = 0
  integer :: ipiv(n),indxr(n),indxc(n)
!
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
            print *,'ERR: mqb_module -- error in matrix inversion'
            stop
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
      print *,'ERR: mqb_module -- singular matrix encountered'
      stop
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

end module ensrf_general
