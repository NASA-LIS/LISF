MODULE cable_sli_roots

  USE cable_dimensions, ONLY: r_2, i_d
  USE cable_sli_numbers,       ONLY: zero, half, one, two, e3

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: setroots, getrex ! routines

  ! b1, b2 - params used to get root distribution param b (Li et al., J Hydrolo 2001).
  REAL(r_2), PARAMETER :: b1     = 24.66
  REAL(r_2), PARAMETER :: b2     = 1.59
  REAL(r_2), PARAMETER :: lambda = 1.0

  INTERFACE setroots
     MODULE PROCEDURE setroots_1d
     MODULE PROCEDURE setroots_2d
  END INTERFACE setroots

  INTERFACE getrex
     MODULE PROCEDURE getrex_1d
     MODULE PROCEDURE getrex_2d
  END INTERFACE getrex

  ! Definitions:
  ! setroots - subroutine to set current root distribution based on Li, K.Y., De Jong, R. and J.B. Boisvert (2001).
  !            An exponential root-water-uptake model with water stress compensation. J. Hydrol. 252:189-204.
  ! getrex   - subroutine to get rate of water extraction from layers.

  !**********************************************************************************************************************

CONTAINS

  !**********************************************************************************************************************

  SUBROUTINE setroots_1d(x, F10, Zr, Fs)
 
    IMPLICIT NONE

    REAL(r_2), DIMENSION(:), INTENT(IN)  :: x
    REAL(r_2),               INTENT(IN)  :: F10
    REAL(r_2),               INTENT(IN)  :: Zr
    REAL(r_2), DIMENSION(:), INTENT(OUT) :: Fs
    ! Sets current weighted root length density distribn (Fs).
    !      Li et al. (J Hydrolo, 2001)
    ! Definitions of arguments:
    ! x(:) - depths to bottom of layers (cm).
    ! F10  - fraction of root length density in top 10% of the root zone
    ! Zr   - rooting depth (cm).
    INTEGER(i_d)                    :: ms
    REAL(r_2)                       :: b, extr
    REAL(r_2), DIMENSION(1:size(x)) :: ext0, ext1, xend, Fi
    REAL(r_2)                       :: tmp

    ms = size(x) ! soil layers

    xend(1)    = zero
    xend(2:ms) = x(1:ms-1)
    b          = b1 * F10**b2 / Zr ! root distrib param
    extr       = exp(-b*Zr)
    ext1(:)    = exp(-b*x(:))
    ext0(1)    = one
    ext0(2:ms) = ext1(1:ms-1)
    ! get fraction of root length density in layer i
    Fi(:) = (log((one+ext0(:))/(one+ext1(:))) + half*(ext0(:)-ext1(:))) / &
         (log(two/(one+extr)) + half*(one-extr))
    Fs(:) = exp(lambda*log(Fi(:))) ! weighted Fi
    where (xend(:) >= Zr) Fs(:) = zero

    ! ensure Fs sums to one
    tmp = sum(Fs(:))
    if (tmp > zero) then
       Fs(:) = Fs(:)/tmp
    else
       Fs(:) = zero
    endif

  END SUBROUTINE setroots_1d

  SUBROUTINE setroots_2d(x, F10, Zr, Fs)
 
    IMPLICIT NONE

    REAL(r_2), DIMENSION(:,:), INTENT(IN)  :: x
    REAL(r_2), DIMENSION(:),   INTENT(IN)  :: F10
    REAL(r_2), DIMENSION(:),   INTENT(IN)  :: Zr
    REAL(r_2), DIMENSION(:,:), INTENT(OUT) :: Fs
    ! Sets current weighted root length density distribn (Fs).
    !      Li et al. (J Hydrolo, 2001)
    ! Definitions of arguments:
    ! x(:) - depths to bottom of layers (cm).
    ! F10  - fraction of root length density in top 10% of the root zone
    ! Zr   - rooting depth (cm).
    INTEGER(i_d)                                  :: mp, ms
    REAL(r_2), DIMENSION(1:size(x,1))             :: b, extr
    REAL(r_2), DIMENSION(1:size(x,1),1:size(x,2)) :: ext0, ext1, xend, Fi
    REAL(r_2), DIMENSION(1:size(x,1),1:size(x,2)) :: tmp2d

    mp = size(x,1) ! landpoints
    ms = size(x,2) ! soil layers

    xend(1:mp,1)    = zero
    xend(1:mp,2:ms) = x(1:mp,1:ms-1)
    b(1:mp)         = b1 * F10(1:mp)**b2 / Zr(1:mp) ! root distrib param
    extr(1:mp)      = exp(-b(1:mp)*Zr(1:mp))
    tmp2d(:,:)      = spread(b(1:mp),2,ms)
    ext1(:,:)       = exp(-tmp2d(:,:)*x(:,:))
    ext0(1:mp,1)    = one
    ext0(1:mp,2:ms) = ext1(1:mp,1:ms-1)
    ! get fraction of root length density in layer i
    tmp2d(:,:)      = spread(extr(1:mp),2,ms)
    Fi(:,:)         = (log((one+ext0(:,:))/(one+ext1(:,:))) + half*(ext0(:,:)-ext1(:,:))) / &
         (log(two/(one+tmp2d(:,:))) + half*(one-tmp2d(:,:)))
    Fs(:,:)         = exp(lambda*log(Fi(:,:))) ! weighted Fi
    tmp2d(:,:)      = spread(Zr(1:mp),2,ms)
    where (xend(:,:) >= tmp2d(:,:)) Fs(:,:) = zero

    ! ensure Fs sums to one
    tmp2d(:,:) = spread(sum(Fs(1:mp,1:ms),2),2,ms)
    where (tmp2d(:,:) > zero)
       Fs(:,:) = Fs(:,:)/tmp2d(:,:)
    elsewhere
       Fs(:,:) = zero
    endwhere

  END SUBROUTINE setroots_2d

  !**********************************************************************************************************************

  SUBROUTINE getrex_1d(S, rex, fws, Fs, thetaS, thetaw, Etrans, gamma)
 
    ! Lai and Katul formulation for root efficiency function vh 17/07/09
    ! changed to MCs maximum formulation

    IMPLICIT NONE
    
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: S      ! relative saturation
    REAL(r_2), DIMENSION(:), INTENT(OUT)   :: rex    ! water extraction per layer
    REAL(r_2),               INTENT(INOUT) :: fws    ! stomatal limitation factor
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: Fs     ! root length density
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: thetaS ! saturation soil moisture
    REAL(r_2), DIMENSION(:), INTENT(IN)    :: thetaw ! soil moisture at wiliting point
    REAL(r_2),               INTENT(IN)    :: Etrans ! total transpiration
    REAL(r_2),               INTENT(IN)    :: gamma  ! skew of Li & Katul alpha2 function

    ! Gets rate of water extraction compatible with CABLE stomatal conductance model 
    ! theta(:) - soil moisture(m3 m-3)
    ! rex(:)   - rate of water extraction by roots from layers (cm/h).
    INTEGER(i_d)                                  :: ms
    REAL(r_2), DIMENSION(1:size(S)) :: theta, lthetar, alpha_root, delta_root
    REAL(r_2)                       :: trex

    ms = size(S) ! soil layers

    theta(:)   = S(:)*thetaS(:)
    lthetar(:) = log(max(theta(:)-thetaw(:),e3)/thetaS(:))

    where ((theta(:)-thetaw(:)) > e3)
       alpha_root(:) = exp( gamma/max(theta(:)-thetaw(:), e3) * lthetar(:) )
    elsewhere
       alpha_root(:) = zero
    endwhere

    where (Fs(:) > zero)
        delta_root(:) = one
    elsewhere
        delta_root(:) = zero
    endwhere

    rex(:) = alpha_root(:)*Fs(:)

    trex = sum(rex(:))
    if (trex > zero) then
       rex(:) = rex(:)/trex
    else
       rex(:) = zero
    endif
    rex(:) = Etrans*rex(:)

    fws    = maxval(alpha_root(:)*delta_root(:))

  END SUBROUTINE getrex_1d

  SUBROUTINE getrex_2d(S, rex, fws, Fs, thetaS, thetaw, Etrans, gamma)
 
    ! Lai and Katul formulation for root efficiency function vh 17/07/09
    ! changed to MCs maximum formulation

    IMPLICIT NONE
    
    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: S      ! relative saturation
    REAL(r_2), DIMENSION(:,:), INTENT(OUT)   :: rex    ! water extraction per layer
    REAL(r_2), DIMENSION(:),   INTENT(INOUT) :: fws    ! stomatal limitation factor
    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: Fs     ! root length density
    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: thetaS ! saturation soil moisture
    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: thetaw ! soil moisture at wiliting point
    REAL(r_2), DIMENSION(:),   INTENT(IN)    :: Etrans ! total transpiration
    REAL(r_2), DIMENSION(:),   INTENT(IN)    :: gamma  ! skew of Li & Katul alpha2 function

    ! Gets rate of water extraction compatible with CABLE stomatal conductance model 
    ! theta(:) - soil moisture(m3 m-3)
    ! rex(:)   - rate of water extraction by roots from layers (cm/h).
    INTEGER(i_d)                                  :: mp, ms
    REAL(r_2), DIMENSION(1:size(S,1),1:size(S,2)) :: theta, lthetar, alpha_root, delta_root
    REAL(r_2), DIMENSION(1:size(S,1))             :: trex
    REAL(r_2), DIMENSION(1:size(S,1),1:size(S,2)) :: tmp2d

    mp = size(S,1) ! landpoints
    ms = size(S,2) ! soil layers

    theta(:,:)   = S(:,:)*thetaS(:,:)
    lthetar(:,:) = log(max(theta(:,:)-thetaw(:,:),e3)/thetaS(:,:))

    tmp2d = spread(gamma(1:mp),2,ms)
    where ((theta(:,:)-thetaw(:,:)) > e3)
       alpha_root(:,:) = exp( tmp2d(:,:)/max(theta(:,:)-thetaw(:,:), e3) * lthetar(:,:) )
    elsewhere
       alpha_root(:,:) = zero
    endwhere

    where (Fs(:,:) > zero)
        delta_root(:,:) = one
    elsewhere
        delta_root(:,:) = zero
    endwhere

    rex(:,:) = alpha_root(:,:)*Fs(:,:)

    trex  = sum(rex(:,:),2)
    tmp2d = spread(trex(1:mp),2,ms)
    where (tmp2d(:,:) > zero)
       rex(:,:) = rex(:,:)/tmp2d(:,:)
    elsewhere
       rex(:,:) = zero
    endwhere
    tmp2d    = spread(Etrans(1:mp),2,ms)
    rex(:,:) = tmp2d(:,:)*rex(:,:)

    fws(:)   = maxval(alpha_root(:,:)*delta_root(:,:),dim=2)

  END SUBROUTINE getrex_2d

!**********************************************************************************************************************

END MODULE cable_sli_roots
