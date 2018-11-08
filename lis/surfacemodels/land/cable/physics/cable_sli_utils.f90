MODULE cable_sli_utils

  USE cable_dimensions, ONLY: r_2, i_d, SP
  USE cable_types,      ONLY: soil_parameter_type
  USE cable_sli_numbers,       ONLY: zero, half, one, two, four, e3, &
                               Tzero, gravity, Rgas, thousand, Mw, rlambda, Dva, &
                               params, vars_aquifer, vars, rapointer, &
                               rhow, rhoi, lambdaf, lambdas, csice, cswat, &
                               rhocp, &
                               dpmaxr, solve_type, &
                               gf, hmin, csol, rhmin, dsmmax

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dx, dxL, par, plit, sol, x ! soil water parameters
  PUBLIC :: bd, dis, isopar, isotype   ! soil solute parameters
  PUBLIC :: aquifer_props, flux, generic_thomas, getfluxes_vp, getheatfluxes, hyofh, hyofS, isosub ! subroutines
  PUBLIC :: litter_props, massman_sparse, potential_evap, setlitterpar, setpar, setsol, setx, tri
  PUBLIC :: csat, csoil, dthetalmaxdT, dthetalmaxdTh, esat, gammln, igamma, phi, rh0_sol, rtbis_rh0 ! functions
  PUBLIC :: slope_csat, slope_esat, Sofh, Tfrz, thetalmax, weight, zerovars

  REAL(r_2),        DIMENSION(:,:),   ALLOCATABLE :: dx
  REAL(r_2),        DIMENSION(:),     ALLOCATABLE :: dxL
  TYPE(params),     DIMENSION(:,:),   ALLOCATABLE :: par
  TYPE(params),     DIMENSION(:),     ALLOCATABLE :: plit
  TYPE(solve_type), DIMENSION(:),     ALLOCATABLE :: sol
  REAL(r_2),        DIMENSION(:,:),   ALLOCATABLE :: x

  !MC! solute not done yet
  !REAL(r_2),        DIMENSION(:,:),   ALLOCATABLE         :: bd
  !REAL(r_2),        DIMENSION(:,:),   ALLOCATABLE         :: dis
  !TYPE(rapointer),  DIMENSION(:,:,:), ALLOCATABLE         :: isopar
  !CHARACTER(LEN=2), DIMENSION(:,:,:), ALLOCATABLE         :: isotype
  REAL(r_2),        DIMENSION(:),   ALLOCATABLE         :: bd ! function solute not done yet
  REAL(r_2),        DIMENSION(:),   ALLOCATABLE         :: dis
  TYPE(rapointer),  DIMENSION(:,:), ALLOCATABLE         :: isopar
  CHARACTER(LEN=2), DIMENSION(:,:), ALLOCATABLE         :: isotype

  ! Subroutine interfaces

  INTERFACE generic_thomas
     MODULE PROCEDURE generic_thomas_1d
     MODULE PROCEDURE generic_thomas_2d
  END INTERFACE generic_thomas

  INTERFACE getfluxes_vp
     MODULE PROCEDURE getfluxes_vp_1d
     MODULE PROCEDURE getfluxes_vp_2d
  END INTERFACE getfluxes_vp

  INTERFACE getheatfluxes
     MODULE PROCEDURE getheatfluxes_1d
     MODULE PROCEDURE getheatfluxes_2d
  END INTERFACE getheatfluxes

  INTERFACE massman_sparse
     MODULE PROCEDURE massman_sparse_1d
     MODULE PROCEDURE massman_sparse_2d
  END INTERFACE massman_sparse

  INTERFACE tri
     MODULE PROCEDURE tri_1d
     MODULE PROCEDURE tri_2d
  END INTERFACE tri

  ! Function interfaces

  INTERFACE gammln
     MODULE PROCEDURE sgammln, dgammln
  END INTERFACE gammln

  INTERFACE gcf
     MODULE PROCEDURE sgcf, dgcf
  END INTERFACE gcf

  INTERFACE gser
     MODULE PROCEDURE sgser, dgser
  END INTERFACE gser

  INTERFACE igamma
     MODULE PROCEDURE sigamma, digamma
  END INTERFACE igamma

  !P.J. Ross 2005-2007:
  ! This module implements Brooks-Corey (BC) soil water retention and conductivity
  ! functions. It illustrates the structure required without the complicating
  ! detail of more flexible hydraulic property functions.
  ! Definitions of public entities (see above for default values):
  ! params  - type for water parameters. Params the, thre (=the-thr), he, lam,
  !           Ke and eta are for the BC functions. Params he, Ke, KSe, phie
  !           and phiSe are values of variables h, K, KS, phi and phiS at
  !           saturation (denoted by the "e" for "air entry"), needed by
  !           module flow (MF).
  ! vars    - type for water variables used by MF and returned by subroutine
  !           hyofS (except for isat, which is 0 for unsaturated layers and 1
  !           for saturated layers).
  ! gf      - gravity factor for flow direction (usually 1 for straight down).
  ! hmin    - minimum matric head h (used by MF).
  ! par(:)  - hydraulic property params for soil types (used by MF).
  ! allo    - subroutine to allocate parameter storage.
  ! hyofS   - subroutine to get water variable from saturation S (where S<1).
  ! hyofh   - subroutine to get some water variables from h.
  ! Sofh    - subroutine to get S from h.
  ! hypar   - subroutine to set soil hydraulic params.
  ! weight  - subroutine to get gravity flow conductivity weight w.

  !V. Haverd 2008:
  ! modified by vh to include isothermal vapour conductivity in phi
  ! change units from cm, h to m, s
  ! includes thermal properties in soil parameters
  ! litter_props - subroutine to define litter properties

  !**********************************************************************************************************************

CONTAINS

  !**********************************************************************************************************************
  ! SUBROUTINES - alphabetical
  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE aquifer_props(v_aquifer)

    IMPLICIT NONE

    TYPE(vars_aquifer), INTENT(INOUT) :: v_aquifer

    v_aquifer%zzero     = 53.43_r_2  ! water table depth corresponding to Wa = zero
    v_aquifer%Sy        = 0.2_r_2  ! specific yield of aquifer
    ! initialise water content of aquifer
    v_aquifer%Wa        = v_aquifer%Sy*(v_aquifer%zzero-max(v_aquifer%zdelta,v_aquifer%zsoil))
    v_aquifer%isat      = 0
    if (v_aquifer%zdelta <= v_aquifer%zsoil) v_aquifer%isat = 1
    v_aquifer%f         = 1.25_r_2 ! multiplier in exponent of Rs (m-1)
    v_aquifer%Rsmax     = 4.5e-7_r_2  ! maximum discharge rate from aquifer (ms-1)
    v_aquifer%discharge = v_aquifer%Rsmax*exp(-v_aquifer%f*v_aquifer%zdelta)

  END SUBROUTINE aquifer_props

  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE flux(parin, ns1j0, v1, v2, dz, q, qya, qyb, qTa, qTb)
    ! VH modified 25/05/10 to include t-dep component of liquid flux (in frozen soil)

    IMPLICIT NONE

    TYPE(params), INTENT(IN)  :: parin
    LOGICAL,      INTENT(IN)  :: ns1j0
    REAL(r_2),    INTENT(IN)  :: dz
    REAL(r_2),    INTENT(OUT) :: q, qya, qyb, qTa, qTb
    TYPE(vars),   INTENT(IN)  :: v1, v2
    ! Gets flux and partial derivs for specified flow path.
    ! Definitions of arguments:
    ! j   - soil type no.
    ! v1  - water vars at upper end of path.
    ! v2  - ditto at lower end.
    ! dz  - length of path.
    ! q   - flux.
    ! qya - partial deriv of flux wrt S (if unsat) or phi (if sat) at upper end.
    ! qyb - ditto at lower end.
    REAL(r_2) :: w, rdz

    ! gf is gravity factor (0 to 1) assumed available in module
    if (gf < zero) then
       if ((v1%isat /= 0 .and. v2%isat /= 0) .or. v1%h-gf*(-dz) >= parin%he) then
          w = zero
       else
          w = weight(parin, v1%h, v1%K, v1%phi, -dz)
          w = one-w
       end if
    else
       if ((v1%isat /= 0 .and. v2%isat /= 0) .or. v2%h-gf*dz >= parin%he) then
          w = zero
       else
          w = weight(parin, v2%h, v2%K, v2%phi, dz)
       end if
    end if

    rdz = one/dz
    q   = (v1%phi-v2%phi)*rdz + gf*(w*v1%K+(one-w)*v2%K)
    if (v1%isat==0) then
       qya = v1%phiS*rdz + gf*w*v1%KS
       qTa = v1%phiT*rdz + gf*w*v1%KT
    else
       qya = rdz
       qTa = zero
    end if

    if (v2%isat==0) then
       qyb = -v2%phiS*rdz + gf*(one-w)*v2%KS
       qTb = -v2%phiT*rdz + gf*(one-w)*v2%KT
    else
       qyb = -rdz
       qTb = zero
    end if

    !MC! should gravitational be included?
    if (ns1j0) then          ! ignore gravitational flux in top half of top layer (non-ponding)
       q = ((v1%phi)-(v2%phi))*rdz !+gf*(w*v1%K+(one-w)*v2%K)
       if (v1%isat==0) then
          qya = (v1%phiS)*rdz   !+gf*w*v1%KS
          qTa = zero
       else
          qya = rdz
          qTa = zero
       end if
       if (v2%isat==0) then
          qyb = -(v2%phiS)*rdz ! +gf*(1.-w)*v2%KS
          qTb = zero
       else
          qyb = -rdz
          qTb = zero
       end if
    endif

    !MC! delete ?
    ! special for frozen soil
    !if (v1%K.eq.zero.or.v2%K.eq.zero) then
    !    q = zero
    !    qya = zero
    !    qTa = zero
    !    qyb = zero
    !    qTb = zero
    !endif

  END SUBROUTINE flux

  !**********************************************************************************************************************

  SUBROUTINE generic_thomas_1d(n,A,B,C,r,u)

    USE cable_sli_numbers,       ONLY: one

    IMPLICIT NONE

    ! in/out
    INTEGER(i_d),                      INTENT(IN)  :: n
    REAL(r_2), DIMENSION(1:n,1:2,1:2), INTENT(IN)  :: A, B, C
    REAL(r_2), DIMENSION(1:n,1:2),     INTENT(IN)  :: r
    REAL(r_2), DIMENSION(1:n,1:2),     INTENT(OUT) :: u
    ! local
    REAL(r_2), DIMENSION(1:n,1:2,1:2) :: G
    REAL(r_2), DIMENSION(1:2,1:2)     :: bet
    REAL(r_2), DIMENSION(1:2)         :: d
    REAL(r_2)                         :: detbet, detbet1
    INTEGER(i_d)                      :: j

    ! j=1
    bet(1:2,1:2) = B(1,1:2,1:2)
    detbet       = bet(1,1)*bet(2,2) - bet(1,2)*bet(2,1)
    if (abs(detbet) < epsilon(detbet)) then
       write(2,*) 'generic_thomas_1d error1: det = 0'
       stop 'program terminated by generic_thomas_1d'
    endif
    detbet1 = one/detbet
    d(1:2)  = r(1,1:2)
    u(1,1)  = (bet(2,2)*d(1) - bet(1,2)*d(2))*detbet1
    u(1,2)  = (bet(1,1)*d(2) - bet(2,1)*d(1))*detbet1
    ! j=2, n
    do j=2, n
       d(1:2)       = C(j-1,1:2,1)
       G(j,1,1)     = (bet(2,2)*d(1) - bet(1,2)*d(2))*detbet1
       G(j,2,1)     = (bet(1,1)*d(2) - bet(2,1)*d(1))*detbet1
       d(1:2)       = C(j-1,1:2,2)
       G(j,1,2)     = (bet(2,2)*d(1) - bet(1,2)*d(2))*detbet1
       G(j,2,2)     = (bet(1,1)*d(2) - bet(2,1)*d(1))*detbet1
       bet(1:2,1:2) = B(j,1:2,1:2) - matmul(A(j,1:2,1:2),G(j,1:2,1:2))
       detbet       = bet(1,1)*bet(2,2) - bet(1,2)*bet(2,1)
       if (abs(detbet) < epsilon(detbet)) then
          write(2,*) 'generic_thomas_1d error2: det = 0 at j=', j
          stop 'program terminated by generic_thomas_1d'
       endif
       detbet1      = one/detbet
       d(1:2)       = r(j,1:2) - matmul(A(j,1:2,1:2),u(j-1,1:2))
       u(j,1)       = (bet(2,2)*d(1) - bet(1,2)*d(2))*detbet1
       u(j,2)       = (bet(1,1)*d(2) - bet(2,1)*d(1))*detbet1
    end do
    ! back substitution
    do j=n-1, 1, -1
       u(j,1:2) = u(j,1:2) - matmul(G(j+1,1:2,1:2),u(j+1,1:2))
    end do
    !
  END SUBROUTINE generic_thomas_1d

  SUBROUTINE generic_thomas_2d(mp, n,A,B,C,r,u)

    USE cable_sli_numbers,       ONLY: one

    IMPLICIT NONE

    ! in/out
    INTEGER(i_d),                           INTENT(IN)  :: mp
    INTEGER(i_d),                           INTENT(IN)  :: n
    REAL(r_2), DIMENSION(1:mp,1:n,1:2,1:2), INTENT(IN)  :: A, B, C
    REAL(r_2), DIMENSION(1:mp,1:n,1:2),     INTENT(IN)  :: r
    REAL(r_2), DIMENSION(1:mp,1:n,1:2),     INTENT(OUT) :: u
    ! local
    REAL(r_2), DIMENSION(1:mp,1:n,1:2,1:2) :: G
    REAL(r_2), DIMENSION(1:mp,1:2,1:2)     :: bet
    REAL(r_2), DIMENSION(1:mp,1:2)         :: d
    REAL(r_2), DIMENSION(1:mp)             :: detbet, detbet1
    INTEGER(i_d)                           :: j
    REAL(r_2), DIMENSION(1:mp,1:2)         :: tmp1d
    REAL(r_2), DIMENSION(1:mp,1:2,1:2)     :: tmp2d

    ! j=1
    bet(1:mp,1:2,1:2) = B(1:mp,1,1:2,1:2)
    detbet(1:mp)      = bet(1:mp,1,1)*bet(1:mp,2,2) - bet(1:mp,1,2)*bet(1:mp,2,1)
    if (any(abs(detbet(1:mp)) < epsilon(detbet))) then
       write(2,*) 'generic_thomas_2d error1: det = 0'
       stop 'program terminated by generic_thomas_2d'
    endif
    detbet1(1:mp) = one/detbet(1:mp)
    d(1:mp,1:2)   = r(1:mp,1,1:2)
    u(1:mp,1,1)   = (bet(1:mp,2,2)*d(1:mp,1) - bet(1:mp,1,2)*d(1:mp,2))*detbet1(1:mp)
    u(1:mp,1,2)   = (bet(1:mp,1,1)*d(1:mp,2) - bet(1:mp,2,1)*d(1:mp,1))*detbet1(1:mp)
    ! j=2, n
    do j=2, n
       d(1:mp,1:2)       = C(1:mp,j-1,1:2,1)
       G(1:mp,j,1,1)     = (bet(1:mp,2,2)*d(1:mp,1) - bet(1:mp,1,2)*d(1:mp,2))*detbet1(1:mp)
       G(1:mp,j,2,1)     = (bet(1:mp,1,1)*d(1:mp,2) - bet(1:mp,2,1)*d(1:mp,1))*detbet1(1:mp)
       d(1:mp,1:2)       = C(1:mp,j-1,1:2,2)
       G(1:mp,j,1,2)     = (bet(1:mp,2,2)*d(1:mp,1) - bet(1:mp,1,2)*d(1:mp,2))*detbet1(1:mp)
       G(1:mp,j,2,2)     = (bet(1:mp,1,1)*d(1:mp,2) - bet(1:mp,2,1)*d(1:mp,1))*detbet1(1:mp)
       tmp2d(1:mp,1,1)   = A(1:mp,j,1,1)*G(1:mp,j,1,1) + A(1:mp,j,1,2)*G(1:mp,j,2,1)
       tmp2d(1:mp,1,2)   = A(1:mp,j,1,1)*G(1:mp,j,1,2) + A(1:mp,j,1,2)*G(1:mp,j,2,2)
       tmp2d(1:mp,2,1)   = A(1:mp,j,2,1)*G(1:mp,j,1,1) + A(1:mp,j,2,2)*G(1:mp,j,2,1)
       tmp2d(1:mp,2,2)   = A(1:mp,j,2,1)*G(1:mp,j,1,2) + A(1:mp,j,2,2)*G(1:mp,j,2,2)
       bet(1:mp,1:2,1:2) = B(1:mp,j,1:2,1:2) - tmp2d(1:mp,1:2,1:2)
       detbet(1:mp)      = bet(1:mp,1,1)*bet(1:mp,2,2) - bet(1:mp,1,2)*bet(1:mp,2,1)
       if (any(abs(detbet(1:mp)) < epsilon(detbet))) then
          write(2,*) 'generic_thomas_2d error2: det = 0 at j=', j
          stop 'program terminated by generic_thomas_2d'
       endif
       detbet1(1:mp) = one/detbet(1:mp)
       tmp1d(1:mp,1) = A(1:mp,j,1,1)*u(1:mp,j-1,1) + A(1:mp,j,1,2)*u(1:mp,j-1,2)
       tmp1d(1:mp,2) = A(1:mp,j,2,1)*u(1:mp,j-1,1) + A(1:mp,j,2,2)*u(1:mp,j-1,2)
       d(1:mp,1:2)   = r(1:mp,j,1:2) - tmp1d(1:mp,1:2)
       u(1:mp,j,1)   = (bet(1:mp,2,2)*d(1:mp,1) - bet(1:mp,1,2)*d(1:mp,2))*detbet1(1:mp)
       u(1:mp,j,2)   = (bet(1:mp,1,1)*d(1:mp,2) - bet(1:mp,2,1)*d(1:mp,1))*detbet1(1:mp)
    end do
    ! back substitution
    do j=n-1, 1, -1
       tmp1d(1:mp,1) = G(1:mp,j+1,1,1)*u(1:mp,j+1,1) + G(1:mp,j+1,1,2)*u(1:mp,j+1,2)
       tmp1d(1:mp,2) = G(1:mp,j+1,2,1)*u(1:mp,j+1,1) + G(1:mp,j+1,2,2)*u(1:mp,j+1,2)
       u(1:mp,j,1:2) = u(1:mp,j,1:2) - tmp1d(1:mp,1:2)
    end do
    !
  END SUBROUTINE generic_thomas_2d

  !**********************************************************************************************************************

  SUBROUTINE getfluxes_vp_1d(n, ns, dx, vtop, vbot, parin, var, hint, phimin, q, qya, qyb, qTa, qTb, &
       ql, qlya, qlyb, qv, qvT, qvh, qvya, qvyb, iflux, init, getq0, getqn, Tsoil, T0, nsat, nsatlast)

    IMPLICIT NONE

    INTEGER(i_d),                 INTENT(IN)    :: n
    INTEGER(i_d),                 INTENT(IN)    :: ns
    REAL(r_2),    DIMENSION(1:n), INTENT(IN)    :: dx
    TYPE(vars),                   INTENT(IN)    :: vtop
    TYPE(vars),                   INTENT(IN)    :: vbot
    TYPE(params), DIMENSION(1:n), INTENT(IN)    :: parin
    TYPE(vars),   DIMENSION(1:n), INTENT(IN)    :: var
    REAL(r_2),    DIMENSION(1:n), INTENT(INOUT) :: hint
    REAL(r_2),    DIMENSION(1:n), INTENT(INOUT) :: phimin
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: q
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: qya
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: qyb
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: qTa
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: qTb
    REAL(r_2),    DIMENSION(0:n), INTENT(OUT)   :: ql
    REAL(r_2),    DIMENSION(0:n), INTENT(OUT)   :: qlya
    REAL(r_2),    DIMENSION(0:n), INTENT(OUT)   :: qlyb
    REAL(r_2),    DIMENSION(0:n), INTENT(OUT)   :: qv
    REAL(r_2),    DIMENSION(0:n), INTENT(OUT)   :: qvT
    REAL(r_2),    DIMENSION(0:n), INTENT(OUT)   :: qvh
    REAL(r_2),    DIMENSION(0:n), INTENT(OUT)   :: qvya
    REAL(r_2),    DIMENSION(0:n), INTENT(OUT)   :: qvyb
    INTEGER(i_d),                 INTENT(IN)    :: iflux
    LOGICAL,                      INTENT(IN)    :: init
    LOGICAL,                      INTENT(IN)    :: getq0
    LOGICAL,                      INTENT(IN)    :: getqn
    REAL(r_2),    DIMENSION(1:n), INTENT(IN)    :: Tsoil
    REAL(r_2),                    INTENT(IN)    :: T0
    INTEGER(i_d),                 INTENT(IN)    :: nsat
    INTEGER(i_d),                 INTENT(IN)    :: nsatlast
    ! Gets fluxes q and partial derivs qya, qyb wrt S (if unsat) or phi (if sat).
    ! Fluxes at top and bottom of profile, and fluxes due to plant extraction of
    ! water are included.
    ! Definitions of arguments:
    ! k     - land point
    ! n     - no. of soil layers.
    ! jt(1:n)   - layer soil type nos.
    ! dx(1:n)   - layer thicknesses.
    ! dz(1:n-1)   - distances between layer centres.
    ! vtop    - water vars at soil surface.
    ! vbot    - water vars at bottom of profile.
    ! var(1:n)   - water vars at layer centres.
    ! hint(1:n)   - values of h at interfaces are stored sequentially in hint.
    ! phimin(1:n) - similarly for phi at hmin in layers above interfaces.
    ! q(0:n)   - fluxes; q(i), i=1,...,n-1 is flux from layer i to layer i+1.
    !    q(0) is surface flux and q(n) is flux at bottom of profile.
    ! qya(0:n)   - partial deriv of q(i), i=0,...,n, wrt the variable to be solved
    !    for (S, phi or h) at upper end of flow path.
    ! qyb(0:n)   - ditto for var at lower end.
    ! iflux    - if iflux/=1, get only fluxes involving sat layers.
    ! init    - true if hint and phimin to be initialised.
    ! getq0    - true if q(0) required.
    ! getqn    - true if q(n) required.
    LOGICAL               :: flag, limit
    INTEGER(i_d)          :: i, itmp, l
    REAL(r_2)             :: dphii1, dhi, h1, h2, hi, Khi1, Khi2, phii1, q2, qya2, qyb2, y, y1, y2
    REAL(r_2)             :: qTa2, qTb2
    TYPE(vars)            :: vi1, vi2
    REAL(r_2), DIMENSION(1:n-1) :: dz

    dz(:) = half*(dx(1:n-1)+dx(2:n))

    if ((iflux==1) .or. (var(1)%isat /= 0)) then ! get top flux if required
       if (getq0) then
          call flux(parin(1), ns==1, vtop, var(1), half*dx(1), q(0), qya(0), qyb(0), qTa(0), qTb(0))

          q(0)  = q(0)+(T0-Tsoil(1))*(var(1)%kE)/thousand/var(1)%lambdav/dx(1)*two

          qTa(0) = zero
          qTb(0) = -(var(1)%kE)/thousand/var(1)%lambdav/dx(1)*two
          qv(0)  = (vtop%phiv-var(1)%phiv)/dx(1)*two +(T0-Tsoil(1))*(var(1)%kE)/thousand/var(1)%lambdav/dx(1)*two
          ql(0)  = q(0) - qv(0)

          qvya(0) = zero

          if (vtop%isat==0) then
             qvyb(0) = vtop%phivS/dx(1)*two
          else
             qvyb(0) = zero
          end if

          qlya(0) = qya(0) - qvya(0)
          qlyb(0) = qyb(0) - qvyb(0)
       end if
    end if
    ! otherwise undefined
    qvh(0) = zero
    qvT(0) = zero

    ! get other fluxes
    l = 0
    do i=1, n-1
       if (iflux==1 .or. var(i)%isat/=0 .or. var(i+1)%isat/=0 .or. nsat/=nsatlast) then ! get flux
          if (parin(i)%isahorizon == parin(i+1)%isahorizon) then ! same soil type, no interface
             call flux(parin(i), .false., var(i), var(i+1), dz(i), q(i), qya(i), qyb(i), qTa(i), qTb(i))
          else ! interface
             l  = l+1
             if (init) then ! initialise
                call hyofh(hmin, parin(i), vi1%K, Khi1, phimin(l)) ! get phi at hmin
                h1 = var(i)%h
                h2 = var(i+1)%h
                y1 = var(i)%K*dx(i+1)
                y2 = var(i+1)%K*dx(i)
                ! equate fluxes (K constant) to get initial estimate of h at interface
                hint(l) = (y1*h1+y2*h2+half*gf*(var(i)%K-var(i+1)%K)*dx(i)*dx(i+1))/(y1+y2)
             end if
             hi   = hint(l)
             flag = .true.
             itmp = 0
             ! iterate to get hi at interface for equal fluxes using Newton's method
             ! get dphii1 at interface in upper layer, because of better linearity,
             ! then convert to dhi
             do while (flag)
                itmp = itmp+1
                if (itmp>100) then
                   !write(2,*) "getfluxes: too many iterations finding interface h"
                   !stop
                   call flux(parin(i), .false., var(i), var(i+1), dz(i), q(i), qya(i), qyb(i), qTa(i), qTb(i))
                   goto 111
                end if

                if (hi<parin(i)%he) then
                   vi1%isat = 0
                   call hyofh(hi, parin(i), vi1%K, Khi1, phii1)
                   vi1%KS = Khi1/vi1%K ! use dK/dphi, not dK/dS
                else
                   vi1%isat = 1
                   vi1%K    = parin(i)%Ke
                   phii1    = parin(i)%phie+(hi-parin(i)%he)*parin(i)%Ke
                   vi1%KS   = zero
                end if

                vi1%h    = hi
                vi1%phi  = phii1
                vi1%phiS = one ! use dphi/dphi not dphi/dS
                !MC! define phiT=0, KT=0 to be consitent with undefined version
                vi1%phiT = zero
                vi1%KT   = zero
                call flux(parin(i), .false., var(i), vi1, half*dx(i), q(i), qya(i), qyb(i), qTa(i), qTb(i))

                if (hi<parin(i+1)%he) then
                   vi2%isat = 0
                   call hyofh(hi, parin(i+1), vi2%K, Khi2, vi2%phi)
                   vi2%KS = Khi2/vi2%K ! dK/dphi
                else
                   vi2%isat = 1
                   vi2%K    = parin(i+1)%Ke
                   vi2%phi  = parin(i+1)%phie+(hi-parin(i+1)%he)*parin(i+1)%Ke
                end if

                vi2%h    = hi
                vi2%phiS = one ! dphi/dphi
                !MC! define phiT=0, KT=0 to be consitent with undefined version
                vi2%phiT = zero
                vi2%KT   = zero
                call flux(parin(i+1), .false., vi2, var(i+1), half*dx(i+1), q2, qya2, qyb2, qTa2, qTb2)
                qya2   = qya2*vi2%K/vi1%K ! partial deriv wrt phii1
                ! adjust for equal fluxes
                dphii1 = -(q(i)-q2)/(qyb(i)-qya2)
                limit  = .false.
                if (phii1+dphii1<=phimin(l)) then ! out of range
                   limit  = .true.
                   dphii1 = -half*(phii1-phimin(l))
                end if
                phii1 = phii1+dphii1
                dhi   = dphii1/(vi1%K+half*vi1%KS*dphii1) ! 2nd order Pade approx
                if (-vi1%KS*dphii1 > 1.5_r_2*vi1%K) then ! use 1st order approx for dhi
                   dhi = dphii1/vi1%K
                end if
                hi = hi+dhi

                ! check for convergence - dphi/(mean phi)<=dpmaxr
                if (.not.(limit .or. abs(dphii1/(phii1-half*dphii1))>dpmaxr)) then
                   flag = .false.
                end if
             end do ! while flag

             q(i)    = q(i) + qyb(i)*dphii1
             hint(l) = hi
             ! adjust derivs
             y      = one/(qya2-qyb(i))
             qya(i) = qya(i)*qya2*y
             qyb(i) = -qyb2*qyb(i)*y
          end if
       end if
   
111    ql(i)  = q(i)
       qTa(i) = qTa(i)+(var(i)%kE+var(i+1)%kE)/thousand/var(1)%lambdav/two/dz(i)
       qTb(i) = qTb(i)-(var(i)%kE+var(i+1)%kE)/thousand/var(1)%lambdav/two/dz(i)
       qvT(i) = (Tsoil(i)-Tsoil(i+1))*(var(i)%kE+var(i+1)%kE)/thousand/var(1)%lambdav/two/dz(i)
       !!! MC !!! The full description (next two lines) gave problems before ->  third line
       !!! MC !!!     Try again original
       !!! VH !!! Do both formulations come to the same?? I found that I couldn't reproduce 
       !!! VH !!! Barnes and Allison semi-analytic sol'n with original
       !MC! this should be re-checked
       qvh(i) = ((((Tsoil(i)+Tzero)/Tzero)**1.88+((Tsoil(i+1)+Tzero)/Tzero)**1.88)/two) &
            * ((var(i)%cvsat+var(i+1)%cvsat)/two)*(var(i)%phiv-var(i+1)%phiv)/dz(i)
       ! qvh(i) = ((var(i)%Dv+var(i+1)%Dv)/two)* ((var(i)%cvsat+var(i+1)%cvsat)/two)*(var(i)%rh-var(i+1)%rh)/dz(i)
       qv(i)  = qvh(i) + qvT(i) ! whole vapour flux has one part from humidity (qvh) and one part from temp diff (qvT)
       q(i)   = qv(i) + ql(i)

       if (var(i)%isat==0) then
          qvya(i) = var(i)%phivS/dz(i) *((((Tsoil(i)+Tzero)/Tzero)**1.88_r_2+ &
               ((Tsoil(i+1)+Tzero)/Tzero)**1.88_r_2)/two) &
               * ((var(i)%cvsat+var(i+1)%cvsat)/two)
       else
          qvya(i) = zero
       end if

       if (var(i)%isat==0) then
          qvyb(i) = -var(i+1)%phivS/dz(i) *((((Tsoil(i)+Tzero)/Tzero)**1.88_r_2+ &
               ((Tsoil(i+1)+Tzero)/Tzero)**1.88_r_2)/two) &
               * ((var(i)%cvsat+var(i+1)%cvsat)/two)
       else
          qvyb(i) = zero
       end if

       qlya(i) = qya(i)
       qlyb(i) = qyb(i)
       qya(i)  = qya(i) + qvya(i)
       qyb(i)  = qyb(i) + qvyb(i)
    end do

    if (iflux==1 .or. var(n)%isat/=0) then ! get bottom flux if required
       if (getqn) then
          !MC! what is i here? It is after the do i=1, n-1 loop. -> take n
          !call flux(jt(n), i, ns, v, vbot, half*dx(n), q(n), qya(n), qyb(n), qTa(n), qTb(n))
          call flux(parin(n), .false., var(n), vbot, half*dx(n), q(n), qya(n), qyb(n), qTa(n), qTb(n))
          qvya(n) = zero
          qvyb(n) = zero
          qlya(n) = qya(n)
          qlyb(n) = zero
       else
          qvya(n) = zero
          qvyb(n) = zero
          qlya(n) = qya(n)
          qlyb(n) = zero
       end if
    else
       qvya(n) = zero
       qvyb(n) = zero
       qlya(n) = qya(n)
       qlyb(n) = zero
    end if
    ! otherwise undefined
    ql(n)  = q(n)
    qv(n)  = zero
    qvh(n) = zero
    qvT(n) = zero

    do i=1, n-1
       if (var(i)%Dv == zero .or. var(i+1)%Dv == zero) then
          q(i)    = q(i) - qv(i)
          qya(i)  = qya(i) - qvya(i)
          qyb(i)  = qyb(i) - qvyb(i)
          qv(i)   = zero
          !qTa(i)  = zero
          !qTb(i)  = zero
          qvya(i) = zero
          qvyb(i) = zero
       endif
    enddo

  END SUBROUTINE getfluxes_vp_1d

  SUBROUTINE getfluxes_vp_2d(ns, dx, vtop, vbot, parin, var, hint, phimin, i_q, i_qya, i_qyb, i_qTa, i_qTb, &
       i_ql, i_qlya, i_qlyb, i_qv, i_qvT, i_qvh, i_qvya, i_qvyb, iflux, init, getq0, getqn, Tsoil, T0, nsat, nsatlast)

    IMPLICIT NONE

    INTEGER(i_d), DIMENSION(:),   INTENT(IN)    :: ns
    REAL(r_2),    DIMENSION(:,:), INTENT(IN)    :: dx      ! 1:n
    TYPE(vars),   DIMENSION(:),   INTENT(IN)    :: vtop
    TYPE(vars),   DIMENSION(:),   INTENT(IN)    :: vbot
    TYPE(params), DIMENSION(:,:), INTENT(IN)    :: parin   ! 1:n
    TYPE(vars),   DIMENSION(:,:), INTENT(IN)    :: var     ! 1:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: hint    ! 1:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: phimin  ! 1:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_q       ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qya     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qyb     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qTa     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qTb     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_ql      ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qlya    ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qlyb    ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qv      ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qvT     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qvh     ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qvya    ! 0:n
    REAL(r_2),    DIMENSION(:,:), INTENT(OUT)   :: i_qvyb    ! 0:n
    INTEGER(i_d), DIMENSION(:),   INTENT(IN)    :: iflux
    LOGICAL,      DIMENSION(:),   INTENT(IN)    :: init
    LOGICAL,      DIMENSION(:),   INTENT(IN)    :: getq0
    LOGICAL,      DIMENSION(:),   INTENT(IN)    :: getqn
    REAL(r_2),    DIMENSION(:,:), INTENT(IN)    :: Tsoil   ! 1:n
    REAL(r_2),    DIMENSION(:),   INTENT(IN)    :: T0
    INTEGER(i_d), DIMENSION(:),   INTENT(IN)    :: nsat
    INTEGER(i_d), DIMENSION(:),   INTENT(IN)    :: nsatlast
    ! Gets fluxes q and partial derivs qya, qyb wrt S (if unsat) or phi (if sat).
    ! Fluxes at top and bottom of profile, and fluxes due to plant extraction of
    ! water are included.
    ! Definitions of arguments:
    ! k     - land point
    ! n     - no. of soil layers.
    ! jt(1:n)   - layer soil type nos.
    ! dx(1:n)   - layer thicknesses.
    ! dz(1:n-1)   - distances between layer centres.
    ! vtop    - water vars at soil surface.
    ! vbot    - water vars at bottom of profile.
    ! var(1:n)   - water vars at layer centres.
    ! hint(1:n)   - values of h at interfaces are stored sequentially in hint.
    ! phimin(1:n) - similarly for phi at hmin in layers above interfaces.
    ! q(0:n)   - fluxes; q(i), i=1,...,n-1 is flux from layer i to layer i+1.
    !    q(0) is surface flux and q(n) is flux at bottom of profile.
    ! qya(0:n)   - partial deriv of q(i), i=0,...,n, wrt the variable to be solved
    !    for (S, phi or h) at upper end of flow path.
    ! qyb(0:n)   - ditto for var at lower end.
    ! iflux    - if iflux/=1, get only fluxes involving sat layers.
    ! init    - true if hint and phimin to be initialised.
    ! getq0    - true if q(0) required.
    ! getqn    - true if q(n) required.
    LOGICAL,      DIMENSION(1:size(dx,1))                :: flag, limit, ltmp, l1, l2, l3
    REAL(r_2),    DIMENSION(1:size(dx,1))                :: dphii1, dhi, h1, h2, hi, Khi1, Khi2, phii1, q2, qya2, qyb2, y, y1, y2
    REAL(r_2),    DIMENSION(1:size(dx,1))                :: qTa2, qTb2, htmp
    REAL(r_2),    DIMENSION(1:size(dx,1))                :: ztmp1, ztmp2, ztmp3, ztmp4, ztmp5
    TYPE(vars),   DIMENSION(1:size(dx,1))                :: vi1, vi2
    REAL(r_2),    DIMENSION(1:size(dx,1),1:size(dx,2)-1) :: dz
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: q
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qya
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qyb
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qTa
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qTb
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: ql
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qlya
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qlyb
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qv
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qvT
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qvh
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qvya
    REAL(r_2),    DIMENSION(1:size(dx,1),0:size(dx,2))   :: qvyb
    TYPE(vars)    :: vtmp
    INTEGER(i_d)  :: i, n, mp, itmp

    mp = size(dx,1)
    n  = size(dx,2)
    dz(:,1:n-1) = half*(dx(:,1:n-1)+dx(:,2:n))
    q(:,0:n)    = i_q(:,1:n+1)
    qya(:,0:n)  = i_qya(:,1:n+1)
    qyb(:,0:n)  = i_qyb(:,1:n+1)
    qTa(:,0:n)  = i_qTa(:,1:n+1)
    qTb(:,0:n)  = i_qTb(:,1:n+1)
    ql(:,0:n)   = zero
    qlya(:,0:n) = zero
    qlyb(:,0:n) = zero
    qv(:,0:n)   = zero
    qvT(:,0:n)  = zero
    qvh(:,0:n)  = zero
    qvya(:,0:n) = zero
    qvyb(:,0:n) = zero

    vtmp = zerovars()
    vi1  = spread(vtmp,1,mp)
    vi2  = spread(vtmp,1,mp)
    ztmp1(:) = zero
    ztmp2(:) = zero
    ztmp3(:) = zero
    ztmp4(:) = zero
    ztmp5(:) = zero

    ltmp(:) = ns(:)==1 ! sets logical ltmp to .true. where ns =1 
    l1(:)   = ((iflux(:)==1) .or. (var(:,1)%isat /= 0)) .and. getq0(:)
    if (any(l1(:))) &
         call flux(parin(:,1), ltmp(:), vtop(:), var(:,1), half*dx(:,1), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
    where (l1(:)) ! get top flux if required
       q(:,0)   = ztmp1(:)
       qya(:,0) = ztmp2(:)
       qyb(:,0) = ztmp3(:)
       qTa(:,0) = ztmp4(:)
       qTb(:,0) = ztmp5(:)

       q(:,0)    = q(:,0)+(T0(:)-Tsoil(:,1))*(var(:,1)%kE)/thousand/var(:,1)%lambdav/dx(:,1)*two
       qTa(:,0)  = zero
       qTb(:,0)  = -(var(:,1)%kE)/thousand/var(:,1)%lambdav/dx(:,1)*two
       qv(:,0)   = (vtop(:)%phiv-var(:,1)%phiv)/dx(:,1)*two +(T0(:)-Tsoil(:,1))*(var(:,1)%kE)/thousand/var(:,1)%lambdav/dx(:,1)*two
       ql(:,0)   = q(:,0) - qv(:,0)
       qvya(:,0) = zero
       qvyb(:,0) = zero ! isat==0 below
       qlya(:,0) = qya(:,0) - qvya(:,0)
       qlyb(:,0) = qyb(:,0) - qvyb(:,0)
    endwhere
    where (l1(:) .and. (vtop(:)%isat==0))
       qvyb(:,0) = vtop(:)%phivS/dx(:,1)*two
       qlyb(:,0) = qyb(:,0) - qvyb(:,0)
    endwhere
    ! otherwise undefined
    qvh(:,0) = zero
    qvT(:,0) = zero

    ! get other fluxes
    ltmp(:) = .false.
    do i=1, n-1
       l1(:) = (iflux(:)==1 .or. var(:,i)%isat/=0 .or. var(:,i+1)%isat/=0 .or. nsat(:)/=nsatlast(:))
       if (any(l1(:) .and. (parin(:,i)%isahorizon==parin(:,i+1)%isahorizon))) &
            call flux(parin(:,i), ltmp(:), var(:,i), var(:,i+1), dz(:,i), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
       where (l1(:) .and. (parin(:,i)%isahorizon==parin(:,i+1)%isahorizon)) ! same soil type, no interface
          q(:,i)   = ztmp1(:)
          qya(:,i) = ztmp2(:)
          qyb(:,i) = ztmp3(:)
          qTa(:,i) = ztmp4(:)
          qTb(:,i) = ztmp5(:)
       endwhere
       l2(:) = l1(:) .and. (parin(:,i)%isahorizon /= parin(:,i+1)%isahorizon) ! interface
       htmp(:) = hmin
       if (any(l2(:) .and. init(:))) &
            call hyofh(htmp(:), parin(:,i), ztmp1(:), ztmp2(:), ztmp3(:))
       where (l2(:) .and. init(:)) ! interface & get phi at hmin
          vi1(:)%K    = ztmp1(:)
          Khi1(:)     = ztmp2(:)
          phimin(:,i) = ztmp3(:)
          h1(:)       = var(:,i)%h
          h2(:)       = var(:,i+1)%h
          y1(:)       = var(:,i)%K*dx(:,i+1)
          y2(:)       = var(:,i+1)%K*dx(:,i)
          ! equate fluxes (K constant) to get initial estimate of h at interface
          hint(:,i) = (y1(:)*h1(:) + y2(:)*h2(:) + half*gf*(var(:,i)%K-var(:,i+1)%K)*dx(:,i)*dx(:,i+1)) / (y1(:)+y2(:))
       endwhere
       !where ((.not. l2(:)) .and. init(:)) hint(:,i) = zero
       hi(:)   = hint(:,i)
       flag(:) = .true.
       ! iterate to get hi at interface for equal fluxes using Newton's method
       ! get dphii1 at interface in upper layer, because of better linearity,
       ! then convert to dhi
       if (any(l2(:))) then
          l3(:)    = l2(:)
          limit(:) = .false.
          do itmp=1, 100
             if (any(l3(:) .and. (hi(:)<parin(:,i)%he))) &
                  call hyofh(hi(:), parin(:,i), ztmp1(:), ztmp2(:), ztmp3(:))
             where (l3(:) .and. (hi(:)<parin(:,i)%he))
                vi1(:)%isat = 0
                vi1(:)%K    = ztmp1(:)
                Khi1(:)     = ztmp2(:)
                phii1(:)    = ztmp3(:)
                vi1(:)%KS   = Khi1(:)/vi1(:)%K ! use dK/dphi, not dK/dS
             endwhere
             where (l3(:) .and. (hi(:)>=parin(:,i)%he))
                vi1(:)%isat = 1
                vi1(:)%K    = parin(:,i)%Ke
                phii1(:)    = parin(:,i)%phie+(hi(:)-parin(:,i)%he)*parin(:,i)%Ke
                vi1(:)%KS   = zero
             endwhere

             where (l3(:))
                vi1(:)%h    = hi(:)
                vi1(:)%phi  = phii1(:)
                vi1(:)%phiS = one ! use dphi/dphi not dphi/dS
                !MC! define phiT=0, KT=0 to be consitent with undefined version
                vi1(:)%phiT = zero
                vi1(:)%KT   = zero
             endwhere

             if (any(l3(:))) &
                  call flux(parin(:,i), ltmp(:), var(:,i), vi1(:), half*dx(:,i), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
             where (l3(:))
                q(:,i)      = ztmp1(:)
                qya(:,i)    = ztmp2(:)
                qyb(:,i)    = ztmp3(:)
                qTa(:,i)    = ztmp4(:)
                qTb(:,i)    = ztmp5(:)
             endwhere

             if (any(l3(:) .and. (hi(:)<parin(:,i+1)%he))) &
                  call hyofh(hi(:), parin(:,i+1), ztmp1(:), ztmp2(:), ztmp3(:))
             where (l3(:) .and. (hi(:)<parin(:,i+1)%he))
                vi2(:)%K    = ztmp1(:)
                Khi2(:)     = ztmp2(:)
                vi2(:)%phi  = ztmp3(:)
                vi2(:)%isat = 0
                vi2(:)%KS   = Khi2(:)/vi2(:)%K ! dK/dphi
             endwhere
             where (l3(:) .and. (hi(:)>=parin(:,i+1)%he))
                vi2(:)%isat = 1
                vi2(:)%K    = parin(:,i+1)%Ke
                vi2(:)%phi  = parin(:,i+1)%phie+(hi(:)-parin(:,i+1)%he)*parin(:,i+1)%Ke
             endwhere

             where (l3(:))
                vi2(:)%h    = hi(:)
                vi2(:)%phiS = one ! dphi/dphi
                !MC! define phiT=0, KT=0 to be consitent with undefined version
                vi2(:)%phiT = zero
                vi2(:)%KT   = zero
             endwhere
             if (any(l3(:))) &
                  call flux(parin(:,i+1), ltmp(:), vi2(:), var(:,i+1), half*dx(:,i+1), &
                  ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
             where (l3(:))
                q2(:)     = ztmp1(:)
                qya2(:)   = ztmp2(:)
                qyb2(:)   = ztmp3(:)
                qTa2(:)   = ztmp4(:)
                qTb2(:)   = ztmp5(:)
                qya2(:)   = qya2(:)*vi2(:)%K/vi1(:)%K ! partial deriv wrt phii1
                ! adjust for equal fluxes
                dphii1(:) = -(q(:,i)-q2(:))/(qyb(:,i)-qya2(:))
                limit(:)  = .false.
             endwhere
             where (l3(:) .and. (phii1(:)+dphii1(:)<=phimin(:,i))) ! out of range
                limit(:)  = .true.
                dphii1(:) = -half*(phii1(:)-phimin(:,i))
             endwhere
             where (l3(:))
                phii1(:) = phii1(:)+dphii1(:)
                dhi(:)   = dphii1(:)/(vi1(:)%K+half*vi1(:)%KS*dphii1(:)) ! 2nd order Pade approx
             endwhere
             where (l3(:) .and. (-vi1%KS*dphii1 > 1.5_r_2*vi1%K)) ! use 1st order approx for dhi
                dhi(:) = dphii1(:)/vi1(:)%K
             endwhere
             where (l3(:))
                hi(:) = hi(:)+dhi(:)
             endwhere

             ! check for convergence - dphi/(mean phi)<=dpmaxr
             where (l3(:) .and. &
                  .not. (limit(:) .or. (abs(dphii1(:)/(phii1(:)-half*dphii1(:)))>dpmaxr))) l3(:) = .false.
             if (.not. any(l3(:))) exit
          end do ! do itmp=1, 100
          if (itmp>=100) then
             !write(2,*) "getfluxes: too many iterations finding interface h"
             !stop
             if (any(l2(:) .and. l3(:))) &
                call flux(parin(:,i), ltmp(:), var(:,i), var(:,i+1), dz(:,i), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
             where (l2(:) .and. l3(:))
                q(:,i)   = ztmp1(:)
                qya(:,i) = ztmp2(:)
                qyb(:,i) = ztmp3(:)
                qTa(:,i) = ztmp4(:)
                qTb(:,i) = ztmp5(:)
             endwhere
          else
             where (l2(:) .and. (.not. l3(:)))
                q(:,i)    = q(:,i) + qyb(:,i)*dphii1(:)
                hint(:,i) = hi(:)
                ! adjust derivs
                y(:)      = one/(qya2(:)-qyb(:,i))
                qya(:,i) = qya(:,i)*qya2(:)*y(:)
                qyb(:,i) = -qyb2(:)*qyb(:,i)*y(:)
             endwhere
          end if
       end if

       ql(:,i)  = q(:,i)
       qTa(:,i) = qTa(:,i)+(var(:,i)%kE+var(:,i+1)%kE)/thousand/var(:,1)%lambdav/two/dz(:,i)
       qTb(:,i) = qTb(:,i)-(var(:,i)%kE+var(:,i+1)%kE)/thousand/var(:,1)%lambdav/two/dz(:,i)
       qvT(:,i) = (Tsoil(:,i)-Tsoil(:,i+1))*(var(:,i)%kE+var(:,i+1)%kE)/thousand/var(:,1)%lambdav/two/dz(:,i)
       !!! MC !!! The full description (next two lines) gave problems before ->  third line
       !!! MC !!!     Try again original
       !!! VH !!! Do both formulations come to the same?? I found that I couldn't reproduce 
       !!! VH !!! Barnes and Allison semi-analytic sol'n with original
       !MC! this should be re-checked
       qvh(:,i) = ((((Tsoil(:,i)+Tzero)/Tzero)**1.88+((Tsoil(:,i+1)+Tzero)/Tzero)**1.88)/two) &
            * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)*(var(:,i)%phiv-var(:,i+1)%phiv)/dz(:,i)
       ! qvh(:,i) = ((var(:,i)%Dv+var(:,i+1)%Dv)/two)* ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)*(var(:,i)%rh-var(:,i+1)%rh)/dz(:,i)
       qv(:,i)  = qvh(:,i) + qvT(:,i) ! whole vapour flux has one part from humidity (qvh) and one part from temp diff (qvT)
       q(:,i)   = qv(:,i) + ql(:,i)

       where (var(:,i)%isat==0)
          qvya(:,i) = var(:,i)%phivS/dz(:,i) *((((Tsoil(:,i)+Tzero)/Tzero)**1.88_r_2+ &
               ((Tsoil(:,i+1)+Tzero)/Tzero)**1.88_r_2)/two) &
               * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
       elsewhere
          qvya(:,i) = zero
       endwhere

       where (var(:,i)%isat==0)
          qvyb(:,i) = -var(:,i+1)%phivS/dz(:,i) *((((Tsoil(:,i)+Tzero)/Tzero)**1.88_r_2+ &
               ((Tsoil(:,i+1)+Tzero)/Tzero)**1.88_r_2)/two) &
               * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
       elsewhere
          qvyb(:,i) = zero
       endwhere

       qlya(:,i) = qya(:,i)
       qlyb(:,i) = qyb(:,i)
       qya(:,i)  = qya(:,i) + qvya(:,i)
       qyb(:,i)  = qyb(:,i) + qvyb(:,i)
    end do

    l1(:) = (iflux(:)==1) .or. (var(:,n)%isat/=0)
    !MC! what is i here? It is after the do i=1, n-1 loop. -> take n
    !call flux(jt(n), i, ns, v, vbot, half*dx(n), q(n), qya(n), qyb(n), qTa(n), qTb(n))
    if (any(l1(:) .and. getqn(:))) &  ! get bottom flux if required
         call flux(parin(:,n), ltmp(:), var(:,n), vbot(:), half*dx(:,n), ztmp1(:), ztmp2(:), ztmp3(:), ztmp4(:), ztmp5(:))
    where (l1(:) .and. getqn(:))
       q(:,n)    = ztmp1(:)
       qya(:,n)  = ztmp2(:)
       qyb(:,n)  = ztmp3(:)
       qTa(:,n)  = ztmp4(:)
       qTb(:,n)  = ztmp5(:)
       qvya(:,n) = zero
       qvyb(:,n) = zero
       qlya(:,n) = qya(:,n)
       qlyb(:,n) = zero
    endwhere

    where (l1(:) .and. (.not. getqn(:)))
       qvya(:,n) = zero
       qvyb(:,n) = zero
       qlya(:,n) = qya(:,n)
       qlyb(:,n) = zero
    endwhere
    where (.not. l1(:))
       qvya(:,n) = zero
       qvyb(:,n) = zero
       qlya(:,n) = qya(:,n)
       qlyb(:,n) = zero
    endwhere
    ! otherwise undefined
    ql(:,n)  = q(:,n)
    qv(:,n)  = zero
    qvh(:,n) = zero
    qvT(:,n) = zero

    do i=1, n-1
       where (var(:,i)%Dv == zero .or. var(:,i+1)%Dv == zero)
          q(:,i)    = q(:,i) - qv(:,i)
          qya(:,i)  = qya(:,i) - qvya(:,i)
          qyb(:,i)  = qyb(:,i) - qvyb(:,i)
          qv(:,i)   = zero
          !qTa(:,i)  = zero
          !qTb(:,i)  = zero
          qvya(:,i) = zero
          qvyb(:,i) = zero
       endwhere
    enddo

    i_q(:,1:n+1)    = q(:,0:n)    
    i_qya(:,1:n+1)  = qya(:,0:n)  
    i_qyb(:,1:n+1)  = qyb(:,0:n)  
    i_qTa(:,1:n+1)  = qTa(:,0:n)  
    i_qTb(:,1:n+1)  = qTb(:,0:n)  
    i_ql(:,1:n+1)   = ql(:,0:n)   
    i_qlya(:,1:n+1) = qlya(:,0:n) 
    i_qlyb(:,1:n+1) = qlyb(:,0:n) 
    i_qv(:,1:n+1)   = qv(:,0:n)   
    i_qvT(:,1:n+1)  = qvT(:,0:n)  
    i_qvh(:,1:n+1)  = qvh(:,0:n)  
    i_qvya(:,1:n+1) = qvya(:,0:n) 
    i_qvyb(:,1:n+1) = qvyb(:,0:n) 

  END SUBROUTINE getfluxes_vp_2d

  !**********************************************************************************************************************

  SUBROUTINE getheatfluxes_1d(n, ns, h0, dx, dxL, qh, qhya, qhyb, qhTa, qhTb, var, vlit, T, TL, T0, litter, &
								q, qya, qyb, qTa, qTb,qadv, qadvya, qadvyb, qadvTa, qadvTb,advection)
    ! modified 25/05/10 to include contribution to heat flux from liquid water flux in the presence of ice
    IMPLICIT NONE

    INTEGER(i_d),               INTENT(IN)    :: n
    INTEGER(i_d),               INTENT(IN)    :: ns
    REAL(r_2),                  INTENT(IN)    :: h0
    REAL(r_2),  DIMENSION(1:n), INTENT(IN)    :: dx
    REAL(r_2),                  INTENT(IN)    :: dxL
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qh, q, qadv
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qhya, qya, qadvya
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qhyb, qyb, qadvyb
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qhTa, qTa, qadvTa
    REAL(r_2),  DIMENSION(0:n), INTENT(INOUT) :: qhTb, qTb, qadvTb
    TYPE(vars), DIMENSION(1:n), INTENT(IN)    :: var
    TYPE(vars),                 INTENT(IN)    :: vlit
    REAL(r_2),  DIMENSION(1:n), INTENT(IN)    :: T
    REAL(r_2),                  INTENT(IN)    :: TL
    REAL(r_2),                  INTENT(IN)    :: T0
    LOGICAL,                    INTENT(IN)    :: litter
	INTEGER(i_d),               INTENT(IN)    :: advection
    ! Gets heat fluxes qh and partial derivs qhya, qhyb wrt T and S (if unsat) or phi (if sat).

    REAL(r_2),  PARAMETER :: kw=0.57
    INTEGER(i_d)          :: i
    REAL(r_2)             :: rdz, w, keff
    REAL(r_2), DIMENSION(1:n-1) :: dz

    dz(:) = half*(dx(1:n-1)+dx(2:n))


    do i=1, n-1
       rdz = one/dz(i)
	   keff = 2_r_2*(var(i)%kth*var(i+1)%kth)/(var(i)%kth*dx(i)+var(i+1)%kth*dx(i+1))
       qh(i) = keff*(T(i)-T(i+1)) +(var(i)%phiv-var(i+1)%phiv)*var(i)%lambdav*thousand*rdz &
            * ((((T(i)+Tzero) /Tzero)**1.88_r_2+((T(i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
            *((var(i)%cvsat+var(i+1)%cvsat)/two)
       if (var(i)%isat==0) then
          qhya(i) = rdz*var(i)%lambdav*thousand*var(i)%phivS*((((T(i)+Tzero) /Tzero)**1.88_r_2+ &
               ((T(i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
               * ((var(i)%cvsat+var(i+1)%cvsat)/two)
       else
          qhya(i) = zero
       end if
       if (var(i+1)%isat==0) then
          qhyb(i) = -rdz*var(i)%lambdav*thousand*var(i+1)%phivS*((((T(i)+Tzero) /Tzero)**1.88_r_2+ &
               ((T(i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
               * ((var(i)%cvsat+var(i+1)%cvsat)/two)
       else
          qhyb(i) = zero
       end if
       qhTa(i) = keff
       qhTb(i) = -keff

	   	   ! add advective terms
	   if (advection==1) then
		   if (q(i).gt.zero) then
			w = (var(i)%kth/dx(i))/(var(i)%kth/dx(i)+var(i+1)%kth/dx(i+1))
		   else
			w = (var(i)%kth/dx(i))/(var(i)%kth/dx(i)+var(i+1)%kth/dx(i+1))
		   endif
		   qadv(i) = rhow*cswat*q(i)*(w*(T(i)+Tzero)+(1.-w)*(T(i+1)+Tzero))
		   qadvya(i) =  rhow*cswat*qya(i)*(w*(T(i)+Tzero)+(1.-w)*(T(i+1)+Tzero))
		   qadvyb(i) =  rhow*cswat*qyb(i)*(w*(T(i)+Tzero)+(1.-w)*(T(i+1)+Tzero))
		   qadvTa(i) =  rhow*cswat*q(i)*w
		   qadvTb(i) =  rhow*cswat*q(i)*(1.-w)

		   qh(i) = qh(i) + qadv(i)
		   qhya(i) = qhya(i) + qadvya(i)
		   qhyb(i) = qhyb(i) + qadvyb(i)
		   qhTa(i) = qhTa(i) + qadvTa(i)
		   qhTb(i) = qhTb(i) + qadvTb(i)
	   endif

    enddo

    qh(n)   = zero
    qhya(n) = zero
    qhyb(n) = zero
    qhTa(n) = zero
    qhTb(n) = zero

  END SUBROUTINE getheatfluxes_1d

  SUBROUTINE getheatfluxes_2d(ns, h0, dx, dxL, i_qh, i_qhya, i_qhyb, i_qhTa, i_qhTb, var, vlit, T, TL, T0, &
                              litter, i_q,i_qya,i_qyb,i_qTa,i_qTb,&
							  i_qadv,i_qadvya, i_qadvyb, i_qadvTa, i_qadvTb, advection)
    ! modified 25/05/10 to include contribution to heat flux from liquid water flux in the presence of ice
    IMPLICIT NONE

    INTEGER(i_d), DIMENSION(:),   INTENT(IN)    :: ns
    REAL(r_2),    DIMENSION(:),   INTENT(IN)    :: h0
    REAL(r_2),    DIMENSION(:,:), INTENT(IN)    :: dx      ! :,1:n
    REAL(r_2),    DIMENSION(:),   INTENT(IN)    :: dxL
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qh    ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qhya  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qhyb  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qhTa  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qhTb  ! :,0:n => :,1:n+1
	REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_q    ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qya  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qyb  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qTa  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qTb  ! :,0:n => :,1:n+1
	REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadv    ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadvya  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadvyb  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadvTa  ! :,0:n => :,1:n+1
    REAL(r_2),    DIMENSION(:,:), INTENT(INOUT) :: i_qadvTb  ! :,0:n => :,1:n+1
    TYPE(vars),   DIMENSION(:,:), INTENT(IN)    :: var
    TYPE(vars),   DIMENSION(:),   INTENT(IN)    :: vlit
    REAL(r_2),    DIMENSION(:,:), INTENT(IN)    :: T       ! :,1:n
    REAL(r_2),    DIMENSION(:),   INTENT(IN)    :: TL
    REAL(r_2),    DIMENSION(:),   INTENT(IN)    :: T0
    LOGICAL,                      INTENT(IN)    :: litter
	INTEGER(i_d),   INTENT(IN)    :: advection
    ! Gets heat fluxes qh and partial derivs qhya, qhyb wrt T and S (if unsat) or phi (if sat).

    REAL(r_2),  PARAMETER :: kw=0.57
    INTEGER(i_d)          :: i, n
    REAL(r_2),  DIMENSION(1:size(dx,1))                :: rdz
	REAL(r_2),  DIMENSION(1:size(dx,1))                :: keff
	REAL(r_2),  DIMENSION(1:size(dx,1))                :: w
    REAL(r_2),  DIMENSION(1:size(dx,1),1:size(dx,2)-1) :: dz
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qh   ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qhya ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qhyb ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qhTa ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qhTb ! :,0:n
	REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: q   ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qya ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qyb ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qTa ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qTb ! :,0:n
	REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadv   ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadvya ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadvyb ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadvTa ! :,0:n
    REAL(r_2),  DIMENSION(1:size(dx,1),0:size(dx,2))   :: qadvTb ! :,0:n

    n = size(dx,2)
    dz(:,1:n-1) = half*(dx(:,1:n-1)+dx(:,2:n))
    qh(:,0:n)   = i_qh(:,1:n+1)
    qhya(:,0:n) = i_qhya(:,1:n+1)
    qhyb(:,0:n) = i_qhyb(:,1:n+1)
    qhTa(:,0:n) = i_qhTa(:,1:n+1)
    qhTb(:,0:n) = i_qhTb(:,1:n+1)

	q(:,0:n)   = i_q(:,1:n+1)
    qya(:,0:n) = i_qya(:,1:n+1)
    qyb(:,0:n) = i_qyb(:,1:n+1)
    qTa(:,0:n) = i_qTa(:,1:n+1)
    qTb(:,0:n) = i_qTb(:,1:n+1)

    qadv(:,0:n)   = i_qadv(:,1:n+1)
    qadvya(:,0:n) = i_qadvya(:,1:n+1)
    qadvyb(:,0:n) = i_qadvyb(:,1:n+1)
    qadvTa(:,0:n) = i_qadvTa(:,1:n+1)
    qadvTb(:,0:n) = i_qadvTb(:,1:n+1)

   
    do i=1, n-1
       rdz =  one/dz(:,i)
	   keff = 2_r_2*(var(:,i)%kth*var(:,i+1)%kth)/(var(:,i)%kth*dx(:,i)+var(:,i+1)%kth*dx(:,i+1))
       qh(:,i) = keff*(T(:,i)-T(:,i+1)) &
            + (var(:,i)%phiv-var(:,i+1)%phiv)*var(:,i)%lambdav*thousand*rdz(:) &
            * ((((T(:,i)+Tzero) /Tzero)**1.88_r_2+((T(:,i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
            *((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
       where (var(:,i)%isat == 0)
          qhya(:,i) = rdz(:)*var(:,i)%lambdav*thousand*var(:,i)%phivS*((((T(:,i)+Tzero) /Tzero)**1.88_r_2+ &
               ((T(:,i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
               * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
       elsewhere
          qhya(:,i) = zero
       endwhere
       where (var(:,i+1)%isat == 0)
          qhyb(:,i) = -rdz(:)*var(:,i)%lambdav*thousand*var(:,i+1)%phivS*((((T(:,i)+Tzero) /Tzero)**1.88_r_2+ &
               ((T(:,i+1)+Tzero) /Tzero)**1.88_r_2)/two) &
               * ((var(:,i)%cvsat+var(:,i+1)%cvsat)/two)
       elsewhere
          qhyb(:,i) = zero
       endwhere
       qhTa(:,i) = keff
       qhTb(:,i) = -keff

	   	   	   ! add advective terms
		if (advection==1) then
			w = (var(:,i)%kth/dx(:,i))/(var(:,i)%kth/dx(:,i)+var(:,i+1)%kth/dx(:,i+1))

		   qadv(:,i) = rhow*cswat*q(:,i)*(w*(T(:,i)+Tzero)+(1.-w)*(T(:,i+1)+Tzero))
		   qadvya(:,i) =  rhow*cswat*qya(:,i)*(w*(T(:,i)+Tzero)+(1.-w)*(T(:,i+1)+Tzero))
		   qadvyb(:,i) =  rhow*cswat*qyb(:,i)*(w*(T(:,i)+Tzero)+(1.-w)*(T(:,i+1)+Tzero))
		   qadvTa(:,i) =  rhow*cswat*q(:,i)*w
		   qadvTb(:,i) =  rhow*cswat*q(:,i)*(1.-w)




		   qh(:,i) = qh(:,i) + qadv(:,i)
		   qhya(:,i) = qhya(:,i) + qadvya(:,i)
		   qhyb(:,i) = qhyb(:,i) + qadvyb(:,i)
		   qhTa(:,i) = qhTa(:,i) + qadvTa(:,i)
		   qhTb(:,i) = qhTb(:,i) + qadvTb(:,i)
	   endif
    enddo

    qh(:,n)   = zero
    qhya(:,n) = zero
    qhyb(:,n) = zero
    qhTa(:,n) = zero
    qhTb(:,n) = zero
    
    i_qh(:,1:n+1)   = qh(:,0:n)
    i_qhya(:,1:n+1) = qhya(:,0:n)
    i_qhyb(:,1:n+1) = qhyb(:,0:n)
    i_qhTa(:,1:n+1) = qhTa(:,0:n)
    i_qhTb(:,1:n+1) = qhTb(:,0:n)

	if (advection==1) then
		i_qadv(:,1:n+1)   = qadv(:,0:n)
		i_qadvya(:,1:n+1) = qadvya(:,0:n)
		i_qadvyb(:,1:n+1) = qadvyb(:,0:n)
		i_qadvTa(:,1:n+1) = qadvTa(:,0:n)
		i_qadvTb(:,1:n+1) = qadvTb(:,0:n)
	endif

  END SUBROUTINE getheatfluxes_2d

  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE hyofh(h, parin, K, Kh, phi)

    IMPLICIT NONE

    REAL(r_2),    INTENT(IN)  :: h
    TYPE(params), INTENT(IN)  :: parin
    REAL(r_2),    INTENT(OUT) :: K
    REAL(r_2),    INTENT(OUT) :: Kh
    REAL(r_2),    INTENT(OUT) :: phi
    ! Get soil water variables from h.
    ! Definitions of arguments:
    ! h   - matric head.
    ! K   - hydraulic conductivity.
    ! Kh  - derivative dK/dh.
    ! phi - matric flux potential (MFP).
    REAL(r_2) :: a

    a    =  -parin%lam * parin%eta
    K    =  parin%Ke * exp(a*log(h/parin%he))
    Kh   =  a * K / h ! exp and log may be faster than **
    phi  =  K * h / (one+a)

  END SUBROUTINE hyofh

  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE hyofS(S, Tsoil, parin, var)

    IMPLICIT NONE

    REAL(r_2),    INTENT(IN)    :: S
    REAL(r_2),    INTENT(IN)    :: Tsoil
    TYPE(params), INTENT(IN)    :: parin
    TYPE(vars),   INTENT(INOUT) :: var
    ! Get soil water variables from S.
    ! Definitions of arguments:
    ! S(1:ms)   - degree of saturation ("effective satn") of layers.
    ! ms        - no. of soil layers.
    ! jt(1:ms)  - layer soil type nos.
    ! var(1:ms) - other water vars of layers.
    REAL(r_2) :: lnS, v1, v2, theta, c, lnSat, v3, v4
    REAL(r_2) :: dhdS, lambda, crh, int
    REAL(r_2) :: dhdT, thetal_max
    REAL(r_2) :: Sliq, Sliqsat
    REAL(r_2) :: A, B, D, E, C1

    theta         = S*(parin%thre) + (parin%the - parin%thre)
    var%lambdav   = rlambda       ! latent heat of vaporisation
    var%lambdav   = 1.91846e6_r_2*((Tsoil+Tzero)/((Tsoil+Tzero)-33.91_r_2))**2  ! Henderson-Sellers, QJRMS, 1984
    var%lambdaf   = lambdaf        ! latent heat of fusion
    var%Tfrz      = Tfrz(S,parin%he,one/parin%lam)
    thetal_max    = thetalmax(Tsoil,S,parin%he,one/parin%lam,parin%thre,parin%the)
    var%dthetaldT = dthetalmaxdT(Tsoil,S,parin%he,1/parin%lam,parin%thre,parin%the)

    if (Tsoil < var%Tfrz) then ! ice
       var%iice   = 1
       var%thetai = (theta - thetal_max)*rhow/rhoi ! volumetric ice content
       var%thetal = thetal_max
       ! liquid water content, relative to saturation
       Sliq    = (var%thetal - (parin%the-parin%thre))/parin%thre
       ! saturated liquid water content (< 1 for frozen soil)
       Sliqsat = (thetal_max - (parin%the-parin%thre))/parin%thre
       lnSat      = log(Sliqsat)
       v1         = exp(-lnSat/parin%lam)
       ! air entry potential, flux matric potential and hydraulic conductivity for saturated frozen soil
       v2     = exp(parin%eta*lnSat)
       var%he   = parin%he*v1
       var%phie = parin%phie*v1*v2
       var%Ksat = parin%Ke*v2

       lnS    = log(Sliq)
       v3     = exp(-lnS/parin%lam)
       v4     = exp(parin%eta*lnS)
       var%h    = parin%he*v3  ! matric potential
       dhdT   =  (lambdaf*273.15_r_2 + csol*Rgas*(Tsoil + Tzero)**2)/(gravity*(Tsoil + Tzero)**2)
       dhdS   = zero
       var%K    = parin%Ke*v4
       var%KS   = zero
       var%KT   = var%dthetaldT * parin%Ke * parin%eta * (Sliq**(parin%eta-one))/parin%thre
       !MC! strange in original: used dthetaldT instead of v%dthetaldT -> KT=0
       !MC! set KT=0 to be consistent
       var%KT   = zero
       var%phi  = parin%phie*v3*v4
       var%phiS = zero
       var%phiT = parin%Ke * (var%h/parin%he)**(-parin%eta*parin%lam) * dhdT
       var%rh   = max(exp(Mw*gravity*var%h/Rgas/(Tsoil+Tzero)),rhmin)
    else ! no ice
       var%he     = parin%he 
       var%phie   = parin%phie
       var%Ksat   = parin%Ke
       var%iice   = 0                                 
       var%thetai = zero
       var%thetal = S*(parin%thre) + (parin%the - parin%thre)
       var%dthetaldT = zero
       dhdS   = zero
       Sliq    = S        ! liquid water content, relative to saturation
       Sliqsat = one
       lnS     = log(Sliq)
       v3      = exp(-lnS/parin%lam)
       v4      = exp(parin%eta*lnS)
    endif
    if (var%thetal < 1.e-12_r_2) then ! completely frozen
       var%dthetaldT = zero
       var%lambdav   = lambdas ! latent heat of sublimation
       var%KT        = zero
    endif
    if ((Tsoil >= var%Tfrz) .and. (S < one)) then ! unsaturated
       var%h    = parin%he*v3  ! matric potential
       dhdT     = zero
       dhdS     = -parin%he/parin%lam*S**(-one/parin%lam-one)
       var%K    = parin%Ke*v4
       var%phi  = parin%phie*v3*v4
       var%KS   = parin%eta*var%K/S
       var%KT   = zero
       var%phiS = (parin%eta-one/parin%lam)*var%phi/S
       var%phiT = zero
       var%rh   = max(exp(Mw*gravity*var%h/Rgas/(Tsoil+Tzero)),rhmin)
    endif
    if ((Tsoil >= var%Tfrz) .and. (S >= one)) then ! saturated
       var%h    = parin%he
       dhdT     = zero
       dhdS     = zero
       var%K    = parin%Ke
       var%phi  = parin%phie
       var%KS   = parin%eta*parin%Ke
       var%phiS = (parin%eta-one/parin%lam)*parin%phie
       var%phiT = zero
       var%rh   = one
       !MC! otherwise undefined
       var%KT   = zero
    endif

    !  variables required for vapour phase transfer
    theta  =  S*(parin%thre) + (parin%the - parin%thre)
    var%sl = slope_esat(Tsoil) * Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3 K-1
    c      = Mw*gravity/Rgas/(Tsoil+Tzero)
    lambda = parin%lam
    if (S < one) then
       var%rh = max(exp(Mw*gravity*var%h/Rgas/(Tsoil+Tzero)),rhmin)
    else
       var%rh = one
    endif
    crh        = c*var%rh
    var%hS     = dhdS
    var%rhS    = crh*dhdS
    var%cvsat  = esat(Tsoil)*Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3
    var%cv     = var%rh*var%cvsat
    var%cvS    = var%rhS *var%cvsat
    var%cvsatT = slope_esat(Tsoil)*Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3 K-1

    var%Dv    = Dva*parin%tortuosity*(parin%the-theta)  * ((Tsoil+Tzero)/Tzero)**1.88_r_2 ! m2 s-1
    int       = (-c*parin%he)**lambda * igamma(one-lambda,-c*var%h)
    var%phiv  = Dva*parin%tortuosity * (parin%thre*exp(c*var%h) -parin%thre*int)
    var%phivS = dhdS * Dva*parin%tortuosity * &
         ( (parin%thre)*c*exp(c*var%h) - parin%thre*c*(var%h/parin%he)**(-parin%lam)*exp(c*var%h) )
    var%kv    = var%Dv * var%cvsat *c * var%rh

    ! calculate v%kH as in Campbell (1985) p.32 eq. 4.20
    A          = 0.65_r_2 - 0.78_r_2*parin%rho/thousand + 0.60_r_2*(parin%rho/thousand)**2 ! need to substitute 1600 for rhob
    C1         = one + 2.6_r_2 * parin%clay**(-0.5_r_2)
    B          = 2.8_r_2 * (one-parin%thre)*theta
    D          = 0.03_r_2 + 0.7_r_2*(one-parin%thre)**2
    E          = four
    var%kH     = (A + B*theta-(A-D)*exp(-(C1*theta)**E))
    var%eta_th = one
    var%kE     = var%Dv*var%rh*var%sl*thousand*var%lambdav*var%eta_th
    var%kth    = var%kE + var%kH ! thermal conductivity of soil (includes contribution from vapour phase)

    var%csoil    = parin%css*parin%rho + rhow*cswat*var%thetal + rhoi*csice*var%thetai
    var%csoileff = var%csoil + rhow*lambdaf*var%dthetaldT ! increase effective heat capacity due to presence of ice

  END SUBROUTINE hyofS

  !**********************************************************************************************************************

  SUBROUTINE isosub(iso, c, p, f, fd)

    IMPLICIT NONE

    CHARACTER(LEN=2),       INTENT(IN)    :: iso
    REAL(r_2),               INTENT(IN)    :: c
    REAL(r_2), DIMENSION(:), INTENT(INOUT) :: p
    REAL(r_2),               INTENT(OUT)   :: f, fd
    ! Subroutine to get adsorbed solute (units/g soil) from concn in soil water
    ! according to chosen isotherm code ("Fr" for Freundlich, "La" for Langmuir
    ! and "Ll" for Langmuir-linear).
    ! Definitions of arguments:
    ! iso  - 2 character code.
    ! c    - concn in soil water.
    ! p(:) - isotherm parameters.
    ! f    - adsorbed mass/g soil.
    ! fc   - deriv of f wrt c (slope of isotherm curve).
    REAL(r_2) :: x

    select case (iso)
    case ("Fr")
       if (p(3)==zero) then ! linearise near zero
          p(3) = (0.01_r_2*dsmmax/p(1))**(one/p(2)) ! concn at 0.01*dsmmax
          p(4) = p(1)*p(3)**(p(2)-one) ! slope
       end if
       if (c < p(3)) then
          fd = p(4)
          f  = fd*c
       else
          x  = p(1)*exp((p(2)-one)*log(c))
          f  = x*c
          fd = p(2)*x
       end if
    case ("La")
       x  = one/(one+p(2)*c)
       f  = p(1)*c*x
       fd = p(1)*(x-p(2)*c*x**2)
    case ("Ll")
       x  = one/(one+p(2)*c)
       f  = p(1)*c*x+p(3)*c
       fd = p(1)*(x-p(2)*c*x**2)+p(3)
    case default
       write(2,*) "isosub: illegal isotherm type"
       stop
    end select

  END SUBROUTINE isosub

  !**********************************************************************************************************************

  SUBROUTINE massman_sparse_1d(aa, aah, bb, bbh, cc, cch, dd, ddh, ee, eeh, ff, ffh, gg, ggh, dy, dT, condition)

    USE cable_dimensions, ONLY: r_2, i_d
    USE cable_sli_numbers,       ONLY: zero, one

    IMPLICIT NONE

    ! in/out
    REAL(r_2), DIMENSION(:), INTENT(IN)  :: aa, aah, bb, bbh, ee, eeh, ff, ffh
    REAL(r_2), DIMENSION(:), INTENT(IN)  :: cc, cch, dd, ddh, gg, ggh
    REAL(r_2), DIMENSION(:), INTENT(OUT) :: dy, dT
    INTEGER(i_d),  OPTIONAL, INTENT(IN)  :: condition
    ! local
    INTEGER(i_d)                         :: n, n2
    REAL(r_2), DIMENSION(size(cc),2,2)   :: A, B, C
    REAL(r_2), DIMENSION(size(cc),2)     :: d, x
    ! for conditioning
    REAL(r_2), DIMENSION(2*size(cc))       :: lST, cST
    REAL(r_2), DIMENSION(size(cc))         :: lS, lT, cS, cT
    REAL(r_2), DIMENSION(2*size(cc)*2*size(cc)) :: allvec
    REAL(r_2), DIMENSION(2*size(cc),2*size(cc)) :: allmat
    REAL(r_2)    :: eps
    INTEGER(i_d) :: docond ! 0: no conditioning, 1: columns, 2: lines, 3: both
    ! CHARACTER(LEN=20) :: form1
    ! integer :: i, nn
    !
    ! check input sizes
    if (.not. all((/size(aa)+1,size(bb)+1,size(dd),size(ee)+1,size(ff)+1,size(gg)/) == size(cc))) then
       write(2,*) 'massman_sparse_1d error1: unequal humidity coeffs.'
       stop 'program terminated by massman_sparse_1d'
    end if
    if (.not. all((/size(aah)+1,size(bbh)+1,size(ddh),size(eeh)+1,size(ffh)+1,size(ggh)/) == size(cch))) then
       write(2,*) 'massman_sparse_1d error2: unequal temperature coeffs.'
       stop 'program terminated by massman_sparse_1d'
    end if
    if (size(cc) /= size(cch)) then
       write(2,*) 'massman_sparse_1d error3: unequal temperature and humidity coeffs.'
       stop 'program terminated by massman_sparse_1d'
    end if
    n = size(cc)
    if (present(condition)) then
       docond = condition
    else
       docond = 0
    endif
    !
    ! Overall matrix
    if (docond >= 1 .and. docond <= 3) then
       eps = epsilon(one)
       n2  = 2*n*2*n
       allvec(:) = zero
       allvec(1:n2:4*n+2)         = cc(1:n)
       allvec(2:n2:4*n+2)         = dd(1:n)
       allvec(3:n2-4*n:4*n+2)     = ee(1:n-1)
       allvec(4:n2-4*n:4*n+2)     = ff(1:n-1)
       allvec(2*n+1:n2:4*n+2)     = cch(1:n)
       allvec(2*n+2:n2:4*n+2)     = ddh(1:n)
       allvec(2*n+3:n2-4*n:4*n+2) = eeh(1:n-1)
       allvec(2*n+4:n2-4*n:4*n+2) = ffh(1:n-1)
       allvec(4*n+1:n2:4*n+2)     = aa(1:n-1)
       allvec(4*n+2:n2:4*n+2)     = bb(1:n-1)
       allvec(6*n+1:n2:4*n+2)     = aah(1:n-1)
       allvec(6*n+1:n2:4*n+2)     = bbh(1:n-1)
       allmat = reshape(allvec,shape=(/2*n,2*n/),order=(/2,1/))
    endif
    ! Get conditioning numbers
    select case (docond)
    case (1)
       cST = maxval(abs(allmat),1)
       where (cST < eps) cST = one
       cST     = one / cST
       cS(1:n) = cST(1:2*n-1:2)
       cT(1:n) = cST(2:2*n:2)
       lS(1:n) = one
       lT(1:n) = one
    case (2)
       lST = maxval(abs(allmat),2)
       where (lST < eps) lST = one
       lST     = one / lST
       lS(1:n) = lST(1:2*n-1:2)
       lT(1:n) = lST(2:2*n:2)
       cS(1:n) = one
       cT(1:n) = one
    case (3)
       lST = maxval(abs(allmat),2)
       where (lST < eps) lST = one
       lST     = one / lST
       lS(1:n) = lST(1:2*n-1:2)
       lT(1:n) = lST(2:2*n:2)
       allmat  = spread(lST,2,2*n)*allmat
       cST = maxval(abs(allmat),1)
       where (cST < eps) cST = one
       cST     = one / cST
       cS(1:n) = cST(1:2*n-1:2)
       cT(1:n) = cST(2:2*n:2)
    case default
       cS(1:n) = one
       cT(1:n) = one
       lS(1:n) = one
       lT(1:n) = one
    end select
    !
    ! fill matrices of generic thomas algorithm
    A(1,1:2,1:2) = zero
    A(2:n,1,1)   = aa(1:n-1)  * lS(2:n)   * cS(1:n-1)
    A(2:n,1,2)   = bb(1:n-1)  * lS(2:n)   * cT(1:n-1)
    A(2:n,2,1)   = aah(1:n-1) * lT(2:n)   * cS(1:n-1)
    A(2:n,2,2)   = bbh(1:n-1) * lT(2:n)   * cT(1:n-1)
    B(1:n,1,1)   = cc(1:n)    * lS(1:n)   * cS(1:n)
    B(1:n,1,2)   = dd(1:n)    * lS(1:n)   * cT(1:n)
    B(1:n,2,1)   = cch(1:n)   * lT(1:n)   * cS(1:n)
    B(1:n,2,2)   = ddh(1:n)   * lT(1:n)   * cT(1:n)
    C(1:n-1,1,1) = ee(1:n-1)  * lS(1:n-1) * cS(2:n)
    C(1:n-1,1,2) = ff(1:n-1)  * lS(1:n-1) * cT(2:n)
    C(1:n-1,2,1) = eeh(1:n-1) * lT(1:n-1) * cS(2:n)
    C(1:n-1,2,2) = ffh(1:n-1) * lT(1:n-1) * cT(2:n)
    C(n,1:2,1:2) = zero
    d(1:n,1)     = gg(1:n)    * lS(1:n)
    d(1:n,2)     = ggh(1:n)   * lT(1:n)
    !
    ! Call Generic Thomas algorithm
    call generic_thomas(n,A,B,C,d,x)
    dy(1:n) = x(1:n,1) * cS(1:n)
    dT(1:n) = x(1:n,2) * cT(1:n)
    !
  END SUBROUTINE massman_sparse_1d

  SUBROUTINE massman_sparse_2d(aa, aah, bb, bbh, cc, cch, dd, ddh, ee, eeh, ff, ffh, gg, ggh, dy, dT, condition)

    USE cable_dimensions, ONLY: r_2, i_d
    USE cable_sli_numbers,       ONLY: zero, one

    IMPLICIT NONE

    ! in/out
    REAL(r_2), DIMENSION(:,:), INTENT(IN)  :: aa, aah, bb, bbh, ee, eeh, ff, ffh
    REAL(r_2), DIMENSION(:,:), INTENT(IN)  :: cc, cch, dd, ddh, gg, ggh
    REAL(r_2), DIMENSION(:,:), INTENT(OUT) :: dy, dT
    INTEGER(i_d),    OPTIONAL, INTENT(IN)  :: condition
    ! local
    INTEGER(i_d)                                        :: n, mp
    INTEGER(i_d)                                        :: n2
    REAL(r_2), DIMENSION(1:size(cc,1),1:size(cc,2),2,2) :: A, B, C
    REAL(r_2), DIMENSION(1:size(cc,1),1:size(cc,2),2)   :: d, x
    ! for conditioning
    REAL(r_2), DIMENSION(1:size(cc,1),2*size(cc,2))     :: lST, cST
    REAL(r_2), DIMENSION(1:size(cc,1),size(cc,2))       :: lS, lT, cS, cT
    REAL(r_2), DIMENSION(1:size(cc,1),2*size(cc,2)*2*size(cc,2)) :: allvec
    REAL(r_2), DIMENSION(1:size(cc,1),2*size(cc,2),2*size(cc,2)) :: allmat
    REAL(r_2)    :: eps
    INTEGER(i_d) :: docond ! 0: no conditioning, 1: columns, 2: lines, 3: both
    ! CHARACTER(LEN=20) :: form1
    ! integer :: i, k, nn
    !
    ! check input sizes
    if (.not. all((/size(aa,1)+1,size(bb,1)+1,size(dd,1),size(ee,1)+1,size(ff,1)+1,size(gg,1)/) == size(cc,1))) then
       write(2,*) 'massman_sparse_2d error1: unequal humidity coeffs (1st dim).'
       stop 'program terminated by massman_sparse_2d'
    end if
    if (.not. all((/size(aah,1)+1,size(bbh,1)+1,size(ddh,1),size(eeh,1)+1,size(ffh,1)+1,size(ggh,1)/) == size(cch,1))) then
       write(2,*) 'massman_sparse_2d error2: unequal temperature coeffs (1st dim).'
       stop 'program terminated by massman_sparse_2d'
    end if
    if (size(cc,1) /= size(cch,1)) then
       write(2,*) 'massman_sparse_2d error3: unequal temperature and humidity coeffs (1st dim).'
       stop 'program terminated by massman_sparse_2d'
    end if
    if (.not. all((/size(aa,2)+1,size(bb,2)+1,size(dd,2),size(ee,2)+1,size(ff,2)+1,size(gg,2)/) == size(cc,2))) then
       write(2,*) 'massman_sparse_2d error4: unequal humidity coeffs (2nd dim).'
       stop 'program terminated by massman_sparse_2d'
    end if
    if (.not. all((/size(aah,2)+1,size(bbh,2)+1,size(ddh,2),size(eeh,2)+1,size(ffh,2)+1,size(ggh,2)/) == size(cch,2))) then
       write(2,*) 'massman_sparse_2d error5: unequal temperature coeffs (2nd dim).'
       stop 'program terminated by massman_sparse_2d'
    end if
    if (size(cc,2) /= size(cch,2)) then
       write(2,*) 'massman_sparse_2d error6: unequal temperature and humidity coeffs (2nd dim).'
       stop 'program terminated by massman_sparse_2d'
    end if

    mp   = size(cc,1)
    n    = size(cc,2)
    if (present(condition)) then
       docond = condition
    else
       docond = 0
    endif
    !
    ! Overall matrix
    if (docond >= 1 .and. docond <= 3) then
       eps = epsilon(one)
       n2  = 2*n*2*n
       allvec(1:mp,:) = zero
       allvec(1:mp,1:n2:4*n+2)         = cc(1:mp,1:n)
       allvec(1:mp,2:n2:4*n+2)         = dd(1:mp,1:n)
       allvec(1:mp,3:n2-4*n:4*n+2)     = ee(1:mp,1:n-1)
       allvec(1:mp,4:n2-4*n:4*n+2)     = ff(1:mp,1:n-1)
       allvec(1:mp,2*n+1:n2:4*n+2)     = cch(1:mp,1:n)
       allvec(1:mp,2*n+2:n2:4*n+2)     = ddh(1:mp,1:n)
       allvec(1:mp,2*n+3:n2-4*n:4*n+2) = eeh(1:mp,1:n-1)
       allvec(1:mp,2*n+4:n2-4*n:4*n+2) = ffh(1:mp,1:n-1)
       allvec(1:mp,4*n+1:n2:4*n+2)     = aa(1:mp,1:n-1)
       allvec(1:mp,4*n+2:n2:4*n+2)     = bb(1:mp,1:n-1)
       allvec(1:mp,6*n+1:n2:4*n+2)     = aah(1:mp,1:n-1)
       allvec(1:mp,6*n+1:n2:4*n+2)     = bbh(1:mp,1:n-1)
       allmat = reshape(allvec,shape=(/mp,2*n,2*n/),order=(/1,3,2/))
    endif
    ! Get conditioning numbers
    select case (docond)
    case (1)
       cST = maxval(abs(allmat),2)
       where (cST < eps) cST = one
       cST     = one / cST
       cS(1:mp,1:n) = cST(1:mp,1:2*n-1:2)
       cT(1:mp,1:n) = cST(1:mp,2:2*n:2)
       lS(1:mp,1:n) = one
       lT(1:mp,1:n) = one
    case (2)
       lST = maxval(abs(allmat),3)
       where (lST < eps) lST = one
       lST  = one / lST
       lS(1:mp,1:n) = lST(1:mp,1:2*n-1:2)
       lT(1:mp,1:n) = lST(1:mp,2:2*n:2)
       cS(1:mp,1:n) = one
       cT(1:mp,1:n) = one
    case (3)
       lST = maxval(abs(allmat),3)
       where (lST < eps) lST = one
       lST  = one / lST
       lS(1:mp,1:n) = lST(1:mp,1:2*n-1:2)
       lT(1:mp,1:n) = lST(1:mp,2:2*n:2)
       allmat  = spread(lST,2,2*n)*allmat
       cST = maxval(abs(allmat),2)
       where (cST < eps) cST = one
       cST = one / cST
       cS(1:mp,1:n) = cST(1:mp,1:2*n-1:2)
       cT(1:mp,1:n) = cST(1:mp,2:2*n:2)
    case default
       cS(1:mp,1:n) = one
       cT(1:mp,1:n) = one
       lS(1:mp,1:n) = one
       lT(1:mp,1:n) = one
    end select
    !
    ! fill matrices of generic thomas algorithm
    A(1:mp,1,1:2,1:2) = zero
    A(1:mp,2:n,1,1)   = aa(1:mp,1:n-1)  * lS(1:mp,2:n)   * cS(1:mp,1:n-1)
    A(1:mp,2:n,1,2)   = bb(1:mp,1:n-1)  * lS(1:mp,2:n)   * cT(1:mp,1:n-1)
    A(1:mp,2:n,2,1)   = aah(1:mp,1:n-1) * lT(1:mp,2:n)   * cS(1:mp,1:n-1)
    A(1:mp,2:n,2,2)   = bbh(1:mp,1:n-1) * lT(1:mp,2:n)   * cT(1:mp,1:n-1)
    B(1:mp,1:n,1,1)   = cc(1:mp,1:n)    * lS(1:mp,1:n)   * cS(1:mp,1:n)
    B(1:mp,1:n,1,2)   = dd(1:mp,1:n)    * lS(1:mp,1:n)   * cT(1:mp,1:n)
    B(1:mp,1:n,2,1)   = cch(1:mp,1:n)   * lT(1:mp,1:n)   * cS(1:mp,1:n)
    B(1:mp,1:n,2,2)   = ddh(1:mp,1:n)   * lT(1:mp,1:n)   * cT(1:mp,1:n)
    C(1:mp,1:n-1,1,1) = ee(1:mp,1:n-1)  * lS(1:mp,1:n-1) * cS(1:mp,2:n)
    C(1:mp,1:n-1,1,2) = ff(1:mp,1:n-1)  * lS(1:mp,1:n-1) * cT(1:mp,2:n)
    C(1:mp,1:n-1,2,1) = eeh(1:mp,1:n-1) * lT(1:mp,1:n-1) * cS(1:mp,2:n)
    C(1:mp,1:n-1,2,2) = ffh(1:mp,1:n-1) * lT(1:mp,1:n-1) * cT(1:mp,2:n)
    C(1:mp,n,1:2,1:2) = zero
    d(1:mp,1:n,1)     = gg(1:mp,1:n)    * lS(1:mp,1:n)
    d(1:mp,1:n,2)     = ggh(1:mp,1:n)   * lT(1:mp,1:n)
    !
    ! Call Generic Thomas algorithm
    call generic_thomas(mp,n,A,B,C,d,x)
    dy(1:mp,1:n) = x(1:mp,1:n,1) * cS(1:mp,1:n)
    dT(1:mp,1:n) = x(1:mp,1:n,2) * cT(1:mp,1:n)
    !
  END SUBROUTINE massman_sparse_2d

  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE litter_props(S, Tsoil, vlit, plit, h0)

    IMPLICIT NONE

    REAL(r_2),    INTENT(IN)    :: S
    REAL(r_2),    INTENT(IN)    :: Tsoil
    TYPE(vars),   INTENT(INOUT) :: vlit
    TYPE(params), INTENT(IN)    :: plit
    REAL(r_2),    INTENT(IN)    :: h0
    ! Get soil water variables from S.
    ! Definitions of arguments:
    ! S(1:ms)    - degree of saturation ("effective satn") of layers.
    ! ms        - no. of soil layers.
    ! jt(1:ms)    - layer soil type nos.
    ! var(1:ms) - other water vars of layers.
    REAL(r_2), PARAMETER :: u = 1.0 ! wind speed at litter surface
    REAL(r_2) :: lnS, theta, c
    REAL(r_2) :: dhdS, sl, lambda
    REAL(r_2) :: rhoL, copo, c_w ! params for thermal vapour transfer enhancement factor (Campbell 1985)
    REAL(r_2) :: chi, DT0
    
    lnS     = log(S)
    sl      = slope_esat(Tsoil)*Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3 K-1
    vlit%sl = sl
    theta   = S*(plit%thre) + (plit%the - plit%thre)
    c       = Mw*gravity/Rgas/(Tsoil+Tzero)
    lambda  = plit%lam
    rhoL    = plit%rho
    dhdS    = -plit%he/plit%lam*S**(-one/plit%lam-one)*(thousand/rhoL*plit%the)**(-one/plit%lam)
    vlit%hS = dhdS
    ! Mathews (2006), A process-based model of offine fuel moisture,
    !                 International Journal of Wildland Fire 15,155–168
    chi     = 2.08_r_2+u*2.38_r_2 ! (Eq. 45, Tab. 1)
    DT0     = Dva*exp(u*2.6_r_2) ! (Eq. 46, Tab. 1)
    vlit%Dv = DT0*exp(-half*chi) ! heat and vapour diffusivity half-way through layer (Eq. 11)
    if (S < one) then
       vlit%h = plit%he*(thousand/rhoL*S*plit%the)**(-one/plit%lam)
       vlit%K = zero
    else
       vlit%K = vlit%Ke
       vlit%h = plit%he
    endif
    vlit%rh     = max(exp(Mw*gravity*vlit%h/Rgas/(Tsoil+Tzero)),rhmin)
    vlit%cvsat  = esat(Tsoil)*Mw/thousand/Rgas/(Tsoil+Tzero) ! m3 m-3
    vlit%cv     = vlit%cvsat*vlit%rh
    vlit%cvS    = vlit%rhS *vlit%cvsat
    vlit%KS     = zero
    vlit%phiv   = vlit%Dv*vlit%cvsat*vlit%rh
    vlit%phivS  = dhdS*vlit%phiv*c
    vlit%phi    = vlit%phiv
    vlit%phiS   = vlit%phivS
    vlit%rhS    = dhdS*c*vlit%rh
    ! sensible heat conductivity of wet soil W m-1 K-1 (Matthews 2006)
    vlit%kH     = 0.2_r_2 + 0.14_r_2*theta*thousand/rhoL
    copo        = 1932._r_2
    c_w         = 4.185e6_r_2
    vlit%csoil  = copo*rhoL + c_w*theta + c_w*h0 !volumetric heat capacity of wet soil
    vlit%eta_th = one !enhancement factor for transport of water vapour due to a temperature gradient
    vlit%kE     = vlit%Dv*vlit%rh*sl*thousand*rlambda
    vlit%kth    = vlit%kE + vlit%kH ! thermal conductivity of soil 

  END SUBROUTINE litter_props

  !**********************************************************************************************************************

  ELEMENTAL PURE SUBROUTINE potential_evap(Rn, rbh, rbw, Ta, rha, Tsoil, k, dz,  Ts, E, H, G, &
       dEdrha, dEdTa, dEdTsoil, dGdTa, dGdTsoil)

    ! Pennman-Monteith equation, with additional account for heat flux into the surface

    IMPLICIT NONE

    REAL(r_2), INTENT(IN)  :: Rn, rbh, rbw, Ta, rha, Tsoil, k, dz
    REAL(r_2), INTENT(OUT) :: Ts, E, H, G, dEdrha, dEdTa, dEdTsoil, dGdTa, dGdTsoil
 
    REAL(r_2) :: s, ea, dEdea, dEdesat, dTsdTa, dEdDa, Da
    REAL(r_2), PARAMETER :: gamma = 62.0 ! psychrometric constant

    s  = slope_esat(Ta)
    ea = esat(Ta) * rha
    Da = esat(Ta) - ea

    E  = (rhocp*(Da*(k*rbh + dz*rhocp) + rbh*s*(dz*Rn + k*(-Ta + Tsoil)))) / &
         (gamma*rbw*(k*rbh + dz*rhocp) + dz*rbh*rhocp*s)
    Ts = Ta + E*gamma*rbw/s/rhocp - Da/s
    H  = rhocp*(Ts - Ta)/rbh
    G  = k*(Ts-Tsoil)/dz

    dEdDa    = (-(k*rbh*rhocp) - dz*rhocp**2)/(gamma*k*rbh*rbw + dz*gamma*rbw*rhocp + dz*rbh*rhocp*s)
    dEdea    = -dEdDa
    dEdesat  = dEdea
    dEdrha   = dEdea *esat(Ta)
    dEdTa    = (k*rbh*rhocp*s)/(gamma*k*rbh*rbw + dz*gamma*rbw*rhocp + dz*rbh*rhocp*s) + dEdesat *slope_esat(Ta)
    dEdTsoil = -((k*rbh*rhocp*s)/(gamma*k*rbh*rbw + dz*gamma*rbw*rhocp + dz*rbh*rhocp*s))

    dTsdTa   = (-(dz*gamma*rbw*rhocp) - dz*rbh*rhocp*s)/(gamma*k*rbh*rbw + dz*gamma*rbw*rhocp + dz*rbh*rhocp*s)

    dGdTa    = k/dz * dTsdTa
    dGdTsoil = -k/dz

  END SUBROUTINE potential_evap

  !**********************************************************************************************************************

  SUBROUTINE setlitterpar(mp, soil)

    IMPLICIT NONE

    INTEGER(i_d),              INTENT(IN) :: mp
    TYPE(soil_parameter_type), INTENT(IN) :: soil

    allocate(plit(mp))
    allocate(dxL(mp))

    plit%the   = 0.09_r_2
    plit%thre  = 0.09_r_2
    plit%he    = -35.0_r_2
    plit%lam   = 1._r_2/2.4_r_2
    plit%Ke    = zero
    plit%eta   = zero
    plit%KSe   = zero
    plit%phie  = zero
    plit%phiSe = zero
    plit%rho   = 63.5_r_2
    dxL        = zero            ! litter params
    dxL        = real(soil%clitt,r_2)*two/plit%rho*0.1_r_2

  END SUBROUTINE setlitterpar

  !**********************************************************************************************************************

  SUBROUTINE setpar(mp, ms, x2dx, soil)

    IMPLICIT NONE

    INTEGER(i_d),              INTENT(IN)    :: mp
    INTEGER(i_d),              INTENT(IN)    :: ms
    REAL(r_2), DIMENSION(:,:), INTENT(IN)    :: x2dx
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    REAL(r_2), DIMENSION(1:mp) :: depthA
    INTEGER(i_d) :: i

    allocate(par(mp,ms))

    depthA(:) = real(soil%DepthA(:),r_2)
    do i=1, ms
       where (x2dx(:,i) < depthA(:))
          par(:,i)%isahorizon = one
          par(:,i)%isbhorizon = zero
          par(:,i)%thw        = real(soil%swilt(:),r_2)
          par(:,i)%the        = real(soil%ssat(:),r_2)
          par(:,i)%thr        = zero    ! Assume theta_r is zero
          par(:,i)%thre       = real(soil%ssat(:),r_2) - par(:,i)%thr  
          par(:,i)%he         = real(soil%sucs(:),r_2)
          par(:,i)%Ke         = real(soil%hyds(:),r_2)
          par(:,i)%lam        = one/real(soil%bch(:),r_2)
          !p%eta                 = 9.14 !(special for Braud test-cases)
          par(:,i)%eta        = two/par(:,i)%lam + two + one
          par(:,i)%KSe        = par(:,i)%eta * par(:,i)%Ke    ! dK/dS at he
          par(:,i)%phie       = par(:,i)%Ke * par(:,i)%he / (one - par(:,i)%lam * par(:,i)%eta) ! MFP at he
          par(:,i)%phiSe      = (par(:,i)%eta - one/par(:,i)%lam) * par(:,i)%phie    ! dphi/dS at he
          par(:,i)%kd         = real(soil%cnsd(:),r_2)
          par(:,i)%css        = real(soil%css(:),r_2)
          par(:,i)%rho        = real(soil%rhosoil(:),r_2)
          par(:,i)%tortuosity = 0.67_r_2
          par(:,i)%clay       = real(soil%clay(:),r_2)
       elsewhere
          par(:,i)%isahorizon = zero
          par(:,i)%isbhorizon = one
          par(:,i)%thw        = real(soil%swiltB(:),r_2)
          par(:,i)%the        = real(soil%ssatB(:),r_2)
          par(:,i)%thr        = zero    ! Assume theta_r is zero
          par(:,i)%thre       = real(soil%ssatB(:),r_2) - par(:,i)%thr 
          par(:,i)%he         = real(soil%sucsB(:),r_2)
          par(:,i)%Ke         = real(soil%hydsB(:),r_2)
          par(:,i)%lam        = one/real(soil%bchB(:),r_2)
          !p%eta                  = 9.14 !(special for Braud test-cases)
          par(:,i)%eta        = two/par(:,i)%lam + two + one
          par(:,i)%KSe        = par(:,i)%eta * par(:,i)%Ke    ! dK/dS at he
          par(:,i)%phie       = par(:,i)%Ke * par(:,i)%he / (one - par(:,i)%lam * par(:,i)%eta) ! MFP at he
          par(:,i)%phiSe      = (par(:,i)%eta - one/par(:,i)%lam) * par(:,i)%phie    ! dphi/dS at he
          par(:,i)%kd         = real(soil%cnsdB(:),r_2)
          par(:,i)%css        = real(soil%cssB(:),r_2)
          par(:,i)%rho        = real(soil%rhosoilB(:),r_2)
          par(:,i)%tortuosity = 0.67_r_2
          par(:,i)%clay       = real(soil%clayB(:),r_2)
       endwhere
    enddo

  END SUBROUTINE setpar

  !**********************************************************************************************************************

  SUBROUTINE setsol(mp)

    IMPLICIT NONE

    INTEGER(i_d), INTENT(IN) :: mp

    allocate(sol(mp))

    sol%T1      = zero
    sol%Ta      = zero
    sol%cva     = zero
    sol%Rnet    = zero
    sol%hr1     = zero
    sol%hra     = zero
    sol%Dv      = zero
    sol%gv      = zero
    sol%gh      = zero
    sol%Dh      = zero
    sol%dz      = zero
    sol%phie    = zero
    sol%he      = zero
    sol%K1      = zero
    sol%eta     = zero
    sol%lambda  = zero
    sol%Ks      = zero
    sol%lambdav = zero

  END SUBROUTINE setsol

  !**********************************************************************************************************************

  SUBROUTINE setx(mp, ms, soil)

    IMPLICIT NONE

    INTEGER(i_d),              INTENT(IN) :: mp
    INTEGER(i_d),              INTENT(IN) :: ms
    TYPE(soil_parameter_type), INTENT(IN) :: soil
    
    REAL(r_2), DIMENSION(ms) :: tmp1d
    INTEGER :: i

    allocate(x(mp,ms))
    allocate(dx(mp,ms))
    
    ! cumulative soil layer depths = bottom of soil layers
    tmp1d(1) = real(soil%zse(1),r_2)
    do i=2, ms
       tmp1d(i) = tmp1d(i-1) + real(soil%zse(i),r_2)
    end do

    ! soil layer depth
    dx(1:mp,1:ms) = spread(real(soil%zse(1:ms),r_2),1,mp)
    ! bottom of each soil layer
    x(1:mp,1:ms)  = spread(tmp1d(1:ms),1,mp)

  END SUBROUTINE setx

  !**********************************************************************************************************************

  SUBROUTINE tri_1d(ns, n, aa, bb, cc, dd, dy)

    IMPLICIT NONE

    INTEGER(i_d),                 INTENT(IN)    :: ns
    INTEGER(i_d),                 INTENT(IN)    :: n
    REAL(r_2),    DIMENSION(0:n), INTENT(IN)    :: aa
    REAL(r_2),    DIMENSION(0:n), INTENT(IN)    :: cc
    REAL(r_2),    DIMENSION(0:n), INTENT(IN)    :: dd
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: bb
    REAL(r_2),    DIMENSION(0:n), INTENT(INOUT) :: dy
    ! Solves tridiag set of linear eqns. Coeff arrays aa and cc left intact.
    ! Definitions of arguments:
    ! ns      - start index for eqns.
    ! n       - end index.
    ! aa(0:n) - coeffs below diagonal; ns+1:n used.
    ! bb(0:n) - coeffs on diagonal; ns:n used.
    ! cc(0:n) - coeffs above diagonal; ns:n-1 used.
    ! dd(0:n) - rhs coeffs; ns:n used.
    ! ee(0:n) - work space.
    ! dy(0:n) - solution in ns:n.
    REAL(r_2),   DIMENSION(0:n) :: ee
    INTEGER(i_d)                :: i

    dy(ns) = dd(ns) ! decomposition and forward substitution
    do i=ns, n-1
       ee(i)   = cc(i)/bb(i)
       dy(i)   = dy(i)/bb(i)
       bb(i+1) = bb(i+1)-aa(i+1)*ee(i)
       dy(i+1) = dd(i+1)-aa(i+1)*dy(i)
    end do

    dy(n) = dy(n)/bb(n) ! back substitution
    do i=n-1, ns, -1
       dy(i) = dy(i)-ee(i)*dy(i+1)
    end do

  END SUBROUTINE tri_1d

  SUBROUTINE tri_2d(mp, ns, n, aa, bb, cc, dd, dy)

    IMPLICIT NONE

    INTEGER(i_d),                      INTENT(IN)    :: mp
    INTEGER(i_d),                      INTENT(IN)    :: ns
    INTEGER(i_d),                      INTENT(IN)    :: n
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(IN)    :: aa
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(IN)    :: cc
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(IN)    :: dd
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(INOUT) :: bb
    REAL(r_2),    DIMENSION(1:mp,0:n), INTENT(INOUT) :: dy
    ! Solves tridiag set of linear eqns. Coeff arrays aa and cc left intact.
    ! Definitions of arguments:
    ! mp      - number of land patches
    ! ns      - start index for eqns.
    ! n       - end index.
    ! aa(0:n) - coeffs below diagonal; ns+1:n used.
    ! bb(0:n) - coeffs on diagonal; ns:n used.
    ! cc(0:n) - coeffs above diagonal; ns:n-1 used.
    ! dd(0:n) - rhs coeffs; ns:n used.
    ! ee(0:n) - work space.
    ! dy(0:n) - solution in ns:n.
    REAL(r_2),   DIMENSION(1:mp,0:n) :: ee
    INTEGER(i_d)                     :: i

    dy(:,ns) = dd(:,ns) ! decomposition and forward substitution
    do i=ns, n-1
       ee(:,i)   = cc(:,i)/bb(:,i)
       dy(:,i)   = dy(:,i)/bb(:,i)
       bb(:,i+1) = bb(:,i+1)-aa(:,i+1)*ee(:,i)
       dy(:,i+1) = dd(:,i+1)-aa(:,i+1)*dy(:,i)
    end do

    dy(:,n) = dy(:,n)/bb(:,n) ! back substitution
    do i=n-1, ns, -1
       dy(:,i) = dy(:,i)-ee(:,i)*dy(:,i+1)
    end do

  END SUBROUTINE tri_2d

  !**********************************************************************************************************************
  ! FUNCTIONS - alphabetical
  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION csat(T)
    !returns  sat vapour pressure curve in kg m-3
    USE cable_sli_numbers, ONLY: Tzero, Rgas, Mw, esata, esatb, esatc

    IMPLICIT NONE

    real(r_2), intent(in) :: T

    csat = esata * exp(esatb*T/(T+esatc)) * Mw/Rgas/(T+Tzero)

  END FUNCTION csat

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION csoil(thetal,thetai,css,rho)
    ! determines heat capacity of soil (Jkg-1K-1)
    USE cable_sli_numbers, ONLY: rhow, cswat, rhoi, csice

    IMPLICIT NONE

    real(r_2), intent(in) ::  thetal, thetai, css, rho

    csoil = css*rho + rhow*cswat*thetal + rhoi*csice*thetai

  END FUNCTION csoil

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION dthetalmaxdT(Tin,S,he,b,thre,the)
    ! determines derivative of thetalmax wrt T
    USE cable_sli_numbers, ONLY: gravity, lambdaf, csol, Rgas, Tzero

    IMPLICIT NONE

    real(r_2), intent(in) :: Tin,S,he,b,thre,the
    real(r_2)             :: dthetaldh, h, psi, PI, T, thetalmax

    T = min(Tfrz(S,he,b),Tin)
    PI  = -csol*Rgas*(T+Tzero)/gravity ! osmotic potential (m)
    psi = lambdaf*T/(gravity*(T+Tzero))! total matric potential in presence of ice
    h   = psi-PI       ! moisture potential in presence of ice

    thetalmax = thre*(h/he)**(-1/b) + (the-thre)
    dthetaldh = -thetalmax/b/h

    dthetalmaxdT = dthetaldh*(lambdaf/(T+Tzero)/gravity-lambdaf*T/gravity/(T+Tzero)**2+csol*Rgas/gravity)

  END FUNCTION dthetalmaxdT

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION dthetalmaxdTh(Tin,S,he,b,thre,the)
    ! determines derivative of thetalmax wrt T
    ! (uses dthetaldh defined in terms of S ...required for crossing freezing point)
    USE cable_sli_numbers, ONLY: gravity, lambdaf, csol, Rgas, Tzero

    IMPLICIT NONE

    real(r_2), INTENT(in) :: Tin,S,he,b,thre,the
    real(r_2)             :: dthetaldh, h, theta, T

    T     = min(Tfrz(S,he,b),Tin)
    h     = he * S**(-b)
    theta = S*thre + (the-thre)

    dthetaldh = -theta/b/h    

    dthetalmaxdTh = dthetaldh*(lambdaf/(T+Tzero)/gravity-lambdaf*T/gravity/(T+Tzero)**2+csol*Rgas/gravity)

  END FUNCTION dthetalmaxdTh

  !*********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION esat(T)
    !returns sat vapour pressure curve in Pa
    USE cable_sli_numbers, ONLY: esata, esatb, esatc

    IMPLICIT NONE

    real(r_2), intent(in) :: T 

    esat = esata * exp(esatb*T/(T+esatc))

  END FUNCTION esat

  !**********************************************************************************************************************

  REAL(SP) ELEMENTAL PURE FUNCTION sgammln(z)
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !  Reference:
    !       Lanczos, C. 'A precision approximation of the gamma
    !               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996
    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: z
    ! Local variables
    REAL(r_2), PARAMETER :: a(9) = (/ &
         0.9999999999995183_r_2, 676.5203681218835_r_2, -1259.139216722289_r_2, &
         771.3234287757674_r_2, -176.6150291498386_r_2, 12.50734324009056_r_2, &
         -0.1385710331296526_r_2, 0.9934937113930748E-05_r_2, 0.1659470187408462E-06_r_2 /)
    REAL(r_2), PARAMETER :: zero = 0.0_r_2
    REAL(r_2), PARAMETER :: one = 1.0_r_2
    REAL(r_2), PARAMETER :: lnsqrt2pi =  0.9189385332046727_r_2
    REAL(r_2), PARAMETER :: half = 0.5_r_2
    REAL(r_2), PARAMETER :: sixpt5 = 6.5_r_2
    REAL(r_2), PARAMETER :: seven = 7.0_r_2
    REAL(r_2)    :: tmp, tmpgammln, ztmp
    INTEGER(i_d) :: j

    ztmp = real(z,r_2)
    tmpgammln = zero
    tmp = ztmp + seven
    DO j = 9, 2, -1
       tmpgammln = tmpgammln + a(j)/tmp
       tmp = tmp - one
    END DO
    tmpgammln = tmpgammln + a(1)
    tmpgammln = LOG(tmpgammln) + lnsqrt2pi - (ztmp + sixpt5) + (ztmp - half)*LOG(ztmp + sixpt5)
    sgammln   = real(tmpgammln,SP)
    RETURN

  END FUNCTION sgammln

  REAL(r_2) ELEMENTAL PURE FUNCTION dgammln(z)
    !  Uses Lanczos-type approximation to ln(gamma) for z > 0.
    !  Reference:
    !       Lanczos, C. 'A precision approximation of the gamma
    !               function', J. SIAM Numer. Anal., B, 1, 86-96, 1964.
    !  Accuracy: About 14 significant digits except for small regions
    !            in the vicinity of 1 and 2.
    !  Programmer: Alan Miller
    !              1 Creswick Street, Brighton, Vic. 3187, Australia
    !  e-mail: amiller @ bigpond.net.au
    !  Latest revision - 14 October 1996
    IMPLICIT NONE

    REAL(r_2), INTENT(IN) :: z
    ! Local variables
    REAL(r_2), PARAMETER :: a(9) = (/ &
         0.9999999999995183_r_2, 676.5203681218835_r_2, -1259.139216722289_r_2, &
         771.3234287757674_r_2, -176.6150291498386_r_2, 12.50734324009056_r_2, &
         -0.1385710331296526_r_2, 0.9934937113930748E-05_r_2, 0.1659470187408462E-06_r_2 /)
    REAL(r_2), PARAMETER :: zero = 0.0_r_2
    REAL(r_2), PARAMETER :: one = 1.0_r_2
    REAL(r_2), PARAMETER :: lnsqrt2pi =  0.9189385332046727_r_2
    REAL(r_2), PARAMETER :: half = 0.5_r_2
    REAL(r_2), PARAMETER :: sixpt5 = 6.5_r_2
    REAL(r_2), PARAMETER :: seven = 7.0_r_2
    REAL(r_2)    :: tmp, tmpgammln
    INTEGER(i_d) :: j

    tmpgammln = zero
    tmp = z + seven
    DO j = 9, 2, -1
       tmpgammln = tmpgammln + a(j)/tmp
       tmp = tmp - one
    END DO
    tmpgammln = tmpgammln + a(1)
    dgammln = LOG(tmpgammln) + lnsqrt2pi - (z + sixpt5) + (z - half)*LOG(z + sixpt5)
    RETURN

  END FUNCTION dgammln

  !**********************************************************************************************************************

  REAL(SP) ELEMENTAL PURE FUNCTION sgcf(a,x)

    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a,x
    INTEGER(i_d), PARAMETER :: ITMAX=100
    REAL(SP),     PARAMETER :: EPS=epsilon(x)
    REAL(SP),     PARAMETER :: FPMIN=tiny(x)/EPS
    INTEGER(i_d) :: i
    REAL(SP)     :: an, b, c, d, del, h

    if (x == 0.0) then
       sgcf=1.0
       RETURN
    end if
    b = x + 1.0_sp - a
    c = 1.0_sp/FPMIN
    d = 1.0_sp/b
    h = d
    do i=1, ITMAX
       an  = -i*(i-a)
       b   = b + 2.0_sp
       d   = an*d + b
       if (abs(d)<FPMIN) d=FPMIN
       c   = b + an/c
       if (abs(c)<FPMIN) c=FPMIN
       d   = 1.0_sp/d
       del = d*c
       h   = h*del
       if (abs(del-1.0_sp) <= EPS) exit
    end do
    sgcf = exp(-x + a*log(x) - gammln(a)) * h

  END FUNCTION sgcf


  REAL(r_2) ELEMENTAL PURE FUNCTION dgcf(a,x)

    IMPLICIT NONE

    REAL(r_2), INTENT(IN) :: a,x
    INTEGER(i_d), PARAMETER :: ITMAX=100
    REAL(r_2),    PARAMETER :: EPS=epsilon(x)
    REAL(r_2),    PARAMETER :: FPMIN=tiny(x)/EPS
    INTEGER(i_d) :: i
    REAL(r_2)     :: an, b, c, d, del, h

    if (x == 0.0) then
       dgcf=1.0
       RETURN
    end if
    b = x + 1.0_r_2 - a
    c = 1.0_r_2/FPMIN
    d = 1.0_r_2/b
    h = d
    do i=1, ITMAX
       an  = -i*(i-a)
       b   = b + 2.0_r_2
       d   = an*d + b
       if (abs(d)<FPMIN) d=FPMIN
       c   = b + an/c
       if (abs(c)<FPMIN) c=FPMIN
       d   = 1.0_r_2/d
       del = d*c
       h   = h*del
       if (abs(del-1.0_r_2) <= EPS) exit
    end do
    dgcf = exp(-x + a*log(x) - gammln(a)) * h

  END FUNCTION dgcf

  !**********************************************************************************************************************

  REAL(SP) ELEMENTAL PURE FUNCTION sgser(a,x)

    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a, x
    INTEGER(i_d), PARAMETER :: ITMAX=100
    REAL(SP),     PARAMETER :: EPS=epsilon(x)
    INTEGER(i_d) :: n
    REAL(SP)     :: ap, del, summ

    if (x == 0.0) then
       sgser = 0.0
       RETURN
    end if
    ap   = a
    summ = 1.0_sp/a
    del  = summ
    do n=1, ITMAX
       ap   = ap + 1.0_sp
       del  = del*x/ap
       summ = summ + del
       if (abs(del) < abs(summ)*EPS) exit
    end do
    sgser=summ*exp(-x+a*log(x)-gammln(a))

  END FUNCTION sgser


  REAL(r_2) ELEMENTAL PURE FUNCTION dgser(a,x)

    IMPLICIT NONE

    REAL(r_2), INTENT(IN) :: a, x
    INTEGER(i_d), PARAMETER :: ITMAX=100
    REAL(r_2),     PARAMETER :: EPS=epsilon(x)
    INTEGER(i_d) :: n
    REAL(r_2)     :: ap, del, summ

    if (x == 0.0) then
       dgser = 0.0
       RETURN
    end if
    ap   = a
    summ = 1.0_r_2/a
    del  = summ
    do n=1, ITMAX
       ap   = ap + 1.0_r_2
       del  = del*x/ap
       summ = summ + del
       if (abs(del) < abs(summ)*EPS) exit
    end do
    dgser=summ*exp(-x+a*log(x)-gammln(a))

  END FUNCTION dgser

  !**********************************************************************************************************************

  REAL(SP) ELEMENTAL PURE FUNCTION sigamma(a,x)

    IMPLICIT NONE

    REAL(SP), INTENT(IN) :: a, x
    REAL(SP) :: gln
    gln = gammln(a)
    if (x < a+1.0_sp) then
       sigamma = 1.0_sp - gser(a,x)
    else
       sigamma = gcf(a,x)
    end if
    sigamma = sigamma * exp(gln)

  END FUNCTION sigamma


  REAL(r_2) ELEMENTAL PURE FUNCTION digamma(a,x)
    USE cable_dimensions, ONLY: r_2
    IMPLICIT NONE
    REAL(r_2), INTENT(IN) :: a, x
    REAL(r_2) :: gln
    gln = gammln(a)
    if (x < a+1.0_r_2) then
       digamma = 1.0_r_2 - gser(a,x)
    else
       digamma = gcf(a,x)
    end if
    digamma = digamma * exp(gln)

  END FUNCTION digamma

  !**********************************************************************************************************************
  !MC! this routines has to be adjusted for csol in freezing, probably
  REAL(r_2) ELEMENTAL PURE FUNCTION phi(hr0, lambda, eta, phie, he, T, Ksat)

    USE cable_sli_numbers, ONLY: zero, one, gravity, lambdaf, Rgas, Tzero, csol

    IMPLICIT NONE

    REAL(r_2), INTENT(IN)           :: hr0, lambda, eta, phie, he, T
    REAL(r_2), INTENT(IN), OPTIONAL :: Ksat

    REAL(r_2) :: h, csol1

    !MC! freezing point? csol?
    csol1 = zero ! use zero instead of csol for the moment
    if (present(Ksat)) then
       if (T < zero) then     ! frozen soil  !! need to adjust freezing point for csol and use global csol
          h   = (lambdaf*T/gravity/(T+Tzero)) + csol1*Rgas* (T + Tzero)/gravity
          phi = Ksat * he/(one-eta*lambda) * (h/he)**(one-eta*lambda)
       else
          phi = zero
       endif
    else
       if (hr0 < one) then
          phi = phie * (Rgas *(T+Tzero) *log(hr0) /(gravity * he * Mw))**(one-eta*lambda)
       else
          phi = phie
       endif
    endif

  END FUNCTION phi

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION rh0_sol(hr0, solin)

    IMPLICIT NONE

    REAL(r_2),        INTENT(in) :: hr0
    TYPE(solve_type), INTENT(in) :: solin

    !    Rnet = Rnet + rhow*lambdaf*var(1)%iice*((phi(hr1,lambda,eta,phie,he,T1) &
    !         - phi(hr0,lambda,eta,phie,he,T1))/dz - K1)
    rh0_sol = (solin%Dv*(-hr0 + solin%hr1)*csat(solin%T1))/solin%dz - &
         (thousand*(solin%dz*solin%K1 + phi(hr0,solin%lambda,solin%eta,solin%phie,solin%he,solin%T1) &
         - phi(solin%hr1,solin%lambda,solin%eta,solin%phie,solin%he,solin%T1)))/solin%dz - &
         (solin%Dv*solin%hr1*(solin%cva*solin%gv*solin%lambdav + solin%Rnet + solin%gh*(-solin%T1 + solin%Ta) &
         - solin%gv*hr0*solin%lambdav*csat(solin%T1))*slope_csat(solin%T1)) / &
         (solin%DH + solin%dz*solin%gh + solin%dz*solin%gv*hr0*solin%lambdav*slope_csat(solin%T1)) + &
         ((solin%DH + solin%dz*solin%gh)*solin%gv*(solin%cva - hr0*csat(solin%T1)) &
         - solin%dz*solin%gv*hr0*(solin%Rnet + solin%gh*(-solin%T1 + solin%Ta))*slope_csat(solin%T1)) / &
         (solin%DH + solin%dz*solin%gh + solin%dz*solin%gv*hr0*solin%lambdav*slope_csat(solin%T1))

  END FUNCTION rh0_sol

  !**********************************************************************************************************************

  ! Using an elemental subroutine as a function argument is not Fortran90 standard.
  ! This would have been rh0_sol in function rtbis. intel's ifort and nag's f95 accept it but gfortran does not.
  ! Write rh0_sol into rtbis -> rtbis_rh0
  ! It does not check that f(x1) and f(x2) have different signs and not if iteration > MAXIT

  REAL(r_2) ELEMENTAL PURE FUNCTION rtbis_rh0(sol, x1, x2, xacc)

    IMPLICIT NONE

    TYPE(solve_type), INTENT(IN) :: sol
    REAL(r_2),        INTENT(IN) :: x1, x2, xacc

    INTEGER(i_d), PARAMETER :: MAXIT=40
    INTEGER(i_d) :: j
    REAL(r_2)    :: dx, f, fmid, xmid

    fmid = rh0_sol(x2,sol)
    f    = rh0_sol(x1,sol)
    if (f < 0.0_r_2) then
       rtbis_rh0 = x1
       dx        = x2-x1
    else
       rtbis_rh0 = x2
       dx        = x1-x2
    end if
    do j=1, MAXIT
       dx   = dx*0.5_r_2
       xmid = rtbis_rh0+dx
       fmid = rh0_sol(xmid,sol)
       if (fmid <= 0.0_r_2) rtbis_rh0 = xmid
       if (abs(dx) < xacc .or. fmid == 0.0_r_2) RETURN
    end do

  END FUNCTION rtbis_rh0

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION slope_csat(T)
    !returns slope of sat vapour pressure curve in kg m-3 K-1
    USE cable_sli_numbers, ONLY: Tzero, Rgas, Mw, esata, esatb, esatc

    IMPLICIT NONE

    real(r_2), intent(in) :: T
    real(r_2) :: csat

    csat       = esata * exp(esatb*T/(T+esatc)) * Mw/Rgas/(T+Tzero)
    slope_csat = csat * esatb*esatc/(T+esatc)**2

  END FUNCTION slope_csat

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION slope_esat(T)
    !returns slope of sat vapour pressure curve in Pa K^-1
    USE cable_sli_numbers, ONLY: esata, esatb, esatc

    IMPLICIT NONE

    real(r_2), intent(in) :: T 
    real(r_2) :: esat

    esat       = esata * exp(esatb*T/(T+esatc))
    slope_esat = esat * esatb*esatc/(T+esatc)**2

  END FUNCTION slope_esat

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION Sofh(h,parin)

    IMPLICIT NONE

    REAL(r_2),    INTENT(IN) :: h
    TYPE(params), INTENT(IN) :: parin
    
    ! Get saturation S from matric head h.
    ! Definitions of arguments:
    ! h   - matric head.

    Sofh = (h/parin%he)**(-parin%lam) ! Sofh not used much so ** not an issue

  END FUNCTION Sofh

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION Tfrz(S,he,b)
    ! determines freezing point temperature for a given moisture content
    USE cable_sli_numbers, ONLY: gravity, lambdaf, csol, Rgas, Tzero

    IMPLICIT NONE

    real(r_2), intent(in) :: S, he, b

    if (csol > e3) then
       Tfrz = (gravity*he - (S)**b* (lambdaf + 2.*csol*Rgas*Tzero) + &
            sqrt(gravity**2*he**2 - 2.*gravity*he*lambdaf*(S)**b + &
            lambdaf*(S)**(2.*b)*(lambdaf + 4.*csol*Rgas*Tzero)))/ &
            (2.*csol*Rgas*(S)**b)
    else
       Tfrz = (gravity*he*Tzero) / (-(gravity*he) + lambdaf*(S)**b)
    endif

  END FUNCTION Tfrz

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION thetalmax(Tin,S,he,b,thre,the)
    ! determines maximum liquid water content, given T
    USE cable_sli_numbers, ONLY: gravity, lambdaf, csol, Rgas, Tzero

    IMPLICIT NONE

    real(r_2), intent(in) :: Tin,S,he,b,thre,the
    real(r_2)             :: PI, psi, h, T

    T   = min(Tfrz(S,he,b),Tin)
    PI  = -csol *Rgas *(T+Tzero)/gravity ! osmotic potential (m)
    psi = lambdaf*T/(gravity*(T+Tzero)) ! matric potential in presence of ice
    h   = psi-PI       ! moisture potential in presence of ice

    thetalmax = thre*(h/he)**(-1/b) + (the-thre)

  END FUNCTION thetalmax

  !**********************************************************************************************************************

  REAL(r_2) ELEMENTAL PURE FUNCTION weight(parin, h, K, phi, dz)

    IMPLICIT NONE

    TYPE(params), INTENT(IN) :: parin
    REAL(r_2),    INTENT(IN) :: h
    REAL(r_2),    INTENT(IN) :: K
    REAL(r_2),    INTENT(IN) :: phi
    REAL(r_2),    INTENT(IN) :: dz
    ! Get conductivity weighting for gravity flux calculations.
    ! Definitions of arguments:
    ! l   - land point
    ! j   - soil type no.
    ! h   - matric head.
    ! K   - conductivity.
    ! phi - MFP.
    ! dz  - flow path length.
    LOGICAL   :: done
    REAL(r_2) :: a, hz, Khz, Kz, phiz, w, x

    done = .false.
    hz   = h-gf*dz ! gf is gravity fac in direction of dz
    if (h<parin%he) then
       a = parin%lam*parin%eta
       x = -gf*dz/h
       if (a <= 3.0_r_2 .or. x*(a-3.0_r_2) <= 4.0_r_2) then ! use predetermined approx.
          w    = (60.0_r_2+x*(70.0_r_2+10.0_r_2*a+x*(16.0_r_2+a*(5.0_r_2+a))))/ &
               (120.0_r_2+x*(120.0_r_2+x*(22.0_r_2+2.0_r_2*a**2)))
          done = .true.
       end if
    end if
    if (.not. done) then
       call hyofh(hz, parin, Kz, Khz, phiz) ! accurate but slower
       w = -((phiz-phi)/(gf*dz)+K)/(Kz-K)
    end if
    weight = min(max(w,zero),one)

  END FUNCTION weight

  !**********************************************************************************************************************

  FUNCTION zerovars()

    ! Sets all fields of type vars to zero

    IMPLICIT NONE
    
    TYPE(vars) :: zerovars
    
    zerovars%isat      = 0
    zerovars%h         = zero
    zerovars%phi       = zero
    zerovars%phiS      = zero
    zerovars%K         = zero
    zerovars%KS        = zero
    zerovars%Dv        = zero
    zerovars%cvsat     = zero
    zerovars%rh        = zero
    zerovars%phiv      = zero
    zerovars%phivS     = zero
    zerovars%kH        = zero
    zerovars%kE        = zero
    zerovars%kth       = zero
    zerovars%csoil     = zero
    zerovars%eta_th    = zero
    zerovars%hS        = zero
    zerovars%rhS       = zero
    zerovars%sl        = zero
    zerovars%cv        = zero
    zerovars%cvsatT    = zero
    zerovars%cvS       = zero
    zerovars%kv        = zero
    zerovars%iice      = 0
    zerovars%thetai    = zero
    zerovars%thetal    = zero
    zerovars%phiT      = zero
    zerovars%KT        = zero
    zerovars%lambdav   = zero
    zerovars%lambdaf   = zero
    zerovars%he        = zero
    zerovars%phie      = zero
    zerovars%Ksat      = zero
    zerovars%dthetaldT = zero
    zerovars%Tfrz      = zero
    zerovars%csoileff  = zero

  END FUNCTION zerovars

  !**********************************************************************************************************************

END MODULE cable_sli_utils
