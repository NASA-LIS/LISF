!
! Original code by P.J. Ross 2001-2007  See: Ross, P.J. (2003)
!     Modeling soil water and solute transport - fast, simplified numerical solutions. Agron. J. 95:1352-1361.
!
! Modified by V. Haverd 2008: matrix expanded by factor of 2 to allow solution of coupled heat/moisture fluxes
! Explicit calculation of heat/moisture fluxes at surface
! Switchable option for litter
! Isotope subroutine (V. Haverd and M. Cuntz, September 2008)
!
! Frozen Soil (V. Haverd Sep 2010 - Jan 2011)
! Pond lumped with top soil layer (V. Haverd Jan 2011)
! Snow included in layers -2:0 (V. Haverd Jan 2011)
! Convergence test for dthetaldT (V.Haverd Feb 2011)
! Include heat advection by liquid water flux (V.Haverd Feb 2011)

!
MODULE cable_sli_solve

  USE cable_dimensions, ONLY: r_2, i_d
  USE cable_sli_numbers,       ONLY: &
       zero, one, two, half, thousand, e1, e2, e3, e4, e5, e7, &  ! numbers
       Tzero, rlambda, lambdaf, Dva, rhocp, rhow, gf, hmin, & ! parameters
	   lambdas, csice , cswat, rhoi, &
       params, vars_aquifer, vars_met, vars, solve_type,  & ! types
       dSfac, h0min, Smax, dh0max, h0max, dSmax, dSmaxr, dtmax, dSmax, dSmaxr, & ! numerical limits
       dtmax, dtmin, dTsoilmax, dTLmax, nsteps_ice_max, tol_dthetaldT, &
       hbot, botbc ! boundary condition

  USE cable_sli_utils,         ONLY: &
       Sofh, hyofh, hyofS, litter_props, massman_sparse, tri, &
       aquifer_props, Tfrz, csoil, thetalmax, dthetalmaxdTh, &
	   dthetalmaxdT, &
       getfluxes_vp, getheatfluxes, flux, sol, rtbis_rh0, phi, &
       csat, slope_csat, potential_evap, isosub, setsol, zerovars
!       csat, slope_csat, potential_evap, tri, isosub, setsol, zerovars
  USE cable_sli_roots,         ONLY: getrex

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: solve ! solution routine

  INTEGER(i_d), DIMENSION(:), ALLOCATABLE :: nless, n_noconverge ! global counters

  ! Definitions of public entities and private parameters (see above for default
  ! values):
  ! botbc    - bottom boundary condn for water; "constant head", "free drainage",
  !    "seepage", or "zero flux". Constant head means that matric head h
  !    is specified. Free drainage means zero gradient of matric head,
  !    i.e. unit hydraulic gradient. Seepage means zero flux when the
  !    matric head is below zero and an upper limit of zero for the head.
  ! h0max    - max pond depth allowed before runoff.
  ! hbot    - matric head at bottom of profile when botbc set to "constant head".
  ! dSmax    - max change in S (the "effective saturation") of any unsaturated
  !    layer to aim for each time step; controls time step size.
  ! dSmaxr   - maximum negative relative change in S each time step. This
  !    parameter helps avoid very small or negative S.
  ! dtmax    - max time step allowed.
  ! dsmmax   - max solute change per time step (see dSmax); user should set this
  !    according to solute units used. Units for different solutes can be
  !    scaled by the user (e.g. to an expected max of around 1.0).
  ! dSfac    - a change in S of up to dSfac*dSmax is accepted.
  ! dpmaxr   - relative change in matric flux potential (MFP) phi that is
  !    accepted for convergence when finding head h at soil interfaces.
  ! h0min    - min (negative) value for surface pond when it empties.
  ! Smax    - max value for layer saturation to allow some overshoot.
  ! dh0max   - allowable overshoot when pond reaches max allowed depth.
  ! solve    - sub to call to solve RE
  !

CONTAINS

  !**********************************************************************************************************************

  SUBROUTINE solve(ts, tfin, irec, mp, qprec, n, nsol, dx, h0, S,thetai,Jsensible,Tsoil, evap, evap_pot, runoff, &
       infil, drainage, discharge, qh, nsteps, vmet, vlit,csoil, T0, rh0, Tsurface, rhsurface, Hcum, lEcum, &
       Gcum,Qadvcum,Jcol_sensible,Jcol_latent_S,Jcol_latent_T, deltaice_cum_T, deltaice_cum_S, dxL, zdelta, &
       SL, TL, plit, par, Etrans, fwsoil, wex, heads, FS, gamma, &
       ciso, cisos, ciso0, cisoL, cprec, cali, qiso_in, qiso_out, qiso_evap_cum, qiso_trans_cum, qiso_liq_adv, &
       qiso_vap_adv, qiso_liq_diff, qiso_vap_diff, qvsig, qlsig, qvTsig, qvh, deltaTa, lE_old, TL_test, &
       dolitter, doisotopologue, dotestcase, dosepts, docondition, doadvection)

    IMPLICIT NONE

    REAL(r_2),                             INTENT(IN)              :: ts, tfin
    INTEGER(i_d),                          INTENT(IN)              :: irec, mp
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN)              :: qprec
    INTEGER(i_d),                          INTENT(IN)              :: n, nsol
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(IN)              :: dx
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: h0
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT)           :: S
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT)             :: thetai
	REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT)             :: Jsensible
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT)           :: Tsoil
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: evap, evap_pot, runoff, infil
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: drainage, discharge
    REAL(r_2),      DIMENSION(1:mp,0:n),   INTENT(OUT)             :: qh
    INTEGER(i_d),   DIMENSION(1:mp),       INTENT(INOUT)           :: nsteps
    TYPE(vars_met), DIMENSION(1:mp),       INTENT(INOUT)           :: vmet
    TYPE(vars),     DIMENSION(1:mp),       INTENT(INOUT)           :: vlit
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: T0, rh0, Tsurface, rhsurface
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: Hcum, lEcum, Gcum,Qadvcum, deltaice_cum_T, deltaice_cum_S
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT)             :: Jcol_sensible, Jcol_latent_S, Jcol_latent_T
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT)             :: csoil
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN)              :: dxL
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT)           :: zdelta, SL, Tl
    TYPE(params),   DIMENSION(1:mp),       INTENT(IN)              :: plit
    TYPE(params),   DIMENSION(1:mp,1:n),   INTENT(IN)              :: par
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN)              :: Etrans
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT),   OPTIONAL :: fwsoil
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: gamma
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(INOUT), OPTIONAL :: wex
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT),   OPTIONAL :: heads
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(IN),    OPTIONAL :: FS
    REAL(r_2),      DIMENSION(1:mp,0:n),   INTENT(INOUT), OPTIONAL :: ciso
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT), OPTIONAL :: cisos, ciso0, cisoL
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: cprec, cali
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT),   OPTIONAL :: qiso_in, qiso_out
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT),   OPTIONAL :: qiso_evap_cum, qiso_trans_cum
    REAL(r_2),      DIMENSION(1:mp,1:n),   INTENT(OUT),   OPTIONAL :: qiso_liq_adv, qiso_vap_adv
    REAL(r_2),      DIMENSION(1:mp,1:n-1), INTENT(OUT),   OPTIONAL :: qiso_liq_diff, qiso_vap_diff
    REAL(r_2),      DIMENSION(1:mp,0:n),   INTENT(OUT),   OPTIONAL :: qvsig, qlsig, qvTsig, qvh
    REAL(r_2),      DIMENSION(1:mp),       INTENT(INOUT), OPTIONAL :: deltaTa
    REAL(r_2),      DIMENSION(1:mp),       INTENT(IN),    OPTIONAL :: lE_old
    REAL(r_2),      DIMENSION(1:mp),       INTENT(OUT),   OPTIONAL :: TL_test
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: dolitter       ! 0: no; 1: normal; 2: resistance
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: doisotopologue ! 0: no isotope; 1: HDO; 2: H218O
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: dotestcase     ! isotopic tests
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: dosepts        ! 0: normal; 1: uncouple T & S
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: docondition    ! 0: no cond., 1: columns, 2: lines, 3: both
    INTEGER(i_d),                          INTENT(IN),    OPTIONAL :: doadvection       ! 0: off; 1: onn

    ! Solves the RE and, optionally, the ADE from time ts to tfin.
    ! Definitions of arguments:
    ! Required args:
    ! ts   - start time (h).
    ! tfin   - finish time.
    ! qprec   - precipitation (or water input) rate (fluxes are in cm/h).
    ! qevap   - potl evaporation rate from soil surface.
    ! n    - no. of soil layers.
    ! nsol   - no. of solutes.
    ! dx(1:n) - layer thicknesses.
    ! h0   - surface head, equal to depth of surface pond.
    ! S(1:n)  - degree of saturation ("effective satn") of layers.
    ! evap   - cumulative evaporation from soil surface (cm, not initialised).
    ! runoff  - cumulative runoff.
    ! infil   - cumulative net infiltration (time integral of flux across surface).
    ! drn   - cumulative net drainage (time integral of flux across bottom).
    ! nsteps  - cumulative no. of time steps for RE soln.
    ! Optional args:
    ! heads(1:n)   - matric heads h of layers at finish.
    ! qexsub    - subroutine to get layer water extraction rates (cm/h) by
    !     plants. Note that there is no solute extraction and osmotic
    !     effects due to solute are ignored. Arguments:
    !     qex(1:n) - layer extraction rates; qexh(1:n) - partial
    !     derivs of qex wrt h.
    ! wex(1:n)    - cumulative water extraction from layers.
    ! cin(1:nsol)   - solute concns in water input (user's units/cc).
    ! c0(1:nsol)   - solute concns in surface pond.
    ! sm(1:n,1:nsol)  - solute (mass) concns in layers.
    ! soff(1:nsol)   - cumulative solute runoff (user's units).
    ! sinfil(1:nsol)  - cumulative solute infiltration.
    ! sdrn(1:nsol)   - cumulative solute drainage.
    ! nssteps(1:nsol) - cumulative no. of time steps for ADE soln.
    ! isosub    - subroutine to get adsorbed solute (units/g soil) from concn
    !     in soil water according to chosen isotherm code.
    !     Arguments: iso - 2 character code; c - concn in soil water;
    !     p(:) - isotherm parameters; f - adsorbed mass/g soil;
    !     fc - deriv of f wrt c (slope of isotherm curve). Note that
    !     linear adsorption does not require a sub, and other types
    !     are available in sub isosub.

    REAL(r_2),    DIMENSION(1:mp)       :: Ebal0
    REAL(r_2),    DIMENSION(1:mp)       :: precip, qevap
    REAL(r_2),    DIMENSION(1:mp)       :: fws
    REAL(r_2),    DIMENSION(1:mp)       :: qL, qhL, qybL, qTbL, qhTbL, qhybL, rexcol, wcol, ql0, qv0, qvT0, tol
    LOGICAL,      DIMENSION(1:mp)       :: again, getq0,getqn,init, maxpond
    LOGICAL,      DIMENSION(1:mp)       :: lns
	LOGICAL,      DIMENSION(1:mp,1:n)       :: again_ice
    INTEGER(i_d), DIMENSION(1:mp)       :: ih0, iok, itmp, ns, nsat, nsatlast, nsteps0
    REAL(r_2),    DIMENSION(1:mp)       :: accel, dmax, dt, dwinfil, dwoff, fac, Khmin1, Kmin1, phimin1, phip
    REAL(r_2),    DIMENSION(1:mp)       :: qpme, rsig, rsigdt, sig, t
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: Sbot, Tbot
    REAL(r_2),    DIMENSION(1:mp,1:n-1) :: dz
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: hint, phimin, qex, qexd
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: aa, bb, cc, dd, dy, ee, ff, gg, q, qya, qyb, qTa, qTb
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: aah, bbh, cch, ddh, eeh, ffh, ggh, dTsoil, qhya, qhyb, qhTa, qhTb
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: qadv, qadvya, qadvyb, qadvTa, qadvTb
    TYPE(vars)                          :: vtmp
    TYPE(vars),   DIMENSION(1:mp,1:n)   :: var
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: qsig, qhsig, qadvsig
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: qliq, qv, qvT, qlya, qlyb, qvya, qvyb
    TYPE(vars),   DIMENSION(1:mp,1:n)   :: vcall
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: deltaJ, deltaJ1, deltaS, Ssig
    REAL(r_2),    DIMENSION(1:mp,0:n)   :: tmp2d1, tmp2d2
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: S0, Sliq0, Sliq, deltaSliq, cv0, deltacv
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: tmp_thetasat, tmp_thetar, tmp_tortuosity
    REAL(r_2),    DIMENSION(1:mp,1:n)   :: delthetai, delJsensible
	REAL(r_2),    DIMENSION(1:mp)       :: delthetai_col, delJsensible_col
    INTEGER(i_d), DIMENSION(1:mp,1:n)   :: isave, nsteps_ice
    TYPE(vars),         DIMENSION(1:mp) :: vtop, vbot
    TYPE(vars_aquifer), DIMENSION(1:mp) :: v_aquifer
    REAL(r_2),          DIMENSION(1:mp) :: qd, dwcol, dwdrainage, drn,inlit, dwinlit, drexcol, dwdischarge
    REAL(r_2),          DIMENSION(1:mp) :: dJcol, Jcol,dicecol_S,dicecol_T
    REAL(r_2),          DIMENSION(1:mp) :: dJcol_latent_S, dJcol_latent_T, dJcol_sensible ,dT0dTsoil,w
    REAL(r_2),          DIMENSION(1:mp,1:n):: deltaJ_latent_S, deltaJ_latent_T, deltaJ_sensible_S, deltaJ_sensible_T
    REAL(r_2),          DIMENSION(1:mp) :: qevapsig
    REAL(r_2),          DIMENSION(1:mp) :: qrunoff, deltatheta, h0sig
    REAL(r_2),          DIMENSION(1:mp) :: tmp1d1, tmp1d2, tmp1d3,  tmp1d4
    REAL(r_2),          DIMENSION(1:mp) :: ponding, deltah0
    REAL(r_2),          DIMENSION(1:mp) :: Tsurface_pot, Epot, Hpot, Gpot, dEdrha, dEdTa, dEdTsoil, dGdTa, dGdTsoil
    REAL(r_2),          DIMENSION(1:mp) :: SL0, deltaSL, cvL0, SLliq0, deltacvL, SLliq, deltaSLliq
    REAL(r_2),          DIMENSION(1:mp) :: qiso_evap, qiso_trans
    REAL(r_2),          DIMENSION(1:mp) :: H, lE0, G0
    REAL(r_2),          DIMENSION(1:mp) :: T0_test, rh0_test, qevap_test, G0_test
    REAL(r_2),          DIMENSION(1:mp) :: Tfreezing
    REAL(r_2),          DIMENSION(1:mp) :: dtdT
    REAL(r_2),          DIMENSION(1:mp) :: e41d, one1d
    REAL(r_2),          DIMENSION(1:mp,0:n) :: LHS, RHS,LHS_h, RHS_h
    INTEGER(i_d),       DIMENSION(1:mp) :: pondcase
    INTEGER(i_d),       DIMENSION(1:mp) :: nns, iflux
    LOGICAL      :: litter
    INTEGER(i_d) :: i, j, k, kk, condition
    INTEGER(i_d) :: littercase, isotopologue,advection, testcase, septs ! switches
    REAL(r_2)    :: ztmp, tmp1, tmp2, hr0, TT0, Qqevap, LlE0, GG0, HH0
    REAL(r_2)    :: dTsoil_test, dS_test

    !open (unit=7, file="Test.out", status="replace", position="rewind")
    ! The derived types params and vars hold soil water parameters and variables.
    ! Parameter names often end in e, which loosely denotes "air entry", i.e.,
    ! values at h=he. While values of water content th and hydraulic conductivity K
    ! at h=he are equal to those at saturation, the derivs wrt S are nonzero. The
    ! MFP phi for h>he is given by phi=phie+Ke*(h-he). The saturation status of a
    ! layer is stored as 0 or 1 in isat since S may be >1 (because of previous
    ! overshoot) when a layer desaturates. Fluxes at the beginning of a time step
    ! and their partial derivs wrt S or phi of upper and lower layers or boundaries
    ! are stored in q, qya and qyb.

    ! set switches
    if (present(dolitter)) then
       littercase = dolitter
    else
       littercase = 0
    endif
    if (present(doadvection)) then
       advection = doadvection
    else
       advection = 0
    endif
    if (littercase > 2) then
       write(2,*) 'dolitter not in [0-2]: ', littercase
       stop
    endif

    if (present(doisotopologue)) then
       isotopologue = doisotopologue
    else
       isotopologue = 0
    endif
    if (isotopologue > 2) then
       write(2,*) 'doisotopologue not in [0-2]: ', isotopologue
       stop
    endif
    if (isotopologue /= 0 .and. (.not. present(ciso))) then
       write(2,*) 'doisotopologue /= 0 but no ciso present.'
       stop
    endif

    if (present(dotestcase)) then
       testcase = dotestcase
    else
       testcase = 0
    endif
    if (testcase > 8) then
       write(2,*) 'dotestcase not in [0-8]: ', testcase
       stop
    endif

    if (present(dosepts)) then
       septs = dosepts
    else
       septs = 0
    endif
    if (septs > 1) then
       write(2,*) 'dosepts not in [0-1]: ', septs
       stop
    endif

    if (present(docondition)) then
       condition = docondition
    else
       condition = 0
    endif
    if (condition < 0 .or. condition > 3) then
       write(2,*) 'docondition not in [0-3]: ', condition
       stop
    endif

    ! global counters
    if (.not. allocated(nless)) allocate(nless(mp))
    nless(:) = 0
    if (.not. allocated(n_noconverge)) allocate(n_noconverge(mp))
    n_noconverge(:) = 0

    ! set solve_type for numerical derivatives
    if (.not. allocated(sol)) call setsol(mp)

    ! initialise cumulative variables
    wcol(:)      = zero
    Jcol(:)      = zero
    Jcol_sensible(:) = zero
    Jcol_latent_S(:) = zero
    Jcol_latent_T(:) = zero
    deltaice_cum_T(:) = zero
    deltaice_cum_S(:) = zero
    deltaJ_sensible_S(:,:) = zero
    deltaJ_sensible_T(:,:) = zero
    deltaJ_latent_S(:,:) = zero
    deltaJ_latent_T(:,:) = zero
    drainage(:)  = zero
    discharge(:) = zero
    infil(:)     = zero
    inlit(:)     = zero
    evap(:)      = zero
    evap_pot(:)  = zero
    runoff(:)    = zero
    rexcol(:)    = zero
    ponding(:)   = zero
    Hcum(:)      = zero
    Gcum(:)      = zero
    lEcum(:)     = zero
    Qadvcum(:)  = zero
    wex(:,:)          = zero
    precip(:)         = zero
    drn(:)            = zero
    if (isotopologue /= 0) then
       qiso_evap_cum(:)  = zero
       qiso_trans_cum(:) = zero
    endif
    maxpond(:) = .false.
    deltah0(:) = zero
    ! zero var-structure that contains all the hydrological variables
    vtmp = zerovars()
    vtmp%h       = one
    vtmp%lambdav = rlambda
    vtmp%lambdaf = lambdaf
    var  = spread(spread(vtmp,1,mp),2,n)
    ! zero vars at the bottom and top of the soil column
    vtop = spread(vtmp,1,mp)
    vbot = spread(vtmp,1,mp)
    vlit = spread(vtmp,1,mp)
    hint(:,:)   = zero
    phimin(:,:) = zero
    q(:,:)      = zero
    qya(:,:)    = zero
    qyb(:,:)    = zero
    qTa(:,:)    = zero
    qTb(:,:)    = zero
    qhya(:,:)   = zero
    qhyb(:,:)   = zero
    qhTa(:,:)   = zero
    qhTb(:,:)   = zero
    aa(:,:)     = zero
    aah(:,:)    = zero
    bb(:,:)     = zero
    bbh(:,:)    = zero
    cc(:,:)     = zero
    cch(:,:)    = zero
    dd(:,:)     = zero
    ddh(:,:)    = zero
    ee(:,:)     = zero
    eeh(:,:)    = zero
    ff(:,:)     = zero
    ffh(:,:)    = zero
    gg(:,:)     = zero
    ggh(:,:)    = zero
    dy(:,:)     = zero
    dTsoil(:,:) = zero

    e41d(:)  = e4
    one1d(:) = one

    litter = .false.
    if (littercase == 1) litter=.true. ! full litter model

    qex(:,:)    = zero
    qexd(:,:)   = zero
    fwsoil(:)   = zero

    phip(:) = max(par(:,1)%phie-par(:,1)%he*par(:,1)%Ke, 1.00001_r_2*par(:,1)%phie) ! phi at h=0

    ! get K, Kh and phi at hmin (hmin is smallest h, stored in hy-props)
    call hyofh(spread(hmin,1,mp), par(:,1), Kmin1(:), Khmin1(:), phimin1(:))

    dz(:,:) = half*(dx(:,1:n-1)+dx(:,2:n)) ! flow paths

    !----- set up for boundary conditions
    getq0(:) = .false.
    getqn(:) = .false.
    if (botbc == "constant head") then ! h at bottom bdry specified
       getqn(:)  = .true.
       tmp1d1(:)  = hbot
       ! for hbot < he
       Sbot(:,:) = spread(Sofh(tmp1d1,par(:,n)),1,n)
       !MC! strange in original code: all calc with bottom n except Tsoil(1)
       Tbot(:,:) = spread(Tsoil(:,1),1,n)
       !Tbot(:,:) = spread(Tsoil(:,n),1,n)
       call hyofS(Sbot, Tbot, par, vcall)
       ! for hbot >= he
       vtmp = zerovars()
       vtmp%isat    = 1
       vtmp%h       = hbot
       vtmp%rh      = one
       vtmp%lambdav = rlambda
       vtmp%lambdaf = lambdaf
       vbot = spread(vtmp,1,mp)
       vbot(:)%phi = (hbot-par(:,n)%he)*par(:,n)%Ke+par(:,n)%phie
       vbot(:)%K   = par(:,n)%Ke
       where (par(:,n)%he > hbot)
          vbot(:)      = vcall(:,n)
          vbot(:)%isat = 0
       endwhere
    end if
    !----- end set up for boundary conditions

    !----- initialise
    t(:)       = ts
    nsteps0(:) = nsteps
    nsat(:)    = 0
    qd(:)      = zero
    ! initialise saturated regions
    var(:,:)%isat  = 0
    where (S(:,:) >= one)
       var(:,:)%phi  = par(:,:)%phie
       var(:,:)%K    = par(:,:)%Ke
       var(:,:)%isat = 1
    end where
    vlit(:)%isat = 0
    where (SL(:) >= one) vlit(:)%isat = 1

    ! initialise acquifer
    v_aquifer(:)%zsoil  = sum(dx(:,:),2)
    v_aquifer(:)%zdelta = zdelta(:)
    call aquifer_props(v_aquifer(:))

    ! initialise litter
    if (littercase == 1 .or. littercase == 2) then
       call litter_props(Sl(:), Tl(:), vlit(:), plit(:), h0(:))
    endif
    ! Add resistance through litter for simple litter model
    if (littercase == 2) then
       vmet(:)%rbw = vmet(:)%rbw + dxL(:)/vlit(:)%Dv
       ztmp        = one/rhocp
       vmet(:)%rbh = vmet(:)%rbh + dxL(:)/(vlit(:)%kth*ztmp)
    endif

    lE0(:) = lE_old(:) ! used for initial guess of litter temperature
    !----- end initialise

    !----- solve until tfin
    init(:) = .true. ! flag to initialise h at soil interfaces
    do kk=1, mp
	   call getrex(S(kk,:), qex(kk,:), fws(kk), FS(kk,:), par(kk,:)%the, par(kk,:)%thw, Etrans(kk), gamma(kk))
       do while (t(kk) < tfin)

          !----- take next time step
          iflux(kk)=1
          again(kk)  = .true. ! flag for recalcn of fluxes (default=false)
          do while (again(kk)) ! sometimes need twice to adjust phi at satn

             nsatlast(kk) = nsat(kk) ! for detecting onset of profile saturation
             !nsat(kk)     = sum(var(:,:)%isat,2) ! no. of sat layers
             nsat(kk)     = sum(var(kk,:)%isat,1) ! no. of sat layers
             sig(kk)      = half
             !where (nsat(kk) /= 0) sig(kk) = one ! time weighting sigma
             if (nsat(kk) /= 0) sig(kk) = one ! time weighting sigma
             rsig(kk)     = one/sig(kk)

             ! update variables
             if (iflux(kk)==1) then
                ! Calc flux matric potentials (and derivatives) from S
                ! this set var-structure
                isave(kk,:) = var(kk,:)%isat
                call hyofS(S(kk,:), Tsoil(kk,:), par(kk,:), var(kk,:)) ! for layers where S<1
                ! phi at h=0 (special for frozen soil)
                !MC! strange in original: should be var(1) but is var(n)
                ! v is pointing to var(n) after first loop; to var(1) in first loop
                phip(kk) = max(var(kk,n)%phie-var(kk,n)%he*var(kk,n)%Ksat, (one+e5)*var(kk,n)%phie)
                !phip(kk) = max(var(kk,1)%phie-var(kk,1)%he*var(kk,1)%Ksat, (one+e5)*var(kk,1)%phie)
                var(kk,:)%isat = isave(kk,:)
                thetai(kk,:)     = var(kk,:)%thetai
				Jsensible(kk,1) = (var(kk,1)%csoil* dx(kk,1)+h0(kk)*cswat)*(Tsoil(kk,1)+Tzero)	 
				Jsensible(kk,2:n) = var(kk,2:n)%csoil*(Tsoil(kk,2:n)+Tzero)* dx(kk,2:n)
             endif

             ! initialise litter vars
             if (iflux(kk)==1 .and. (littercase==1 .or. littercase == 2)) then
                call litter_props(Sl(kk), Tl(kk), vlit(kk), plit(kk), h0(kk))
             endif

             ! phi is solution var at satn, so h calc from phi where S>=1 - done in hyofS above for S<1
             !MC! one or the other
             where (S(kk,:) >= one) &
                  var(kk,:)%h = par(kk,:)%he + (var(kk,:)%phi-par(kk,:)%phie)/par(kk,:)%Ke
             where (S(kk,:) >= one) & ! (special for frozen soil)
                  var(kk,:)%h = var(kk,:)%he + (var(kk,:)%phi-var(kk,:)%phie)/var(kk,:)%Ksat
             !MC! one or the other

             if (botbc=="aquifer") then
                !where (abs(v_aquifer(kk)%zdelta-v_aquifer(kk)%zsoil) < e2)
                !   v_aquifer(kk)%K = var(kk,n)%K
                !elsewhere
                !   v_aquifer(kk)%K = var(kk,n)%K * &
                !        (one - exp(-v_aquifer(kk)%f*(v_aquifer(kk)%zdelta-v_aquifer(kk)%zsoil))) / &
                !        (v_aquifer(kk)%f*(v_aquifer(kk)%zdelta-v_aquifer(kk)%zsoil))
                !endwhere
                if (abs(v_aquifer(kk)%zdelta-v_aquifer(kk)%zsoil) < e2) then
                   v_aquifer(kk)%K = var(kk,n)%K
                else
                   v_aquifer(kk)%K = var(kk,n)%K * &
                        (one - exp(-v_aquifer(kk)%f*(v_aquifer(kk)%zdelta-v_aquifer(kk)%zsoil))) / &
                        (v_aquifer(kk)%f*(v_aquifer(kk)%zdelta-v_aquifer(kk)%zsoil))
                endif
             endif

             !----- get fluxes and derivs
             ! get surface condition
             ! ns==1 if no pond or full pond, i.e. do not solve for change in pond height
             ! ns==0 then change for pond height
             if ((var(kk,1)%phi <= phip(kk) .and. h0(kk) <= zero .and. nsat(kk) < n).or.(var(kk,1)%isat.eq.0)) then ! no ponding
                ns(kk)    = 1 ! start index for eqns
                lns(kk) = ns(kk)==1 
             else ! ponding
                ns(kk)    = 0
                lns(kk) = ns(kk)==1 
                !getq0(kk) = .true.
                TL(kk)    = vmet(kk)%Ta ! initialise pond T
                vtop(kk)%isat    = 1
                vtop(kk)%h       = h0(kk)
                vtop(kk)%phi     = (h0(kk)-var(kk,1)%he)*var(kk,1)%Ksat+var(kk,1)%phie 
                vtop(kk)%K       = var(kk,1)%Ksat
                vtop(kk)%rh      = one
                vtop(kk)%lambdav = rlambda
                vtop(kk)%lambdaf = lambdaf
                !call flux(jt(1), 0, ns, vtop, var(1), half*dx(1), q(0), qya(0), qyb(0),qTa(0), qTb(0))
                call flux(par(kk,1), lns(kk), vtop(kk), var(kk,1), half*dx(kk,1), &
                     q(kk,0), qya(kk,0), qyb(kk,0), qTa(kk,0), qTb(kk,0))
             endif

             ! saturated litter
             !where (Sl(kk) >= one)
             !   vlit(kk)%isat = 1
             !   vlit(kk)%h    = plit(kk)%he
             !endwhere
             if (Sl(kk) >= one) then
                vlit(kk)%isat = 1
                vlit(kk)%h    = plit(kk)%he
             endif

             ! get bottom boundary condn
             if (botbc=="seepage") then
                vtmp = zerovars()
                vbot = spread(vtmp,1,mp)
                if (var(kk,n)%h > -half*gf*dx(kk,n)) then
                   getqn(kk) = .true.
                   vbot(kk)%isat    = 1                 
                   vbot(kk)%phi     = (zero-var(kk,n)%he)*var(kk,n)%Ksat+var(kk,n)%phie ! (special for frozen soil)
                   vbot(kk)%K       = var(kk,n)%Ksat
                   vbot(kk)%rh      = one
                   vbot(kk)%lambdav = rlambda
                   vbot(kk)%lambdaf = lambdaf
                else
                   getqn(kk) = .false.
                endif
             end if



             ! Fluxes and derivatives at the surfaces (air/litter, litter/soil, and similar)
             ! three cases
10           pondcase(kk) = 0
             !where (ns(kk)==1 .and. (.not. litter) .and. (.not. maxpond(kk))) pondcase(kk) = 1 ! no litter and no pond
             !where (ns(kk)==1 .and.        litter  .and. (.not. maxpond(kk))) pondcase(kk) = 2 ! litter and no pond
             !where (ns(kk)==0                      .or.         maxpond(kk) ) pondcase(kk) = 3 ! pond with and without litter
             if (ns(kk)==1 .and. (.not. litter) .and. (.not. maxpond(kk))) pondcase(kk) = 1 ! no litter and no pond
             if (ns(kk)==1 .and.        litter  .and. (.not. maxpond(kk))) pondcase(kk) = 2 ! litter and no pond
             if (ns(kk)==0                      .or.         maxpond(kk) ) pondcase(kk) = 3 ! pond with and without litter
             sol(kk)%Ks = par(kk,1)%Ke
             select case (pondcase(kk))
             case (1) ! no litter and no pond
                ! structure that can be seen by rh0_sol_0d etc. which is used in rtbis (bisection)
                ! Look in paper Appendix A.6: upper boundary condition: no ponding or litter
                ! it just equates fluxes towards and away from the interface; 1 for water, 1 for heat
                ! have to solve numerically because of non-linearity in phi(h)
                sol(kk)%T1      = Tsoil(kk,1)
                sol(kk)%Ta      = vmet(kk)%Ta
                sol(kk)%cva     = vmet(kk)%cva*thousand
                sol(kk)%Rnet    = vmet(kk)%Rn
                sol(kk)%hr1     = var(kk,1)%rh
                sol(kk)%hra     = vmet(kk)%rha
                sol(kk)%Dv      = var(kk,1)%Dv
                sol(kk)%gv      = one/vmet(kk)%rbw
                sol(kk)%gh      = one/vmet(kk)%rrc * rhocp
                sol(kk)%Dh      = var(kk,1)%kth
                sol(kk)%dz      = half*dx(kk,1)
                sol(kk)%phie    = par(kk,1)%phie
                sol(kk)%he      = par(kk,1)%he
                sol(kk)%K1      = var(kk,1)%K
                sol(kk)%eta     = par(kk,1)%eta
                sol(kk)%lambda  = par(kk,1)%lam
                sol(kk)%lambdav = var(kk,1)%lambdav

                ! make sure there is a solution
                ! get rh_sol at upper and lower bounds of humidity
                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e4
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e3,e7)
                ! there is a zero in between
                if ((tmp1*tmp2) < zero) then
                   rh0(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else ! tests in Mathematica showed that if no zero between tmp1 and tmp2 then rh=1.0
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (01). Stop."
                      stop
                   endif
                endif

                ! Get temperature at interface
                ! bisection gave you rh0 as a return, now calc T0 as well
                !    Rnet = Rnet + rhow*lambdaf*var(1)%iice*((phi(hr1,lambda,eta,phie,he,T1) &
                !           - phi(hr0,lambda,eta,phie,he,T1))/Dz - K1)
                !T0 = T0_sol(rh0)
                hr0 = rh0(kk)
                T0(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                rhsurface(kk) = rh0(kk)

                ! get all interfacial fluxes
                ! plug the solution (rho, To) into the conservation equations
                hr0 = rh0(kk)
                TT0 = T0(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                LlE0   = Qqevap*thousand*sol(kk)%lambdav
                GG0    = sol(kk)%Dh/sol(kk)%dz*(TT0-sol(kk)%T1)
                HH0    = sol(kk)%gh*(TT0-sol(kk)%Ta)
                qevap(kk) = Qqevap
                lE0(kk)   = LlE0
                G0(kk)    = GG0
                H(kk)     = HH0

                Tsurface(kk) = T0(kk)

                qpme(kk)   = qprec(kk) - qevap(kk)
                q(kk,0)  = qprec(kk) - qevap(kk)
                qh(kk,0) = G0(kk)

                ! q refers to moisture in numerator
                ! qh refers to heat in numerator
                ! y refers to moisture in denominator
                ! T refers to temperature in denominator
                ! a refers to the layer above
                ! b refers to the layer below
                ! L flux into litter from above
                ! e.g. qhya = d(qH_k)/dS_k-1
                !      qhyb = d(qH_k)/dS_k
                !      qTa  = d(q_k)/dT_k-1
                !      qTb  = d(q_k)/dT_k
                !      qTbL = d(qL)/dTL
                !      qybL = d(qL)/dSL

                ! initialise derivatives to zero
                qya(kk,0)  = zero
                qTa(kk,0)  = zero
                qhya(kk,0) = zero
                qhyb(kk,0) = zero
                qhTa(kk,0) = zero
                dT0dTsoil(kk) = zero
                Ebal0(kk)  = vmet(kk)%Rn - lE0(kk) - H(kk) - G0(kk)

                qybL(kk) = zero
                qL(kk)   = zero
                qTbL(kk) = zero

                ! liquid and vapour fluxes
                qliq(kk,0) = zero
                qlya(kk,0) = zero
                qlyb(kk,0) = zero
                qv(kk,0)   = q(kk,0)
                qvya(kk,0) = qya(kk,0)
                qvyb(kk,0) = qyb(kk,0)

                ! perturb S and T to calculate partial derivs qyb(0), qTb(0), qhTb(0)
                dTsoil_test = e1
                dS_test     = e2
                sol(kk)%hr1  = min((var(kk,1)%rh + var(kk,1)%rhS*dS_test), one)
                sol(kk)%K1   = min((var(kk,1)%K + var(kk,1)%KS*dS_test), par(kk,1)%Ke)

                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e4
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e3,e7)
                if ((tmp1*tmp2) < zero) then
                   rh0_test(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0_test(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (02). Stop."
                      stop
                   endif
                endif
                !T0_test = T0_sol(rh0_test)
                hr0 = rh0_test(kk)
                T0_test(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))

                hr0 = rh0_test(kk)
                TT0 = T0_test(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                qevap_test(kk) = Qqevap

                ztmp = one/dS_test
                qyb(kk,0) = -(qevap_test(kk)-qevap(kk))*ztmp

                sol(kk)%T1  = Tsoil(kk,1) + dTsoil_test
                sol(kk)%hr1 = var(kk,1)%rh
                sol(kk)%K1  = var(kk,1)%K
                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e4
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e3,e7)
                if ((tmp1*tmp2) < zero) then
                   rh0_test(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0_test(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (03). Stop."
                      stop
                   endif
                endif
                !T0_test = T0_sol(rh0_test)
                hr0 = rh0_test(kk)
                T0_test(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))

                hr0 = rh0_test(kk)
                TT0 = T0_test(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                GG0    = sol(kk)%Dh/sol(kk)%dz*(TT0-sol(kk)%T1)
                qevap_test(kk) = Qqevap
                G0_test(kk) = GG0

                ztmp = one/dTsoil_test
                qTb(kk,0)  = -(qevap_test(kk)-qevap(kk))*ztmp
                qhTb(kk,0) = (G0_test(kk)-G0(kk))*ztmp
                dT0dTsoil(kk) = (T0_test(kk)-T0(kk))*ztmp
                ! end of partial derivative evaluation

             case (2) ! litter and no pond
                ! irec=1 then first external time step
                ! At first two iterations get an estimate of litter temperature
                !   based of current met but previous surface fluxes
                if (n_noconverge(kk) < 3 .and. irec > 1) then
                   !TL_test = (vmet%Rn - lE_old + vmet%Ta*(one/vmet%rrc * rhocp) + &
                   TL_test(kk) = (vmet(kk)%Rn - lE0(kk) + vmet(kk)%Ta*(one/vmet(kk)%rrc * rhocp) + &
                        Tsoil(kk,1)*(1./(half*dxL(kk)/vlit(kk)%kth + half*dx(kk,1)/var(kk,1)%kth)))/ &
                        (one/(half*dxL(kk)/vlit(kk)%kth + half*dx(kk,1)/var(kk,1)%kth) + one/vmet(kk)%rrc * rhocp)
                   deltaTa(kk) = TL_test(kk)-TL(kk)
                else
                   deltaTa(kk) = zero ! take normal approach, which is Tl0 is value from previous time step
                endif

                ! litter/air interface
                ! 1. get temperature and humidity at interface
                ! 2. get partial derivatives at interface
                sol(kk)%T1     = TL(kk) + deltaTa(kk)
                sol(kk)%Ta     = vmet(kk)%Ta
                sol(kk)%Rnet   = vmet(kk)%Rn
                sol(kk)%cva    = vmet(kk)%cva*thousand
                sol(kk)%hr1    = vlit(kk)%rh
                sol(kk)%hra    = vmet(kk)%rha
                sol(kk)%Dv     = vlit(kk)%Dv
                sol(kk)%gv     = one/vmet(kk)%rbw
                sol(kk)%gh     = one/vmet(kk)%rrc * rhocp
                sol(kk)%Dh     = vlit(kk)%kth
                sol(kk)%dz     = half*dxL(kk)
                sol(kk)%phie   = zero
                sol(kk)%he     = plit(kk)%he
                sol(kk)%K1     = zero
                sol(kk)%eta    = par(kk,1)%eta
                sol(kk)%lambda = par(kk,1)%lam
                sol(kk)%lambdav = var(kk,1)%lambdav

                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e5
                if (rhsurface(kk) > 0.9_r_2) tol(kk) = max((one-rhsurface(kk))*e4,e7)
                if ((tmp1*tmp2) < zero) then
                   rhsurface(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rhsurface(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (04). Stop."
                      stop
                   endif
                endif
                !Tsurface = T0_sol(rhsurface)
                hr0 = rhsurface(kk)
                Tsurface(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))

                hr0 = rhsurface(kk)
                TT0 = Tsurface(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                LlE0   = Qqevap*thousand*sol(kk)%lambdav
                GG0    = sol(kk)%Dh/sol(kk)%dz*(TT0-sol(kk)%T1)
                HH0    = sol(kk)%gh*(TT0-sol(kk)%Ta)
                qevap(kk) = Qqevap
                lE0(kk)   = LlE0
                G0(kk)    = GG0
                H(kk)     = HH0

                qpme(kk)  = qprec(kk) - qevap(kk)
                qL(kk)    = qprec(kk) - qevap(kk)
                qybL(kk)  = -vlit(kk)%rhS*sol(kk)%cva*sol(kk)%Dv/sol(kk)%dz
                qTbL(kk)  = -vlit(kk)%cvsatT *sol(kk)%gv*vlit(kk)%rh
                qhL(kk)   = G0(kk)
                qhTbL(kk) = -sol(kk)%Dh /sol(kk)%dz
                qhybL(kk) = zero
                Ebal0(kk) = vmet(kk)%Rn - lE0(kk) - H(kk) - G0(kk)

                ! perturb S and T to calcuate partial derivs qybL, qTbL, qhTbL, qhybL
                dTsoil_test = e1
                dS_test     = e2
                sol(kk)%hr1  = min((vlit(kk)%rh + vlit(kk)%rhS*dS_test), one)

                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e5
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e4,e7)
                if ((tmp1*tmp2) < zero) then
                   rh0_test(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0_test(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (05). Stop."
                      stop
                   endif
                endif

                !T0_test = T0_sol(rh0_test)
                hr0 = rh0_test(kk)
                T0_test(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))

                hr0 = rh0_test(kk)
                TT0 = T0_test(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                GG0    = sol(kk)%Dh/sol(kk)%dz*(TT0-sol(kk)%T1)
                qevap_test(kk) = Qqevap
                G0_test(kk)    = GG0
                ztmp = one/dS_test
                qybL(kk)  = -(qevap_test(kk)-qevap(kk))*ztmp
                qhybL(kk) = (G0_test(kk)-G0(kk))*ztmp

                sol(kk)%T1  = TL(kk) + deltaTa(kk) + dTsoil_test
                sol(kk)%hr1 = vlit(kk)%rh

                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e5
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e4,e7)
                if ((tmp1*tmp2) < zero) then
                   rh0_test(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0_test(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (06). Stop."
                      stop
                   endif
                endif

                !T0_test = T0_sol(rh0_test)
                hr0 = rh0_test(kk)
                T0_test(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))

                hr0 = rh0_test(kk)
                TT0 = T0_test(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                GG0    = sol(kk)%Dh/sol(kk)%dz*(TT0-sol(kk)%T1)
                qevap_test(kk) = Qqevap
                G0_test(kk)    = GG0
                ztmp = one/dTsoil_test
                qTbL(kk)  = -(qevap_test(kk)-qevap(kk))*ztmp
                qhTbL(kk) = (G0_test(kk)-G0(kk))*ztmp
                ! end of partial derivative evaluation for litter/air interface

                ! end litter/air interface

                ! soil/litter interface
                ! as before: 1. temp, humidity; 2. derivatives
                sol(kk)%T1      = Tsoil(kk,1)
                sol(kk)%Ta      = TL(kk) + deltaTa(kk)
                sol(kk)%cva     = vlit(kk)%cv*thousand
                sol(kk)%Rnet    = zero
                sol(kk)%hr1     = var(kk,1)%rh
                sol(kk)%hra     = vlit(kk)%rh
                sol(kk)%Dv      = var(kk,1)%Dv
                sol(kk)%gv      = (vlit(kk)%Dv/(half*dxL(kk)))
                sol(kk)%gh      = (vlit(kk)%kth/(half*dxL(kk)))
                sol(kk)%Dh      = var(kk,1)%kth
                sol(kk)%dz      = half*dx(kk,1)
                sol(kk)%phie    = par(kk,1)%phie
                sol(kk)%he      = par(kk,1)%he
                sol(kk)%K1      = var(kk,1)%K
                sol(kk)%eta     = par(kk,1)%eta
                sol(kk)%lambda  = par(kk,1)%lam
                sol(kk)%lambdav = var(kk,1)%lambdav

                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e5
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e4,e7)
                if ((tmp1*tmp2) < zero) then
                   rh0(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (07). Stop."
                      stop
                   endif
                endif
                !T0 = T0_sol(rh0)
                hr0 = rh0(kk)
                T0(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0 = rh0(kk)
                TT0 = T0(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                q(kk,0)    = -Qqevap
                qhya(kk,0) = zero
                qhyb(kk,0) = zero
                qhTa(kk,0) = sol(kk)%gv
                qhTb(kk,0) = -sol(kk)%Dh/sol(kk)%dz

                ! perturb S1, T1, SL, TL to calcuate partial derivs
                !    qyb(0),qhyb(0), qTb(0), qhTb(0), qya(0), qhya(0), qTa(0), qhTa(0)
                dTsoil_test = e1
                dS_test     = e2
                sol(kk)%hr1 = min((var(kk,1)%rh + var(kk,1)%rhS*dS_test), one)
                sol(kk)%K1  = min((var(kk,1)%K + var(kk,1)%KS*dS_test),par(kk,1)%Ke)
                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e4
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e3,e7)
                if ((tmp1*tmp2) < zero) then
                   rh0_test(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0_test(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (08). Stop."
                      stop
                   endif
                endif
                !T0_test = T0_sol(rh0_test)
                hr0 = rh0_test(kk)
                T0_test(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0 = rh0_test(kk)
                TT0 = T0_test(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                GG0    = sol(kk)%Dh/sol(kk)%dz*(TT0-sol(kk)%T1)
                qevap_test(kk) = Qqevap
                G0_test(kk)    = GG0
                ztmp = one/dS_test
                qyb(kk,0)  = -(qevap_test(kk)-(-q(kk,0)))*ztmp
                qhyb(kk,0) = (G0_test(kk)-qh(kk,0))*ztmp

                sol(kk)%T1  = Tsoil(kk,1) + dTsoil_test
                sol(kk)%hr1 = var(kk,1)%rh
                sol(kk)%K1  = var(kk,1)%K
                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e4
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e3,e7)
                if ((tmp1*tmp2) < zero) then
                   rh0_test(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0_test(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (09). Stop."
                      stop
                   endif
                endif
                !T0_test = T0_sol(rh0_test)
                hr0 = rh0_test(kk)
                T0_test(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0 = rh0_test(kk)
                TT0 = T0_test(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                GG0    = sol(kk)%Dh/sol(kk)%dz*(TT0-sol(kk)%T1)
                qevap_test(kk) = Qqevap
                G0_test(kk)    = GG0
                ztmp = one/dTsoil_test
                qTb(kk,0)  = -(qevap_test(kk)-(-q(kk,0)))*ztmp
                qhTb(kk,0) = (G0_test(kk)-qh(kk,0))*ztmp

                sol(kk)%T1  = Tsoil(kk,1)
                sol(kk)%cva = (vlit(kk)%cv+ vlit(kk)%cvS*dS_test)*thousand
                sol(kk)%hra = vlit(kk)%rh + vlit(kk)%rhS*dS_test
                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e4
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e3,e7)
                if ((tmp1*tmp2) < zero) then
                   rh0_test(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0_test(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (10). Stop."
                      stop
                   endif
                endif
                !T0_test = T0_sol(rh0_test)
                hr0 = rh0_test(kk)
                T0_test(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0 = rh0_test(kk)
                TT0 = T0_test(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                GG0    = sol(kk)%Dh/sol(kk)%dz*(TT0-sol(kk)%T1)
                qevap_test(kk) = Qqevap
                G0_test(kk)    = GG0
                ztmp = one/dS_test
                qya(kk,0)  = -(qevap_test(kk)-(-q(kk,0)))*ztmp
                qhya(kk,0) = (G0_test(kk)-qh(kk,0))*ztmp

                sol(kk)%cva = vlit(kk)%cv*thousand
                sol(kk)%hra = vlit(kk)%rh
                sol(kk)%Ta  = TL(kk) + deltaTa(kk) + dTsoil_test
                hr0  = e4 
                tmp1 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0  = one
                tmp2 = (sol(kk)%Dv*(-hr0 + sol(kk)%hr1)*csat(sol(kk)%T1))/sol(kk)%dz - &
                     (thousand*(sol(kk)%dz*sol(kk)%K1 + phi(hr0,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1) &
                     - phi(sol(kk)%hr1,sol(kk)%lambda,sol(kk)%eta,sol(kk)%phie,sol(kk)%he,sol(kk)%T1)))/sol(kk)%dz - &
                     (sol(kk)%Dv*sol(kk)%hr1*(sol(kk)%cva*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%Rnet &
                     + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta) &
                     - sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1)) + &
                     ((sol(kk)%DH + sol(kk)%dz*sol(kk)%gh)*sol(kk)%gv*(sol(kk)%cva - hr0*csat(sol(kk)%T1)) &
                     - sol(kk)%dz*sol(kk)%gv*hr0*(sol(kk)%Rnet + sol(kk)%gh*(-sol(kk)%T1 + sol(kk)%Ta))*slope_csat(sol(kk)%T1)) / &
                     (sol(kk)%DH + sol(kk)%dz*sol(kk)%gh + sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                tol(kk)  = e4
                if (rh0(kk) > 0.9_r_2) tol(kk) = max((one-rh0(kk))*e3,e7)
                if ((tmp1*tmp2) < zero) then
                   rh0_test(kk) = rtbis_rh0(sol(kk),e41d(kk),one1d(kk),tol(kk))
                else
                   if (tmp1 > zero .and. tmp2 > zero) then
                      rh0_test(kk) = one
                   else
                      write(2,*) "Found no solution for rh at surface (11). Stop."
                      stop
                   endif
                endif
                !T0_test = T0_sol(rh0_test)
                hr0 = rh0_test(kk)
                T0_test(kk) = ( sol(kk)%cva*sol(kk)%dz*sol(kk)%gv*sol(kk)%lambdav + sol(kk)%dz*sol(kk)%Rnet &
                     + sol(kk)%dh*sol(kk)%T1 + sol(kk)%dz*sol(kk)%gh*sol(kk)%Ta - &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*csat(sol(kk)%T1) + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*sol(kk)%T1*slope_csat(sol(kk)%T1) ) / &
                     (sol(kk)%dh + sol(kk)%dz*sol(kk)%gh + &
                     sol(kk)%dz*sol(kk)%gv*hr0*sol(kk)%lambdav*slope_csat(sol(kk)%T1))
                hr0 = rh0_test(kk)
                TT0 = T0_test(kk)
                Qqevap = sol(kk)%gv*(-sol(kk)%cva + hr0*(csat(sol(kk)%T1) + (TT0 - sol(kk)%T1)*slope_csat(sol(kk)%T1)))/thousand
                GG0    = sol(kk)%Dh/sol(kk)%dz*(TT0-sol(kk)%T1)
                qevap_test(kk) = Qqevap
                G0_test(kk)    = GG0
                ztmp = one/dTsoil_test
                qTa(kk,0)  = -(qevap_test(kk)-(-q(kk,0)))*ztmp
                qhTa(kk,0) = (G0_test(kk)-qh(kk,0))*ztmp

                ! end of partial derivative evaluation

                qd(kk) = zero ! drainage flux from litter to soil
                if (vlit(kk)%isat==1) then
                   qd(kk) = qprec(kk)
                   if (Sl(kk)>1) then
                      qd(kk) = qprec(kk) !+ (Sl-one)*dxL*plit%thre/dt
                   endif
                endif
                q(kk,0) = q(kk,0) + qd(kk)
                ! end litter/soil interface

             case (3) ! pond with and without litter
                deltaTa(kk) = zero
                dT0dTsoil(kk) = zero
                ! if ponding then evaporates with potential evaporation
                ! if (h0 > zero) then
                !    call potential_evap(vmet%Rn, vmet%rrc, vmet%rbw, vmet%Ta, vmet%rha,TL, kw, half*h0, &
                !         Tsurface_pot, Epot, Hpot, Gpot, dEdrha, dEdTa, dEdTsoil, dGdTa, dGdTsoil)
                ! elseif (h0<=zero) then
                call potential_evap(vmet(kk)%Rn, vmet(kk)%rrc, vmet(kk)%rbw, vmet(kk)%Ta, vmet(kk)%rha, &
                     Tsoil(kk,1), var(kk,1)%kth, half*dx(kk,1)+h0(kk), Tsurface_pot(kk), Epot(kk), Hpot(kk), &
                     Gpot(kk), dEdrha(kk), dEdTa(kk), dEdTsoil(kk), dGdTa(kk), dGdTsoil(kk))
                ! endif
                Tsurface(kk)  = Tsurface_pot(kk)
                T0(kk)        = Tsurface_pot(kk)
                rh0(kk)       = one
                rhsurface(kk) = one
                G0(kk)        = Gpot(kk)
                H(kk)         = Hpot(kk)
                lE0(kk)       = Epot(kk)
                Ebal0(kk)     = vmet(kk)%Rn - lE0(kk) - H(kk) - G0(kk)
                qevap(kk)     = Epot(kk)*e3/var(kk,1)%lambdav
                qpme(kk)        = qprec(kk) - qevap(kk)
                if (ns(kk)==0) then
                   ! if above maximum pond and precip-evap exceeds infiltration then do not solve for pond height
                   if (h0(kk)>=h0max .and. qpme(kk)>q(kk,0)) then
                      maxpond(kk) = .true.
                      ns(kk)        = 1
                  ! elseif (h0(kk)>=h0max .and. nsat(kk).eq.n.and.qpme(kk)>minval(var(kk,1:n)%Ksat).and.iflux(kk).eq.1) then ! steady state flow
				      !maxpond(kk) = .true.
                      !ns(kk)        = 1
				   else
				    
                      maxpond(kk) = .false.
                      ! pond included in top soil layer
                      var(kk,1)%phi = (h0(kk)-var(kk,1)%he)*var(kk,1)%Ksat+var(kk,1)%phie        
                      q(kk,0) = qpme(kk)
                      qya(kk,0) = zero
                      qyb(kk,0) = zero
                      qTa(kk,0) = zero
                      qTb(kk,0) = zero
                   end if
                endif

                qh(kk,0)   = G0(kk)
                qhya(kk,0) = zero
                qhyb(kk,0) = zero
                qhTa(kk,0) = zero
                qhTb(kk,0) = dGdTsoil(kk)

                qL(kk)    = zero
                qybL(kk)  = zero
                qTbL(kk)  = zero
                qhL(kk)   = zero
                qhTbL(kk) = zero
                qhybL(kk) = zero

                ql0(kk)   = -qevap(kk)
                qvT0(kk)  = zero
                qv0(kk)   = zero
                if (litter) then
                   qL(kk)  = qpme(kk)
                   qhL(kk) = G0(kk)
                endif

                ! liquid and vapour fluxes
                qliq(kk,0) = q(kk,0)
                qlya(kk,0) = qya(kk,0)
                qlyb(kk,0) = qyb(kk,0)
                qv(kk,0)   = zero
                qvya(kk,0) = zero
                qvyb(kk,0) = zero

             case default
                write(2,*) "solve: illegal pond case."
                stop
             end select ! pondcase
             ! finished all the surfaces

             ! get fluxes and derivatives (at time t=0, i.e. q0 etc.)
             call getfluxes_vp(n, ns(kk), dx(kk,1:n), vtop(kk), vbot(kk), par(kk,1:n), var(kk,1:n), & ! moisture fluxes
                  hint(kk,1:n), phimin(kk,1:n), q(kk,0:n), qya(kk,0:n), qyb(kk,0:n), qTa(kk,0:n), qTb(kk,0:n), &
                  qliq(kk,0:n), qlya(kk,0:n), qlyb(kk,0:n), qv(kk,0:n), qvT(kk,0:n), qvh(kk,0:n), qvya(kk,0:n), &
                  qvyb(kk,0:n), iflux(kk), init(kk), getq0(kk), getqn(kk), Tsoil(kk,1:n), T0(kk), nsat(kk), nsatlast(kk))
             qTa(kk,n) = zero
             qTb(kk,n) = zero
             call getheatfluxes(ns(kk:kk), h0(kk:kk), dx(kk:kk,1:n), dxL(kk:kk), &
                  qh(kk:kk,0:n), qhya(kk:kk,0:n), qhyb(kk:kk,0:n), qhTa(kk:kk,0:n), qhTb(kk:kk,0:n), &
                  var(kk:kk,1:n), vlit(kk:kk), Tsoil(kk:kk,1:n), TL(kk:kk), T0(kk:kk), litter, &
                  q(kk:kk,0:n), qya(kk:kk,0:n), qyb(kk:kk,0:n), qTa(kk:kk,0:n), qTb(kk:kk,0:n), &
                  qadv(kk:kk,0:n),qadvya(kk:kk,0:n), qadvyb(kk:kk,0:n), qadvTa(kk:kk,0:n), qadvTb(kk:kk,0:n), &
				  advection) ! heat fluxes

             if (ns(kk).eq.0) then ! pond included in top soil layer
                ! change qya(1) from dq/dphi (returned by getfluxes) to dq/dh
                qya(kk,1) = var(kk,1)%Ksat*qya(kk,1)
                if (advection==1) then
                   qhya(kk,1) = qhya(kk,1) - qadvya(kk,1)
                   w(kk) = (var(kk,1)%kth/dx(kk,1))/(var(kk,1)%kth/dx(kk,1)+var(kk,2)%kth/dx(kk,2))
                   qadvya(kk,1) =  rhow*cswat*qya(kk,1)*(w(kk)*(Tsoil(kk,1)+Tzero)+(1.-w(kk))*(Tsoil(kk,2)+Tzero))
                   qhya(kk,1) = qhya(kk,1) + qadvya(kk,1)
                endif
             endif



             ! adjust for top and bottom bdry condns
             if (botbc=="zero flux") then
                qliq(kk,n) = zero
                qv(kk,n)   = zero
                q(kk,n)    = zero
                qya(kk,n)  = zero
                qlya(kk,n) = zero
                qvya(kk,n) = zero
             endif

             ! specify mositure flux at bottom of soil column (heat flux set to zero)
             if (botbc /= "constant head") then
                select case (botbc)
                case ("zero flux")
                   q(kk,n)   = zero
                   qya(kk,n) = zero
                case ("free drainage")
                   q(kk,n) = gf*var(kk,n)%K
                   if (var(kk,n)%isat == 0) then
                      qya(kk,n) = gf*var(kk,n)%KS
                   else
                      qya(kk,n) = zero
                   end if
                case ("seepage")
                   if (var(kk,n)%h <= -half*gf*dx(kk,n)) then
                      q(kk,n)   = zero
                      qya(kk,n) = zero
                   end if
                case ("aquifer")
                   if (v_aquifer(kk)%isat==0) then
                      q(kk,n)   = v_aquifer(kk)%K*(one+var(kk,n)%h / &
                           (v_aquifer(kk)%zdelta-(v_aquifer(kk)%zsoil - half*gf*dx(kk,n))))
                      qya(kk,n) = zero
                   else
                      q(kk,n) = gf*var(kk,n)%K
                      if (var(kk,n)%isat == 0) then
                         qya(kk,n) = gf*var(kk,n)%KS
                      else
                         qya(kk,n) = zero
                      end if
                   end if
                case default
                   write(2,*) "solve: illegal bottom boundary condition."
                   stop
                end select
             end if

             if (testcase==7 .or. testcase==8) then
                q(kk,n)   = q(kk,0)
                qya(kk,n) = zero
             endif
             if (testcase==8) qh(kk,n) = G0(kk)

             ! adjust upper and lower heat flux for advection
             if (advection==1) then
                qadv(kk,n) = rhow*cswat*(Tsoil(kk,n)+Tzero)*q(kk,n)
                qadvya(kk,n) = rhow*cswat*(Tsoil(kk,n)+Tzero)*qya(kk,n)
                qadvTa(kk,n)= rhow*cswat*q(kk,n)
                qh(kk,n) = qh(kk,n) + qadv(kk,n)
                qhya(kk,n) = qhya(kk,n) + qadvya(kk,n) 
                qhTa(kk,n) = qhTa(kk,n) + qadvTa(kk,n) 

                !qadv(kk,0) = rhow*cswat*q(kk,0)*(T0(kk)+Tzero)
                qadv(kk,0) = rhow*cswat*q(kk,0)*(max(vmet(kk)%Ta,0.01)+Tzero)
                qh(kk,0) = qh(kk,0) + qadv(kk,0)
                !qadvyb(kk,0) =  rhow*cswat*qyb(kk,0)*(T0(kk)+Tzero)
                qadvyb(kk,0) =  rhow*cswat*qyb(kk,0)*(max(vmet(kk)%Ta,0.01)+Tzero)
                qadvTb(kk,0) = dT0dTsoil(kk)*rhow*cswat*q(kk,0)
                qhyb(kk,0) = qhyb(kk,0) +  qadvyb(kk,0)
                qhTb(kk,0) = qhTb(kk,0) + qadvTb(kk,0)
             else
                qadv(:,:) = zero; qadvya(:,:)=zero; qadvyb(:,:)=zero
                qadvTa(:,:) = zero; qadvTb(:,:) = zero
             endif

			 !if (nsat(kk).eq.n.and.advection.eq.1) then
			!	tmp1d1(kk) = minval(var(kk,1:n)%Ksat)
			!	do i=1,n-1
			!	    qh(kk,i) = qh(kk,i)-qadv(kk,i)
			!		qhya(kk,i) = qhya(kk,i) - qadvya(kk,i)
			!		qhyb(kk,i) = qhyb(kk,i) - qadvyb(kk,i)
			!		qhTa(kk,i) = qhTa(kk,i) - qadvTa(kk,i)
			!		qhTb(kk,i) = qhTb(kk,i) - qadvTb(kk,i)

			!		tmp1d2(kk) = (var(kk,i)%kth/dx(kk,i))/(var(kk,i)%kth/dx(kk,i)+var(kk,i+1)%kth/dx(kk,i+1))
			!		qadv(kk,i) = rhow*cswat*tmp1d1(kk)*(tmp1d2(kk)*(Tsoil(kk,i)+Tzero)+(1.-tmp1d2(kk))*(Tsoil(kk,i+1)+Tzero))
			!		qadvya(kk,i) =  zero
			!		qadvyb(kk,i) =  zero
			!		qadvTa(kk,i) =  rhow*cswat*tmp1d1(kk)*tmp1d2(kk)
			!		qadvTb(kk,i) =   rhow*cswat*tmp1d1(kk)*(1.-tmp1d2(kk))
			!		qh(kk,i) = qh(kk,i)+qadv(kk,i)
			!		qhya(kk,i) = qhya(kk,i) + qadvya(kk,i)
			!		qhyb(kk,i) = qhyb(kk,i) + qadvyb(kk,i)
			!		qhTa(kk,i) = qhTa(kk,i) + qadvTa(kk,i)
			!		qhTb(kk,i) = qhTb(kk,i) + qadvTb(kk,i)
			!	enddo
			!	qh(kk,0) = qh(kk,0) - qadv(kk,0)
			!	qadv(kk,0) = rhow*cswat*tmp1d1(kk)*(max(vmet(kk)%Ta,0.01)+Tzero)
             !   qh(kk,0) = qh(kk,0) + qadv(kk,0)
			!	qhyb(kk,0) = qhyb(kk,0) -  qadvyb(kk,0)
            !    qhTb(kk,0) = qhTb(kk,0) - qadvTb(kk,0)
             !   qadvyb(kk,0) =  zero
             !   qadvTb(kk,0) = zero
             !   qhyb(kk,0) = qhyb(kk,0) +  qadvyb(kk,0)
             !   qhTb(kk,0) = qhTb(kk,0) + qadvTb(kk,0)

			!	qh(kk,n) = qh(kk,n) - qadv(kk,n)
			!	qhya(kk,n) = qhya(kk,n) - qadvya(kk,n) 
            !    qhTa(kk,n) = qhTa(kk,n) - qadvTa(kk,n) 
		!		qadv(kk,n) = rhow*cswat*(Tsoil(kk,n)+Tzero)*tmp1d1(kk)
            !    qadvya(kk,n) = zero
            !    qadvTa(kk,n)= rhow*cswat*tmp1d1(kk)
            !    qh(kk,n) = qh(kk,n) + qadv(kk,n)
            !    qhya(kk,n) = qhya(kk,n) + qadvya(kk,n) 
            !    qhTa(kk,n) = qhTa(kk,n) + qadvTa(kk,n) 
			! endif

             ! get rate of extraction (roots) - ignore solute
             !  if (Etrans>zero) then
             !call getrex(S(kk,:), qex(kk,:), fws(kk), FS(kk,:), par(kk,:)%the, par(kk,:)%thw, Etrans(kk), gamma(kk))
             ! Ross had a derivative. Set it to zero for the moment
             !qexd = qexd*var%phiS/var%K ! to get deriv qexS
!!! MC !!! might have to change with more sophisticated soil water uptake
!!! MC !!! investigate or leave it for the moment?
!!! VH !!! leave for the moment: would only change for different root uptake
             qexd(kk,1:n) = zero
             fwsoil(kk)   = fws(kk)
             !  end if
             again(kk)  = .false. ! flag for recalcn of fluxes (default=false)
             !----- end get fluxes and derivs

             !----- first estimate of time step dt before the calculation
             !      gets revised after the calculation
             dmax(kk)     = zero
             tmp2d1(kk,:) = zero
             tmp2d2(kk,:) = zero !  temp storage
             ! estimate rate of change of moisture storage [m/s]
             where (var(kk,1:n)%isat==0) tmp2d1(kk,1:n) = &
                  abs(q(kk,1:n)-q(kk,0:n-1)-qex(kk,1:n))/(par(kk,1:n)%thre*dx(kk,1:n))
             ! estimate rate of change of heat storage [K/s]
             tmp2d2(kk,1:n) = abs(qh(kk,1:n)-qh(kk,0:n-1))/(var(kk,1:n)%csoileff*dx(kk,1:n))
			 if (advection.eq.1) then
				!tmp2d2(kk,1:n) = abs(qh(kk,1:n)-qh(kk,0:n-1)+rhow*cswat*qex(kk,1:n)*(Tsoil(kk,1:n)+Tzero))/(var(kk,1:n)%csoileff*dx(kk,1:n))
               tmp2d2(kk,1:n) = abs((qh(kk,1:n)-qadv(kk,1:n))-(qh(kk,0:n-1)-qadv(kk,0:n-1)))/(var(kk,1:n)%csoileff*dx(kk,1:n))
			 endif
			 if (litter .and. ns(kk)==1 .and. (.not. maxpond(kk))) then ! litter , no pond
                if (vlit(kk)%isat==0) then ! estimate rate of change of moisture storage [m/s]
                   tmp2d1(kk,0) = abs(q(kk,0) - qL(kk))/(plit(kk)%thre*dxL(kk))
                endif
                tmp2d2(kk,0) = abs(qh(kk,0) - qhL(kk))/(vlit(kk)%csoileff*dxL(kk)) ! estimate rate of change of heat storage [K/s]
                tmp1d3(kk) = (dTLmax-abs(deltaTa(kk))) / tmp2d2(kk,0)
             else
                tmp2d1(kk,0) = zero
                tmp2d2(kk,0) = zero
                tmp1d3(kk)   = dtmax
             endif

             dmax(kk) = maxval(tmp2d1(kk,:),1) ! max derivative |dS/dt|

             if (abs(minval(tmp2d1(kk,:),1)) > maxval(tmp2d1(kk,:),1)) then
                write(2,*) 'Should not be here (01)'
                stop
                dmax(kk) = minval(tmp2d1(kk,:),1)
             endif
             tmp1d1(kk) = maxval(tmp2d2(kk,0:n),1) ! max derivative |dTsoil/dt|
             if (dmax(kk) > zero) then
                dt(kk) = min(dSmax/dmax(kk), dTLmax/tmp1d1(kk), tmp1d3(kk)) ! constrained either by moisture or temp
                ! if pond, overwrite dt
                if (h0(kk)>zero .and. (q(kk,0)-qpme(kk))*dt(kk)>h0(kk).and.ns(kk)==1) &
                     dt(kk) = (h0(kk)-half*h0min)/(q(kk,0)-qpme(kk))
                if (h0(kk)>zero .and. (q(kk,1)-qpme(kk))*dt(kk)>h0(kk).and.ns(kk)==0) &
                     dt(kk) = (h0(kk)-half*h0min)/(q(kk,1)-qpme(kk))
             else ! steady state flow
                if (qpme(kk)>=q(kk,n)) then ! if saturated soil columnn and more precip then drainige -> finish
                   dt(kk) = tfin-t(kk) ! step to finish
                else ! otherwise adjust dt because change of pond height
                   dt(kk) = -(h0(kk)-half*h0min)/(qpme(kk)-q(kk,n))
                end if
             end if
             if (dt(kk)>dtmax) dt(kk) = dtmax ! user's limit

             ! if initial step, improve phi where S>=1
             ! might be that you get better derivatives, especially at sat/non-sat interfaces
             if (nsteps(kk)==nsteps0(kk) .and. nsat(kk)>0 .and. iflux(kk)==1) then
                again(kk) = .true.
                dt(kk)    = dtmin
             end if
             ! if fully saturated but was not fully saturated before, adjust even within each time step iteration
             if (nsat(kk)==n .and. nsatlast(kk)<n .and. iflux(kk)==1) then
                again(kk) = .true. ! profile has just become saturated so adjust phi values
                dt(kk)    = dtmin
             end if
             ! sprint to the end
             if (t(kk)+1.1_r_2*dt(kk)>tfin) then ! step to finish
                dt(kk) = tfin-t(kk)
                t(kk)  = tfin
             else
                t(kk) = t(kk)+dt(kk) ! tentative update
                if (again(kk)) t(kk) = t(kk)-dt(kk)
             end if

			 !write(339,"(i4,i4,14f15.6)") irec, nsteps(kk), dt(kk), dSmax/max(dmax(kk),1e-6), dTLmax/tmp1d1(kk)
             !----- end estimate time step dt

             !----- get and solve eqns
             rsigdt(kk) = one/(sig(kk)*dt(kk))

             ! aa, bb, cc and dd hold coeffs and rhs of linear equation eqn set
             if (septs == 1) then ! uncoupling of T and S
                qTa(kk,:)   = zero
                qTb(kk,:)   = zero
                qhya(kk,:)  = zero
                qhyb(kk,:)  = zero
                qTbL(kk)    = zero
                qhybL(kk)   = zero
             endif
             aa(kk,1:n)   =  qya(kk,0:n-1)
             ee(kk,0:n-1) = -qyb(kk,0:n-1)
             bb(kk,1:n)   =  qTa(kk,0:n-1)
             ff(kk,0:n-1) = -qTb(kk,0:n-1)

             aah(kk,1:n)   =  qhya(kk,0:n-1)
             eeh(kk,0:n-1) = -qhyb(kk,0:n-1)
             bbh(kk,1:n)   =  qhTa(kk,0:n-1)
             ffh(kk,0:n-1) = -qhTb(kk,0:n-1)


             gg(kk,1:n) = -(q(kk,0:n-1)-q(kk,1:n)-qex(kk,1:n))*rsig(kk)
			

             iok(kk)            = 0 ! flag for time step test
             itmp(kk)           = 0 ! counter to abort if not getting solution
             again_ice(kk,:)      = .false.
             nsteps_ice(kk,1:n) = 0
             do while (iok(kk)==0) ! keep reducing time step until all ok
                itmp(kk)  = itmp(kk)+1
                accel(kk) = one - 0.05_r_2*real(min(10,max(0,itmp(kk)-4)),r_2) ! acceleration [0.5,1], start with 1
                if (itmp(kk) > 1000) then
                   write(2,*) "solve: too many iterations of equation solution"
                   !write(2,*) " irec, k, SL, TL, deltaTa, h0"
                   !write(2,*) irec, k, SL(kk), TL(kk), deltaTa(kk), h0(kk)
                   stop
                end if
				 ggh(kk,1:n) = -(qh(kk,0:n-1)-qh(kk,1:n))*rsig(kk)
				if (advection==1) then
					ggh(kk,1:n) = ggh(kk,1:n) + qex(kk,1:n)*rsig(kk)*(Tsoil(kk,1:n)+Tzero)*cswat*rhow
				endif

                if (ns(kk)<1) then ! ns=0 then pond height otherwise litter moisture
                   cc(kk,0) = -qya(kk,0) - rsigdt(kk)
                else
                   cc(kk,0) = -qya(kk,0) - rsigdt(kk)*plit(kk)%thre*dxL(kk) + qybL(kk)
                end if
                gg(kk,0)  = -(qprec(kk)-qevap(kk)-q(kk,0))*rsig(kk)
                ggh(kk,0) = -(G0(kk)-qh(kk,0))*rsig(kk)
                dd(kk,0)  = -qTa(kk,0)

                where (var(kk,1:n)%isat==0)
                   cc(kk,1:n) = qyb(kk,0:n-1) - qya(kk,1:n) - par(kk,1:n)%thre*dx(kk,1:n)*rsigdt(kk) - qexd(kk,1:n)
                elsewhere
                   cc(kk,1:n) = qyb(kk,0:n-1) - qya(kk,1:n) - qexd(kk,1:n)
                endwhere

                if (ns(kk)<1) then ! pond included in top soil layer
                   cc(kk,1) = -qya(kk,1)-rsigdt(kk) -qexd(kk,1)
                endif

                dd(kk,1:n)  = qTb(kk,0:n-1)-qTa(kk,1:n)
                ddh(kk,1:n) = qhTb(kk,0:n-1)-qhTa(kk,1:n)-var(kk,1:n)%csoileff*dx(kk,1:n)*rsigdt(kk)

				
                cch(kk,1:n) = qhyb(kk,0:n-1)-qhya(kk,1:n)
                ddh(kk,1) = ddh(kk,1) - h0(kk)*rhow*cswat*rsigdt(kk)*real(-var(kk,1)%iice+1) &
                     - h0(kk)*rhow*csice*rsigdt(kk)*real(var(kk,1)%iice)
                where (var(kk,1:n)%isat==0) cch(kk,1:n) = cch(kk,1:n) + &
                     real(var(kk,1:n)%iice,r_2)*rhow*lambdaf*par(kk,1:n)%thre*dx(kk,1:n)*rsigdt(kk)
                aa(kk,0)  = qya(kk,0)
                aah(kk,0) = zero
                bbh(kk,0) = zero

                if (advection==1) then
                   if (ns(kk).eq.0) then  ! pond included in top soil layer
                      cch(kk,1) = cch(kk,1) -rhow*cswat*(Tsoil(kk,1)+Tzero)*rsigdt(kk)*real(-var(kk,1)%iice+1) &
                           -rhow*csice*(Tsoil(kk,1)+Tzero)*rsigdt(kk)*real(var(kk,1)%iice) 
                   else
                      cch(kk,1) = cch(kk,1) -rhow*cswat*(Tsoil(kk,1)+Tzero)*par(kk,1)%thre*dx(kk,1)*rsigdt(kk)* &
                           real(-var(kk,1)%isat+1)*real(-var(kk,1)%iice+1) &
                           -rhow*csice*(Tsoil(kk,1)+Tzero)*par(kk,1)%thre*dx(kk,1)*rsigdt(kk)* &
                           real(-var(kk,1)%isat+1)*real(var(kk,1)%iice)

                   endif
                   cch(kk,2:n) = cch(kk,2:n) -rhow*cswat*(Tsoil(kk,2:n)+Tzero)*par(kk,2:n)%thre*dx(kk,2:n)*rsigdt(kk)* &
                        real(-var(kk,2:n)%isat+1)*real(-var(kk,1)%iice+1) &
                        -rhow*csice*(Tsoil(kk,2:n)+Tzero)*par(kk,2:n)%thre*dx(kk,2:n)*rsigdt(kk)*real(-var(kk,2:n)%isat+1)* &
                        real(var(kk,1)%iice)
                endif

                if (litter .and. ns(kk)==1) then ! litter and no pond
                   ! watch for deltaTa
                   cc(kk,0)  = qybL(kk) -qya(kk,0) -rsigdt(kk)*plit(kk)%thre*dxL(kk)
                   dd(kk,0)  = qTbL(kk) - qTa(kk,0)
                   ee(kk,0)  = -qyb(kk,0)
                   ff(kk,0)  = -qTb(kk,0)
                   gg(kk,0)  = (q(kk,0) - qL(kk))*rsig(kk) + deltaTa(kk)*(qTbL(kk)-qTa(kk,0))
                   gg(kk,1)  = -(q(kk,0)-q(kk,1)-qex(kk,1))*rsig(kk) + deltaTa(kk)*qTa(kk,0)
                   cch(kk,0) = qhybL(kk) - qhya(kk,0)
                   ddh(kk,0) = -qhTa(kk,0)-vlit(kk)%csoil*dxL(kk)*rsigdt(kk)+qhTbL(kk)
                   eeh(kk,0) = -qhyb(kk,0)
                   ffh(kk,0) = -qhTb(kk,0)
                   ggh(kk,0) = (qh(kk,0)-qhL(kk))*rsig(kk) + deltaTa(kk)*(qhTbL(kk) - qhTa(kk,0))
                   ggh(kk,1) = -(qh(kk,0)-qh(kk,1))*rsig(kk) + deltaTa(kk)*qhTa(kk,0)
                endif

                if (litter .and. ns(kk)==0) then ! litter and pond
                   ddh(kk,0) = -qhTa(kk,0)-cswat*rhow*h0(kk)*rsigdt(kk) -vlit(kk)%csoil*dxL(kk)*rsigdt(kk)
                endif

                if ((.not. litter) .and. ns(kk)==0) then ! no litter and pond
                   cch(kk,0) = zero
                   !ddh(kk,0) = -qhTa(kk,0)-cw*max(h0(kk),e3)*rsigdt(kk)
                   ddh(kk,0) = one
                   eeh(kk,0) = zero
                   ffh(kk,0) = zero
                   ggh(kk,0) = one ! these assignments have the effect of decoupling Tsoil(0) from other layers
                   dd(kk,0)  = zero
                   bb(kk,1)  = zero
                   bbh(kk,1) = zero
                endif

                if (maxpond(kk) .and. ns(kk)==1) then ! no litter and pond
                   cc(kk,0) = zero
                   dd(kk,0) = zero
                   ee(kk,0) = zero
                   ff(kk,0) = zero
                   gg(kk,0) = zero
                   cch(kk,0) = zero
                   ddh(kk,0) = zero
                   eeh(kk,0) = zero
                   ffh(kk,0) = zero
                   ggh(kk,0) = zero
                endif

                if (septs == 1) then ! uncoupled of T and S
                   ! Solve two tridiagonal matrices instead of 1 block matrix
                   ! Could try two different time steps for moisture and temperature
!!! MC !!! This works now and gives about 10% speed increase.
!!! MC !!!    Do you think it is worth to rewrite the code
!!! MC !!!       to separate moisture and temperature?
!!! MC !!!    Might speed up things further but first we should see
!!! MC !!!       if vapour transport is needed somewhere in Australia or not.
!!! MC !!!    This was part of your Australia project, isn't it.
!!! VH !!!    leave for the moment.
!!! VH !!!    Could check to see how often T is limiting time-step and by how much
                   if ((ns(kk) == 1 .or. maxpond(kk)) .and. (.not. litter)) then ! no pond or maxpond, no litter
                      nns(kk) = 1
                   else
                      nns(kk) = 0
                   endif
                   !nn(kk) = n-nns(kk)+1
                   !call dgtsv(nn, 1, aa(nns+1:n), cc(nns:n), ee(nns:n-1), gg(nns:n), nn, info)
                   !dy(nns:n) = gg(nns:n)
                   call tri(nns(kk), n, aa(kk,0:n), cc(kk,0:n), ee(kk,0:n), gg(kk,0:n), dy(kk,0:n))
                   if (nns(kk) == 1) dy(kk,0) = zero
                   !call dgtsv(nn, 1, bbh(nns+1:n), ddh(nns:n), ffh(nns:n-1), ggh(nns:n), nn, info)
                   !dTsoil(nns:n) = ggh(nns:n)
                   call tri(nns(kk), n, bbh(kk,0:n), ddh(kk,0:n), ffh(kk,0:n), gg(kk,0:n), dTsoil(kk,0:n))
                   if (nns(kk) == 1) dTsoil(kk,0) = zero
                   if (nns(kk)==0 .and. h0(kk)<e3 .and. (.not. litter)) dTsoil(kk,0) = zero
                else ! coupled of T and S
                   nns(kk) = 1  ! pond included in top soil layer
                   dy(kk,0)     = zero ! to be safe if nns=1
                   dTsoil(kk,0) = zero
                   call massman_sparse(aa(kk,nns(kk)+1:n), aah(kk,nns(kk)+1:n), bb(kk,nns(kk)+1:n), bbh(kk,nns(kk)+1:n), &
                        cc(kk,nns(kk):n), cch(kk,nns(kk):n), dd(kk,nns(kk):n), ddh(kk,nns(kk):n), ee(kk,nns(kk):n-1), &
                        eeh(kk,nns(kk):n-1), &
                        ff(kk,nns(kk):n-1), ffh(kk,nns(kk):n-1), gg(kk,nns(kk):n), ggh(kk,nns(kk):n), &
                        dy(kk,nns(kk):n), dTsoil(kk,nns(kk):n),condition=condition)
                   if (nns(kk)==1) dy(kk,0) = zero
                   if (nns(kk)==0 .and. h0(kk)<one .and. (.not. litter)) dTsoil(kk,0) = zero


                   ! evaluate soil fluxes at sigma of time step for use in isotope routine

				  
				   qsig(kk,0)  = q(kk,0)+sig(kk)*qyb(kk,0)*dy(kk,1) + sig(kk)*qya(kk,0)*dy(kk,0) + &
                        sig(kk)*qTb(kk,0)*dTsoil(kk,1) + sig(kk)*qTa(kk,0)*dTsoil(kk,0)

                   qhsig(kk,0) = qh(kk,0) + sig(kk)*qhyb(kk,0)*dy(kk,1) + sig(kk)*qhya(kk,0)*dy(kk,0) + &
                        sig(kk)*qhTb(kk,0)*dTsoil(kk,1) + sig(kk)*qhTa(kk,0)*dTsoil(kk,0)


                   if (ns(kk)<1) then
                      qlsig(kk,0) = qsig(kk,0)
                      qvsig(kk,0) = zero
                   else
                      qvsig(kk,0) = qsig(kk,0)
                      qlsig(kk,0) = zero
                   endif

                   qadvsig(kk,0) = qadv(kk,0) + sig(kk)*qadvyb(kk,0)*dy(kk,1) &
                        + sig(kk)*qadvya(kk,0)*dy(kk,0) + sig(kk)*qadvTb(kk,0)*dTsoil(kk,1) + sig(kk)*qadvTa(kk,0)*dTsoil(kk,0)

                   qsig(kk,1:n-1) = q(kk,1:n-1) + sig(kk)*(qya(kk,1:n-1)*dy(kk,1:n-1)+qyb(kk,1:n-1)*dy(kk,2:n) &
                        +qTa(kk,1:n-1)*dTsoil(kk,1:n-1)+qTb(kk,1:n-1)*dTsoil(kk,2:n))
                   qsig(kk,n)     = q(kk,n) + sig(kk)*(qya(kk,n)*dy(kk,n)+qTa(kk,n)*dTsoil(kk,n))

                   qhsig(kk,1:n-1) = qh(kk,1:n-1) + sig(kk)*(qhya(kk,1:n-1)*dy(kk,1:n-1)+qhyb(kk,1:n-1)*dy(kk,2:n) &
                        +qhTa(kk,1:n-1)*dTsoil(kk,1:n-1)+qhTb(kk,1:n-1)*dTsoil(kk,2:n))
                   qhsig(kk,n) = qh(kk,n) + sig(kk)*(qhya(kk,n)*dy(kk,n)+qhTa(kk,n)*dTsoil(kk,n))

                   qadvsig(kk,1:n-1) = qadv(kk,1:n-1) + sig(kk)*(qadvya(kk,1:n-1)*dy(kk,1:n-1)+qadvyb(kk,1:n-1)*dy(kk,2:n) &
                        +qadvTa(kk,1:n-1)*dTsoil(kk,1:n-1)+qadvTb(kk,1:n-1)*dTsoil(kk,2:n))
                   qadvsig(kk,n) = qadv(kk,n) + sig(kk)*(qadvya(kk,n)*dy(kk,n)+qadvTa(kk,n)*dTsoil(kk,n))


                   LHS_h(kk,1) = (dx(kk,1)*var(kk,1)%csoileff+h0(kk)*rhow*cswat)*dTsoil(kk,1)/dt(kk) - &
                        dx(kk,1)*dy(kk,1)/dt(kk)*par(kk,1)%thre*var(kk,1)%iice*rhow*lambdaf*real(-var(kk,1)%isat+1)
                   RHS_h(kk,1) = qhsig(kk,0) -qhsig(kk,1) 
                   do i=1,n
                      if(i ==1) then
                         LHS_h(kk,i) = (dx(kk,i)*var(kk,i)%csoileff+h0(kk)*rhow*cswat)*dTsoil(kk,i)/dt(kk) &
                              - dx(kk,i)*dy(kk,i)/dt(kk)*par(kk,i)%thre*var(kk,i)%iice*rhow*lambdaf*real(-var(kk,i)%isat+1)
                      else
                         LHS_h(kk,i) = (dx(kk,i)*var(kk,i)%csoileff)*dTsoil(kk,i)/dt(kk) &
                              - dx(kk,i)*dy(kk,i)/dt(kk)*par(kk,i)%thre*var(kk,i)%iice*rhow*lambdaf*real(-var(kk,i)%isat+1)
                      endif

                      if (advection==1) then
                         if (ns(kk)==0 .and. i == 1) then
                            LHS_h(kk,i) = LHS_h(kk,i) + rhow*cswat*(Tsoil(kk,i)+Tzero)*dy(kk,i)*real(-var(kk,1)%iice+1)/dt(kk) &
                                 + rhow*csice*(Tsoil(kk,i)+Tzero)*dy(kk,i)*real(var(kk,1)%iice)/dt(kk) 
                         else
                            LHS_h(kk,i) = LHS_h(kk,i) + rhow*cswat*(Tsoil(kk,i)+Tzero)*dx(kk,i)*dy(kk,i) &
                                 *par(kk,i)%thre*real(-var(kk,i)%isat+1)/dt(kk)
                         endif
                      endif

					  if (advection==1) then
						RHS_h(kk,i) = RHS_h(kk,i) - qex(kk,i)*(Tsoil(kk,i)+Tzero)*cswat*rhow
					  endif
                      RHS_h(kk,i) = qhsig(kk,i-1) -qhsig(kk,i) 
					  
					  if (advection==1) then
						RHS_h(kk,i) = RHS_h(kk,i) - qex(kk,i)*(Tsoil(kk,i)+Tzero)*cswat*rhow
					  endif
                   enddo



                   dJcol = sum(LHS_h(:,1:n),2)*dt(:)

                   ! check mass balance on top soil layer
                   if (ns(kk).eq.0) then  ! pond included in top soil layer
                      LHS(kk,1) = dy(kk,1)/dt(kk)
                   else
                      LHS(kk,1) = dy(kk,1)*par(kk,1)%thre*dx(kk,1)*real(-var(kk,1)%isat+1)/dt(kk)
                   endif
                   RHS(kk,1) = qsig(kk,0) - qsig(kk,1)- qex(kk,1)


                endif ! septs==1

                ! dy contains dS or, for sat layers, dphi values
                iok(kk) = 1
                fac(kk) = one

                if (.not. again(kk)) then

                   ! check if time step ok, if not then set fac to make it less
                   iok(kk) = 1
                   do i=1, n
                      if (var(kk,i)%isat==0) then ! check change in S in initially unsaturated layers
                         if (abs(dy(kk,i)) > dSfac*dSmax) then
                            fac(kk) = max(half,accel(kk)*abs(dSmax/dy(kk,i)))
                            iok(kk) = 0
                            exit
                         end if
                         if (-dy(kk,i) > dSmaxr*S(kk,i)) then ! Check relative moisture change
                            fac(kk) = max(half,accel(kk)*dSmaxr*S(kk,i)/(-dSfac*dy(kk,i)))
                            iok(kk) = 0
                            exit
                         end if
                         if (S(kk,i) < one .and. S(kk,i)+dy(kk,i) > Smax) then ! Check for oversaturating
                            fac(kk) = accel(kk)*(half*(one+Smax)-S(kk,i))/dy(kk,i)
                            iok(kk) = 0
                            exit
                         end if
                         if (S(kk,i) >= one .and. dy(kk,i) > half*(Smax-one)) then ! Check for limit at oversaturation
                            fac(kk) = 0.25_r_2*(Smax-one)/dy(kk,i)
                            iok(kk) = 0
                            exit
                         end if
                      end if
                   end do

                   do i=1, n
                      if (abs(dTsoil(kk,i)) > dTsoilmax) then ! Check absolute soil temperature change
                         fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTsoilmax/dTsoil(kk,i))))
                         iok(kk) = 0
                         exit
                      end if
                   enddo

                   if (litter .and. ns(kk)==0) then ! litter and pond
                      if (abs(dTsoil(kk,0)) > dTLmax) then ! Check absolute litter temperature change
                         fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTLmax/dTsoil(kk,0))))
                         iok(kk) = 0
                      end if
                   endif

                   if (litter .and. ns(kk)==1) then ! litter, no pond
                      if (-dy(kk,0) > SL(kk)*dSmaxr) then ! Check relative litter moisture change
                         fac(kk) = min(fac(kk),max(half,accel(kk)*dSmaxr*SL(kk)/(-dSfac*dy(kk,0))))
                         iok(kk) = 0
                      endif

                      if (abs(deltaTa(kk)) > zero) then
                         tmp1d1(kk) = dTsoil(kk,0) - deltaTa(kk)
                         if (abs(tmp1d1(kk)) > dTLmax) then ! Check absolute litter temperature change
                            fac(kk) = min(fac(kk),max(half,accel(kk)*dTLmax/abs(tmp1d1(kk))))
                            iok(kk) = 0
                            if (itmp(kk) > 10) then
                               iok(kk) = 1
                               t(kk)   = t(kk)-dt(kk)
                               n_noconverge(kk) = n_noconverge(kk)+1
                               goto 10
                            endif
                         endif
                      else
                         if (abs(dTsoil(kk,0)) > dTLmax) then ! Check absolute litter temperature change
                            fac(kk) = min(fac(kk),max(half,accel(kk)*abs(dTLmax/dTsoil(kk,0))))
                            iok(kk) = 0
                         end if
                      endif
                   endif ! litter, no pond

                   ! pond decreasing or runoff starts
                   if (iok(kk)==1 .and. ns(kk)<1 .and. h0(kk)<h0max .and. h0(kk)+dy(kk,1)>h0max+dh0max) then ! start of runoff
                      fac(kk) = (h0max+half*dh0max-h0(kk))/dy(kk,1)
                      iok(kk) = 0
                   end if
                   if (iok(kk)==1 .and. ns(kk)<1 .and. h0(kk)+dy(kk,1)<h0min) then ! pond going
                      fac(kk) = -(h0(kk)-half*h0min)/dy(kk,1)
                      iok(kk) = 0
                   end if

				   	!reduce time-step at onset of melting
				!	do i=1, n  
				!		Tfreezing(kk) = Tfrz((S(kk,i)+dy(kk,i)*real(1-var(kk,i)%isat,r_2)),par(kk,i)%he,one/par(kk,i)%lam)
				!		var(kk,i)%Tfrz = Tfreezing(kk)
				!		if (var(kk,i)%iice.eq.1.and.(Tsoil(kk,i)+dTsoil(kk,i))>var(kk,i)%Tfrz+0.01) then	! check for onset of freezing in top soil layer	

				!			fac(kk)=min(max(half,accel(kk)*abs((var(kk,i)%Tfrz-Tsoil(kk,i)+0.005)/dTsoil(kk,i))),fac(kk))
				!			iok(kk)=0
				!		  exit
				!		endif
				!	enddo

                   ! repeat calculation if iice==1 and dthetamaxdT > 0.1 and no time step reduction

				   if (iok(kk)==1) then
						do i=1, n
							if ((.not.again_ice(kk,i))) then
                     
								if (var(kk,i)%iice.eq.1) then
									Tfreezing(kk) = Tfrz((S(kk,i)+dy(kk,i)*real(1-var(kk,i)%isat,r_2)),par(kk,i)%he,one/par(kk,i)%lam)
									var(kk,i)%Tfrz = Tfreezing(kk)
									tmp1d1(kk)    = ( thetalmax(min(Tsoil(kk,i)+dTsoil(kk,i), Tfreezing(kk)), &
									(S(kk,i)+dy(kk,i)*real(1-var(kk,i)%isat,r_2)), par(kk,i)%he, one/par(kk,i)%lam, &
									par(kk,i)%thre,par(kk,i)%the) - &
									thetalmax(Tsoil(kk,i), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, par(kk,i)%thre, par(kk,i)%the) &
									) / (min(Tsoil(kk,i)+dTsoil(kk,i),Tfreezing(kk))-Tsoil(kk,i))
									var(kk,i)%dthetaldT = half*(var(kk,i)%dthetaldT+tmp1d1(kk))
									nsteps_ice(kk,i) = nsteps_ice(kk,i) + 1

									tmp1d2(kk) = var(kk,i)%dthetaldT*(min(Tsoil(kk,i)+dTsoil(kk,i), Tfreezing(kk))-Tsoil(kk,i))
									if (i.eq.1) then
										tmp1d3(kk) = dy(kk,i)*real(1-var(kk,i)%isat,r_2)*par(kk,i)%thre*real(ns(kk))
									else
										tmp1d3(kk) = dy(kk,i)*real(1-var(kk,i)%isat,r_2)*par(kk,i)%thre
									endif

									tmp1d4(kk) = thetalmax(min(Tsoil(kk,i)+dTsoil(kk,i), Tfreezing(kk)), &
									(S(kk,i)+dy(kk,i)*real(1-var(kk,i)%isat,r_2)), par(kk,i)%he, one/par(kk,i)%lam, &
									par(kk,i)%thre,par(kk,i)%the)- &
									thetalmax(Tsoil(kk,i), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, par(kk,i)%thre, par(kk,i)%the)
									

									


									if (abs(tmp1d4(kk)-tmp1d2(kk))<tol_dthetaldT) then
										again_ice(kk,i) = .true.
									endif
									    

									if (nsteps_ice(kk,i).gt.nsteps_ice_max) then
										fac(kk) = 0.5
										nsteps_ice(kk,i) = 0
									endif
									var(kk,i)%csoileff = var(kk,i)%csoil + var(kk,i)%dthetaldT*rhow*lambdaf
									iok(kk) = 0
								endif
							endif
                      enddo
                      if (fac(kk).lt.1.0) then
                         again_ice(kk,1:n) = .false.  ! reset all again_ice if calc is to be repated with smaller time-step
                      endif
                   endif

                   ! reduce time step
                   if (iok(kk)==0) then
                      t(kk)      = t(kk)-dt(kk)
                      dt(kk)     = fac(kk)*dt(kk)
                      t(kk)      = t(kk)+dt(kk)
                      rsigdt(kk) = one/(sig(kk)*dt(kk))
                      nless(kk)  = nless(kk)+1 ! count step size reductions
					 ! write(339,"(i4,i4,i4,14f15.6)") irec, nsteps(kk),nless(kk), dt(kk), fac 
                   end if
                   if (var(kk,1)%isat/=0 .and. iflux(kk)==1 .and. var(kk,1)%phi<phip(kk) .and. &
                        var(kk,1)%phi+dy(kk,1)>phip(kk)) then
                      ! incipient (onset of) ponding - adjust state of saturated regions
                      t(kk)      = t(kk)-dt(kk)
                      dt(kk)     = dtmin
                      rsigdt(kk) = one/(sig(kk)*dt(kk))
                      again(kk)  = .true.
                      iok(kk)    = 0
                   end if
                end if
             end do ! while (iok==0) ----- end get and solve eqns

             !----- update unknowns
             ! i.e. update state variables to end of time step
             !      cumulate some surface fluxes
             ih0(kk) = 0
             if (.not. again(kk)) then
                ! get surface fluxes at sigma of the time step
                dwoff(kk)  = zero
                deltah0(kk) = zero
                precip(kk) = precip(kk) + qprec(kk)*dt(kk)
                if (ns(kk)<1) then ! pond
                   ! adjust pond height
                   h0sig(kk)   = h0(kk) + sig(kk)*dy(kk,0)
                   h0(kk)      = h0(kk)+dy(kk,0)
                   deltah0(kk) = dy(kk,0)
                   if (h0(kk)<zero .and. dy(kk,0)<zero) ih0(kk)=1 ! pond gone


                   ! note that fluxes required are q at sigma of time step
                   ! infiltration
                   dwinfil(kk)  = (q(kk,0)+sig(kk)*(qya(kk,0)*dy(kk,0)+qyb(kk,0)*dy(kk,1)+qTb(kk,0)*dTsoil(kk,1)))*dt(kk) &
                        - deltah0(kk) ! pond included in top soil layer
                   dwcol(kk)    = zero                   
                   do j=2, n
                      dwcol(kk) = sum(par(kk,1:n)%thre*dx(kk,1:n)*dy(kk,1:n)*(abs(var(kk,1:n)%isat-1)),1)
                   enddo

                   ! change in heat stored in soil column
                   do j=1,n

                      if (j.eq.1) then
                         deltaJ_latent_S(kk,j) = -dx(kk,j)*dy(kk,j)*par(kk,j)%thre*real(var(kk,j)%iice)* &
                              real(-var(kk,j)%isat+1)*rhow*lambdaf*real(ns(kk)) ! pond included in top soil layer
                         deltaJ_sensible_T(kk,j) = var(kk,j)%csoil*dx(kk,j)*dTsoil(kk,j) &
                              + h0(kk)*rhow*cswat*dTsoil(kk,j)        

                         if (advection==1) then
                            deltaJ_sensible_S(kk,j) = rhow*cswat*(Tsoil(kk,j)+Tzero)*dy(kk,j)*real(-var(kk,j)%iice+1) &
                                 + rhow*csice*(Tsoil(kk,j)+Tzero)*dy(kk,j)*real(var(kk,j)%iice)
                         endif
                      else
                         deltaJ_latent_S(kk,j) = -dx(kk,j)*dy(kk,j)*par(kk,j)%thre*real(var(kk,j)%iice)* &
                              real(-var(kk,j)%isat+1)*rhow*lambdaf
                         deltaJ_sensible_T(kk,j) = var(kk,j)%csoil*dx(kk,j)*dTsoil(kk,j) 

                         if (advection==1) then
                            deltaJ_sensible_S(kk,j) = rhow*cswat*(Tsoil(kk,j)+Tzero)*dx(kk,j)*dy(kk,j)*par(kk,j)%thre* &
                                 real(-var(kk,j)%isat+1)*real(-var(kk,j)%iice+1) &
                                 + rhow*csice*(Tsoil(kk,j)+Tzero)*dx(kk,j)*dy(kk,j)*par(kk,j)%thre &
                                 *real(-var(kk,j)%isat+1)*real(var(kk,j)%iice)
                         endif
                      endif
                      deltaJ_latent_T(kk,j) = var(kk,j)%dthetaldT*dx(kk,j)*dTsoil(kk,j)*rhow*lambdaf
                   end do
                   h0sig(kk) = h0(kk) + sig(kk)*dy(kk,1)
                   h0(kk) = h0(kk) + dy(kk,1)
                   deltah0(kk) = dy(kk,1)
                   if (h0(kk)<zero .and. dy(kk,1)<zero) then
                      ih0(kk)=1 ! pond gone
                   endif

                else   ! no pond


                   dwcol(kk) = sum(par(kk,1:n)%thre*dx(kk,1:n)*dy(kk,1:n)*(abs(var(kk,1:n)%isat-1)),1)
                   dwcol(kk)   = dwcol(kk) + plit(kk)%thre*dy(kk,0)*dxL(kk)
                   dwinfil(kk) = (q(kk,0) + &
                        sig(kk)*(qya(kk,0)*dy(kk,0)+qTa(kk,0)*dTsoil(kk,0)+qyb(kk,0)*dy(kk,1)+qTb(kk,0)*dTsoil(kk,1))) * dt(kk)
                   ! change in heat stored in soil column
                   do j=1,n

                      deltaJ_latent_S(kk,j) = -dx(kk,j)*dy(kk,j)*par(kk,j)%thre*real(var(kk,j)%iice)* &
                           real(-var(kk,j)%isat+1)*rhow*lambdaf
                      deltaJ_latent_T(kk,j) = var(kk,j)%dthetaldT*dx(kk,j)*dTsoil(kk,j)*rhow*lambdaf
                      if (j.eq.1) then
                         deltaJ_sensible_T(kk,j) = var(kk,j)%csoil*dx(kk,j)*dTsoil(kk,j) + h0(kk)*rhow*cswat*dTsoil(kk,j)   
                         if (advection==1) then
                            deltaJ_sensible_S(kk,j) =	rhow*cswat*(Tsoil(kk,j)+Tzero)*dx(kk,j)*dy(kk,j)*par(kk,j)%thre* &
                                 real(-var(kk,j)%isat+1)*real(-var(kk,j)%iice+1) &			
                                 + rhow*csice*(Tsoil(kk,j)+Tzero)*dx(kk,j)*dy(kk,j)*par(kk,j)%thre* &
                                 real(-var(kk,j)%isat+1)*real(var(kk,j)%iice)
                         endif
                      else
                         deltaJ_sensible_T(kk,j) = var(kk,j)%csoil*dx(kk,j)*dTsoil(kk,j) 
                         if (advection==1) then
                            deltaJ_sensible_S(kk,j) =   rhow*cswat*(Tsoil(kk,j)+Tzero)*dx(kk,j)*dy(kk,j)*par(kk,j)%thre* &
                                 real(-var(kk,j)%isat+1)*real(-var(kk,j)%iice+1) &					
                                 + rhow*csice*(Tsoil(kk,j)+Tzero)*dx(kk,j)*dy(kk,j)*par(kk,j)%thre*real(-var(kk,j)%isat+1)* &
                                 real(var(kk,j)%iice)
                         endif
                      endif
                   end do

                end if ! end no pond

                ! cumulate evaporation from top of soil column or top of litter/pond or top of snow pack
                select case (pondcase(kk))
                case(1)  ! no pond or litter
                   evap(kk)     = evap(kk) +qevap(kk)*dt(kk) - sig(kk)*(qyb(kk,0)*dy(kk,1)+qTb(kk,0)*DTsoil(kk,1))*dt(kk)
                   qevapsig(kk) = qevap(kk) - sig(kk)*(qyb(kk,0)*dy(kk,1)+qTb(kk,0)*DTsoil(kk,1))
                   Gcum(kk) = Gcum(kk)+(qhsig(kk,0)-qadvsig(kk,0))*dt(kk) 
                   Qadvcum(kk) = Qadvcum(kk) + qadvsig(kk,0)*dt(kk) - qadvsig(kk,n)*dt(kk)
                   lEcum(kk)    = lEcum(kk) + qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk)
                   Hcum(kk) = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-(qhsig(kk,0)-qadvsig(kk,0))*dt(kk)- &
                        qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk))
                   dwinfil(kk) = (qprec(kk)-qevapsig(kk))*dt(kk)
                case(2)  ! litter, no pond
                   evap(kk)     = evap(kk) - (qL(kk)-qprec(kk) + sig(kk)*(qybL(kk)*dy(kk,0)+qTbL(kk)*(dTsoil(kk,0) &
                        -deltaTa(kk))))*dt(kk)
                   qevapsig(kk) = - (qL(kk)-qprec(kk) + sig(kk)*(qybL(kk)*dy(kk,0)+qTbL(kk)*(dTsoil(kk,0)-deltaTa(kk))))
                   Gcum(kk) = Gcum(kk)+qhsig(kk,0)*dt(kk) 
                   lEcum(kk)    = lEcum(kk) + qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk)
                   Hcum(kk) = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-qhsig(kk,0)*dt(kk)- qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk))

                   dwinlit(kk)    = (qL(kk) + sig(kk)*(qybL(kk)*dy(kk,0)+qTbL(kk)*(dTsoil(kk,0)-deltaTa(kk))))*dt(kk)
                   dwinfil(kk) = qsig(kk,0)*dt(kk)
                case(3)! pond
				if (maxpond(kk)) then

					evap(kk)     = evap(kk) +qevap(kk)*dt(kk)
					qevapsig(kk) = qevap(kk)
					dwinfil(kk)  = (q(kk,0)+sig(kk)*(qya(kk,0)*dy(kk,0)+qyb(kk,0)*dy(kk,1)+qTb(kk,0)*dTsoil(kk,1)))*dt(kk) - &
                        deltah0(kk) ! pond included in top soil layer
					dwoff(kk)    = (qprec(kk)-qevapsig(kk))*dt(kk) -dwinfil(kk)
					if (dwoff(kk)<0) then
					!write(*,*) 'Negative runoff', dwoff, irec

					deltah0 = dwoff
					h0 = h0+deltah0
					dwoff = zero

					! adjust sensible heat content
					endif
				else
                   evap(kk)     = evap(kk) +qevap(kk)*dt(kk) - sig(kk)*(qyb(kk,0)*dy(kk,1)+qTb(kk,0)*DTsoil(kk,1))*dt(kk)
                   qevapsig(kk) = qevap(kk) - sig(kk)*(qyb(kk,0)*dy(kk,1)+qTb(kk,0)*DTsoil(kk,1))
                   dwinfil(kk)  = (q(kk,0)+sig(kk)*(qya(kk,0)*dy(kk,0)+qyb(kk,0)*dy(kk,1)+qTb(kk,0)*dTsoil(kk,1)))*dt(kk) - &
                        deltah0(kk) ! pond included in top soil layer
				
				 endif
				   Gcum(kk) = Gcum(kk)+(qhsig(kk,0)-qadvsig(kk,0))*dt(kk) 
                   Qadvcum(kk) = Qadvcum(kk) + qadvsig(kk,0)*dt(kk) - qadvsig(kk,n)*dt(kk)
                   lEcum(kk)    = lEcum(kk) + qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk)
                   Hcum(kk) = Hcum(kk) + (vmet(kk)%Rn*dt(kk)-(qhsig(kk,0)-qadvsig(kk,0))*dt(kk)- &
                        qevapsig(kk)*thousand*var(kk,1)%lambdav*dt(kk))

				
                   
                end select ! pondcase


                dwdrainage(kk)     = q(kk,n)*dt(kk) +sig(kk)*dt(kk)*(qya(kk,n)*dy(kk,n)+qTa(kk,n)*dTsoil(kk,n))

                if (botbc=="aquifer" .and. v_aquifer(kk)%isat==0) then
                   dwdischarge(kk) = v_aquifer(kk)%discharge*dt(kk)
                else
                   dwdischarge(kk) = dwdrainage(kk)
                endif

                drexcol(kk) = sum(qex(kk,:),1)*dt(kk)
                qrunoff(kk) = dwoff(kk)/dt(kk)

                ! cumulative fluxes
                infil(kk)     = infil(kk) + dwinfil(kk)
                inlit(kk)     = inlit(kk) + dwinlit(kk)
                wcol(kk)      = wcol(kk) + dwcol(kk)
                drainage(kk)  = drainage(kk) + dwdrainage(kk)
                discharge(kk) = discharge(kk) + dwdischarge(kk)
                runoff(kk)    = runoff(kk) + dwoff(kk)
                rexcol(kk)    = rexcol(kk) + drexcol(kk)
                ponding(kk)   = ponding(kk) + deltah0(kk)
				Qadvcum(kk) = Qadvcum(kk) - sum(qex(kk,:)*(Tsoil(kk,:)+Tzero),1)*dt(kk)*rhow*cswat

                ! evaluate soil fluxes at sigma of time step for use in isotope routine
                if (ns(kk)<1) then
                   qsig(kk,0)  = q(kk,0)  + sig(kk)*(qyb(kk,0)*dy(kk,1)  + qya(kk,0)*dy(kk,0)  + qTb(kk,0)*dTsoil(kk,1) &
                        + qTa(kk,0)*dTsoil(kk,0))
                   qhsig(kk,0) = qh(kk,0) + sig(kk)*(qhyb(kk,0)*dy(kk,1) + qhya(kk,0)*dy(kk,0) + qhTb(kk,0)*dTsoil(kk,1) &
                        + qhTa(kk,0)*dTsoil(kk,0))
                   qlsig(kk,0) = qsig(kk,0)
                   qvsig(kk,0) = zero
                else
                   qsig(kk,0)  = q(kk,0)  + sig(kk)*(qyb(kk,0)*dy(kk,1)  + qya(kk,0)*dy(kk,0)  + qTb(kk,0)*dTsoil(kk,1) &
                        + qTa(kk,0)*dTsoil(kk,0))
                   qhsig(kk,0) = qh(kk,0) + sig(kk)*(qhyb(kk,0)*dy(kk,1) + qhya(kk,0)*dy(kk,0) + qhTb(kk,0)*dTsoil(kk,1) &
                        + qhTa(kk,0)*dTsoil(kk,0))
                   qvsig(kk,0) = qsig(kk,0)
                   qlsig(kk,0) = zero
                endif
                qsig(kk,1:n-1)   = q(kk,1:n-1) + sig(kk)*(qya(kk,1:n-1)*dy(kk,1:n-1) + qyb(kk,1:n-1)*dy(kk,2:n) &
                     + qTa(kk,1:n-1)*dTsoil(kk,1:n-1) + qTb(kk,1:n-1)*dTsoil(kk,2:n))
                qsig(kk,n)       = q(kk,n) + sig(kk)*(qya(kk,n)*dy(kk,n) + qTa(kk,n)*dTsoil(kk,n))
                qhsig(kk,1:n-1)  = qh(kk,1:n-1) + sig(kk)*(qhya(kk,1:n-1)*dy(kk,1:n-1) + qhyb(kk,1:n-1)*dy(kk,2:n) &
                     + qhTa(kk,1:n-1)*dTsoil(kk,1:n-1) + qhTb(kk,1:n-1)*dTsoil(kk,2:n))
                qhsig(kk,n)      = qh(kk,n) + sig(kk)*(qhya(kk,n)*dy(kk,n) + qhTa(kk,n)*dTsoil(kk,n))
                qvsig(kk,1:n-1)  = qv(kk,1:n-1) + sig(kk)*(qvya(kk,1:n-1)*dy(kk,1:n-1) + qvyb(kk,1:n-1)*dy(kk,2:n) &
                     + qTa(kk,1:n-1)*dTsoil(kk,1:n-1) + qTb(kk,1:n-1)*dTsoil(kk,2:n))
                qvsig(kk,n)      = zero
                qvTsig(kk,1:n-1) = qvT(kk,1:n-1) + sig(kk)*(qTa(kk,1:n-1)*dTsoil(kk,1:n-1) + qTb(kk,1:n-1)*dTsoil(kk,2:n))
                qvTsig(kk,n)     = zero
                qlsig(kk,1:n-1)  = qsig(kk,1:n-1) - qvsig(kk,1:n-1)
                qlsig(kk,n)      = qsig(kk,n)


                if (botbc=="aquifer") then ! update aquifer props
                   if (v_aquifer(kk)%isat==0) then
                      if (v_aquifer(kk)%WA+(dwdrainage(kk)-dwdischarge(kk))*dt(kk) > &
                           (v_aquifer(kk)%zzero-v_aquifer(kk)%zsoil)*v_aquifer(kk)%Sy) then
                         v_aquifer(kk)%isat = 1  ! aquifer saturated
                         S(kk,n) = S(kk,n) + (-(v_aquifer(kk)%WA+(dwdrainage(kk)-dwdischarge(kk))*dt(kk)) &
                              +(v_aquifer(kk)%zzero-v_aquifer(kk)%zsoil)*v_aquifer(kk)%Sy)/(dx(kk,n)*par(kk,n)%thre)
                      endif
                      v_aquifer(kk)%WA        = min(v_aquifer(kk)%WA+(dwdrainage(kk)-dwdischarge(kk))*dt(kk), &
                           (v_aquifer(kk)%zzero-v_aquifer(kk)%zsoil)*v_aquifer(kk)%Sy)
                      v_aquifer(kk)%zdelta    = v_aquifer(kk)%zzero - v_aquifer(kk)%Wa/v_aquifer(kk)%Sy
                      ! new discharge rate
                      v_aquifer(kk)%discharge = v_aquifer(kk)%Rsmax*exp(-v_aquifer(kk)%f*v_aquifer(kk)%zdelta)
                   elseif (v_aquifer(kk)%isat==1) then
                      ! check for desat of aquifer
                      if (dwdrainage(kk) < v_aquifer(kk)%Rsmax*exp(-v_aquifer(kk)%f*v_aquifer(kk)%zsoil)) then
                         v_aquifer(kk)%isat      = 0
                         v_aquifer(kk)%zdelta    = v_aquifer(kk)%zsoil
                         ! new discharge rate
                         v_aquifer(kk)%discharge = v_aquifer(kk)%Rsmax*exp(-v_aquifer(kk)%f*v_aquifer(kk)%zdelta)
                      endif
                   endif
                   zdelta(kk) = v_aquifer(kk)%zdelta
                endif       ! end update aquifer props

                if (botbc=="constant head") then
                   drn(kk) = drn(kk)+(q(kk,n)+sig(kk)*qya(kk,n)*dy(kk,n))*dt(kk)
                else
                   drn(kk) = drn(kk)+(q(kk,n)+sig(kk)*qya(kk,n)*dy(kk,n))*dt(kk)
                end if

                if (present(wex)) then
                   ! where (var%isat==0)!
                   !wex(kk,1:n) = wex(kk,1:n)+(qex(kk,1:n)+spread(sig(kk),2,n)*qexd(kk,1:n)*dy(kk,1:n))*spread(dt(kk),2,n)
                   wex(kk,1:n) = wex(kk,1:n)+(qex(kk,1:n)+spread(sig(kk),1,n)*qexd(kk,1:n)*dy(kk,1:n))*spread(dt(kk),1,n)
                end if

                if (litter .and. ns(kk)==1 .and. (.not. maxpond(kk))) then ! adjust litter moisture content if no ponding
                   SL(kk)      = SL(kk) + dy(kk,0)
                   deltaSL(kk) = dy(kk,0)
                else
                   deltaSL(kk) = zero
                endif

                TL(kk) = TL(kk) + dTsoil(kk,0)   ! Tlitter assumed the same as pond T
                if (SL(kk) >= one) vlit(kk)%isat = 1   ! check for litter saturation

                csoil(kk,:) = var(kk,1:n)%csoil


             endif ! .not. again

             ! update variables (S,T) to end of time step
             ! update variables to sig for use in isotope routine
             do i=1, n
                if (var(kk,i)%isat==0) then
                   if (.not.again(kk)) then
                      deltaS(kk,i) = dy(kk,i)  !required for isotope subroutine
                      Ssig(kk,i)   = S(kk,i) + two*sig(kk)*dy(kk,i)  !required for isotope subroutine
                      if (i.eq.1) then
                         S(kk,i) = S(kk,i) + dy(kk,i)*real(ns(kk)) ! pond included in top soil layer
                      else
                         S(kk,i)      = S(kk,i)+dy(kk,i)
                      endif
                      Tsoil(kk,i)  = Tsoil(kk,i) + dTsoil(kk,i)  ! update soil temperature
                      if (S(kk,i)>one .and. dy(kk,i)>zero) then ! saturation of layer
                         var(kk,i)%isat = 1
                         var(kk,i)%K    = par(kk,i)%Ke
                         var(kk,i)%phi  = par(kk,i)%phie
                         if (var(kk,i)%iice==1 ) then
                            var(kk,i)%K    = var(kk,i)%Ksat
                            var(kk,i)%phi  = var(kk,i)%phie
                         endif
                      end if
                   end if
                else
                   if (.not.again(kk)) Tsoil(kk,i)  = Tsoil(kk,i) + dTsoil(kk,i)
                   deltaS(kk,i)  = zero    !required for isotope subroutine
                   Ssig(kk,i)    = S(kk,i)    !required for isotope subroutine
                   if (i.eq.1) then
                      var(kk,i)%phi = var(kk,i)%phi + dy(kk,i)*real(ns(kk)) ! pond included in top soil layer
                   else
                      var(kk,i)%phi = var(kk,i)%phi + dy(kk,i)
                   endif

                   if (i==1 .and. ih0(kk)/=0 .and. var(kk,i)%phi>=par(kk,i)%phie) var(kk,i)%phi = zero ! pond gone
                   if (i==1 .and. ih0(kk)/=0 .and. var(kk,i)%phi>=var(kk,i)%phie) var(kk,i)%phi = zero ! pond gone (frozen soil)
                   if (var(kk,i)%phi<var(kk,i)%phie .and. var(kk,i)%iice==0) then ! desaturation of layer
                      var(kk,i)%isat = 0
                      var(kk,i)%K    = par(kk,i)%Ke
                      var(kk,i)%phi  = par(kk,i)%phie
                      var(kk,i)%KS   = par(kk,i)%KSe
                      var(kk,i)%phiS = par(kk,i)%phiSe
                   elseif (var(kk,i)%phi<var(kk,i)%phie.and. var(kk,i)%iice==1) then
                      var(kk,i)%isat = 0
                      var(kk,i)%K    = var(kk,i)%Ksat
                      var(kk,i)%phi  = var(kk,i)%phie
                      var(kk,i)%KS   = var(kk,i)%KS
                      var(kk,i)%phiS = var(kk,i)%phiS
                   endif
                end if
                
				if (.not. again(kk)) then
				! calculate freezing point temperature at new S
                Tfreezing(kk) = Tfrz(S(kk,i),par(kk,i)%he,one/par(kk,i)%lam)
                ! check for onset of freezing and adjust (increase) temperature to account for latent heat
                ! release by ice formation

                if (Tsoil(kk,i) < Tfreezing(kk) .and. var(kk,i)%iice==0) then
                   dtdT(kk)      = dthetalmaxdTh(Tfreezing(kk), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, &
                        par(kk,i)%thre,par(kk,i)%the)
				   tmp1d3(kk) = 1000.
				    k = 1
                   if (i.eq.1.and.h0(kk)>zero) then
                      h0(kk) = h0(kk) - deltah0(kk)

                      do while ((k.lt.nsteps_ice_max).and.(tmp1d3(kk).gt.tol_dthetaldT))
					  tmp1d1(kk) = ((var(kk,i)%csoil*Tsoil(kk,i) + &
                           rhow*lambdaf*dtdT(kk)*Tfreezing(kk))*dx(kk,i)+h0(kk)* &
                           (Tfreezing(kk)*rhow*csice+cswat*rhow*(Tsoil(kk,i)-Tfreezing(kk)))) &
                           / ((var(kk,1)%csoil + rhow*lambdaf*dtdT(kk))*dx(kk,i)+ csice*rhow*h0(kk))

					  !tmp1d2(kk) = dthetalmaxdT(tmp1d1(kk), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, &  !  estimate of dthetaldT at coreected temperature
                       ! par(kk,i)%thre,par(kk,i)%the)

					  !dtdT(kk) = (dtdT(kk) + tmp1d2(kk))/2.                ! updated estiamte of mean dthetaldT

                      dtdT(kk) = 	(thetalmax(tmp1d1(kk), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, par(kk,i)%thre, par(kk,i)%the) &
									 - (S(kk,i)*par(kk,i)%thre + (par(kk,i)%the-par(kk,i)%thre)))/ &
									(tmp1d1(kk) - Tfreezing(kk))

					  tmp1d3(kk) = abs((tmp1d1(kk)- Tfreezing(kk))*dtdT(kk))

					  k = k+1
					 
					 enddo  ! end of do while

                      dTsoil(kk,i) = tmp1d1(kk) - (Tsoil(kk,i)-dTsoil(kk,i))
                      Tsoil(kk,i) = tmp1d1(kk)
                      deltaJ_sensible_T(kk,i) = (var(kk,i)%csoil*dx(kk,i)+ cswat*rhow*h0(kk))*dTsoil(kk,i)
                      h0(kk) = h0(kk) + deltah0(kk)
                   else



					 do while ((k.lt.nsteps_ice_max).and.(tmp1d3(kk).gt.tol_dthetaldT))
					  tmp1d1(kk) = (var(kk,i)%csoil*Tsoil(kk,i) + rhow*lambdaf*dtdT(kk)*Tfreezing(kk)) &  ! zeroth estimate of corrected temperature
                           / (var(kk,i)%csoil + rhow*lambdaf*dtdT(kk))

					  
					  !tmp1d2(kk) = dthetalmaxdT(tmp1d1(kk), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, &  !  estimate of dthetaldT at coreected temperature
                      !  par(kk,i)%thre,par(kk,i)%the)

					  dtdT(kk) = 	(thetalmax(tmp1d1(kk), S(kk,i), par(kk,i)%he, one/par(kk,i)%lam, par(kk,i)%thre, par(kk,i)%the) &
									 - (S(kk,i)*par(kk,i)%thre + (par(kk,i)%the-par(kk,i)%thre)))/ &
									(tmp1d1(kk) - Tfreezing(kk))

					  tmp1d3(kk) = abs((tmp1d1(kk)- Tfreezing(kk))*dtdT(kk))

					  k = k+1
					 
					 enddo  ! end of do while
					 !tmp1d1(kk) = (var(kk,i)%csoil*Tsoil(kk,i) + rhow*lambdaf*dtdT(kk)*Tfreezing(kk)) &  ! first estimate of corrected temperature
                      !     / (var(kk,i)%csoil + rhow*lambdaf*dtdT(kk))
					  
					  						   
					   dTsoil(kk,i) = tmp1d1(kk)- (Tsoil(kk,i)-dTsoil(kk,i))
					  
					   Tsoil(kk,i) = tmp1d1(kk)
				

					  
					  			
                      deltaJ_sensible_T(kk,i) = var(kk,i)%csoil*dTsoil(kk,i)*dx(kk,i)
                   endif
                   deltaJ_latent_T(kk,i) = deltaJ_latent_T(kk,i) - dtdT(kk)*(Tfreezing(kk)- &
                        Tsoil(kk,i))*dx(kk,i)*rhow*lambdaf
                   var(kk,i)%csoileff = var(kk,i)%csoil + dtdT(kk)*rhow*lambdaf
                endif
                ! check for onset of thawing and decrease temperature to account for latent heat required to melt ice 
                if (Tsoil(kk,i)>Tfreezing(kk) .and. var(kk,i)%iice==1) then
                   ! Correct for "overthawing". This extra energy is used to heat the soil even more
                   ! The correction comes from equating the energy balances before and after correction.
                   if (i.eq.1.and.h0(kk)>zero) then
                      h0(kk) = h0(kk) - deltah0(kk)
                      deltaJ_latent_T(kk,i) = var(kk,i)%dthetaldT*(Tfreezing(kk)-(Tsoil(kk,i) &
                           -dTsoil(kk,i)))*dx(kk,i)*rhow*lambdaf

                      tmp1d1(kk) =  (Tsoil(kk,i)*var(kk,i)%csoil*dx(kk,i) + &
                           rhow*lambdaf*var(kk,i)%dthetaldT*(Tsoil(kk,i)-Tfreezing(kk))*dx(kk,i) &
                           + rhow*cswat*h0(kk)*Tfreezing(kk) + rhow*csice*h0(kk)*(Tsoil(kk,i)-Tfreezing(kk))) &
                           /(var(kk,i)%csoil*dx(kk,i) + cswat*rhow*h0(kk))	
                      dTsoil(kk,i) = tmp1d1(kk) - (Tsoil(kk,i)-dTsoil(kk,i))
                      Tsoil(kk,i) = tmp1d1(kk)
                      deltaJ_sensible_T(kk,i) = (var(kk,i)%csoil*dx(kk,i)+ cswat*rhow*h0(kk))*dTsoil(kk,i)
                      h0(kk) = h0(kk) + deltah0(kk)

                   else
                      deltaJ_latent_T(kk,i) = var(kk,i)%dthetaldT*(Tfreezing(kk)-(Tsoil(kk,i)-dTsoil(kk,i)))*dx(kk,i)*rhow*lambdaf
                      dTsoil(kk,i) = Tsoil(kk,i) + rhow*lambdaf/var(kk,i)%csoil*var(kk,i)%dthetaldT*(Tsoil(kk,i)-Tfreezing(kk)) &
                           - (Tsoil(kk,i)-dTsoil(kk,i))
                      Tsoil(kk,i) = Tsoil(kk,i) + rhow*lambdaf/var(kk,i)%csoil*var(kk,i)%dthetaldT*(Tsoil(kk,i)-Tfreezing(kk))
                      deltaJ_sensible_T(kk,i) = var(kk,i)%csoil*dTsoil(kk,i)*dx(kk,i)
                   endif
                endif
				endif
             end do ! i=1, n => update variables S,T

			if (.not. again(kk)) then
              ! update T0, rh0 and var, needed for isotopes
                cv0(kk,1:n)   = var(kk,1:n)%cv
                S0(kk,1:n)    = S(kk,1:n) - deltaS(kk,1:n)
                Sliq0(kk,1:n) = (S0(kk,1:n) - cv0(kk,1:n))/(one-cv0(kk,1:n))
                isave(kk,1:n) = var(kk,1:n)%isat
                call hyofS(S(kk,1:n), Tsoil(kk,1:n), par(kk,1:n), var(kk,1:n))
                var(kk,1:n)%isat  = isave(kk,1:n)
                deltacv(kk,1:n)   = var(kk,1:n)%cv - cv0(kk,1:n)
                Sliq(kk,1:n)      = (S(kk,1:n) - var(kk,1:n)%cv)/(one-var(kk,1:n)%cv)
                deltaSliq(kk,1:n) = Sliq(kk,1:n) - Sliq0(kk,1:n)
                if (littercase==1) then
                   SL0(kk)    = SL(kk) - deltaSL(kk)
                   cvL0(kk)   = vlit(kk)%cv
                   SLliq0(kk) = (SL0(kk) - cvL0(kk))/(one-cvL0(kk))
                   call litter_props(Sl(kk), Tl(kk), vlit(kk), plit(kk), h0(kk))
                   deltacvL(kk)   = vlit(kk)%cv - cvL0(kk)
                   SLliq(kk)      = (SL(kk) - vlit(kk)%cv)/(one-vlit(kk)%cv)
                   deltaSLliq(kk) = SLliq(kk) - SLliq0(kk)
                endif
             endif

             iflux(kk) = iflux(kk)+1
             !if (.not. again(kk)) exit
          end do ! do while (again(kk)) ! iflux loop

          ! just checking
          if (dt(kk) <= dtmin) then
             write(2,*) "solve: time step =  (irec)", dt(kk), irec
             !stop
          end if



          if (isotopologue /= 0) then

             ! remove negative h0 (optional)
             if (h0(kk)<zero .and. var(kk,1)%isat==0) then
                infil(kk)     = infil(kk)+h0(kk)
                ciso(kk,1)   = (ciso(kk,1)*S(kk,1)*par(kk,1)%thre*dx(kk,1) + h0(kk)*ciso(kk,0)) / &
                     (S(kk,1)*par(kk,1)%thre*dx(kk,1) + h0(kk))
                S(kk,1)      = S(kk,1)+h0(kk)/(par(kk,1)%thre*dx(kk,1))
                deltaS(kk,1) = deltaS(kk,1) + h0(kk)/(par(kk,1)%thre*dx(kk,1))  ! required for isotope
                Ssig(kk,1)   = S(kk,1) + (two*sig(kk) - one) *deltaS(kk,1)           ! required for isotope
                qsig(kk,0)   = qsig(kk,0) + h0(kk)/dt(kk)
                qlsig(kk,0)  = qlsig(kk,0) + h0(kk)/dt(kk)
                dy(kk,0)     = dy(kk,0) - h0(kk)
                h0(kk)        = zero
                h0sig(kk)     = h0(kk) + (two*sig(kk) - one) *dy(kk,0)
                deltah0(kk)   = dy(kk,0)

                S0(kk,1)     = S(kk,1) - deltaS(kk,1)
                Sliq0(kk,1)  = (S0(kk,1) - cv0(kk,1))/(one-cv0(kk,1))

                isave(kk,1:n)    = var(kk,1:n)%isat
                call hyofS(S(kk,1:n), Tsoil(kk,1:n), par(kk,1:n), var(kk,1:n))
                var(kk,1:n)%isat = isave(kk,1:n)
                rh0(kk)          = min(rh0(kk), var(kk,1)%rh)
                deltacv(kk,1)    = var(kk,1)%cv - cv0(kk,1)
                Sliq(kk,1)       = (S(kk,1) - var(kk,1)%cv)/(one-var(kk,1)%cv)
                deltaSliq(kk,1)  = Sliq(kk,1) - Sliq0(kk,1)
             end if

             tmp_thetasat(kk,1:n)   = par(kk,1:n)%thre
             tmp_tortuosity(kk,1:n) = par(kk,1:n)%tortuosity
             tmp_thetar(kk,1:n)     = par(kk,1:n)%the - par(kk,1:n)%thre
             qvsig(kk,0)  = qv0(kk)
             qvTsig(kk,0) = qvT0(kk)
             qlsig(kk,0)  = ql0(kk)

             call isotope_vap(isotopologue, testcase, litter, n, &
                  ns(kk), dx(kk,1:n), dz(kk,1:n-1), sig(kk), dt(kk), dxL(kk), maxpond(kk), &
                  Tsoil(kk,1:n), dTsoil(kk,0:n), Sliq(kk,1:n), deltaSliq(kk,1:n), &
                  Tsurface(kk), TL(kk), T0(kk), h0(kk), deltah0(kk), SLliq(kk), deltaSLliq(kk), &
                  qsig(kk,0:n), qlsig(kk,0:n), qvsig(kk,0:n), &
                  qprec(kk), qevapsig(kk), qrunoff(kk), qex(kk,1:n), qd(kk), &
                  var(kk,1:n), tmp_thetasat(kk,1:n), tmp_thetar(kk,1:n), tmp_tortuosity(kk,1:n), &
                  deltacv(kk,1:n), vmet(kk)%ra, vmet(kk)%rs, vlit(kk), vmet(kk)%cva, vmet(kk)%civa, &
                  plit(kk)%the, deltacvL(kk), &
                  cprec(kk), cali(kk), &
                  ql0(kk), qv0(kk), &
                  ciso(kk,0:n), ciso0(kk), cisoL(kk), cisos(kk), &
                  qiso_in(kk), qiso_out(kk), qiso_evap(kk), qiso_trans(kk), &
                  qiso_liq_adv(kk,1:n), qiso_vap_adv(kk,1:n), qiso_liq_diff(kk,1:n-1), qiso_vap_diff(kk,1:n-1))

             qiso_evap_cum(kk)  = qiso_evap_cum(kk)  + qiso_evap(kk)*dt(kk)
             qiso_trans_cum(kk) = qiso_trans_cum(kk) + qiso_trans(kk)*dt(kk)

          else                                                  ! (if iso==0)
             ! remove negative h0 (optional)
             if (h0(kk)<zero .and. var(kk,1)%isat==0) then
                infil(kk)     = infil(kk) + h0(kk)
                S(kk,1)       = S(kk,1)  + h0(kk)/(par(kk,1)%thre*dx(kk,1))
                dy(kk,1)      = dy(kk,1) - h0(kk)
                h0(kk)        = zero
                deltah0(kk)   = dy(kk,1)
                ! for frozen soil, adjust change in ice storage and temperature
                tmp1d2(kk) = deltaJ_latent_T(kk,1)+deltaJ_sensible_T(kk,1) + deltaJ_latent_S(kk,1) ! energy change in top soil layer
                deltaJ_latent_S(kk,1) = dx(kk,1)*deltaS(kk,1)*par(kk,1)%thre*real(var(kk,1)%iice)*rhow*lambdaf
                tmp1d1(kk) = Tsoil(kk,1)-dTsoil(kk,1)  ! old soil temp

                Tsoil(kk,1) = (tmp1d2(kk) - deltaJ_latent_S(kk,1) + var(kk,1)%csoileff*dx(kk,1)*tmp1d1(kk))/ &
                     (var(kk,1)%csoileff*dx(kk,1))
                dTsoil(kk,1) = Tsoil(kk,1)-tmp1d1(kk)
                deltaJ_sensible_T(kk,1) = dTsoil(kk,1)*dx(kk,1)*var(kk,1)%csoil
                deltaJ_latent_T(kk,1) = dTsoil(kk,1)*dx(kk,1)*(var(kk,1)%csoileff - var(kk,1)%csoil)

                isave(kk,1:n) = var(kk,1:n)%isat
                call hyofS(S(kk,1:n), Tsoil(kk,1:n), par(kk,1:n), var(kk,1:n))
                var(kk,1:n)%isat = isave(kk,1:n)
                rh0(kk)          = min(rh0(kk),var(kk,1)%rh)
                ns(kk)           = 1
             end if
          endif ! isotopologue/=0

          nsteps(kk)        = nsteps(kk) + 1
          delthetai(kk,1:n) = (var(kk,1:n)%thetai-thetai(kk,1:n)) * dx(kk,1:n)*rhoi/rhow 
		  delthetai_col(kk) = sum(delthetai(kk,1:n))
          thetai(kk,1:n)    = var(kk,1:n)%thetai
		  
		  delJsensible(kk,1) = (var(kk,1)%csoil* dx(kk,1)+h0(kk)*cswat)*(Tsoil(kk,1)+Tzero)  - Jsensible(kk,1)
		  Jsensible(kk,1) = (var(kk,1)%csoil* dx(kk,1)+h0(kk)*cswat)*(Tsoil(kk,1)+Tzero)
		  delJsensible(kk,2:n) = var(kk,2:n)%csoil*(Tsoil(kk,2:n)+Tzero)* dx(kk,2:n) - Jsensible(kk,2:n)
		  Jsensible(kk,2:n) = var(kk,2:n)%csoil*(Tsoil(kk,2:n)+Tzero)* dx(kk,2:n)
		  delJsensible_col(kk) = sum(delJsensible(kk,1:n))

          deltaTa(kk)       = zero
          init(kk)          = .false.


          ! change in heat stored in soil column
          dJcol_latent_S(kk) = sum(deltaJ_latent_S(kk,1:n))
          Jcol_latent_S(kk) = Jcol_latent_S(kk) + dJcol_latent_S(kk)
          dJcol_latent_T(kk) = sum(deltaJ_latent_T(kk,1:n))
          Jcol_latent_T(kk) = Jcol_latent_T(kk) + dJcol_latent_T(kk)
          dJcol_sensible(kk) = sum(deltaJ_sensible_T(kk,1:n)) + sum(deltaJ_sensible_S(kk,1:n))
          Jcol_sensible(kk) = Jcol_sensible(kk) + dJcol_sensible(kk)
          Jcol(kk) = Jcol_sensible(kk) + Jcol_latent_S(kk) + Jcol_latent_T(kk)
          dicecol_S(kk) = -dJcol_latent_S(kk)/rhow/lambdaf
          dicecol_T(kk) = -dJcol_latent_T(kk)/rhow/lambdaf
		  tmp2d1(kk,1:n) = -(deltaJ_latent_S(kk,1:n)+deltaJ_latent_T(kk,1:n))/rhow/lambdaf

          deltaice_cum_S(kk) = -Jcol_latent_S(kk)/rhow/lambdaf
          deltaice_cum_T(kk) = -Jcol_latent_T(kk)/rhow/lambdaf
		  

       end do ! while (t<tfin)
    end do ! kk=1, mp

    ! get heads if required
    if (present(heads)) then
       isave    = var%isat
       !MC! strange in original code: all calc with bottom n except Tsoil(1)
       !MC! Tbot from above
       call hyofS(Sbot, Tbot, par, var)
       var%isat = isave
       heads    = var%h
       where (S(:,:) >= one) heads(:,:) = par(:,:)%he + (var(:,:)%phi-par(:,:)%phie)/par(:,:)%Ke
    end if

    ! write(*,*) 'not15: ', h0
    ! write(*,*) 'not16: ', S
    ! write(*,*) 'not17: ', thetai
    ! write(*,*) 'not18: ', Tsoil
    ! write(*,*) 'not19: ', evap
    ! write(*,*) 'not20: ', qh
    ! write(*,*) 'not21: ', Hcum
    ! write(*,*) 'not22: ', LEcum
    ! write(*,*) 'not23: ', Gcum
    ! write(*,*) 'not24: ', fwsoil

  END SUBROUTINE solve

  !**********************************************************************************************************************

  ! SUBROUTINE solute(ti,tf,thi,thf,win,cin,n,ns,dx,jt,dsmmax,sm,sdrn,nssteps,c, &
  !      isosub)

  !   USE sli_utils, ONLY: dis, isotype, bd, isopar

  !   IMPLICIT NONE

  !   INTEGER(i_d),INTENT(IN)::n,ns,jt(n)
  !   REAL(r_2),INTENT(IN)::ti,tf,thi(n),thf(n),win,cin(ns),dx(n),dsmmax
  !   INTEGER(i_d),INTENT(INOUT)::nssteps(ns)
  !   REAL(r_2),INTENT(INOUT)::sm(n,ns),sdrn(ns),c(n,ns)
  !   OPTIONAL::isosub
  !   INTERFACE
  !      SUBROUTINE isosub(iso,c,p,f,fc)
  !        USE sli_numbers, ONLY: r_2
  !        CHARACTER(LEN=2),INTENT(IN)::iso
  !        REAL(r_2),INTENT(IN)::c
  !        REAL(r_2),DIMENSION(:),INTENT(INOUT)::p
  !        REAL(r_2),INTENT(OUT)::f, fc
  !      END SUBROUTINE isosub
  !   END INTERFACE
  !   ! Solves the ADE from time ti to tf. Diffusion of solute ignored - dispersion
  !   ! coeff = dispersivity * abs(pore water velocity).
  !   ! Definitions of arguments:
  !   ! Required args:
  !   ! ti   - start time (h).
  !   ! tf   - finish time.
  !   ! thi(1:n)  - initial layer water contents.
  !   ! thf(1:n)  - final layer water contents.
  !   ! win   - water in at top of profile.
  !   ! cin(1:ns)  - solute concn in win.
  !   ! n    - no. of soil layers.
  !   ! ns   - no. of solutes.
  !   ! dx(1:n)  - layer thicknesses.
  !   ! jt(1:n)  - layer soil type nos.
  !   ! dsmmax(1:ns) - max change in sm of any layer to aim for each time step;
  !   !      controls time step size.
  !   ! sm(1:n,1:ns) - layer masses of solute per cc.
  !   ! sdrn(1:ns) - cumulative solute drainage.
  !   ! nssteps(1:ns) - cumulative no. of time steps for ADE soln.
  !   ! Optional args:
  !   ! isosub  - subroutine to get adsorbed solute (units/g soil) from concn
  !   !      in soil water according to chosen isotherm code.
  !   !      Arguments: iso - 2 character code; c - concn in soil water;
  !   !      p(:) - isotherm parameters; f - adsorbed mass/g soil;
  !   !      fc - deriv of f wrt c (slope of isotherm curve).
  !   INTEGER(i_d),PARAMETER::itmax=20 ! max iterations for finding c from sm
  !   REAL(r_2),PARAMETER::eps=0.00001 ! for stopping
  !   INTEGER(i_d)::i,it,j,k
  !   REAL(r_2)::dc,dm,dmax,dt,dz(n-1),f,fc,r,rsig,rsigdt,sig,sigdt,t,tfin,th,v1,v2
  !   REAL(r_2),DIMENSION(n-1)::coef1,coef2
  !   REAL(r_2),DIMENSION(n)::csm,tht
  !   REAL(r_2),DIMENSION(0:n)::aa,bb,cc,dd,dy,q,qw,qya,qyb
  !   INTEGER(i_d) :: info

  !   sig=half
  !   rsig=one/sig
  !   tfin=tf
  !   dz=half*(dx(1:n-1)+dx(2:n))
  !   !get average water fluxes
  !   r=one/(tf-ti)
  !   qw(0)=r*win
  !   tht=r*(thf-thi)
  !   do i=1,n
  !      qw(i)=qw(i-1)-dx(i)*tht(i)
  !   end do
  !   !get constant coefficients
  !   do i=1,n-1
  !      v1=half*qw(i)
  !      v2=half*(dis(jt(i))+dis(jt(i+1)))*abs(qw(i))/dz(i)
  !      coef1(i)=v1+v2
  !      coef2(i)=v1-v2
  !   end do
  !   do j=1,ns
  !      t=ti
  !      if (qw(0)>zero) then
  !         q(0)=qw(0)*cin(j)
  !      else
  !         q(0)=zero
  !      end if
  !      qyb(0)=zero
  !      do while (t<tfin)
  !         ! get fluxes
  !         do i=1,n
  !            ! get c and csm=dc/dsm (with theta constant)
  !            k=jt(i)
  !            th=thi(i)+(t-ti)*tht(i)
  !            if (isotype(k,j)=="no" .or. sm(i,j)<zero) then ! handle sm<0 here
  !               csm(i)=one/th
  !               c(i,j)=csm(i)*sm(i,j)
  !            else if (isotype(k,j)=="li") then
  !               csm(i)=one/(th+bd(k)*isopar(k,j)%p(1))
  !               c(i,j)=csm(i)*sm(i,j)
  !            else
  !               do it=1,itmax ! get c from sm using Newton's method and bisection
  !                  if (c(i,j)<zero) c(i,j)=zero ! c and sm are >=0
  !                  call isosub(isotype(k,j), c(i,j), isopar(k,j)%p(:), f, fc)
  !                  csm(i)=one/(th+bd(k)*fc)
  !                  dm=sm(i,j)-(bd(k)*f+th*c(i,j))
  !                  dc=dm*csm(i)
  !                  if (sm(i,j)>=zero .and. c(i,j)+dc<zero) then
  !                     c(i,j)=half*c(i,j)
  !                  else
  !                     c(i,j)=c(i,j)+dc
  !                  end if
  !                  if (abs(dm)<eps*(sm(i,j)+10.0_r_2*dsmmax)) exit
  !                  if (it==itmax) then
  !                     write(2,*) "solute: too many iterations getting c"
  !                     stop
  !                  end if
  !               end do
  !            end if
  !         end do
  !         q(1:n-1)=coef1*c(1:n-1,j)+coef2*c(2:n,j)
  !         qya(1:n-1)=coef1*csm(1:n-1)
  !         qyb(1:n-1)=coef2*csm(2:n)
  !         q(n)=qw(n)*c(n,j)
  !         qya(n)=qw(n)*csm(n)
  !         ! get time step
  !         dmax=maxval(abs(q(1:n)-q(0:n-1))/dx)
  !         if (dmax==zero) then
  !            dt=tfin-t
  !         elseif (dmax<zero) then
  !            write(2,*) "solute: errors in fluxes prevent continuation"
  !            stop
  !         else
  !            dt=dsmmax/dmax
  !         end if
  !         if (t+1.1_r_2*dt>tfin) then
  !            dt=tfin-t
  !            t=tfin
  !         else
  !            t=t+dt
  !         end if
  !         sigdt=sig*dt
  !         rsigdt=one/sigdt
  !         ! adjust q for change in theta
  !         q(1:n-1)=q(1:n-1)-sigdt*(qya(1:n-1)*tht(1:n-1)*c(1:n-1,j)+ &
  !              qyb(1:n-1)*tht(2:n)*c(2:n,j))
  !         q(n)=q(n)-sigdt*qya(n)*tht(n)*c(n,j)
  !         ! get and solve eqns
  !         aa(2:n)=qya(1:n-1)
  !         cc(1:n-1)=-qyb(1:n-1)
  !         bb(1:n)=qyb(0:n-1)-qya(1:n)-dx*rsigdt
  !         dd(1:n)=-(q(0:n-1)-q(1:n))*rsig
  !         call dgtsv(n, 1, aa(2:n), bb(1:n), cc(1:n-1), dd(1:n), n, info)
  !         if (info > 0) then
  !            write(2,*) 'solute: singular matrix (01).'
  !            stop
  !         endif
  !         dy(1:n) = dd(1:n)

  !         ! update unknowns
  !         sdrn(j)=sdrn(j)+(q(n)+sig*qya(n)*dy(n))*dt
  !         sm(:,j)=sm(:,j)+dy(1:n)
  !         nssteps(j)=nssteps(j)+1
  !      end do
  !   end do
  ! END SUBROUTINE solute

  !**********************************************************************************************************************

  SUBROUTINE isotope_vap(isotopologue, testcase, litter, n, & ! scalar in
       ns, dx, deltaz, sig, dt, dxL, maxpond, &     ! in soil
       Tsoil0, deltaT, Sliq, deltaSliq, & ! in variables
       Ts, TL, T0, h0new, dh0, SLliq, deltaSLliq, &
       qsig, qlsig, qvsig, qprec, qevap, qrunoff, qex, qd, & ! in fluxes
       var, thetasat, thetar, tortuosity, deltacv, ram, rbh, vlit, cva, civa, & ! in parameter
       thetasatL, deltacvL, &
       cprec, cali, & ! in iso
       ql0, qv0, & ! inout water
       ciso, ciso0, cisoL, cisos, & ! inout iso
       qiso_in, qiso_out, qiso_evap, qiso_trans, qiso_liq_adv, qiso_vap_adv, qiso_liq_diff, qiso_vap_diff) ! out iso

    IMPLICIT NONE

    INTEGER(i_d),                  INTENT(IN)    :: isotopologue ! which isotope
    INTEGER(i_d),                  INTENT(IN)    :: testcase     ! Braud test cases
    LOGICAL,                       INTENT(IN)    :: litter       ! litter or not
    INTEGER(i_d),                  INTENT(IN)    :: n            ! # of soil layers
    INTEGER(i_d),                  INTENT(IN)    :: ns           ! pond or litter (0), nothing (1)
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: dx           ! soil depths
    REAL(r_2),   DIMENSION(1:n-1), INTENT(IN)    :: deltaz       ! soil layer thickness
    REAL(r_2),                     INTENT(IN)    :: sig          ! implicit/explicit time steping constant
    REAL(r_2),                     INTENT(IN)    :: dt           ! time step
    REAL(r_2),                     INTENT(IN)    :: dxL          ! litter layer thickness
    LOGICAL,                       INTENT(IN)    :: maxpond      ! reach maximum pond height
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: Tsoil0       ! soil temperatures soil
    REAL(r_2),   DIMENSION(0:n),   INTENT(IN)    :: deltaT       ! soil temperatures change
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: Sliq         ! soil saturation
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: deltaSliq    ! soil sturation change
    REAL(r_2),                     INTENT(IN)    :: Ts           ! surface temperature
    REAL(r_2),                     INTENT(IN)    :: TL           ! litter/pond temperature
    REAL(r_2),                     INTENT(IN)    :: T0           ! air temperature
    REAL(r_2),                     INTENT(IN)    :: h0new        ! pond height
    REAL(r_2),                     INTENT(IN)    :: dh0          ! pond height change
    REAL(r_2),                     INTENT(IN)    :: SLliq        ! litter saturation
    REAL(r_2),                     INTENT(IN)    :: deltaSLliq   ! litter saturation change
    REAL(r_2),   DIMENSION(0:n),   INTENT(IN)    :: qsig         ! water flux
    REAL(r_2),   DIMENSION(0:n),   INTENT(IN)    :: qlsig        ! liquid water flux
    REAL(r_2),   DIMENSION(0:n),   INTENT(IN)    :: qvsig        ! vapour flux
    REAL(r_2),                     INTENT(IN)    :: qprec        ! precip
    REAL(r_2),                     INTENT(IN)    :: qevap        ! evaporation
    REAL(r_2),                     INTENT(IN)    :: qrunoff      ! runoff
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: qex          ! root extraction
    REAL(r_2),                     INTENT(IN)    :: qd           ! litter to soil drainage
    TYPE(vars),  DIMENSION(1:n),   INTENT(IN)    :: var          ! soil variables
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: thetasat     ! saturation moisture
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: thetar       ! residual moisture
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: tortuosity   ! soil tortuosity
    REAL(r_2),   DIMENSION(1:n),   INTENT(IN)    :: deltacv      !
    REAL(r_2),                     INTENT(IN)    :: ram          ! aerodynamic resistance
    REAL(r_2),                     INTENT(IN)    :: rbh          ! boundary layer resistance
    TYPE(vars),                    INTENT(IN)    :: vlit         ! litter variables
    REAL(r_2),                     INTENT(IN)    :: cva          !
    REAL(r_2),                     INTENT(IN)    :: civa         !
    REAL(r_2),                     INTENT(IN)    :: thetasatL    ! litter saturation moisture
    REAL(r_2),                     INTENT(IN)    :: deltacvL     ! 
    REAL(r_2),                     INTENT(IN)    :: cprec        ! iso conc in precip
    REAL(r_2),                     INTENT(IN)    :: cali         ! iso conc in alimenation water (from below)
    REAL(r_2),                     INTENT(INOUT) :: ql0          ! liquid flux into soil
    REAL(r_2),                     INTENT(INOUT) :: qv0          ! vapour flux into soil
    REAL(r_2),   DIMENSION(0:n),   INTENT(INOUT) :: ciso         ! iso conc in soil
    REAL(r_2),                     INTENT(INOUT) :: ciso0        ! iso conc in air
    REAL(r_2),                     INTENT(INOUT) :: cisoL        ! iso conc in litter/pond
    REAL(r_2),                     INTENT(INOUT) :: cisos        ! iso conc on surface
    REAL(r_2),                     INTENT(OUT)   :: qiso_in      ! iso flux into soil
    REAL(r_2),                     INTENT(OUT)   :: qiso_out     ! iso flux out of soil
    REAL(r_2),                     INTENT(OUT)   :: qiso_evap    ! iso flux of evaporation
    REAL(r_2),                     INTENT(OUT)   :: qiso_trans   ! iso flux of transpiration
    REAL(r_2),   DIMENSION(1:n),   INTENT(OUT)   :: qiso_liq_adv ! liquid iso flux in soil due to advection
    REAL(r_2),   DIMENSION(1:n),   INTENT(OUT)   :: qiso_vap_adv ! vapour iso flux in soil due to advection
    REAL(r_2),   DIMENSION(1:n-1), INTENT(OUT)   :: qiso_liq_diff ! liquid iso flux in soil due to diffusion
    REAL(r_2),   DIMENSION(1:n-1), INTENT(OUT)   :: qiso_vap_diff ! vapour iso flux in soil due to diffusion

    ! Local variables

    REAL(r_2), DIMENSION(0:n)   :: aa, bb, cc, dd, dc,  LHS, RHS
    REAL(r_2), DIMENSION(0:n)   :: alphaplus, dalphaplusdT
    ! diffusivities (liquid, liq-vap, coefft for D, surface H2O vapour, surface minor isotopologue vapour)
    REAL(r_2), DIMENSION(0:n)   :: Dl, Dv
    REAL(r_2)                   :: Dvs, Divs
    REAL(r_2)                   :: patm, nk, alphak, alphak_vdiff, alphak_ldiff
    REAL(r_2)                   :: cevapin, cevapout, qevapin, qevapout, dcevapoutdciso
    ! concentrations of advective fluxes and corresponding partial derivs wrt ciso
    REAL(r_2), DIMENSION(0:n)   :: cql, dcqldca, dcqldcb, cqv, dcqvdca, dcqvdcb
    REAL(r_2), DIMENSION(0:n-1) :: wcql, wcqv
    REAL(r_2), DIMENSION(0:n)   :: beta, deltabeta, betaqv, dbetaqv
    REAL(r_2), DIMENSION(0:n)   :: Dlmean, Dvmean, wl, wv
    REAL(r_2)                   :: coefA, coefB, coefC
    REAL(r_2), DIMENSION(0:n)   :: Seff, deltaSeff, S, Tsoil, cvsig
    REAL(r_2)                   :: h0
    INTEGER(i_d)                :: ns_ciso, j
    REAL(r_2)                   :: num, den, cvL, cv1
    REAL(r_2)                   :: alphaplus_s, alphaplus_0, dalphaplusdT_s
    REAL(r_2)                   :: deltaz0, qevapL, cv0,cvs, qevapoutL, qevapinL
    REAL(r_2)                   :: cevapinL, cevapoutL, dcevapoutdcisoL, dcevapindcisoL
    REAL(r_2), DIMENSION(0:n)   :: Dveff
    REAL(r_2)                   :: w1
    !REAL(r_2), DIMENSION(1:n)   :: ctmp1, ctmp2
    INTEGER(i_d)                :: info
    !CHARACTER(LEN=20)          :: form1
    INTEGER(i_d), PARAMETER        :: formulation=1
    REAL(r_2),    PARAMETER        :: qmin=1.e-15

    aa(0) = zero
    bb(0) = zero
    cc(0) = zero
    dd(0) = zero
    ns_ciso = ns ! index of top layer (0 (pond or litter) or 1 (soil))
    if (maxpond .or. litter) ns_ciso = 0
    if (litter .and. (.not. maxpond) .and. (ns==1)) ciso(0) = cisoL
    deltaz0 = half*dxL + half*dx(1)
    patm    = one

    Tsoil(1:n) = Tsoil0(1:n) + (sig-one)*deltaT(1:n)
    Tsoil(0)   = TL + (sig-one)*deltaT(0)         ! litter or pond temperature

    ! equilibrium fractionation factors at the surface and in the soil
    coefA = zero
    coefB = zero
    coefC = zero
    if (isotopologue==1) then
       coefA = 24844.0_r_2
       coefB = -76.248_r_2
       coefC = 0.052612_r_2
    endif
    if (isotopologue==2) then
       coefA = 1137.0_r_2
       coefB = -0.4156_r_2
       coefC = -0.0020667_r_2
    endif
    alphaplus_s       = one/exp(coefA/((Ts+Tzero)**2)+coefB/(Ts+Tzero)+coefC)       ! at soil or litter surface
    alphaplus_0       = one/exp(coefA/((T0+Tzero)**2)+coefB/(T0+Tzero)+coefC)       ! at soil/litter interface
    alphaplus(0:n)    = one/exp(coefA/((Tsoil(0:n)+Tzero)**2)+coefB/(Tsoil(0:n)+Tzero)+coefC)
    dalphaplusdT_s    = (two*coefA/(Ts+Tzero)**3 + coefB/(Ts+Tzero)**2) &
         / exp(coefA/((Ts+Tzero)**2)+coefB/(Ts+Tzero)+coefC)
    dalphaplusdT(0:n) = (two*coefA/(Tsoil(0:n)+Tzero)**3 + coefB/(Tsoil(0:n)+Tzero)**2) &
         / exp(coefA/((Tsoil(0:n)+Tzero)**2)+coefB/(Tsoil(0:n)+Tzero)+coefC)
    if (testcase==1 .or. testcase==2) then
       alphaplus_s    = one
       alphaplus      = one
       dalphaplusdT   = zero
       dalphaplusdT_s = zero
       alphaplus_0    = one
    endif

    !beta = cv/ cl
    beta(0:n) = alphaplus(0:n)      ! dimensionless
    if (ns==0 .or. maxpond) beta(0) = zero
    ! delta beta
    deltabeta(0:n) =  dalphaplusdT(0:n)*deltaT(0:n)
    if (ns==0 .or. maxpond) deltabeta(0) = zero
    beta = beta + sig*deltabeta         !beta_sig

    ! adjust S and h and cv for sig of time-step
    if (litter) then
       S(0)         = SLliq + deltaSLliq*(sig-one)
       cvsig(0)     =  vlit%cv + deltacvL*(sig-one)
       Seff(0)      = (S(0) + sig*deltaSLliq) + (cvsig(0)*beta(0) + sig*(beta(0)*deltacvL + cvsig(0)*deltabeta(0))) &
            - ( cvsig(0)*S(0)*beta(0) + sig*(S(0)*beta(0)*deltacvL &
            + cvsig(0)*beta(0)*deltaSLliq + cvsig(0)*S(0)*deltabeta(0)))

       deltaSeff(0) = deltaSLliq + beta(0)*deltacvL + cvsig(0)*deltabeta(0) &
            - (S(0)*beta(0)*deltacvL + cvsig(0)*beta(0)*deltaSLliq + cvsig(0)*S(0)*deltabeta(0))
    else
       S(0)         = one
       cvsig(0)     = zero
       Seff(0)      = one
       deltaSeff(0) = one
    endif

    S(1:n)     = Sliq(1:n) + deltaSliq(1:n)*(sig-one)      ! set S to Ssig. N.B. S = Sliq
    h0         = h0new +dh0*(sig-one)    ! set h0 to h0sig
    cvsig(1:n) = var(1:n)%cv + deltacv(1:n)*(sig-one)

    Seff(1:n) = S(1:n) + cvsig(1:n)*beta(1:n)- cvsig(1:n)*S(1:n)*beta(1:n) &
         + thetar(1:n)/thetasat(1:n)

    deltaSeff(1:n) = deltaSliq(1:n) + beta(1:n)*deltacv(1:n) + cvsig(1:n)*deltabeta(1:n) &
         - (S(1:n)*beta(1:n)*deltacv(1:n) + cvsig(1:n)*beta(1:n)*deltaSliq(1:n) + cvsig(1:n)*S(1:n)*deltabeta(1:n))

    ! diffusional fractionation factor in air
    alphak_vdiff =  one
    if (isotopologue==1) alphak_vdiff = one / 1.0251_r_2   ! HDO diffusivity in air (Merlivat 1978)
    if (isotopologue==2) alphak_vdiff = one / 1.0285_r_2   ! H218O diffusivity in air (Merlivat 1978)
    if (testcase >= 1 .and. testcase <= 5) alphak_vdiff = one

    ! kinetic fractionation factor at the surface
    Dvs  = Dva*1.e5_r_2/patm*((Ts+Tzero)/Tzero)**1.88_r_2 ! vapour diffuxivity of water in air (m2s-1)
    Divs =  Dvs
    if (isotopologue==1) Divs = Dvs * alphak_vdiff  ! HDO diffusivity in air
    if (isotopologue==2) Divs = Dvs * alphak_vdiff
    nk = ((thetasat(1)*S(1)-zero)*half + (thetasat(1)*(one-S(1))))/(thetasat(1)-zero)
    if (testcase==7 .or. testcase==8) nk = one
    alphak = one/((Dvs/Divs)**nk)
    !alphak = one/(((Dvs/Divs)**nk + ram/rbh)/(one + ram/rbh)) ! kinetic fractionation factor (< 1)
    if (testcase >= 1 .and. testcase <= 5) alphak = one

    ! liquid diffusivity in the pond and soil
    alphak_ldiff = one
    if (isotopologue==1) alphak_ldiff = one / 1.013_r_2
    if (isotopologue==2) alphak_ldiff = one / 1.026_r_2
    !molecular diffusion of HDO in normal liquid water (m2s-1)
    Dl(1:n) = tortuosity(1:n) * alphak_ldiff * 1.0e-7_r_2*exp(-577.0_r_2/((Tsoil(1:n)+Tzero)-145._r_2))
    !molecular diffusion of HDO in normal liquid water (m2s-1) (pond)
    Dl(0)   = alphak_ldiff * 1.0e-7_r_2*exp(-577.0_r_2/((T0+Tzero)-145._r_2))
    Dl(1:n) = Dl(1:n) * (S(1:n) * (thetasat(1:n)-thetar(1:n)) + thetar(1:n))
    if (testcase >= 1 .and. testcase <= 4) Dl = zero

    ! vapour diffusivity in the soil
    Dv(1:n) = var%Dv * alphak_vdiff  ! isotope diffusivity in soil air spaces
    if (testcase >= 1 .and. testcase <= 5) Dv(1:n) = var%Dv
    Dv(0)   = zero
    Dv(1:n) = Dv(1:n) * cvsig(1:n) !* thetasat(1:n) * (1. - S(1:n))

    do j=1,n-1
       if (abs(cvsig(j)-cvsig(j+1)) > 1.e-8) then
          Dveff(j) = qvsig(j)/(cvsig(j)-cvsig(j+1))*deltaz(j) * alphak_vdiff
          if (testcase >= 1 .and. testcase <= 5) then
             Dveff(j) = qvsig(j)/(cvsig(j)-cvsig(j+1))*deltaz(j)
          endif
       else
          Dveff(j) = var(j)%Dv*alphak_vdiff
          if (testcase >= 1 .and. testcase <= 5) then
             Dveff(j) = var(j)%Dv
          endif
       endif

       if (Dveff(j) <= zero) then
          Dveff(j) = var(j)%Dv*alphak_vdiff
          if (testcase >= 1 .and. testcase <= 5) then
             Dveff(j) = var(j)%Dv
          endif
       endif
    enddo

    if (litter .and. (.not. maxpond) .and. (ns==1)) then
       Dv(0) = vlit%Dv*alphak_vdiff
       if (testcase >= 1 .and. testcase <= 5) Dv(0) = vlit%Dv
       Dv(0) = Dv(0) * cvsig(0)
    endif

    ! upper boundary condition
    cvs = cva + qevap*(ram+rbh) ! concentration of water vapour at soil/air interface (m3 (H2O liq)/ m3 (air))

    ! denominator can be zero, check
    if ((cvs-(var(1)%cv-deltacv(1))) /= zero) then
       Dveff(0) = qv0/(cvs-(var(1)%cv-deltacv(1)))*(half*dx(1))
    else
       ! What should it be in this case?
       Dveff(0) = var(1)%Dv
    endif

    if (var(1)%Dv /= zero) then
       cv1 = -qv0*(half*dx(1))/var(1)%Dv + cvs
    else
       cv1 = var(1)%cv
    endif

    ! do as before
    if (litter) then
       Dveff(0) = -qevap/(cvs-vlit%cv)*(half*dxL)
       if (vlit%Dv /= zero) then
          cvL = qevap*(half*dxL)/vlit%Dv + cvs
       else
          cvL = vlit%cv
       endif
    endif

    if (ql0 > zero) then
       w1 = zero
    else
       w1 = one
    endif

    if (litter .and. (.not. maxpond) .and. (ns==1)) then
       num            = alphak_vdiff*vlit%Dv/(half*dxL)*cvL*alphaplus(0)*ciso(0)  +civa*alphak/(ram+rbh)
       den            = alphak*cvs*alphaplus_s/(ram+rbh) + alphaplus(0)*alphak_vdiff*cvs*vlit%Dv/(half*dxL)
       cisos          = num/den
       dcevapoutdciso = alphak*alphaplus_s*alphak_vdiff*vlit%Dv/(half*dxL)*cvL*alphaplus(0)/den

       qevapL = -(qsig(0) - qd)                     ! vapour flux from top soil layer to litter
       cv0    = qevapL*(half*dxL)/vlit%Dv + (vlit%cv -deltacvL)   ! vapour conc at soil/liter interface
       num    = alphak_vdiff*vlit%Dv/(half*dxL)*(vlit%cv -deltacvL)*alphaplus(0)*ciso(0) &
            +alphak_vdiff*var(1)%Dv/(half*dx(1))*(var(1)%cv-deltacv(1))*alphaplus(1)*ciso(1) &
            - ql0*ciso(1)*w1  + Dl(1)*ciso(1)/(half*dx(1))

       den             = alphak_vdiff*vlit%Dv/(half*dxL)*cv0*alphaplus(0)  + ql0*(one-w1) + &
            alphaplus(1)*alphak_vdiff*cv0*var(1)%Dv/(half*dx(1)) +Dl(1)/(half*dx(1))
       ciso0           = num/den
       dcevapoutdcisoL = alphak*alphaplus_0/den*alphak_vdiff*vlit%Dv/(half*dxL)*vlit%cv*alphaplus(0)*ciso(0)


    elseif ((.not. litter) .and. (.not. maxpond) .and. (ns==1)) then
       Dveff(0) = var(1)%Dv
       num            = alphak_vdiff*Dveff(0)/(half*dx(1))*cv1*alphaplus(1)*ciso(1) &
            - ql0*ciso(1)*w1 +civa*alphak/(ram+rbh) + Dl(1)*ciso(1)/(half*dx(1))
       den            = alphak*cvs*alphaplus_s/(ram+rbh) + ql0*(one-w1) + &
            alphaplus(1)*alphak_vdiff*cvs*Dveff(0)/(half*dx(1)) +Dl(1)/(half*dx(1))
       cisos          = num/den
       dcevapoutdciso = alphak*alphaplus_s
       dcevapoutdciso = dcevapoutdciso*(alphak_vdiff*Dveff(0)/(half*dx(1))*var(1)%cv*alphaplus(1) - &
            ql0*w1  + Dl(1)/(half*dx(1)))/den
    else
       cisos          = ciso(0)
       num            = zero
       den            = zero
       dcevapoutdciso = alphak*alphaplus_s
    endif

    qevapin  = cva/(ram+rbh)
    qevapout = cvs/(ram+rbh)
    cevapin  = civa/cva * alphak
    cevapout = alphak*alphaplus_s * cisos

    if (litter .and. (.not. maxpond) .and. (ns==1)) then
       qevapL          = -(qsig(0) - qd)                        ! vapour flux from top soil layer to litter
       cv0             = qevapL*(half*dxL)/vlit%Dv + vlit%cv    ! vapour conc at soil/liter interface
       qevapoutL       = cv0*vlit%Dv/(half*dxL)
       qevapinL        = vlit%cv*vlit%Dv/(half*dxL)
       cevapinL        = alphak*alphaplus(0)*ciso(0)
       cevapoutL       = alphak*alphaplus_0*ciso(1)
       dcevapoutdcisoL = alphak*alphaplus_0*ciso0
       dcevapindcisoL  = alphak*alphaplus(0)
    endif

    ! concentrations of advective fluxes and corresponding partial derivs wrt ciso
    select case (formulation)
    case (1)
       cql(0:n-1)     = merge(ciso(0:n-1), ciso(1:n), qlsig(0:n-1)>zero)
       dcqldca(0:n-1) = merge(one, zero, qlsig(0:n-1)>zero)
       dcqldcb(0:n-1) = merge(zero,one, qlsig(0:n-1)>zero)
       cql(n)         = ciso(n)
       dcqldca(n)     = one
       dcqldcb(n)     = zero

       cqv(0:n-1)     = merge(ciso(0:n-1), ciso(1:n), qvsig(0:n-1)>zero)
       dcqvdca(0:n-1) = merge(one, zero, qvsig(0:n-1)>zero)
       dcqvdcb(0:n-1) = merge(zero, one, qvsig(0:n-1)>zero)
       cqv(n)         = ciso(n)
       dcqvdca(n)     = one
       dcqvdcb(n)     = zero

       betaqv(0:n-1)  = merge(beta(0:n-1),beta(1:n),qvsig(0:n-1)>0) * alphak_vdiff
       betaqv(n)      = beta(n) * alphak_vdiff

       dbetaqv(1:n-1) = (beta(1:n-1) - beta(2:n))/deltaz(1:n-1)
       dbetaqv(0)     = zero
       dbetaqv(n)     = zero
    case (2)
       cql(0:n-1)     = (ciso(0:n-1)+ciso(1:n))/two
       dcqldca(0:n-1) = half
       dcqldcb(0:n-1) = half
       cql(n)         = ciso(n)
       dcqldca(n)     = one
       dcqldcb(n)     = zero

       cqv(0:n-1)     = (ciso(0:n-1)*beta(0:n-1)+ciso(1:n)*beta(1:n))/two
       dcqvdca(0:n-1) = half * beta(0:n-1)
       dcqvdcb(0:n-1) = half * beta(1:n)
       cqv(n)         = ciso(n)*beta(n)
       dcqvdca(n)     = one*beta(n)
       dcqvdcb(n)     = zero

       betaqv(0:n-1)  = alphak_vdiff
       betaqv(n)      = alphak_vdiff
    case (3)
       wcql(0) = one
       where (qlsig(1:n-1) > qmin)
          wcql(1:n-1) = var(1:n-1)%K*(1.-var(1:n-1)%h)/(var(1:n-1)%K*(1.-var(1:n-1)%h) +var(2:n)%K*var(2:n)%h)
          wcql(1:n-1) = dx(1:n-1)
       elsewhere
          wcql(1:n-1) = half
       endwhere
       cql(0:n-1)     = (wcql(0:n-1)*ciso(0:n-1)+(one-wcql(0:n-1))*ciso(1:n))
       dcqldca(0:n-1) = wcql(0:n-1)
       dcqldcb(0:n-1) = (one-wcql(0:n-1))
       cql(n)         = ciso(n)
       dcqldca(n)     = one
       dcqldcb(n)     = zero

       where (qvsig(1:n-1) > qmin)
          wcqv(1:n-1) = (var(1:n-1)%Kv*var(1:n-1)%h - var(1:n-1)%kE/rhow/rlambda*Tsoil(1:n-1))/ &
               (var(1:n-1)%Kv*var(1:n-1)%h - var(1:n-1)%kE/rhow/rlambda*Tsoil(1:n-1) &
               + var(2:n)%Kv*var(2:n)%h - var(2:n)%kE/rhow/rlambda*Tsoil(2:n))
       elsewhere
          wcqv(1:n-1) = half
       endwhere

       cqv(0)         = zero
       dcqvdca(0)     = zero
       dcqvdcb(0)     = zero
       cqv(1:n-1)     = (wcqv(1:n-1)*ciso(1:n-1)*beta(1:n-1)+(one-wcqv(1:n-1))*ciso(2:n)*beta(2:n))
       dcqvdca(1:n-1) = wcqv(1:n-1)* beta(1:n-1)
       dcqvdcb(1:n-1) = (one-wcqv(1:n-1)) * beta(2:n)
       cqv(n)         = zero
       dcqvdca(n)     = zero
       dcqvdcb(n)     = zero

       betaqv(0:n-1)  = alphak_vdiff
       betaqv(n)      = alphak_vdiff
    case default
       write(2,*) "isotope_vap: illegal formulation [1-3]: ", formulation
       stop
    end select

    ! mean diffusivities
    wl(0)         = h0
    wl(1:n)       =  dx(1:n)
    Dlmean(0:n-1) = (Dl(0:n-1)*wl(0:n-1) + Dl(1:n)*wl(1:n))/(wl(0:n-1)+wl(1:n))
    DLmean(n)     = zero
    where (Dv(1:n) > 1.e-16_r_2)
       wv(1:n) = dx(1:n)*Dv(1:n)
    elsewhere
       wv(1:n) = dx(1:n)
    endwhere
    wv(1:n)       = dx(1:n)
    wv(0)         = dxL
    Dvmean(0)     = zero
    Dvmean(1:n-1) = (Dv(1:n-1)*wv(1:n-1) + Dv(2:n)*wv(2:n))/(wv(1:n-1)+wv(2:n))
    Dvmean(n)     = zero
    if (litter .and. (.not. maxpond) .and. (ns==1)) then
       Dvmean(0) = (Dv(0)*wv(0) + Dv(1)*wv(1))/(wv(0)+wv(1))
    endif

    ! coefficients of tridiagonal matrix
    aa(1) = zero

    aa(2:n) = qlsig(1:n-1)*dcqldca(1:n-1) +qvsig(1:n-1)*betaqv(1:n-1)*dcqvdca(1:n-1) &
         + qvsig(1:n-1)*dbetaqv(1:n-1)*dcqvdca(1:n-1)*Dvmean(1:n-1) &
         + Dlmean(1:n-1)/deltaz(1:n-1) &
         + Dvmean(1:n-1)/deltaz(1:n-1)*beta(1:n-1)

    bb(1)  = -Seff(1)*thetasat(1)*dx(1)/sig/dt &
         - qevapout*dcevapoutdciso &
         - qlsig(1)*dcqldca(1) - qvsig(1)*betaqv(1)*dcqvdca(1) &
         - qvsig(1)*dbetaqv(1)*dcqvdca(1)*Dvmean(1) &
         - Dlmean(1)/deltaz(1) &
         - Dvmean(1)/deltaz(1)*beta(1) &
         - qex(1)

    bb(2:n-1) = -Seff(2:n-1)*thetasat(2:n-1)*dx(2:n-1)/sig/dt &
         + qlsig(1:n-2)*dcqldcb(1:n-2)  +qvsig(1:n-2)*betaqv(1:n-2)*dcqvdcb(1:n-2) &
         + qvsig(1:n-2)*dbetaqv(1:n-2)*dcqvdcb(1:n-2)* Dvmean(1:n-2) &
         - qlsig(2:n-1)*dcqldca(2:n-1) - qvsig(2:n-1)*betaqv(2:n-1)*dcqvdca(2:n-1) &
         - qvsig(2:n-1)*dbetaqv(2:n-1)*dcqvdca(2:n-1)*Dvmean(2:n-1) &
         - (Dlmean(1:n-2)/deltaz(1:n-2) + Dlmean(2:n-1)/deltaz(2:n-1) ) &
         - (Dvmean(1:n-2)/deltaz(1:n-2) + Dvmean(2:n-1)/deltaz(2:n-1))*beta(2:n-1) &
         - qex(2:n-1)

    bb(n) = -Seff(n)*thetasat(n)*dx(n)/sig/dt &
         + qlsig(n-1)*dcqldcb(n-1)  +qvsig(n-1)*betaqv(n-1)*dcqvdcb(n-1) &
         + qvsig(n-1)*dbetaqv(n-1)*dcqvdcb(n-1) *Dvmean(n-1) &
         - qlsig(n)*dcqldca(n) &
         - Dlmean(n-1)/deltaz(n-1) &
         - Dvmean(n-1)/deltaz(n-1)*beta(n) &
         - qex(n)

    cc(1:n-1) = -qlsig(1:n-1)*dcqldcb(1:n-1)  - qvsig(1:n-1)*betaqv(1:n-1)*dcqvdcb(1:n-1) &
         - qvsig(1:n-1)*dbetaqv(1:n-1)*dcqvdcb(1:n-1)*Dvmean(1:n-1) &
         + Dlmean(1:n-1)/deltaz(1:n-1) &
         + Dvmean(1:n-1)/deltaz(1:n-1)*beta(2:n)

    cc(n) = zero

    dd(1) = thetasat(1)*dx(1)/sig/dt*ciso(1)*deltaSeff(1) &
         - qprec*cprec/sig - qevapin*cevapin/sig + qevapout*cevapout/sig &
         + qlsig(1)*cql(1)/sig + qvsig(1)*betaqv(1)*cqv(1)/sig &
         + qvsig(1)*dbetaqv(1)*cqv(1)*Dvmean(1)/sig &
         + Dlmean(1)/deltaz(1)/sig*(ciso(1) - ciso(2)) &
         + Dvmean(1)/deltaz(1)/sig*(ciso(1)*beta(1) - ciso(2)*beta(2) ) &
         + qex(1)*ciso(1)/sig

    dd(2:n-1) = thetasat(2:n-1)*dx(2:n-1)/sig/dt*ciso(2:n-1)*deltaSeff(2:n-1) &
         - qlsig(1:n-2)*cql(1:n-2)/sig &
         - qvsig(1:n-2)*betaqv(1:n-2)*cqv(1:n-2)/sig &
         - qvsig(1:n-2)*dbetaqv(1:n-2)*cqv(1:n-2)*Dvmean(1:n-2)/sig &
         + qlsig(2:n-1)*cql(2:n-1)/sig &
         + qvsig(2:n-1)*betaqv(2:n-1)*cqv(2:n-1)/sig &
         + qvsig(2:n-1)*dbetaqv(2:n-1)*cqv(2:n-1)*Dvmean(2:n-1)/sig &
         - Dlmean(1:n-2)/deltaz(1:n-2)/sig*(ciso(1:n-2) - ciso(2:n-1)) &
         - Dvmean(1:n-2)/deltaz(1:n-2)/sig*(ciso(1:n-2)*beta(1:n-2) - ciso(2:n-1)*beta(2:n-1)) &
         + Dlmean(2:n-1)/deltaz(2:n-1)/sig*(ciso(2:n-1) - ciso(3:n)) &
         + Dvmean(2:n-1)/deltaz(2:n-1)/sig*(ciso(2:n-1)*beta(2:n-1) - ciso(3:n)*beta(3:n)) &
         + qex(2:n-1)*ciso(2:n-1)/sig

    dd(n) = thetasat(n)*dx(n)/sig/dt*ciso(n)*deltaSeff(n) &
         - qlsig(n-1)*cql(n-1)/sig &
         - qvsig(n-1)*betaqv(n-1)*cqv(n-1)/sig &
         - qvsig(n-1)*dbetaqv(n-1)*cqv(n-1)*Dvmean(n-1)/sig &
         + qlsig(n)*cql(n)/sig &
         - Dlmean(n-1)/deltaz(n-1)/sig*(ciso(n-1) - ciso(n)) &
         - Dvmean(n-1)/deltaz(n-1)/sig*(ciso(n-1)*beta(n-1) - ciso(n)*beta(n) ) &
         + qex(n)*ciso(n)/sig

    if (testcase==7 .or. testcase==8) then
       bb(n) = -Seff(n)*thetasat(n)*dx(n)/sig/dt &
            + qlsig(n-1)*dcqldcb(n-1)  +qvsig(n-1)*betaqv(n-1)*dcqvdcb(n-1) &
            + qvsig(n-1)*dbetaqv(n-1)*dcqvdcb(n-1)* Dvmean(n-1) &
            - Dlmean(n-1)/deltaz(n-1) &
            - Dvmean(n-1)/deltaz(n-1)*beta(n) &
            - qex(n)

       dd(n) = thetasat(n)*dx(n)/sig/dt*ciso(n)*deltaSeff(n) &
            - qlsig(n-1)*cql(n-1)/sig &
            - qvsig(n-1)*betaqv(n-1)*cqv(n-1)/sig &
            - qvsig(n-1)*dbetaqv(n-1)*cqv(n-1)*Dvmean(n-1)/sig &
            + qlsig(n)*cali/sig &
            - Dlmean(n-1)/deltaz(n-1)/sig*(ciso(n-1) - ciso(n)) &
            - Dvmean(n-1)/deltaz(n-1)/sig*(ciso(n-1)*beta(n-1) - ciso(n)*beta(n)) &
            + qex(n)*ciso(n)/sig
    endif

    if (maxpond .or. litter) ns_ciso = 0

    if (ns_ciso==0) then
       aa(0) = zero

       bb(0) = -qevapout*dcevapoutdciso -qsig(0) - h0/sig/dt &
            - Dlmean(0)/(half*h0 + half*dx(1))

       cc(0) = Dlmean(0)/(half*h0 + half*dx(1))

       dd(0) = ciso(0)*dh0/dt/sig - cprec*qprec/sig + qevapout*cevapout/sig - qevapin*cevapin/sig + &
            qsig(0)*ciso(0)/sig + Dlmean(0)/(half*h0 + half*dx(1))/sig*(ciso(0) - ciso(1))

       aa(1) = qsig(0) + Dlmean(0)/(h0 + half*dx(1))

       bb(1)  = - Seff(1)*thetasat(1)*dx(1)/sig/dt &
            - qlsig(1)*dcqldca(1) - qvsig(1)*betaqv(1)*dcqvdca(1) &
            - (Dlmean(0)/(half*h0 + half*dx(1)) + Dlmean(1)/deltaz(1)) &
            - ( Dvmean(1)/deltaz(1))*beta(1) &
            - qex(1)

       cc(1) = -qlsig(1)*dcqldcb(1)  - qvsig(1)*betaqv(1)*dcqvdcb(1) &
            + Dlmean(1)/deltaz(1) &
            + Dvmean(1)/deltaz(1)*beta(2)

       dd(1) = thetasat(1)*dx(1)/sig/dt*ciso(1)*deltaSeff(1) &
            - qsig(0)*ciso(0)/sig &
            + qlsig(1)*cql(1)/sig + qvsig(1)*betaqv(1)*cqv(1)/sig &
            - Dlmean(0)/(half*h0 + half*dx(1))/sig*(ciso(0) - ciso(1)*S(1)*thetasat(1)) &
            + Dlmean(1)/deltaz(1)/sig*(ciso(1) - ciso(2)) &
            + Dvmean(1)/deltaz(1)/sig*(ciso(1)*beta(1) - ciso(2)*beta(2)) &
            + qex(1)*ciso(1)/sig

       if (maxpond) then
          bb(0) = bb(0) - qrunoff
          dd(0) = dd(0) + qrunoff*ciso(0)/sig
       endif

       if (litter .and. (.not. maxpond) .and. (ns==1)) then               ! litter and no ponding
          bb(0) = -qevapout*dcevapoutdciso -qevapinL*dcevapindcisoL &
               - Seff(0)*thetasatL*dxL/sig/dt  - ( Dvmean(0)/deltaz0*beta(0) )

          cc(0) = qevapoutL*dcevapoutdcisoL + Dvmean(0)/deltaz0*beta(1)

          dd(0) = thetasatL*dxL/sig/dt*ciso(0)*deltaSeff(0) &
               - cprec*qprec/sig + qevapout*cevapout/sig - qevapin*cevapin/sig &
               - qevapoutL*cevapoutL/sig  + qevapinL*cevapinL/sig + qd*cprec/sig &
               + Dvmean(0)/deltaz0*(ciso(0)*beta(0) - ciso(1)*beta(1))/sig

          aa(1) = qevapinL*dcevapindcisoL + Dvmean(0)/deltaz0*beta(0)

          bb(1) = -Seff(1)*thetasat(1)*dx(1)/sig/dt &
               - qevapoutL*dcevapoutdcisoL &
               - qlsig(1)*dcqldca(1) - qvsig(1)*betaqv(1)*dcqvdca(1) &
               - Dvmean(0)/deltaz0*beta(1) &
               - Dlmean(1)/deltaz(1) &
               - (Dvmean(1)/deltaz(1))*beta(1) &
               - qex(1)

          cc(1) = -qlsig(1)*dcqldcb(1)  -qvsig(1)*betaqv(1)*dcqvdcb(1) &
               + Dlmean(1)/deltaz(1) &
               + Dvmean(1)/deltaz(1)*beta(2)

          dd(1) = thetasat(1)*dx(1)/sig/dt*ciso(1)*deltaSeff(1) &
               - qd*cprec/sig + qevapoutL*cevapoutL/sig -qevapinL*cevapinL/sig &
               + qlsig(1)*cql(1)/sig + qvsig(1)*betaqv(1)*cqv(1)/sig &
               - Dvmean(0)/(deltaz0)/sig*(ciso(0)*beta(0) - ciso(1)*beta(1)) &
               + Dlmean(1)/deltaz(1)/sig*(ciso(1) - ciso(2)) &
               + Dvmean(1)/deltaz(1)/sig*(ciso(1)*beta(1) - ciso(2)*beta(2)) &
               + qex(1)*ciso(1)/sig

       endif                      ! end litter

    endif ! ns_ciso = 0

    ! Solve tri-diagonal matrix
    dc(0) = zero
    !call dgtsv(n-ns_ciso+1, 1, aa(ns_ciso+1:n), bb(ns_ciso:n), cc(ns_ciso:n-1), dd(ns_ciso:n), n-ns_ciso+1, info)
    !if (info > 0) then
     !  write(2,*) 'isotope_vap: singular matrix (01).'
      ! stop
    !endif
    !dc(ns_ciso:n) = dd(ns_ciso:n)
    !call tri(ns_ciso,n,aa,bb,cc,dd,ee,dc)

    ! check for mass balance
    if (1 == 0) then
       LHS (1:n) = thetasat(1:n)*dx(1:n)/dt*(dc(1:n)*Seff(1:n) + ciso(1:n)*deltaSeff(1:n))

       RHS(2:n-1) = qlsig(1:n-2)*(cql(1:n-2) + sig*dc(1:n-2)*dcqldca(1:n-2) + sig*dc(2:n-1)*dcqldcb(1:n-2)) &
            + qvsig(1:n-2)*betaqv(1:n-2)*(cqv(1:n-2) + sig*dc(1:n-2)*dcqvdca(1:n-2) &
            + sig*dc(2:n-1)*dcqvdcb(1:n-2)) &
            -qlsig(2:n-1)*(cql(2:n-1) + sig*dc(2:n-1)*dcqldca(2:n-1) + sig*dc(3:n)*dcqldcb(2:n-1) ) &
            - qvsig(2:n-1)*betaqv(2:n-1)*(cqv(2:n-1) + sig*dc(2:n-1)*dcqvdca(2:n-1) + sig*dc(3:n)*dcqvdcb(2:n-1)) &
            + Dlmean(1:n-2)*((ciso(1:n-2) + sig*dc(1:n-2))-(ciso(2:n-1)+sig*dc(2:n-1)))/deltaz(1:n-2) &
            - Dlmean(2:n-1)*((ciso(2:n-1) + sig*dc(2:n-1))-(ciso(3:n)+sig*dc(3:n)))/deltaz(2:n-1) &
            + Dvmean(1:n-2)*((ciso(1:n-2) + sig*dc(1:n-2))*beta(1:n-2) - (ciso(2:n-1) &
            + sig*dc(2:n-1))*beta(2:n-1))/deltaz(1:n-2) &
            - Dvmean(2:n-1)*((ciso(2:n-1) + sig*dc(2:n-1))*beta(2:n-1) - (ciso(3:n) &
            + sig*dc(3:n))*beta(3:n))/deltaz(2:n-1) &
            - qex(2:n-1)*(ciso(2:n-1) + sig*dc(2:n-1))

       RHS(1) = qprec*cprec - qevapout*(cevapout + sig*dc(1)*dcevapoutdciso) +qevapin*cevapin &
            -qlsig(1)*(cql(1) + sig*dc(1)*dcqldca(1) + sig*dc(2)*dcqldcb(1)) &
            - qvsig(1)*betaqv(1)*(cqv(1) + sig*dc(1)*dcqvdca(1) + sig*dc(2)*dcqvdcb(1)) &
            - Dlmean(1)*((ciso(1) + sig*dc(1))-(ciso(2)+sig*dc(2)))/deltaz(1) &
            - Dvmean(1)*((ciso(1) + sig*dc(1))*beta(1) - (ciso(2) + sig*dc(2))*beta(2))/deltaz(1) &
            - qex(1)*(ciso(1) + sig*dc(1))

       RHS(n) =  qlsig(n-1)*(cql(n-1) + sig*dc(n-1)*dcqldca(n-1) + sig*dc(n)*dcqldcb(n-1)) &
            + qvsig(n-1)*betaqv(n-1)*(cqv(n-1) + sig*dc(n-1)*dcqvdca(n-1) + sig*dc(n)*dcqvdcb(n-1)) &
            -qsig(n)*(ciso(n) + sig*dc(n)) &
            + Dlmean(n-1)*((ciso(n-1) + sig*dc(n-1))-(ciso(n)+sig*dc(n)))/deltaz(n-1) &
            + Dvmean(n-1)*((ciso(n-1) + sig*dc(n-1))*beta(n-1) - (ciso(n) + sig*dc(n))*beta(n))/deltaz(n-1) &
            -qex(n)*(ciso(n)+sig*dc(n))

       if (ns_ciso==0) then
          LHS(0) = (ciso(0)*dh0 + h0*dc(0))/dt

          RHS(0) = qprec*cprec - qevapout*(cevapout + sig*dc(0)*dcevapoutdciso) &
               + qevapin*cevapin - qsig(0)*(ciso(0) + sig*dc(0)) &
               - Dlmean(0)*((ciso(0) + sig*dc(0))-(ciso(1)+sig*dc(1)))/(half*h0 + half*dx(1))

          RHS(1) = qsig(0)*(ciso(0)+sig*dc(0)) &
               -qlsig(1)*(cql(1) + sig*dc(1)*dcqldca(1) + sig*dc(2)*dcqldcb(1)) &
               - qvsig(1)*betaqv(1)*(cqv(1) + sig*dc(1)*dcqvdca(1) + sig*dc(2)*dcqvdcb(1)) &
               + Dlmean(0)*((ciso(0) + sig*dc(0))-(ciso(1)+sig*dc(1)))/(half*h0 + half*dx(1)) &
               - Dlmean(1)*((ciso(1) + sig*dc(1))-(ciso(2)+sig*dc(2)))/deltaz(1) &
               - Dvmean(1)*((ciso(1) + sig*dc(1))*beta(1) - (ciso(2) + sig*dc(2))*beta(2))/deltaz(1) &
               - qex(1)*(ciso(1) + sig*dc(1))

          if (maxpond) then
             RHS(0) = RHS(0) - qrunoff*(ciso(0)+sig*dc(0))
          endif

          if (litter .and. (.not. maxpond) .and. (ns==1)) then               ! litter and no ponding
             RHS(0) = thetasatL*dxL/dt*(ciso(0)*deltaSeff(0)+ dc(0)*Seff(0) )

             LHS(0) = qprec*cprec-qd*cprec - qevapout*(cevapout+sig*dcevapoutdciso*dc(0)) &
                  + qevapin*cevapin &
                  + qevapoutL*(cevapoutL+sig*dcevapoutdcisoL*dc(1)) &
                  -qevapinL*(cevapinL + sig*dcevapindcisoL *dc(0)) &
                  -Dvmean(0)/(deltaz0)*(beta(0)*(ciso(0)+dc(0)) - beta(1)*(ciso(1)+sig*dc(1)) )

             RHS(1) = qd*cprec - qevapoutL*(cevapoutL+ sig*dcevapoutdcisoL*dc(1)) &
                  +qevapinL*(cevapinL + sig*dcevapindcisoL*dc(0)) &
                  -qlsig(1)*(cql(1) + sig*dc(1)*dcqldca(1) + sig*dc(2)*dcqldcb(1)) &
                  - qvsig(1)*betaqv(1)*(cqv(1) + sig*dc(1)*dcqvdca(1) + sig*dc(2)*dcqvdcb(1)) &
                  + Dvmean(0)*(beta(0)*(ciso(0) + sig*dc(0))-beta(1)*(ciso(1)+sig*dc(1)))/(deltaz0) &
                  - Dlmean(1)*((ciso(1) + sig*dc(1))-(ciso(2)+sig*dc(2)))/deltaz(1) &
                  - Dvmean(1)*((ciso(1) + sig*dc(1))*beta(1) - (ciso(2) + sig*dc(2))*beta(2))/deltaz(1) &
                  - qex(1)*(ciso(1) + sig*dc(1))
          endif                      ! end litter and no ponding
       end if ! ns_ciso==0

       if (testcase==7 .or. testcase==8) then
          RHS(n) =   qlsig(n-1)*(cql(n-1) + sig*dc(n-1)*dcqldca(n-1) + sig*dc(n)*dcqldcb(n-1)) &
               + qvsig(n-1)*betaqv(n-1)*(cqv(n-1) + sig*dc(n-1)*dcqvdca(n-1) +  sig*dc(n)*dcqvdcb(n-1)) &
               -qsig(n)*cali &
               + Dlmean(n-1)*((ciso(n-1) + sig*dc(n-1))-(ciso(n)+sig*dc(n)))/deltaz(n-1) &
               + Dvmean(n-1)*((ciso(n-1) + sig*dc(n-1))*beta(n-1) - (ciso(n) + sig*dc(n))*beta(n))/deltaz(n-1) &
               -qex(n)*(ciso(n)+sig*dc(n))
       endif
    endif ! 1==0

    ! isotopic fluxes
    qiso_in    = qprec*cprec +qevapin*cevapin
    qiso_out   = qevapout*(cevapout + sig*dc(1)*dcevapoutdciso)
    qiso_evap  = qevapout*(cevapout + sig*dc(1)*dcevapoutdciso) - qevapin*cevapin
    qiso_trans = sum(qex(1:n)*(ciso(1:n)+sig*dc(1:n)),1)

    qiso_liq_diff(1:n-1) = Dlmean(1:n-1)*((ciso(1:n-1) + sig*dc(1:n-1))-(ciso(2:n)+sig*dc(2:n)))/deltaz(1:n-1)
    qiso_vap_diff(1:n-1) = Dvmean(1:n-1)*((ciso(1:n-1) + sig*dc(1:n-1))*beta(1:n-1) &
         - (ciso(2:n) + sig*dc(2:n))*beta(2:n))/deltaz(1 :n-1)
    qiso_liq_adv(1:n-1)  = qlsig(1:n-1)*(cql(1:n-1) + sig*dc(1:n-1)*dcqldca(1:n-1) + sig*dc(2:n)*dcqldcb(1:n-1))
    qiso_liq_adv(n)      = qlsig(n)*(cql(n) + sig*dc(n))
    qiso_vap_adv(1:n-1)  = qvsig(1:n-1) * betaqv(1:n-1) &
         * (cqv(1:n-1) + sig*dc(1:n-1)*dcqvdca(1:n-1) + sig*dc(2:n)*dcqvdcb(1:n-1))
    qiso_vap_adv(n)      = zero

    ciso = ciso + dc
    if (litter .and. (.not. maxpond) .and. (ns==1)) cisoL = ciso(0)

    ! if (testcase==8) then
    !      ctmp1(1:n) = var(1:n)%cvsat
    !      ctmp2(1:n) = var(1:n)%rh
    !      Dveff(n) = zero ! otherwise undefined
    !      write(form1,"(A,I3,A)") "(", 6*(n+1)+6, "f17.7)"
    !      write(31,form1) Dl(1), Dl(1:n), Dveff(0)/alphak_vdiff, Dveff(1:n)/alphak_vdiff, &
    !           ctmp1(1:n), ctmp2(1:n), alphaplus_s, alphaplus(1:n), Ts, Tsoil(1:n), alphak, qevap
    ! endif

  END SUBROUTINE isotope_vap

  !**********************************************************************************************************************

END MODULE cable_sli_solve
