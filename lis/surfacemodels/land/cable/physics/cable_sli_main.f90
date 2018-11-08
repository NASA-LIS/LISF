SUBROUTINE cable_sli_main(irec, dt, veg, soil, ssoil, met, canopy, air) 

  ! Main subroutine for Soil-litter-iso soil model
  ! Vanessa Haverd, CSIRO Marine and Atmospheric Research, 2010
  USE cable_dimensions,  ONLY: r_1, r_2, i_d, mp_patch, ms ! mp_patch=#landpoints*patches, ms=# soil layers
  USE cable_types,       ONLY: veg_parameter_type, soil_parameter_type, soil_snow_type, met_type, &
                                canopy_type, air_type
  USE cable_physical_constants, ONLY: tfrz
  USE cable_sli_numbers,        ONLY: zero, half, one, four, thousand, & ! numbers
                                Tzero, rhow, cswat, &                           ! variables
                                vars_met, vars , params                    ! types
  USE cable_sli_utils,          ONLY: x, dx, par, setpar, setx, plit, dxL, setlitterpar, esat
  USE cable_sli_roots,          ONLY: setroots
  USE cable_sli_solve,          ONLY: solve

  IMPLICIT NONE

  INTEGER(i_d),              INTENT(IN)    :: irec
  REAL,                      INTENT(IN)    :: dt
  TYPE(veg_parameter_type),  INTENT(INOUT) :: veg     ! all r_1
  TYPE(soil_parameter_type), INTENT(INOUT) :: soil    ! all r_1
  TYPE(soil_snow_type),      INTENT(INOUT) :: ssoil   ! r_1, r_2 desaster
  TYPE(met_type),            INTENT(INOUT) :: met     ! all r_1
  TYPE(canopy_type),         INTENT(INOUT) :: canopy  ! all r_1
  TYPE(air_type),            INTENT(INOUT) :: air     ! all r_1

  INTEGER(i_d), PARAMETER :: ns=0 ! soil horizons, solute
  REAL(r_2),    PARAMETER :: emsoil=0.97
  REAL(r_2),    PARAMETER :: rhocp=1.1822e3
  REAL(r_2),    PARAMETER :: Dva = 2.17e-5
  INTEGER(i_d) :: j, setroot
  REAL(r_2)    :: ti, tf
  TYPE(vars_met), DIMENSION(1:mp_patch)      :: vmet ! Meteorology above soil
  TYPE(vars),     DIMENSION(1:mp_patch)      :: vlit
  INTEGER(i_d),   DIMENSION(1:mp_patch)      :: nsteps
  REAL(r_2),      DIMENSION(1:mp_patch,1:ms) :: Tsoil, S, thetai, Jsensible
  REAL(r_2),      DIMENSION(1:mp_patch)      :: SL, TL, T0
  REAL(r_2),      DIMENSION(1:mp_patch)      :: drn, evap, infil, qprec, runoff
  REAL(r_2),      DIMENSION(1:mp_patch)      :: win, wp, wpi,hp, hpi, deltah0, h0old, discharge
  REAL(r_2),      DIMENSION(1:mp_patch)      :: ip, ipi ! volumetric ice content of profile (final and initial)
  REAL(r_2),      DIMENSION(1:mp_patch,1:ms) :: FS, wex, csoil
  REAL(r_2),      DIMENSION(1:mp_patch)      :: rh0, rhsurface
  ! surface temperature (top of top soil layer or top of litter layer)
  REAL(r_2),      DIMENSION(1:mp_patch)      :: Tsurface
  REAL(r_2),      DIMENSION(1:mp_patch)      :: gr, grc
  REAL(r_2),      DIMENSION(1:mp_patch)      :: Etrans, gamma
  REAL(r_2),      DIMENSION(1:mp_patch)      :: G0, H, lE
  REAL(r_2),      DIMENSION(1:mp_patch,0:ms) :: qh, qvsig, qlsig, qvTsig, qvh
  REAL(r_2),      DIMENSION(1:mp_patch)      :: esatvar
  REAL(r_2),      DIMENSION(1:mp_patch)      :: deltaTa, lE_old, TL_test, SA, SB, wpAi, wpBi, wpA, wpB
  REAL(r_2),      DIMENSION(1:mp_patch)      :: evap_pot, deltaice_cum_T, deltaice_cum_S, zdelta, fws
  REAL(r_2),      DIMENSION(1:mp_patch)      :: Qadvcum,Jcol_sensible,Jcol_latent_S,Jcol_latent_T
  REAL(r_2),      DIMENSION(1:mp_patch)      :: tmp1d1, tmp1d2
    ! Model switches
  INTEGER(i_d), PARAMETER :: litter       = 2 ! which litter model
  ! 0: no litter
  ! 1: full litter
  ! 2: litter resistance
    ! Model switches
  INTEGER(i_d), PARAMETER :: advection      = 0 ! heat advection by water
  INTEGER(i_d), PARAMETER :: isotopologue = 0 ! which isotope
  ! 0: no isotope calculations
  ! 1: HDO
  ! 2: H218O
  ! 3: HDO & H218O
  INTEGER(i_d), PARAMETER :: testcase     = 0 ! isotopic test cases 1-8
  ! 0: normal run
  INTEGER(i_d), PARAMETER :: septs        = 0 ! coupled or uncoupled energy and water calculation
  ! 0: coupled calc
  ! 1: uncoupled energy (T) and moisture (S)
  INTEGER(i_d), PARAMETER :: condition    = 3 ! condition matrix before solving
  ! 0: no conditioning
  ! 1: condition columns
  ! 2: condition lines
  ! 3: condition first lines then columns


  ! output files for testing purposes
  if (irec .eq. 1) then
     open (unit=332,file="vh08.out",status="replace",position="rewind")
     open (unit=334,file="S.out",status="replace",position="rewind")
     open (unit=336,file="Tsoil.out",status="replace",position="rewind")
     open (unit=335,file="SEB.out",status="replace",position="rewind")
     !open(unit=37, file = "c:\soil_model\cable_met_test.inp",status="replace",position="rewind")
     !open (unit=337,file="soil_log.out",status="replace",position="rewind")
     open(unit=338, file="thetai.out", status="replace", position="rewind")
     open(unit=339, file="test_dt.out",status="replace", position="rewind")
  endif

  ! Save soil / snow surface temperature from last time step:
  ssoil%tssold(:) = ssoil%tss(:)

  ! set layer thicknesses
  if (.not. allocated(x)) call setx(mp_patch, ms, soil)

  ! Set root density distribution (leave in for sli offline)
  setroot = 0  ! reset rooting depths
  if (setroot == 1) then
     call setroots(x(:,:)*100.0_r_2, real(veg%F10(:),r_2), real(veg%ZR(:),r_2)*100.0_r_2, FS(:,:))
  else
     FS(:,:) = real(veg%froot(:,:),r_2)
  endif

  ! set required soil hydraulic params
  if (.not. allocated(par)) then
     call setpar(mp_patch, ms, x(:,:)-half*dx(:,:), soil)
     soil%swilt_vec(:,:) = merge(spread(soil%swiltB(1:mp_patch),2,ms),soil%swilt_vec(:,:),par(:,:)%isbhorizon==one)
  endif

  ! If we want solutes:
!!$     if (.not. allocated(bd))      allocate(bd(soil%nhorizons(k)))
  
  ! Litter parameters:
  if (.not. allocated(plit)) call setlitterpar(mp_patch, soil)
  
  ! Met data above soil:
  vmet(:)%Ta    = real(met%Tvair(:)-tfrz,r_2)
  vmet(:)%Da    = real(met%dva(:),r_2)
  vmet(:)%Rn    = real(canopy%fns(:),r_2)
  vmet(:)%rbh   = real(ssoil%rtsoil(:),r_2)
  vmet(:)%rbw   = real(ssoil%rtsoil(:),r_2)
  esatvar(:)    = esat(vmet(:)%Ta)
  vmet(:)%rha   = max(min((esat(vmet(:)%Ta)-vmet(:)%Da)/esat(vmet(:)%Ta),one),0.1_r_2)
  vmet(:)%cva   = vmet(:)%rha * esat(vmet(:)%Ta)*0.018_r_2/thousand/8.314_r_2/(vmet(:)%Ta+Tzero) ! m3 H2O (liq) m-3 (air)
  vmet(:)%phiva = Dva * vmet(:)%cva
  Etrans(:)     = real(canopy%fevc(:)/air%rlam(:),r_2)/thousand ! m s-1
  
  ! Initialisations:
  if (irec == 1) then
     S(:,:)         = min((real(ssoil%wb(:,:),r_2)-par(:,:)%thr) / par(:,:)%thre, 1.0_r_2) ! Assume theta_r is zero
     rh0(:)         = vmet(:)%rha
     rhsurface(:)   = rh0(:) ! initialise rel. humidity at surface
     Tsurface(:)    = vmet(:)%Ta
     T0(:)          = Tsurface(:)
     deltaTa(:)     = zero
     lE_old(:)      = zero
     ssoil%h0(:)    = zero
     ssoil%thetai(:,:) = zero
     met%tk_old(:)     = met%tk(:)
     ! do not seem to be initialised elsewhere
     ssoil%snowd(:) = 0.0_r_1
     zdelta(:)      = x(:,ms)
	 ssoil%gammzz = zero
	 !ssoil%tgg(:,6)=5.0+Tfrz
  else
     rhsurface(:)   = ssoil%rhsurface(:)
     Tsurface(:)    = ssoil%Tsurface(:)
     T0(:)          = ssoil%Tsurface(:)
     rh0(:)         = ssoil%rh0(:)
     deltaTa(:)     = zero
     lE_old(:)      = ssoil%lE(:)
     S(:,:)         = ssoil%S(:,:)                ! degree of soil saturation
     zdelta(:)      = ssoil%zdelta(:)
  endif
  
  SL(:)       = ssoil%SL(:)   ! degree of litter saturation
  TL(:)       = ssoil%TL(:)   ! litter T
  Tsoil(:,:)  = real(ssoil%tgg(:,:)-tfrz,r_2)
  ssoil%smelt = zero
  
  !MC! one or the other
  deltaTa(:)  = real(met%tk(:)-met%tk_old(:),r_2)
  deltaTa(:)  = zero
  !MC! one or the other
  lE_old(:)   = ssoil%lE(:)
  gamma(:)    = real(veg%gamma(:),r_2)
  
  gr(:)       = four * emsoil * (vmet(:)%Ta+tfrz)**3 *5.67e-8_r_2 ! radiation conductance Wm-2K-1
  grc(:)      = one/vmet(:)%rbh + gr(:)/rhocp
  vmet(:)%rrc = one/grc(:)                            ! resistance to radiative and convective heat transfer
  qprec(:)    = real(canopy%through(:),r_2)/thousand/dt              ! precip rate (m s-1)
  h0old(:)    = ssoil%h0(:) ! pond height 
  
  ! Heat balance variables
 
  
  ! Water balance variables:
  ipi(:)   = sum(ssoil%thetai(:,:)*dx(:,:),2)     ! ice in profile initially
  ! water in profile initially
  wpi(:)   = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*dx(:,:),2) + plit(:)%thre*SL(:)*dxL(:)
  wpAi(:)  = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*par(:,:)%isahorizon*dx(:,:),2) + &
       plit(:)%thre*SL(:)*dxL(:) ! water in profile initially
  wpBi(:)  = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*par(:,:)%isbhorizon*dx(:,:),2)

  nsteps    = 0
  ti        = zero
  tf        = dt ! initial and final times
  win(:)    = zero ! water input (total precip)
  evap(:)   = zero
  runoff(:) = zero
  infil(:)  = zero
  drn(:)    = zero
  
  ! rewind(337)
  ! write(337,*) irec,  ' call_soil_litter'
  
  if (irec.eq.429) then ! for debugging
    ! write (*,*)
   endif

  call solve(ti, tf, irec, mp_patch, qprec(:), ms, ns, dx(:,:), &
       ssoil%h0(:), S(:,:), thetai(:,:),Jsensible(:,:), Tsoil(:,:), evap(:), &
       evap_pot(:), runoff(:), infil(:), drn(:), discharge(:), qh(:,:), &
       nsteps, vmet(:), vlit(:),csoil(:,:),T0(:), rh0(:), Tsurface(:), &
       rhsurface(:), H(:), lE(:), G0(:),Qadvcum(:),Jcol_sensible(:), &
	   Jcol_latent_S(:),Jcol_latent_T(:), deltaice_cum_T(:), &
       deltaice_cum_S(:), dxL(:), zdelta(:), SL(:), TL(:), &
       plit(:), par(:,:), Etrans(:), fws(:), gamma=gamma(:), &
       wex=wex(:,:), FS=FS(:,:), qvsig=qvsig(:,:), qlsig=qlsig(:,:), qvTsig=qvTsig(:,:), qvh=qvh(:,:), &
       deltaTa=deltaTa(:), lE_old=lE_old(:), TL_test=TL_test(:), &
       dolitter=litter, doisotopologue=isotopologue, dotestcase=testcase, dosepts=septs, docondition=condition, &
	   doadvection =advection)
  
  H(:)      = H(:)/(tf-ti)
  lE(:)     = lE(:)/(tf-ti)
  G0(:)     = G0(:)/(tf-ti)
  Jcol_latent_S(:) = Jcol_latent_S(:)/(tf-ti)
  Jcol_latent_T(:) = Jcol_latent_T(:)/(tf-ti)
  Jcol_sensible(:) = Jcol_sensible(:)/(tf-ti)
  Qadvcum(:) = Qadvcum(:)/(tf-ti)

  
  deltah0(:) = ssoil%h0(:)-h0old(:)
  
  ip(:)  = sum(thetai(:,:)*dx(:,:),2)  ! ice in profile at tf
  ! water at tf
  wp(:)  = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*dx(:,:),2) + plit(:)%thre*SL(:)*dxL(:)
  win(:) = win(:) + qprec(:)*(tf-ti)
  wpA(:) = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*par(:,:)%isahorizon*dx(:,:),2) + &
       plit(:)%thre*SL(:)*dxL(:) ! water in profile initially
  wpB(:) = sum((par(:,:)%thr + (par(:,:)%the-par(:,:)%thr)*S(:,:))*par(:,:)%isbhorizon*dx(:,:),2)
  SA(:)  = (wpA(:)/sum(dx(:,:)*par(:,:)%isahorizon,2) - par(:,1)%thr)/(par(:,1)%the - par(:,1)%thr)
  SB(:)  = (wpB(:)/sum(dx(:,:)*par(:,:)%isbhorizon,2) - par(:,ms)%thr)/(par(:,ms)%the - par(:,ms)%thr)
  where (all(par(:,:)%isbhorizon==zero,2)) SB(:) = zero
       
  !if (mp_patch == 1) then
  if (1 == 1) then
     write(332,"(i8,i8,14e15.6)") irec,nsteps(1),wp(1)-wpi(1),infil(1)-drn(1),runoff(1),&
          win(1)-(wp(1)-wpi(1)+deltah0(1)+runoff(1)+evap(1)+drn(1))-Etrans(1)*dt,wp(1),evap(1),evap_pot(1),infil(1),&
          drn(1),ssoil%h0(1),Etrans(1)*dt,discharge(1), fws(1), (ip(1)-ipi(1))*920./1000.
     write(334,"(100f15.6)") S(1,:)
     write(336,"(100f15.6)") Tsoil(1,:)
     write(335,"(100f20.6)") vmet(1)%Ta, T0(1), rh0(1), H(1), lE(1), &
	                         G0(1),Jcol_sensible(1),Jcol_latent_S(1), Jcol_latent_T(1), &
	                         vmet(1)%Rn, TL(1), SL(1), deltaice_cum_T(1), &
                              deltaice_cum_S(1), rhsurface(1), Tsurface(1), vmet(1)%rha, &
							   csoil(1,:)*dx(1,:), Qadvcum(1), sum((Jsensible(1,:)-ssoil%gammzz(1,:)),1)
     write(338,"(100f18.6)") thetai(1,:)
  endif

  tmp1d1 = sum((Jsensible(1,:)-ssoil%gammzz(1,:)),1)

  
  ! Update variables for output:
  ssoil%tss(:)       = real(Tsurface(:),r_1) + tfrz
  ssoil%tgg(:,:)     = real(Tsoil(:,:),r_1) + tfrz
  ssoil%wb(:,:)      = real(S(:,:)*(par(:,:)%thr+(par(:,:)%the-par(:,:)%thr)),r_1)
  ssoil%wbtot(:)     = real(wp(:)*thousand,r_1)
  canopy%ga(:)       = real(G0(:),r_1)
  canopy%fhs(:)      = real(H(:),r_1)
  canopy%fes(:)      = real(lE(:),r_1)
  ssoil%hflux(:,:)   = qh(:,:)
  ssoil%rnof1(:)     = real(runoff(:)*thousand,r_1)
  ssoil%rnof2(:)     = real(discharge(:)*thousand,r_1)
  ssoil%zdelta(:)    = zdelta(:)
  ssoil%S(:,:)       = S(:,:)
  ssoil%SL(:)        = SL(:)
  ssoil%TL(:)        = TL(:)
  ssoil%delwcol(:)   = (wp(:)-wpi(:)+deltah0(:))*thousand
  ssoil%Tsurface(:)  = Tsurface(:)
  ssoil%rh0(:)       = rh0(:)
  ssoil%rhsurface(:) = rhsurface(:)
  ssoil%lE(:)        = lE(:)
  ssoil%evap(:)      = evap(:)*thousand
  ssoil%SA(:)        = SA(:)
  ssoil%SB(:)        = SB(:)
  ssoil%delwcolA(:)  = (wpA(:)-wpAi(:))*thousand
  ssoil%delwcolB(:)  = (wpB(:)-wpBi(:))*thousand
  ssoil%rex(:,:)     = wex(:,:)*thousand
  ssoil%rlitt(:)     = dxL(:)/vlit(:)%Dv
  ssoil%rexA(:)      = sum(wex(:,:)*par(:,:)%isahorizon,2)*thousand
  ssoil%rexB(:)      = sum(wex(:,:)*par(:,:)%isbhorizon,2)*thousand
  ssoil%thetai(:,:)  = thetai(:,:)
  ssoil%gammzz(:,:) = Jsensible(:,:)
  do j=1, ms
     where(par(:,j)%isahorizon==one) ssoil%leachAB(:) = (qvsig(:,j)+qlsig(:,j))*thousand
  enddo

  ! Update total latent heat to reflect updated soil component:
  canopy%fe(:) = canopy%fev(:) + canopy%fes(:)
  ! Update total sensible heat to reflect updated soil component:
  canopy%fh(:) = canopy%fhv(:) + canopy%fhs(:)

  ! store air temprature
  met%tk_old(:)  = met%tk(:)

END SUBROUTINE cable_sli_main
