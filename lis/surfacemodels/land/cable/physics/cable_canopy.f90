! cable_canopy.f90
!
! Source file containing canopy code for CABLE
!
! Development by Ying-Ping Wang, Eva Kowalczyk, Ray Leuning
! Gab Abramowitz, Martin Dix, Harvey Davies, Mike Raupach, 
!
! bugs to bernard.pak@csiro.au
!
! This file contains modules:
!   canopy_module
! The subroutines included are:
!   define_canopy
! The functions included are:
!   qsatf,
!   ej3x,
!   ej4x,
!   xvcmxt4,
!   xvcmxt3,
!   xejmxt3,
!   psim,
!   psis,
!   rplant, and
!   rsoil
!
! Most user-defined types (e.g. met%tk) are defined in cable_types module
! in cable_variables.f90

MODULE cable_canopy
  USE cable_dimensions, ONLY: r_1,r_2,mp_patch,mf,ms,i_d
  USE cable_types
  USE cable_photosynthetic_constants
  USE cable_radiation, ONLY: radiation, sinbet
  USE cable_air, ONLY: define_air
  USE cable_physical_constants
  IMPLICIT NONE
  PRIVATE
  PUBLIC define_canopy
CONTAINS
  
  !--------------------------------------------------------------------------
  SUBROUTINE define_canopy(ktau,bal,rad,rough,air,met,dels,ssoil,soil,&
       veg,bgc,canopy,model_structure)
    TYPE (balances_type),INTENT(INOUT)   :: bal
    TYPE (radiation_type), INTENT(INOUT) :: rad
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (air_type), INTENT(INOUT)       :: air
    TYPE (met_type), INTENT(INOUT)       :: met
    REAL(r_1), INTENT(IN)                :: dels ! integration time setp (s)
    TYPE (soil_snow_type), INTENT(INOUT) :: ssoil
    TYPE (bgc_pool_type),INTENT(IN)      :: bgc
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE (canopy_type), INTENT(INOUT)    :: canopy
    TYPE (model_structure_type), INTENT(IN)  :: model_structure
    INTEGER(i_d), INTENT(IN)             :: ktau ! integration step number
    REAL(r_1), DIMENSION(mp_patch,mf):: a1c3 ! Spatially varying a1c3
    REAL(r_1) :: a1c3_380 = 4.34 ! fitted value from Hawkwsbury data
    REAL(r_1) :: a1c3_620 = 5.23 ! fitted value from Hawkwsbury data
    REAL(r_1), DIMENSION(mp_patch,mf)   :: abs_deltlf ! ABS(deltlf)
    REAL(r_2), DIMENSION(mp_patch,mf,3) :: ancj ! soln to quad eqn
    REAL(r_2), DIMENSION(mp_patch,mf) :: anx ! net photos. prev iteration
    REAL(r_2), DIMENSION(mp_patch,mf) :: an_y ! net photosynthesis soln
  !  REAL(r_1), DIMENSION(mp_patch) :: avgtrs !root weighted mean soil temperature
  !  REAL(r_1), DIMENSION(mp_patch) :: avgwrs !root weighted mean soil moisture
    REAL(r_1), DIMENSION(mp_patch,mf) :: ca2  ! 2D CO2 concentration
    REAL(r_1), DIMENSION(mp_patch) :: cansat ! max canopy intercept. (mm)
    REAL(r_2), DIMENSION(mp_patch,mf,3) :: ci ! intercellular CO2 conc.
    REAL(r_1), PARAMETER  :: co2cp3=0.0 ! CO2 compensation pt C3
    REAL(r_2), DIMENSION(mp_patch,mf,3) :: coef0 ! CO2 comp. pt coeff 1
    REAL(r_2), DIMENSION(mp_patch,mf,3) :: coef1 ! " 2
    REAL(r_2), DIMENSION(mp_patch,mf,3) :: coef2 ! " 3
    REAL(r_1), DIMENSION(mp_patch,mf) :: conkct ! Michaelis Menton const.
    REAL(r_1), DIMENSION(mp_patch,mf) :: conkot ! Michaelis Menton const.
    REAL(r_1), DIMENSION(mp_patch):: cscale ! scaling between 2 hawkesbury elev co2 vals
    REAL(r_2), DIMENSION(mp_patch,mf) :: csx ! leaf surface CO2 concentration
    REAL(r_2), DIMENSION(mp_patch,mf,3) :: cx  ! "d_{3}" in Wang and Leuning, 1998, appendix E
    REAL(r_1), DIMENSION(mp_patch,mf) :: d0c3 ! Empirical coef for vpd sensitivity of stomata
    REAL(r_1), DIMENSION(mp_patch,mf) :: da2 ! 2D sat vap pres deficit
    REAL(r_1), DIMENSION(mp_patch,mf) :: dva2 ! 2D in canopy sat vap pres deficit
    REAL(r_2), DIMENSION(mp_patch,mf,3) :: delcx ! discriminant  in quadratic in eq. E7 Wang and Leuning, 1998
    REAL(r_1), DIMENSION(mp_patch,mf) :: deltlf ! deltlfy of prev iter.
    REAL(r_1), DIMENSION(mp_patch,mf) :: deltlfy ! del temp successive iteration
    REAL(r_1), DIMENSION(mp_patch) :: dq ! sat spec hum diff.
    REAL(r_1), DIMENSION(mp_patch,mf) :: dsatdk2 ! 2D dsatdk
    REAL(r_1), DIMENSION(mp_patch,mf) :: dsx ! leaf surface vpd
    REAL(r_2), DIMENSION(mp_patch,mf) :: ecx ! lat. hflux big leaf
    REAL(r_1), DIMENSION(mp_patch,mf) :: ejmax2 ! jmax of big leaf
    REAL(r_1), DIMENSION(mp_patch,mf) :: ejmxt3 ! jmax big leaf C3 plants
    REAL(r_2), DIMENSION(mp_patch,mf) :: ecy ! lat heat fl dry big leaf
    REAL(r_1), DIMENSION(mp_patch,mf) :: frac42 ! 2D frac4
    REAL(r_2), DIMENSION(mp_patch) :: fwsoil ! soil water modifier of stom. cond.
    REAL(r_2), DIMENSION(mp_patch) :: gaw ! aerodynamic conduct. for water
    REAL(r_2), DIMENSION(mp_patch,mf) :: gaw2 ! 2D gaw
    REAL(r_2), DIMENSION(mp_patch,mf) :: gbhf ! freeConvectionBndLayerConductance mol/m2/s
    REAL(r_2), DIMENSION(mp_patch,mf) :: gbhu ! forcedConvectionBoundaryLayerConductance
    REAL(r_2), DIMENSION(mp_patch) :: gbvtop ! bnd layer cond. top leaf
    REAL(r_1), DIMENSION(mp_patch,mf) :: gras ! Grashof coeff
    REAL(r_1), DIMENSION(mp_patch,mf) :: gswmin ! min stomatal conductance
    REAL(r_2), DIMENSION(mp_patch,mf) :: gswx ! stom cond for water
    REAL(r_2), DIMENSION(mp_patch,mf) :: gw  ! cond for water for a dry canopy
    REAL(r_2), DIMENSION(mp_patch,mf) :: gh  ! cond for heat for a dry canopy
    REAL(r_2), DIMENSION(mp_patch,mf) :: ghr ! dry canopy cond for heat & thermal radiat'n
    REAL(r_2), DIMENSION(mp_patch) :: gwwet  ! cond for water for a wet canopy
    REAL(r_2), DIMENSION(mp_patch) :: ghwet  ! cond for heat for a wet canopy
    REAL(r_2), DIMENSION(mp_patch) :: ghrwet ! wet canopy cond: heat & thermal radiat'n
    REAL(r_2), DIMENSION(mp_patch,mf) :: hcx ! sens heat fl big leaf prev iteration
    REAL(r_2), DIMENSION(mp_patch,mf) :: hcy ! veg. sens heat
    INTEGER(i_d) :: iter ! iteration #
    INTEGER(i_d) :: iterplus !
    INTEGER(i_d) :: k  ! interation count
    REAL(r_1), DIMENSION(mp_patch,mf) :: psycst ! modified pych. constant
    REAL(r_1), DIMENSION(mp_patch,mf) :: rdx ! daytime leaf resp rate, prev iteration
    REAL(r_1), DIMENSION(mp_patch,mf) :: rdy ! daytime leaf resp rate
    REAL(r_2), DIMENSION(mp_patch,mf) :: rnx ! net rad prev timestep
    REAL(r_2), DIMENSION(mp_patch,mf) :: rny ! net rad
    REAL(r_1), DIMENSION(mp_patch) :: rt0 ! turbulent resistance
    REAL(r_1), DIMENSION(mp_patch) :: ortsoil ! turbulent resistance, prev time step
    REAL(r_1), DIMENSION(mp_patch) :: rt1usc ! eq. 3.53, SCAM manual, 1997
    REAL(r_1), DIMENSION(mp_patch) :: rwater ! soil water availability
    REAL(r_1), DIMENSION(mp_patch,mf) :: tair2 ! 2D tair
    REAL(r_1), DIMENSION(mp_patch,mf) :: tvair2 ! 2D tair
    REAL(r_1), DIMENSION(mp_patch,mf) :: tdiff ! leaf air temp diff.
    REAL(r_1), DIMENSION(mp_patch,mf) :: tlfx ! leaf temp prev. iteration
    REAL(r_1), DIMENSION(mp_patch,mf) :: tlfxx ! leaf temperature of current iteration
    REAL(r_1), DIMENSION(mp_patch,mf) :: tlfy ! leaf temp
    REAL(r_1), DIMENSION(mp_patch,mf) :: vcmax2 ! vcmax big leaf
    REAL(r_1), DIMENSION(mp_patch,mf) :: vcmxt3 ! vcmax big leaf C3
    REAL(r_1), DIMENSION(mp_patch,mf) :: vcmxt4 ! vcmax big leaf C4
    REAL(r_2), DIMENSION(mp_patch,mf,2):: vx3 ! carboxylation C3 plants
    REAL(r_2), DIMENSION(mp_patch,mf,2):: vx4 ! carboxylation C4 plants
    REAL(r_1), DIMENSION(mp_patch) :: wetfac ! degree of soil water limitation on stage 2 soil evaporation
    REAL(r_1), DIMENSION(mp_patch,mf) :: xdleaf2 ! 2D dleaf
    REAL(r_2), DIMENSION(mp_patch,mf) :: xleuning ! leuning stomatal coeff
    REAL(r_1), DIMENSION(mp_patch,niter):: zetar ! stability correction
    REAL(r_1), PARAMETER :: jtomol = 4.6e-6 ! Conversion from Joule to Mol for light
    REAL(r_1), PARAMETER :: EHaVc  = 73637.0  !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER :: EHdVc  = 149252.0 !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER :: EntropVc = 486.0  !J/mol/K (Leuning 2002)
    REAL(r_1), PARAMETER :: xVccoef = 1.17461 !derived parameter
    !	xvccoef=1.0+exp((Entropjx*TrefK-EHdjx)/(Rconst*TrefK))
    REAL(r_1), PARAMETER :: EHaJx  = 50300.0  !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER :: EHdJx  = 152044.0 !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER :: EntropJx = 495.0 !J/mol/K (Leuning 2002)
    REAL(r_1), PARAMETER :: xjxcoef = 1.16715 !derived parameter
    REAL(r_1), PARAMETER :: effc4 = 4000.0  !Vc=effc4*Ci*Vcmax (see
    ! Bonan,LSM version 1.0, p106)
    REAL(r_1), DIMENSION(mp_patch) :: oldcansto ! prev t step canopy storage
    REAL(r_1), DIMENSION(mp_patch) :: cc ! limitation term for canopy interception per timestep		   
    REAL(r_2), DIMENSION(mp_patch) :: ccfevw ! limitation term for wet canopy evaporation rate  
    REAL(r_1), DIMENSION(mp_patch) :: denom ! denominator in calculating screen temperature, humidity etc
    REAL(r_1), DIMENSION(mp_patch) :: tstar ! 
    REAL(r_1), DIMENSION(mp_patch) :: zscrn !
    REAL(r_1), DIMENSION(mp_patch) :: qstar !
    REAL(r_1), DIMENSION(mp_patch) :: rsts  !
    REAL(r_1), DIMENSION(mp_patch) :: qsurf !
    REAL(r_1), DIMENSION(mp_patch) :: qtgnet !
    REAL(r_1), DIMENSION(mp_patch) :: evapfb !
    REAL(r_1), DIMENSION(mp_patch,ms) :: evapfbl !
    REAL(r_1), DIMENSION(mp_patch) :: phenps ! Leaf phenology influence on vcmax and jmax
    REAL(r_1), DIMENSION(mp_patch) :: poolcoef1 ! non-leaf carbon turnover rate * non-leaf pool size
    REAL(r_1), DIMENSION(mp_patch) :: poolcoef1w ! wood carbon turnover rate * wood pool size
    REAL(r_1), DIMENSION(mp_patch) :: poolcoef1r ! root carbon turnover rate * root pool size
    REAL(r_1), DIMENSION(mp_patch) :: rbw ! leaf boundary layer resistance for water
    REAL(r_1), DIMENSION(mp_patch) :: rsw ! stomatal resistance for water
    REAL(r_1), DIMENSION(mp_patch) :: dmah ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch) :: dmbh ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch) :: dmch ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch) :: dmae ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch) :: dmbe ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch) :: dmce ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch) :: tss4 ! soil/snow temperature**4
!    REAL(r_1), DIMENSION(mp_patch) :: sss ! variable for Penman-Monteith evap for soil
!    REAL(r_1), DIMENSION(mp_patch) :: cc1 ! variable for Penman-Monteith evap for soil
!    REAL(r_1), DIMENSION(mp_patch) :: cc2 ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp_patch) :: qstvair ! sat spec humidity at leaf temperature
    REAL(r_1), DIMENSION(mp_patch) :: qstss ! sat spec humidity at soil/snow temperature
    REAL(r_1), DIMENSION(mp_patch) :: xx ! delta-type function for sparse canopy limit, p20 SCAM manual
    REAL(r_1), DIMENSION(mp_patch,mf) :: temp ! vcmax big leaf C3
    REAL(r_2), DIMENSION(mp_patch,mf) :: deltecy ! YP & Mao (jun08)
    !
    !	xjxcoef=1.0+exp((Entropjx*TrefK-EHdjx)/(Rconst*TrefK))
    ! 1-oct-2002 ypw: to keep the unit consistent for resistance or conductance
    ! s/m for r; mol/m2/s for g, and always use g where appropriate
    ! replace rh, rhr, rw  with ghdry/ghwet,ghrdry, ghrwet, gwdry, gwwet
    
    LOGICAL, DIMENSION(mp_patch) :: mdb_mask ! Needs to be set

    WHERE(model_structure%canopy=='hawkesbury') 
       mdb_mask = .TRUE. ! this needs to be set elsewhere!
       ! Note how far between 380ppm and 620ppm each site is:
       cscale = MIN(MAX((met%ca-380.0e-6)/240.0e-6,0.0),1.0) ! 0<cscale<1
       WHERE(mdb_mask) ! In the Murray Darling
          a1c3(:,1) = cscale*a1c3_620 + (1-cscale)*a1c3_380 ! scale linearly
          veg%vcmax = veg%vcmax * (1 - 0.25*cscale)
       ELSEWHERE ! use current C02 level value of a1c3
          a1c3(:,1) = a1c3_380
       END WHERE
       ! Set the same values for the shaded leaf:
       a1c3(:,2) = a1c3(:,1)
       ! Empirical coef for vpd sensitivity of stomata
       d0c3(:,1) = d0c3_hawkesbury ! set in photosynthetic constants
       d0c3(:,2) = d0c3(:,1)
    ELSEWHERE
       a1c3(:,1) = a1c3_default ! set in photosynthetic constants
       a1c3(:,2) = a1c3(:,1)
       d0c3(:,1) = d0c3_default ! set in photosynthetic constants
       d0c3(:,2) = d0c3(:,1)
    END WHERE

    ! Set surface water vapour pressure deficit:
    met%da = (qsatf(met%tc,met%pmb) - met%qv ) * rmair/rmh2o * met%pmb * 100.0
    ! Soil water limitation on stomatal conductance:
    rwater = MAX(1.0e-4, &
         SUM(veg%froot * MIN(1.0,REAL(ssoil%wb,r_1) - &
         SPREAD(soil%swilt, 2, ms)),2) / (soil%sfc-soil%swilt))
    ! construct function to limit stom cond for soil water availability
    fwsoil = MAX(1.0e-4,MIN(1.0, veg%vbeta * rwater))
    ! BATS-type canopy saturation proportional to LAI:
    cansat = veg%canst1 * canopy%vlaiw
    ! Leaf phenology influence on vcmax and jmax
    ! rml 22/10/07 only apply to deciduous types
    WHERE (veg%deciduous)
       phenps = MAX (1.0e-4, MIN(1.0, 1.0 - ( (veg%tmaxvj - ssoil%tgg(:,4)+tfrz) &
            / (veg%tmaxvj - veg%tminvj) )**2 ) )
       WHERE ( ssoil%tgg(:,4) < (veg%tminvj + tfrz) ) phenps = 0.0
       WHERE ( ssoil%tgg(:,4) > (veg%tmaxvj + tfrz) ) phenps = 1.0
    ELSEWHERE
       phenps = 1.0
    END WHERE
    ! Set previous time step canopy water storage:
    oldcansto=canopy%cansto
    ! Rainfall variable is limited so canopy interception is limited,
    ! used to stabilise latent fluxes.
    ! to avoid excessive direct canopy evaporation (EK nov2007, snow scheme)
    !      cc =MIN(met%precip-met%precip_s, 4./(1440./(dels/60.)))
    ! modified further by timestep requirement (EAK aug08)
    cc =MIN(met%precip-met%precip_s, 4.0 * MIN(dels,1800.0) / 60.0 /1440.0 )
    !    ! to avoid canopy temp. oscillations
    !    cc =MIN(met%precip, 4./(1440./(dels/60.)))
    ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
    canopy%wcint = MERGE(MIN(MAX(cansat - REAL(canopy%cansto,r_1),0.0), cc), 0.0, &
         cc > 0.0  )  ! EK nov2007, snow scheme
    !         cc > 0.0  .AND. met%tk > tfrz)
    ! Define canopy throughfall (100% of precip if temp < 0C, see above):
    canopy%through = met%precip_s + MIN( met%precip - met%precip_s , &
         MAX(0.0, met%precip - met%precip_s - canopy%wcint) )  ! EK nov2007
    ! Delete line below in case it affects snow sites (precip_s) (EK Jul08)
    !      canopy%through = MIN(met%precip,MAX(0.0, met%precip - canopy%wcint))
    ! Add canopy interception to canopy storage term:
    canopy%cansto = canopy%cansto + canopy%wcint
    wetfac = MAX(0.0, MIN(1.0, &
         (REAL(ssoil%wb(:,1),r_1) - soil%swilt) / (soil%sfc - soil%swilt)))
    ! Temporay fix to account for reduction of soil evaporation when freezing
    wetfac = wetfac * REAL(((1.0-ssoil%wbice(:,1)/ssoil%wb(:,1)))**2,r_1)
    zetar(:,1) = zeta0 ! stability correction terms
    zetar(:,2) = zetpos + 1
    xdleaf2 = SPREAD(veg%dleaf, 2, mf) ! characteristic leaf length
    dsatdk2 = SPREAD(air%dsatdk, 2, mf)! deriv of sat vap pressure wrt temp
    ca2 = SPREAD(met%ca, 2, mf)        ! CO2 concentration
    csx = ca2                     ! initialise leaf surface CO2 concentration
    da2 = SPREAD(met%da, 2, mf)   ! water vapour pressure deficit
    dsx = da2                     ! init. leaf surface vpd
    tair2 = SPREAD(met%tc, 2, mf) ! air temp in C
    ejmax2 = SPREAD(veg%ejmax*phenps, 2,mf) ! max. pot. electr transp. rate of top leaf(mol/m2s)
    vcmax2 = SPREAD(veg%vcmax*phenps, 2,mf) ! max. RuBP carboxylsation rate of top leaf(mol/m2s)
    tlfx = tair2  ! initialise leaf temp iteration memory variable
    tlfy = tair2  ! initialise current leaf temp
    frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants
    ! weight min stomatal conductance by C3 an C4 plant fractions
    rdy = 0.0       ! init daytime leaf respiration rate
    rdx = 0.0       ! init daytime leaf respiration rate
    an_y = 0.0      ! init current estimate net photos.
    gswx = 1e-3     ! default stomatal conuctance 
    gbhf = 1e-3     ! default free convection boundary layer conductance
    gbhu = 1e-3     ! default forced convection boundary layer conductance
    anx = 0.0       ! init net photos. iteration memory variable
    psycst = SPREAD(air%psyc, 2, mf) ! modified psyc. constant
    ! add on by ypw 1-oct-2002
    gw = 1.0e-3 ! default values of conductance
    gh = 1.0e-3
    ghr= 1.0e-3
    ancj = 0.0
    ghwet = 1.0e-3
    gwwet = 1.0e-3
    ghrwet= 1.0e-3
    ! Initialise in-canopy temperatures and humidity:
    met%tvair = met%tk
    met%tvrad = met%tk
    met%qvair = met%qv
    ortsoil = ssoil%rtsoil
    ssoil%tss = (1-ssoil%isflag)*ssoil%tgg(:,1)+ssoil%isflag*ssoil%tggsn(:,1)
    tss4 = ssoil%tss**4
    DO iter = 1, niter
       CALL define_air (met, air)
       psycst = SPREAD(air%psyc, 2, mf)
       dsatdk2 = SPREAD(air%dsatdk, 2, mf)
       CALL radiation(ssoil, veg, air, met, rad, canopy)
       hcx = 0.0       ! init sens heat iteration memory variable
       ecx = rad%rniso ! init lat heat iteration memory variable
       rnx = rad%rniso ! init net rad iteration memory variable
       rny = rad%rniso ! init current estimate net rad
       hcy = 0.0       ! init current estimate lat heat
       ecy = rny - hcy ! init current estimate lat heat
       gswmin = rad%scalex * (gsw03 * (1.0 - frac42) + gsw04 * frac42)
       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       canopy%us = MAX( 1.0e-6, vonk * MAX(met%ua,umin) &
            / ( LOG(rough%zref/rough%z0m) - psim(zetar(:,iter)) &
            + psim(zetar(:,iter)*rough%z0m/rough%zref) ) )
       ! Turbulent aerodynamic resistance from roughness sublayer depth to
       ! reference height,
       ! x=1 if zref+disp>zruffs, 0 otherwise:
       ! thus rt1usc = 0 if zref+disp<zruffs
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
       xx = 0.5 + SIGN(0.5,rough%zref+rough%disp-rough%zruffs)
       !              correction  by Ian Harman to the 2nd psis term
       rt1usc = xx * ( LOG( rough%zref/MAX(rough%zruffs-rough%disp, &
            rough%z0soilsn) ) - psis(zetar(:,iter)) &
            + psis( zetar(:,iter)*MAX(rough%zruffs-rough%disp, &
            rough%z0soilsn)/rough%zref ) )/vonk
       ! rt0 = turbulent resistance from soil to canopy:
!!$     ! correction  by Ian Harman to rough%rt0us = f( zetar )
!!$     WHERE (canopy%vlaiw.LT.0.01 .OR. rough%hruff.LT. rough%z0soilsn)
!!$       rough%rt0us  = 0.0
!!$       rt0old  = 0.0
!!$     ELSEWHERE
!!$!      rough%term6=exp(3.*rough%coexp*(rough%disp/rough%hruff -1.))
!!$       rt0old  = rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) &
!!$            + (1-zdlin))*(EXP(2*csw*canopy%vlaiw) - rough%term2)/rough%term3
!!$       rough%rt0us=rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) &
!!$!          - psis( zetar(:,iter) * rough%disp/rough%zref/rough%term6)  &
!!$!          + psis( zetar(:,iter) * rough%z0soilsn/rough%zref/rough%term6) &
!!$           + (1-zdlin))*(EXP(2*csw*canopy%vlaiw) - rough%term2)/rough%term3 &
!!$              / rough%term6
!!$     ENDWHERE
!!$     rt0old = rt0old / canopy%us
!!$     rt0 = MAX(5.,rough%rt0us / canopy%us)
       rt0 = rough%rt0us / canopy%us
       ! Aerodynamic resistance (sum 3 height integrals)/us
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
       rough%rt1 = MAX(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
       WHERE (ssoil%snowd > 0.1)
          wetfac = 1.
       END WHERE
       ssoil%rtsoil = rt0 + rough%rt1*(0.5 + SIGN(0.5,0.02-canopy%vlaiw)) 
       ssoil%rtsoil = MAX(25.,ssoil%rtsoil)   
       WHERE (ssoil%rtsoil > 2.*ortsoil .OR. ssoil%rtsoil < 0.5*ortsoil)
          ssoil%rtsoil = MAX(25.,0.5*(ssoil%rtsoil + ortsoil))
       END WHERE

       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! prandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       gbvtop = air%cmolar*apol * air%visc / prandt / veg%dleaf &
            * (canopy%us / MIN(rough%usuh, 0.2) &
            * veg%dleaf / air%visc)**0.5 * prandt**(1.0/3.0) / veg%shelrb
       ! Forced convection boundary layer conductance
       ! (see Wang & Leuning 1998, AFM):
       ! gbhu corrected by F.Cruz & A.Pitman on 13/03/07
       gbhu(:,1) = gbvtop*(1.0-EXP(-canopy%vlaiw*(0.5*rough%coexp+rad%extkb)))&
            / (rad%extkb+0.5*rough%coexp)
       gbhu(:,2) = (2.0/rough%coexp)*gbvtop*  &
            (1.0 - EXP(-0.5*rough%coexp*canopy%vlaiw))-gbhu(:,1)

       ! Aerodynamic conductance:
       gaw = air%cmolar / rough%rt1
       WHERE (veg%meth > 0 )
          gaw=100000.0
       END WHERE
       gaw2 = SPREAD(gaw, 2, mf)
       abs_deltlf = 999.0
       deltlfy = 999.0
       WHERE(rad%fvlai <= 1.0e-2)
          abs_deltlf=0.0
          hcx = 0.0 ! intialise
          ecx = 0.0 ! intialise
          anx = 0.0 ! intialise
          rnx = 0.0 ! intialise
          rny = rnx ! store initial values
          hcy = hcx ! store initial values
          ecy = ecx ! store initial values
          rdy = rdx ! store initial values
          an_y = anx ! store initial values
       END WHERE
       deltlfy = abs_deltlf
       k = 0
       DO WHILE (ANY(abs_deltlf > 0.1)  .AND.  k < maxiter)
          k = k + 1
          WHERE (rad%fvlai > 1e-2 .AND. abs_deltlf > 0.1)
             ! Grashof number (Leuning et al, 1995) eq E4:
             gras = MAX(1.0e-6,1.595E8*ABS(tlfx-tair2)*(xdleaf2**3.))
             ! See Appendix E in (Leuning et al, 1995):
             gbhf = rad%fvlai * SPREAD(air%cmolar, 2, mf) * 0.5 * dheat &
                  * (gras**0.25) / xdleaf2
             ! Conductance for heat:
             gh = 1.0/(MIN(1.0e3_r_2, SPREAD(1.0/gaw, 2, mf) + 0.5/(gbhu+gbhf)))
             ! Conductance for heat and longwave radiation:
             ghr = rad%gradis+gh
             temp =  xvcmxt3(tlfx+tfrz)
             !  Leuning 2002 (P C & E) equation for temperature response
             !  used for Vcmax for C3 plants:
             vcmxt3 = (1.0-frac42)*vcmax2*rad%scalex *temp
             temp=  xvcmxt4(tlfx)
             ! Temperature of Vcmax for C4 plants (Collatz et al 1989):
             vcmxt4 = frac42 * vcmax2 * rad%scalex * temp
             temp= xejmxt3(tlfx+tfrz)
             !  Leuning 2002 (P C & E) equation for temperature response
             !  used for Jmax for C3 plants:
             ejmxt3 = (1.0-frac42) * ejmax2 * rad%scalex * temp
             ! Difference between leaf temperature and reference temperature:
             tdiff  =  tlfx+tfrz-trefk
             ! Michaelis menten constant of Rubisco for CO2:
             conkct = conkc0*EXP((ekc/(rgas* trefk)) *(1.-trefk/(tlfx+tfrz)))
             ! Michaelis menten constant of Rubisco for oxygen:
             conkot = conko0*EXP((eko/(rgas* trefk)) *(1.-trefk/(tlfx+tfrz)))
             ! "d_{3}" in Wang and Leuning, 1998, appendix E:
             cx(:,:,1) = conkct*(1.0+0.21/conkot)
             cx(:,:,2) = 2.0* gam0*(1.+gam1*tdiff + gam2*tdiff*tdiff) !gamma*
             ! All equations below in appendix E in Wang and Leuning 1998 are
             ! for calculating anx, csx and gswx for Rubisco limited,
             ! RuBP limited, sink limited
             vx3(:,:,1) = vcmxt3
             vx4(:,:,1) = vcmxt4
             temp = rad%qcan(:,:,1)*jtomol*(1.0-frac42)
             vx3(:,:,2) = ej3x(temp,ejmxt3)
             temp = frac42*rad%qcan(:,:,1)*jtomol
             vx4(:,:,2) = ej4x(temp,vcmxt4)
             rdx = (cfrd3 * vcmxt3 + cfrd4 * vcmxt4) * SPREAD(fwsoil, 2, mf)
             xleuning = (1.0-frac42) * a1c3 / (1.0 + dsx/d0c3) &
                  + frac42 * a1c4 /(1.0 + dsx/d0c4)
             xleuning = xleuning * SPREAD(fwsoil, 2, mf) / (csx-co2cp3)
             ! Rubisco limited:
             coef2(:,:,1) = gswmin/rgswc+xleuning *(vx3(:,:,1)-(rdx-vx4(:,:,1)))
             coef1(:,:,1) = (1.0-csx*xleuning) *(vx3(:,:,1)+vx4(:,:,1)-rdx) &
                  +(gswmin/rgswc)*(cx(:,:,1)-csx) -xleuning*(vx3(:,:,1) &
                  *cx(:,:,2)/2.0+cx(:,:,1)*(rdx-vx4(:,:,1)))
             coef0(:,:,1) = -(1.0-csx*xleuning) *(vx3(:,:,1)*cx(:,:,2)/2.0  &
                  +cx(:,:,1)*(rdx-vx4(:,:,1))) -(gswmin/rgswc)*cx(:,:,1)*csx
             ! Discriminant in quadratic in eq. E7 Wang and Leuning, 1998
             delcx(:,:,1) = coef1(:,:,1)**2 -4.0*coef0(:,:,1)*coef2(:,:,1)
             ci(:,:,1) = (-coef1(:,:,1)+SQRT(MAX(0.0_r_2,delcx(:,:,1)))) &
                  / (2.0*coef2(:,:,1))
             ci(:,:,1) = MAX(0.0_r_2,ci(:,:,1))
             ancj(:,:,1) = vx3(:,:,1)*(ci(:,:,1)-cx(:,:,2)/2.0) &
                  / (ci(:,:,1) + cx(:,:,1)) + vx4(:,:,1) - rdx
             ! RuBP limited:
             coef2(:,:,2) = gswmin/rgswc+xleuning *(vx3(:,:,2)-(rdx-vx4(:,:,2)))
             coef1(:,:,2) = (1.0-csx*xleuning) *(vx3(:,:,2)+vx4(:,:,2)-rdx) &
                  +(gswmin/rgswc)*(cx(:,:,2)-csx) -xleuning*(vx3(:,:,2) &
                  *cx(:,:,2)/2.0+cx(:,:,2)*(rdx-vx4(:,:,2)))
             coef0(:,:,2) = -(1.0-csx*xleuning) *(vx3(:,:,2)*cx(:,:,2)/2.0  &
                  +cx(:,:,2)*(rdx-vx4(:,:,2))) -(gswmin/rgswc)*cx(:,:,2)*csx
             delcx(:,:,2) = coef1(:,:,2)**2 -4.0*coef0(:,:,2)*coef2(:,:,2)
             ci(:,:,2) = (-coef1(:,:,2) + SQRT(MAX(0.0_r_2,delcx(:,:,2)))) &
                  / (2.0*coef2(:,:,2))
             ci(:,:,2) = MAX(0.0_r_2,ci(:,:,2))
             ancj(:,:,2) = vx3(:,:,2) * (ci(:,:,2)-cx(:,:,2)/2.0) &
                  / (ci(:,:,2)+cx(:,:,2)) + vx4(:,:,2) - rdx
             ! Sink limited:
             coef2(:,:,3) = xleuning
             coef1(:,:,3) = gswmin/rgswc + xleuning * (rdx - 0.5*vcmxt3)  &
                  + effc4 * vcmxt4 - xleuning * csx * effc4 *vcmxt4
             coef0(:,:,3) = -(gswmin/rgswc)*csx *effc4*vcmxt4 &
                  + (rdx -0.5*vcmxt3)*gswmin/rgswc
             delcx(:,:,3) = coef1(:,:,3)**2 -4.0*coef0(:,:,3)*coef2(:,:,3)
             ancj(:,:,3)  = (-coef1(:,:,3)+SQRT(MAX(0.0_r_2,delcx(:,:,3)))) &
                  / (2.0*coef2(:,:,3))
             anx = MIN(ancj(:,:,1),ancj(:,:,2),ancj(:,:,3))
             csx = ca2 - anx * (1.0/gaw2+rgbwc/(gbhu + gbhf))
             gswx = gswmin+MAX(0.0_r_2,rgswc*xleuning *anx)
             ! Recalculate conductance for water:
             gw = 1.0/(1.0/gswx+1.0/(1.075*(gbhu+gbhf)) + SPREAD(1.0/gaw, 2, mf))
             ! Modified psychrometric constant (Monteith and Unsworth, 1990)
             psycst = SPREAD(air%psyc, 2, mf) * REAL( ghr / gw ,r_1)
             ! Store leaf temperature:
             tlfxx = tlfx
             ! Update canopy latent heat flux:
             ecx = (dsatdk2*rad%rniso +capp*rmair*da2*ghr) /(dsatdk2+psycst)
             ! Update canopy sensible heat flux:
             hcx = (rad%rniso-ecx)*gh/ghr
             ! Update leaf temperature:
             tlfx=tair2+REAL(hcx,r_1)/(capp*rmair*REAL(gh,r_1))
             ! Update net radiation for canopy:
             rnx = rad%rniso-capp*rmair*(tlfx -tair2)*rad%gradis
             ! Update leaf surface vapour pressure deficit:
             ! dsx = ecx*100.0* SPREAD(met%pmb, 2, mf) &
             !      /(gswx*rmh2o*SPREAD(air%rlam, 2, mf))
             dsx = da2 + dsatdk2 * (tlfx-tair2)
             ! Store change in leaf temperature between successive iterations:
             deltlf = tlfxx-tlfx
             abs_deltlf = ABS(deltlf)
          END WHERE
          ! Where leaf temp change b/w iterations is significant, and difference
          ! is smaller than the previous iteration, store results:
          WHERE (abs_deltlf > 0.1 .AND. abs_deltlf < ABS(deltlfy) )
             deltlfy = deltlf
             tlfy = tlfx
             rny = rnx
             hcy = hcx
             ecy = ecx
             rdy = rdx
             an_y = anx
          END WHERE
          WHERE (abs_deltlf > 0.1)
             ! after four iterations, take the mean value of current and previous
             ! estimates as the next estimate of leaf temperature, to avoid
             ! oscillation
             tlfx = (0.5*(MAX(0,k-5)/(k-4.9999))) *tlfxx &
                  & + (1.0- (0.5*(MAX(0,k-5)/(k-4.9999))))*tlfx
          END WHERE
          IF(k==1) THEN
             !        taken the first iterated estimates as the defaults
             tlfy = tlfx
             rny = rnx
             hcy = hcx
             ecy = ecx
             rdy = rdx
             an_y = anx
          END IF
       END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND. k < maxiter)
       
       ! VEGETATION SENSIBLE AND LATENT HEAT FLUXES fevw, fhvw (W/m2) for a 
       ! wet canopy
       ! calculate total thermal resistance, rthv in s/m
       ghwet=gaw*SUM(gbhu,2)/(0.5*gaw+SUM(gbhu,2))
       gwwet=1.075*SUM(gbhu,2)*gaw/(1.075*SUM(gbhu,2)+gaw)
       WHERE (veg%meth > 0 )
          ghwet=2.*SUM(gbhu,2)
          gwwet=1.075*SUM(gbhu,2)
       END WHERE
       ghrwet=SUM(rad%gradis,2)+ghwet
       ! Calculate fraction of canopy which is wet:
       canopy%fwet   = MAX(0.0,MIN(1.0,0.8*canopy%cansto/MAX(cansat,0.01)))
       ! Calculate lat heat from wet canopy, may be neg. if dew onto wet canopy
       ! to avoid excessive evaporation:
       ccfevw = MIN(canopy%cansto*air%rlam/dels,2./(1440./(dels/60.))*air%rlam)
       canopy%fevw = MIN(canopy%fwet*((air%dsatdk*SUM(rad%rniso,2)+ &
            capp*rmair*met%da*ghrwet) &
            /(air%dsatdk+air%psyc*ghrwet/gwwet)), ccfevw)
       ! canopy potential evapotranspiratn for output purposes (YP & Mao jun08)
       canopy%potev_c = MAX(0.0_r_2,(air%dsatdk*SUM(rad%rniso,2)+ &
            capp*rmair*met%da*ghrwet)/(air%dsatdk+air%psyc*ghrwet/gwwet))
       ! Calculate sens heat from wet canopy:
       canopy%fhvw = (canopy%fwet*SUM(rad%rniso,2)-canopy%fevw)*ghwet/ghrwet
       ! Calculate (dry) canopy transpiration, may be negative if dew
       canopy%fevc = (1.0 - canopy%fwet) * SUM(ecy,2)
       evapfbl = 0.
       DO k = 1,ms
          WHERE(canopy%fevc > 0.0)
             evapfb = REAL(canopy%fevc,r_1) * dels/air%rlam ! convert to mm/dt
             evapfbl(:,k) = MIN(evapfb*veg%froot(:,k), &
                  MAX(0.,MIN(REAL(ssoil%wb(:,k),r_1)-soil%swilt, &
                  REAL(ssoil%wb(:,k)-1.05*ssoil%wbice(:,k),r_1))) &
                  * soil%zse(k)*1000.0)
          END WHERE
       END DO
       WHERE(SUM(ecy,2)>0.0)
          canopy%fevc=SUM(evapfbl,2)*air%rlam/dels
       END WHERE
       
      ! replaced 6 lines below with 3 lines above to keep negative fevc values
      ! (YP & Mao jun08)
!      canopy%fevc = 0.
!      DO k = 1,ms
!        WHERE(activepatch) ! where vegetation/soil patch > 0%
!          canopy%fevc=canopy%fevc+ evapfbl(:,k)*air%rlam/dels
!        END WHERE
!      END DO
      tvair2 = SPREAD(met%tvair-tfrz, 2, mf) ! within-canopy air temp in C
      ! Calculate latent heat from vegetation:
      canopy%fev = REAL(canopy%fevc + canopy%fevw,r_1)
      ! recalculate for checking energy balance (YP & Mao, jun08)
      deltecy(:,1) = (canopy%fevc/(MAX((1.0-canopy%fwet),1.0e-10)))*rad%fvlai(:,1) &
           / (rad%fvlai(:,1)+rad%fvlai(:,2)+1.0e-10)
      deltecy(:,2) = (canopy%fevc/(MAX((1.0-canopy%fwet),1.0e-10)))*rad%fvlai(:,2) &
           / (rad%fvlai(:,1)+rad%fvlai(:,2)+1.0e-10)
      ecy(:,1)     = deltecy(:,1)
      ecy(:,2)     = deltecy(:,2)
      hcy(:,1)     = (rad%rniso(:,1) - ecy(:,1))*gh(:,1)/ghr(:,1)
      hcy(:,2)     = (rad%rniso(:,2) - ecy(:,2))*gh(:,2)/ghr(:,2)
      tlfy(:,1)    = tvair2(:,1)+REAL(hcy(:,1),r_1)/(capp*rmair*REAL(gh(:,1),r_1))
      tlfy(:,2)    = tvair2(:,2)+REAL(hcy(:,2),r_1)/(capp*rmair*REAL(gh(:,2),r_1))
      rny(:,1)     = rad%rniso(:,1) - capp*rmair * (tlfy(:,1) &
           - tvair2(:,1)) * rad%gradis(:,1)
      rny(:,2)     = rad%rniso(:,2) - capp*rmair * (tlfy(:,2) &
           - tvair2(:,2)) * rad%gradis(:,2)
      ! end recalculation (YP & Mao, jun08)
      ! Calculate sensible heat from vegetation:
      canopy%fhv = (1.0 - canopy%fwet) * REAL(SUM(hcy,2) + canopy%fhvw, r_1)
      ! Calculate net rad absorbed by canopy:
      canopy%fnv = (1.0 - canopy%fwet) * REAL(SUM(rny,2) &
           + canopy%fevw + canopy%fhvw , r_1)
      ! canopy radiative temperature is calculated based on long-wave
      ! radiation balance
      !   Q_lw=Q_lw_iso - (1.0-fwet)*SUM(capp*rmair*(tlfy-tair)*gri &
      !                 - canopy%fhvw*gr/ghw
      !   Q_lw=(1-transd)*(L_soil+L_sky-2.0*L_can)
      ! therefore
      !   Q_lw_iso-Q_lw=2(1-transd)*emleaf*(Tv^4-Tc^4)
      
      !rad%lwabv = (1.0-canopy%fwet)*(capp*rmair*(tlfy(:,1) &
      !          - met%tc)*rad%gradis(:,1) &
      !          + capp*rmair*(tlfy(:,2) - met%tc)*rad%gradis(:,2)) &
      !          + canopy%fhvw*SUM(rad%gradis,2)/ghwet
      rad%lwabv = (1.0-canopy%fwet)*(capp*rmair*(tlfy(:,1) &
           - tvair2(:,1))*REAL(rad%gradis(:,1),r_1) &
           + capp*rmair*(tlfy(:,2) - tvair2(:,2))*REAL(rad%gradis(:,2),r_1)) &
           !            YP & Mao (jun08) replaced met%tk with tvair2
           !             - (met%tk-tfrz))*rad%gradis(:,1) &
           !             + capp*rmair*(tlfy(:,2) - (met%tk-tfrz))*rad%gradis(:,2)) &
           + REAL(canopy%fhvw*SUM(rad%gradis,2)/ghwet,r_1)
      ! add if condition here to avoid dividing by zero
      ! ie when rad%transd=1.0 Ypw:24-02-2003
      WHERE (canopy%vlaiw > 0.01 .AND. rough%hruff > rough%z0soilsn)
         canopy%tv = (rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf) &
              + met%tvair**4)**0.25
         !                 YP & Mao (jun08) replaced met%tk with tvair
         !                  & + met%tk**4)**0.25
      ELSEWHERE ! sparse canopy
         canopy%tv = met%tvair
         !         YP & Mao (jun08) replaced met%tk with tvair
         !          canopy%tv = met%tk
      END WHERE

      ! Calculate ground heat flux:
      canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux &
           + ssoil%isflag*canopy%sghflux
      ! Saturation specific humidity at soil/snow surface temperature:
      qstss = qsatf((ssoil%tss - tfrz),met%pmb)
      ! Spec hum deficit at soil/snow surface:
      dq = qstss - met%qv
      ! excessive dew over snow area
      WHERE (ssoil%snowd > 0.1)
         dq = MAX( -0.1e-3, dq)
      END WHERE
      ! Calculate net rad to soil:
      canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd) &
           *emleaf*sboltz*canopy%tv**4 - emsoil*sboltz* tss4
      ! Penman-Monteith formula
      !        sss=air%dsatdk
      !        cc1=sss/(sss+air%psyc )
      !        cc2=air%psyc /(sss+air%psyc )
      !        ssoil%potev = cc1 * (canopy%fns - canopy%ghflux) + cc2 * air%rho &
      !          & * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
      ! method alternative to P-M formula above
      ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
      ! Soil latent heat:
      canopy%fes= wetfac * ssoil%potev
      WHERE (ssoil%snowd < 0.1 .AND. canopy%fes > 0.0)
         ! Reduce for wilting point limitation:
         canopy%fes= MIN( canopy%fes, &
              MAX( 0.0, (REAL(ssoil%wb(:,1),r_1)-soil%swilt) ) &
              *soil%zse(1)*1000.0*air%rlam/dels )
         ! Reduce for soil ice limitation:
         canopy%fes = MIN(canopy%fes,REAL(ssoil%wb(:,1)-ssoil%wbice(:,1),r_1) &
              * soil%zse(1) * 1000.0 * air%rlam / dels)
      END WHERE
      ssoil%cls=1.
      WHERE (ssoil%snowd >= 0.1)
         ssoil%cls =1.1335
         canopy%fes=MIN(wetfac*ssoil%potev,ssoil%snowd/dels*air%rlam*ssoil%cls)
      END WHERE
      ! Calculate soil sensible heat:
      canopy%fhs = air%rho*capp*(ssoil%tss - met%tk) /ssoil%rtsoil
      ! Calculate ground heat flux:
      canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
      ! Calculate total latent heat:
      canopy%fe = canopy%fev + canopy%fes
      ! Calculate total sensible heat:
      canopy%fh = canopy%fhv + canopy%fhs

      ! Initialise in-canopy temperature and humidity:
      met%tvair = met%tk
      met%qvair = met%qv
     
      WHERE ( veg%meth > 0 .AND. canopy%vlaiw > 0.01 .AND. &
           rough%hruff > rough%z0soilsn ) 
         ! use the dispersion matrix (DM) to find the air temperature 
         ! and specific humidity (Raupach, Finkele and Zhang 1997, pp 17)
         ! leaf boundary layer resistance for water
         rbw = air%cmolar/REAL(SUM(gbhu+gbhf,2),r_1)
         ! leaf stomatal resistance for water
         rsw = air%cmolar/REAL(SUM(gswx,2),r_1)
         ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmah = (rt0+rough%rt1) * ((1.+air%epsi)/rsw + 1.0/rbw) &
              + air%epsi * (rt0*rough%rt1) / (rbw*rsw)
         ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmbh = (-air%rlam/capp)*(rt0*rough%rt1)/(rbw*rsw)
         ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmch = ((1.+air%epsi)/rsw + 1.0/rbw) * rt0 * rough%rt1 &
              * (canopy%fhv+canopy%fhs) / (air%rho*capp)
         ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmae = (-air%epsi*capp/air%rlam)*(rt0*rough%rt1)/(rbw*rsw)
         ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmbe = (rt0 + wetfac*rough%rt1) * ((1.+air%epsi)/rsw + 1.0/rbw) &
              + (rt0*rough%rt1) / (rbw*rsw)
         ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmce = ((1.+air%epsi)/rsw + 1.0/rbw) * rt0 * rough%rt1 &
              * (canopy%fev + canopy%fes) / (air%rho*air%rlam)
!ccc Vanessa says this iterative update of the temperature doesn't always converge.
! She thinks it is because the in canopy temperature should be updated at the same time
! as the soil temperature. Now, it's using old soil temperature in the canopy iterations.
! Comment update for the moement
!         ! Within canopy air temperature:
!         met%tvair = met%tk &
!              + (dmbe*dmch-dmbh*dmce)/(dmah*dmbe-dmae*dmbh+1.0e-12)
!         ! Within canopy specific humidity:
!         met%qvair = met%qv &
!              + (dmah*dmce-dmae*dmch)/(dmah*dmbe-dmae*dmbh+1.0e-12)
!         met%qvair = MAX(0.0,met%qvair)
      END WHERE

      ! Saturated specific humidity in canopy:
      qstvair = qsatf((met%tvair-tfrz),met%pmb)
      ! Saturated vapour pressure deficit in canopy:
      met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.
      ! Set radiative temperature as within canopy air temp:
      met%tvrad = met%tvair
      ! recalculate air%dsatdk (EAK aug08)
      CALL define_air (met, air)
      dsatdk2 = SPREAD(air%dsatdk, 2, mf)
      ! 2 Dim saturated vapour pressure deficit in canopy:
      dva2 = SPREAD(met%dva, 2, mf)
      ! 2 dim Within canopy air temperature in degrees C:
      tvair2 = SPREAD(met%tvair-tfrz, 2, mf)
      ! recalculate using canopy within temperature
      !     where (veg%meth .eq. 0 )
      WHERE (rad%fvlai > 1e-2) ! where LAI of sunlit or shaded leaf
         ! is significant:
         ! Recalculate fluxes and leaf temperature using within canopy air vpd:
         ecy = (dsatdk2*rad%rniso +capp*rmair*dva2*ghr) /(dsatdk2+psycst)
         hcy = (rad%rniso-ecy)*gh/ghr
         !tlfx=tvair+hcx/(capp*rmair*gh)
         tlfy=tvair2+REAL(hcy,r_1)/(capp*rmair*REAL(gh,r_1))
         ! YP & Mao (jun08) added rny calculation here
         rny = rad%rniso - capp*rmair * (tlfy - tvair2) * rad%gradis
      END WHERE

      WHERE(veg%meth > 0)
         canopy%fevc = (1.0 - canopy%fwet) * SUM(ecy,2)
         WHERE (SUM(ecy,2) > 0.0)  ! +ve ecy == +ve fevc (YP & Mao jun08)
            evapfb = REAL(canopy%fevc,r_1)*dels/air%rlam ! convert to mm/dt
            ! Calcualte contribution by different soil layers
            ! to canopy transpiration:
            evapfbl(:,1) = MIN( evapfb*veg%froot(:,1), &
                 MAX(0.,REAL(MIN(ssoil%wb(:,1)-soil%swilt, &
                 ssoil%wb(:,1)-1.05*ssoil%wbice(:,1)),r_1))*soil%zse(1)*1000.)
            evapfbl(:,2) = MIN( evapfb*veg%froot(:,2), &
                 MAX(0.,REAL(MIN(ssoil%wb(:,2)-soil%swilt, &
                 ssoil%wb(:,2)-1.05*ssoil%wbice(:,2)),r_1))*soil%zse(2)*1000.)
            evapfbl(:,3) = MIN( evapfb*veg%froot(:,3), &
                 MAX(0.,REAL(MIN(ssoil%wb(:,3)-soil%swilt, &
                 ssoil%wb(:,3)-1.05*ssoil%wbice(:,3)),r_1))*soil%zse(3)*1000.)
            evapfbl(:,4) = MIN( evapfb*veg%froot(:,4), &
                 MAX(0.,REAL(MIN(ssoil%wb(:,4)-soil%swilt, &
                 ssoil%wb(:,4)-1.05*ssoil%wbice(:,4)),r_1))*soil%zse(4)*1000.)
            evapfbl(:,5) = MIN( evapfb*veg%froot(:,5), &
                 MAX(0.,REAL(MIN(ssoil%wb(:,5)-soil%swilt, &
                 ssoil%wb(:,5)-1.05*ssoil%wbice(:,5)),r_1))*soil%zse(5)*1000.)
            evapfbl(:,6) = MIN( evapfb*veg%froot(:,6), &
                 MAX(0.,REAL(MIN(ssoil%wb(:,6)-soil%swilt, &
                 ssoil%wb(:,6)-1.05*ssoil%wbice(:,6)),r_1))*soil%zse(6)*1000.)
            ! fevc recalculated within WHERE construct (YP & Mao jun08)
            ! so that negative fevc values can be retained
            canopy%fevc = ( evapfbl(:,1)+evapfbl(:,2)+evapfbl(:,3)+evapfbl(:,4) &
                 +evapfbl(:,5)+evapfbl(:,6) ) * air%rlam/dels
         END WHERE
         ! Set total vegetation latent heat:
         canopy%fev = REAL(canopy%fevc + canopy%fevw, r_1)
         ! recalculate for checking energy balance (YP & Mao, jun08)
         deltecy(:,1) = (canopy%fevc/(MAX((1.0-canopy%fwet),1.0e-10))) &
              * rad%fvlai(:,1) / (rad%fvlai(:,1)+rad%fvlai(:,2)+1.0e-10)
         deltecy(:,2) = (canopy%fevc/(MAX((1.0-canopy%fwet),1.0e-10))) &
              * rad%fvlai(:,2) / (rad%fvlai(:,1)+rad%fvlai(:,2)+1.0e-10)
         ecy(:,1)     = deltecy(:,1)
         ecy(:,2)     = deltecy(:,2)
         hcy(:,1)     = (rad%rniso(:,1) - ecy(:,1))*gh(:,1)/ghr(:,1)
         hcy(:,2)     = (rad%rniso(:,2) - ecy(:,2))*gh(:,2)/ghr(:,2)
         tlfy(:,1)    = tvair2(:,1)+REAL(hcy(:,1),r_1)/(capp*rmair*REAL(gh(:,1),r_1))
         tlfy(:,2)    = tvair2(:,2)+REAL(hcy(:,2),r_1)/(capp*rmair*REAL(gh(:,2),r_1))
         rny(:,1)     = rad%rniso(:,1) - capp*rmair * (tlfy(:,1) &
              - tvair2(:,1)) * rad%gradis(:,1)
         rny(:,2)     = rad%rniso(:,2) - capp*rmair * (tlfy(:,2) &
              - tvair2(:,2)) * rad%gradis(:,2)
         ! end recalculation (YP & Mao, jun08)
         ! Set total vegetation sensible heat:
         canopy%fhv = (1.0 - canopy%fwet) * REAL(SUM(hcy,2) + canopy%fhvw, r_1)
         ! YP & Mao (jun08) added fnv calculation here
         canopy%fnv = canopy%fev + canopy%fhv
         ! Longwave absorbed by vegetation:
         ! YP & Mao (jun08) replaced tvair with tvair2
         rad%lwabv = (1.0-canopy%fwet)*(capp*rmair*(tlfy(:,1) &
              - tvair2(:,1))*REAL(rad%gradis(:,1),r_1) &
              + capp*rmair*(tlfy(:,2) - tvair2(:,2))*REAL(rad%gradis(:,2),r_1)) &
              !           & - (met%tvair-tfrz))*rad%gradis(:,1) &
              !           & + capp*rmair*(tlfy(:,2) - (met%tvair-tfrz))*rad%gradis(:,2)) &
              + REAL(canopy%fhvw*SUM(rad%gradis,2)/ghwet,r_1)
         ! Set canopy temperature:
         WHERE (rad%transd <= 0.98)
            canopy%tv = ( rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf) &
                 + met%tvair**4 )**0.25
         ELSEWHERE
            ! sparse canopy 
            canopy%tv = met%tvair
         END WHERE
         ! Ground heat flux:
         canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux &
              + ssoil%isflag*canopy%sghflux
         dq = qstss - met%qvair
         WHERE (ssoil%snowd > 0.1)
            dq = MAX( -0.1e-3, dq)
         END WHERE
         ! Net radiation absorbed by soil: 
         canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd) &
              *emleaf*sboltz*canopy%tv**4 - emsoil*sboltz* tss4
         ! YP & Mao (jun08) reinstating Penman-Monteith formula
         ! Penman-Monteith formula
         !        sss=air%dsatdk
         !        cc1=sss/(sss+air%psyc )
         !        cc2=air%psyc /(sss+air%psyc )
         !        ssoil%potev = cc1 * (canopy%fns - canopy%ghflux) + cc2 * air%rho &
         !            & * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
         ! method alternative to P-M formula above
         ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
         ! Soil evaporation:
         canopy%fes= wetfac * ssoil%potev
         WHERE (ssoil%snowd < 0.1 .AND. canopy%fes > 0.0 )
            ! Reduce for wilting point limitation:
            canopy%fes = MIN( canopy%fes , &
                 MAX( 0., (REAL(ssoil%wb(:,1),r_1)-soil%swilt)) &
                 * soil%zse(1) * 1000.0 * air%rlam / dels )
            ! Reduce for soil ice limitation:
            canopy%fes = MIN( canopy%fes , &
                 REAL( ssoil%wb(:,1)-ssoil%wbice(:,1) , r_1 ) &
                 * soil%zse(1) * 1000.0 * air%rlam / dels )
         END WHERE
         ssoil%cls=1.
         WHERE (ssoil%snowd >= 0.1)
            ssoil%cls  = 1.1335
            canopy%fes = MIN( wetfac*ssoil%potev , &
                 ssoil%snowd/dels*air%rlam*ssoil%cls)
         END WHERE
         ! Soil sensible heat:
         canopy%fhs = air%rho*capp*(ssoil%tss - met%tvair) /ssoil%rtsoil
         canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
         ! Set total latent heat:
         canopy%fe = canopy%fev + canopy%fes
         ! Set total sensible heat:
         canopy%fh = canopy%fhv + canopy%fhs
      END WHERE ! where veg%meth > 0
      
      ! monin-obukhov stability parameter zetar=zref/l
      !	recompute zetar for the next iteration, except on last iteration
      IF (iter < niter) THEN ! dont compute zetar on the last iter
         iterplus = MAX(iter+1,2)
         zetar(:,iterplus) = -( vonk*grav*rough%zref &
              * (canopy%fh+0.07*canopy%fe) ) &
              / (air%rho*capp*met%tk*canopy%us**3)
         ! case niter=2: final zetar=zetmul*zetar(2) (compute only when iter=1)
         IF (niter == 2) THEN
            zetar(:,2) = zetmul * zetar(:,2)
            WHERE (met%fsd ==  0.0)
               zetar(:,2) = 0.5 * zetar(:,2)
            END WHERE
         END IF
         ! constrain zeta to zetpos and zetneg (set in param0)
         zetar(:,iterplus) = MIN(zetpos,zetar(:,iterplus))  ! zetar too +
         zetar(:,iterplus) = MAX(zetneg,zetar(:,iterplus))  ! zetar too -
      END IF
   END DO ! DO iter = 1, niter
 
   ! screen temp., windspeed and relative humidity at 1.8m
   tstar = - canopy%fh / ( air%rho*capp*canopy%us)
   qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
   zscrn = MAX(rough%z0m,1.8-rough%disp)
   denom = ( LOG(rough%zref/zscrn)- psim(zetar(:,iterplus)) &
        + psim(zetar(:,iterplus) * zscrn / rough%zref) ) /vonk
   ! Calculate screen temperature:
   canopy%tscrn = met%tc - tstar * denom
   rsts = qsatf(canopy%tscrn, met%pmb)
   qtgnet = rsts * wetfac - met%qv
   canopy%cduv = canopy%us * canopy%us / (MAX(met%ua,umin))**2 ! EK jun08
   !      canopy%cduv = canopy%us * canopy%us / MAX(met%ua,umin)
   WHERE (qtgnet > 0.0)
      qsurf = rsts * wetfac
   ELSEWHERE
      qsurf = 0.1*rsts*wetfac + 0.9*met%qv
   END WHERE
   canopy%qscrn = qsurf + qstar * denom
   canopy%uscrn = MAX(0.0, MAX(met%ua,umin) - canopy%us*denom ) ! at present
   ! incorrect
   ! avgwrs = REAL(SUM(veg%froot * ssoil%wb,2),r_1)
   ! avgtrs = MAX(0.0,SUM(veg%froot * ssoil%tgg,2)-tfrz)
   poolcoef1=(SUM(spread(bgc%ratecp,1,mp_patch)*bgc%cplant,2) - &
        bgc%ratecp(1)*bgc%cplant(:,1))
   poolcoef1w=(SUM(spread(bgc%ratecp,1,mp_patch)*bgc%cplant,2) -  &
        bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(3)*bgc%cplant(:,3))
   poolcoef1r=(SUM(spread(bgc%ratecp,1,mp_patch)*bgc%cplant,2) -  &
        bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(2)*bgc%cplant(:,2))
   ! Carbon uptake from photosynthesis: 
   canopy%frp = veg%rp20*((3.22-0.046*met%tc)**(0.1*(met%tc-20.))) &
        * poolcoef1 /(365.0*24.0*3600.0)   ! 24/05
   canopy%frpw = veg%rp20*((3.22-0.046*met%tc)**(0.1*(met%tc-20.))) &
        * poolcoef1w /(365.0*24.0*3600.0)    ! 24/05
   canopy%frpr = veg%rp20*((3.22-0.046*met%tc)**(0.1*(met%tc-20.))) &
        * poolcoef1r /(365.0*24.0*3600.0)    ! 24/05
   
   ! This section to be updated as part of carbon module upgrade;
   ! frs is currently calculated in carbon module.
   !canopy%frs  = rsoil(soil%rs20, avgwrs, avgtrs)
   !canopy%frs  = canopy%frs &
   !    * SUM(spread(bgc%ratecs,1, mp_patch) * bgc%csoil,2) &
   !     /(365.0*24.0*3600.0)     !convert 1/year to 1/second
   !WHERE (ssoil%snowd > 1.)
   !   canopy%frs	= canopy%frs / MIN(100.,ssoil%snowd)
   !END WHERE
   canopy%frday = 12.0 * SUM(rdy, 2)
   canopy%fpn = -12.0 * REAL(SUM(an_y, 2),r_1)
   ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
   canopy%dewmm = - REAL(MIN(0.0_r_2,canopy%fevw) + &
        MIN(0.0_r_2,canopy%fevc),r_1) * dels * 1.0e3 / (rhow*air%rlam)
   ! Add dewfall to canopy water storage:
   canopy%cansto = canopy%cansto + canopy%dewmm
   ! Calculate canopy water storage excess:
   canopy%spill=MAX(0.,MIN(0.2*canopy%cansto,MAX(0.0, canopy%cansto-cansat)))
   ! Move excess canopy water to throughfall:
   canopy%through = canopy%through + canopy%spill
   ! Initialise 'throughfall to soil' as 'throughfall from canopy';
   ! snow may absorb
   canopy%precis = canopy%through
   ! Update canopy storage term:
   canopy%cansto=canopy%cansto - canopy%spill
   ! Modify canopy water storage for evaporation:
   canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw,r_1))*dels*1.0e3/ &
        (rhow*air%rlam), 0.0)
   ! Calculate the total change in canopy water store (mm/dels):
   canopy%delwc = canopy%cansto-oldcansto
   ! calculate dgdtg, derivative of ghflux
   ssoil%dfn_dtg = (-1.)*4.*emsoil*sboltz*tss4/ssoil%tss ! d(canopy%fns)
   ! /d(ssoil%tgg)
   ssoil%dfh_dtg = air%rho*capp/ssoil%rtsoil  ! d(canopy%fhs)/d(ssoil%tgg)
   ssoil%dfe_ddq = wetfac*air%rho*air%rlam/ssoil%rtsoil !d(canopy%fes)/d(dq)
   ssoil%ddq_dtg = (rmh2o/rmair)/met%pmb *tetena*tetenb*tetenc &
        / ((tetenc+ssoil%tss-tfrz)**2) &
        * EXP(tetenb*(ssoil%tss-tfrz)/(tetenc+ssoil%tss-tfrz))
   canopy%dgdtg = ssoil%dfn_dtg - ssoil%dfh_dtg &
        - ssoil%dfe_ddq * ssoil%ddq_dtg
   !ypw: energy balance of the dry canopy
   !    bal%drybal=ecy(:,1)+ecy(:,2)+hcy(:,1)+hcy(:,2) &
   !         -rad%rniso(:,1)-rad%rniso(:,2) &
   !         +capp*rmair*(tlfy(:,1)-(met%tvair-tfrz))*rad%gradis(:,1) &
   !         +capp*rmair*(tlfy(:,2)-(met%tvair-tfrz))*rad%gradis(:,2)

   bal%canopy_drybal=REAL(ecy(:,1)+hcy(:,1)-rad%rniso(:,1),r_1) &
        +capp*rmair*(tlfy(:,1)-(met%tvair-tfrz))*REAL(rad%gradis(:,1),r_1)
   bal%canopy_drybal=REAL(ecy(:,2)+hcy(:,2) - rad%rniso(:,2),r_1) &
        +capp*rmair*(tlfy(:,2)-(met%tvair-tfrz))*REAL(rad%gradis(:,2),r_1)
   !ypw: energy balance of the wet canopy
   bal%canopy_wetbal=REAL(canopy%fevw+canopy%fhvw &
        -(rad%rniso(:,1)+rad%rniso(:,2))*REAL(canopy%fwet,r_2) &
        +canopy%fhvw*(rad%gradis(:,1)+rad%gradis(:,2))/ghwet,r_1)

  CONTAINS
    !--------------------------------------------------------------------------
    ELEMENTAL FUNCTION qsatf(tair,pmb) RESULT(r)
      ! MRR, 1987
      ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
      ! HUMIDITY (KG/KG) FROM TETEN FORMULA
      REAL(r_1), INTENT(IN) :: tair ! air temperature (C)
      REAL(r_1), INTENT(IN) :: pmb  ! pressure PMB (mb)
      REAL(r_1)             :: r    ! result; sat sp humidity
      r = (rmh2o/rmair) * (tetena*EXP(tetenb*tair/(tetenc+tair))) / pmb
    END FUNCTION qsatf
    !---------------------------------------------------------
    ELEMENTAL FUNCTION ej3x(parx,x) result(z)
      REAL(r_1), INTENT(IN) :: parx
      REAL(r_1), INTENT(IN) :: x
      REAL(r_1)             :: z
      z = MAX( 0.0 , &
           0.25*((alpha3*parx+x-SQRT((alpha3*parx+x)**2 - &
           4.0*convx3*alpha3*parx*x)) /(2.0*convx3)) )
    END FUNCTION ej3x
    !---------------------------------------------------------
    ELEMENTAL FUNCTION ej4x(parx,x) result(z)
      REAL(r_1), INTENT(IN) :: parx
      REAL(r_1), INTENT(IN) :: x
      REAL(r_1)             :: z
      z = MAX( 0.0 , &
           (alpha4*parx+x-SQRT((alpha4*parx+x)**2 - &
           4.0*convx4*alpha4*parx*x))/(2.0*convx4))
    END FUNCTION ej4x
    !---------------------------------------------------------
    ! Explicit array dimension as temporary work around for NEC inlining problem
    FUNCTION xvcmxt4(x) result(z)
      REAL(r_1), PARAMETER :: q10c4 = 2.0
      REAL(r_1), DIMENSION(mp_patch,mf), INTENT(IN) :: x
      REAL(r_1), DIMENSION(mp_patch,mf)             :: z
      z = q10c4 ** (0.1 * x - 2.5) / &
           ((1.0 + EXP(0.3 * (13.0 - x))) * (1.0 + EXP(0.3 * (x - 36.0))))
    END FUNCTION xvcmxt4
    !---------------------------------------------------------
    FUNCTION xvcmxt3(x) result(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for vcmax for c3 plants
      REAL(r_1), DIMENSION(mp_patch,mf), INTENT(IN) :: x
      REAL(r_1), DIMENSION(mp_patch,mf)             :: xvcnum
      REAL(r_1), DIMENSION(mp_patch,mf)             :: xvcden
      REAL(r_1), DIMENSION(mp_patch,mf)             :: z
      xvcnum = xvccoef * EXP((ehavc/(rgas*trefk))*(1.-trefk/x))
      xvcden = 1.0 + EXP((entropvc*x-ehdvc)/(rgas*x))
      z = MAX(0.0,xvcnum/xvcden)
    END FUNCTION xvcmxt3
    !---------------------------------------------------------
    FUNCTION xejmxt3(x) result(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for jmax for c3 plants
      REAL(r_1), DIMENSION(mp_patch,mf), INTENT(IN) :: x
      REAL(r_1), DIMENSION(mp_patch,mf)             :: xjxnum
      REAL(r_1), DIMENSION(mp_patch,mf)             :: xjxden
      REAL(r_1), DIMENSION(mp_patch,mf)             :: z
      xjxnum = xjxcoef * EXP((ehajx/(rgas*trefk))*(1.-trefk/x))
      xjxden = 1.0 + EXP((entropjx*x-ehdjx)/(rgas*x))
      z = MAX(0.0, xjxnum/xjxden)
    END FUNCTION xejmxt3
    !---------------------------------------------------------
    ELEMENTAL FUNCTION psim(zeta) result(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psim(z/l) (z/l=zeta)
      ! for momentum, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).
      USE cable_math_constants
      REAL(r_1), INTENT(IN) :: zeta
      REAL(r_1)             :: r
      REAL(r_1)             :: x
      REAL(r_1), PARAMETER  :: gu = 16.0
      REAL(r_1), PARAMETER  :: gs = 5.0
      x = (1.0 + gu * ABS(zeta))**0.25
      r = MERGE(LOG((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*ATAN(x) &
           + pi*0.5, -1.0*gs*zeta, zeta < 0.0)
    END FUNCTION psim
    !---------------------------------------------------------
    ELEMENTAL FUNCTION psis(zeta) result(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psis(z/l) (z/l=zeta)
      ! for scalars, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).
      REAL(r_1), INTENT(IN) :: zeta
      REAL(r_1)             :: r
      REAL(r_1), PARAMETER :: gu = 16.0
      REAL(r_1), PARAMETER :: gs = 5.0
      r = MERGE(2.0 * LOG((1.0 + SQRT(1.0 + gu * ABS(zeta))) * 0.5), &
           -1.0 * gs * zeta, zeta < 0.0)
    END FUNCTION psis
    !---------------------------------------------------------
    ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) result(z)
      REAL(r_1), INTENT(IN) :: rpconst
      REAL(r_1), INTENT(IN) :: rpcoef
      REAL(r_1), INTENT(IN) :: tair
      REAL(r_1)             :: z
      z = rpconst * EXP(rpcoef * tair)
    END FUNCTION rplant
    !---------------------------------------------------------
    ELEMENTAL FUNCTION rsoil(rsconst, avgwrs, avgtrs) result(z)
      REAL(r_1), INTENT(IN) :: rsconst
      REAL(r_1), INTENT(IN) :: avgwrs
      REAL(r_1), INTENT(IN) :: avgtrs
      REAL(r_1)             :: z
      z = rsconst * MIN( 1.0 , MAX( 0.0 , MIN( &
           - 0.0178 + 0.2883*avgwrs + 5.0176*avgwrs*avgwrs &
           - 4.5128*avgwrs*avgwrs*avgwrs, &
           0.3320+22.6726*EXP(-5.8184*avgwrs) ) ) ) &
           * MIN( 1.0 , MAX( 0.0 , MIN( 0.0104*(avgtrs**1.3053), &
           5.5956 - 0.1189*avgtrs ) ) )
    END FUNCTION rsoil
    !---------------------------------------------------------
  END SUBROUTINE define_canopy
END MODULE cable_canopy
