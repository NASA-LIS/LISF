! cable_canopy_vh.f90
!********************************************************************************
! Edited by VH and MC 16/10/09: Moved incanopy T to top of iteration loop; 
! moved stability calc to top of iteration loop;
! removed restriction on fevc due to extractible water;
! introduced tolerances for qvair and Tvair to reduce number of iterations
! merge wet canopy calc into leaf T iteration loop
! recalc wetbal and drybal
! fwsoil modified to max(alpha2)
! initialised all local variables
!********************************************************************************

! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
! Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 coding by Harvey Davies, Gab Abramowitz and Martin Dix
! bugs to gabsun@gmail.com.
!
! This file contains module canopy_module and subroutine define_canopy.
! The functions included are:
!   qsatf,
!   ej3
!   ej4x,
!   xvcmxt4,
!   xvcmxt3,
!   xejmxt3,
!   psim,
!   psis,
!   rplant, and
!   rsoil
!
MODULE cable_canopy_vh
  USE cable_photosynthetic_constants
  USE cable_radiation
  USE cable_roughness
  USE cable_air
  USE cable_types
  USE cable_physical_constants
  USE cable_dimensions
  IMPLICIT NONE
  PRIVATE
  PUBLIC define_canopy_vh, sinbet
CONTAINS
  
  !--------------------------------------------------------------------------
  SUBROUTINE define_canopy_vh(ktau,bal,rad,rough,air,met,dels,ssoil,soil, &
       veg,bgc,canopy,model_structure)
    TYPE (balances_type),INTENT(INOUT)  :: bal
    TYPE (radiation_type), INTENT(INOUT):: rad
    TYPE (roughness_type), INTENT(INOUT):: rough
    TYPE (air_type), INTENT(INOUT)	:: air
    TYPE (met_type), INTENT(INOUT)	:: met
    REAL(r_1), INTENT(IN)		:: dels ! integration time setp (s)
    TYPE (soil_snow_type), INTENT(INOUT):: ssoil
    TYPE (bgc_pool_type),INTENT(IN)	:: bgc
    TYPE (soil_parameter_type), INTENT(INOUT)	:: soil
    TYPE (veg_parameter_type), INTENT(INOUT)	:: veg
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    TYPE (model_structure_type), INTENT(IN)  :: model_structure
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    REAL(r_1), DIMENSION(mp_patch,mf):: a1c3 ! Spatially varying a1c3
    REAL(r_1) :: a1c3_380 = 4.34 ! fitted value from Hawkwsbury data
    REAL(r_1) :: a1c3_620 = 5.23 ! fitted value from Hawkwsbury data
    REAL(r_1), DIMENSION(mp_patch,mf)	        :: abs_deltlf ! ABS(deltlf)
    REAL(r_2), DIMENSION(mp_patch,mf,3)	:: ancj ! soln to quad eqn
    REAL(r_1), DIMENSION(mp_patch,mf)		:: anx ! net photos. prev iteration
    REAL(r_1), DIMENSION(mp_patch,mf)		:: an_y ! net photosynthesis soln
    !  REAL(r_1), DIMENSION(mp_patch)		:: avgtrs !root weighted mean soil temperature
    !  REAL(r_1), DIMENSION(mp_patch)		:: avgwrs !root weighted mean soil moisture
    REAL(r_1), DIMENSION(mp_patch,mf)		:: ca2	 ! 2D CO2 concentration
    REAL(r_1), DIMENSION(mp_patch)		:: cansat ! max canopy intercept. (mm)
    REAL(r_2), DIMENSION(mp_patch,mf,3)	:: ci ! intercellular CO2 conc.
    REAL(r_1), PARAMETER		:: co2cp3=0.0 ! CO2 compensation pt C3
    REAL(r_2), DIMENSION(mp_patch,mf,3)	:: coef0 ! CO2 comp. pt coeff 1
    REAL(r_2), DIMENSION(mp_patch,mf,3)	:: coef1 ! " 2
    REAL(r_2), DIMENSION(mp_patch,mf,3)	:: coef2 ! " 3
    REAL(r_1), DIMENSION(mp_patch,mf)		:: conkct ! Michaelis Menton const.
    REAL(r_1), DIMENSION(mp_patch,mf)		:: conkot ! Michaelis Menton const.
    REAL(r_1), DIMENSION(mp_patch):: cscale ! scaling between 2 hawkesbury elev co2 vals
    REAL(r_1), DIMENSION(mp_patch,mf)		:: csx ! leaf surface CO2 concentration
    REAL(r_2), DIMENSION(mp_patch,mf,3)	:: cx  ! "d_{3}" in Wang and Leuning, 1998, appendix E
    REAL(r_1), DIMENSION(mp_patch,mf)		:: da2 ! 2D sat vap pres deficit
    REAL(r_1), DIMENSION(mp_patch,mf)		:: dva2 ! 2D in canopy sat vap pres deficit
    REAL(r_2), DIMENSION(mp_patch,mf,3)	:: delcx ! discriminant  in quadratic in eq. E7 Wang and Leuning, 1998
    REAL(r_1), DIMENSION(mp_patch,mf)		:: deltlf ! deltlfy of prev iter.
    REAL(r_1), DIMENSION(mp_patch,mf)		:: deltlfy ! del temp successive iteration
    REAL(r_1), DIMENSION(mp_patch)		:: deltvair ! deltvair 
    REAL(r_1), DIMENSION(mp_patch)		:: delqvair ! delqvair 
    REAL(r_1), DIMENSION(mp_patch)		:: tvair_old ! tvair of prev iter.
    REAL(r_1), DIMENSION(mp_patch)		:: qvair_old ! qvair of prev iteration

    !! variables for method alternative to P-M formula
    !    REAL(r_1), DIMENSION(mp_patch)		:: dq ! sat spec hum diff.
    !    REAL(r_1), DIMENSION(mp_patch)            :: qstss ! sat spec hunidity at soil/snow temperature
    !! end variables for alternative method

    REAL(r_1), DIMENSION(mp_patch,mf)		:: dsatdk2	! 2D dsatdk
    REAL(r_1), DIMENSION(mp_patch,mf)		:: dsx ! leaf surface vpd
    REAL(r_2), DIMENSION(mp_patch,mf)		:: ecx ! lat. hflux big leaf
    REAL(r_1), DIMENSION(mp_patch,mf)		:: ejmax2 ! jmax of big leaf
    REAL(r_1), DIMENSION(mp_patch,mf)		:: ejmxt3 ! jmax big leaf C3 plants
    REAL(r_2), DIMENSION(mp_patch,mf)		:: ecy ! lat heat fl dry big leaf
    REAL(r_1), DIMENSION(mp_patch,mf) :: d0c3 ! Empirical coef for vpd sensitivity of stomata
    REAL(r_1), DIMENSION(mp_patch,mf)		:: frac42	! 2D frac4
    REAL(r_1), DIMENSION(mp_patch)		:: fwsoil ! soil water modifier of stom. cond.
    REAL(r_1), DIMENSION(mp_patch,mf)		:: fwsoil2 ! soil water modifier of stom. cond.
    REAL(r_1), DIMENSION(mp_patch)		:: gaw ! aerodynamic conduct. for water
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gaw2	! 2D gaw
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gbhf ! freeConvectionBndLayerConductance mol/m2/s
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gbhu ! forcedConvectionBoundaryLayerConductance
    REAL(r_1), DIMENSION(mp_patch)		:: gbvtop ! bnd layer cond. top leaf
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gras ! Grashof coeff
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gswmin ! min stomatal conductance
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gswx ! stom cond for water
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gw  ! cond for water for a dry canopy
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gh  ! cond for heat for a dry canopy
    REAL(r_1), DIMENSION(mp_patch,mf)		:: ghr ! dry canopy cond for heat & thermal radiat'n
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gwwet  ! cond for water for a wet canopy
    REAL(r_1), DIMENSION(mp_patch,mf)		:: ghwet  ! cond for heat for a wet canopy
    REAL(r_1), DIMENSION(mp_patch,mf)		:: gbw ! cond for water transfer (leaf boundary layer)
    REAL(r_2), DIMENSION(mp_patch,mf)		:: hcx ! sens heat fl big leaf prev iteration
    REAL(r_2), DIMENSION(mp_patch,mf)		:: hcy ! veg. sens heat
    INTEGER(i_d)			:: iter ! iteration #
    INTEGER(i_d)			:: iterplus !
    INTEGER(i_d)			:: k		! interation count
    INTEGER(i_d)			:: kk		! interation count
    REAL(r_1), DIMENSION(mp_patch,mf)		:: psycst ! modified pych. constant
    REAL(r_1), DIMENSION(mp_patch,mf)		:: rdx ! daytime leaf resp rate, prev iteration
    REAL(r_1), DIMENSION(mp_patch,mf)		:: rdy ! daytime leaf resp rate
    REAL(r_2), DIMENSION(mp_patch,mf)		:: rnx ! net rad prev timestep
    REAL(r_2), DIMENSION(mp_patch,mf)		:: rny ! net rad
    REAL(r_1), DIMENSION(mp_patch)		:: rt0 ! turbulent resistance
    REAL(r_1), DIMENSION(mp_patch)		:: ortsoil ! turbulent resistance, prev time step
    REAL(r_1), DIMENSION(mp_patch)		:: rt1usc ! eq. 3.53, SCAM manual, 1997
    REAL(r_1), DIMENSION(mp_patch)		:: rwater ! soil water availability
    REAL(r_1), DIMENSION(mp_patch,mf)		:: tair2 ! 2D tair
    REAL(r_1), DIMENSION(mp_patch,mf)		:: tvair2 ! 2D tair
    REAL(r_1), DIMENSION(mp_patch,mf)		:: tdiff ! leaf air temp diff.
    REAL(r_1), DIMENSION(mp_patch,mf)		:: tlfx ! leaf temp prev. iteration
    REAL(r_1), DIMENSION(mp_patch,mf)		:: tlfxx ! leaf temperature of current iteration
    REAL(r_1), DIMENSION(mp_patch,mf)		:: tlfy ! leaf temp
    REAL(r_1), DIMENSION(mp_patch,mf)		:: vcmax2 ! vcmax big leaf
    REAL(r_1), DIMENSION(mp_patch,mf)		:: vcmxt3 ! vcmax big leaf C3
    REAL(r_1), DIMENSION(mp_patch,mf)		:: vcmxt4 ! vcmax big leaf C4
    REAL(r_2), DIMENSION(mp_patch,mf,2)	:: vx3 ! carboxylation C3 plants
    REAL(r_2), DIMENSION(mp_patch,mf,2)	:: vx4 ! carboxylation C4 plants
    REAL(r_1), DIMENSION(mp_patch,mf)		:: xdleaf2	! 2D dleaf
    REAL(r_1), DIMENSION(mp_patch,mf)		:: xleuning ! leuning stomatal coeff
    REAL(r_1), DIMENSION(mp_patch,niter)	:: zetar ! stability correction
    REAL(r_1), PARAMETER		:: jtomol = 4.6e-6 ! Conversion from Joule to Mol for light
    REAL(r_1), PARAMETER		:: EHaVc  = 73637.0  !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER		:: EHdVc  = 149252.0 !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER		:: EntropVc = 486.0  !J/mol/K (Leuning 2002)
    REAL(r_1), PARAMETER		:: xVccoef = 1.17461 !derived parameter
    !	xvccoef=1.0+exp((Entropjx*TrefK-EHdjx)/(Rconst*TrefK))
    REAL(r_1), PARAMETER		:: EHaJx  = 50300.0  !J/mol (Leuning 2002)
    REAL(r_1), PARAMETER		:: EHdJx  = 152044.0	!J/mol (Leuning 2002)
    REAL(r_1), PARAMETER		:: EntropJx = 495.0	!J/mol/K (Leuning 2002)
    REAL(r_1), PARAMETER		:: xjxcoef = 1.16715	!derived parameter
    REAL(r_1), PARAMETER		:: effc4 = 4000.0  !Vc=effc4*Ci*Vcmax (see
    ! Bonan,LSM version 1.0, p106)
    REAL(r_1), DIMENSION(mp_patch)		:: oldcansto ! prev t step canopy storage
    REAL(r_1), DIMENSION(mp_patch)		:: cc ! limitation term for canopy interception per timestep		   
    REAL(r_1), DIMENSION(mp_patch)		:: ccfevw ! limitation term for wet canopy evaporation rate  
    REAL(r_1), DIMENSION(mp_patch)		:: denom ! denominator in calculating screen temperature, humidity etc
    REAL(r_1), DIMENSION(mp_patch)		:: tstar ! 
    REAL(r_1), DIMENSION(mp_patch)		:: zscrn !
    REAL(r_1), DIMENSION(mp_patch)		:: qstar !
    REAL(r_1), DIMENSION(mp_patch)		:: rsts  !
    REAL(r_1), DIMENSION(mp_patch)		:: qsurf !
    REAL(r_1), DIMENSION(mp_patch)		:: qtgnet !
    REAL(r_1), DIMENSION(mp_patch)		:: evapfb !
    REAL(r_1), DIMENSION(mp_patch,ms)		:: evapfbl, temp1, temp2, temp3, temp4, temp5 !
    REAL(r_1), DIMENSION(mp_patch)		:: phenps ! Leaf phenology influence on vcmax and jmax
    REAL(r_1), DIMENSION(mp_patch)		:: poolcoef1 ! leaf carbon turnover rate * leaf pool size
    REAL(r_1), DIMENSION(mp_patch)		:: poolcoef1w ! wood carbon turnover rate * wood pool size
    REAL(r_1), DIMENSION(mp_patch)		:: poolcoef1r ! root carbon turnover rate * root pool size
    REAL(r_1), DIMENSION(mp_patch)		:: rbw ! leaf boundary layer resistance for water
    REAL(r_1), DIMENSION(mp_patch)            :: rrbw ! recipr. leaf boundary layer resistance for water
    REAL(r_1), DIMENSION(mp_patch)		:: rsw ! stomatal resistance for water
    REAL(r_1), DIMENSION(mp_patch)            :: rrsw ! recipr. stomatal resistance for water
    REAL(r_1), DIMENSION(mp_patch)		:: dmah ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch)		:: dmbh ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch)		:: dmch ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch)		:: dmae ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch)		:: dmbe ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch)		:: dmce ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp_patch)		:: tss4 ! soil/snow temperature**4
    REAL(r_1), DIMENSION(mp_patch)		:: sss ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp_patch)		:: cc1 ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp_patch)		:: cc2 ! variable for Penman-Monteith evap for soil
    REAL(r_1), DIMENSION(mp_patch)		:: qstvair ! sat spec hunidity at leaf temperature
    REAL(r_1), DIMENSION(mp_patch)		:: xx ! delta-type function for sparse canopy limit, p20 SCAM manual
    REAL(r_1), DIMENSION(mp_patch,mf)		:: temp ! vcmax big leaf C3
    REAL(r_1), DIMENSION(mp_patch,mf)         :: deltecy ! YP & Mao (jun08)
    REAL(r_1), DIMENSION(mp_patch,mf)         :: fwet ! fraction wet canopy
    REAL(r_1), DIMENSION(mp_patch,mf)         :: Ecansto  ! supply limited evap from wet leaves
    REAL(r_1), DIMENSION(mp_patch)  :: dummy1, dummy2, test

    !%% changes by Ashok Luhar (low wind speed)
    REAL(r_1), PARAMETER                :: alpha1=4.0
    REAL(r_1), PARAMETER                :: beta1=0.5
    REAL(r_1), PARAMETER                :: gamma1=0.3
    REAL(r_1), DIMENSION(mp_patch)            :: zeta1
    REAL(r_1), DIMENSION(mp_patch)            :: zeta2
    !%%
    !**************************************************************************************** ! vh 17/07/09     
    REAL(r_1), DIMENSION(mp_patch,ms)         :: alpha1a_root, alpha1b_root, &
         alpha1_root, alpha2a_root,  &
         alpha2_root, alpha_root, delta_root, pwcol
    INTEGER(i_d)			:: n
    !**************************************************************************************** ! vh 17/07/09     

    REAL(r_1), DIMENSION(mp_patch)            :: flpwb, emair  
    INTEGER,   DIMENSION(mp_patch,mf)            :: Flag_fwet   
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

    abs_deltlf = 0
    ancj= 0
    anx= 0
    an_y= 0
    ca2	= 0
    cansat = 0
    ci= 0
    coef0= 0
    coef1= 0
    coef2= 0
    conkct= 0
    conkot= 0
    csx= 0
    cx= 0
    da2= 0
    dva2= 0
    delcx= 0
    deltlf= 0
    deltlfy= 0
    deltvair = 0
    delqvair= 0
    tvair_old = 0
    qvair_old= 0
    dsatdk2	= 0
    dsx= 0
    ecx= 0
    ejmax2= 0
    ejmxt3= 0
    ecy= 0
    frac42	= 0
    fwsoil= 0
    fwsoil2= 0
    gaw= 0
    gaw2 = 0
    gbhf= 0
    gbhu= 0
    gbvtop= 0
    gras= 0
    gswmin= 0
    gswx= 0
    gw= 0
    gh= 0
    ghr= 0
    gwwet = 0
    ghwet= 0
    gbw= 0
    hcx= 0
    hcy= 0
    iter= 0
    iterplus= 0
    k = 0
    kk = 0
    psycst= 0
    rdx= 0
    rdy= 0
    rnx= 0
    rny= 0
    rt0= 0
    ortsoil= 0
    rt1usc= 0
    rwater= 0
    tair2= 0
    tvair2= 0
    tdiff= 0
    tlfx= 0
    tlfxx= 0
    tlfy= 0
    vcmax2= 0
    vcmxt3= 0
    vcmxt4= 0
    vx3= 0
    vx4= 0
    xdleaf2 = 0
    xleuning= 0
    zetar= 0
    oldcansto= 0
    cc= 0
    ccfevw= 0
    denom= 0
    tstar= 0
    zscrn= 0
    qstar= 0
    rsts= 0
    qsurf= 0
    qtgnet= 0
    evapfb= 0
    evapfbl = 0
    temp1 = 0
    temp2 = 0
    temp3 = 0
    temp4 = 0
    temp5= 0
    phenps= 0
    poolcoef1= 0
    poolcoef1w= 0
    poolcoef1r= 0
    rbw= 0
    rrbw= 0
    rsw= 0
    rrsw= 0
    dmah= 0
    dmbh= 0
    dmch= 0
    dmae= 0
    dmbe= 0
    dmce= 0
    tss4= 0
    sss= 0
    cc1= 0
    cc2= 0
    qstvair= 0
    xx= 0
    temp= 0
    deltecy= 0
    fwet = 0
    dummy1 = 0
    dummy2 = 0
    zeta1 = 0
    zeta2 = 0
    alpha1a_root = 0
    alpha1b_root = 0
    alpha1_root = 0
    alpha2_root = 0
    alpha_root = 0
    pwcol = 0
    n = 0
    flpwb = 0
    emair= 0


    !
    !	xjxcoef=1.0+exp((Entropjx*TrefK-EHdjx)/(Rconst*TrefK))
    ! 1-oct-2002 ypw: to keep the unit consistent for resistance or conductance
    ! s/m for r; mol/m2/s for g, and always use g where appropriate
    ! replace rh, rhr, rw  with ghdry/ghwet,ghrdry, ghrwet, gwdry, gwwet

    ! Set surface water vapour pressure deficit:
    met%da = (qsatf(met%tk-tfrz,met%pmb) - met%qv )*rmair/rmh2o*met%pmb*100.0
    ! Soil water limitation on stomatal conductance:
    rwater = MAX(1.0e-4, &
         SUM(veg%froot * MIN(1.0,REAL(ssoil%wb,r_1) - &
         SPREAD(soil%swilt, 2, ms)),2) / (soil%sfc-soil%swilt))
    ! construct function to limit stom cond for soil water availability
    !fwsoil = MAX(1.0e-4,MIN(1.0, veg%vbeta * rwater))

    ! test Lai and Katul formulation for root efficiency function  vh 17/07/09
    !alpha1a_root = MIN(1.0,REAL(ssoil%wb,r_1)/(SPREAD(soil%ssat, 2, ms)-SPREAD(soil%swilt, 2, ms)))
    alpha1a_root = MIN(1.0,REAL(ssoil%wb,r_1)/(SPREAD(soil%ssat, 2, ms)-soil%swilt_vec))
    FORALL (k=1:ms)
       pwcol(:,k) =  SUM(ssoil%wb(:,1:k)*SPREAD(soil%zse (1:k),1,mp_patch),2)                ! partial water column
    END FORALL

    alpha1b_root = REAL(pwcol/SPREAD(pwcol(:,ms),2,ms))

    alpha1_root = MAX(REAL(alpha1a_root,r_1), REAL(alpha1b_root,r_1))

    !alpha2_root =EXP( (SPREAD(veg%gamma,2,ms)/(REAL(ssoil%wb,r_1) - SPREAD(soil%swilt, 2, ms))) *log(alpha1a_root))

    alpha2a_root = MAX((REAL(ssoil%wb,r_1)-soil%swilt_vec),0.001)/(SPREAD(soil%ssat, 2, ms))
    alpha2_root = EXP( (SPREAD(veg%gamma,2,ms)/MAX((REAL(ssoil%wb,r_1) -soil%swilt_vec),1.0e-10)) *log(alpha2a_root))
    WHERE ((ssoil%wb-soil%swilt_vec).gt.0.001)
       alpha2_root =EXP( (SPREAD(veg%gamma,2,ms)/((REAL(ssoil%wb,r_1) -soil%swilt_vec))) *log(alpha2a_root))
    ELSEWHERE
       alpha2_root = 0.0
    ENDWHERE

    WHERE (veg%froot>0.0)
       delta_root = 1.0
    ELSEWHERE
       delta_root = 0.0
    ENDWHERE

    alpha_root = MIN(1.0, alpha1_root*alpha2_root)

    fwsoil = maxval(alpha2_root*delta_root,2)
    fwsoil = MAX(0.0,fwsoil)
    fwsoil2 = SPREAD(fwsoil, 2, mf)             

    ! BATS-type canopy saturation proportional to LAI:
    cansat = veg%canst1 * canopy%vlaiw
    ! Leaf phenology influence on vcmax and jmax
    ! rml 22/10/07 only apply to deciduous types
    WHERE (veg%deciduous)
       phenps = max (1.0e-4, MIN(1.,1. - ( (veg%tmaxvj - ssoil%tgg(:,4)+tfrz)/ &
            (veg%tmaxvj - veg%tminvj) )**2 ) )
       WHERE ( ssoil%tgg(:,4) < (veg%tminvj + tfrz) ) phenps = 0.0
       WHERE ( ssoil%tgg(:,4) > (veg%tmaxvj + tfrz) ) phenps = 1.0
    ELSEWHERE
       phenps = 1.0
    END WHERE
    ! Set previous time step canopy water storage:
    oldcansto=canopy%cansto
    ! to avoid excessive direct canopy evaporation, rainfall rate is limited,
    ! hence canopy interception is limited (EK nov2007, snow scheme)
    ! modified further by timestep requirement to avoid canopy temperature
    ! oscillations (EAK aug08)
    cc =MIN(met%precip-met%precip_s, 4. * MIN(dels,1800.) / 60. /1440. )
    ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
    canopy%wcint = MERGE(MIN(MAX(cansat - REAL(canopy%cansto,r_1),0.0), cc), 0.0, &
         & cc > 0.0  )  ! EK nov2007, snow scheme
    !         cc > 0.0  .AND. met%tk > tfrz)
    ! Define canopy throughfall (100% of precip if temp < 0C, see above):
    canopy%through = met%precip_s + MIN( met%precip - met%precip_s , &
         & MAX(0.0, met%precip - met%precip_s - canopy%wcint) )  ! EK nov2007
    ! Delete line below in case it affects snow sites (precip_s) (EK Jul08)
    !    canopy%through = MIN(met%precip,MAX(0.0, met%precip - canopy%wcint))
    ! Add canopy interception to canopy storage term:
    canopy%cansto = canopy%cansto + canopy%wcint
    ssoil%wetfac = MAX(0.0, MIN(1.0, &
         (REAL(ssoil%wb(:,1),r_1) - soil%swilt) / (soil%sfc - soil%swilt)))
    ! owetfac introduced to reduce sharp changes in dry regions,
    ! especially in offline runs where there may be discrepancies between
    ! timing of precip and temperature change (EAK apr2009)
    !ssoil%wetfac = 0.5*(ssoil%wetfac + ssoil%owetfac)
    ! Temporay fixer for accounting the reduction of soil evap due to freezing
    WHERE ( ssoil%wbice(:,1) > 0. ) ! Prevents divide by zero at glaciated
       ! points where wb and wbice=0.
       ssoil%wetfac = ssoil%wetfac &
            * (1.0 - REAL(ssoil%wbice(:,1)/ssoil%wb(:,1),r_1))**2
    END WHERE
    zetar(:,1) = zeta0 ! stability correction terms
    zetar(:,2) = zetpos + 1 
    xdleaf2 = SPREAD(veg%dleaf, 2, mf) ! characteristic leaf length
    dsatdk2 = SPREAD(air%dsatdk, 2, mf)! deriv of sat vap pressure wrt temp
    ca2 = SPREAD(met%ca, 2, mf)        ! CO2 concentration
    csx = ca2                     ! initialise leaf surface CO2 concentration
    da2 = SPREAD(met%da, 2, mf)   ! water vapour pressure deficit
    dsx = da2                     ! init. leaf surface vpd
    tair2 = SPREAD(met%tvair-tfrz, 2, mf) ! air temp in C
    ejmax2 = SPREAD(veg%ejmax*phenps, 2,mf) !max. pot. electr transp. rate top leaf(mol/m2s)
    vcmax2 = SPREAD(veg%vcmax*phenps, 2,mf) !max. RuBP carboxylsation rate top leaf(mol/m2s)
    tlfx = tair2  ! initialise leaf temp iteration memory variable
    tlfy = tair2  ! initialise current leaf temp
    frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants
    ! weight min stomatal conductance by C3 an C4 plant fractions
    rdy = 0.0       ! init daytime leaf respiration rate
    rdx = 0.0       ! init daytime leaf respiration rate
    an_y = 0.0      ! init current estimate net photos.
    gswx = 0.0     ! default stomatal conuctance 
    gbhf = 1e-3     ! default free convection boundary layer conductance
    gbhu = 1e-3     ! default forced convection boundary layer conductance
    anx = 0.0       ! init net photos. iteration memory variable
    ancj = 0.0    
    !    psycst = SPREAD(air%psyc, 2, mf) ! modified psyc. constant
    ! add on by ypw 1-oct-2002
    gw = 1.0e-3 ! default values of conductance
    gh = 1.0e-3
    ghr= 1.0e-3
    gwwet = 1.0e-3
    ghwet = 1.0e-3
    ! Initialise in-canopy temperatures and humidity:
    met%tvair = met%tk
    met%tvrad = met%tk
    met%qvair = met%qv
    ortsoil = ssoil%rtsoil
    tss4 = ssoil%tss**4	
    deltvair = 999.0
    delqvair = 999.0
    tvair_old = met%tvair
    qvair_old = met%qvair	
    ! Calculate fraction of canopy which is wet:
    canopy%fwet = MAX(0.0,MIN(1.0,0.8*canopy%cansto/MAX(cansat,0.01)))
    fwet = SPREAD(canopy%fwet,2,mf)
    iter = 0		
    Ecansto = SPREAD(canopy%fwet,2,mf)*rhow*SPREAD(air%rlam,2,mf)/(dels*1.0e3)*rad%fvlai/SPREAD(canopy%vlaiw,2,mf) ! supply limited evap from wet leaves
  !  CALL define_air(met, air)
    CALL radiation(ssoil, veg, air, met, rad, canopy)
    
    DO WHILE ((ANY(ABS(deltvair) > 0.1).OR.ANY(ABS(delqvair) > 0.005))  .AND.  iter < niter) ! iterate over in-canopy T and q
       iter = iter+1
       CALL define_air(met, air)
       dsatdk2 = SPREAD(air%dsatdk, 2, mf)

      
       ! monin-obukhov stability parameter zetar=zref/l
       if (iter.eq.1.and.ktau.eq.1) then
          zetar(:,iter) = zeta0
       else
          zetar(:,iter) = -(vonk*grav*rough%zref*(canopy%fh+0.07*canopy%fe))/ &
               max( (air%rho*capp*met%tk*canopy%us**3), 1.e-12)
       endif

       !	constrain zeta to zetpos and zetneg (set in param0)
       zetar(:,iter) = min(zetpos,zetar(:,iter))	 ! zetar too +
       zetar(:,iter) = max(zetneg,zetar(:,iter))	 ! zetar too -

       gswmin = rad%scalex * (gsw03 * (1. - frac42) + gsw04 * frac42)
       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       canopy%us = MAX(1.e-6, &
            vonk * MAX(met%ua,umin) / ( &
            LOG(rough%zref / rough%z0m) - &
            psim(zetar(:,iter)) + &
            psim(zetar(:,iter) * rough%z0m / rough%zref) ))
       !%%change by Ashok Luhar - low wind formulation
       where (zetar(:,iter) > 0.7)
          zeta1=zetar(:,iter) * rough%z0m / rough%zref
          canopy%us = MAX(1.e-6, &
               vonk * MAX(met%ua,umin) / ( &
               alpha1* ((zetar(:,iter)**beta1*  &
               (1.0+gamma1*zetar(:,iter)**(1.0-beta1)))  &
               - (zeta1**beta1*(1.0+gamma1*zeta1**(1.0-beta1))))))
       endwhere
       !%%
       ! Turbulent aerodynamic resistance from roughness sublayer depth to reference height,
       ! x=1 if zref+disp>zruffs, 0 otherwise: thus rt1usc = 0 if zref+disp<zruffs
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
       xx = 0.5 + sign(0.5,rough%zref+rough%disp-rough%zruffs)
       !              correction  by Ian Harman to the 2nd psis term
       rt1usc = xx * (LOG(rough%zref/MAX(rough%zruffs-rough%disp, rough%z0soilsn)) &
            - psis( zetar(:,iter) ) &
            + psis( zetar(:,iter)*(MAX(rough%zruffs-rough%disp,rough%z0soilsn))/rough%zref ) &
            )/vonk

       ! rt0 = turbulent resistance from soil to canopy:
       rt0 = rough%rt0us / canopy%us
       ! Aerodynamic resistance (sum 3 height integrals)/us
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
       rough%rt1 = max(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
       WHERE (ssoil%snowd > 0.1)
          ssoil%wetfac = 1.0
       END WHERE
       ! change canopy%vlaiw requirement to 0.01 for conformity (BP may 2009)
       WHERE (canopy%vlaiw > 0.01)
          ssoil%rtsoil = rt0
       ELSEWHERE
          ssoil%rtsoil = rt0 + rough%rt1
       END WHERE
       ssoil%rtsoil = max(25.,ssoil%rtsoil)   
       WHERE ( ssoil%rtsoil .GT. 2.* ortsoil .OR. ssoil%rtsoil .LT. 0.5*ortsoil )
          ssoil%rtsoil = MAX(25.,0.5*(ssoil%rtsoil + ortsoil))
       END WHERE

       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! prandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       gbvtop = air%cmolar*apol * air%visc / prandt / veg%dleaf *	&
            (canopy%us / MIN(rough%usuh, 0.2) * &
            veg%dleaf / air%visc)**0.5 * prandt**(1.0/3.0) / veg%shelrb
       ! Forced convection boundary layer conductance (see Wang & Leuning 1998, AFM):
       !                                gbhu corrected by F.Cruz & A.Pitman on 13/03/07
       gbhu(:,1) = gbvtop*(1.0-EXP(-canopy%vlaiw*(0.5*rough%coexp+rad%extkb))) / &
            (rad%extkb+0.5*rough%coexp)
       gbhu(:,2) = (2.0/rough%coexp)*gbvtop*  &
            (1.0-EXP(-0.5*rough%coexp*canopy%vlaiw))-gbhu(:,1)
       ! Aerodynamic conductance:
       gaw = air%cmolar / rough%rt1
       WHERE (veg%meth > 0 )
          gaw=100000.0
       END WHERE
       gaw2 = SPREAD(gaw, 2, mf)

       WHERE (veg%meth > 0 .and. canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn) 
          !      use the dispersion matrix (DM) to find the air temperature and specific humidity 
          !      (Raupach, Finkele and Zhang 1997, pp 17)
          ! leaf boundary layer resistance for water
          rbw = air%cmolar/sum(gbhu+gbhf,2)    ! gbhf initially set to 1e-3
          rrbw = sum(gbhu+gbhf,2)/air%cmolar  ! MJT
          ! leaf stomatal resistance for water
          rsw = air%cmolar/sum(gswx,2)
          rrsw = sum(gswx,2)/air%cmolar ! MJT
          ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmah = (rt0+rough%rt1)*((1.+air%epsi)*rrsw +rrbw) &
               + air%epsi * (rt0*rough%rt1)*(rrbw*rrsw)
          ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmbh = (-air%rlam/capp)*(rt0*rough%rt1)*(rrbw*rrsw)
          ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmch = ((1.+air%epsi)*rrsw +rrbw)*rt0*rough%rt1* &
               (canopy%fhv + canopy%fhs)/(air%rho*capp)
          dmch = ((1.+air%epsi)*rrsw +rrbw)*rt0*rough%rt1* &
               (canopy%fhv )/(air%rho*capp)


          ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmae = (-air%epsi*capp/air%rlam)*(rt0*rough%rt1)*(rrbw*rrsw)
          ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmbe = (rt0+ssoil%wetfac*rough%rt1)*((1.+air%epsi)*rrsw +rrbw)+(rt0*rough%rt1)*(rrbw*rrsw)
          ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
          dmce = ((1.+air%epsi)*rrsw +rrbw)*rt0*rough%rt1*(canopy%fev + canopy%fes)/ &
               (air%rho*air%rlam)

          dmce = ((1.+air%epsi)*rrsw +rrbw)*rt0*rough%rt1*(canopy%fev )/ &
               (air%rho*air%rlam)
          tvair_old = met%tvair
          qvair_old = met%qvair
          ! Within canopy air temperature:
          where (veg%vlai.gt.0.1)
             met%tvair = met%tk  + (dmbe*dmch-dmbh*dmce)/(dmah*dmbe-dmae*dmbh+1.0e-12)
          elsewhere
             met%tvair = met%tk 
          endwhere

         ! where (abs(met%tvair-met%tk).gt.5.0)
         !    met%tvair = met%tk
         ! endwhere

          ! Within canopy specific humidity:
          where (veg%vlai.gt.0.1)
             met%qvair = met%qv + (dmah*dmce-dmae*dmch)/(dmah*dmbe-dmae*dmbh+1.0e-12)
          elsewhere
             met%qvair = met%qv
          endwhere
          met%qvair = max(0.0,met%qvair)

       END WHERE
       ! Saturated specific humidity in canopy:
       qstvair = qsatf((met%tvair-tfrz),met%pmb)
       met%qvair = min(qstvair,met%qvair)          ! avoid -ve dva
       ! Saturated vapour pressure deficit in canopy:
       met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.
       ! 2 Dim saturated vapour pressure deficit in canopy:
       dva2 = SPREAD(met%dva, 2, mf)
       ! 2 dim Within canopy air temperature in degrees C:
       tvair2 = SPREAD(met%tvair-tfrz, 2, mf)  ! N.B. tvair2 and tair2 ned to be equal because 
	                                           !the temperature difference (Ta-Tl) need to be 
											   !the same for both the sensible heat flux and non-isothermal radiation flux at the leaf 
	   tair2 = tvair2
       ! Set radiative temperature as within canopy air temp:
       met%tvrad = met%tvair
	   ! call radiation here: longwave isothermal radiation absorption
	   ! and radiation conductance depends on tvrad
	   CALL radiation(ssoil, veg, air, met, rad, canopy)
       ! Store change in canopy air temperature and humidity between successive iterations
       deltvair = met%tvair - tvair_old
       delqvair = met%qvair - qvair_old

       CALL define_air(met, air)
       dsatdk2 = SPREAD(air%dsatdk, 2, mf)

       hcx = 0.0       ! init sens heat iteration memory variable
       ecx = rad%rniso ! init lat heat iteration memory variable
       rnx = rad%rniso ! init net rad iteration memory variable
       rny = rad%rniso ! init current estimate net rad
       hcy = 0.0       ! init current estimate lat heat
       ecy = rny - hcy ! init current estimate lat heat
       abs_deltlf = 999.0
       deltlfy = 999.0
       ! Initialise, over each gridpoint, sunlit and shaded leaves:
       DO k=1,mp_patch
          DO kk=1,mf
             IF(rad%fvlai(k,kk) <=1.0e-7) THEN
                abs_deltlf(k,kk)=0.0
                hcx(k,kk) = 0.0 ! intialise
                ecx(k,kk) = 0.0 ! intialise
                anx(k,kk) = 0.0 ! intialise
                rnx(k,kk) = 0.0 ! intialise
                rny(k,kk) = rnx(k,kk) ! store initial values
                hcy(k,kk) = hcx(k,kk) ! store initial values
                ecy(k,kk) = ecx(k,kk) ! store initial values
                rdy(k,kk) = rdx(k,kk) ! store initial values
                an_y(k,kk) = anx(k,kk) ! store initial values
             END IF
          ENDDO
       ENDDO
       deltlfy = abs_deltlf
       k = 0
       Flag_fwet = 0

       DO WHILE ((ANY(abs_deltlf > 0.1) .OR. ANY(Flag_fwet.eq.1))  .AND.  k < maxiter)
          k = k + 1
          ! Where vegetation and no convergence...
          WHERE (rad%fvlai > 1e-7 .and. abs_deltlf > 0.1 .or. Flag_fwet.eq.1)
             ! Grashof number (Leuning et al, 1995) eq E4:
             gras = max(1.0e-6,1.595E8*ABS(tlfx-tair2)*(xdleaf2**3.))
             ! See Appendix E in (Leuning et al, 1995):
             gbhf = rad%fvlai*SPREAD(air%cmolar, 2, mf)*0.5*dheat *(gras**0.25)/xdleaf2
             ! Conductance for heat:
             !             YP removing double counting of gh and gw (Feb 2009)
             !             gh = 1.0/(MIN(1e3, SPREAD(1.0/gaw, 2, mf) + 0.5/(gbhu+gbhf)))
             gh = 2.0 * (gbhu + gbhf)
             ghwet=2.0*(gbhu + gbhf) ! NB changed wet leaf bl conductance to include forced component VH 16/10/09
             gh = fwet*ghwet + (1.-fwet)*gh
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
             ! for calculating anx, csx and gswx for Rubisco limited, RuBP limited,
             ! sink limited
             vx3(:,:,1) = vcmxt3
             vx4(:,:,1) = vcmxt4
             temp = rad%qcan(:,:,1)*jtomol*(1.0-frac42)
             vx3(:,:,2) = ej3x(temp,ejmxt3)
             temp = frac42*rad%qcan(:,:,1)*jtomol
             vx4(:,:,2) = ej4x(temp,vcmxt4)
             rdx = (cfrd3*vcmxt3+cfrd4*vcmxt4)*fwsoil2
             xleuning = (1.0-frac42)*a1c3/(1.0+dsx/d0c3) +frac42*a1c4/(1.0+dsx/d0c4)
             xleuning = xleuning * fwsoil2 / (csx-co2cp3)
             ! Rubisco limited:
             coef2(:,:,1) = gswmin*fwsoil2/rgswc+xleuning *(vx3(:,:,1)-(rdx-vx4(:,:,1)))
             coef1(:,:,1) = (1.0-csx*xleuning) *(vx3(:,:,1)+vx4(:,:,1)-rdx)	&
                  +(gswmin*fwsoil2/rgswc)*(cx(:,:,1)-csx) -xleuning*(vx3(:,:,1)*cx(:,:,2)/2.0 &
                  +cx(:,:,1)*(rdx-vx4(:,:,1)))
             coef0(:,:,1) = -(1.0-csx*xleuning) *(vx3(:,:,1)*cx(:,:,2)/2.0  &
                  +cx(:,:,1)*(rdx-vx4(:,:,1))) -(gswmin*fwsoil2/rgswc)*cx(:,:,1)*csx
             ! Discriminant in quadratic in eq. E7 Wang and Leuning, 1998
             delcx(:,:,1) = coef1(:,:,1)**2 -4.0*coef0(:,:,1)*coef2(:,:,1)
             ci(:,:,1) = (-coef1(:,:,1)+SQRT(MAX(0.0_r_2,delcx(:,:,1)))) /(2.0*coef2(:,:,1))
             ci(:,:,1) = MAX(0.0_r_2,ci(:,:,1))
             ancj(:,:,1) = vx3(:,:,1)*(ci(:,:,1)-cx(:,:,2)/2.0)	&
                  / (ci(:,:,1) + cx(:,:,1)) + vx4(:,:,1) - rdx
             ! RuBP limited:
             coef2(:,:,2) = gswmin*fwsoil2/rgswc+xleuning *(vx3(:,:,2)-(rdx-vx4(:,:,2)))
             coef1(:,:,2) = (1.0-csx*xleuning) *(vx3(:,:,2)+vx4(:,:,2)-rdx)	&
                  +(gswmin*fwsoil2/rgswc)*(cx(:,:,2)-csx) -xleuning*(vx3(:,:,2)*cx(:,:,2)/2.0 &
                  +cx(:,:,2)*(rdx-vx4(:,:,2)))
             coef0(:,:,2) = -(1.0-csx*xleuning) *(vx3(:,:,2)*cx(:,:,2)/2.0  &
                  +cx(:,:,2)*(rdx-vx4(:,:,2))) -(gswmin*fwsoil2/rgswc)*cx(:,:,2)*csx
             delcx(:,:,2) = coef1(:,:,2)**2 -4.0*coef0(:,:,2)*coef2(:,:,2)
             ci(:,:,2) = (-coef1(:,:,2)+SQRT(MAX(0.0_r_2,delcx(:,:,2)))) /(2.0*coef2(:,:,2))
             ci(:,:,2) = MAX(0.0_r_2,ci(:,:,2))
             ancj(:,:,2) = vx3(:,:,2)*(ci(:,:,2)-cx(:,:,2)/2.0)	&
                  /(ci(:,:,2)+cx(:,:,2)) +vx4(:,:,2)-rdx
             ! Sink limited:
             coef2(:,:,3) = xleuning
             coef1(:,:,3) = gswmin*fwsoil2/rgswc + xleuning * (rdx - 0.5*vcmxt3)  +  &
                  effc4 * vcmxt4 - xleuning * csx * effc4 *vcmxt4
             coef0(:,:,3) = -(gswmin*fwsoil2/rgswc)*csx *effc4*vcmxt4 +	&
                  (rdx -0.5*vcmxt3)*gswmin*fwsoil2/rgswc
             delcx(:,:,3) = coef1(:,:,3)**2 -4.0*coef0(:,:,3)*coef2(:,:,3)
             ancj(:,:,3)  = (-coef1(:,:,3)+SQRT(MAX(0.0_r_2,delcx(:,:,3)))) &
                  /(2.0*coef2(:,:,3))
             anx = REAL(MIN(ancj(:,:,1),ancj(:,:,2),ancj(:,:,3)),r_1)
             !             YP removing double counting of gh and gw (Feb 2009)
             !             csx = ca2 - anx * (1.0/gaw2+rgbwc/(gbhu + gbhf))
             csx = ca2 - anx * (gbhu + gbhf) / rgbwc
             gswx = gswmin*fwsoil2+MAX(0.0,rgswc*xleuning *anx)
             ! Recalculate conductance for water:
             gbw  = 1.075*(gbhu+gbhf)
             gwwet= 1.075*(gbhu+gbhf)  ! NB changed wet leaf bl conductance to include forced component VH 16/10/09
             !             YP removing double counting of gh and gw (Feb 2009)
             !             gw = 1.0/(1.0/gswx + 1.0/(1.075*(gbhu+gbhf)) + SPREAD(1.0/gaw, 2, mf))
			 where (gswx>1.e-15)
             gw = 1.0/(1.0/gswx + 1.0/(gbw))
			 elsewhere
			 gw = 0.0
			 endwhere
             gw = (1.-fwet)*gw + fwet*gbw  ! corrected vapour conductance for influence of wet part of leaf. Still need to reduce carbon conducatnce by factor of (1-fwet)
            ! gw = 1.0/(1.0/gswx + 1.0/(1.075*(gbhu+gbhf)))
             ! Modified psychrometric constant (Monteith and Unsworth, 1990)
             psycst = SPREAD(air%psyc, 2, mf) *ghr/gw
             ! Store leaf temperature:
             tlfxx = tlfx
             ! Update canopy latent heat flux:
             ecx = (dsatdk2*rad%rniso +capp*rmair*da2*ghr) /(dsatdk2+psycst)
             ! Update canopy sensible heat flux:
             hcx = (rad%rniso-ecx)*gh/ghr
             ! Update leaf temperature:
             tlfx=tair2+REAL(hcx,r_1)/(capp*rmair*gh)
             ! Update net radiation for canopy:
             rnx = rad%rniso-capp*rmair*(tlfx -tair2)*rad%gradis
             ! Update leaf surface vapour pressure deficit:
             ! dsx = ecx*100.0* SPREAD(met%pmb, 2, mf) /(gswx*rmh2o*SPREAD(air%rlam, 2, mf))
             dsx = da2 + dsatdk2 * (tlfx-tair2)
             ! Store change in leaf temperature between successive iterations:
             deltlf = tlfxx-tlfx
             abs_deltlf = ABS(deltlf)
          END WHERE
          ! Where leaf temp change b/w iterations is significant, and difference is 
          ! smaller than the previous iteration, store results:
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
             !        after four iteration, take the mean value of current and previous estimates
             !        as the next estimate of leaf temperature, to avoid oscillation
             tlfx = (0.5*(MAX(0,k-5)/(k-4.9999))) *tlfxx + &
                  (1.0- (0.5*(MAX(0,k-5)/(k-4.9999))))*tlfx
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
       END DO  ! DO WHILE (ANY(abs_deltlf > 0.1)	.AND.  k < maxiter)

       canopy%fev = sum(ecy,2)
       canopy%fhv = sum(hcy,2)
       canopy%fnv = sum(rny,2)
       an_y = an_y*(1.-fwet)  ! only allow photosynthesis on dry part of leaf


       IF(ANY(model_structure%soil=='soilsnow')) THEN ! evaulate canopy%fes for use in dispersion matrix calc (only if Soil-Snow is being used)
          ! Penman-Monteith formula
          sss=air%dsatdk
          cc1=sss/(sss+air%psyc )
          cc2=air%psyc /(sss+air%psyc )

		  canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*emleaf* &
          sboltz*met%tvair**4 - emsoil*sboltz* tss4

          ssoil%potev = cc1 * (canopy%fns - canopy%ga) + cc2 * air%rho  &
               & * air%rlam * (qsatf((met%tk-tfrz),met%pmb) - met%qv) / ssoil%rtsoil 

          ! Soil latent heat:
          canopy%fes= ssoil%wetfac * ssoil%potev
          WHERE (ssoil%snowd < 0.1 .AND. canopy%fes > 0.0)
             ! Reduce for wilting point limitation:
             canopy%fes= MIN( canopy%fes, MAX(0.0, &
                  (REAL(ssoil%wb(:,1),r_1)-soil%swilt)) *soil%zse(1)*1000.0*air%rlam/dels)
             ! Reduce for soil ice limitation:
             canopy%fes = MIN(canopy%fes,REAL(ssoil%wb(:,1)-ssoil%wbice(:,1),r_1) &
                  * soil%zse(1) * 1000. * air%rlam / dels)
          END WHERE
          ssoil%cls=1.
          WHERE (ssoil%snowd >= 0.1)
             ssoil%cls = 1.1335
             canopy%fes= MIN(ssoil%wetfac * ssoil%potev,ssoil%snowd/dels*air%rlam*ssoil%cls)
          END WHERE

          ! Calculate soil sensible heat:
          canopy%fhs = air%rho*capp*(ssoil%tss - met%tk) /ssoil%rtsoil
          ! Calculate ground heat flux:
          canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
       END IF

       ! Calculate total latent heat:
       canopy%fe = canopy%fev + canopy%fes
       ! Calculate total sensible heat:
       canopy%fh = canopy%fhv + canopy%fhs


       where (fwet*ecy>Ecansto)  ! move this inside Tleaf loop?
          fwet = Ecansto/ecy
          Flag_fwet = 1
       elsewhere (ecy.lt.0.0_r_2) ! dew formation
		fwet = 1.0
		Flag_fwet = 1
	   elsewhere
          Flag_fwet = 0
       endwhere


    END DO	     ! do iter = 1, niter


	where (sum(ecy,2).gt.0.0_r_2)  ! evaporation
!		canopy%fevw = min(canopy%fwet*REAL(sum(ecy,2),r_1),&
!					max(0.0,canopy%cansto)*rhow*air%rlam/(dels*1.0e3))
!		canopy%fevw = min(sum(REAL(ecy,r_1)*fwet,2),&
!					max(0.0,canopy%cansto)*rhow*air%rlam/(dels*1.0e3))
		canopy%fevw = min(sum(REAL(ecy,r_1)*fwet*gbw/max(gw,1.e-12),2),&
					max(0.0,canopy%cansto)*rhow*air%rlam/(dels*1.0e3))
		
		canopy%fwet = REAL(canopy%fevw/sum(ecy,2),r_1)
	elsewhere ! condensation
		canopy%fevw = sum(ecy*fwet*gbw/max(gw,1.e-12),2)
		canopy%fwet = 1.0
	endwhere

    where (sum(gswx,2)>0.0)
		canopy%fevc = canopy%fev - canopy%fevw 
	elsewhere
		canopy%fevc = 0.0
	endwhere
	if (any(canopy%fevc<-0.01)) then
		write (*,*) 'negative trans ', ktau, canopy%fev, canopy%fevc
	endif 
	where (canopy%fevc<0.0.and.canopy%fevc>-1.e-5)  ! negative values of fevc due to precision
		canopy%fevc = 0.0
		canopy%fevw =canopy%fev
	endwhere
    canopy%fhvw = canopy%fhv*canopy%fwet 

    rad%lwabv = (capp*rmair*(tlfy(:,1) - &
         tvair2(:,1))*rad%gradis(:,1) &
         +capp*rmair*(tlfy(:,2) - tvair2(:,2))*rad%gradis(:,2))      ! non-isothermal emitted long-wave radiation

	
	test = -rad%lwabv + sum(rad%rniso,2)-canopy%fnv

    WHERE (canopy%vlaiw > 0.01 .and. rough%hruff > rough%z0soilsn)
       canopy%tv = (rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf)+met%tvrad**4)**0.25
     ELSEWHERE ! sparse canopy
       canopy%tv = met%tvair
    END WHERE

    ! Calculate net radiation absorbed by soil:
    canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd)*emleaf* &
         sboltz*canopy%tv**4 - emsoil*sboltz* tss4


	 ! Calculate radiative/skin temperature:
    rad%trad = ( (1.-rad%transd)*canopy%tv**4 &
         & + rad%transd * ssoil%tss**4 )**0.25

    if (1.eq.0) then              ! comment out calc of variables at screen height vh 15/07/09
       ! screen temp., windspeed and relative humidity at 1.8m
       tstar = - canopy%fh / ( air%rho*capp*canopy%us)
       qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
       zscrn = max(rough%z0m,1.8-rough%disp)
       !    denom = ( log(rough%zref/zscrn)- psim(zetar(:,iterplus)) + &
       !         psim(zetar(:,iterplus) * zscrn / rough%zref) ) /vonk
       denom = ( log(rough%zref/zscrn)- psis(zetar(:,iterplus)) + &
            psis(zetar(:,iterplus) * zscrn / rough%zref) ) /vonk

       !%% change by Ashok Luhar
       where (zetar(:,iterplus) > 0.7)
          zeta2=zetar(:,iterplus) * zscrn / rough%zref
          denom =alpha1* ((zetar(:,iterplus)**beta1* &
               (1.0+gamma1*zetar(:,iterplus)**(1.0-beta1)))  &
               - (zeta2*beta1*(1.0+gamma1*zeta2**(1.0-beta1)))) /vonk
       endwhere
       !%%

       ! Calculate screen temperature:
       canopy%tscrn = met%tk-tfrz - tstar * denom
       rsts = qsatf(canopy%tscrn, met%pmb)
       qtgnet = rsts * ssoil%wetfac - met%qv
       canopy%cduv = canopy%us * canopy%us / (MAX(met%ua,umin))**2 ! EK jun08
       !    canopy%cduv = canopy%us * canopy%us / max(met%ua,umin)
       WHERE (qtgnet > 0.0)
          qsurf = rsts * ssoil%wetfac
       ELSEWHERE
          qsurf = 0.1*rsts*ssoil%wetfac + 0.9*met%qv
       END WHERE
       canopy%qscrn = qsurf + qstar * denom
       canopy%uscrn = max(0.0, max(met%ua,umin) - canopy%us * denom )	 ! at present incorrect
    endif              ! comment out calc of variables at screen height vh 15/07/09
    !  avgwrs = REAL(SUM(veg%froot * ssoil%wb,2),r_1)
    !  avgtrs = max(0.0,sum(veg%froot * ssoil%tgg,2)-tfrz)
    poolcoef1=(sum(spread(bgc%ratecp,1,mp_patch)*bgc%cplant,2) - &
         bgc%ratecp(1)*bgc%cplant(:,1))
    poolcoef1w=(sum(spread(bgc%ratecp,1,mp_patch)*bgc%cplant,2) -  &
         bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(3)*bgc%cplant(:,3))
    poolcoef1r=(sum(spread(bgc%ratecp,1,mp_patch)*bgc%cplant,2) -  &
         bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(2)*bgc%cplant(:,2))
    ! YP & Mao (jun08) capped the tc values
    !    WHERE (met%tk-tfrz>=70.)
    !       met_tc=69.  ! use a temporary variable below to replace met%tk-tfrz
    !    END WHERE
    ! Carbon uptake from photosynthesis: 
    canopy%frp = veg%rp20*((3.22-0.046*(met%tk-tfrz)) &
         **(0.1*(met%tk-tfrz-20.0))) * poolcoef1 /(365.0*24.0*3600.0)
    canopy%frpw = veg%rp20*((3.22-0.046*(met%tk-tfrz)) &
         **(0.1*(met%tk-tfrz-20.0))) * poolcoef1w /(365.0*24.0*3600.0)
    canopy%frpr = veg%rp20*((3.22-0.046*(met%tk-tfrz)) &
         **(0.1*(met%tk-tfrz-20.0))) * poolcoef1r /(365.0*24.0*3600.0)

    ! This section to be updated as part of carbon module upgrade;
    ! frs is currently calculated in carbon module.
    !canopy%frs  = rsoil(soil%rs20, avgwrs, avgtrs)
    !canopy%frs  = canopy%frs &
    !     * sum(spread(bgc%ratecs,1, mp_patch) * bgc%csoil,2)	&
    !     /(365.0*24.0*3600.0)		     !convert 1/year to 1/second
    !WHERE (ssoil%snowd > 1.)
    !   canopy%frs	= canopy%frs / min(100.,ssoil%snowd)
    !END WHERE
    canopy%frday = 12.0 * sum(rdy, 2)
    canopy%fpn = -12.0 * sum(an_y, 2)
    ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
    canopy%dewmm = - REAL(min(0.0_r_2,canopy%fevw) + min(0.0_r_2,canopy%fevc),r_1) * &
         dels * 1.0e3 / (rhow*air%rlam)
    ! Add dewfall to canopy water storage:
    canopy%cansto = canopy%cansto + canopy%dewmm
    ! Calculate canopy water storage excess:
    canopy%spill=max(0.,min(0.2*canopy%cansto,max(0.0, canopy%cansto-cansat)))
    ! Move excess canopy water to throughfall:
    canopy%through = canopy%through + canopy%spill
    ! Initialise 'throughfall to soil' as 'throughfall from canopy'; snow may absorb
    canopy%precis = canopy%through
    ! Update canopy storage term:
    canopy%cansto=canopy%cansto - canopy%spill
    ! Modify canopy water storage for evaporation:
    canopy%cansto = max(canopy%cansto-max(0.0,REAL(canopy%fevw,r_1))*dels*1.0e3/ &
         (rhow*air%rlam), 0.0)
    ! Calculate the total change in canopy water store (mm/dels):
    canopy%delwc = canopy%cansto-oldcansto
    ! calculate dgdtg, derivative of ga
    ssoil%dfn_dtg = (-1.)*4.*emsoil*sboltz*tss4/ssoil%tss	     ! d(canopy%fns)/d(ssoil%tgg)
    ssoil%dfh_dtg = air%rho*capp/ssoil%rtsoil	   ! d(canopy%fhs)/d(ssoil%tgg)
    ssoil%dfe_ddq = ssoil%wetfac*air%rho*air%rlam/ssoil%rtsoil	! d(canopy%fes)/d(dq)
    ssoil%ddq_dtg = (rmh2o/rmair)/met%pmb *tetena*tetenb*tetenc &
         /((tetenc+ssoil%tss-tfrz)**2)*exp(tetenb*(ssoil%tss-tfrz)/(tetenc+ssoil%tss-tfrz))
    !MC Copy cls from soil_snow part above (2 lines)
    ssoil%cls=1.
    WHERE (ssoil%snowd >= 0.1) ssoil%cls = 1.1335
    canopy%dgdtg = ssoil%dfn_dtg - ssoil%dfh_dtg &
         - ssoil%cls * ssoil%dfe_ddq * ssoil%ddq_dtg

    bal%canopy_wetbal = canopy%fwet*(canopy%fev+canopy%fhv-canopy%fnv)
    bal%canopy_drybal = (1.-canopy%fwet)*(canopy%fev+canopy%fhv-canopy%fnv)
    
    ! owetfac will need to be outside define_canopy 
    ! because the UM driver will call define_canopy twice
    ssoil%owetfac = ssoil%wetfac

    !******************************************************************************************
    canopy%gw = gw   ! edit vh 6/7/09
    canopy%ancj = ancj*(-12.0) ! edit vh 7/7/09
    canopy%gswx = gswx ! edit vh 7/7/09
    canopy%tlfy = tlfy ! edit vh 7/7/09
    canopy%ecy =ecy ! edit vh 7/7/09
    canopy%ecx =ecx ! edit vh 7/7/09
    canopy%ci = ci ! edit vh 7/7/09
    canopy%fwsoil = fwsoil ! edit vh 7/7/09
    !******************************************************************************************

  CONTAINS
    !--------------------------------------------------------------------------
    ELEMENTAL FUNCTION qsatf(tair,pmb) RESULT(r)
      ! MRR, 1987
      ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
      ! HUMIDITY (KG/KG) FROM TETEN FORMULA
      REAL(r_1), INTENT(IN) :: tair ! air temperature (C)
      REAL(r_1), INTENT(IN) :: pmb  ! pressure PMB (mb)
      REAL(r_1)		  :: r    ! result; sat sp humidity
      r = (rmh2o/rmair) * (tetena*EXP(tetenb*tair/(tetenc+tair))) / pmb
    END FUNCTION qsatf
    !---------------------------------------------------------
    ELEMENTAL FUNCTION ej3x(parx,x) result(z)
      REAL(r_1), INTENT(IN)	:: parx
      REAL(r_1), INTENT(IN)	:: x
      REAL(r_1)			:: z
      z = max(0.0, &
           0.25*((alpha3*parx+x-sqrt((alpha3*parx+x)**2 - &
           4.0*convx3*alpha3*parx*x)) /(2.0*convx3)) )
    END FUNCTION ej3x
    !---------------------------------------------------------
    ELEMENTAL FUNCTION ej4x(parx,x) result(z)
      REAL(r_1), INTENT(IN)	:: parx
      REAL(r_1), INTENT(IN)	:: x
      REAL(r_1)			:: z
      z = max(0.0, &
           (alpha4*parx+x-sqrt((alpha4*parx+x)**2 - &
           4.0*convx4*alpha4*parx*x))/(2.0*convx4))
    END FUNCTION ej4x
    !---------------------------------------------------------
    ! Explicit array dimensions as temporary work around for NEC inlining problem
    FUNCTION xvcmxt4(x) result(z)
      REAL(r_1), PARAMETER	:: q10c4 = 2.0
      REAL(r_1), DIMENSION(mp_patch,mf), INTENT(IN)	:: x
      REAL(r_1), DIMENSION(mp_patch,mf)			:: z
      z = q10c4 ** (0.1 * x - 2.5) / &
           ((1.0 + exp(0.3 * (13.0 - x))) * (1.0 + exp(0.3 * (x - 36.0))))
    END FUNCTION xvcmxt4
    !---------------------------------------------------------
    FUNCTION xvcmxt3(x) result(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for vcmax for c3 plants
      REAL(r_1), DIMENSION(mp_patch,mf), INTENT(IN)	:: x
      REAL(r_1), DIMENSION(mp_patch,mf)		:: xvcnum
      REAL(r_1), DIMENSION(mp_patch,mf)		:: xvcden
      REAL(r_1), DIMENSION(mp_patch,mf)		:: z
      xvcnum=xvccoef*exp((ehavc/(rgas*trefk))*(1.-trefk/x))
      xvcden=1.0+exp((entropvc*x-ehdvc)/(rgas*x))
      z = max(0.0,xvcnum/xvcden)
    END FUNCTION xvcmxt3
    !---------------------------------------------------------
    FUNCTION xejmxt3(x) result(z)
      !  leuning 2002 (p c & e) equation for temperature response
      !  used for jmax for c3 plants
      REAL(r_1), DIMENSION(mp_patch,mf), INTENT(IN)	:: x
      REAL(r_1), DIMENSION(mp_patch,mf)		:: xjxnum
      REAL(r_1), DIMENSION(mp_patch,mf)		:: xjxden
      REAL(r_1), DIMENSION(mp_patch,mf)		:: z
      xjxnum=xjxcoef*exp((ehajx/(rgas*trefk))*(1.-trefk/x))
      xjxden=1.0+exp((entropjx*x-ehdjx)/(rgas*x))
      z = max(0.0, xjxnum/xjxden)
    END FUNCTION xejmxt3
    !---------------------------------------------------------
    ELEMENTAL FUNCTION psim(zeta) result(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psim(z/l) (z/l=zeta)
      ! for momentum, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).
      USE cable_math_constants
      REAL(r_1), INTENT(IN)	:: zeta
      REAL(r_1)			:: r
      REAL(r_1)			:: x
      REAL(r_1), PARAMETER	:: gu = 16.0
      REAL(r_1), PARAMETER	:: gs = 5.0
      REAL(r_1)                 :: z
      REAL(r_1)                 :: stable
      REAL(r_1)                 :: unstable
      !      x = (1.0 + gu*abs(zeta))**0.25
      !      r = merge(log((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*atan(x) &
      !           + pi*0.5, -gs*zeta, zeta < 0.0)
      z        = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable
      stable   = -gs*zeta
      x        = (1.0 + gu*abs(zeta))**0.25
      unstable = log((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*atan(x) + pi*0.5
      r        = z*stable + (1.0-z)*unstable
    END FUNCTION psim
    !---------------------------------------------------------
    ELEMENTAL FUNCTION psis(zeta) result(r)
      ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
      ! computes integrated stability function psis(z/l) (z/l=zeta)
      ! for scalars, using the businger-dyer form for unstable cases
      ! and the webb form for stable cases. see paulson (1970).
      REAL(r_1), INTENT(IN)	:: zeta
      REAL(r_1)			:: r
      REAL(r_1), PARAMETER	:: gu = 16.0
      REAL(r_1), PARAMETER	:: gs = 5.0
      REAL(r_1)                 :: z
      REAL(r_1)                 :: y
      REAL(r_1)                 :: stable
      REAL(r_1)                 :: unstable
      !      r = merge(2.0 * log((1.0 + sqrt(1.0 + gu * abs(zeta))) * 0.5), &
      !           - gs * zeta, zeta < 0.0)
      z        = 0.5 + sign(0.5,zeta)    ! z=1 in stable, 0 in unstable 
      stable   = -gs*zeta
      y        = (1.0 + gu*abs(zeta))**0.5
      unstable = 2.0 * log((1+y)*0.5)
      r        = z*stable + (1.0-z)*unstable
    END FUNCTION psis
    !---------------------------------------------------------
    ELEMENTAL FUNCTION rplant(rpconst, rpcoef, tair) result(z)
      REAL(r_1), INTENT(IN)	:: rpconst
      REAL(r_1), INTENT(IN)	:: rpcoef
      REAL(r_1), INTENT(IN)	:: tair
      REAL(r_1)			:: z
      z = rpconst * exp(rpcoef * tair)
    END FUNCTION rplant
    !---------------------------------------------------------
    ELEMENTAL FUNCTION rsoil(rsconst, avgwrs, avgtrs) result(z)
      REAL(r_1), INTENT(IN)	:: rsconst
      REAL(r_1), INTENT(IN)	:: avgwrs
      REAL(r_1), INTENT(IN)	:: avgtrs
      REAL(r_1)			:: z
      z = rsconst * min(1.0, max(0.0, min(&
           -0.0178+0.2883*avgwrs+5.0176*avgwrs*avgwrs-4.5128*avgwrs*avgwrs*avgwrs, &
           0.3320+22.6726*exp(-5.8184*avgwrs)))) &
           * min(1.0, max(0.0, min( 0.0104*(avgtrs**1.3053), 5.5956-0.1189*avgtrs)))
    END FUNCTION rsoil
    !---------------------------------------------------------
  END SUBROUTINE define_canopy_vh
END MODULE cable_canopy_vh

