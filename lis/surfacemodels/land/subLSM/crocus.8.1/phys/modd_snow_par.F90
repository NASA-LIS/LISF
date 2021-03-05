!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ######################
      MODULE MODD_SNOW_PAR
!     ######################
!
!!****  *MODD_SNOW_PAR* - declaration of parameters related
!!                          to the snow parameterization
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     parameters related to the surface parameterization of snow.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004                    
!! P. Samuelsson  10/2014   MEB complements
!! P. Hagenmuller 07/2014   Mepra/Crocus complements
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!--------------------------------------------------------------------------------
! Snow on the ground: Given in ini_surf_csts and/or in NAM_SURF_CSTS
!--------------------------------------------------------------------------------
!
! Snow emissivity:
!
REAL, SAVE       :: XEMISSN
!
! Minimum and maximum values of the albedo of snow:
!
REAL, SAVE       :: XANSMIN
REAL, SAVE       :: XANSMAX 
!
! Minimum and maximum values of the albedo of permanet snow/ice:
!
REAL, SAVE       :: XAGLAMIN
REAL, SAVE       :: XAGLAMAX
!
! Use recommended settings for snow albedo (FALSE = ISBA default)
! 
LOGICAL,SAVE     :: LMEBREC
!
! Fraction of maximum value of the albedo of snow that is reached for melting
! snow
!
REAL, SAVE       :: XANSFRACMEL
!
! Threeshold temperature above which the snow albedo starts to decrease 
!
REAL, SAVE       :: XTEMPANS
!
! Minimum value of the albedo of snow reached under canopy vegetation:
!
REAL, SAVE       :: XANSMINMEB
! 
! Prescribed ice albedo in 3 spectral bands for glacier simulation with CROCUS scheme.
REAL, SAVE       :: XALBICE1,XALBICE2,XALBICE3
!

! Density threshold for ice detection in CROCUS scheme.
REAL, SAVE       :: XRHOTHRESHOLD_ICE

!for ageing effects
REAL, SAVE      :: XVAGING_NOGLACIER, XVAGING_GLACIER

! percentage of the total pore volume to compute the max liquid water holding capacity
REAL, SAVE      :: XPERCENTAGEPORE

! Height (m) of aged snow in glacier case (allows Pn=1)
!
REAL, SAVE       :: XHGLA
! 
! Coefficient for calculation of snow fraction over vegetation
!
REAL, SAVE       :: XWSNV
!
! Roughness length of pure snow surface (m)
!
REAL, SAVE       :: XZ0SN  
!
! Roughness length for heat of pure snow surface (m)
!
REAL, SAVE       :: XZ0HSN
!
! Roughness length ratio between ice and snow
REAL, SAVE       :: XZ0ICEZ0SNOW
!
! Snow Melt timescale with D95 (s): needed to prevent time step 
! dependence of melt when snow fraction < unity.
!
REAL, SAVE       :: XTAU_SMELT
! Snow impurity deposition rates
REAL,DIMENSION(5), SAVE :: XIMPUR_DRY !(g m-2 s-1)
REAL,DIMENSION(5), SAVE     :: XIMPUR_WET !(g m-2 s-1)

!
!	 Grooming and Snowmaking option by P.Spandre 20160211
REAL,SAVE                        :: XPSR_SNOWMAK
REAL,SAVE                       :: XRHO_SNOWMAK
REAL, SAVE                      :: XTIMESNOWMAK
REAL, SAVE                      :: XPTA_SEUIL
REAL, DIMENSION(5), SAVE        :: XPROD_SCHEME
REAL, DIMENSION(9500), SAVE     :: XPROD_COUNT
REAL, DIMENSION(4), SAVE        :: XSM_END
INTEGER,SAVE                    :: XFREQ_GRO
!
!--------------------------------------------------------------------------------
! Snow on the ground: PARAMETER
!--------------------------------------------------------------------------------
!
! Critical value of the equivalent water content
! of the snow reservoir for snow fractional coverage and albedo computations
!
REAL, PARAMETER       :: XWCRN      = 10.0   ! (kg m-2) Veg (default value)
REAL, PARAMETER       :: XWCRN_EXPL =  1.0   ! (kg m-2) Veg explicit
REAL, PARAMETER       :: XWCRN_ROOF =  1.0   ! (kg m-2)  Roofs 
REAL, PARAMETER       :: XWCRN_ROAD =  1.0   ! (kg m-2)  Roads
REAL, PARAMETER       :: XWCRN_VEG  =  1.0   ! (kg m-2)  Urban veg
!
! Critical value of the total snow depth for ground snow fractional coverage
!
REAL, PARAMETER       :: XDCRN_EXPL = 0.01  ! (m) Veg explicit
!
! Critical value of snow emissivity
!
REAL, PARAMETER       :: XEMCRIN = 0.98
!
! Minimum and maximum values of the albedo of snow:
!
REAL, PARAMETER       :: XANSMIN_ROOF = 0.30 ! (-)   Roofs
REAL, PARAMETER       :: XANSMIN_ROAD = 0.15 ! (-)   Roads
!
REAL, PARAMETER       :: XANSMAX_ROOF = 0.85 ! (-)   Roofs
REAL, PARAMETER       :: XANSMAX_ROAD = 0.85 ! (-)   Roads
!
! Snow aging coefficients (albedo and Force-Restore density):
!
REAL, PARAMETER       :: XANS_TODRY    = 0.008     ! (-) Veg (default value)
REAL, PARAMETER       :: XANS_TODRY_ROOF = 0.008   ! (-)  Roofs
REAL, PARAMETER       :: XANS_TODRY_ROAD = 0.008   ! (-)  Roads
REAL, PARAMETER       :: XANS_TODRY_MEB  = 0.016   ! (-) Surface under canopy vegetation
!
REAL, PARAMETER       :: XANS_T        = 0.240     ! (-) Veg (default value)
REAL, PARAMETER       :: XANS_T_ROOF     = 0.174   ! (-)  Roofs
REAL, PARAMETER       :: XANS_T_ROAD     = 0.174   ! (-)  Roads (alley simul)
REAL, PARAMETER       :: XANS_T_MEB    = 0.480     ! (-) Surface under canopy vegetation
!
! Minimum and maximum values of the density of snow 
! for Force-Restore snow option
!
REAL, PARAMETER       :: XRHOSMIN = 100.       ! (kg m-3)   Veg (Default value)
REAL, PARAMETER       :: XRHOSMIN_ROOF = 100.  ! (kg m-3)   Roofs
REAL, PARAMETER       :: XRHOSMIN_ROAD = 100.  ! (kg m-3)   Roads
!
REAL, PARAMETER       :: XRHOSMAX = 300.       ! (kg m-3)   Veg (Default value)
REAL, PARAMETER       :: XRHOSMAX_ROOF = 300.  ! (kg m-3)   Roofs
REAL, PARAMETER       :: XRHOSMAX_ROAD = 350.  ! (kg m-3)   Roads
!
! Minimum and maximum values of the density of snow 
! for ISBA-ES snow option
!
REAL, PARAMETER       :: XRHOSMIN_ES =  50.  ! (kg m-3)
REAL, PARAMETER       :: XRHOSMAX_ES = 750.  ! (kg m-3)
!
! ISBA-ES Critical snow depth at which snow grid thicknesses constant
!
REAL, PARAMETER       :: XSNOWCRITD = 0.03  ! (m)
!                                       
! ISBA-ES Minimum total snow depth for thermal calculations. 
! Used to prevent numerical problems as snow becomes vanishingly thin. 
!
REAL, PARAMETER      :: XSNOWDMIN = 0.000001  ! (m)
!
!Coefficients for Morin impurities model
! (unitless)

REAL,PARAMETER        :: XIMPUR_EFOLD = 0.005 !(m) e-folding of the exponential decay rate with depth below the surface of the middle of the considered snow layer (0.5*PSNOWDZ(JJ,1)) for the deposition of snow impurities
!  
!Cluzet et al 2016 liquid water content options parameters

! percentage of the total pore volume to compute the max liquid water holding capacity
REAL, PARAMETER      :: XPERCENTAGEPORE_B92 = 0.05 !(%) original parameter value from Crocus, according to Brun et al. 1992 
REAL, PARAMETER      :: XPERCENTAGEPORE_O04 = 0.033!(%) different value used in CLM from Oleson et al. 2004
!
REAL, PARAMETER      :: XWHOLDMAX_S02 = 0.08 !(-)        fixed value for the maximum liquid water mass fracton in SNOWPACK (Lehning et al. 2002)
!                                   
! Maximum Richardson number limit for very stable conditions using the ISBA-ES 'RIL' option
!
REAL, SAVE       :: X_RI_MAX

!                                       
! ISBA-ES Maximum snow liquid water holding capacity (fraction by mass) parameters:
!
REAL, PARAMETER       :: XWSNOWHOLDMAX2   = 0.10  ! (-) 
REAL, PARAMETER       :: XWSNOWHOLDMAX1   = 0.03  ! (-)
REAL, PARAMETER       :: XSNOWRHOHOLD     = 200.0 ! (kg/m3)
!                                       
! ISBA-ES arameters for grain size computation :
!
REAL, PARAMETER       :: XSNOW_AGRAIN = 1.6e-4   ! (m)
REAL, PARAMETER       :: XSNOW_BGRAIN = 1.1e-13  ! (m13/kg4)
REAL, PARAMETER       :: XSNOW_CGRAIN = 0.5e-4   ! (m)
REAL, PARAMETER       :: XDSGRAIN_MAX = 2.796e-3 ! m
!
!--------------------------------------------------------------------------------
! Calibration coefficients for CROCUS and ES albedo computation
!--------------------------------------------------------------------------------
!
REAL, PARAMETER :: XD1 = 1., XD2 = 3., XD3 = 4., XX = 99., &
                   XVALB2 = .96, XVALB3 = 1.58, XVALB4 = .92, XVALB5 = .90, &
                   XVALB6 = 15.4, XVALB7 = 346.3, XVALB8 = 32.31, XVALB9 = .88, &
                   XVALB10 = .200, XVALB11 = .6, XVDIOP1 = 2.3E-3, XVRPRE1 = .5, &
                   XVRPRE2=1.5
!
! for ageing effects:
REAL, PARAMETER :: XVPRES1 = 87000.
!
! spectral bands
!
INTEGER, PARAMETER :: NSPEC_BAND_SNOW = 3
!
! for spectral distribution and thickness effects
REAL, PARAMETER :: XVSPEC1 = .71, XVSPEC2 = .21, XVSPEC3 = .08
!
! for thickness effects
REAL, PARAMETER :: XVW1 = .80, XVW2 = .20 , XVD1 = .02, XVD2 = .01
!
!--------------------------------------------------------------------------------
! calibration coefficients for exctinction computation
REAL, PARAMETER :: XVBETA1 = 1.92E-3, XVBETA2 = 40., XVBETA3 = 1.098E-2, &
                   XVBETA4 = 100.,  XVBETA5 = 2000.
!
! ISBA-ES minimum cosinus of zenithal angle
REAL, PARAMETER :: XMINCOSZEN = 0.01
!
!--------------------------------------------------------------------------------
! ISBA-ES Thermal conductivity coefficients from Anderson (1976):
! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
!
REAL, PARAMETER :: XSNOWTHRMCOND1 = 0.02    ! [W/(m K)]
REAL, PARAMETER :: XSNOWTHRMCOND2 = 2.5E-6  ! [W m5/(kg2 K)]
!
! ISBA-ES Thermal conductivity: Implicit vapor diffn effects
! (sig only for new snow OR high altitudes)
! from Sun et al. (1999): based on data from Jordan (1991)
! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
!
REAL, PARAMETER :: XSNOWTHRMCOND_AVAP = -0.06023 ! [W/(m K)]
REAL, PARAMETER :: XSNOWTHRMCOND_BVAP = -2.5425  ! (W/m)
REAL, PARAMETER :: XSNOWTHRMCOND_CVAP = -289.99  ! (K)
!
! Crocus thermal conducitivity coefficient from Yen (1981)
REAL, PARAMETER :: XVRKZ6 = 1.88
! Cluzet et al 2016, Crocus thermal conductivity coefficients from Calonne et al. 2011
!
REAL, PARAMETER :: XSNOWTHRMCOND_C11_1 = 2.5E-6   ! (W m5 K-1 kg-2)
REAL, PARAMETER :: XSNOWTHRMCOND_C11_2 = -1.23E-4 ! (W m2 K-1 km-1)
REAL, PARAMETER :: XSNOWTHRMCOND_C11_3 = 0.024    ! (W m-1 K-1) 
!--------------------------------------------------------------------------------
! ISBA-ES CROCUS (Pahaut 1976): snowfall density coefficients:
!
REAL, PARAMETER :: XSNOWFALL_A_SN = 109.0  ! kg/m3
REAL, PARAMETER :: XSNOWFALL_B_SN =   6.0  ! kg/(m3 K)
REAL, PARAMETER :: XSNOWFALL_C_SN =  26.0  ! kg/(m7/2 s1/2)
!
! Cluzet et al 2016, Coefficients for A76 fresh snow density option (from Anderson 76)
REAL, PARAMETER       :: XRHOS_A76_1 = 50.  ! (kg m-3)
REAL, PARAMETER       :: XRHOS_A76_2 = 1.7  ! (K-1)
REAL, PARAMETER       :: XRHOS_A76_3 = 15.   ! (K)
! Coefficents for S02 fresh snow density option (from Schmucki and al. 2014)
REAL, PARAMETER       :: XRHOS_S02_1 = 3.28 !
REAL, PARAMETER       :: XRHOS_S02_2 = 0.03 !  
REAL, PARAMETER       :: XRHOS_S02_3 = -0.36 !
REAL, PARAMETER       :: XRHOS_S02_4 = -0.75 ! 
REAL, PARAMETER       :: XRHOS_S02_5 = 0.8  !  (fixed relative humidity value when snowing (RH = 0.8)
REAL, PARAMETER       :: XRHOS_S02_6 = 0.3  !  


! Coefficients for P75 (Pahaut 1975 original law quotesd by Brun 1989 but with different values
REAL, PARAMETER :: XSNOWFALL_A_SN_P75 = 109.0  ! kg/m3
REAL, PARAMETER :: XSNOWFALL_B_SN_P75 =   8.0  ! kg/(m3 K) this one is different from Brun 89
REAL, PARAMETER :: XSNOWFALL_C_SN_P75 =  26.0  ! kg/(m7/2 s1/2)

!
!
! Coefficients for the optimal vertical grid calculation
REAL, PARAMETER :: XDZ1 = 0.01
REAL, PARAMETER :: XDZ2 = 0.0125
REAL, PARAMETER :: XDZ3 = 0.015
REAL, PARAMETER :: XDZ3_BIS = 0.03
REAL, PARAMETER :: XDZ4 = 0.04
REAL, PARAMETER :: XDZ5 = 0.05
REAL, PARAMETER :: XDZ_BASE = 0.02
REAL, PARAMETER :: XDZ_INTERNAL = 0.07
REAL, PARAMETER :: XSCALE_CM = 100.
REAL,DIMENSION(5), PARAMETER :: XDZMAX_INTERNAL = (/0.5,1.,2.,4.,10./)
REAL, PARAMETER :: XDZMIN_TOP_EXTREM = 0.0001
!
! Below this threshold of snowfall, new snowfall are aggregated with surface layer to avoid numerical problems
! (0.03 mm/h)
REAL,PARAMETER :: XSNOWFALL_THRESHOLD = 0.0333/3600.

! The ratio between a new surface layer thickness and the second layer surface thickness is limited to 1/10
REAL,PARAMETER :: XRATIO_NEWLAYER = 0.1

! Coefficients for cases with very thick snowpacks
REAL, PARAMETER :: XDEPTH_THRESHOLD1 = 3.
REAL, PARAMETER :: XDEPTH_THRESHOLD2 = 20.
REAL, PARAMETER :: XDEPTH_SURFACE = 3.
!
! Coefficients for computing the difference in 2 snow layer characteristics
REAL, PARAMETER :: XDIFF_1 = 20.
REAL, PARAMETER :: XDIFF_MAX = 200.
REAL, PARAMETER :: XSCALE_DIFF = 25.
!
! Coeefficients for snow layer splitting
REAL, PARAMETER :: XDZMIN_TOP = 0.01
REAL, PARAMETER :: XDZMIN_TOP_BIS = 0.005
REAL, PARAMETER :: XDZMIN_BOT = 0.02
REAL, PARAMETER :: XSPLIT_COEF = 8.
!
! Coeefficients for snow layer agregation 
REAL, PARAMETER :: XAGREG_COEF_1 = 5.
REAL, PARAMETER :: XAGREG_COEF_2 = 4.5
!
!--------------------------------------------------------------------------------
!
! Calibration coefficients
REAL, PARAMETER :: XVTIME = 48*3600. ! characteristic time for
!compaction and metamorphism by wind drift
!
REAL, PARAMETER :: XVROMAX = 350. !  maximum density for
! drift compaction     UNIT : kg m-3
REAL, PARAMETER :: XVROMIN = 50.  !  minimum density for
! mobility computation UNIT : kg m-3
REAL, PARAMETER :: XVMOB1 = 0.295  !  coefficient for computing
! the mobility index
REAL, PARAMETER :: XVMOB2 = 0.833  !  coefficient for computing
! the mobility index
REAL, PARAMETER :: XVMOB3 = 0.583  !  coefficient for computing
! the mobility index
REAL, PARAMETER :: XVMOB4 = -0.0583 !  coefficient for computing
! the mobility index
REAL, PARAMETER :: XVDRIFT1 = 2.868 !  coefficient for computing
! the drift index
REAL, PARAMETER :: XVDRIFT2 = 0.085 !  coefficient for computing
! the drift index
REAL, PARAMETER :: XVDRIFT3 = 3.25  !  coefficient for computing
! the drift index
REAL, PARAMETER :: XVSIZEMIN = 3.E-4 !  minimum size decrease 
! by drift  UNIT = m
!
! modif_EB pour sublim 
! a pour but de tenir compte du fait que le vent moyen est > rafales
! on en tient compte egalement pour diminuer la duree de l'effet
REAL, PARAMETER :: XCOEF_FF = 1.25 ! coefficient for gust diagnosis from average wind 
REAL, PARAMETER :: XCOEF_EFFECT = 1.0 ! coefficient for impact on density du drift
REAL, PARAMETER :: XQS_REF = 2.E-5 ! valeur de reference de ZQS pour effet neige
!
!--------------------------------------------------------------------------------
!
! ISBA-ES snow grid parameters
!
REAL, PARAMETER, DIMENSION(3)     :: XSGCOEF1  = (/0.25, 0.50, 0.25/) 
REAL, PARAMETER, DIMENSION(2)     :: XSGCOEF2  = (/0.05, 0.34/)       
REAL, PARAMETER, DIMENSION(10)    :: XSGCOEF3  = (/0.025, 0.033, 0.043, &
                                     0.055, 0.071, 0.091, 0.117, 0.150, &
                                     0.193, 0.247/) 
!
! Minimum total snow depth at which surface layer thickness is constant:
!
REAL, PARAMETER                   :: XSNOWTRANS = 0.20                ! (m)
REAL, PARAMETER                   :: XSNOWTRANS1 = 0.40                ! (m)
REAL, PARAMETER                   :: XSNOWTRANS2 = 0.6061                ! (m)
REAL, PARAMETER                   :: XSNOWTRANS3 = 0.7143               ! (m)
REAL, PARAMETER                   :: XSNOWTRANS4 = 0.9259                ! (m)
REAL, PARAMETER                   :: XSNOWTRANS5 = 1.4493                ! (m)
!
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
! Parameters for MEPRA grain type
! Note that mixed types are symetric, i.e. PP_DF is PP/DF or DF/PP. There is
! no distinction between primary and secondary snow types.
! dendritic snow
!
INTEGER, PARAMETER :: JP_PP_PP = 0    ! new snow (PP)                                     fr
INTEGER, PARAMETER :: JP_PP_DF = 1    ! new snow (PP) and decomposed snow (DF)            fr_lb
INTEGER, PARAMETER :: JP_DF_DF = 2    ! decomposed snow (DF)                              lb
INTEGER, PARAMETER :: JP_DF_RG = 3    ! decomposed snow (DF) and rounded grains (RG)      lb_fin
INTEGER, PARAMETER :: JP_DF_FC = 4    ! decomposed snow (DF) and faceted crystals (FC)    lb_ang
! non dendritic snow
INTEGER, PARAMETER :: JP_PP_GP = 5    ! graupel (PPgp) !!!Never used.                     roul
INTEGER, PARAMETER :: JP_RG_RG = 6    ! rounded grains (RG)                               fin
INTEGER, PARAMETER :: JP_RG_MF = 7    ! rounded grains (RG) and melt forms (MF)           fin_ar
INTEGER, PARAMETER :: JP_RG_FC = 8    ! rounded grains (RG) and faceted crystals (FC)     fin_ang
INTEGER, PARAMETER :: JP_FC_FC = 9    ! faceted crystals (FC)                             pl
INTEGER, PARAMETER :: JP_FC_DH = 10   ! faceted crystals (FC) and depth hoar (DH)         pl_gob
INTEGER, PARAMETER :: JP_DH_DH = 11   ! depth hoar (DH)                                   gob
INTEGER, PARAMETER :: JP_MF_MF = 12   ! melt forms (MF)                                   gel
INTEGER, PARAMETER :: JP_MF_DH = 13   ! depth hoar (DH) and melt forms (MF)               gob_fon
INTEGER, PARAMETER :: JP_MF_FC = 14   ! melt forms (MF) and faceted crystals (FC)         ron_ang
!
!
! Correspondance between classes of (dendricity,sphericity) and snow type for
! dendritic snow.
INTEGER, DIMENSION(100), PARAMETER :: JPTAB_DEND = (/&
JP_DF_FC,JP_DF_FC,JP_DF_FC,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP,&
JP_DF_FC,JP_DF_FC,JP_DF_FC,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP,&
JP_DF_FC,JP_DF_FC,JP_DF_FC,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP,&
JP_DF_FC,JP_DF_FC,JP_DF_FC,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP,&
JP_DF_FC,JP_DF_FC,JP_DF_FC,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP,&
JP_DF_RG,JP_DF_RG,JP_DF_RG,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP,&
JP_DF_RG,JP_DF_RG,JP_DF_RG,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP,&
JP_DF_RG,JP_DF_RG,JP_DF_RG,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP,&
JP_DF_RG,JP_DF_RG,JP_DF_RG,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP,&
JP_DF_RG,JP_DF_RG,JP_DF_RG,JP_DF_DF,JP_DF_DF,JP_DF_DF,JP_PP_DF,JP_PP_DF,JP_PP_PP,JP_PP_PP/)
!
!
! Correspondance between classes of (historic,grain size ,sphericity) and snow
! type for non dendritic snow.
INTEGER, DIMENSION(180), PARAMETER :: JPTAB_NODEND = (/&
! Historic = 0, other cases
JP_FC_FC,JP_FC_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_RG,JP_RG_RG,&
JP_FC_FC,JP_FC_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_RG,JP_RG_RG,&
JP_FC_FC,JP_FC_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_RG,JP_RG_RG,&
! Historic = 1, has been angular and was never in contact with liquid water
JP_FC_FC,JP_FC_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_FC,JP_RG_RG,JP_RG_RG,&
JP_FC_DH,JP_FC_DH,JP_FC_DH,JP_FC_DH,JP_FC_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,&
JP_DH_DH,JP_DH_DH,JP_DH_DH,JP_DH_DH,JP_DH_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,&
! Historic = 2, has been in contact with liquid water but was never angular
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_RG_MF,JP_RG_MF,JP_RG_MF,JP_RG_MF,JP_RG_MF,&
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_MF,JP_MF_MF,JP_MF_MF,JP_MF_MF,JP_MF_MF,&
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_MF,JP_MF_MF,JP_MF_MF,JP_MF_MF,JP_MF_MF,&
! Historic = 3, has been in contact with liquid water but has been angular
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,&
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,&
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,&
! Historic = 4, has undergone several melt-freeze cycles and was never angular
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_RG_MF,JP_RG_MF,JP_RG_MF,JP_RG_MF,JP_RG_MF,&
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_MF,JP_MF_MF,JP_MF_MF,JP_MF_MF,JP_MF_MF,&
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_MF,JP_MF_MF,JP_MF_MF,JP_MF_MF,JP_MF_MF,&
! Historic = 5, has undergone several melt-freeze cycles and has been angular
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,&
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,&
JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_FC,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH,JP_MF_DH/)
!
!------------------------------------------------------------------------------
! Parameter of accidental risk calculations
!
! Accidental risk indices
INTEGER*1, PARAMETER :: JPACC_NUL = 0   ! null
INTEGER*1, PARAMETER :: JPACC_LOW = 1   ! low
INTEGER*1, PARAMETER :: JPACC_MOD = 2   ! moderate
INTEGER*1, PARAMETER :: JPACC_HIG = 3   ! high
INTEGER*1, PARAMETER :: JPACC_NAN = 4   ! empty

! Thresholds on strength/stress ratio (-)
REAL, PARAMETER :: XACC_RAT_HIG = 1.5   ! threshold on strength/stress ratio (-)
REAL, PARAMETER :: XACC_RAT_MOD = 2.5   ! threshold on strength/stress ratio (-)

! Threshold on slab shear strength (0.981 kPa)
REAL, PARAMETER :: XACC_SLA_STR = 1.3


!--------------------------------------------------------------------------------------------------
! Parameter of natural risk calculations

! Natural risks indices
INTEGER*1, PARAMETER :: JPNAT_VLO = 0    !very low
INTEGER*1, PARAMETER :: JPNAT_LOW = 1    !low
INTEGER*1, PARAMETER :: JPNAT_MOA = 2    !moderate ascending
INTEGER*1, PARAMETER :: JPNAT_MOD = 3    !moderate descending
INTEGER*1, PARAMETER :: JPNAT_HIG = 4    !high
INTEGER*1, PARAMETER :: JPNAT_VHI = 5    !very high
INTEGER*1, PARAMETER :: JPNAT_NAN = 6    !empty

! Thresholds on strength/stress ratio (-)
REAL, PARAMETER :: XNAT_RAT_HIG  = 2.0
REAL, PARAMETER :: XNAT_RAT_MOD  = 3.0

! Thresholds on release thickness (m)
REAL, PARAMETER :: XNAT_HEI_MIN  = 0.1
REAL, PARAMETER :: XNAT_HEI_HIG  = 0.8!0.2
REAL, PARAMETER :: XNAT_HEI_MOD  = 0.4!0.4
REAL, PARAMETER :: XNAT_HEI_LOW  = 0.2!0.8

! Correspondance between classes of (sup. profil type, strenght/stress ratio, release depth) and
! natural risk indices
INTEGER*1,DIMENSION(45),PARAMETER :: JPNAT_TAB = (/&
!profil sup. NEW
JPNAT_VLO,JPNAT_VLO,JPNAT_VLO,&
JPNAT_VLO,JPNAT_LOW,JPNAT_LOW,&
JPNAT_LOW,JPNAT_MOD,JPNAT_MOA,&
JPNAT_LOW,JPNAT_MOD,JPNAT_HIG,&
JPNAT_LOW,JPNAT_MOD,JPNAT_VHI,&
!profil sup. WET
JPNAT_VLO,JPNAT_VLO,JPNAT_VLO,&
JPNAT_LOW,JPNAT_LOW,JPNAT_LOW,&
JPNAT_LOW,JPNAT_MOD,JPNAT_MOA,&
JPNAT_LOW,JPNAT_MOD,JPNAT_HIG,&
JPNAT_LOW,JPNAT_MOD,JPNAT_VHI,&
!profil sup. FRO
JPNAT_VLO,JPNAT_VLO,JPNAT_VLO,&
JPNAT_VLO,JPNAT_LOW,JPNAT_LOW,&
JPNAT_VLO,JPNAT_LOW,JPNAT_MOD,&
JPNAT_VLO,JPNAT_LOW,JPNAT_MOD,&
JPNAT_VLO,JPNAT_LOW,JPNAT_MOD/)!fake line to be uniform

! Table indicating how to change the current risk as a function of the current and previous
! risk indices for profil sup. NEW
!!!!!!!TO DO what if JPNAT_NAN = NI1UNDEF
! current = column, previous = line
INTEGER*1,DIMENSION(49),PARAMETER :: JPNAT_ACT = (/&
!vlow    ,low      ,moa      ,mod      ,hig      ,vhi      ,nan
JPNAT_VLO,JPNAT_LOW,JPNAT_MOA,JPNAT_MOD,JPNAT_HIG,JPNAT_VHI,JPNAT_NAN,&!vlo
JPNAT_VLO,JPNAT_LOW,JPNAT_LOW,JPNAT_LOW,JPNAT_HIG,JPNAT_VHI,JPNAT_NAN,&!low
JPNAT_VLO,JPNAT_LOW,JPNAT_LOW,JPNAT_LOW,JPNAT_HIG,JPNAT_VHI,JPNAT_NAN,&!moa
JPNAT_VLO,JPNAT_LOW,JPNAT_LOW,JPNAT_LOW,JPNAT_MOD,JPNAT_HIG,JPNAT_NAN,&!mod
JPNAT_VLO,JPNAT_LOW,JPNAT_MOA,JPNAT_MOD,JPNAT_MOD,JPNAT_HIG,JPNAT_NAN,&!hig
JPNAT_VLO,JPNAT_LOW,JPNAT_MOA,JPNAT_MOD,JPNAT_MOD,JPNAT_HIG,JPNAT_NAN,&!vhi
JPNAT_VLO,JPNAT_LOW,JPNAT_MOA,JPNAT_MOD,JPNAT_HIG,JPNAT_VHI,JPNAT_NAN/)!nan

!------------------------------------------------------------------------------
! Parameter of superior/inferior profiles characterization
! (used for natural risk calculation)

! Superior profil types
INTEGER*1, PARAMETER :: JPPRO_SUP_NAN = 6   ! ps_vid in profil.f
INTEGER*1, PARAMETER :: JPPRO_SUP_NEW = 0   ! aggregation of ps_rec_ins, ps_rec_sta, ps_rec_het, ps_rec_pla
INTEGER*1, PARAMETER :: JPPRO_SUP_WET = 4   ! ps_fon_deg
INTEGER*1, PARAMETER :: JPPRO_SUP_FRO = 5   ! ps_fon_gel

! Inferior profil types
INTEGER*1, PARAMETER :: JPPRO_INF_NAN = 6   ! pi_vid (6) in profil.f
INTEGER*1, PARAMETER :: JPPRO_INF_HAR = 0   ! aggregation of pi_hetero (0), pi_fili_pla (3), pi_homo (4), pi_homo_fb (5)
INTEGER*1, PARAMETER :: JPPRO_INF_SOF = 1   ! aggregation of pi_hete_squ (1), pi_fili (2)
!
! Thresholds to define the type of profile
REAL, PARAMETER :: XPRO_SUP_CRU = 0.01      ! threshold on crust height for prosup type NEW (m)
REAL, PARAMETER :: XPRO_SUP_DEP = 0.03      ! threshold on depth for prosup type WET and FRO (m)
REAL, PARAMETER :: XPRO_INF_RAM = 8
REAL, PARAMETER :: XPRO_INF_COE = 0.25
!
!
!--------------------------------------------------------------------------------------------------
! Avalanche types
! (used for natural risk update)
!
INTEGER*1, PARAMETER :: JPAVA_NEW_DRY = 0    ! recente seche
INTEGER*1, PARAMETER :: JPAVA_NEW_WET = 1    ! recente humide
INTEGER*1, PARAMETER :: JPAVA_NEW_MIX = 3    ! recente mixte
INTEGER*1, PARAMETER :: JPAVA_SLA_SUR = 2    ! plaque surface (never used)
INTEGER*1, PARAMETER :: JPAVA_MEL_SUR = 4    ! fonte de surface             (fon_sur)
INTEGER*1, PARAMETER :: JPAVA_MEL_GRO = 5    ! fonte de fond                (fon_fon)
INTEGER*1, PARAMETER :: JPAVA_NAN     = 6    ! vide
!
!
END MODULE MODD_SNOW_PAR
