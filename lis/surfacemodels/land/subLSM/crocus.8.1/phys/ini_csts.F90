!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE INI_CSTS 
!     ##################
!
!!****  *INI_CSTS * - routine to initialize the module MODD_CST
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the physical constants
!     stored in  module MODD_CST.
!      
!
!!**  METHOD
!!    ------
!!      The physical constants are set to their numerical values 
!!     
!!
!!    EXTERNAL
!!    --------
!!      FMLOOK : to retrieve logical unit number associated to a file
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST     : contains physical constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module MODD_CST, routine INI_CSTS)
!!      
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    18/05/94 
!!      J. Stein    02/01/95  add the volumic mass of liquid water
!!      J.-P. Pinty 13/12/95  add the water vapor pressure over solid ice
!!      J. Stein    29/06/97  add XTH00
!!      V. Masson   05/10/98  add XRHOLI
!!      C. Mari     31/10/00  add NDAYSEC
!!      V. Masson   01/03/03  add XCONDI
!!      A. Voldoire 01/12/09  add XTTSI, XICEC, XTTS for ESM
!!      J. Escobar  28/03/2014 for pb with emissivity/aerosol reset XSURF_TINY=1.0e-80 in real8 case
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!USE MODI_INI_CTURBS
!
!USE MODI_INI_OCEAN_CSTS
!
!USE MODI_INI_SURF_CSTS
!
IMPLICIT NONE
!  
!-------------------------------------------------------------------------------
!
!*       1.     FUNDAMENTAL CONSTANTS
!               ---------------------
!

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('INI_CSTS',0,ZHOOK_HANDLE)

#ifdef SFX_MNH
#ifdef MNH_MPI_DOUBLE_PRECISION
XSURF_TINY    = 1.0e-80
#else
XSURF_TINY    = TINY    (XSURF_TINY    )
#endif
#else
XSURF_TINY    = TINY    (XSURF_TINY    )
#endif
XSURF_TINY_12 = SQRT    (XSURF_TINY    )
XSURF_EPSILON = EPSILON (XSURF_EPSILON ) * 10.0

XPI         = 2.*ASIN(1.)
XKARMAN     = 0.4
XBOLTZ      = 1.380658E-23
XLIGHTSPEED = 299792458.
XPLANCK     = 6.6260755E-34
XAVOGADRO   = 6.0221367E+23
!
!-------------------------------------------------------------------------------
!
!*       2.     ASTRONOMICAL CONSTANTS
!               ----------------------
!
XDAY   = 86400.
XSIYEA = 365.25*XDAY*2.*XPI/ 6.283076
XSIDAY = XDAY/(1.+XDAY/XSIYEA)
XOMEGA = 2.*XPI/XSIDAY
NDAYSEC = 24*3600 ! Number of seconds in a day
!
!-------------------------------------------------------------------------------!
!
!
!*       3.     TERRESTRIAL GEOIDE CONSTANTS
!               ----------------------------
!
XRADIUS = 6371229.
XG      = 9.80665
!
!-------------------------------------------------------------------------------
!
!*       4.     REFERENCE PRESSURE
!               -------------------
!
XP00 = 1.E5
XTH00 = 300.
!-------------------------------------------------------------------------------
!
!*       5.     RADIATION CONSTANTS
!               -------------------
!
!JUAN OVERFLOW XSTEFAN = 2.* XPI**5 * XBOLTZ**4 / (15.* XLIGHTSPEED**2 * XPLANCK**3)
XSTEFAN = ( 2.* XPI**5 / 15. ) * ( (XBOLTZ / XPLANCK)* XBOLTZ ) * (XBOLTZ/(XLIGHTSPEED*XPLANCK))**2 
XI0     = 1370.
!
!-------------------------------------------------------------------------------
!
!*       6.     THERMODYNAMIC CONSTANTS
!               -----------------------
!
XMD    = 28.9644E-3
XMV    = 18.0153E-3
XRD    = XAVOGADRO * XBOLTZ / XMD
XRV    = XAVOGADRO * XBOLTZ / XMV
XCPD   = 7.* XRD /2.
XCPV   = 4.* XRV
XRHOLW = 1000.
XRHOLI = 917.
XCONDI = 2.22
XCL    = 4.218E+3
XCI    = 2.106E+3
XTT    = 273.16
XTTSI  = XTT - 1.8
XICEC  = 0.5
XTTS   = XTT*(1-XICEC) + XTTSI*XICEC
XLVTT  = 2.5008E+6
XLSTT  = 2.8345E+6
XLMTT  = XLSTT - XLVTT
XESTT  = 611.14
XGAMW  = (XCL - XCPV) / XRV
XBETAW = (XLVTT/XRV) + (XGAMW * XTT)
XALPW  = LOG(XESTT) + (XBETAW /XTT) + (XGAMW *LOG(XTT))
XGAMI  = (XCI - XCPV) / XRV
XBETAI = (XLSTT/XRV) + (XGAMI * XTT)
XALPI  = LOG(XESTT) + (XBETAI /XTT) + (XGAMI *LOG(XTT))
!
!-------------------------------------------------------------------------------
!
!*       7.     TURBULENCE CONSTANTS
!               --------------------
!
! CALL INI_CTURBS
!-------------------------------------------------------------------------------
!
!*       8.     OCEAN CONSTANTS
!               ---------------
!
! CALL INI_OCEAN_CSTS
!
!*       9.     SURFACE CONSTANTS
!               -----------------
!
! CALL INI_SURF_CSTS
IF (LHOOK) CALL DR_HOOK('INI_CSTS',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_CSTS 
