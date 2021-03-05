!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     ##################
      MODULE MODE_SNOW3L
!     ##################
!
!!****  *MODE_SNOW * - contains explicit snow (ISBA-ES) characteristics functions
!!                     for total liquid water holding capacity of a snow layer (m)
!!                     and the thermal heat capacity of a layer (J K-1 m-3)
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!    direct calculation
!!
!!    EXTERNAL
!!    --------
!!       
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!    Boone and Etchevers, J. HydroMeteor., 2001
!!      
!!
!!    AUTHOR
!!    ------
!!	A. Boone       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        01/08/02
!!      V. Masson       01/2004  add snow grid computations
!!      V. Vionnet      06/2008 -Introduction of Crocus formulation to
!                       calculate maximum liquid water holding capacity
!!                              - New algorithm to compute snow grid :
!                       10 layers
!!                              - Routine to aggregate snow grain type
!                       from 2 layers    
!!                              _ Routine to compute average grain
!                       type when snow depth< 3 cm. 
!     S. Morin          02/2011 - Add routines for Crocus
!     A. Boone          02/2012 - Add optimization of do-loops.
!     C. Carmagnola     12/2012 - Add the case CSNOWMETAMO!='B92' in subroutine SNOW3LAVGRAIN and in function SNOW3LDIFTYP
!     M. Lafaysse       01/2013 - Remove SNOWCROWLIQMAX routines (not used)
!     M. Lafaysse       08/2013 - simplification of routine SNOW3LAVGRAIN (logical GDENDRITIC)
!     B. Decharme       07/2013 - SNOW3LGRID cleanning 
!                                 New algorithm to compute snow grid for 6-L or 12-L
!     A. Boone          10/2014 - Added snow thermal conductivity routines
!     B. Decharme       01/2015 - Added optical snow grain size diameter
!     B. Cluzet         08/2015 - deleted unused procedures SNOWCROHOLD_3,2,1D
!                               - added lwc options (functions SNOWO04HOLD_0D, SNOWS02HOLD_0D) 
!----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
INTERFACE SNOW3LWLIQMAX
  MODULE PROCEDURE SNOW3LWLIQMAX_3D
  MODULE PROCEDURE SNOW3LWLIQMAX_2D
  MODULE PROCEDURE SNOW3LWLIQMAX_1D
END INTERFACE
!
INTERFACE SNOW3LHOLD
  MODULE PROCEDURE SNOW3LHOLD_3D
  MODULE PROCEDURE SNOW3LHOLD_2D
  MODULE PROCEDURE SNOW3LHOLD_1D
  MODULE PROCEDURE SNOW3LHOLD_0D
END INTERFACE
!
INTERFACE SNOWCROHOLD
  MODULE PROCEDURE SNOWCROHOLD_3D
  MODULE PROCEDURE SNOWCROHOLD_2D
  MODULE PROCEDURE SNOWCROHOLD_1D
  MODULE PROCEDURE SNOWCROHOLD_0D
END INTERFACE
! Cluzet et al 2016
INTERFACE SNOWSPKHOLD
  MODULE PROCEDURE SNOWSPKHOLD_0D
END INTERFACE
!
INTERFACE SNOWO04HOLD
  MODULE PROCEDURE SNOWO04HOLD_0D
END INTERFACE
! fin Cluzet et al 2016
INTERFACE SNOW3LSCAP
  MODULE PROCEDURE SNOW3LSCAP_3D
  MODULE PROCEDURE SNOW3LSCAP_2D
  MODULE PROCEDURE SNOW3LSCAP_1D
  MODULE PROCEDURE SNOW3LSCAP_0D
END INTERFACE
!
INTERFACE SNOW3L_MARBOUTY
  MODULE PROCEDURE SNOW3L_MARBOUTY
END INTERFACE
!
INTERFACE SNOW3LGRID
  MODULE PROCEDURE SNOW3LGRID_2D
  MODULE PROCEDURE SNOW3LGRID_1D
END INTERFACE
!
INTERFACE SNOW3LAGREG
  MODULE PROCEDURE SNOW3LAGREG
END INTERFACE
!
INTERFACE SNOW3LAVGRAIN
  MODULE PROCEDURE SNOW3LAVGRAIN
END INTERFACE
!
INTERFACE SNOW3LDIFTYP
  MODULE PROCEDURE SNOW3LDIFTYP
END INTERFACE
!
INTERFACE GET_MASS_HEAT
  MODULE PROCEDURE GET_MASS_HEAT
END INTERFACE
!
INTERFACE GET_DIAM
  MODULE PROCEDURE GET_DIAM
END INTERFACE
!
INTERFACE SNOW3LRADABS
  MODULE PROCEDURE SNOW3LRADABS_2D
  MODULE PROCEDURE SNOW3LRADABS_1D
  MODULE PROCEDURE SNOW3LRADABS_0D
END INTERFACE
!
INTERFACE SNOW3LRADABS_SFC
  MODULE PROCEDURE SNOW3LRADABS_SFC
END INTERFACE 
!
INTERFACE SNOW3LTHRM
  MODULE PROCEDURE SNOW3LTHRM
END INTERFACE
!
INTERFACE SNOW3LDOPT
  MODULE PROCEDURE SNOW3LDOPT_2D
  MODULE PROCEDURE SNOW3LDOPT_1D
  MODULE PROCEDURE SNOW3LDOPT_0D
END INTERFACE
!
INTERFACE SNOW3LALB
  MODULE PROCEDURE SNOW3LALB
END INTERFACE
!
INTERFACE SYVAGRE
  MODULE PROCEDURE SYVAGRE
END INTERFACE
!
INTERFACE SNOW3LFALL
  MODULE PROCEDURE SNOW3LFALL
END INTERFACE
!
INTERFACE SNOW3LTRANSF
  MODULE PROCEDURE SNOW3LTRANSF
END INTERFACE
!
INTERFACE SNOW3LCOMPACTN
  MODULE PROCEDURE SNOW3LCOMPACTN
END INTERFACE
!-------------------------------------------------------------------------------
CONTAINS
!
!####################################################################
      FUNCTION SNOW3LWLIQMAX_3D(PSNOWRHO) RESULT(PWLIQMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES, XSNOWRHOHOLD,      &
                          XWSNOWHOLDMAX2, XWSNOWHOLDMAX1
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN)                                  :: PSNOWRHO ! (kg/m3)
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2),SIZE(PSNOWRHO,3)) :: PWLIQMAX ! (kg/m3)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2),SIZE(PSNOWRHO,3)) :: ZHOLDMAXR, ZSNOWRHO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
! Evaluate capacity using upper density limit:
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LWLIQMAX_3D',0,ZHOOK_HANDLE)
ZSNOWRHO(:,:,:) = MIN(XRHOSMAX_ES, PSNOWRHO(:,:,:))
!
! Maximum ratio of liquid to SWE:
!
ZHOLDMAXR(:,:,:) = XWSNOWHOLDMAX1 + (XWSNOWHOLDMAX2-XWSNOWHOLDMAX1)*    &
                  MAX(0.,XSNOWRHOHOLD-ZSNOWRHO(:,:,:))/XSNOWRHOHOLD 
!
! Maximum liquid water holding capacity of the snow (kg/m3):
!
PWLIQMAX(:,:,:) = ZHOLDMAXR(:,:,:)*ZSNOWRHO(:,:,:)           
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LWLIQMAX_3D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LWLIQMAX_3D
!####################################################################
      FUNCTION SNOW3LWLIQMAX_2D(PSNOWRHO) RESULT(PWLIQMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES, XSNOWRHOHOLD,      &
                          XWSNOWHOLDMAX2, XWSNOWHOLDMAX1
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)                   :: PSNOWRHO ! (kg/m3)
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: PWLIQMAX ! (kg/m3)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZHOLDMAXR, ZSNOWRHO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
! Evaluate capacity using upper density limit:
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LWLIQMAX_2D',0,ZHOOK_HANDLE)
ZSNOWRHO(:,:) = MIN(XRHOSMAX_ES, PSNOWRHO(:,:))
!
! Maximum ratio of liquid to SWE:
!
ZHOLDMAXR(:,:) = XWSNOWHOLDMAX1 + (XWSNOWHOLDMAX2-XWSNOWHOLDMAX1)*    &
                  MAX(0.,XSNOWRHOHOLD-ZSNOWRHO(:,:))/XSNOWRHOHOLD 
!
! Maximum liquid water holding capacity of the snow (kg/m3):
!
PWLIQMAX(:,:) = ZHOLDMAXR(:,:)*ZSNOWRHO(:,:)               
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LWLIQMAX_2D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LWLIQMAX_2D
!####################################################################
      FUNCTION SNOW3LWLIQMAX_1D(PSNOWRHO) RESULT(PWLIQMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES, XSNOWRHOHOLD,      &
                           XWSNOWHOLDMAX2, XWSNOWHOLDMAX1
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)  :: PSNOWRHO ! (kg/m3)
!
REAL, DIMENSION(SIZE(PSNOWRHO)) :: PWLIQMAX ! (kg/m3)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO)) :: ZHOLDMAXR, ZSNOWRHO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
! Evaluate capacity using upper density limit:
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LWLIQMAX_1D',0,ZHOOK_HANDLE)
ZSNOWRHO(:) = MIN(XRHOSMAX_ES, PSNOWRHO(:))
!
! Maximum ratio of liquid to SWE:
!
ZHOLDMAXR(:) = XWSNOWHOLDMAX1 + (XWSNOWHOLDMAX2-XWSNOWHOLDMAX1)*    &
                  MAX(0.,XSNOWRHOHOLD-ZSNOWRHO(:))/XSNOWRHOHOLD 
!
! Maximum liquid water holding capacity of the snow (kg/m3):
!
PWLIQMAX(:) = ZHOLDMAXR(:)*ZSNOWRHO(:)                
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LWLIQMAX_1D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LWLIQMAX_1D

!####################################################################
!####################################################################
!####################################################################
      FUNCTION SNOW3LHOLD_3D(PSNOWRHO,PSNOWDZ) RESULT(PWHOLDMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_CSTS,     ONLY : XRHOLW
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES, XSNOWRHOHOLD,      &
                          XWSNOWHOLDMAX2, XWSNOWHOLDMAX1
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN)                 :: PSNOWDZ, PSNOWRHO
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2),SIZE(PSNOWRHO,3)) :: PWHOLDMAX
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2),SIZE(PSNOWRHO,3)) :: ZHOLDMAXR, ZSNOWRHO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------
! Evaluate capacity using upper density limit:
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LHOLD_3D',0,ZHOOK_HANDLE)
ZSNOWRHO(:,:,:) = MIN(XRHOSMAX_ES, PSNOWRHO(:,:,:))
!
! Maximum ratio of liquid to SWE:
!
ZHOLDMAXR(:,:,:) = XWSNOWHOLDMAX1 + (XWSNOWHOLDMAX2-XWSNOWHOLDMAX1)*    &
                  MAX(0.,XSNOWRHOHOLD-ZSNOWRHO(:,:,:))/XSNOWRHOHOLD 
!
! Maximum liquid water holding capacity of the snow (m):
!
PWHOLDMAX(:,:,:) = ZHOLDMAXR(:,:,:)*PSNOWDZ(:,:,:)*ZSNOWRHO(:,:,:)/XRHOLW
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LHOLD_3D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LHOLD_3D
!####################################################################
      FUNCTION SNOW3LHOLD_2D(PSNOWRHO,PSNOWDZ) RESULT(PWHOLDMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_CSTS,     ONLY : XRHOLW
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES, XSNOWRHOHOLD,      &
                          XWSNOWHOLDMAX2, XWSNOWHOLDMAX1
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)                   :: PSNOWDZ, PSNOWRHO
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: PWHOLDMAX
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZHOLDMAXR, ZSNOWRHO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
! Evaluate capacity using upper density limit:
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LHOLD_2D',0,ZHOOK_HANDLE)
ZSNOWRHO(:,:) = MIN(XRHOSMAX_ES, PSNOWRHO(:,:))
!
! Maximum ratio of liquid to SWE:
!
ZHOLDMAXR(:,:) = XWSNOWHOLDMAX1 + (XWSNOWHOLDMAX2-XWSNOWHOLDMAX1)*    &
                  MAX(0.,XSNOWRHOHOLD-ZSNOWRHO(:,:))/XSNOWRHOHOLD 
!
! Maximum liquid water holding capacity of the snow (m):
!
PWHOLDMAX(:,:) = ZHOLDMAXR(:,:)*PSNOWDZ(:,:)*ZSNOWRHO(:,:)/XRHOLW
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LHOLD_2D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LHOLD_2D
!####################################################################
      FUNCTION SNOW3LHOLD_1D(PSNOWRHO,PSNOWDZ) RESULT(PWHOLDMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_CSTS,     ONLY : XRHOLW
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES, XSNOWRHOHOLD,     &
                          XWSNOWHOLDMAX2, XWSNOWHOLDMAX1 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)                     :: PSNOWDZ, PSNOWRHO
!
REAL, DIMENSION(SIZE(PSNOWRHO))                    :: PWHOLDMAX
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO))                    :: ZHOLDMAXR, ZSNOWRHO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------------
! Evaluate capacity using upper density limit:
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LHOLD_1D',0,ZHOOK_HANDLE)
ZSNOWRHO(:) = MIN(XRHOSMAX_ES, PSNOWRHO(:))
!
! Maximum ratio of liquid to SWE:
!
ZHOLDMAXR(:) = XWSNOWHOLDMAX1 + (XWSNOWHOLDMAX2-XWSNOWHOLDMAX1)*     &
                  MAX(0.,XSNOWRHOHOLD-ZSNOWRHO(:))/XSNOWRHOHOLD 
!
! Maximum liquid water holding capacity of the snow (m):
!
PWHOLDMAX(:) = ZHOLDMAXR(:)*PSNOWDZ(:)*ZSNOWRHO(:)/XRHOLW               
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LHOLD_1D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LHOLD_1D
!####################################################################
      FUNCTION SNOW3LHOLD_0D(PSNOWRHO,PSNOWDZ) RESULT(PWHOLDMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_CSTS,     ONLY : XRHOLW
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES, XSNOWRHOHOLD,     &
                          XWSNOWHOLDMAX2, XWSNOWHOLDMAX1
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)        :: PSNOWDZ, PSNOWRHO
!
REAL                    :: PWHOLDMAX
!
!*      0.2    declarations of local variables
!
REAL                    :: ZHOLDMAXR, ZSNOWRHO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
! Evaluate capacity using upper density limit:
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LHOLD_0D',0,ZHOOK_HANDLE)
ZSNOWRHO = MIN(XRHOSMAX_ES, PSNOWRHO)
!
! Maximum ratio of liquid to SWE:
!
ZHOLDMAXR = XWSNOWHOLDMAX1 + (XWSNOWHOLDMAX2-XWSNOWHOLDMAX1) *     &
                             MAX(0.,XSNOWRHOHOLD-ZSNOWRHO)/XSNOWRHOHOLD 
!
! Maximum liquid water holding capacity of the snow (m):
!
PWHOLDMAX = ZHOLDMAXR*PSNOWDZ*ZSNOWRHO/XRHOLW              
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LHOLD_0D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LHOLD_0D

!####################################################################
      FUNCTION SNOWCROHOLD_3D(PSNOWRHO,PSNOWLIQ,PSNOWDZ) RESULT(PWHOLDMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_CSTS,     ONLY : XRHOLW,XRHOLI
USE MODD_SNOW_PAR, ONLY : XPERCENTAGEPORE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN)                 :: PSNOWDZ, PSNOWLIQ, PSNOWRHO
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2),SIZE(PSNOWRHO,3)) :: PWHOLDMAX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
! computation of water holding capacity based on Crocus, 
!taking into account the conversion between wet and dry density - 
!S. Morin/V. Vionnet 2010 12 09

! PWHOLDMAX is expressed in m water for each layer
! In short, PWHOLDMAX = XPERCENTAGEPORE_B92 * porosity * PSNOWDZ .
! The porosity is computed as (rho_ice - (rho_snow - lwc))/(rho_ice)
! where everything has to be in kg m-3 units. In practice, since
! PSNOWLIQ is expressed in m water, expressing the lwc in kg m-3
! is achieved as PSNOWLIQ*XRHOLW/PSNOWDZ. After some rearranging one
! obtains the equation given above.
! Note that equation (19) in Vionnet et al., GMD 2012, is wrong,
! because it does not take into account the fact that liquid water
! content has to be substracted from total density to compute the
! porosity.


IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_3D',0,ZHOOK_HANDLE)
PWHOLDMAX(:,:,:) = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ * (XRHOLI-PSNOWRHO) + PSNOWLIQ*XRHOLW)
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_3D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOWCROHOLD_3D
!####################################################################
      FUNCTION SNOWCROHOLD_2D(PSNOWRHO,PSNOWLIQ,PSNOWDZ) RESULT(PWHOLDMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_CSTS,     ONLY : XRHOLW,XRHOLI
USE MODD_SNOW_PAR, ONLY : XPERCENTAGEPORE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)                   :: PSNOWDZ, PSNOWRHO, PSNOWLIQ
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: PWHOLDMAX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
! computation of water holding capacity based on Crocus, 
!taking into account the conversion between wet and dry density - 
!S. Morin/V. Vionnet 2010 12 09

! PWHOLDMAX is expressed in m water for each layer
! In short, PWHOLDMAX = XPERCENTAGEPORE * porosity * PSNOWDZ .
! The porosity is computed as (rho_ice - (rho_snow - lwc))/(rho_ice)
! where everything has to be in kg m-3 units. In practice, since
! PSNOWLIQ is expressed in m water, expressing the lwc in kg m-3
! is achieved as PSNOWLIQ*XRHOLW/PSNOWDZ. After some rearranging one
! obtains the equation given above.
! Note that equation (19) in Vionnet et al., GMD 2012, is wrong,
! because it does not take into account the fact that liquid water
! content has to be substracted from total density to compute the
! porosity.

IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_2D',0,ZHOOK_HANDLE)
PWHOLDMAX(:,:) = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ * (XRHOLI-PSNOWRHO) + PSNOWLIQ*XRHOLW)
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_2D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOWCROHOLD_2D
!####################################################################
!####################################################################
!####################################################################
      FUNCTION SNOWCROHOLD_1D(PSNOWRHO,PSNOWLIQ,PSNOWDZ) RESULT(PWHOLDMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_CSTS,     ONLY : XRHOLW,XRHOLI
USE MODD_SNOW_PAR, ONLY : XPERCENTAGEPORE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)                     :: PSNOWDZ, PSNOWRHO, PSNOWLIQ
!
REAL, DIMENSION(SIZE(PSNOWRHO))                    :: PWHOLDMAX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
! computation of water holding capacity based on Crocus, 
!taking into account the conversion between wet and dry density -
!S. Morin/V. Vionnet 2010 12 09

! PWHOLDMAX is expressed in m water for each layer
! In short, PWHOLDMAX = XPERCENTAGEPORE * porosity * PSNOWDZ .
! The porosity is computed as (rho_ice - (rho_snow - lwc))/(rho_ice)
! where everything has to be in kg m-3 units. In practice, since
! PSNOWLIQ is expressed in m water, expressing the lwc in kg m-3
! is achieved as PSNOWLIQ*XRHOLW/PSNOWDZ. After some rearranging one
! obtains the equation given above.
! Note that equation (19) in Vionnet et al., GMD 2012, is wrong,
! because it does not take into account the fact that liquid water
! content has to be substracted from total density to compute the
! porosity.

IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_1D',0,ZHOOK_HANDLE)
PWHOLDMAX(:) = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ * (XRHOLI-PSNOWRHO) + PSNOWLIQ*XRHOLW)
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_1D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOWCROHOLD_1D
!####################################################################
      FUNCTION SNOWCROHOLD_0D(PSNOWRHO,PSNOWLIQ,PSNOWDZ) RESULT(PWHOLDMAX)
!
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s).
!
USE MODD_CSTS,     ONLY : XRHOLW,XRHOLI
USE MODD_SNOW_PAR, ONLY : XPERCENTAGEPORE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)        :: PSNOWDZ, PSNOWRHO, PSNOWLIQ
!
REAL                    :: PWHOLDMAX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
! computation of water holding capacity based on Crocus, 
!taking into account the conversion between wet and dry density - 
!S. Morin/V. Vionnet 2010 12 09

! PWHOLDMAX is expressed in m water for each layer
! In short, PWHOLDMAX = XPERCENTAGEPORE * porosity * PSNOWDZ .
! The porosity is computed as (rho_ice - (rho_snow - lwc))/(rho_ice)
! where everything has to be in kg m-3 units. In practice, since
! PSNOWLIQ is expressed in m water, expressing the lwc in kg m-3
! is achieved as PSNOWLIQ*XRHOLW/PSNOWDZ. After some rearranging one
! obtains the equation given above.
! Note that equation (19) in Vionnet et al., GMD 2012, is wrong,
! because it does not take into account the fact that liquid water
! content has to be substracted from total density to compute the
! porosity.

IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_0D',0,ZHOOK_HANDLE)
PWHOLDMAX = XPERCENTAGEPORE/XRHOLI * (PSNOWDZ * (XRHOLI-PSNOWRHO) + PSNOWLIQ*XRHOLW)
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_0D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOWCROHOLD_0D
!####################################################################
      FUNCTION SNOWO04HOLD_0D(PSNOWRHO,PSNOWLIQ,PSNOWDZ) RESULT(PWHOLDMAX)
!     Cluzet et al 2016
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s), with the CLM model for max lwc (similar SNOWCROHOLD_0D) 
!     see Oleson et al. 2004
!
USE MODD_CSTS,     ONLY : XRHOLW,XRHOLI
USE MODD_SNOW_PAR, ONLY : XPERCENTAGEPORE_O04
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)        :: PSNOWDZ, PSNOWRHO, PSNOWLIQ
!
REAL                    :: PWHOLDMAX
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_0D',0,ZHOOK_HANDLE)
PWHOLDMAX = XPERCENTAGEPORE_O04/XRHOLI * (PSNOWDZ * (XRHOLI-PSNOWRHO) + PSNOWLIQ*XRHOLW)
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWCROHOLD_0D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOWO04HOLD_0D
!####################################################################
        FUNCTION SNOWSPKHOLD_0D(PSNOWRHO,PSNOWLIQ,PSNOWDZ) RESULT(PWHOLDMAX)
!     Lafaysse et al 2017
!!    PURPOSE
!!    -------
!     Calculate the maximum liquid water holding capacity of
!     snow layer(s), with the SNOWPACK model for maximum volumetric water content
!     see Wever 
!
USE MODD_CSTS,     ONLY: XRHOLW,XRHOLI
USE MODD_SNOW_METAMO, ONLY : XUEPSI


USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB

IMPLICIT NONE

REAL, INTENT(IN)        :: PSNOWDZ, PSNOWRHO, PSNOWLIQ
REAL                    :: PWHOLDMAX
REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL                    :: ZTHETAI

IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWSPKHOLD_0D',0,ZHOOK_HANDLE)
! PSNOWLIQ in m
! PSNOWLIQ*XRHOLW/PSNOWDZ liquid water content in kg/m3
IF (PSNOWDZ>XUEPSI) THEN
  ZTHETAI=(PSNOWRHO-PSNOWLIQ*XRHOLW/PSNOWDZ)/XRHOLI
ELSE
  ZTHETAI=PSNOWRHO/XRHOLI
END IF
! In equation 12 Lafaysse et al 2017, capacity in kg m-3
! Here capacity in m, so eq 12 is multiplied by PSNOWDZ/XRHOLW
! 
IF (ZTHETAI<=0.23) THEN
  PWHOLDMAX = PSNOWDZ * ( 0.08- 0.1023 * (ZTHETAI - 0.03))
ELSEIF(ZTHETAI>0.812) THEN
  PWHOLDMAX = 0.
ELSE
  PWHOLDMAX = PSNOWDZ * ( 0.0264+0.0099*(1-ZTHETAI)/ZTHETAI )

ENDIF

IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOWSPKHOLD_0D',1,ZHOOK_HANDLE)

END FUNCTION SNOWSPKHOLD_0D
!####################################################################
      FUNCTION SNOW3LSCAP_3D(PSNOWRHO) RESULT(PSCAP)
!
!!    PURPOSE
!!    -------
!     Calculate the heat capacity of a snow layer.
!
USE MODD_CSTS,ONLY : XCI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN)                                  :: PSNOWRHO
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2),SIZE(PSNOWRHO,3)) :: PSCAP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!     The method of Verseghy (1991), Int. J. Climat., 11, 111-133.
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LSCAP_3D',0,ZHOOK_HANDLE)
PSCAP(:,:,:) = PSNOWRHO(:,:,:)*XCI      ! (J K-1 m-3)
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LSCAP_3D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LSCAP_3D
!####################################################################
      FUNCTION SNOW3LSCAP_2D(PSNOWRHO) RESULT(PSCAP)
!
!!    PURPOSE
!!    -------
!     Calculate the heat capacity of a snow layer.
!
USE MODD_CSTS,ONLY : XCI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)                   :: PSNOWRHO
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: PSCAP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!     The method of Verseghy (1991), Int. J. Climat., 11, 111-133.
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LSCAP_2D',0,ZHOOK_HANDLE)
PSCAP(:,:) = PSNOWRHO(:,:)*XCI      ! (J K-1 m-3)
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LSCAP_2D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LSCAP_2D
!####################################################################
      FUNCTION SNOW3LSCAP_1D(PSNOWRHO) RESULT(PSCAP)
!
!!    PURPOSE
!!    -------
!     Calculate the heat capacity of a snow layer.
!
USE MODD_CSTS,ONLY : XCI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)                   :: PSNOWRHO
!
REAL, DIMENSION(SIZE(PSNOWRHO))                  :: PSCAP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!     The method of Verseghy (1991), Int. J. Climat., 11, 111-133.
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LSCAP_1D',0,ZHOOK_HANDLE)
PSCAP(:) = PSNOWRHO(:)*XCI      ! (J K-1 m-3)
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LSCAP_1D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LSCAP_1D
!####################################################################
      FUNCTION SNOW3LSCAP_0D(PSNOWRHO) RESULT(PSCAP)
!
!!    PURPOSE
!!    -------
!     Calculate the heat capacity of a snow layer.
!
USE MODD_CSTS,ONLY : XCI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)       :: PSNOWRHO
!
REAL                   :: PSCAP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!     The method of Verseghy (1991), Int. J. Climat., 11, 111-133.
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LSCAP_0D',0,ZHOOK_HANDLE)
PSCAP = PSNOWRHO*XCI      ! (J K-1 m-3)
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LSCAP_0D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LSCAP_0D
!
!####################################################################
!####################################################################
!####################################################################
      FUNCTION SNOW3L_MARBOUTY(PSNOWRHO,PSNOWTEMP,PGRADT) RESULT(PDANGL) 
!**** *ZDANGL* - CROISSANCE DES GRAINS NON DENDRITIQUES ET ANGULEUX .
!              - GROWTH RATES FOR NON DENDRITIC GRAINS WITH SPHERICITY=0 


!     OBJET.
!     ------

!**   INTERFACE.
!     ----------

!     *ZDANGL*(PST,PSRO,PGRADT)*

!        *PST*     -  TEMPERATURE DE LA STRATE DE NEIGE.
!        *PSRO*    -  MASSE VOLUMIQUE DE LA STRATE DE NEIGE.
!        *PGRADT*  -  GRADIENT DE TEMPERATURE AFFECTANT LA STRATE DE NEIGE.

!     METHODE.
!     --------
!     THE FUNCTION RETURN A VALUE BETWEEN 0 AND 1 WHICH IS USED IN THE DETERMINATION OF THE 
!     GROWTH RATE FOR THE CONSIDERED LAYER.
!     THIS VALUE (WITHOUT UNIT) IS THE PRODUCT OF 3 FUNCTIONS (WHICH HAVE THEIR SOLUTIONS BETWEEN 0 AND 1) :
!     F(TEMPERATURE) * H(DENSITY) * G(TEMPERATURE GRADIENT)

!     EXTERNES.
!     ---------

!     REFERENCES.
!     -----------
!        MARBOUTY D (1980) AN EXPERIMENTAL STUDY OF TEMPERATURE GRADIENT 
!                          METAMORPHISM J GLACIOLOGY

!     AUTEURS.
!     --------
!        DOMINIQUE MARBOUTY (FMARBO/GMARBO/HMARBO).

!     MODIFICATIONS.
!     --------------
!        08/95: YANNICK DANIELOU - CODAGE A LA NORME DOCTOR.
!        03/06: JM WILLEMET      - F90 AND SI UNITS
!        01/08: JM WILLEMET      - ERROR ON THE FORMULATION OF G FUNCTION (WITH GRADIENT) IS CORRECTED 

USE MODD_CSTS, ONLY : XTT
USE MODD_SNOW_METAMO  
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!     DECLARATIONS.
!     -------------
!
REAL ,INTENT(IN) :: PSNOWTEMP, PSNOWRHO, PGRADT
!
REAL             :: PDANGL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3L_MARBOUTY',0,ZHOOK_HANDLE)
!
PDANGL = 0.0
!
! INFLUENCE DE LA TEMPERATURE /TEMPERATURE INFLUENCE.
IF( PSNOWTEMP>=XTT-XVTANG1 ) THEN
  !
  IF ( PSNOWTEMP>=XTT-XVTANG2 ) THEN
    PDANGL = XVTANG4 + XVTANG5 * (XTT-PSNOWTEMP) / XVTANG6
  ELSEIF( PSNOWTEMP>=XTT-XVTANG3 ) THEN
    PDANGL = XVTANG7 - XVTANG8 * (XTT-XVTANG2-PSNOWTEMP) / XVTANG9
  ELSE
    PDANGL = XVTANGA - XVTANGB * (XTT-XVTANG3-PSNOWTEMP) / XVTANGC
  ENDIF
  !
  ! INFLUENCE DE LA MASSE VOLUMIQUE / DENSITY INFLUENCE.
  IF ( PSNOWRHO<=XVRANG1 ) THEN
    !
    IF ( PSNOWRHO>XVRANG2 ) THEN
      PDANGL = PDANGL * ( 1. - (PSNOWRHO-XVRANG2)/(XVRANG1-XVRANG2) )
    ENDIF
    !
    ! INFLUENCE DU GRADIENT DE TEMPERATURE / TEMPERATURE GRADIENT INFLUENCE.
    IF ( PGRADT<=XVGANG1 ) THEN
      !
      IF ( PGRADT<=XVGANG2 ) THEN
        PDANGL = PDANGL * XVGANG5 * (PGRADT-XVGANG6)/(XVGANG2-XVGANG6)
      ELSEIF( PGRADT<=XVGANG3 ) THEN
        PDANGL = PDANGL * ( XVGANG7 + XVGANG8 * (PGRADT-XVGANG2)/(XVGANG3-XVGANG2) )
      ELSEIF( PGRADT<=XVGANG4 )THEN
        PDANGL = PDANGL * ( XVGANG9 + XVGANGA * (PGRADT-XVGANG3)/(XVGANG4-XVGANG3) )
      ELSE
        PDANGL = PDANGL * ( XVGANGB + XVGANGC * (PGRADT-XVGANG4)/(XVGANG1-XVGANG4) )
      ENDIF
      !
    ENDIF
    !
  ELSE
    !
    PDANGL = 0.
    !
  ENDIF
  !
ELSE
  !
  PDANGL = 0.
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3L_MARBOUTY',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3L_MARBOUTY     
!       
!####################################################################
!####################################################################
!####################################################################
!
      SUBROUTINE SNOW3LGRID_2D(PSNOWDZ,PSNOW,PSNOWDZ_OLD)
!
!!    PURPOSE
!!    -------
!     Once during each time step, update grid to maintain
!     grid proportions. Similar to approach of Lynch-Steiglitz,
!     1994, J. Clim., 7, 1842-1855. Corresponding mass and
!     heat adjustments made directly after the call to this
!     routine. 3 grid configurations:
!     1) for very thin snow, constant grid spacing
!     2) for intermediate thicknesses, highest resolution at soil/snow
!        interface and at the snow/atmosphere interface
!     3) for deep snow, vertical resoution finest at snow/atmosphere
!        interface (set to a constant value) and increases with snow depth.
!        Second layer can't be more than an order of magnitude thicker
!        than surface layer.
!
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_SNOW_PAR,   ONLY : XSNOWCRITD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:  ), INTENT(IN )           :: PSNOW
REAL, DIMENSION(:,:), INTENT(OUT)           :: PSNOWDZ
REAL, DIMENSION(:,:), INTENT(IN ), OPTIONAL :: PSNOWDZ_OLD
!
!*      0.1    declarations of local variables
!
INTEGER                           :: JJ, JI
!
INTEGER                           :: INLVLS, INI
!   
REAL,     DIMENSION(SIZE(PSNOW))  :: ZWORK
!
LOGICAL , DIMENSION(SIZE(PSNOW))  :: GREGRID

! ISBA-ES snow grid parameters
!
REAL, PARAMETER, DIMENSION(3)     :: ZSGCOEF1  = (/0.25, 0.50, 0.25/) 
REAL, PARAMETER, DIMENSION(2)     :: ZSGCOEF2  = (/0.05, 0.34/)       
!      
REAL, PARAMETER, DIMENSION(3)     :: ZSGCOEF   = (/0.3, 0.4, 0.3/) 
!
! Minimum total snow depth at which surface layer thickness is constant:
!
REAL, PARAMETER                   :: ZSNOWTRANS = 0.20                ! (m)
!      
! Minimum snow depth by layer for 6-L or 12-L configuration :
!
REAL, PARAMETER                   ::  ZDZ1=0.01
REAL, PARAMETER                   ::  ZDZ2=0.05
REAL, PARAMETER                   ::  ZDZ3=0.15
REAL, PARAMETER                   ::  ZDZ4=0.50
REAL, PARAMETER                   ::  ZDZ5=1.00
REAL, PARAMETER                   ::  ZDZN0=0.02
REAL, PARAMETER                   ::  ZDZN1=0.1
REAL, PARAMETER                   ::  ZDZN2=0.5
REAL, PARAMETER                   ::  ZDZN3=1.0
!
REAL, PARAMETER                   ::  ZCOEF1 = 0.5
REAL, PARAMETER                   ::  ZCOEF2 = 1.5
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LGRID_2D',0,ZHOOK_HANDLE)
!
INLVLS = SIZE(PSNOWDZ(:,:),2)
INI    = SIZE(PSNOWDZ(:,:),1)
!
ZWORK  (:) = 0.0
GREGRID(:) = .TRUE.
!
! 1. Calculate current grid for 3-layer (default) configuration):
! ---------------------------------------------------------------
! Based on formulation of Lynch-Stieglitz (1994)
! except for 3 modifications: 
! i) smooth transition here at ZSNOWTRANS
! ii) constant ratio for very thin snow:
! iii) ratio of layer 2 to surface layer <= 10
!
IF(INLVLS == 1)THEN
!
  DO JI=1,INI
     PSNOWDZ(JI,1)  = PSNOW(JI)
  ENDDO
!
ELSEIF(INLVLS == 3)THEN
!
   WHERE(PSNOW <= XSNOWCRITD+0.01)
      PSNOWDZ(:,1) = MIN(0.01, PSNOW(:)/INLVLS)
      PSNOWDZ(:,3) = MIN(0.01, PSNOW(:)/INLVLS)
      PSNOWDZ(:,2) = PSNOW(:) - PSNOWDZ(:,1) - PSNOWDZ(:,3)
   END WHERE
!
   WHERE(PSNOW <= ZSNOWTRANS .AND. PSNOW > XSNOWCRITD+0.01)
      PSNOWDZ(:,1) = PSNOW(:)*ZSGCOEF1(1)
      PSNOWDZ(:,2) = PSNOW(:)*ZSGCOEF1(2)
      PSNOWDZ(:,3) = PSNOW(:)*ZSGCOEF1(3)
   END WHERE
!
   WHERE(PSNOW > ZSNOWTRANS)
      PSNOWDZ(:,1) = ZSGCOEF2(1)
      PSNOWDZ(:,2) = (PSNOW(:)-ZSGCOEF2(1))*ZSGCOEF2(2) + ZSGCOEF2(1)
!
! When using simple finite differences, limit the thickness
! factor between the top and 2nd layers to at most 10
! 
      PSNOWDZ(:,2) = MIN(10*ZSGCOEF2(1),  PSNOWDZ(:,2))
      PSNOWDZ(:,3) = PSNOW(:) - PSNOWDZ(:,2) - PSNOWDZ(:,1)
   END WHERE
!
!
! 2. Calculate current grid for 6-layer :
! ---------------------------------------------------------------
!
ELSEIF(INLVLS == 6)THEN
!
! critere a satisfaire pour remaillage
!
  IF(PRESENT(PSNOWDZ_OLD))THEN
    GREGRID(:) = PSNOWDZ_OLD(:,1) < ZCOEF1 * MIN(ZDZ1 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,1) > ZCOEF2 * MIN(ZDZ1 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,2) < ZCOEF1 * MIN(ZDZ2 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,2) > ZCOEF2 * MIN(ZDZ2 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,6) < ZCOEF1 * MIN(ZDZN1,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,6) > ZCOEF2 * MIN(ZDZN1,PSNOW(:)/INLVLS)
  ENDIF
!
  WHERE(GREGRID(:))  
!      top layers 
       PSNOWDZ(:,1) = MIN(ZDZ1,PSNOW(:)/INLVLS) 
       PSNOWDZ(:,2) = MIN(ZDZ2,PSNOW(:)/INLVLS) 
!      last layers 
       PSNOWDZ(:,6) = MIN(ZDZN1,PSNOW(:)/INLVLS)
!      remaining snow for remaining layers
       ZWORK(:)     = PSNOW(:) - PSNOWDZ(:,1) - PSNOWDZ(:,2) - PSNOWDZ(:,6)
       PSNOWDZ(:,3) = ZWORK(:)*ZSGCOEF(1)
       PSNOWDZ(:,4) = ZWORK(:)*ZSGCOEF(2)
       PSNOWDZ(:,5) = ZWORK(:)*ZSGCOEF(3)
!      layer 3 tickness >= layer 2 tickness
       ZWORK(:)=MIN(0.0,PSNOWDZ(:,3)-PSNOWDZ(:,2))
       PSNOWDZ(:,3)=PSNOWDZ(:,3)-ZWORK(:)
       PSNOWDZ(:,4)=PSNOWDZ(:,4)+ZWORK(:) 
!      layer 5 tickness >= layer 6 tickness  
       ZWORK(:)=MIN(0.0,PSNOWDZ(:,5)-PSNOWDZ(:,6))
       PSNOWDZ(:,5)=PSNOWDZ(:,5)-ZWORK(:)
       PSNOWDZ(:,4)=PSNOWDZ(:,4)+ZWORK(:)
  ENDWHERE
!
! 3. Calculate current grid for 9-layer :
! ---------------------------------------------------------------
!
ELSEIF(INLVLS == 9)THEN
!
! critere a satisfaire pour remaillage
!
  IF(PRESENT(PSNOWDZ_OLD))THEN
    GREGRID(:) = PSNOWDZ_OLD(:,1) < ZCOEF1 * MIN(ZDZ1 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,1) > ZCOEF2 * MIN(ZDZ1 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,2) < ZCOEF1 * MIN(ZDZ2 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,2) > ZCOEF2 * MIN(ZDZ2 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,9) < ZCOEF1 * MIN(ZDZN0,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,9) > ZCOEF2 * MIN(ZDZN0,PSNOW(:)/INLVLS) 
  ENDIF
!             
  WHERE(GREGRID(:))  
!      top layers 
       PSNOWDZ(:,1) = MIN(ZDZ1,PSNOW(:)/INLVLS) 
       PSNOWDZ(:,2) = MIN(ZDZ2,PSNOW(:)/INLVLS) 
       PSNOWDZ(:,3) = MIN(ZDZ3,PSNOW(:)/INLVLS)
!      last layers 
       PSNOWDZ(:,9)= MIN(ZDZN0,PSNOW(:)/INLVLS)
       PSNOWDZ(:,8)= MIN(ZDZN1,PSNOW(:)/INLVLS)
       PSNOWDZ(:,7)= MIN(ZDZN2,PSNOW(:)/INLVLS)
!      remaining snow for remaining layers
       ZWORK(:) = PSNOW(:) - PSNOWDZ(:, 1) - PSNOWDZ(:, 2) - PSNOWDZ(:, 3) &
                           - PSNOWDZ(:, 7) - PSNOWDZ(:, 8) - PSNOWDZ(:, 9)
       PSNOWDZ(:,4) = ZWORK(:)*ZSGCOEF(1)
       PSNOWDZ(:,5) = ZWORK(:)*ZSGCOEF(2)
       PSNOWDZ(:,6) = ZWORK(:)*ZSGCOEF(3)
!      layer 4 tickness >= layer 3 tickness
       ZWORK(:)=MIN(0.0,PSNOWDZ(:,4)-PSNOWDZ(:,3))
       PSNOWDZ(:,4)=PSNOWDZ(:,4)-ZWORK(:)
       PSNOWDZ(:,5)=PSNOWDZ(:,5)+ZWORK(:) 
!      layer 6 tickness >= layer 7 tickness  
       ZWORK(:)=MIN(0.0,PSNOWDZ(:,6)-PSNOWDZ(:,7))
       PSNOWDZ(:,6)=PSNOWDZ(:,6)-ZWORK(:)
       PSNOWDZ(:,5)=PSNOWDZ(:,5)+ZWORK(:)
  ENDWHERE
!
! 4. Calculate current grid for 12-layer :
! ---------------------------------------------------------------
!
ELSEIF(INLVLS == 12)THEN
!
! critere a satisfaire pour remaillage
!
  IF(PRESENT(PSNOWDZ_OLD))THEN
    GREGRID(:) = PSNOWDZ_OLD(:, 1) < ZCOEF1 * MIN(ZDZ1 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:, 1) > ZCOEF2 * MIN(ZDZ1 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:, 2) < ZCOEF1 * MIN(ZDZ2 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:, 2) > ZCOEF2 * MIN(ZDZ2 ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,12) < ZCOEF1 * MIN(ZDZN0,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,12) > ZCOEF2 * MIN(ZDZN0,PSNOW(:)/INLVLS) 
  ENDIF
!             
  WHERE(GREGRID(:))  
!      top layers 
       PSNOWDZ(:,1) = MIN(ZDZ1,PSNOW(:)/INLVLS) 
       PSNOWDZ(:,2) = MIN(ZDZ2,PSNOW(:)/INLVLS) 
       PSNOWDZ(:,3) = MIN(ZDZ3,PSNOW(:)/INLVLS)
       PSNOWDZ(:,4) = MIN(ZDZ4,PSNOW(:)/INLVLS)
       PSNOWDZ(:,5) = MIN(ZDZ5,PSNOW(:)/INLVLS)
!      last layers 
       PSNOWDZ(:,12)= MIN(ZDZN0,PSNOW(:)/INLVLS)
       PSNOWDZ(:,11)= MIN(ZDZN1,PSNOW(:)/INLVLS)
       PSNOWDZ(:,10)= MIN(ZDZN2,PSNOW(:)/INLVLS)
       PSNOWDZ(:, 9)= MIN(ZDZN3,PSNOW(:)/INLVLS)
!      remaining snow for remaining layers
       ZWORK(:) = PSNOW(:) - PSNOWDZ(:, 1) - PSNOWDZ(:, 2) - PSNOWDZ(:, 3) &
                           - PSNOWDZ(:, 4) - PSNOWDZ(:, 5) - PSNOWDZ(:, 9) &
                           - PSNOWDZ(:,10) - PSNOWDZ(:,11) - PSNOWDZ(:,12)
       PSNOWDZ(:,6) = ZWORK(:)*ZSGCOEF(1)
       PSNOWDZ(:,7) = ZWORK(:)*ZSGCOEF(2)
       PSNOWDZ(:,8) = ZWORK(:)*ZSGCOEF(3)
!      layer 6 tickness >= layer 5 tickness
       ZWORK(:)=MIN(0.0,PSNOWDZ(:,6)-PSNOWDZ(:,5))
       PSNOWDZ(:,6)=PSNOWDZ(:,6)-ZWORK(:)
       PSNOWDZ(:,7)=PSNOWDZ(:,7)+ZWORK(:) 
!      layer 8 tickness >= layer 9 tickness  
       ZWORK(:)=MIN(0.0,PSNOWDZ(:,8)-PSNOWDZ(:,9))
       PSNOWDZ(:,8)=PSNOWDZ(:,8)-ZWORK(:)
       PSNOWDZ(:,7)=PSNOWDZ(:,7)+ZWORK(:)
  ENDWHERE
!
! 4. Calculate other non-optimized grid :
! ---------------------------------------
! 
ELSEIF(INLVLS<10.AND.INLVLS/=3.AND.INLVLS/=6.AND.INLVLS/=9) THEN
!
  DO JJ=1,INLVLS
     DO JI=1,INI
        PSNOWDZ(JI,JJ)  = PSNOW(JI)/INLVLS
     ENDDO
  ENDDO
!
  PSNOWDZ(:,INLVLS) = PSNOWDZ(:,INLVLS) + (PSNOWDZ(:,1) - MIN(0.05, PSNOWDZ(:,1)))
  PSNOWDZ(:,1)      = MIN(0.05, PSNOWDZ(:,1))
!
ELSE !(INLVLS>=10 and /=12)  
!
! critere a satisfaire pour remaillage
!
  IF(PRESENT(PSNOWDZ_OLD))THEN
    GREGRID(:) = PSNOWDZ_OLD(:,     1) < ZCOEF1 * MIN(ZDZ1         ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,     1) > ZCOEF2 * MIN(ZDZ1         ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,     2) < ZCOEF1 * MIN(ZDZ2         ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,     2) > ZCOEF2 * MIN(ZDZ2         ,PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,INLVLS) < ZCOEF1 * MIN(0.05*PSNOW(:),PSNOW(:)/INLVLS) .OR. &
               & PSNOWDZ_OLD(:,INLVLS) > ZCOEF2 * MIN(0.05*PSNOW(:),PSNOW(:)/INLVLS) 
  ENDIF
!
  WHERE(GREGRID(:))  
       PSNOWDZ(:,1     ) = MIN(ZDZ1         ,PSNOW(:)/INLVLS) 
       PSNOWDZ(:,2     ) = MIN(ZDZ2         ,PSNOW(:)/INLVLS) 
       PSNOWDZ(:,3     ) = MIN(ZDZ3         ,PSNOW(:)/INLVLS)
       PSNOWDZ(:,4     ) = MIN(ZDZ4         ,PSNOW(:)/INLVLS)
       PSNOWDZ(:,5     ) = MIN(ZDZ5         ,PSNOW(:)/INLVLS)          
       PSNOWDZ(:,INLVLS) = MIN(0.05*PSNOW(:),PSNOW(:)/INLVLS) 
  ENDWHERE
!
  DO JJ=6,INLVLS-1,1
     DO JI=1,INI
        IF(GREGRID(JI))THEN           
          ZWORK(JI) = PSNOWDZ(JI,1)+PSNOWDZ(JI,2)+PSNOWDZ(JI,3)+PSNOWDZ(JI,4)+PSNOWDZ(JI,5)
          PSNOWDZ(JI,JJ) = (PSNOW(JI)-ZWORK(JI)-PSNOWDZ(JI,INLVLS))/(INLVLS-6) 
        ENDIF
     ENDDO
  ENDDO
!
ENDIF
!
DO JI=1,INI
  IF(PSNOW(JI)==XUNDEF) PSNOWDZ(JI,:) = XUNDEF
ENDDO
!
IF (PRESENT(PSNOWDZ_OLD)) THEN
  DO JI=1,INI
    IF (.NOT.GREGRID(JI)) PSNOWDZ(JI,:)=PSNOWDZ_OLD(JI,:) 
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LGRID_2D',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW3LGRID_2D
!####################################################################
!####################################################################
!####################################################################
!
      SUBROUTINE SNOW3LGRID_1D(PSNOWDZ,PSNOW,PSNOWDZ_OLD)
!
!!    PURPOSE
!!    -------
!     Once during each time step, update grid to maintain
!     grid proportions. Similar to approach of Lynch-Steiglitz,
!     1994, J. Clim., 7, 1842-1855. Corresponding mass and
!     heat adjustments made directly after the call to this
!     routine. 3 grid configurations:
!     1) for very thin snow, constant grid spacing
!     2) for intermediate thicknesses, highest resolution at soil/snow
!        interface and at the snow/atmosphere interface
!     3) for deep snow, vertical resoution finest at snow/atmosphere
!        interface (set to a constant value) and increases with snow depth.
!        Second layer can't be more than an order of magnitude thicker
!        than surface layer.
!
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_SNOW_PAR,   ONLY : XSNOWCRITD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN )           :: PSNOW
REAL, DIMENSION(:), INTENT(OUT)           :: PSNOWDZ
REAL, DIMENSION(:), INTENT(IN ), OPTIONAL :: PSNOWDZ_OLD
!
!*      0.1    declarations of local variables
!
INTEGER JJ
!
INTEGER                           :: INLVLS
!
REAL                              :: ZWORK
!
! modif_EB pour maillage
LOGICAL                           :: GREGRID

! ISBA-ES snow grid parameters
!
REAL, PARAMETER, DIMENSION(3)     :: ZSGCOEF1  = (/0.25, 0.50, 0.25/) 
REAL, PARAMETER, DIMENSION(2)     :: ZSGCOEF2  = (/0.05, 0.34/)       
!      
REAL, PARAMETER, DIMENSION(3)     :: ZSGCOEF   = (/0.3, 0.4, 0.3/) 
!
! Minimum total snow depth at which surface layer thickness is constant:
!
REAL, PARAMETER                   :: ZSNOWTRANS  = 0.20                ! (m)
!      
! Minimum snow depth by layer for 6-L or 12-L configuration :
!
REAL, PARAMETER                   ::  ZDZ1=0.01
REAL, PARAMETER                   ::  ZDZ2=0.05
REAL, PARAMETER                   ::  ZDZ3=0.15
REAL, PARAMETER                   ::  ZDZ4=0.50
REAL, PARAMETER                   ::  ZDZ5=1.00
REAL, PARAMETER                   ::  ZDZN0=0.02
REAL, PARAMETER                   ::  ZDZN1=0.1
REAL, PARAMETER                   ::  ZDZN2=0.5
REAL, PARAMETER                   ::  ZDZN3=1.0
!
REAL, PARAMETER                   ::  ZCOEF1 = 0.5
REAL, PARAMETER                   ::  ZCOEF2 = 1.5
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LGRID_1D',0,ZHOOK_HANDLE)
!
INLVLS = SIZE(PSNOWDZ(:),1)
!
GREGRID = .TRUE.
!
! 1. Calculate current grid for 3-layer (default) configuration):
! ---------------------------------------------------------------
! Based on formulation of Lynch-Stieglitz (1994)
! except for 3 modifications: 
! i) smooth transition here at ZSNOWTRANS
! ii) constant ratio for very thin snow:
! iii) ratio of layer 2 to surface layer <= 10
!
IF(INLVLS == 1)THEN
!
  PSNOWDZ(1) = PSNOW
!
ELSEIF(INLVLS == 3)THEN
!
   IF(PSNOW <= XSNOWCRITD+0.01)THEN
      PSNOWDZ(1) = MIN(0.01, PSNOW/INLVLS)
      PSNOWDZ(3) = MIN(0.01, PSNOW/INLVLS)
      PSNOWDZ(2) = PSNOW - PSNOWDZ(1) - PSNOWDZ(3)
   ENDIF
!
   IF(PSNOW <= ZSNOWTRANS .AND. PSNOW > XSNOWCRITD+0.01)THEN
      PSNOWDZ(1) = PSNOW*ZSGCOEF1(1)
      PSNOWDZ(2) = PSNOW*ZSGCOEF1(2)
      PSNOWDZ(3) = PSNOW*ZSGCOEF1(3)
   ENDIF
!
   IF(PSNOW > ZSNOWTRANS)THEN
      PSNOWDZ(1) = ZSGCOEF2(1)
      PSNOWDZ(2) = (PSNOW-ZSGCOEF2(1))*ZSGCOEF2(2) + ZSGCOEF2(1)
!
! When using simple finite differences, limit the thickness
! factor between the top and 2nd layers to at most 10
! 
      PSNOWDZ(2) = MIN(10*ZSGCOEF2(1),  PSNOWDZ(2))
      PSNOWDZ(3) = PSNOW - PSNOWDZ(2) - PSNOWDZ(1)
   END IF
!
!
! 2. Calculate current grid for 6-layer :
! ---------------------------------------------------------------
!
ELSEIF(INLVLS == 6)THEN
!
! critere a satisfaire pour remaillage
!
  IF(PRESENT(PSNOWDZ_OLD))THEN
    GREGRID    = PSNOWDZ_OLD(1) < ZCOEF1 * MIN(ZDZ1 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(1) > ZCOEF2 * MIN(ZDZ1 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(2) < ZCOEF1 * MIN(ZDZ2 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(2) > ZCOEF2 * MIN(ZDZ2 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(6) < ZCOEF1 * MIN(ZDZN1,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(6) > ZCOEF2 * MIN(ZDZN1,PSNOW/INLVLS)
  ENDIF
!
  IF(GREGRID)THEN
!      top layers 
       PSNOWDZ(1) = MIN(ZDZ1,PSNOW/INLVLS) 
       PSNOWDZ(2) = MIN(ZDZ2,PSNOW/INLVLS) 
!      last layers 
       PSNOWDZ(6) = MIN(ZDZN1,PSNOW/INLVLS)
!      remaining snow for remaining layers
       ZWORK      = PSNOW - PSNOWDZ(1) - PSNOWDZ(2) - PSNOWDZ(6)
       PSNOWDZ(3) = ZWORK*ZSGCOEF(1)
       PSNOWDZ(4) = ZWORK*ZSGCOEF(2)
       PSNOWDZ(5) = ZWORK*ZSGCOEF(3)
!      layer 3 tickness >= layer 2 tickness
       ZWORK=MIN(0.0,PSNOWDZ(3)-PSNOWDZ(2))
       PSNOWDZ(3)=PSNOWDZ(3)-ZWORK
       PSNOWDZ(4)=PSNOWDZ(4)+ZWORK
!      layer 5 tickness >= layer 6 tickness  
       ZWORK=MIN(0.0,PSNOWDZ(5)-PSNOWDZ(6))
       PSNOWDZ(5)=PSNOWDZ(5)-ZWORK
       PSNOWDZ(4)=PSNOWDZ(4)+ZWORK
  ENDIF
!
! 3. Calculate current grid for 9-layer :
! ---------------------------------------------------------------
!
ELSEIF(INLVLS == 9)THEN
!
! critere a satisfaire pour remaillage
!
  IF(PRESENT(PSNOWDZ_OLD))THEN
    GREGRID    = PSNOWDZ_OLD(1) < ZCOEF1 * MIN(ZDZ1 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(1) > ZCOEF2 * MIN(ZDZ1 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(2) < ZCOEF1 * MIN(ZDZ2 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(2) > ZCOEF2 * MIN(ZDZ2 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(9) < ZCOEF1 * MIN(ZDZN0,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(9) > ZCOEF2 * MIN(ZDZN0,PSNOW/INLVLS) 
  ENDIF
!             
  IF(GREGRID)THEN
!      top layers 
       PSNOWDZ(1) = MIN(ZDZ1,PSNOW/INLVLS) 
       PSNOWDZ(2) = MIN(ZDZ2,PSNOW/INLVLS) 
       PSNOWDZ(3) = MIN(ZDZ3,PSNOW/INLVLS)
!      last layers 
       PSNOWDZ(9)= MIN(ZDZN0,PSNOW/INLVLS)
       PSNOWDZ(8)= MIN(ZDZN1,PSNOW/INLVLS)
       PSNOWDZ(7)= MIN(ZDZN2,PSNOW/INLVLS)
!      remaining snow for remaining layers
       ZWORK = PSNOW - PSNOWDZ( 1) - PSNOWDZ( 2) - PSNOWDZ( 3) &
                     - PSNOWDZ( 7) - PSNOWDZ( 8) - PSNOWDZ( 9)
       PSNOWDZ(4) = ZWORK*ZSGCOEF(1)
       PSNOWDZ(5) = ZWORK*ZSGCOEF(2)
       PSNOWDZ(6) = ZWORK*ZSGCOEF(3)
!      layer 4 tickness >= layer 3 tickness
       ZWORK=MIN(0.0,PSNOWDZ(4)-PSNOWDZ(3))
       PSNOWDZ(4)=PSNOWDZ(4)-ZWORK
       PSNOWDZ(5)=PSNOWDZ(5)+ZWORK
!      layer 6 tickness >= layer 7 tickness  
       ZWORK=MIN(0.0,PSNOWDZ(6)-PSNOWDZ(7))
       PSNOWDZ(6)=PSNOWDZ(6)-ZWORK
       PSNOWDZ(5)=PSNOWDZ(5)+ZWORK
  ENDIF
!
! 4. Calculate current grid for 12-layer :
! ---------------------------------------------------------------
!
ELSEIF(INLVLS == 12)THEN
!
! modif_EB pour maillage
! critere a satisfaire pour remaillage
  IF(PRESENT(PSNOWDZ_OLD))THEN
    GREGRID    = PSNOWDZ_OLD(1)  < ZCOEF1 * MIN(ZDZ1 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(1)  > ZCOEF2 * MIN(ZDZ1 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(2)  < ZCOEF1 * MIN(ZDZ2 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(2)  > ZCOEF2 * MIN(ZDZ2 ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(12) < ZCOEF1 * MIN(ZDZN0,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(12) > ZCOEF2 * MIN(ZDZN0,PSNOW/INLVLS) 
  ENDIF
!
  IF (GREGRID)THEN
!    top layers 
     PSNOWDZ(1) = MIN(ZDZ1,PSNOW/INLVLS) 
     PSNOWDZ(2) = MIN(ZDZ2,PSNOW/INLVLS) 
     PSNOWDZ(3) = MIN(ZDZ3,PSNOW/INLVLS)
     PSNOWDZ(4) = MIN(ZDZ4,PSNOW/INLVLS)
     PSNOWDZ(5) = MIN(ZDZ5,PSNOW/INLVLS)
!    last layers 
     PSNOWDZ(12)= MIN(ZDZN0,PSNOW/INLVLS)
     PSNOWDZ(11)= MIN(ZDZN1,PSNOW/INLVLS)
     PSNOWDZ(10)= MIN(ZDZN2,PSNOW/INLVLS)
     PSNOWDZ( 9)= MIN(ZDZN3,PSNOW/INLVLS)
!    remaining snow for remaining layers
     ZWORK = PSNOW - PSNOWDZ( 1) - PSNOWDZ( 2) - PSNOWDZ( 3) &
                   - PSNOWDZ( 4) - PSNOWDZ( 5) - PSNOWDZ( 9) &
                   - PSNOWDZ(10) - PSNOWDZ(11) - PSNOWDZ(12)
     PSNOWDZ(6) = ZWORK*ZSGCOEF(1)
     PSNOWDZ(7) = ZWORK*ZSGCOEF(2)
     PSNOWDZ(8) = ZWORK*ZSGCOEF(3)
!    layer 6 tickness >= layer 5 tickness
     ZWORK=MIN(0.0,PSNOWDZ(6)-PSNOWDZ(5))
     PSNOWDZ(6)=PSNOWDZ(6)-ZWORK
     PSNOWDZ(7)=PSNOWDZ(7)+ZWORK 
!    layer 8 tickness >= layer 9 tickness  
     ZWORK=MIN(0.0,PSNOWDZ(8)-PSNOWDZ(9))
     PSNOWDZ(8)=PSNOWDZ(8)-ZWORK
     PSNOWDZ(7)=PSNOWDZ(7)+ZWORK
  ENDIF
!
! 4. Calculate other non-optimized grid to allow CROCUS PREP :
! ------------------------------------------------------------
! 
ELSE IF(INLVLS<10.AND.INLVLS/=3.AND.INLVLS/=6.AND.INLVLS/=9) THEN
!        
  DO JJ=1,INLVLS
     PSNOWDZ(JJ)  = PSNOW/INLVLS
  ENDDO
!
  PSNOWDZ(INLVLS) = PSNOWDZ(INLVLS) + (PSNOWDZ(1) - MIN(0.05, PSNOWDZ(1)))
  PSNOWDZ(1)      = MIN(0.05, PSNOWDZ(1))
!
ELSE   
!
  IF(PRESENT(PSNOWDZ_OLD))THEN
    GREGRID    = PSNOWDZ_OLD(     1) < ZCOEF1 * MIN(ZDZ1      ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(     1) > ZCOEF2 * MIN(ZDZ1      ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(     2) < ZCOEF1 * MIN(ZDZ2      ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(     2) > ZCOEF2 * MIN(ZDZ2      ,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(INLVLS) < ZCOEF1 * MIN(0.05*PSNOW,PSNOW/INLVLS) .OR. &
               & PSNOWDZ_OLD(INLVLS) > ZCOEF2 * MIN(0.05*PSNOW,PSNOW/INLVLS) 
  ENDIF
!
  IF (GREGRID)THEN
     PSNOWDZ(     1) = MIN(ZDZ1      ,PSNOW/INLVLS) 
     PSNOWDZ(     2) = MIN(ZDZ2      ,PSNOW/INLVLS)
     PSNOWDZ(     3) = MIN(ZDZ3      ,PSNOW/INLVLS)
     PSNOWDZ(     4) = MIN(ZDZ4      ,PSNOW/INLVLS)
     PSNOWDZ(     5) = MIN(ZDZ5      ,PSNOW/INLVLS)
     PSNOWDZ(INLVLS) = MIN(0.05*PSNOW,PSNOW/INLVLS) 
     ZWORK = SUM(PSNOWDZ(1:5))          
     DO JJ=6,INLVLS-1,1
        PSNOWDZ(JJ) = (PSNOW - ZWORK -PSNOWDZ(INLVLS))/(INLVLS-6) 
     END DO
  ENDIF
!
ENDIF
!
IF (PSNOW==XUNDEF) PSNOWDZ(:) = XUNDEF
!
IF (PRESENT(PSNOWDZ_OLD)) THEN
  IF (.NOT.GREGRID) PSNOWDZ(:) = PSNOWDZ_OLD(:)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LGRID_1D',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW3LGRID_1D
!
!###################################################################################
!###################################################################################
SUBROUTINE SNOW3LAGREG(PSNOWDZN,PSNOWDZ,PSNOWRHO,PSNOWGRAN1, PSNOWGRAN2, &
                       PSNOWHIST,PSNOWGRAN1N,PSNOWGRAN2N,PSNOWHISTN,     &
                       KL1,KL2,PSNOWDDZ                        ) 
!
USE MODD_SNOW_METAMO
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!       0.1 declarations of arguments        
!        
REAL, DIMENSION(:), INTENT(IN)  :: PSNOWDZN,PSNOWDZ,PSNOWRHO,PSNOWDDZ  
!                                                    
REAL, DIMENSION(:), INTENT(IN)  :: PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST 
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWGRAN1N,PSNOWGRAN2N,PSNOWHISTN                                              
!
INTEGER, INTENT(IN) :: KL1  ! Indice couche de reference (i)
INTEGER, INTENT(IN) :: KL2 ! Indice de la couche (i-1 ou i+1) dont une 
                                ! partie est aggregee à la couche (i)
!
!       0.2 declaration of local variables
!        
REAL, DIMENSION(SIZE(PSNOWRHO,1)) :: ZSNOWRHO
REAL, DIMENSION(SIZE(PSNOWRHO,1)) :: ZDIAMD,ZDIAMV,ZSPHERD,ZSPHERV,&
                                     ZDIAMN,ZSPHERN,ZDENT 
!
REAL :: ZDELTA, ZCOMP
!
INTEGER :: IDENT, IVIEU, IL
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE 
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LAGREG',0,ZHOOK_HANDLE)
!
IF( KL1<KL2 ) THEN
  ZDELTA = 0.0
  IL = KL1
ELSE
  ZDELTA = 1.0
  IL = KL2
ENDIF
!
! Mean Properties
!
!       1. History
!
IF ( PSNOWHIST(KL1)/=PSNOWHIST(KL2) ) THEN
  PSNOWHISTN(KL1) = 0.0
ENDIF
!
!       2. New grain types
!
!       2.1 Same grain type
!
IF (  PSNOWGRAN1(KL1)*PSNOWGRAN1(KL2)>0.0 .OR. &
    ( PSNOWGRAN1(KL1)==0.0 .AND. PSNOWGRAN1(KL2)>=0.0 ) .OR. &
    ( PSNOWGRAN1(KL2)==0.0 .AND. PSNOWGRAN1(KL1)>=0.0 ) ) THEN 
  !
  !code original vincent          PSNOWGRAN1N(KL1)=(PSNOWGRAN1(KL1)*PSNOWRHO(KL1)&
  !code original vincent        *(PSNOWDZN(KL1)-(1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))-ZDELTA*&
  !code original vincent        ABS(PSNOWDDZ(KL2)))+PSNOWGRAN1(KL2)*                   &
  !code original vincent        PSNOWRHO(KL2)*((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))+        &   
  !code original vincent        ZDELTA*ABS(PSNOWDDZ(KL2))))/((PSNOWDZN(KL1)-(1.0-ZDELTA)&
  !code original vincent        *ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2)))*        &
  !code original vincent        PSNOWRHO(KL1)+((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))+        &   
  !code original vincent        ZDELTA*ABS(PSNOWDDZ(KL2)))*PSNOWRHO(KL2))
  !code original vincent !
  !code original vincent          PSNOWGRAN2N(KL1)=(PSNOWGRAN2(KL1)*PSNOWRHO(KL1) &
  !code original vincent        *(PSNOWDZN(KL1)-(1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))-ZDELTA* &
  !code original vincent        ABS(PSNOWDDZ(KL2)))+PSNOWGRAN2(KL2)*                   &
  !code original vincent        PSNOWRHO(KL2)*((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))        &
  !code original vincent        +ZDELTA*ABS(PSNOWDDZ(KL2))))/((PSNOWDZN(KL1)-(1.0-ZDELTA)&
  !code original vincent        *ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2)))*        &
  !code original vincent        PSNOWRHO(KL1)+((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))+        &   
  !code original vincent        ZDELTA*ABS(PSNOWDDZ(KL2)))*PSNOWRHO(KL2))
  !     
  !plm
  CALL GET_AGREG(KL1,KL2,PSNOWGRAN1(KL1),PSNOWGRAN1(KL2),PSNOWGRAN1N(KL1))
  !
  CALL GET_AGREG(KL1,KL2,PSNOWGRAN2(KL1),PSNOWGRAN2(KL2),PSNOWGRAN2N(KL1))
  !
  !plm
  !     
ELSE
  !
  !       2.2 Different types
  !        
  IF ( PSNOWGRAN1(KL1)<0.0 ) THEN
    IDENT = KL1
    IVIEU = KL2
  ELSE
    IDENT = KL2
    IVIEU = KL1
  ENDIF  
  !                        
  ZDIAMD (KL1) = - PSNOWGRAN1(IDENT)/XGRAN * XDIAET + ( 1.0 + PSNOWGRAN1(IDENT)/XGRAN ) * &
                 ( PSNOWGRAN2(IDENT)/XGRAN * XDIAGF + ( 1.0 - PSNOWGRAN2(IDENT)/XGRAN ) * XDIAFP )
  !
  ZSPHERD(KL1) = PSNOWGRAN2(IDENT)/XGRAN                
  ZDIAMV (KL1) = PSNOWGRAN2(IVIEU)
  ZSPHERV(KL1) = PSNOWGRAN1(IVIEU)/XGRAN
  !IF(KL1==1)THEN
  !write(*,*) 'ZDD1',ZDIAMD(1),'ZSD1',ZSPHERD(1)
  !write(*,*) 'ZDV1',ZDIAMV(1),'ZSV1',ZSPHERV(1)
  !ENDIF       
  !
  IF ( IDENT==KL1 ) THEN
    !code original vincent        ZDIAMN(KL1)= (ZDIAMD(KL1)*PSNOWRHO(IDENT)*&
    !code original vincent            (PSNOWDZN(IDENT)-(1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))-ZDELTA*      &
    !code original vincent            ABS(PSNOWDDZ(KL2)))+ZDIAMV(KL1)*PSNOWRHO(IVIEU)*(       &
    !code original vincent            (1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2))))/&
    !code original vincent            ((PSNOWDZN(KL1)-(1.0-ZDELTA)*                                    &
    !code original vincent            ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2)))*            &
    !code original vincent            PSNOWRHO(KL1)+((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))+          &   
    !code original vincent            ZDELTA*ABS(PSNOWDDZ(KL2)))*PSNOWRHO(KL2))
    !
    !plm
    CALL GET_AGREG(IDENT,IVIEU,ZDIAMD(KL1),ZDIAMV(KL1),ZDIAMN(KL1))
    !
    !plm
    !         
    !code original vincent        ZSPHERN(KL1)= (ZSPHERD(KL1)*PSNOWRHO(IDENT)*&
    !code original vincent            (PSNOWDZN(IDENT)-(1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))-ZDELTA*      &
    !code original vincent            ABS(PSNOWDDZ(KL2)))+ZSPHERV(KL1)*PSNOWRHO(IVIEU)*(       &
    !code original vincent            (1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2))))/&
    !code original vincent            ((PSNOWDZN(KL1)-(1.0-ZDELTA)*                                    &
    !code original vincent            ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2)))*            &
    !code original vincent            PSNOWRHO(KL1)+((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))+          &   
    !code original vincent            ZDELTA*ABS(PSNOWDDZ(KL2)))*PSNOWRHO(KL2))

    !plm
    CALL GET_AGREG(IDENT,IVIEU,ZSPHERD(KL1),ZSPHERV(KL1),ZSPHERN(KL1))   
    !plm
    !
  ELSE
    !code original vincent        ZDIAMN(KL1)= (ZDIAMD(KL1)*PSNOWRHO(IDENT)*&
    !code original vincent            ((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))+ZDELTA*ABS(PSNOWDDZ(KL2)))&
    !code original vincent            +ZDIAMV(KL1)*PSNOWRHO(IVIEU)*(PSNOWDZN(IVIEU)-(1.0-ZDELTA)*  & 
    !code original vincent            ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2))))/&
    !code original vincent            ((PSNOWDZN(KL1)-(1.0-ZDELTA)*                          &
    !code original vincent            ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2)))*    &
    !code original vincent            PSNOWRHO(KL1)+((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))+   &
    !code original vincent            ZDELTA*ABS(PSNOWDDZ(KL2)))*PSNOWRHO(KL2))
    !code original vincent!            
    !code original vincent         ZSPHERN(KL1)= (ZSPHERD(KL1)*PSNOWRHO(IDENT)*&
    !code original vincent            ((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))+ZDELTA*ABS(PSNOWDDZ(KL2)))&
    !code original vincent            +ZSPHERV(KL1)*PSNOWRHO(IVIEU)*(PSNOWDZN(IVIEU)-(1.0-ZDELTA)* & 
    !code original vincent            ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2))))/&
    !code original vincent            ((PSNOWDZN(KL1)-(1.0-ZDELTA)*                                    &
    !code original vincent            ABS(PSNOWDDZ(KL1))-ZDELTA*ABS(PSNOWDDZ(KL2)))*            &
    !code original vincent            PSNOWRHO(KL1)+((1.0-ZDELTA)*ABS(PSNOWDDZ(KL1))+          &   
    !code original vincent            ZDELTA*ABS(PSNOWDDZ(KL2)))*PSNOWRHO(KL2))
    !plm
    !
    CALL GET_AGREG(IVIEU,IDENT,ZDIAMV(KL1),ZDIAMD(KL1),ZDIAMN(KL1))
    !           
    CALL GET_AGREG(IVIEU,IDENT,ZSPHERV(KL1),ZSPHERD(KL1),ZSPHERN(KL1))
    !plm
    !
  ENDIF
  !       
  ZCOMP = ZSPHERN(KL1) * XDIAGF + ( 1.-ZSPHERN(KL1) ) * XDIAFP
  IF( ZDIAMN(KL1) < ZCOMP ) THEN
    !
    ZDENT(KL1) = ( ZDIAMN(KL1) - ZCOMP ) / ( XDIAET - ZCOMP ) 
    !IF(KL1==1) write(*,*) 'ZDENT',ZDENT(1)   
    PSNOWGRAN1N(KL1) = - XGRAN * ZDENT  (KL1)
    PSNOWGRAN2N(KL1) =   XGRAN * ZSPHERN(KL1)
    !
  ELSE
    !
    PSNOWGRAN1N(KL1) = XGRAN * ZSPHERN(KL1)
    PSNOWGRAN2N(KL1) = ZDIAMN(KL1)
    !
  ENDIF
  !  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LAGREG',1,ZHOOK_HANDLE)
!
!       3. Update snow grains parameters : GRAN1, GRAN2        
!        PSNOWGRAN1(KL1)=ZSNOWGRAN1(KL1)
!        PSNOWGRAN2(KL1)=ZSNOWGRAN2(KL1)
!
CONTAINS
!
SUBROUTINE GET_AGREG(KID1,KID2,PFIELD1,PFIELD2,PFIELD)
!
INTEGER, INTENT(IN) :: KID1, KID2
REAL, INTENT(IN) :: PFIELD1, PFIELD2
REAL, INTENT(OUT) :: PFIELD
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LAGREG:GET_AGREG',0,ZHOOK_HANDLE)
!
PFIELD = ( PFIELD1 * PSNOWRHO(KID1) * ( PSNOWDZN(KID1) - ABS(PSNOWDDZ(IL)) ) &
         + PFIELD2 * PSNOWRHO(KID2) * ABS(PSNOWDDZ(IL))                        ) / &
         (           PSNOWRHO (KL1) * ( PSNOWDZN (KL1) - ABS(PSNOWDDZ(IL)) ) + &
                     PSNOWRHO (KL2) * ABS(PSNOWDDZ(IL))                        ) 
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LAGREG:GET_AGREG',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_AGREG
!
END SUBROUTINE SNOW3LAGREG    
!###############################################################################    
!###############################################################################
!
!       
!ajout EB : ajout des arguments "N" pour faire idem variables d'origine
SUBROUTINE SNOW3LAVGRAIN(PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,                &
                         PSNOWGRAN1N,PSNOWGRAN2N,PSNOWHISTN,PNDENT,PNVIEU,&
                         HSNOWMETAMO) 
!
USE MODD_SNOW_METAMO, ONLY : XVDIAM6, XUEPSI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!       0.1 declarations of arguments        
!        
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST 
! ajout EB
REAL, DIMENSION(:,:), INTENT(INOUT)  :: PSNOWGRAN1N,PSNOWGRAN2N,PSNOWHISTN 
!
REAL, DIMENSION(:), INTENT(IN)    :: PNDENT, PNVIEU          
!
CHARACTER(3), INTENT(IN)              :: HSNOWMETAMO
!       0.2 declaration of local variables
!
REAL, DIMENSION(SIZE(PSNOWGRAN1,1)) :: ZGRAN1, ZGRAN2, ZHIST 
!
LOGICAL, DIMENSION(SIZE(PSNOWGRAN1,1),SIZE(PSNOWGRAN1,2)) :: GDENDRITIC
!
INTEGER :: JI, JL
INTEGER :: INLVLS, INI
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!         
!       0.3 initialization         
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LAVGRAIN',0,ZHOOK_HANDLE)
!
INLVLS    = SIZE(PSNOWGRAN1,2) 
INI       = SIZE(PSNOWGRAN1,1)
!
ZGRAN1(:) = 0.0
ZGRAN2(:) = 0.0
ZHIST (:) = 0.0
!     
DO JI = 1,INI
  !
  IF ( PNDENT(JI)==0.0 .AND. PNVIEU(JI)==0.0 ) THEN
    !
    ZGRAN1(JI) = 1.0
    ZGRAN2(JI) = 1.0 
    ZHIST (JI) = 1.0
    !
  ELSE
    !
    DO JL = 1,INLVLS
      IF ( HSNOWMETAMO=='B92' ) THEN
        GDENDRITIC(JI,JL) = ( PSNOWGRAN1(JI,JL) < 0.0 )
      ELSE
        GDENDRITIC(JI,JL) = ( PSNOWGRAN1(JI,JL) < XVDIAM6*(4.-PSNOWGRAN2(JI,JL)) - XUEPSI )
      ENDIF
    ENDDO
    !
    IF ( PNDENT(JI)>=PNVIEU(JI) ) THEN      ! more dendritic than non dendritic snow layer
      !
      DO JL = 1,INLVLS
        IF ( GDENDRITIC(JI,JL) ) THEN
          ZGRAN1(JI) = ZGRAN1(JI) + PSNOWGRAN1(JI,JL)
          ZGRAN2(JI) = ZGRAN2(JI) + PSNOWGRAN2(JI,JL)
        ENDIF
      ENDDO
      !
      PSNOWGRAN1N(JI,:) = ZGRAN1(JI) / PNDENT(JI)   
      PSNOWGRAN2N(JI,:) = ZGRAN2(JI) / PNDENT(JI)
      PSNOWHISTN (JI,:) = 0.0
      !
    ELSE                              ! more non dendritic than dendritic snow layers  
      !
      DO JL = 1,INLVLS
        IF ( .NOT.GDENDRITIC(JI,JL) ) THEN
          ZGRAN1(JI) = ZGRAN1(JI) + PSNOWGRAN1(JI,JL)
          ZGRAN2(JI) = ZGRAN2(JI) + PSNOWGRAN2(JI,JL)
          ZHIST (JI) = ZHIST (JI) + PSNOWHIST (JI,JL) 
        ENDIF
      ENDDO
      !
      PSNOWGRAN1N(JI,:) = ZGRAN1(JI) / PNVIEU(JI)
      PSNOWGRAN2N(JI,:) = ZGRAN2(JI) / PNVIEU(JI)
      PSNOWHISTN (JI,:) = ZHIST (JI) / PNVIEU(JI)
      !    
    ENDIF
    !
  ENDIF
  !
ENDDO


 
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LAVGRAIN',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW3LAVGRAIN         
!        
!####################################################################
!####################################################################
!####################################################################
FUNCTION SNOW3LDIFTYP(PGRAIN1,PGRAIN2,PGRAIN3,PGRAIN4,HSNOWMETAMO,&
                      PSNOWRHO1,PSNOWRHO2,PSNOWAGE1,PSNOWAGE2) RESULT(ZDIFTYPE)
!
! à remplacer sans doute par une routine equivalente du nouveau crocus
!*    CALCUL DE LA DIFFERENCE ENTRE DEUX TYPES DE GRAINS
!     VALEUR ENTRE 200 ET 0
!
USE MODD_SNOW_METAMO, ONLY : XGRAN, XVDIAM6, XUEPSI
USE MODD_SNOW_PAR, ONLY : XRHOTHRESHOLD_ICE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!*      0.1    declarations of arguments
REAL, INTENT(IN) :: PGRAIN1, PGRAIN2, PGRAIN3, PGRAIN4, PSNOWRHO1, PSNOWRHO2, PSNOWAGE1, PSNOWAGE2
CHARACTER(3), INTENT(IN)              :: HSNOWMETAMO
REAL :: ZDIFTYPE, ZCOEF3, ZCOEF4
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!*      0.2    calcul de la difference entre type de grains
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LDIFTYP',0,ZHOOK_HANDLE)
!
IF (ABS(PSNOWAGE1 - PSNOWAGE2) > 600.) THEN
  ZDIFTYPE = 200.
ELSEIF (((PSNOWRHO1 < XRHOTHRESHOLD_ICE).AND.(PSNOWRHO2 >= XRHOTHRESHOLD_ICE)) &
   .OR. ((PSNOWRHO1 >= XRHOTHRESHOLD_ICE).AND.(PSNOWRHO2 < XRHOTHRESHOLD_ICE))) THEN
  ZDIFTYPE = 200.
ELSEIF ( HSNOWMETAMO=='B92' ) THEN 
  !
  IF ( ( PGRAIN1<0.  .AND. PGRAIN2>=0.) .OR. ( PGRAIN1>=0. .AND. PGRAIN2<0. ) ) THEN
    ZDIFTYPE = 200.
  ELSEIF ( PGRAIN1<0. ) THEN
    ZDIFTYPE = ABS( PGRAIN1-PGRAIN2 ) * .5 + ABS( PGRAIN3-PGRAIN4 ) * .5
  ELSE
    ZDIFTYPE = ABS( PGRAIN1-PGRAIN2 )      + ABS( PGRAIN3-PGRAIN4 ) * 5. * 10000.
  ENDIF
  !
ELSE
  !
  ZCOEF3 = XVDIAM6 * (4.-PGRAIN3) - XUEPSI
  ZCOEF4 = XVDIAM6 * (4.-PGRAIN4) - XUEPSI 
  IF ( ( PGRAIN1<ZCOEF3 .AND. PGRAIN2>=ZCOEF4 ) .OR. ( PGRAIN1>=ZCOEF3 .AND. PGRAIN2<ZCOEF4 ) ) THEN
    ZDIFTYPE = 200.
  ELSEIF ( PGRAIN1<ZCOEF3 ) THEN
    ZDIFTYPE = ABS( (PGRAIN3-PGRAIN4)*XGRAN ) * .5 + &
               ABS( ( (PGRAIN1/XVDIAM6 - 4. + PGRAIN3) / (PGRAIN3 - 3.) - &
                      (PGRAIN2/XVDIAM6 - 4. + PGRAIN4) / (PGRAIN4 - 3.) ) * XGRAN ) * .5
             
  ELSE
    ZDIFTYPE = ABS( (PGRAIN3-PGRAIN4)*XGRAN )      + ABS( ZCOEF3-ZCOEF4 ) * 5. * 10000.
  ENDIF  
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LDIFTYP',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LDIFTYP

!####################################################################
SUBROUTINE GET_MASS_HEAT(KJ,KNLVLS_NEW,KNLVLS_OLD,                                &
                         PSNOWZTOP_OLD,PSNOWZTOP_NEW,PSNOWZBOT_OLD,PSNOWZBOT_NEW, &
                         PSNOWRHOO,PSNOWDZO,PSNOWGRAN1O,PSNOWGRAN2O,PSNOWHISTO,   &
                         PSNOWAGEO,PSNOWIMPURO,PSNOWHEATO,                        &
                         PSNOWRHON,PSNOWDZN,PSNOWGRAN1N,PSNOWGRAN2N,PSNOWHISTN,   &
                         PSNOWAGEN, PSNOWIMPURN,PSNOWHEATN,HSNOWMETAMO            )
!
USE MODD_SNOW_PAR, ONLY : XSNOWCRITD, XD1, XD2, XD3, XX, XVALB5, XVALB6
!
USE MODD_SNOW_METAMO, ONLY : XUEPSI
USE MODD_PREP_SNOW,   ONLY : NIMPUR
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER, INTENt(IN) :: KJ
INTEGER, INTENT(IN) :: KNLVLS_NEW, KNLVLS_OLD
REAL, DIMENSION(:), INTENT(IN) :: PSNOWZTOP_OLD, PSNOWZBOT_OLD
REAL, DIMENSION(:), INTENT(IN) :: PSNOWZTOP_NEW, PSNOWZBOT_NEW
REAL, DIMENSION(:), INTENT(IN) :: PSNOWRHOO, PSNOWDZO, PSNOWGRAN1O, PSNOWGRAN2O, &
                                  PSNOWHISTO, PSNOWAGEO, PSNOWHEATO
REAL, DIMENSION(:,:), INTENT(IN) :: PSNOWIMPURO
REAL, DIMENSION(:), INTENT(IN) :: PSNOWDZN
CHARACTER(3), INTENT(IN)       :: HSNOWMETAMO
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWRHON, PSNOWGRAN1N, PSNOWGRAN2N, &
                                   PSNOWHISTN, PSNOWAGEN, PSNOWHEATN
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWIMPURN
!
REAL :: ZPROPOR, ZMASDZ_OLD, ZDIAM, ZMASTOT_T07
REAL :: ZSNOWHEAN, ZMASTOTN, ZDENTMOYN, ZSPHERMOYN, ZALBMOYN, ZHISTMOYN
REAL :: ZAGEMOYN
!
REAL, DIMENSION(NIMPUR) :: ZIMPURMOYN
!
INTEGER :: JST_NEW, JST_OLD, JIMP
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GET_MASS_HEAT',0,ZHOOK_HANDLE)
!
PSNOWRHON  (:) = 0.
PSNOWGRAN1N(:) = 0.
PSNOWGRAN2N(:) = 0.
PSNOWHISTN (:) = 0.
PSNOWAGEN  (:) = 0.
DO JIMP=1,NIMPUR
  PSNOWIMPURN(:,JIMP) = 0.
ENDDO
PSNOWHEATN (:) = 0.
!
DO JST_NEW = 1,KNLVLS_NEW
  !
  ZSNOWHEAN   = 0.
  ZMASTOTN    = 0.
  ZMASTOT_T07 = 0.
  ZDENTMOYN   = 0.
  ZSPHERMOYN  = 0.
  ZALBMOYN    = 0.
  ZDIAM       = 0.
  ZHISTMOYN   = 0.
  ZAGEMOYN    = 0.
  DO JIMP=1,NIMPUR
    ZIMPURMOYN(JIMP) =0.
  ENDDO
  !
  ! lopp over the ols snow layers 
  DO JST_OLD = 1,KNLVLS_OLD
    !
    IF( PSNOWZTOP_OLD(JST_OLD)<=PSNOWZBOT_NEW(JST_NEW) ) THEN
      ! JST_OLD lower than JJ_NEW ==> no contribution
    ELSEIF ( PSNOWZBOT_OLD(JST_OLD)>=PSNOWZTOP_NEW(JST_NEW) ) THEN
      ! JST_OLD higher than JJ_NEW ==> no contribution
    ELSE
      ! old layer contributes to the new one
      ! computation of its contributing ratio and mass/heat 
      !
      !        NEW                                           OLD
      !                                              ------------------- PSNOWZTOP_OLD(JST_OLD)
      !------------------- PSNOWZTOP_NEW(JST_NEW)                                                     |
      !                                                                                               | ZPROPOR
      !                                              ------------------- PSNOWZTOP_OLD(JST_OLD)       |
      !------------------- PSNOWZTOP_NEW(JST_NEW)
      !
      ! The ratio is done in term of layer thickness but it is equivalent to perform this ratio in term of mass (SWE).
      ! Indeed the ratio only concern one old layer at each loop step, and as its SWE is constant it is equivalent to 
      ! think in term of depth or in term os mass (only a constant factor RhoOld between the two approach).
      ZPROPOR = ( MIN( PSNOWZTOP_OLD(JST_OLD), PSNOWZTOP_NEW(JST_NEW) )   &
                - MAX( PSNOWZBOT_OLD(JST_OLD), PSNOWZBOT_NEW(JST_NEW) ) ) &
                 / PSNOWDZO(JST_OLD) 
      ZMASDZ_OLD = ZPROPOR * PSNOWRHOO(JST_OLD) * PSNOWDZO(JST_OLD)
      !
      ! The mass of new snow is incremented with the different old layers contributing
      ZMASTOTN    = ZMASTOTN + ZMASDZ_OLD
      ZMASTOT_T07 = ZMASTOT_T07 + 1.
      !
      ZSNOWHEAN = ZSNOWHEAN + ZPROPOR * PSNOWHEATO(JST_OLD)
      !
      IF ( HSNOWMETAMO=='B92' ) THEN
        !
        ! contribution to the grain optical size and then to the albedo
        IF ( PSNOWGRAN1O(JST_OLD)<0. ) THEN
          ZDIAM = -PSNOWGRAN1O(JST_OLD)*XD1/XX + (1.+PSNOWGRAN1O(JST_OLD)/XX) * &
                 ( PSNOWGRAN2O(JST_OLD)*XD2/XX + (1.-PSNOWGRAN2O(JST_OLD)/XX)*XD3 ) 
          ZDIAM = ZDIAM/10000.      
          ZDENTMOYN  = ZDENTMOYN  - ZMASDZ_OLD * PSNOWGRAN1O(JST_OLD) / XX
          ZSPHERMOYN = ZSPHERMOYN + ZMASDZ_OLD * PSNOWGRAN2O(JST_OLD) / XX
        ELSE
          ZDIAM = PSNOWGRAN2O(JST_OLD)
          ZDENTMOYN  = ZDENTMOYN  + ZMASDZ_OLD * 0.
          ZSPHERMOYN = ZSPHERMOYN + ZMASDZ_OLD * PSNOWGRAN1O(JST_OLD) / XX
        ENDIF
        !
      ELSE
        !
        ZDIAM = PSNOWGRAN1O(JST_OLD)
        ZSPHERMOYN = ZSPHERMOYN + ZMASDZ_OLD * PSNOWGRAN2O(JST_OLD)
        !
      ENDIF
      !
      ZALBMOYN  = ZALBMOYN  + MAX( 0., ZMASDZ_OLD * (XVALB5-XVALB6*SQRT(ZDIAM)) )
      ZHISTMOYN = ZHISTMOYN + ZMASDZ_OLD * PSNOWHISTO(JST_OLD)
      ZAGEMOYN  = ZAGEMOYN  + ZMASDZ_OLD * PSNOWAGEO (JST_OLD)
      ! In fact Zimpurmoyen is adding the contibution in impurity content of all old layers that are used to build a new one.
      !Zpropor is the 
      IF (ZPROPOR>XUEPSI) THEN
          DO JIMP=1,NIMPUR
            ZIMPURMOYN(JIMP)=ZIMPURMOYN(JIMP)+ PSNOWIMPURO(JST_OLD,JIMP) * ZPROPOR
          ENDDO
      ENDIF
      !
    ENDIF
    !
  ENDDO
  ! 
  ! the new layer inherits from the weihted average properties of the old ones
  ! heat and mass
  PSNOWHEATN(JST_NEW) = ZSNOWHEAN
  PSNOWRHON (JST_NEW) = ZMASTOTN / PSNOWDZN(JST_NEW)
  ! grain type and size decuced from the average albedo
  ZALBMOYN   = ZALBMOYN / ZMASTOTN
  ZSPHERMOYN = MAX( 0., ZSPHERMOYN/ZMASTOTN )
  ZDENTMOYN  = MAX( 0., ZDENTMOYN /ZMASTOTN )
  ZDIAM = ( (XVALB5-ZALBMOYN)/XVALB6 )**2
  !
  IF ( HSNOWMETAMO=='B92' ) THEN
    !
    ! size between D2 and D3 and dendricity < 0         
    ! sphericity is firts preserved, if possible. If not,
    ! denditricity =0
    PSNOWGRAN1N(JST_NEW) = -XX * ZDENTMOYN
    !
    IF ( ZDENTMOYN/=1.) THEN
      PSNOWGRAN2N(JST_NEW) = XX * ( ( ZDIAM*10000. + PSNOWGRAN1N(JST_NEW)*XD1/XX ) &
                                 / ( 1. + PSNOWGRAN1N(JST_NEW)/XX ) - XD3 )        &
                             / ( XD2-XD3 )
    ENDIF
    !
    ! dendricity is preserved if possible and sphericity is adjusted
    IF ( ZDIAM < XD2/10000. - 0.0000001 ) THEN
      !
      IF ( ABS( PSNOWGRAN1N(JST_NEW)+XX ) < 0.01 ) THEN
        !
        PSNOWGRAN2N(JST_NEW) = XX * ZSPHERMOYN
        !
      ELSEIF ( ABS( PSNOWGRAN1N(JST_NEW) ) < 0.0001 ) THEN ! dendritic snow
        !
        PSNOWGRAN1N(JST_NEW) = XX * ZSPHERMOYN
        PSNOWGRAN2N(JST_NEW) = ZDIAM
        !
      ELSEIF ( PSNOWGRAN2N(JST_NEW) < 0. ) THEN ! non dendritic
        !
        PSNOWGRAN2N(JST_NEW) = 0.
        !
      ELSEIF ( PSNOWGRAN2N(JST_NEW) > XX + 0.0000001 ) THEN ! non dendritic
        !
        PSNOWGRAN2N(JST_NEW) = XX
        !
      ENDIF
      !
    ELSEIF ( ZDIAM > XD3/10000. .OR. ZDENTMOYN <= 0. + 0.0000001 .OR. &
             PSNOWGRAN2N(JST_NEW) < 0. .OR. PSNOWGRAN2N(JST_NEW) > XX ) THEN
      !
      ! dendritic snow
      ! inconsistency with ZDIAM ==>  dendricity = 0
      ! size between D2 and D3 and dendricity == 0          
      PSNOWGRAN1N(JST_NEW) = XX * ZSPHERMOYN
      PSNOWGRAN2N(JST_NEW) = ZDIAM
      !
    ENDIF
    !
  ELSE
    !
    PSNOWGRAN1N(JST_NEW) = ZDIAM
    PSNOWGRAN2N(JST_NEW) = MIN( 1., ZSPHERMOYN )
    !
  ENDIF
  !
  PSNOWHISTN(JST_NEW) = NINT( ZHISTMOYN/ZMASTOTN )
  PSNOWAGEN (JST_NEW) = ZAGEMOYN / ZMASTOTN
  DO JIMP=1,NIMPUR
    PSNOWIMPURN(JST_NEW,JIMP)= ZIMPURMOYN(JIMP)
  ENDDO
  !
ENDDO   
!print *, '~PSNOWHEATO ', PSNOWHEATO(1:6)  ! MN
!print *, '~~PSNOWHEATN ', PSNOWHEATN(1:4)

!
IF (LHOOK) CALL DR_HOOK('GET_MASS_HEAT',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_MASS_HEAT
!####################################################################
SUBROUTINE GET_DIAM(PSNOWGRAN1,PSNOWGRAN2,PDIAM,HSNOWMETAMO)
!
USE MODD_SNOW_PAR, ONLY : XD1, XD2, XD3, XX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PSNOWGRAN1
REAL, INTENT(IN) :: PSNOWGRAN2
REAL, INTENT(OUT) :: PDIAM
!
CHARACTER(3), INTENT(IN)              :: HSNOWMETAMO
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GET_DIAM',0,ZHOOK_HANDLE)
!
IF ( HSNOWMETAMO=='B92' ) THEN
  !
  IF( PSNOWGRAN1<0. ) THEN
    PDIAM = -PSNOWGRAN1*XD1/XX + (1.+PSNOWGRAN1/XX) * &
           ( PSNOWGRAN2*XD2/XX + (1.-PSNOWGRAN2/XX) * XD3 ) 
    PDIAM = PDIAM/10000.      
  ELSE 
    PDIAM = PSNOWGRAN2*PSNOWGRAN1/XX + &
            MAX( 0.0004, 0.5*PSNOWGRAN2 ) * ( 1.-PSNOWGRAN1/XX )      
  ENDIF
  !
ELSE
  !
  PDIAM = PSNOWGRAN1
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('GET_DIAM',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_DIAM
!####################################################################
SUBROUTINE SYVAGRE(PSNOWGRAN1I,PSNOWGRAN2I,PSNOWGRAN1J,PSNOWGRAN2J,      &
                   PSNOWGRAN1F,PSNOWGRAN2F,PZI,PZJ)
!
!!    PURPOSE
!!    -------
!!    Aggregate snow grain characteristics in 2 layers (I and J) to get the
!!    averaged grain characteristics 
!!
!!    METHOD
!!    -------
!!    Based on the subroutine get_mass_heat in mode_snow3l.f90
!!
!!    AUTHOR
!!    ------
!!    V. Vionnet  * Meteo-France *   Implementation in SURFEX
!!
!
USE MODD_SNOW_PAR, ONLY : XSNOWCRITD, XD1, XD2, XD3, XX, XVALB5, XVALB6
!   
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)      :: PSNOWGRAN1I,PSNOWGRAN2I,PSNOWGRAN1J,PSNOWGRAN2J
REAL, INTENT(IN)      :: PZI,PZJ      !   Weight of layer I and J
REAL, INTENT(OUT)     :: PSNOWGRAN1F,PSNOWGRAN2F

!
!*      0.2    declarations of local variables
!
REAL ZDIAMI,ZDIAMJ,ZDIAMMOY
REAL ZDENTI,ZDENTJ,ZDENTMOY
REAL ZSPHERI,ZSPHERJ,ZSPHERMOY
REAL ZALBI,ZALBJ,ZALBMOY

!
!    1. Compute properties for each layer
!
! For layer I
IF(PSNOWGRAN1I<0.) THEN
    ZDIAMI = -PSNOWGRAN1I*XD1/XX + (1.+PSNOWGRAN1I/XX) * &
                 ( PSNOWGRAN2I*XD2/XX +(1.-PSNOWGRAN2I/XX)*XD3 ) 
    ZDIAMI = ZDIAMI/10000.
    ZDENTI  = -PSNOWGRAN1I/XX
    ZSPHERI = PSNOWGRAN2I/XX
ELSE
   ZDIAMI = PSNOWGRAN2I
   ZDENTI =0.
   ZSPHERI = PSNOWGRAN1I/XX
ENDIF
ZALBI = MAX( 0.,  (XVALB5-XVALB6*SQRT(ZDIAMI)) )
!
! For layer J
IF(PSNOWGRAN1J<0.) THEN
    ZDIAMJ = -PSNOWGRAN1J*XD1/XX + (1.+PSNOWGRAN1J/XX) * &
                 ( PSNOWGRAN2J*XD2/XX +(1.-PSNOWGRAN2J/XX)*XD3 )
    ZDIAMJ = ZDIAMJ/10000.
    ZDENTJ  = -PSNOWGRAN1J/XX
    ZSPHERJ = PSNOWGRAN2J/XX
ELSE
   ZDIAMJ = PSNOWGRAN2J
   ZDENTJ =0.
   ZSPHERJ = PSNOWGRAN1J/XX
ENDIF
ZALBJ = MAX( 0.,  (XVALB5-XVALB6*SQRT(ZDIAMJ)) )
!
!    2. Compute averaged properties
!
ZDENTMOY  = MAX(0.,(ZDENTI*PZI+ZDENTJ*PZJ)/( PZI+PZJ))
ZSPHERMOY = MAX(0.,(ZSPHERI*PZI+ZSPHERJ*PZJ)/( PZI+PZJ))
ZALBMOY   = MAX(0.,(ZALBI*PZI+ZALBJ*PZJ)/( PZI+PZJ))
ZDIAMMOY = ( (XVALB5-ZALBMOY)/XVALB6 )**2
!
!    3. Compute GRAN1 and GRAN2 of averaged layer
!
! size between D2 and D3 and dendricity < 0         
! sphericity is first preserved, if possible. If not,
! dendricity = 0
PSNOWGRAN1F= -XX * ZDENTMOY
!
IF(ZDENTMOY/=1.) THEN
      PSNOWGRAN2F = XX * ( ( ZDIAMMOY*10000. + PSNOWGRAN1F*XD1/XX) &
                     / ( 1. + PSNOWGRAN1F/XX ) - XD3 )/ ( XD2-XD3 )
ENDIF
!
! dendricity is preserved if possible and sphericity is adjusted
IF ( ZDIAMMOY < XD2/10000. - 0.0000001 ) THEN
    !
    IF ( ABS( PSNOWGRAN1F+XX ) < 0.01 ) THEN
    !
       PSNOWGRAN2F = XX * ZSPHERMOY
    !
    ELSEIF ( ABS( PSNOWGRAN1F) < 0.0001 ) THEN ! dendritic snow
    !
      PSNOWGRAN1F = XX * ZSPHERMOY
      PSNOWGRAN2F = ZDIAMMOY
!
    ELSEIF ( PSNOWGRAN2F < 0. ) THEN ! non dendritic
!
      PSNOWGRAN2F = 0.
!
    ELSEIF ( PSNOWGRAN2F > XX + 0.0000001 ) THEN ! non dendritic
!
      PSNOWGRAN2F = XX
!
    ENDIF
!
ELSEIF ( ZDIAMMOY > XD3/10000. .OR. ZDENTMOY <= 0. + 0.0000001 .OR. &
     PSNOWGRAN2F < 0. .OR. PSNOWGRAN2F > XX ) THEN
!
! dendritic snow
! inconsistency with ZDIAM ==>  dendricity = 0
! size between D2 and D3 and dendricity == 0          
    PSNOWGRAN1F = XX * ZSPHERMOY
    PSNOWGRAN2F = ZDIAMMOY
!
ENDIF

END SUBROUTINE SYVAGRE
!####################################################################
!####################################################################
FUNCTION SNOW3LRADABS_0D(PSNOWRHO,PSNOWDZ,PSPECTRALALBEDO,PZENITH,PPERMSNOWFRAC,PDSGRAIN) RESULT(PCOEF)
!
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave radiation within the snowpack
!     (with depth)
!     A. Boone 02/2011
!     A. Boone 11/2014 Updated to use spectral dependence.
!                      NOTE, assumes 3 spectral bands
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_MEB_PAR,  ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
USE MODD_SNOW_PAR, ONLY : XVSPEC1,XVSPEC2,XVSPEC3,XVBETA1,XVBETA2, &
                          XVBETA4,XVBETA3,XVBETA5, XMINCOSZEN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN)                   :: PSNOWRHO        ! snow density    (kg m-3)
REAL,               INTENT(IN)                   :: PSNOWDZ         ! layer thickness (m)
REAL,               INTENT(IN)                   :: PZENITH         ! zenith angle    (rad)
REAL,               INTENT(IN)                   :: PPERMSNOWFRAC   ! permanent snow fraction (-)
REAL, DIMENSION(:), INTENT(IN)                   :: PSPECTRALALBEDO ! spectral albedo (-)
REAL,               INTENT(IN)                   :: PDSGRAIN        ! Snow optical grain diameter (m)
!
REAL                                             :: PCOEF           ! -
!
!*      0.2    declarations of local variables
!
REAL                                             :: ZWORK, ZPROJLAT,                  &
                                                    ZBETA1, ZBETA2, ZBETA3,           &
                                                    ZOPTICALPATH1, ZOPTICALPATH2,     &
                                                    ZOPTICALPATH3
!
REAL(KIND=JPRB)                                  :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LRADABS_0D',0,ZHOOK_HANDLE)
!
! Coefficient for taking into account the increase of path length of rays
! in snow due to zenithal angle
!
ZPROJLAT         = (1.0-PPERMSNOWFRAC)+PPERMSNOWFRAC/ &
                   MAX(XMINCOSZEN,COS(PZENITH))
!
! Extinction coefficient:
!
ZWORK            = SQRT(PDSGRAIN)
ZBETA1           = MAX(XVBETA1*PSNOWRHO/ZWORK,XVBETA2)
ZBETA2           = MAX(XVBETA3*PSNOWRHO/ZWORK,XVBETA4)
ZBETA3           = XVBETA5
!
ZOPTICALPATH1    = ZBETA1*PSNOWDZ
ZOPTICALPATH2    = ZBETA2*PSNOWDZ
ZOPTICALPATH3    = XUNDEF
!
IF(PSPECTRALALBEDO(3)==XUNDEF)THEN 
   PCOEF         = XSW_WGHT_VIS*(1.0-PSPECTRALALBEDO(1))*EXP(-ZOPTICALPATH1*ZPROJLAT) &
                 + XSW_WGHT_NIR*(1.0-PSPECTRALALBEDO(2))*EXP(-ZOPTICALPATH2*ZPROJLAT) 
ELSE
   ZOPTICALPATH3 = ZBETA3*PSNOWDZ
   PCOEF         = XVSPEC1*(1.0-PSPECTRALALBEDO(1))*EXP(-ZOPTICALPATH1*ZPROJLAT) &
                 + XVSPEC2*(1.0-PSPECTRALALBEDO(2))*EXP(-ZOPTICALPATH2*ZPROJLAT) &
                 + XVSPEC3*(1.0-PSPECTRALALBEDO(3))*EXP(-ZOPTICALPATH3*ZPROJLAT)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LRADABS_0D',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION SNOW3LRADABS_0D
!####################################################################
!####################################################################
!####################################################################
FUNCTION SNOW3LRADABS_1D(PSNOWRHO,PSNOWDZ,PSPECTRALALBEDO,PZENITH,PPERMSNOWFRAC,PDSGRAIN) RESULT(PCOEF)
!
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave radiation within the snowpack
!     (with depth)
!     A. Boone 02/2011
!     A. Boone 11/2014 Updated to use spectral dependence
!                      NOTE, assumes 3 spectral bands
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_MEB_PAR,  ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
USE MODD_SNOW_PAR, ONLY : XVSPEC1,XVSPEC2,XVSPEC3,XVBETA1,XVBETA2, &
                          XVBETA4,XVBETA3,XVBETA5, XMINCOSZEN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:),   INTENT(IN)                 :: PSNOWRHO        ! snow density    (kg m-3)
REAL, DIMENSION(:),   INTENT(IN)                 :: PSNOWDZ         ! layer thickness (m)
REAL, DIMENSION(:),   INTENT(IN)                 :: PZENITH         ! zenith angle    (rad)
REAL, DIMENSION(:),   INTENT(IN)                 :: PPERMSNOWFRAC   ! permanent snow fraction (-)
REAL, DIMENSION(:,:), INTENT(IN)                 :: PSPECTRALALBEDO ! spectral albedo (-)
REAL, DIMENSION(:),   INTENT(IN)                 :: PDSGRAIN        ! Snow optical grain diameter (m)
!
REAL, DIMENSION(SIZE(PSNOWRHO))                  :: PCOEF           ! -
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO))                  :: ZWORK, ZPROJLAT,                  &
                                                    ZBETA1, ZBETA2, ZBETA3,           &
                                                    ZOPTICALPATH1, ZOPTICALPATH2,     &
                                                    ZOPTICALPATH3
!
REAL(KIND=JPRB)                                  :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LRADABS_1D',0,ZHOOK_HANDLE)
!
! Coefficient for taking into account the increase of path length of rays
! in snow due to zenithal angle
!
ZPROJLAT(:)         = (1.0-PPERMSNOWFRAC(:))+PPERMSNOWFRAC(:)/ &
                      MAX(XMINCOSZEN,COS(PZENITH(:)))
!
! Extinction coefficient:
!
ZWORK(:)            = SQRT(PDSGRAIN(:))
ZBETA1(:)           = MAX(XVBETA1*PSNOWRHO(:)/ZWORK(:),XVBETA2)
ZBETA2(:)           = MAX(XVBETA3*PSNOWRHO(:)/ZWORK(:),XVBETA4)
ZBETA3(:)           = XVBETA5
!
ZOPTICALPATH1(:)    = ZBETA1(:)*PSNOWDZ(:)
ZOPTICALPATH2(:)    = ZBETA2(:)*PSNOWDZ(:)
ZOPTICALPATH3(:)    = XUNDEF
!
WHERE(PSPECTRALALBEDO(:,3)==XUNDEF)
   PCOEF(:)         = XSW_WGHT_VIS*(1.0-PSPECTRALALBEDO(:,1))*EXP(-ZOPTICALPATH1(:)*ZPROJLAT(:)) &
                    + XSW_WGHT_NIR*(1.0-PSPECTRALALBEDO(:,2))*EXP(-ZOPTICALPATH2(:)*ZPROJLAT(:)) 
ELSEWHERE
   ZOPTICALPATH3(:) = ZBETA3(:)*PSNOWDZ(:)
   PCOEF(:)         = XVSPEC1*(1.0-PSPECTRALALBEDO(:,1))*EXP(-ZOPTICALPATH1(:)*ZPROJLAT(:)) &
                    + XVSPEC2*(1.0-PSPECTRALALBEDO(:,2))*EXP(-ZOPTICALPATH2(:)*ZPROJLAT(:)) &
                    + XVSPEC3*(1.0-PSPECTRALALBEDO(:,3))*EXP(-ZOPTICALPATH3(:)*ZPROJLAT(:))
END WHERE
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LRADABS_1D',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION SNOW3LRADABS_1D
!####################################################################
!####################################################################
!####################################################################
FUNCTION SNOW3LRADABS_2D(PSNOWRHO,PSNOWDZ,PSPECTRALALBEDO,PZENITH,PPERMSNOWFRAC,PDSGRAIN) RESULT(PCOEF)
!
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave radiation within the snowpack
!     (with depth)
!     A. Boone 02/2011
!     A. Boone 11/2014 Updated to use spectral dependence
!                      NOTE, assumes 3 spectral bands
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_MEB_PAR,  ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
USE MODD_SNOW_PAR, ONLY : XVSPEC1,XVSPEC2,XVSPEC3,XVBETA1,XVBETA2, &
                          XVBETA4,XVBETA3,XVBETA5, XMINCOSZEN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:),   INTENT(IN)                 :: PSNOWRHO        ! snow density    (kg m-3)
REAL, DIMENSION(:,:),   INTENT(IN)                 :: PSNOWDZ         ! layer thickness (m)
REAL, DIMENSION(:,:),   INTENT(IN)                 :: PZENITH         ! zenith angle    (rad)
REAL, DIMENSION(:,:),   INTENT(IN)                 :: PPERMSNOWFRAC   ! permanent snow fraction (-)
REAL, DIMENSION(:,:,:), INTENT(IN)                 :: PSPECTRALALBEDO ! spectral albedo (-)
REAL, DIMENSION(:,:),   INTENT(IN)                 :: PDSGRAIN        ! Snow optical grain diameter (m)
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: PCOEF           ! -
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZWORK, ZPROJLAT,                  &
                                                      ZBETA1, ZBETA2, ZBETA3,           &
                                                      ZOPTICALPATH1, ZOPTICALPATH2,     &
                                                      ZOPTICALPATH3
!
REAL(KIND=JPRB)                                    :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LRADABS_2D',0,ZHOOK_HANDLE)
!
! Coefficient for taking into account the increase of path length of rays
! in snow due to zenithal angle
!
ZPROJLAT(:,:)         = (1.0-PPERMSNOWFRAC(:,:))+PPERMSNOWFRAC(:,:)/ &
                        MAX(XMINCOSZEN,COS(PZENITH(:,:)))
!
! Extinction coefficient:
!
ZWORK(:,:)            = SQRT(PDSGRAIN(:,:))
ZBETA1(:,:)           = MAX(XVBETA1*PSNOWRHO(:,:)/ZWORK(:,:),XVBETA2)
ZBETA2(:,:)           = MAX(XVBETA3*PSNOWRHO(:,:)/ZWORK(:,:),XVBETA4)
ZBETA3(:,:)           = XVBETA5
!
ZOPTICALPATH1(:,:)    = ZBETA1(:,:)*PSNOWDZ(:,:)
ZOPTICALPATH2(:,:)    = ZBETA2(:,:)*PSNOWDZ(:,:)
ZOPTICALPATH3(:,:)    = XUNDEF
!
WHERE(PSPECTRALALBEDO(:,:,3)==XUNDEF)
   PCOEF(:,:)         = XSW_WGHT_VIS*(1.0-PSPECTRALALBEDO(:,:,1))*EXP(-ZOPTICALPATH1(:,:)*ZPROJLAT(:,:)) &
                      + XSW_WGHT_NIR*(1.0-PSPECTRALALBEDO(:,:,2))*EXP(-ZOPTICALPATH2(:,:)*ZPROJLAT(:,:)) 
ELSEWHERE
   ZOPTICALPATH3(:,:) = ZBETA3(:,:)*PSNOWDZ(:,:)
   PCOEF(:,:)         = XVSPEC1*(1.0-PSPECTRALALBEDO(:,:,1))*EXP(-ZOPTICALPATH1(:,:)*ZPROJLAT(:,:)) &
                      + XVSPEC2*(1.0-PSPECTRALALBEDO(:,:,2))*EXP(-ZOPTICALPATH2(:,:)*ZPROJLAT(:,:)) &
                      + XVSPEC3*(1.0-PSPECTRALALBEDO(:,:,3))*EXP(-ZOPTICALPATH3(:,:)*ZPROJLAT(:,:))
END WHERE
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LRADABS_2D',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION SNOW3LRADABS_2D
!####################################################################
!####################################################################
!####################################################################
FUNCTION SNOW3LRADABS_SFC(PSNOWRHO,PSNOWDZ,PSPECTRALALBEDO,PZENITH,PPERMSNOWFRAC,PDSGRAIN) RESULT(PCOEF)
!
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave radiation within the snowpack
!     (with depth) for 3 spectral bands
!     A. Boone 02/2011
!     A. Boone 11/2014 Updated to use spectral dependence
!                      NOTE, assumes 3 spectral bands
!     A. Boone 06/2017 ONLY considers albedo of uppermost layer  
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_SNOW_PAR, ONLY : XVSPEC1,XVSPEC2,XVSPEC3,XVBETA1,XVBETA2, &
                          XVBETA4,XVBETA3,XVBETA5, XMINCOSZEN
USE MODD_MEB_PAR,  ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:),   INTENT(IN)                 :: PSNOWRHO        ! snow density    (kg m-3)
REAL, DIMENSION(:,:),   INTENT(IN)                 :: PSNOWDZ         ! layer thickness (m)
REAL, DIMENSION(:),     INTENT(IN)                 :: PZENITH         ! zenith angle    (rad)
REAL, DIMENSION(:),     INTENT(IN)                 :: PPERMSNOWFRAC   ! permanent snow fraction (-)
REAL, DIMENSION(:,:),   INTENT(IN)                 :: PSPECTRALALBEDO ! spectral albedo (-)
REAL, DIMENSION(:,:),   INTENT(IN)                 :: PDSGRAIN        ! Snow optical grain diameter (m)
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: PCOEF           ! -
!
!*      0.2    declarations of local variables
!
INTEGER                                            :: JJ, JI, INLVLS, INI
REAL, DIMENSION(SIZE(PSNOWRHO,1))                  :: ZPROJLAT,  ZOPTICALPATH1,         &
                                                      ZOPTICALPATH2, ZOPTICALPATH3
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZWORK, ZBETA1, ZBETA2, ZBETA3
!
REAL(KIND=JPRB)                                    :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LRADABS_SFC',0,ZHOOK_HANDLE)
!
INI    = SIZE(PSNOWDZ(:,:),1)
INLVLS = SIZE(PSNOWDZ(:,:),2)
!
! Coefficient for taking into account the increase of path length of rays
! in snow due to zenithal angle
!
ZPROJLAT(:)           = (1.0-PPERMSNOWFRAC(:))+PPERMSNOWFRAC(:)/ &
                        MAX(XMINCOSZEN,COS(PZENITH(:)))
!
! Extinction coefficient:
!
ZWORK(:,:)            = SQRT(PDSGRAIN(:,:))
ZBETA1(:,:)           = MAX(XVBETA1*PSNOWRHO(:,:)/ZWORK(:,:),XVBETA2)
ZBETA2(:,:)           = MAX(XVBETA3*PSNOWRHO(:,:)/ZWORK(:,:),XVBETA4)
ZBETA3(:,:)           = XVBETA5
!
ZOPTICALPATH1(:)      = 0.0
ZOPTICALPATH2(:)      = 0.0
ZOPTICALPATH3(:)      = 0.0
!
DO JJ=1,INLVLS
   DO JI=1,INI
      ZOPTICALPATH1(JI) = ZOPTICALPATH1(JI) + ZBETA1(JI,JJ)*PSNOWDZ(JI,JJ)
      ZOPTICALPATH2(JI) = ZOPTICALPATH2(JI) + ZBETA2(JI,JJ)*PSNOWDZ(JI,JJ)

      IF(PSPECTRALALBEDO(JI,3)==XUNDEF)THEN 

         PCOEF (JI,JJ)     = XSW_WGHT_VIS*(1.0-PSPECTRALALBEDO(JI,1))*EXP(-ZOPTICALPATH1(JI)*ZPROJLAT(JI)) &
                           + XSW_WGHT_NIR*(1.0-PSPECTRALALBEDO(JI,2))*EXP(-ZOPTICALPATH2(JI)*ZPROJLAT(JI))  
      ELSE
      
         ZOPTICALPATH3(JI) = ZOPTICALPATH3(JI) + ZBETA3(JI,JJ)*PSNOWDZ(JI,JJ)

         PCOEF (JI,JJ)     = XVSPEC1*(1.0-PSPECTRALALBEDO(JI,1))*EXP(-ZOPTICALPATH1(JI)*ZPROJLAT(JI)) &
                           + XVSPEC2*(1.0-PSPECTRALALBEDO(JI,2))*EXP(-ZOPTICALPATH2(JI)*ZPROJLAT(JI)) &
                           + XVSPEC3*(1.0-PSPECTRALALBEDO(JI,3))*EXP(-ZOPTICALPATH3(JI)*ZPROJLAT(JI))

      ENDIF
      
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LRADABS_SFC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END FUNCTION SNOW3LRADABS_SFC
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LTHRM(PSNOWRHO,PSCOND,PSNOWTEMP,PPS)
!
!!    PURPOSE
!!    -------
!     Calculate snow thermal conductivity from
!     Sun et al. 1999, J. of Geophys. Res., 104, 19587-19579 (vapor) 
!     and Yen, 1981, CRREL Rep 81-10 (snow)
!     or Anderson, 1976, NOAA Tech. Rep. NWS 19 (snow).
!
!
USE MODD_CSTS,     ONLY : XP00, XCONDI, XRHOLW
!
USE MODD_SNOW_PAR, ONLY : XVRKZ6, XSNOWTHRMCOND1, &
                          XSNOWTHRMCOND2,         &
                          XSNOWTHRMCOND_AVAP,     &
                          XSNOWTHRMCOND_BVAP,     &
                          XSNOWTHRMCOND_CVAP 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)      :: PPS
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP, PSNOWRHO
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSCOND
!
!
!*      0.2    declarations of local variables
!
INTEGER                              :: JJ, JI
!
INTEGER                              :: INI
INTEGER                              :: INLVLS
!
CHARACTER(LEN=5)                     :: YSNOWCOND !should be in namelist
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LTHRM',0,ZHOOK_HANDLE)
!
INI    = SIZE(PSNOWRHO(:,:),1)
INLVLS = SIZE(PSNOWRHO(:,:),2)
!
! 1. Snow thermal conductivity
! ----------------------------
!
YSNOWCOND='YEN81' !should be in namelist
!
IF(YSNOWCOND=='AND76')THEN
!  Thermal conductivity coefficients from Anderson (1976)
  PSCOND(:,:) = (XSNOWTHRMCOND1 + XSNOWTHRMCOND2*PSNOWRHO(:,:)*PSNOWRHO(:,:))
ELSE
! Thermal conductivity coefficients from Yen (1981)
  PSCOND(:,:) = XCONDI * EXP(XVRKZ6*LOG(PSNOWRHO(:,:)/XRHOLW))
ENDIF
!
! 2. Implicit vapor diffn effects
! -------------------------------
!
DO JJ=1,INLVLS
   DO JI=1,INI
    PSCOND(JI,JJ) = PSCOND(JI,JJ) + MAX(0.0,(XSNOWTHRMCOND_AVAP+(XSNOWTHRMCOND_BVAP/(PSNOWTEMP(JI,JJ) &
                                  + XSNOWTHRMCOND_CVAP)))*(XP00/PPS(JI)))
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LTHRM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LTHRM
!####################################################################
!####################################################################
!####################################################################
FUNCTION SNOW3LDOPT_2D(PSNOWRHO,PSNOWAGE) RESULT(PDOPT)
!
!!    PURPOSE
!!    -------
!     Calculate the optical grain diameter.
!
USE MODD_SNOW_PAR, ONLY : XDSGRAIN_MAX,XSNOW_AGRAIN, & 
                          XSNOW_BGRAIN,XSNOW_CGRAIN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)                   :: PSNOWRHO,PSNOWAGE
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: PDOPT
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZAGE
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSRHO4
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LDOPT_2D',0,ZHOOK_HANDLE)
!
ZAGE(:,:) = MIN(15.,PSNOWAGE(:,:))
!
ZSRHO4(:,:) = PSNOWRHO(:,:)*PSNOWRHO(:,:)*PSNOWRHO(:,:)*PSNOWRHO(:,:)
!
PDOPT(:,:) = MIN(XDSGRAIN_MAX,XSNOW_AGRAIN+XSNOW_BGRAIN*ZSRHO4(:,:)+XSNOW_CGRAIN*ZAGE(:,:))
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LDOPT_2D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LDOPT_2D
!####################################################################
FUNCTION SNOW3LDOPT_1D(PSNOWRHO,PSNOWAGE) RESULT(PDOPT)
!
!!    PURPOSE
!!    -------
!     Calculate the optical grain diameter.
!
USE MODD_SNOW_PAR, ONLY : XDSGRAIN_MAX,XSNOW_AGRAIN, & 
                          XSNOW_BGRAIN,XSNOW_CGRAIN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)  :: PSNOWRHO,PSNOWAGE
!
REAL, DIMENSION(SIZE(PSNOWRHO)) :: PDOPT
REAL, DIMENSION(SIZE(PSNOWRHO)) :: ZAGE
REAL, DIMENSION(SIZE(PSNOWRHO)) :: ZSRHO4
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LDOPT_1D',0,ZHOOK_HANDLE)
!
ZAGE(:) = MIN(15.,PSNOWAGE(:))
!
ZSRHO4(:) = PSNOWRHO(:)*PSNOWRHO(:)*PSNOWRHO(:)*PSNOWRHO(:)
!
PDOPT(:) = MIN(XDSGRAIN_MAX,XSNOW_AGRAIN+XSNOW_BGRAIN*ZSRHO4(:)+XSNOW_CGRAIN*ZAGE(:))
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LDOPT_1D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LDOPT_1D
!####################################################################
FUNCTION SNOW3LDOPT_0D(PSNOWRHO,PSNOWAGE) RESULT(PDOPT)
!
!!    PURPOSE
!!    -------
!     Calculate the optical grain diameter.
!
USE MODD_SNOW_PAR, ONLY : XDSGRAIN_MAX,XSNOW_AGRAIN, & 
                          XSNOW_BGRAIN,XSNOW_CGRAIN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)  :: PSNOWRHO,PSNOWAGE
!
REAL :: PDOPT
REAL :: ZAGE
REAL :: ZSRHO4
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LDOPT_0D',0,ZHOOK_HANDLE)
!
ZAGE = MIN(15.,PSNOWAGE)
!
ZSRHO4 = PSNOWRHO*PSNOWRHO*PSNOWRHO*PSNOWRHO
!
PDOPT = MIN(XDSGRAIN_MAX,XSNOW_AGRAIN+XSNOW_BGRAIN*ZSRHO4+XSNOW_CGRAIN*ZAGE)
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LDOPT_0D',1,ZHOOK_HANDLE)
!
END FUNCTION SNOW3LDOPT_0D
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LALB(PALBEDOSC,PSPECTRALALBEDO,PSNOWRHO,PSNOWAGE,   &
                     PPERMSNOWFRAC,PPS)  
!
!!    PURPOSE
!!    -------
!     Calculate the snow surface albedo. Use the method of
!     CROCUS with 3 spectral albedo depending on snow density 
!     and age
!
!
USE MODD_SNOW_PAR, ONLY : XVAGING_GLACIER, XVAGING_NOGLACIER,     &
                          XVALB2,XVALB3,XVALB4,XVALB5,XVALB6,     &
                          XVALB7,XVALB8,XVALB9,XVALB10,XVALB11,   &
                          XVDIOP1,XVRPRE1,XVRPRE2,XVPRES1,        &
                          XVW1,XVW2,XVSPEC1,XVSPEC2,XVSPEC3
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWRHO
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWAGE
REAL, DIMENSION(:), INTENT(IN)      :: PPERMSNOWFRAC
REAL, DIMENSION(:), INTENT(IN)      :: PPS
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PALBEDOSC
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSPECTRALALBEDO
!
!*      0.2    declarations of local variables
!
REAL, PARAMETER                 :: ZALBNIR1 = 0.3
REAL, PARAMETER                 :: ZALBNIR2 = 0.0
!
REAL, DIMENSION(SIZE(PSNOWRHO)) :: ZVAGING, ZDIAM, ZAGE,  &
                                   ZWORK, ZPRES_EFFECT
!
REAL, DIMENSION(SIZE(PSNOWRHO)) :: ZALB1, ZALB2, ZALB3
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LALB',0,ZHOOK_HANDLE)
!
! 0. Initialize:
! ------------------
!
!Snow age effect parameter for Visible (small over glacier)
ZVAGING(:)=XVAGING_GLACIER*PPERMSNOWFRAC(:) + XVAGING_NOGLACIER*(1.0-PPERMSNOWFRAC(:))
!
!Atm pression effect parameter on albedo
ZPRES_EFFECT(:) = XVALB10*MIN(MAX(PPS(:)/XVPRES1,XVRPRE1),XVRPRE2)
!
! 1. Snow optical grain diameter :
! --------------------------------
!
!Snow optical diameter do not depend on snow age over glacier or polar regions
ZAGE(:) = (1.0-PPERMSNOWFRAC(:))*PSNOWAGE(:)
!
ZDIAM(:) = SNOW3LDOPT(PSNOWRHO(:),ZAGE(:))
!
! 2. spectral albedo over 3 bands :
! ---------------------------------
!
!Snow age effect limited to 1 year
ZAGE(:) = MIN(365.,PSNOWAGE(:))
!
ZWORK(:)=SQRT(ZDIAM(:))
!
! Visible
ZALB1(:)=MIN(XVALB4,XVALB2-XVALB3*ZWORK(:))
ZALB1(:)=MAX(XVALB11,ZALB1(:)-ZPRES_EFFECT(:)*ZAGE(:)/ZVAGING(:))
!
! near Infra-red 1
ZALB2(:)=XVALB5-XVALB6*ZWORK(:)
ZALB2(:)=MAX(ZALBNIR1,ZALB2(:))
!
! near Infra-red 2
ZDIAM(:)=MIN(XVDIOP1,ZDIAM(:))
ZWORK(:)=SQRT(ZDIAM(:))
ZALB3(:)=XVALB7*ZDIAM(:)-XVALB8*ZWORK(:)+XVALB9
ZALB3(:)=MAX(ZALBNIR2,ZALB3(:))
!
PSPECTRALALBEDO(:,1)=ZALB1(:)
PSPECTRALALBEDO(:,2)=ZALB2(:)
PSPECTRALALBEDO(:,3)=ZALB3(:)
!
! 3. total albedo :
! -----------------
!
PALBEDOSC(:)=XVSPEC1*ZALB1(:)+XVSPEC2*ZALB2(:)+XVSPEC3*ZALB3(:)
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LALB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LALB
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LFALL(PTSTEP,PSR,PTA,PVMOD,PSNOW,PSNOWRHO,PSNOWDZ,        &
                      PSNOWHEAT,PSNOWHMASS,PSNOWAGE,PPERMSNOWFRAC)  
!
!!    PURPOSE
!!    -------
!     Calculate changes to snowpack resulting from snowfall.
!     Update mass and heat content of uppermost layer.
!
!
USE MODD_CSTS,     ONLY : XLMTT, XTT, XCI
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN, &
                          XSNOWFALL_A_SN,         &
                          XSNOWFALL_B_SN,         &
                          XSNOWFALL_C_SN
!                   
USE YOMHOOK   ,    ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,    ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PSR, PTA, PVMOD, PPERMSNOWFRAC
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOW
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ, PSNOWHEAT, PSNOWAGE
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWHMASS
!
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWFALL, ZRHOSNEW,        &
                                       ZSNOW, ZSNOWTEMP,           &
                                       ZSNOWFALL_DELTA, ZSCAP,     &
                                       ZAGENEW
!                               
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LFALL',0,ZHOOK_HANDLE)
!
INI             = SIZE(PSNOWDZ(:,:),1)
INLVLS          = SIZE(PSNOWDZ(:,:),2)
!
ZRHOSNEW(:)     = XRHOSMIN_ES
ZAGENEW (:)     = 0.0
ZSNOWFALL(:)    = 0.0
ZSCAP(:)        = 0.0
ZSNOW(:)        = PSNOW(:)
!
PSNOWHMASS(:)   = 0.0
!
! 1. Incorporate snowfall into snowpack:
! --------------------------------------
!
!
! Heat content of newly fallen snow (J/m2):
! NOTE for now we assume the snowfall has
! the temperature of the snow surface upon reaching the snow.
! This is done as opposed to using the air temperature since
! this flux is quite small and has little to no impact
! on the time scales of interest. If we use the above assumption
! then, then the snowfall advective heat flux is zero.
!
ZSNOWTEMP(:)  = XTT
ZSCAP    (:)  = SNOW3LSCAP(PSNOWRHO(:,1))
!
WHERE (PSR(:) > 0.0 .AND. PSNOWDZ(:,1)>0.)
  ZSNOWTEMP(:)  = XTT + (PSNOWHEAT(:,1) +                              &
                    XLMTT*PSNOWRHO(:,1)*PSNOWDZ(:,1))/                   &
                    (ZSCAP(:)*MAX(XSNOWDMIN/INLVLS,PSNOWDZ(:,1)))  
  ZSNOWTEMP(:)  = MIN(XTT, ZSNOWTEMP(:))
END WHERE
!
WHERE (PSR(:) > 0.0)
!
  PSNOWHMASS(:) = PSR(:)*(XCI*(ZSNOWTEMP(:)-XTT)-XLMTT)*PTSTEP
!
! Snowfall density: Following CROCUS (Pahaut 1976)
!
   ZRHOSNEW(:)   = MAX(XRHOSMIN_ES, XSNOWFALL_A_SN + XSNOWFALL_B_SN*(PTA(:)-XTT)+         &
                     XSNOWFALL_C_SN*SQRT(PVMOD(:)))  
!
!
! Fresh snowfall changes the snowpack age,
! decreasing in uppermost snow layer (mass weighted average):
!
   PSNOWAGE(:,1) = (PSNOWAGE(:,1)*PSNOWDZ(:,1)*PSNOWRHO(:,1)+ZAGENEW(:)*PSR(:)*PTSTEP) / &
                   (PSNOWDZ(:,1)*PSNOWRHO(:,1)+PSR(:)*PTSTEP)
!
! Augment total pack depth:
!
   ZSNOWFALL(:)  = PSR(:)*PTSTEP/ZRHOSNEW(:)    ! snowfall thickness (m)
!
   PSNOW(:)      = PSNOW(:) + ZSNOWFALL(:)
!
! Fresh snowfall changes the snowpack
! density, increases the total liquid water
! equivalent: in uppermost snow layer:
!
   PSNOWRHO(:,1) = (PSNOWDZ(:,1)*PSNOWRHO(:,1) + ZSNOWFALL(:)*ZRHOSNEW(:))/     &
                   (PSNOWDZ(:,1)+ZSNOWFALL(:))  
!
   PSNOWDZ(:,1)  = PSNOWDZ(:,1) + ZSNOWFALL(:)
!
! Add energy of snowfall to snowpack:
! Update heat content (J/m2) (therefore the snow temperature
! and liquid content):
!
   PSNOWHEAT(:,1)  = PSNOWHEAT(:,1) + PSNOWHMASS(:)
!
END WHERE
!
!
! 2. Case of new snowfall on a previously snow-free surface:
! ----------------------------------------------------------
!
! When snow first falls on a surface devoid of snow,
! redistribute the snow mass throughout the 3 layers:
! (temperature already set in the calling routine
! for this case)
!
ZSNOWFALL_DELTA(:)    = 0.0
WHERE(ZSNOW(:) == 0.0 .AND. PSR(:) > 0.0)
   ZSNOWFALL_DELTA(:) = 1.0
END WHERE
!
DO JJ=1,INLVLS
   DO JI=1,INI
!
      PSNOWDZ(JI,JJ)   = ZSNOWFALL_DELTA(JI)*(ZSNOWFALL(JI) /INLVLS) + &
                        (1.0-ZSNOWFALL_DELTA(JI))*PSNOWDZ(JI,JJ)  
!
      PSNOWHEAT(JI,JJ) = ZSNOWFALL_DELTA(JI)*(PSNOWHMASS(JI)/INLVLS) + &
                       (1.0-ZSNOWFALL_DELTA(JI))*PSNOWHEAT(JI,JJ)  
!
      PSNOWRHO(JI,JJ)  = ZSNOWFALL_DELTA(JI)*ZRHOSNEW(JI)            + &
                       (1.0-ZSNOWFALL_DELTA(JI))*PSNOWRHO(JI,JJ)
!
      PSNOWAGE(JI,JJ)  = ZSNOWFALL_DELTA(JI)*(ZAGENEW(JI)/INLVLS)    + &
                       (1.0-ZSNOWFALL_DELTA(JI))*PSNOWAGE(JI,JJ)  
!
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LFALL',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SNOW3LFALL
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LCOMPACTN(PTSTEP,PSNOWDZMIN,PSNOWRHO,PSNOWDZ,PSNOWTEMP,PSNOW,PSNOWLIQ)  
!
!!    PURPOSE
!!    -------
!     Snow compaction due to overburden and settling.
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. Method
!     of Brun et al (1989) and Vionnet et al. (2012)
!
!     
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODD_CSTS,     ONLY : XTT, XG
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES
!
USE MODD_SNOW_METAMO, ONLY : XVVISC1,XVVISC3,XVVISC4, &
                             XVVISC5,XVVISC6,XVRO11
!
USE YOMHOOK   ,    ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,    ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
REAL, INTENT(IN)                    :: PSNOWDZMIN
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP, PSNOWLIQ
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOW
!
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHO2, ZVISCOCITY, ZF1, &
                                                      ZTEMP, ZSMASS, ZSNOWDZ,     &
                                                      ZWSNOWDZ, ZWHOLDMAX
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LCOMPACTN',0,ZHOOK_HANDLE)
!
INI             = SIZE(PSNOWDZ(:,:),1)
INLVLS          = SIZE(PSNOWDZ(:,:),2)
!
ZSNOWRHO2 (:,:) = PSNOWRHO(:,:)
ZSNOWDZ   (:,:) = MAX(PSNOWDZMIN,PSNOWDZ(:,:))
ZVISCOCITY(:,:) = 0.0
ZTEMP     (:,:) = 0.0
!
! 1. Cumulative snow mass (kg/m2):
! --------------------------------
!
ZSMASS(:,:) = 0.0
DO JJ=2,INLVLS
   DO JI=1,INI
      ZSMASS(JI,JJ) = ZSMASS(JI,JJ-1) + PSNOWDZ(JI,JJ-1)*PSNOWRHO(JI,JJ-1)
   ENDDO
ENDDO
! overburden of half the mass of the uppermost layer applied to itself
ZSMASS(:,1) = 0.5 * PSNOWDZ(:,1) * PSNOWRHO(:,1)
!
! 2. Compaction
! -------------
!
!Liquid water effect
!
ZWHOLDMAX(:,:) = SNOW3LHOLD(PSNOWRHO,PSNOWDZ)
ZWHOLDMAX(:,:) = MAX(1.E-10, ZWHOLDMAX(:,:))
ZF1(:,:) = 1.0/(XVVISC5+10.*MIN(1.0,PSNOWLIQ(:,:)/ZWHOLDMAX(:,:)))
!
!Snow viscocity, density and grid thicknesses
!
DO JJ=1,INLVLS
   DO JI=1,INI
!   
      IF(PSNOWRHO(JI,JJ) < XRHOSMAX_ES)THEN
!   
!       temperature dependence limited to 5K: Schleef et al. (2014)
        ZTEMP     (JI,JJ) = XVVISC4*MIN(5.0,ABS(XTT-PSNOWTEMP(JI,JJ)))
!        
!       Calculate snow viscocity: Brun et al. (1989), Vionnet et al. (2012)
        ZVISCOCITY(JI,JJ) = XVVISC1*ZF1(JI,JJ)*EXP(XVVISC3*PSNOWRHO(JI,JJ)+ZTEMP(JI,JJ))*PSNOWRHO(JI,JJ)/XVRO11
!
!       Calculate snow density:
        ZSNOWRHO2(JI,JJ) = PSNOWRHO(JI,JJ) + PSNOWRHO(JI,JJ)*PTSTEP &
                         * ( (XG*ZSMASS(JI,JJ)/ZVISCOCITY(JI,JJ)) )
!         
!       Conserve mass by decreasing grid thicknesses in response to density increases
        PSNOWDZ(JI,JJ) = PSNOWDZ(JI,JJ)*(PSNOWRHO(JI,JJ)/ZSNOWRHO2(JI,JJ))  
!        
      ENDIF
!
   ENDDO
ENDDO
!
! 3. Update total snow depth and density profile:
! -----------------------------------------------
!
! Compaction/augmentation of total snowpack depth
!
PSNOW(:) = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      PSNOW(JI) = PSNOW(JI) + PSNOWDZ(JI,JJ)
   ENDDO
ENDDO
!
! Update density (kg m-3):
!
PSNOWRHO(:,:)  = ZSNOWRHO2(:,:)
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LCOMPACTN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LCOMPACTN
!####################################################################
!####################################################################
!####################################################################
        SUBROUTINE SNOW3LTRANSF(PSNOW,PSNOWDZ,PSNOWDZN,    &
                                PSNOWRHO,PSNOWHEAT,PSNOWAGE)  
!
!!    PURPOSE
!!    -------
!     Snow mass,heat and characteristics redistibution in case of
!     grid resizing. Total mass and heat content of the overall snowpack
!     unchanged/conserved within this routine.
!     Same method as in Crocus
!
USE MODD_SURF_PAR, ONLY : XUNDEF          
USE MODD_SNOW_PAR, ONLY : XSNOWCRITD
!
USE YOMHOOK   ,    ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,    ONLY : JPRB
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:  ), INTENT(IN)    :: PSNOW
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZN  
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWAGE
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JI, JL, JLO
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHON
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEATN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWAGEN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWZTOP_NEW
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWZBOT_NEW                                       
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHOO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEATO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWAGEO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWDZO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWZTOP_OLD
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWZBOT_OLD  
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEAN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWAGN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZMASTOTN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZMASSDZO
!
REAL, DIMENSION(SIZE(PSNOW)) :: ZPSNOW_OLD, ZPSNOW_NEW
REAL, DIMENSION(SIZE(PSNOW)) :: ZSUMHEAT, ZSUMSWE, ZSUMAGE, ZSNOWMIX_DELTA
!
REAL :: ZPROPOR
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LTRANSF',0,ZHOOK_HANDLE)
!
INI        = SIZE(PSNOWRHO,1)
INLVLS     = SIZE(PSNOWRHO,2)
!
ZPSNOW_NEW(:) = 0.0
ZPSNOW_OLD(:) = PSNOW(:)
!
DO JL=1,INLVLS
   DO JI=1,INI
      ZPSNOW_NEW(JI)=ZPSNOW_NEW(JI)+PSNOWDZN(JI,JL)
   ENDDO
ENDDO
!
! initialization of variables describing the initial snowpack 
!
ZSNOWDZO  (:,:) = PSNOWDZ  (:,:)
ZSNOWRHOO (:,:) = PSNOWRHO (:,:)
ZSNOWHEATO(:,:) = PSNOWHEAT(:,:)
ZSNOWAGEO (:,:) = PSNOWAGE (:,:)
ZMASSDZO  (:,:) = XUNDEF
!
! 1. Calculate vertical grid limits (m):
! --------------------------------------
!
ZSNOWZTOP_OLD(:,1) = ZPSNOW_OLD(:)
ZSNOWZTOP_NEW(:,1) = ZPSNOW_NEW(:)
ZSNOWZBOT_OLD(:,1) = ZSNOWZTOP_OLD(:,1)-ZSNOWDZO(:,1)
ZSNOWZBOT_NEW(:,1) = ZSNOWZTOP_NEW(:,1)-PSNOWDZN(:,1)
!
DO JL=2,INLVLS
   DO JI=1,INI
      ZSNOWZTOP_OLD(JI,JL) = ZSNOWZBOT_OLD(JI,JL-1)
      ZSNOWZTOP_NEW(JI,JL) = ZSNOWZBOT_NEW(JI,JL-1)
      ZSNOWZBOT_OLD(JI,JL) = ZSNOWZTOP_OLD(JI,JL  )-ZSNOWDZO(JI,JL)
      ZSNOWZBOT_NEW(JI,JL) = ZSNOWZTOP_NEW(JI,JL  )-PSNOWDZN(JI,JL)
   ENDDO
ENDDO
ZSNOWZBOT_OLD(:,INLVLS)=0.0
ZSNOWZBOT_NEW(:,INLVLS)=0.0
!
! 3. Calculate mass, heat, charcateristics mixing due to vertical grid resizing:
! --------------------------------------------------------------------
!
! loop over the new snow layers
! Summ or avergage of the constituting quantities of the old snow layers
! which are totally or partially inserted in the new snow layer
! For snow age, mass weighted average is used.
!
ZSNOWHEAN(:,:)=0.0
ZMASTOTN (:,:)=0.0
ZSNOWAGN (:,:)=0.0
!
DO JL=1,INLVLS
   DO JLO=1, INLVLS   
      DO JI=1,INI
        IF((ZSNOWZTOP_OLD(JI,JLO)>ZSNOWZBOT_NEW(JI,JL)).AND.(ZSNOWZBOT_OLD(JI,JLO)<ZSNOWZTOP_NEW(JI,JL)))THEN
!                
          ZPROPOR = (MIN(ZSNOWZTOP_OLD(JI,JLO), ZSNOWZTOP_NEW(JI,JL)) &
                  -  MAX(ZSNOWZBOT_OLD(JI,JLO), ZSNOWZBOT_NEW(JI,JL)))&
                  / ZSNOWDZO(JI,JLO) 
!
          ZMASSDZO (JI,JLO)=ZSNOWRHOO(JI,JLO)*ZSNOWDZO(JI,JLO)*ZPROPOR
!
          ZMASTOTN (JI,JL)=ZMASTOTN (JI,JL)+ZMASSDZO  (JI,JLO)
          ZSNOWAGN (JI,JL)=ZSNOWAGN (JI,JL)+ZSNOWAGEO (JI,JLO)*ZMASSDZO(JI,JLO)
!
          ZSNOWHEAN(JI,JL)=ZSNOWHEAN(JI,JL)+ZSNOWHEATO(JI,JLO)*ZPROPOR
!          
        ENDIF
      ENDDO 
    ENDDO 
ENDDO  
!
! the new layer inherits from the weighted average properties of the old ones
! heat and mass
!
ZSNOWHEATN(:,:)= ZSNOWHEAN(:,:)
ZSNOWAGEN (:,:)= ZSNOWAGN (:,:)/ZMASTOTN(:,:)
ZSNOWRHON (:,:)= ZMASTOTN (:,:)/PSNOWDZN(:,:)
!
!
! 4. Vanishing or very thin snowpack check:
! -----------------------------------------
!
! NOTE: ONLY for very shallow snowpacks, mix properties (homogeneous):
! this avoids problems related to heat and mass exchange for
! thin layers during heavy snowfall or signifigant melt: one
! new/old layer can exceed the thickness of several old/new layers.
! Therefore, mix (conservative):
!
ZSUMHEAT(:)       = 0.0
ZSUMSWE(:)        = 0.0
ZSUMAGE(:)        = 0.0
ZSNOWMIX_DELTA(:) = 0.0
!
DO JL=1,INLVLS
   DO JI=1,INI
      IF(PSNOW(JI) < XSNOWCRITD)THEN
         ZSUMHEAT      (JI) = ZSUMHEAT(JI) + PSNOWHEAT(JI,JL)
         ZSUMSWE       (JI) = ZSUMSWE (JI) + PSNOWRHO (JI,JL)*PSNOWDZ(JI,JL)
         ZSUMAGE       (JI) = ZSUMAGE (JI) + PSNOWAGE (JI,JL)
         ZSNOWMIX_DELTA(JI) = 1.0
      ENDIF
   ENDDO
ENDDO
!
! Heat and mass are evenly distributed vertically:
! heat and mass (density and thickness) are constant
! in profile:
!
DO JL=1,INLVLS
   DO JI=1,INI
!         
      ZSNOWHEATN(JI,JL) = ZSNOWMIX_DELTA(JI)*(ZSUMHEAT(JI)/INLVLS)  + &
                         (1.0-ZSNOWMIX_DELTA(JI))*ZSNOWHEATN(JI,JL)  
!
      PSNOWDZN(JI,JL)   = ZSNOWMIX_DELTA(JI)*(PSNOW(JI)/INLVLS)     + &
                         (1.0-ZSNOWMIX_DELTA(JI))*PSNOWDZN(JI,JL)  
!
      ZSNOWRHON(JI,JL)  = ZSNOWMIX_DELTA(JI)*(ZSUMSWE(JI)/PSNOW(JI)) + &
                         (1.0-ZSNOWMIX_DELTA(JI))*ZSNOWRHON(JI,JL)
!
      ZSNOWAGEN(JI,JL)  = ZSNOWMIX_DELTA(JI)*(ZSUMAGE(JI)/INLVLS)  + &
                         (1.0-ZSNOWMIX_DELTA(JI))*ZSNOWAGEN(JI,JL)  
!
   ENDDO
ENDDO
!
! 5. Update mass (density and thickness) and heat:
! ------------------------------------------------
!
PSNOWDZ  (:,:) = PSNOWDZN  (:,:)
PSNOWRHO (:,:) = ZSNOWRHON (:,:)
PSNOWHEAT(:,:) = ZSNOWHEATN(:,:)
PSNOWAGE (:,:) = ZSNOWAGEN (:,:)
!
IF (LHOOK) CALL DR_HOOK('MODE_SNOW3L:SNOW3LTRANSF',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LTRANSF
!####################################################################
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROTHRM(PSNOWRHO,PSCOND,PSNOWTEMP,PPS,PSNOWLIQ, &
                       HSNOWCOND                  )
!
!!    PURPOSE
!!    -------
!     Calculate snow thermal conductivity from
!     Sun et al. 1999, J. of Geophys. Res., 104, 19587-19579
!     (vapor) and Anderson, 1976, NOAA Tech. Rep. NWS 19 (snow).
!
!     Upon activation of flag OCOND_YEN, use the Yen (1981) formula for thermal conductivity
!     This formula was originally used in Crocus.
!
!     05/2016 : Lafaysse/Cluzet : new available options
!
USE MODD_CSTS, ONLY : XP00, XCONDI, XRHOLW
USE MODD_SNOW_PAR, ONLY : XSNOWTHRMCOND1, XSNOWTHRMCOND2, XSNOWTHRMCOND_AVAP, &
                          XSNOWTHRMCOND_BVAP, XSNOWTHRMCOND_CVAP, XVRKZ6, &
                          XSNOWTHRMCOND_C11_1, XSNOWTHRMCOND_C11_2, XSNOWTHRMCOND_C11_3
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)      :: PPS
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP, PSNOWRHO, PSNOWLIQ
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSCOND
!
CHARACTER(3), INTENT(IN)              :: HSNOWCOND ! conductivity option
!
!*      0.2    declarations of local variables
!
INTEGER :: JJ, JST ! looping indexes
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SNOWCROTHRM',0,ZHOOK_HANDLE)
!
! 1. Snow thermal conductivity
! ----------------------------
!
DO JST = 1,SIZE(PSNOWRHO(:,:),2)
  !
  DO JJ = 1,SIZE(PSNOWRHO(:,:),1)
    ! Cluzet et al 2016
    IF ( HSNOWCOND=='Y81') THEN
      PSCOND(JJ,JST) = XCONDI * EXP( XVRKZ6 * LOG( PSNOWRHO(JJ,JST)/XRHOLW ) )
      ! Snow thermal conductivity is set to be above 0.04 W m-1 K-1
      PSCOND(JJ,JST) = MAX( 0.04, PSCOND(JJ,JST) )
    ELSE IF(HSNOWCOND == 'I02') THEN
      PSCOND(JJ,JST) = ( XSNOWTHRMCOND1 + &
                         XSNOWTHRMCOND2 * PSNOWRHO(JJ,JST) * PSNOWRHO(JJ,JST) ) + &
                         MAX( 0.0, ( XSNOWTHRMCOND_AVAP + &
                                    ( XSNOWTHRMCOND_BVAP/(PSNOWTEMP(JJ,JST) + XSNOWTHRMCOND_CVAP) ) ) &
                                   * (XP00/PPS(JJ)) ) 
    ELSE IF(HSNOWCOND == 'C11') THEN
      PSCOND(JJ,JST) = XSNOWTHRMCOND_C11_1 * PSNOWRHO(JJ,JST)* PSNOWRHO(JJ,JST) + &
      		       XSNOWTHRMCOND_C11_2 * PSNOWRHO(JJ,JST) + XSNOWTHRMCOND_C11_3
    ENDIF
    !
    ! In older versions, snow thermal conductivity is annihilated in presence of liquid water.
    ! We decided to remove this incorrect parameterization (May 2016)
     !
   ENDDO ! end loop JST
   !
ENDDO ! end loop JST
!
IF (LHOOK) CALL DR_HOOK('SNOWCROTHRM',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROTHRM
!####################################################################

END MODULE MODE_SNOW3L

