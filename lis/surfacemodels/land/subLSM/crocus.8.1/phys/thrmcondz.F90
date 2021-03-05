!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
    SUBROUTINE THRMCONDZ(PSANDZ,PWSATZ,PCONDDRY,PCONDSLD)
!   ###############################################################
!!****  *THRMCONDZ*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates soil thermal conductivity components
!     using sand fraction and model constants in
!     order to calculate the thermal conductivity
!     following the method of Johansen (1975) as recommended
!     by Farouki (1986) parameterized for SVAT schemes
!     following Peters-Lidard et al. 1998 (JAS). This is
!     used in explicit calculation of CG (soil thermal
!     inertia): it is an option. DEFAULT is method of
!     Noilhan and Planton (1989) (see SOIL.F90).
!              
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!    Peters-Lidard et al. 1998 (JAS)
!!      
!!    AUTHOR
!!    ------
!!
!!      A. Boone           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    25/03/99
!!                  18/02/00    2D for veritcal profiles
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_ISBA_PAR,   ONLY : XDRYWGHT, XSPHSOIL, XCONDQRTZ, XCONDOTH1, XCONDOTH2
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:,:), INTENT(IN) :: PSANDZ     ! soil sand fraction (-)
REAL,   DIMENSION(:,:), INTENT(IN) :: PWSATZ     ! soil porosity (m3 m-3)
!
REAL,   DIMENSION(:,:), INTENT(OUT):: PCONDDRY  ! soil dry thermal conductivity
!                                                 (W m-1 K-1)
REAL,   DIMENSION(:,:), INTENT(OUT):: PCONDSLD  ! soil solids thermal
!                                                 conductivity (W m-1 K-1)
!
!*      0.2    declarations of local variables
!
REAL,    DIMENSION(SIZE(PSANDZ,1),SIZE(PSANDZ,2)) :: ZQUARTZ, ZGAMMAD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('THRMCONDZ',0,ZHOOK_HANDLE)
ZQUARTZ(:,:)   = XUNDEF
ZGAMMAD(:,:)   = XUNDEF
PCONDSLD(:,:)  = XUNDEF
PCONDDRY(:,:)  = XUNDEF
!
!
! Quartz content estimated from sand fraction:
!
WHERE(PSANDZ(:,:)/=XUNDEF)
!
   ZQUARTZ(:,:)   = 0.038 + 0.95*PSANDZ(:,:)
!
! Note, ZGAMMAD (soil dry density) can be supplied from obs, but
! for mesoscale modeling, we use the following approximation
! from Peters-Lidard et al. 1998:
!
   ZGAMMAD(:,:)   = (1.0-PWSATZ(:,:))*XDRYWGHT
!
END WHERE
!
! Soil solids conductivity:
!
WHERE(ZQUARTZ >  0.20 .AND. PSANDZ(:,:)/=XUNDEF)
   PCONDSLD(:,:)  = (XCONDQRTZ**ZQUARTZ(:,:))*                        &
                    (XCONDOTH1**(1.0-ZQUARTZ(:,:)))  
END WHERE
WHERE(ZQUARTZ <= 0.20 .AND. PSANDZ(:,:)/=XUNDEF)
   PCONDSLD(:,:)  = (XCONDQRTZ**ZQUARTZ(:,:))*                        &
                    (XCONDOTH2**(1.0-ZQUARTZ(:,:)))  
ENDWHERE
!
! Soil dry conductivity:
!
WHERE(PSANDZ(:,:)/=XUNDEF)
   PCONDDRY(:,:)     = (0.135*ZGAMMAD(:,:) + 64.7)/                   &
                         (XDRYWGHT - 0.947*ZGAMMAD(:,:))  
END WHERE
IF (LHOOK) CALL DR_HOOK('THRMCONDZ',1,ZHOOK_HANDLE)
!
!
!
END SUBROUTINE THRMCONDZ
