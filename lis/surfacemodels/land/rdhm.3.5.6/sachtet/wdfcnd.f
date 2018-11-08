      SUBROUTINE WDFCND (WDF,WCND,SMC,SMCMAX,BEXP,DKSAT,DWSAT,
     &                   SICEMAX)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE WDFCND
C ----------------------------------------------------------------------
C CALCULATE SOIL WATER DIFFUSIVITY AND SOIL HYDRAULIC CONDUCTIVITY.
C ----------------------------------------------------------------------
        real sh2o
      REAL BEXP
      REAL DKSAT
      REAL DWSAT
      REAL EXPON
      REAL FACTR1
      REAL FACTR2
      REAL SICEMAX
      REAL SMC
      REAL SMCMAX
      REAL VKwgt
      REAL WCND
      REAL WDF

C ----------------------------------------------------------------------
C     CALC THE RATIO OF THE ACTUAL TO THE MAX PSBL SOIL H2O CONTENT
C ----------------------------------------------------------------------
      SMC = SMC
      SMCMAX = SMCMAX
      FACTR1 = 0.2 / SMCMAX
      FACTR2 = SMC / SMCMAX

C ----------------------------------------------------------------------
C PREP AN EXPNTL COEF AND CALC THE SOIL WATER DIFFUSIVITY
C ----------------------------------------------------------------------
      EXPON = BEXP + 2.0
      WDF = DWSAT * FACTR2 ** EXPON

C ----------------------------------------------------------------------
C FROZEN SOIL HYDRAULIC DIFFUSIVITY.  VERY SENSITIVE TO THE VERTICAL
C GRADIENT OF UNFROZEN WATER. THE LATTER GRADIENT CAN BECOME VERY
C EXTREME IN FREEZING/THAWING SITUATIONS, AND GIVEN THE RELATIVELY 
C FEW AND THICK SOIL LAYERS, THIS GRADIENT SUFFERES SERIOUS 
C TRUNCTION ERRORS YIELDING ERRONEOUSLY HIGH VERTICAL TRANSPORTS OF
C UNFROZEN WATER IN BOTH DIRECTIONS FROM HUGE HYDRAULIC DIFFUSIVITY.  
C THEREFORE, WE FOUND WE HAD TO ARBITRARILY CONSTRAIN WDF 
C --
C VERSION D_10CM: ........  FACTR1 = 0.2/SMCMAX
C WEIGHTED APPROACH...................... PABLO GRUNMANN, 28_SEP_1999.
C ----------------------------------------------------------------------
c      IF (SICEMAX .GT. 0.0)  THEN
c        VKWGT = 1./(1.+(500.*SICEMAX)**3.)
c        WDF = VKWGT*WDF + (1.- VKWGT)*DWSAT*FACTR1**EXPON
c      ENDIF
      
C ----------------------------------------------------------------------
C RESET THE EXPNTL COEF AND CALC THE HYDRAULIC CONDUCTIVITY
C ----------------------------------------------------------------------
      EXPON = (2.0 * BEXP) + 3.0
      WCND = DKSAT * FACTR2 ** EXPON
c       write(*,*) smc,bexp,dwsat,dksat,smcmax,sicemax,vkwgt,wdf,wcnd
C ----------------------------------------------------------------------
C END SUBROUTINE WDFCND
C ----------------------------------------------------------------------
      RETURN
      END
