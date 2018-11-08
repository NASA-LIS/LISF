cvk original Noah code
cvk      SUBROUTINE EVAPO (ETA1,SMC,NSOIL,CMC,ETP1,DT,ZSOIL,
cvk     &                  SH2O,
cvk     &                  SMCMAX,BEXP,PC,SMCWLT,DKSAT,DWSAT,
cvk     &                  SMCREF,SHDFAC,CMCMAX,
cvk     &                  SMCDRY,CFACTR,
cvk     &                  EDIR1,EC1,ET1,ETT1,SFCTMP,Q2,NROOT,RTDIS,FXEXP)

      SUBROUTINE EVAPO (ETA1,NSOIL,CMC,ETP1,etpbare,SH2O,SMCMAX,PC,
     &                  SMCWLT,SMCREF,SHDFAC,CMCMAX,SMCDRY,CFACTR,
     &                  EDIR1,EC1,ET1,ETT1,NROOT,RTDIS,FXEXP,avgsmup,
     +                  uztwh,uztwm,bareadj)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE EVAPO
C ----------------------------------------------------------------------
C CALCULATE SOIL MOISTURE FLUX.  THE LIQUID SOIL MOISTURE CONTENT (SH2O - 
C A PER UNIT VOLUME MEASUREMENT) IS A DEPENDENT VARIABLE THAT IS UPDATED 
C WITH PROGNOSTIC EQNS. THE CANOPY MOISTURE CONTENT (CMC) IS ALSO UPDATED.
C FROZEN GROUND VERSION:  NEW STATES ADDED: SH2O, AND FROZEN GROUND
C CORRECTION FACTOR, FRZFACT AND PARAMETER SLOPE.
cvk Notes
c DT in seconds
c SMCWLT and SMCDRY are the same in REDPRM subroutine; however, could bedifferent
c ET1() - transpiration from each layer
c ETT1 - total transpiration
c EDIR1 - bare soil evaporation (mostly from top layer)
c EC1 - canopy water evaporation
c ETA1 - total evapotranspiration 
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER I
      INTEGER K
      INTEGER NSOIL
      INTEGER NROOT

      REAL avgsmup,uztwh,uztwm,etpbare,bareadj
      REAL CFACTR
      REAL CMC
      REAL CMC2MS
      REAL CMCMAX
c      REAL DEVAP
cvk      REAL DKSAT
      REAL DT
cvk      REAL DWSAT
      REAL EC1
      REAL EDIR1
      REAL ET1(NSOIL)
      REAL ETA1
      REAL ETP1
      REAL ETT1
      REAL FXEXP
      REAL PC
cvk      REAL Q2
      REAL RTDIS(NSOIL)
cvk      REAL SFCTMP
      REAL SHDFAC
cvk      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
cvk      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE IF THE POTENTIAL EVAPOTRANSPIRATION IS
C GREATER THAN ZERO.
C ----------------------------------------------------------------------
      EDIR1 = 0.
      EC1 = 0.
      DO K = 1,NSOIL
        ET1(K) = 0.
      END DO
      ETT1 = 0.

      IF (ETP1 .GT. 0.0 .or. etpbare .gt. 0.0) THEN

C ----------------------------------------------------------------------
C RETRIEVE DIRECT EVAPORATION FROM SOIL SURFACE.  CALL THIS FUNCTION
C ONLY IF VEG COVER NOT COMPLETE.
C FROZEN GROUND VERSION:  SH2O STATES REPLACE SMC STATES.
cvk 4/2010  average SM replaces just upper layer SM
C ----------------------------------------------------------------------
        IF (SHDFAC .LT. 1.) THEN
cvk 2010        CALL DEVAP (EDIR1,ETP1,SH2O(1),SHDFAC,SMCMAX,
cvk 2010     &              SMCDRY,SMCREF,SMCWLT,FXEXP)
        CALL DEVAP (EDIR1,etpbare,avgsmup,SHDFAC,SMCMAX,SMCDRY,
     &                 SMCREF,SMCWLT,FXEXP,uztwh,uztwm,bareadj)
        ENDIF

cvk NOTE........................................................................
cvk there is a problem in reduction of EPT1 if canopy water removal rate
cvk is bigger than available CMC. To fix this, first canopy water evaporation
cvk EC1 is estimated, and then transpiration with the use of estimated EC1: 
cvk insted of 1-(CMC/CMCMAX)**CFACTR use 1-EC1/(SHDFAC*ETP1) in subroutine TRANSP
cvk..............................................................................

C ----------------------------------------------------------------------
C INITIALIZE PLANT TOTAL TRANSPIRATION, RETRIEVE PLANT TRANSPIRATION,
C AND ACCUMULATE IT FOR ALL SOIL LAYERS.
C ----------------------------------------------------------------------
cvk TRANSP moved after canopy evaporation
cvk        IF (SHDFAC.GT.0.0) THEN
cvk
cvk          CALL TRANSP (ET1,NSOIL,ETP1,SH2O,CMC,SHDFAC,SMCWLT,
cvk     &                 CMCMAX,PC,CFACTR,SMCREF,NROOT,RTDIS)
cvk
cvk          DO K = 1,NSOIL
cvk            ETT1 = ETT1 + ET1(K)
cvk          END DO

C ----------------------------------------------------------------------
C CALCULATE CANOPY EVAPORATION.
C IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR CMC=0.0.
C ----------------------------------------------------------------------
        IF (SHDFAC.GT.0.0) THEN
          IF (CMC .GT. 0.0) THEN
            EC1 = SHDFAC * ( ( CMC / CMCMAX ) ** CFACTR ) * ETP1
          ELSE
            EC1 = 0.0
          ENDIF

C ----------------------------------------------------------------------
C EC SHOULD BE LIMITED BY THE TOTAL AMOUNT OF AVAILABLE WATER ON THE
C CANOPY.  -F.CHEN, 18-OCT-1994
C ----------------------------------------------------------------------
cvk  Evaporation in SAC will be in mm/dt, so no need to convert by DT
cvk          CMC2MS = CMC / DT
cvk          EC1 = MIN ( CMC2MS, EC1 )
cvk....................................
          EC1 = MIN ( CMC, EC1 )
cvk        ENDIF
cvk      ENDIF

C ----------------------------------------------------------------------
C INITIALIZE PLANT TOTAL TRANSPIRATION, RETRIEVE PLANT TRANSPIRATION,
C AND ACCUMULATE IT FOR ALL SOIL LAYERS.
C ----------------------------------------------------------------------
cvk 5/2009          CALL TRANSP (ET1,NSOIL,ETP1,SH2O,CMC,SHDFAC,SMCWLT,
          CALL TRANSP (ET1,NSOIL,ETP1,SH2O,CMC,SHDFAC,SMCDRY,
     &                 CMCMAX,PC,CFACTR,SMCREF,NROOT,RTDIS,EC1)

          DO K = 1,NSOIL
            ETT1 = ETT1 + ET1(K)
          END DO

        ENDIF
      ENDIF

C ----------------------------------------------------------------------
C TOTAL UP EVAP AND TRANSP TYPES TO OBTAIN ACTUAL EVAPOTRANSP
C ----------------------------------------------------------------------
      ETA1 = EDIR1 + ETT1 + EC1

C ----------------------------------------------------------------------
C END SUBROUTINE EVAPO
C ----------------------------------------------------------------------
      RETURN
      END
