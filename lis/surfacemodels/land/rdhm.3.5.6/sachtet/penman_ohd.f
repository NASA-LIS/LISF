      SUBROUTINE PENMAN_OHD (SFCTMP,SFCPRS,CH,T2V,TH2,PRCP,FDOWN,T24,
     &                   SSOIL, Q2,Q2SAT,ETP,RCH,EPSCA,RR,SNOWNG,FRZGRA,
     &                   DQSDT2,FLX2)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE PENMAN_OHD
C ----------------------------------------------------------------------
C CALCULATE POTENTIAL EVAPORATION FOR THE CURRENT POINT.  VARIOUS
C PARTIAL SUMS/PRODUCTS ARE ALSO CALCULATED AND PASSED BACK TO THE
C CALLING ROUTINE FOR LATER USE.
C ----------------------------------------------------------------------
      LOGICAL SNOWNG
      LOGICAL FRZGRA

      REAL A
      REAL BETA
      REAL CH
      REAL CP
      REAL CPH2O
      REAL CPICE
      REAL DELTA
      REAL DQSDT2
      REAL ELCP
      REAL EPSCA
      REAL ETP
      REAL FDOWN
      REAL FLX2
      REAL FNET
      REAL LSUBC
      REAL LSUBF
      REAL PRCP
      REAL Q2
      REAL Q2SAT
      REAL R
      REAL RAD
      REAL RCH
      REAL RHO
      REAL RR
      REAL SSOIL
      REAL SFCPRS
      REAL SFCTMP
      REAL SIGMA
      REAL T24
      REAL T2V
      REAL TH2

      PARAMETER(CP = 1004.6)
      PARAMETER(CPH2O = 4.218E+3)
      PARAMETER(CPICE = 2.106E+3)
      PARAMETER(R = 287.04)
      PARAMETER(ELCP = 2.4888E+3)
      PARAMETER(LSUBF = 3.335E+5)
      PARAMETER(LSUBC = 2.501000E+6)
      PARAMETER(SIGMA = 5.67E-8)

C ----------------------------------------------------------------------
C EXECUTABLE CODE BEGINS HERE:
C ----------------------------------------------------------------------
      FLX2 = 0.0

C ----------------------------------------------------------------------
C PREPARE PARTIAL QUANTITIES FOR PENMAN_OHD EQUATION.
C ----------------------------------------------------------------------
      DELTA = ELCP * DQSDT2
      T24 = SFCTMP * SFCTMP * SFCTMP * SFCTMP
      RR = T24 * 6.48E-8 /(SFCPRS * CH) + 1.0
      RHO = SFCPRS / (R * T2V)
      RCH = RHO * CP * CH

C ----------------------------------------------------------------------
C ADJUST THE PARTIAL SUMS / PRODUCTS WITH THE LATENT HEAT
C EFFECTS CAUSED BY FALLING PRECIPITATION.
C ----------------------------------------------------------------------
      IF (.NOT. SNOWNG) THEN
        IF (PRCP .GT. 0.0) RR = RR + CPH2O*PRCP/RCH
      ELSE
        RR = RR + CPICE*PRCP/RCH
      ENDIF

      FNET = FDOWN - SIGMA*T24 - SSOIL

C ----------------------------------------------------------------------
C INCLUDE THE LATENT HEAT EFFECTS OF FRZNG RAIN CONVERTING TO ICE ON
C IMPACT IN THE CALCULATION OF FLX2 AND FNET.
C ----------------------------------------------------------------------
      IF (FRZGRA) THEN
        FLX2 = -LSUBF * PRCP
        FNET = FNET - FLX2
      ENDIF

C ----------------------------------------------------------------------
C FINISH PENMAN_OHD EQUATION CALCULATIONS.
C ----------------------------------------------------------------------
      RAD = FNET/RCH + TH2 - SFCTMP
      A = ELCP * (Q2SAT - Q2)
      EPSCA = (A*RR + RAD*DELTA) / (DELTA + RR)
      ETP = EPSCA * RCH / LSUBC

c      write(*,321) 'pen',FNET,SFCPRS,SFCTMP,CH,Q2,Q2SAT,DQSDT2,RR,DELTA
c      write(*,*) rad,rch,th2 ,etp
c321   format(a3,f6.0,f8.0,f9.1,3f7.4,f7.4,2f7.3)

C ----------------------------------------------------------------------
C END SUBROUTINE PENMAN_OHD
C ----------------------------------------------------------------------
      RETURN
      END
      
