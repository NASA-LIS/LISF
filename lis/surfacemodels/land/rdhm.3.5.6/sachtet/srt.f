cv      SUBROUTINE SRT(RHSTT,EDIR,ET,SH2O,SH2OA,NSOIL,
cv    &                ZSOIL,DWSAT,DKSAT,SMCMAX,BEXP, 
cv     &                SICE,AI,BI,CI)

      SUBROUTINE SRT(RHSTT,EDIR,ET,SH2O,SH2OA,NSOIL,
     &                ZSOIL,DWSAT,DKSAT,SMCMAX,BEXP, 
     &                SICE,AI,BI,CI,smc,rdst)
      IMPLICIT NONE
cv  gradient is estimated using total water content instead of liquid
c VK note: Subroutine was changed: All input fluxes except 
c VK       evapotranspiration were excluded
cfews  It includes two options of gradient calculation:
cfews    1. Noah original option (if rdst=1 or noah)that uses an actual SM gradient without 
cfews       regard potential for huge gradient jump between upper/lower zones
cfews    2. OHD option (if rdst=0 or ohd) that uses 'reference' gradient that 
cfews       accounts for the potential jump
 
C ----------------------------------------------------------------------
C SUBROUTINE SRT
C ----------------------------------------------------------------------
C CALCULATE THE RIGHT HAND SIDE OF THE TIME TENDENCY TERM OF THE SOIL
C WATER DIFFUSION EQUATION.  ALSO TO COMPUTE ( PREPARE ) THE MATRIX
C COEFFICIENTS FOR THE TRI-DIAGONAL MATRIX OF THE IMPLICIT TIME SCHEME.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      INTEGER J
      INTEGER JJ      
      INTEGER K
      INTEGER KS
      INTEGER NSOIL
      integer rdst
      real smc(NSOIL),smcx
      REAL AI(NSOLD)
      REAL BEXP
      REAL BI(NSOLD)
      REAL CI(NSOLD)
      REAL DDZ
      REAL DDZ2
      REAL DENOM
      REAL DENOM2
      REAL DICE
      REAL DKSAT
      REAL DMAX(NSOLD)
      REAL DSMDZ
      REAL DSMDZ2
      REAL DWSAT
      REAL EDIR
      REAL ET(NSOIL)
      REAL MXSMC
      REAL MXSMC2
      REAL NUMER
      REAL RHSTT(NSOIL)
      REAL SH2O(NSOIL)
      REAL SH2OA(NSOIL)
      REAL SICE(NSOIL)
      REAL SICEMAX
      REAL SLOPX
      REAL SMCMAX
      REAL WCND
      REAL WCND2
      REAL WDF
      REAL WDF2
      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C LET SICEMAX BE THE GREATEST, IF ANY, FROZEN WATER CONTENT WITHIN SOIL
C LAYERS.
C ----------------------------------------------------------------------
      SICEMAX = 0.0
      DO KS=1,NSOIL
       IF (SICE(KS) .GT. SICEMAX) SICEMAX = SICE(KS)
      END DO
C ----------------------------------------------------------------------
C TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN LINE
C BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
C 'MXSMC = MAX(SH2OA(1), SH2OA(2))'
C ----------------------------------------------------------------------

cvk 8/10  assume that tension only water diffusivity and conductivity 
cvk 8/10  is restricted by lower soil moisture value
cfews  Richards Eq. option: rdst=0 (ohd), rdst=1 (noah)  

      if(rdst .eq. 0) then
       MXSMC = MIN(SH2OA(1), SH2OA(2))
      else 
       MXSMC = SH2OA(1)
      endif 
      CALL WDFCND (WDF,WCND,MXSMC,SMCMAX,BEXP,DKSAT,DWSAT,
     &             SICEMAX)

C ----------------------------------------------------------------------
C CALC THE MATRIX COEFFICIENTS AI, BI, AND CI FOR THE TOP LAYER
C ----------------------------------------------------------------------
      DDZ = 1. / ( -.5 * ZSOIL(2) )
      AI(1) = 0.0
      BI(1) = WDF * DDZ / ( -ZSOIL(1) )
      CI(1) = -BI(1)

C ----------------------------------------------------------------------
C CALC RHSTT FOR THE TOP LAYER AFTER CALC'NG THE VERTICAL SOIL MOISTURE
C GRADIENT BTWN THE TOP AND NEXT TO TOP LAYERS.
C ----------------------------------------------------------------------
cv      DSMDZ = ( SH2O(1) - SH2O(2) ) / ( -.5 * ZSOIL(2) )

cvk 8/10  soil moisture gradient between upper and lower sac zones
cvk 8/10  may be too high during floods when using thick soil layer;
cvk 8/10  to make SM gradient more reasonable, a refernce soil moisture 
cvk 8/10  (smcx) introduced instead of lower layer SM. The reference SM is a 
cvk 8/10  weighted value of upper and lower layers depending on ratio of
cvk 8/10  upper and lower layer thicknesses
cfews
      if(rdst .eq. 0) then
       smcx=smc(1)+(smc(2)-smc(1))*zsoil(1)/(zsoil(2)-zsoil(1))
       dsmdz=(smc(1)-smcx)/( -.5 * ZSOIL(2) )
      else       
       dsmdz=(smc(1)-smc(2))/( -.5 * ZSOIL(2) )
      endif 

      RHSTT(1) = (WDF * DSMDZ + WCND + EDIR + ET(1))/ZSOIL(1)
C ----------------------------------------------------------------------
C INITIALIZE DDZ2
C ----------------------------------------------------------------------
      DDZ2 = 0.0
C ----------------------------------------------------------------------
C LOOP THRU THE REMAINING SOIL LAYERS, REPEATING THE ABV PROCESS
C ----------------------------------------------------------------------
      DO K = 2,NSOIL
        DENOM2 = (ZSOIL(K-1) - ZSOIL(K))
        IF (K .NE. NSOIL) THEN
         SLOPX = 1.
C ----------------------------------------------------------------------
C AGAIN, TO AVOID SPURIOUS DRAINAGE BEHAVIOR, 'UPSTREAM DIFFERENCING' IN
C LINE BELOW REPLACED WITH NEW APPROACH IN 2ND LINE:
C 'MXSMC2 = MAX (SH2OA(K), SH2OA(K+1))'
C ----------------------------------------------------------------------

cvk 8/10 Change similar to the first soil layer above  
cfews
          if(rdst .eq. 0) then
           MXSMC2 = MIN (SH2OA(K), SH2OA(K+1))
          else                   
           MXSMC2 = SH2OA(K)
          endif 
          CALL WDFCND (WDF2,WCND2,MXSMC2,SMCMAX,BEXP,DKSAT,DWSAT,
     &                 SICEMAX)

C ----------------------------------------------------------------------
C CALC SOME PARTIAL PRODUCTS FOR LATER USE IN CALC'NG RHSTT
C ----------------------------------------------------------------------
          DENOM = (ZSOIL(K-1) - ZSOIL(K+1))
cv          DSMDZ2 = (SH2O(K) - SH2O(K+1)) / (DENOM * 0.5)

cvk 8/10 Change similar to the first soil layer above  
cfews
          if(rdst .eq. 0) then
           smcx=smc(k)+(smc(k+1)-smc(k))*(zsoil(k)-zsoil(k-1))/
     +                                   (zsoil(k+1)-zsoil(k))
           dsmdz2=(smc(k)-smcx)/(DENOM * 0.5)
          else            
           dsmdz2=(smc(k)-smc(k+1))/(DENOM * 0.5)
          endif 

C ----------------------------------------------------------------------
C CALC THE MATRIX COEF, CI, AFTER CALC'NG ITS PARTIAL PRODUCT
C ----------------------------------------------------------------------
          DDZ2 = 2.0 / DENOM
          CI(K) = -WDF2 * DDZ2 / DENOM2
        ELSE
         
         SLOPX=0.
C ----------------------------------------------------------------------
C RETRIEVE THE SOIL WATER DIFFUSIVITY AND HYDRAULIC CONDUCTIVITY FOR
C THIS LAYER
C ----------------------------------------------------------------------
          CALL WDFCND (WDF2,WCND2,SH2OA(NSOIL),SMCMAX,BEXP,DKSAT,DWSAT,
     &                 SICEMAX)

C ----------------------------------------------------------------------
C CALC A PARTIAL PRODUCT FOR LATER USE IN CALC'NG RHSTT
C ----------------------------------------------------------------------
          DSMDZ2 = 0.0
C ----------------------------------------------------------------------
C SET MATRIX COEF CI TO ZERO
C ----------------------------------------------------------------------
          CI(K) = 0.0
        ENDIF

C ----------------------------------------------------------------------
C CALC RHSTT FOR THIS LAYER AFTER CALC'NG ITS NUMERATOR
C ----------------------------------------------------------------------
        NUMER = (WDF2 * DSMDZ2) + SLOPX * WCND2 - (WDF * DSMDZ)
     &    - WCND + ET(K)
        RHSTT(K) = NUMER / (-DENOM2)

C ----------------------------------------------------------------------
C CALC MATRIX COEFS, AI, AND BI FOR THIS LAYER
C ----------------------------------------------------------------------
        AI(K) = -WDF * DDZ / DENOM2
        BI(K) = -( AI(K) + CI(K) )

C ----------------------------------------------------------------------
C RESET VALUES OF WDF, WCND, DSMDZ, AND DDZ FOR LOOP TO NEXT LYR
C ----------------------------------------------------------------------
        IF (K .NE. NSOIL) THEN
          WDF = WDF2
          WCND = WCND2
          DSMDZ = DSMDZ2
          DDZ = DDZ2
        ENDIF
      END DO

C ----------------------------------------------------------------------
C END SUBROUTINE SRT
C ----------------------------------------------------------------------
      RETURN
      END
