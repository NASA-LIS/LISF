      SUBROUTINE SSTEP(SH2OOUT,SH2OIN,CMC,RHSTT,RHSCT,DT,
     &                  NSOIL,SMCMAX,CMCMAX,ZSOIL,SMC,SICE,
     &                  SMCDRY,AI,BI,CI,wplus,nups,etles)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE SSTEP
C ----------------------------------------------------------------------
C CALCULATE/UPDATE SOIL MOISTURE CONTENT VALUES AND CANOPY MOISTURE
C CONTENT VALUES. DT in sec.
C ----------------------------------------------------------------------
      INTEGER NSOLD
      PARAMETER(NSOLD = 20)

      integer nups
      INTEGER I
      INTEGER K 
      INTEGER KK11
      INTEGER NSOIL

      REAL AI(NSOLD)
      REAL BI(NSOLD)
      REAL CI(NSOLD)
      REAL CIin(NSOLD)
      REAL CMC
      REAL CMCMAX
      REAL DDZ
      REAL DT
      REAL RHSCT
      REAL RHSTT(NSOIL)
      REAL RHSTTin(NSOIL)
      REAL SH2OIN(NSOIL)
      REAL SH2OOUT(NSOIL)
      REAL SICE(NSOIL)
      REAL SMC(NSOIL)
      REAL SMCMAX
      REAL SMCDRY
      REAL STOT
      REAL WPLUS
      REAL ZSOIL(NSOIL)
      real etles, diff

C ----------------------------------------------------------------------
C CREATE 'AMOUNT' VALUES OF VARIABLES TO BE INPUT TO THE
C TRI-DIAGONAL MATRIX ROUTINE.
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        RHSTT(K) = RHSTT(K) * DT
        AI(K) = AI(K) * DT
        BI(K) = 1. + BI(K) * DT
        CI(K) = CI(K) * DT
      END DO

C ----------------------------------------------------------------------
C COPY VALUES FOR INPUT VARIABLES BEFORE CALL TO ROSR12
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        RHSTTin(K) = RHSTT(K)
      END DO
      DO K = 1,NSOLD
        CIin(K) = CI(K)
      END DO

C ----------------------------------------------------------------------
C CALL ROSR12 TO SOLVE THE TRI-DIAGONAL MATRIX
C ----------------------------------------------------------------------
      CALL ROSR12 (CI,AI,BI,CIin,RHSTTin,RHSTT,NSOIL)

C ----------------------------------------------------------------------
C SUM THE PREVIOUS SMC VALUE AND THE MATRIX SOLUTION TO GET A
C NEW VALUE.  MIN ALLOWABLE VALUE OF SMC WILL BE 0.02.
C RUNOFF3: RUNOFF WITHIN SOIL LAYERS
C ----------------------------------------------------------------------
 
      diff=0.0
      etles=0.0
      WPLUS = 0.0
      DDZ = -ZSOIL(1)
      DO K = 1,NSOIL
        IF (K .NE. 1) DDZ = ZSOIL(K - 1) - ZSOIL(K)
        SH2OOUT(K) = SH2OIN(K) + CI(K) + (diff+WPLUS) / DDZ

cvk 4/2010 check if sh2oout is less than smcdry
        if(k .le. nups) then
         diff= (SH2OOUT(K)-smcdry)*ddz
         if(diff .lt. 0.0) then
          SH2OOUT(K)=smcdry
         else
          diff=0.0
         endif
         etles=diff*1000.
        endif   

        STOT = SH2OOUT(K) + SICE(K)
        IF (STOT .GT. SMCMAX) THEN
          IF (K .EQ. 1) THEN
            DDZ = -ZSOIL(1)
          ELSE
            KK11 = K - 1
            DDZ = -ZSOIL(K) + ZSOIL(KK11)
          ENDIF
          WPLUS = (STOT-SMCMAX) * DDZ
        ELSE
          WPLUS = 0.
        ENDIF
cvk 5/2009  changed threshold from 0.02 to SMSDRY
cvk 5/2009        SMC(K) = MAX ( MIN(STOT,SMCMAX),0.02 )
cvk 5/2009        SH2OOUT(K) = MAX((SMC(K)-SICE(K)),0.02)
        SMC(K) = MAX ( MIN(STOT,SMCMAX),SMCDRY )
        SH2OOUT(K) = MAX((SMC(K)-SICE(K)),SMCDRY)
      END DO
      if(WPLUS .gt. 0.0) then
c       write(*,*) ' WARNING: Excess runoff from soil column',WPLUS
      endif 

C ----------------------------------------------------------------------
C UPDATE CANOPY WATER CONTENT/INTERCEPTION (CMC).  CONVERT RHSCT TO 
C AN 'AMOUNT' VALUE AND ADD TO PREVIOUS CMC VALUE TO GET NEW CMC.
C ----------------------------------------------------------------------
cvk there is no need for DT because in SAC all variables in mm/dt
cvk      CMC = CMC + DT * RHSCT
      CMC = CMC + RHSCT
      IF (CMC .LT. 1.E-20) CMC=0.0
      CMC = MIN(CMC,CMCMAX)

C ----------------------------------------------------------------------
C END SUBROUTINE SSTEP
C ----------------------------------------------------------------------
      RETURN
      END
