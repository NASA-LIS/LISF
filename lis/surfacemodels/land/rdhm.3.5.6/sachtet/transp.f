c      program test
c      
c      parameter (nsoil=4)
c      parameter (nup=2)
c      parameter (nroot=4)
c      
c      integer nsoil, nroot
c      real et1(nsoil),smc(nsoil),zsoil(nsoil),rtdis(nsoil)
c      real etp1,cmc,shdfac,smcwlt,cmcmax,pc,cfactr,smcref,sfctmp,q2
c      
ccvk potential evaporation in kg/(m**2*s)
ccvk if PE=12 mm/day=0.5 mm/hr, it's equal 1.39*10**-4
c      etp1=0.000139
c      cmc=0.0
c      shdfac=1.0
c      smcwlt=0.1
ccvk cmcmax value of canopy retention in m; cfactr=0.5 is power
c      cmcmax=0.0005
c      cfactr=0.5
c      pc=0.8
ccvk field capacity
c      smcref=0.35
c      sfctmp=293
ccvk mixing ratio in kg/kg; not in use in this subroutine
c      q2=0.0
c      
ccvk zsoil array not in use in this version
c      zsoil(1)=-0.1
c      zsoil(2)=-0.3
c      zsoil(3)=-0.8
c      zsoil(4)=-1.5
c      
c      smc(1)=0.30
c      smc(2)=0.30
c      smc(3)=0.30
c      smc(4)=0.30
c
c      rtdis(1)=0.25
c      rtdis(2)=0.25
c      rtdis(3)=0.25
c      rtdis(4)=0.25
c      
c      call TRANSP (ET1,NSOIL,ETP1,SMC,CMC,ZSOIL,SHDFAC,SMCWLT,
c     &                   CMCMAX,PC,CFACTR,SMCREF,SFCTMP,Q2,NROOT,RTDIS)
c     
c      sup=0.
c      slo=0.0
c      do i=1,nup
c       sup=sup+et1(i)*3600
c      enddo
c      do i=nup+1,nsoil
c       slo=slo+et1(i)*3600
c      enddo 
c        
c      write(*,*) pc,etp1*3600,sup+slo,sup,slo
c      write(*,*) zsoil
c      write(*,*) smc
c      write(*,*) rtdis
c      write(*,*) et1*3600
c      
c      stop
c      end

cvk      SUBROUTINE TRANSP (ET1,NSOIL,ETP1,SMC,CMC,ZSOIL,SHDFAC,SMCWLT,
cvk     &                   CMCMAX,PC,CFACTR,SMCREF,SFCTMP,Q2,NROOT,RTDIS)

cvk removed variables not in use: ZSOIL, Q2, SFCTMP
cvk and added EC1 to correctly reduce ETP1

cvk 5/2009 SMCDRY instead SMCWLT      SUBROUTINE TRANSP (ET1,NSOIL,ETP1,SH2O,CMC,SHDFAC,SMCWLT,
      SUBROUTINE TRANSP (ET1,NSOIL,ETP1,SH2O,CMC,SHDFAC,SMCDRY,
     &                   CMCMAX,PC,CFACTR,SMCREF,NROOT,RTDIS,EC1)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE TRANSP
C ----------------------------------------------------------------------
C CALCULATE TRANSPIRATION FOR THE VEG CLASS.
C ----------------------------------------------------------------------
      INTEGER I
      INTEGER K
      INTEGER NSOIL
      INTEGER NROOT

      REAL CFACTR
      REAL CMC
      REAL CMCMAX
      REAL DENOM
      real EC1
      REAL ET1(NSOIL)
      REAL ETP1
      REAL ETP1A
      REAL GX (5)
C.....REAL PART(NSOIL)
      REAL PC
cvk      REAL Q2
      REAL RTDIS(NSOIL)
      REAL RTX
      REAL SFCTMP
      REAL SGX
      REAL SHDFAC
cvk      REAL SMC(NSOIL)
      REAL SH2O(NSOIL)
      REAL SMCREF
      REAL SMCDRY
cvk      REAL ZSOIL(NSOIL)

C ----------------------------------------------------------------------
C INITIALIZE PLANT TRANSP TO ZERO FOR ALL SOIL LAYERS.
C ----------------------------------------------------------------------
      DO K = 1,NSOIL
        ET1(K) = 0.
      END DO

C ----------------------------------------------------------------------
C CALCULATE AN 'ADJUSTED' POTENTIAL TRANSPIRATION
C IF STATEMENT BELOW TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
C NOTE: GX AND OTHER TERMS BELOW REDISTRIBUTE TRANSPIRATION BY LAYER,
C ET(K), AS A FUNCTION OF SOIL MOISTURE AVAILABILITY, WHILE PRESERVING
C TOTAL ETP1A.
C ----------------------------------------------------------------------
cvk....................................................................
cvk changed next two statements to correctly recalculate ETP1
cvk      IF (CMC .NE. 0.0) THEN
cvk        ETP1A = SHDFAC * PC * ETP1 * (1.0 - (CMC /CMCMAX) ** CFACTR)
cvk....................................................................
      IF (EC1 .NE. 0.0 .AND. ETP1 .NE. 0.0) THEN
        ETP1A = SHDFAC * PC * ETP1 * (1.0 - EC1/(SHDFAC*ETP1))
      ELSE
        ETP1A = SHDFAC * PC * ETP1
      ENDIF
      
      SGX = 0.0
      DO I = 1,NROOT
cvk 5/2009 SMCDRY instead SMCWLT GX(I) = ( SH2O(I) - SMCWLT ) / ( SMCREF - SMCWLT )
        GX(I) = ( SH2O(I) - SMCDRY ) / ( SMCREF - SMCDRY )
        GX(I) = MAX ( MIN ( GX(I), 1. ), 0. )
        SGX = SGX + GX (I)
      END DO
      SGX = SGX / NROOT
      
      DENOM = 0.
      DO I = 1,NROOT
        RTX = RTDIS(I) + GX(I) - SGX
        GX(I) = GX(I) * MAX ( RTX, 0. )
        DENOM = DENOM + GX(I)
      END DO
      IF (DENOM .LE. 0.0) DENOM = 1.

      DO I = 1,NROOT
        ET1(I) = ETP1A * GX(I) / DENOM
      END DO

C ----------------------------------------------------------------------
C ABOVE CODE ASSUMES A VERTICALLY UNIFORM ROOT DISTRIBUTION
C CODE BELOW TESTS A VARIABLE ROOT DISTRIBUTION
C ----------------------------------------------------------------------
C      ET(1) = ( ZSOIL(1) / ZSOIL(NROOT) ) * GX * ETP1A
C      ET(1) = ( ZSOIL(1) / ZSOIL(NROOT) ) * ETP1A
C ----------------------------------------------------------------------
C USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
C ----------------------------------------------------------------------
C      ET(1) = RTDIS(1) * ETP1A
C      ET(1) = ETP1A * PART(1)
C ----------------------------------------------------------------------
C LOOP DOWN THRU THE SOIL LAYERS REPEATING THE OPERATION ABOVE,
C BUT USING THE THICKNESS OF THE SOIL LAYER (RATHER THAN THE
C ABSOLUTE DEPTH OF EACH LAYER) IN THE FINAL CALCULATION.
C ----------------------------------------------------------------------
C      DO K = 2,NROOT
C        GX = ( SH2O(K) - SMCWLT ) / ( SMCREF - SMCWLT )
C        GX = MAX ( MIN ( GX, 1. ), 0. )
C TEST CANOPY RESISTANCE
C        GX = 1.0
C        ET(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT))*GX*ETP1A
C        ET(K) = ((ZSOIL(K)-ZSOIL(K-1))/ZSOIL(NROOT))*ETP1A
C ----------------------------------------------------------------------
C USING ROOT DISTRIBUTION AS WEIGHTING FACTOR
C ----------------------------------------------------------------------
C        ET(K) = RTDIS(K) * ETP1A
C        ET(K) = ETP1A*PART(K)
C      END DO      
C ----------------------------------------------------------------------
C END SUBROUTINE TRANSP
C ----------------------------------------------------------------------
      RETURN
      END
