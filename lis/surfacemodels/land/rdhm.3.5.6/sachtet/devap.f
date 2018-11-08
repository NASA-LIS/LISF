cvk original Noah code
cvk      SUBROUTINE DEVAP (EDIR1,ETP1,SMC,ZSOIL,SHDFAC,SMCMAX,BEXP,
cvk     &                DKSAT,DWSAT,SMCDRY,SMCREF,SMCWLT,FXEXP)

      SUBROUTINE DEVAP (EDIR1,ETP1,SH2O,SHDFAC,SMCMAX,SMCDRY,
     &                  SMCREF,SMCWLT,FXEXP,uztwh,uztwm,bareadj)

      IMPLICIT NONE

C ----------------------------------------------------------------------
C SUBROUTINE DEVAP
C ----------------------------------------------------------------------
C CALCULATE DIRECT SOIL EVAPORATION
C ----------------------------------------------------------------------

cvk  Direct evaporation option: if FXEXP=-1, Chen option; other value, Ek
cvk  Chen option: Mahfouf and Noilhan (1991) method, (S-Swlt)/(Sfld-Swlt)
cvk  M. Ek option: non-linear ((S-Swlt)/(smax-Swlt))**FXEXP
cvk                                     FXEXP usually = 2.0 

cfews  excluded options: just Ek-Chen switch
cfews      real xfact,uztwh,uztwm,bareadj,bareadjx,fxexpx
      real uztwh,uztwm,bareadj      
      
C      REAL DEVAP
      REAL EDIR1
      REAL ETP1
      REAL FX
      REAL FXEXP
      REAL SHDFAC
      REAL SH2O
      REAL SMCDRY
      REAL SMCMAX
      REAL SMCREF
      REAL SMCWLT
      REAL SRATIO

cvk 5/2010  added option to switch Chen-Ek option duo to greenness
cvk
cfews      if(bareadj .lt. 0.0) then
cvk  selected switch option
cfews       bareadjx=-bareadj
cfews       if(shdfac .gt. bareadjx) then

       if(shdfac .gt. bareadj) then       
cvk  switch to Chen option
cfews        if(FXEXP .le. 0.0) then
cfews         fxexpx=fxexp
cfews        else
cfews         fxexpx=0.0
cfews        endif  
cvk 3/10  new Chen option (VK extention) 
        SRATIO = (SH2O - SMCDRY) / (SMCREF - SMCDRY)
cfews        if(sratio .lt. 1.0) then
cfews         sratio=sratio*((SMCREF - SMCDRY)/(SMCMAX - SMCDRY))**(-FXEXPX)
cfews        else
cfews         sratio= ((SH2O - SMCDRY)/(SMCMAX - SMCDRY))**(-FXEXPX)
cfews        endif
        if(sratio .gt. 1.0) sratio = 1.0
        EDIR1 = sratio * ( 1.0 - SHDFAC ) * ETP1

       else 
cvk  switch to Ek option
cfews        if(FXEXP .le. 0.0) then
cfews         fxexpx=2.0
cfews        else
cfews         fxexpx=fxexp
cfews        endif
        SRATIO = (SH2O - SMCDRY) / (SMCMAX - SMCDRY)
        IF (SRATIO .GT. 0.) THEN
cfews          FX = SRATIO**FXEXPX
          FX = SRATIO**FXEXP          
          FX = MAX ( MIN ( FX, 1. ) ,0. )
        ELSE
          FX = 0.
        ENDIF
        EDIR1 = FX * ( 1.0 - SHDFAC ) * ETP1         
       endif
cfews      else
       
c  not switch option; use one defined in input card option
cfews       if(fxexp .le. 0.0) then
c  Chen option only
cfews        SRATIO = (SH2O - SMCDRY) / (SMCREF - SMCDRY)
cfews        if(sratio .lt. 1.0) then
cfews         sratio=sratio*((SMCREF - SMCDRY)/(SMCMAX - SMCDRY))**(-FXEXP)
cfews        else
cfews         sratio= ((SH2O - SMCDRY)/(SMCMAX - SMCDRY))**(-FXEXP)
cfews        endif
cfews        if(sratio .gt. 1.0) sratio = 1.0
cfews        EDIR1 = sratio * ( 1.0 - SHDFAC ) * ETP1
cfews       else         
cfews        if(fxexp .ge. 10.) then
c  SAC option only
cfews         sratio=uztwh/uztwm
cfews         EDIR1 = sratio * ( 1.0 - SHDFAC ) * ETP1
cfews         if(edir1 .gt. uztwh) edir1=uztwh
cfews        else       
c  Ek option only
C ----------------------------------------------------------------------
C DIRECT EVAP A FUNCTION OF RELATIVE SOIL MOISTURE AVAILABILITY, LINEAR
C WHEN FXEXP=1.
C FX > 1 REPRESENTS DEMAND CONTROL
C FX < 1 REPRESENTS FLUX CONTROL
C ----------------------------------------------------------------------
cfews         SRATIO = (SH2O - SMCDRY) / (SMCMAX - SMCDRY)
cfews         IF (SRATIO .GT. 0.) THEN
cfews           FX = SRATIO**FXEXP
cfews           FX = MAX ( MIN ( FX, 1. ) ,0. )
cfews         ELSE
cfews           FX = 0.
cfews         ENDIF

C ----------------------------------------------------------------------
C ALLOW FOR THE DIRECT-EVAP-REDUCING EFFECT OF SHADE
C ----------------------------------------------------------------------
cfews         EDIR1 = FX * ( 1.0 - SHDFAC ) * ETP1
cfews        endif 
cfews       endif
cfews      endif
C ----------------------------------------------------------------------
C END SUBROUTINE DEVAP
C ----------------------------------------------------------------------
      RETURN
      END
 
