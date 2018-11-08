
C  FUNCTION CALCULATES FROZEN UPPER AND LOWER BOUNDS 

      SUBROUTINE FRZIND1(SMC,SH2O,TSOIL,ZSOIL,NSOIL,SWLT,FRZDUP,FRZDBT)

#include "config.h"
      
      REAL SMC(*),SH2O(*),ZSOIL(*),TSOIL(*)
      PARAMETER (SICEMN=0.01)
               
      frzdup=0.
      DO II=NSOIL,2,-1
       DICE=SMC(II)-SH2O(II)
       IF(DICE .GT. SICEMN) GOTO 10
      ENDDO
      FRZDBT=0.
      GOTO 21
   10 DZICE=DICE*(ZSOIL(II-1)-ZSOIL(II))/(SMC(II)-SWLT)
      FRZDBT=100.*(DZICE-ZSOIL(II-1))

cc   11 DO JJ=2,NSOIL
cc       DICE=SMC(JJ)-SH2O(JJ)
cc       IF(DICE .GT. SICEMN) GOTO 20
cc      ENDDO  
cc      FRZDUP=0.
cc      GOTO 21
cc   20 IF(JJ .LT. II .AND. JJ .GT. 2) THEN
cc       DZICE=DICE*(ZSOIL(JJ-1)-ZSOIL(JJ))/(SMC(JJ)-SWLT)
cc       FRZDUP=100.*(-DZICE-ZSOIL(JJ))
cc      ELSE
cc       FRZDUP=0.
cc      ENDIF
cc   21 CONTINUE

   11 do jj=2,nsoil 
cc       if(jj .eq. 2) then
       if(tsoil(jj-1) .le. 0.) goto 21
       if(tsoil(jj-1) .gt. 0. .and. tsoil(jj) .le. 0.) then
        if(jj .eq. 2) then
         dz=0.
        else
         dz=zsoil(jj-2)
        endif  
        frzdup=100.*0.5*((dz-zsoil(jj))*tsoil(jj-1)/
     +         (tsoil(jj-1)-tsoil(jj))-(dz+zsoil(jj-1)))
        if(frzdup .gt. frzdbt) frzdup=frzdbt
        goto 21
       endif
      enddo  
cc      endif
Cc    to be consistent with definition of missing values in
C     frzFuncBefore.cpp
      frzdup=-1.                   
   21 continue

      RETURN
      END
