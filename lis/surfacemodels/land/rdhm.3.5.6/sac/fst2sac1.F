C  SUBROUTINE RECALCULATES SOIL MOISTURE STATES INTO SAC-SMA STORAGES
C
      SUBROUTINE FST2SAC1(ZSOIL,SMC,SH2O,NUP,NLW,SWLT,RT,TWH,FWH,
     +                    F1WH,TWC,FWC,F1WC,FROST,PX, HRAPX,HRAPY)
#include "config.h"
      
      REAL ZSOIL(*),SMC(*),SH2O(*)

      FROST=0.         
      STOT=0.
      SLIQ=0.
      SX=0.
      SWH=TWH+FWH+F1WH+1000.*SWLT*(ZSOIL(NUP-1)-ZSOIL(NLW))/RT
      DO I=NUP,NLW

c8  Correction liquid water to be not more than swlt
       if(sh2o(i) .lt. swlt) sh2o(i) = swlt
        
       DZ=ZSOIL(I-1)-ZSOIL(I)
       SX=SX+1000.*SH2O(I)*DZ/RT
       FROST=FROST+1000.*(SMC(I)-SH2O(I))*DZ/RT
       STOT=STOT+1000.*(SMC(I)-SWLT)*DZ/RT
       SLIQ=SLIQ+1000.*(SH2O(I)-SWLT)*DZ/RT
      ENDDO
     
       xxx = abs(stot-TWC-FWC-F1WC) 
#ifdef ENABLE_DEBUGINFO
       if(xxx .gt. 0.5) then
        write(*,*) '** NO TOTAL BALANCE:',HRAPX,HRAPY
        write(*,'(5f12.6,2i2)') stot,twc,fwc,f1wc,rt,nup,nlw
        write(*,'(5f12.6)') stot,twh,fwh,f1wh,swlt
C       write(*,'(5f12.6)') (zsoil(i),i=nup-1,nlw)
        write(*,'(5f12.6)') (smc(i),i=nup,nlw)
        write(*,'(5f12.6)') (sh2o(i),i=nup,nlw)
       endif 
#endif

C  DH2O IS LIQUID WATER CHANGE DUE TO FREEZING/THAWING
      DH2O=SX-SWH
      if(frost .eq. 0. .and. twh. eq. twc .and. fwh .eq. fwc
     +                      .and. f1wh.eq.f1wc) then
       continue
      else 
       IF(FROST .EQ. 0. .AND. DH2O .EQ. 0.) THEN
C  NO FROZEN GROUND, NO LIQUID WATER CHANGE
        TWH=TWC
        FWH=FWC
        F1WH=F1WC
       ELSE 
C  CHANGE LIQUID WATER STORAGES
        SWH=TWH+FWH
        IF(DH2O .LT. 0.) THEN
C  UNFROZEN WATER REDUCTION (FREEZING)       
         IF(TWH .GE. 0.) THEN
          DSH=SWH+DH2O
          IF(DSH .GE. 0.) THEN
           IF(SWH .GT. 1E-04) THEN
c  new
            if(px .gt. 0. .and. sliq .lt. stot) then
             fwh=fwh+dh2o
             if(fwh .lt. 0.) then
              twh=twh+fwh
              fwh=0.
             endif
            else  
             ALP=DH2O*TWH/SWH
             TWH=TWH+ALP
             FWH=FWH+DH2O-ALP
            endif
           ELSE
            TWH=TWH+DH2O
           ENDIF
          ELSE
           TWH=DSH

c VK 2013 fixing lower zone balance
           IF(TWH .LT. 0.0 .AND. F1WH .NE. 0.0) THEN
            FWH=FWH+TWH
            TWH=0.0
            IF(FWH .LT. 0.0) THEN
             fwh=0.0
            ENDIF
           ENDIF
c VK 2013 END fixing lower zone balance

cc           FWH=0.
          ENDIF
         ELSE
          FWH=FWH+DH2O
          IF(FWH .LT. 0.) THEN
           TWH=TWH+FWH
           FWH=0.
          ENDIF
         ENDIF
        ELSE
C  UNFROZEN WATER INCREASE (THAWING)                      
         IF(TWH .LT. 0.) THEN
          TWH=TWH+DH2O
         ELSE
          IF(SWH .GT. 1E-04) THEN
           ALP=DH2O*TWH/SWH
           TWH=TWH+ALP
           FWH=FWH+DH2O-ALP
           IF(FWH .GT. FWC) THEN
            TWH=TWH+FWH-FWC
            FWH=FWC
           ENDIF
          ELSE
           TWH=TWH+DH2O
          ENDIF
         ENDIF
         IF(TWH .GT. TWC) THEN
          FWH=FWH+TWH-TWC
          TWH=TWC
         ENDIF
        ENDIF
       ENDIF        
      endif
      
C  CHECK CONSISTENCY OF THE WATER BALANCE BETWEEN SAC-SMA AND
C  FROZEN GROUND STATES 
      S=TWC+FWC+F1WC

cVK  April 2013 changes  -------------------------------------------
      sx=s-stot
      dshx=sx*0.001/(zsoil(nup-1)-zsoil(nlw))
      do k=nup,nlw
       sh2o(k)=sh2o(k)+dshx
      enddo

#ifdef ENABLE_DEBUGINFO
cc      IF(ABS(STOT-S) .GT. 0.5) then
cc       WRITE(*,*)
cc     + ' ** WARN ** NO TOTAL WATER BALANCE:',NUP,STOT,TWC,FWC,F1WC
cc      endif       
#endif

      S=TWH+FWH+F1WH
      sx=s-sliq
      dshx=sx*0.001/(zsoil(nup-1)-zsoil(nlw))
      do k=nup,nlw
       sh2o(k)=sh2o(k)+dshx
      enddo

#ifdef ENABLE_DEBUGINFO
cc      IF(ABS(SLIQ-S) .GT. 0.5) then
cc       WRITE(*,*)
cc     + ' ** WARN ** NO LIQUID WATER BALANCE:',NUP,SLIQ,TWH,FWH,F1WH
cc      endif 
#endif

      RETURN
      END
            
