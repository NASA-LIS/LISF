C  SUBROUTINE CALLS ESTFST1 TO GENERATE SAC-FROZEN GROUND STATES
c  VK 08/2012 scaling states option added
c  added nsoil_old (soil layer # from old parameters) to pass to estfst1

      SUBROUTINE FSTFG1(NX,NSOIL,NUPL,NPIX,FGCO,STATE,TSOIL,ZSOIL,FGPM,
     +           TBOT,BRT,SMAX,SACST_PRV,SMC,SH2O,DTHR,DTFRZ)

#include "config.h"
      
      REAL FGCO(5),TSOIL(5),FGPM(7),STATE(5),ZSOIL(5)
      REAL SMC(5),SH2O(5),SACST_PRV(5)
      
C--      COMMON/FRDSTFG/SMAX,PSISAT,BRT,SWLT,QUARTZ,STYPE,NUPL,
C--     +               NSAC,RTUP,RTLW,DZUP,DZLOW
cc      COMMON/FRZCNST/ FRST_FACT,CKSOIL,ZBOT
C--      COMMON/FRZCNST/ FRST_FACT,ZBOT      

C--      IF(IVERS .GT. 1) THEN
CVK  ESTIMATE INITIAL SAC-SMA UNFROZEN WATER STORAGES
CVK  IF THEY ARE NOT INPUTED
C--       IF(FGCO(1) .EQ. 999.) THEN

      CALL ESTFST1(NX,NSOIL,NUPL,NPIX,FGCO,STATE,TSOIL,ZSOIL,FGPM,
     +        TBOT,BRT,SMAX,SACST_PRV,SMC,SH2O,DTHR,DTFRZ)

C--       ENDIF   
           
CVK  CHECKING CONSISTENCY OF UNFROZEN WATER STATES
ccc       L=0
ccc       SFI=0.
       DO I=1,5
        IF(FGCO(I) .EQ. -999. .OR. FGCO(I) .GT. STATE(I)) THEN
cvk 2011 change warning message
         delta=FGCO(I)-STATE(I)
         if(delta .gt. 0.0001) then
#if defined(CHPS)
         WRITE (MESSAGESTRING,921) I,FGCO(I),STATE(I)
        call logfromfortran( WARNING_LEVEL, MESSAGESTRING )
#else
             WRITE (*,921) I,FGCO(I),STATE(I)
#endif
         endif
cvk         WRITE (*,921) I,FGCO(I)
            FGCO(I) = STATE(I)
cvk         WRITE (*,922) STATE(I)
ccc         L=1
        ENDIF 
ccc        SFI=SFI+STATE(I)-FGCO(I)
       ENDDO
  921 FORMAT (1H ,49H WARNING* IMPOSSIBLE VALUE OF LIQUID STATE NUMBER,
     1 I2,2F14.7)

  922 FORMAT (1H ,1X,10HCHANGED TO,3X,6F14.8)
          
C--       xx1=ABS(FGCO(1)*FRST_FACT)
C--       xx2=1.1*SFI       
C--       IF(xx1 .GT. xx2) THEN
C--        FGCO(1)=-SFI/FRST_FACT
C--        L=1
C--       ENDIF 
C--       xx1=ABS(FGCO(1)*FRST_FACT)
C--       xx2=0.9*SFI
C--       IF(xx1 .LT. xx2) THEN
C--        FGCO(1)=-SFI/FRST_FACT
C--        L=1
C--       ENDIF

ccc       IF (L.NE.0) THEN
ccc        WRITE (*,920) 
ccc  920 FORMAT (1H0,25X,25H**WARNING FROM VERSION **)
ccc        WRITE (*,921)
ccc  921 FORMAT (1H ,10X,43HINITIAL STATE VARIABLES CONTAIN IMPOSSIBLE ,
ccc     1 7HVALUES.,5X,5HUZTWH,2X,5HUZFWH,2X,5HLZTWH,2X,5HLZFSH,2X,5HLZFPH,
ccc     2 2X,5HFINDX)
ccc        WRITE (*,922) (FGCO(I),I=1,5)
ccc  922 FORMAT (1H ,50X,10HCHANGED TO,3X,6F7.0)
ccc       ENDIF

C--      ENDIF  
      
      RETURN
      END
