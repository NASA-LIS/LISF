c  Estimate sh2o & frz states from smc and tsoil

      subroutine scale_sh2o(smc,tsoil,fgpm,nsoil,nupl,zsoil,
     +                         state,smax,brt,tbot,fgco,sh2o)
      real smc(nsoil),sh2o(5),tsoil(nsoil),fgpm(7)
      real state(5),zsoil(nsoil),fgco(5)
      PARAMETER (T0 = 273.16)

      TBUP=TBND(TSOIL(1)+T0,TSOIL(2)+T0,ZSOIL,FGPM(3),1,NSOIL)
      IF(SH2O(1) .LT. 0.0) SH2O(1)=SMC(1)
      SUZ=STATE(1)+STATE(2)
      SLZ=STATE(3)+STATE(4)+STATE(5)      
      SUP=0.
      SLW=0.
      DO I=2,NSOIL
C  CALCULATE AVERAGE SOIL TEMPERATURE OF I-TH LAYER
       IF(I .NE. NSOIL) THEN
        TBDN=TBND(TSOIL(I)+T0,TSOIL(I+1)+T0,ZSOIL,FGPM(3),
     +            I,NSOIL)
       ELSE
        TBDN=TBND(TSOIL(I)+T0,TBOT,ZSOIL,FGPM(3),I,NSOIL)
       ENDIF
       DZ=ZSOIL(I-1)-ZSOIL(I)
       TS=ST_AVG1(TBUP,TSOIL(I)+T0,TBDN,DZ)
       TBUP=TBDN

C  CALCULATE POTENTIAL UNFROZEN WATER CONTENT
       IF(TS .LE. T0) THEN
        SH2O(I)=FRH2O_356(TS,SMC(I),SMC(I),SMAX,BRT,FGPM(6),FGPM(2))
       ELSE
        SH2O(I)=SMC(I)
       ENDIF
       DSW=1000*(SH2O(I)-FGPM(7))*(ZSOIL(I-1)-ZSOIL(I))/FGPM(4)
       IF(DSW .GT. 0.) THEN
        IF(I .LE. NUPL)THEN
         SUP=SUP+DSW
        ELSE
         SLW=SLW+DSW
        ENDIF
       ENDIF 
      ENDDO

C  ESIMATE FROZEN SAC STATES
      L=0
      DO I=1,NSOIL
       IF(TSOIL(I) .LT. 0.0) L=1
      ENDDO
      IF(L .EQ. 1) THEN
       IF(SUP .GT. SUZ) SUP=SUZ
       IF(SLW .GT. SLZ) SLW=SLZ
       ALP=STATE(1)/SUZ
       FGCO(1)=SUP*ALP
       FGCO(2)=SUP*(1-ALP)
       ALP=STATE(3)/SLZ
       FGCO(3)=SLW*ALP
       ALP1=STATE(4)/SLZ
       FGCO(4)=SLW*ALP1
       FGCO(5)=SLW*(1-ALP-ALP1)
      ELSE
       DO I=1,5
        FGCO(I) = STATE(I)
       ENDDO
      ENDIF

      RETURN
      END


