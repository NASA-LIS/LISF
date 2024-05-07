!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
c     this routine is designed to intialize frozen ground state variables, 
c     including UZTWH, UZFWH, LZTWH, LZFSH, LZFPH, and Noah-like state 
c     variables, including SMC and SH2O. 
c     Author: Brian Cosgrove and Shugong Wang 
c     Date: 03/20/2014 
c     
      SUBROUTINE INIT_SOIL_MOIST(NSOIL, NUPL,
     +            UZTWC, UZFWC, LZTWC, LZFSC, LZFPC,
     +            RSMAX, CKSL, ZBOT, RTUP, RTLW, PSISAT, SWLT,
     +            TSOIL, ZSOIL, TBOT, BRT, SMAX, 
     +            UZTWH, UZFWH, LZTWH, LZFSH, LZFPH, SMC, SH2O)
      PARAMETER (T0=273.16)
      REAL TSOIL(5),ZSOIL(5)
      REAL SMC(5),SH2O(5)
      INTEGER NSOIL, NUPL
      REAL UZTWC, UZFWC, LZTWC, LZFSC, LZFPC
      REAL RSMAX, CKSL, ZBOT, RTUP, RTLW, PSISAT, SWLT
      REAL TBOT, BRT, SMAX
      REAL UZTWH, UZFWH, LZTWH, LZFSH, LZFPH





      SMC(1)=RSMAX*0.15
      SH2O(1)=SMC(1)

c     estimate unfrozen water storages       
      DZUP=ZSOIL(1)-ZSOIL(NUPL)
      DZLOW=ZSOIL(NUPL)-ZSOIL(NSOIL)  
      SUZ=UZTWC+UZFWC
      SLZ=LZTWC+LZFSC+LZFPC
      SMCUZ=0.001*RTUP*SUZ/DZUP+SWLT
      SMCLZ=0.001*RTLW*SLZ/DZLOW+SWLT
      TBUP=TBND_0(TSOIL(1)+T0,TSOIL(2)+T0,ZSOIL,ZBOT,1,NSOIL)
      SUP=0.
      SLW=0.

      DO I=2,NSOIL
c     calcualte average soil temperature of i-th layer 
       IF(I .NE. NSOIL) THEN
        TBDN=TBND_0(TSOIL(I)+T0,TSOIL(I+1)+T0,ZSOIL,ZBOT,
     +           I,NSOIL)
       ELSE
        TBDN=TBND_0(TSOIL(I)+T0,TBOT,ZSOIL,ZBOT,I,NSOIL)
       ENDIF
       DZ=ZSOIL(I-1)-ZSOIL(I)
       TS=ST_AVG1_0(TBUP,TSOIL(I)+T0,TBDN,DZ)
       TBUP=TBDN

c     calculate potential unfrozen water content 
       IF(I .LE. NUPL) THEN
        SMC(I)=SMCUZ
        IF(TS .LE. T0) THEN
         SH2O(I)=FRH2O_0(TS,SMC(I),SMC(I),SMAX,BRT,PSISAT,CKSL)
        ELSE
         SH2O(I)=SMC(I)
        ENDIF
        DSW=1000*(SH2O(I)-SWLT)*(ZSOIL(I-1)-ZSOIL(I))/RTUP
        IF(DSW .GT. 0.) SUP=SUP+DSW
       ELSE
        SMC(I)=SMCLZ
        IF(TS .LE. T0) THEN
         SH2O(I)=FRH2O_0(TS,SMC(I),SMC(I),SMAX,BRT,PSISAT,CKSL)
        ELSE
         SH2O(I)=SMC(I)
        ENDIF
        DSW=1000*(SH2O(I)-SWLT)*(ZSOIL(I-1)-ZSOIL(I))/RTLW
        IF(DSW .GT. 0.) SLW=SLW+DSW
       ENDIF
      ENDDO

c     initialize frozen states
      IF(SUP .GT. SUZ) SUP=SUZ
      IF(SLW .GT. SLZ) SLW=SLZ
      ALP=UZTWC/SUZ
      UZTWH=SUP*ALP
      UZFWH=SUP*(1-ALP)
      ALP=LZTWC/SLZ
      LZTWH=SLW*ALP
      ALP1=LZFSC/SLZ
      LZFSH=SLW*ALP1
      LZFPH=SLW*(1-ALP-ALP1)
      END

      FUNCTION ST_AVG1_0(TUP,TM,TDN,DZ)
      PARAMETER (T0=273.16)
      DZH=DZ*0.5
      IF(TUP .LT. T0) THEN
        IF(TM .LT. T0) THEN
          IF(TDN .LT. T0) THEN
            TAVG=(TUP+2*TM+TDN)/4.
            GOTO 777
          ELSE
            X0=(T0-TM)*DZH/(TDN-TM)
            TAVG=0.5*(TUP*DZH+TM*(DZH+X0)+T0*(2.*DZH-X0))/DZ
            GOTO 777
          ENDIF
        ELSE
          IF(TDN .LT. T0) THEN
            XUP=(T0-TUP)*DZH/(TM-TUP)
            XDN=DZH-(T0-TM)*DZH/(TDN-TM)
            TAVG=0.5*(TUP*XUP+T0*(2.*DZ-XUP-XDN)+TDN*XDN)/DZ
            GOTO 777
          ELSE
            XUP=(T0-TUP)*DZH/(TM-TUP)
            TAVG=0.5*(TUP*XUP+T0*(2.*DZ-XUP))/DZ
            GOTO 777
          ENDIF
        ENDIF
      ELSE
        IF(TM .LT. T0) THEN
          IF(TDN .LT. T0) THEN
            XUP=DZH-(T0-TUP)*DZH/(TM-TUP)
            TAVG=0.5*(T0*(DZ-XUP)+TM*(DZH+XUP)+TDN*DZH)/DZ
            GOTO 777
          ELSE
            XUP=DZH-(T0-TUP)*DZH/(TM-TUP)
            XDN=(T0-TM)*DZH/(TDN-TM)
            TAVG=0.5*(T0*(2.*DZ-XUP-XDN)+TM*(XUP+XDN))/DZ
            GOTO 777
          ENDIF
        ELSE
          IF(TDN .LT. T0) THEN
            XDN=DZH-(T0-TM)*DZH/(TDN-TM)
            TAVG=(T0*(DZ-XDN)+0.5*(T0+TDN)*XDN)/DZ
            GOTO 777
          ELSE
            TAVG=(TUP+2.*TM+TDN)/4.
            GOTO 777
          ENDIF
        ENDIF
      ENDIF
777   ST_AVG1_0=TAVG
      RETURN
      END
      
      FUNCTION TBND_0 (TU, TB, ZSOIL, ZBOT, K, NSOIL)
      PARAMETER (T0=273.16)
      REAL ZSOIL (*)
      IF(K .EQ. 1) THEN
        ZUP=0.
      ELSE
        ZUP=ZSOIL(K-1)
      ENDIF
      IF(K .EQ. NSOIL) THEN
        ZB=2.*ZBOT-ZSOIL(K)
      ELSE
        ZB=ZSOIL(K+1)
      ENDIF
      TBND_0=TU+(TB-TU)*(ZUP-ZSOIL(K))/(ZUP-ZB)
      RETURN
      END
      
      FUNCTION FRH2O_0(T,SMC,SH2O,SMCMAX,B,PSIS,CK)
      PARAMETER (HLICE=3.335 E5)      
      PARAMETER (GS = 9.81)
      PARAMETER (DICE=920.0)
      PARAMETER (DH2O=1000.0)
      PARAMETER (T0=273.16)
      PARAMETER (ERROR=0.005)
      BX = B
      IF ( B .GT. 5.5 ) BX = 5.5
      TR=T
      IF(TR .GE. T0) THEN
        FRH2O_0=SMC
        GOTO 77
      ENDIF  
      SWL=(SMC-SH2O)*100.
      SMCM=SMCMAX*100.
      SMCC=SMC*100.                       
      N=0
   10 DF = ( -PSIS*GS/HLICE ) * ( ( 1+CK*SWL*0.01 )**2 ) *
     +     ( SMCM/(SMCC-SWL ) )**BX
      DENOM = 2 * CK/( 1+CK*SWL*0.01 )+BX/(SMCC - SWL )
      SWLK = SWL - (1-(TR-T0)/(TR*DF))/DENOM
      IF(SWLK .LE. 0.0 .OR. N .GT. 100) THEN
       FK=(((-HLICE/(GS*PSIS))*((TR-T0)/TR))**(-1/B))*SMCMAX
       FRH2O_0 = AMIN1 ( FK, SMC )
       IF(FRH2O_0 .LT. 0.) FRH2O_0 = 0.02
      ELSE         
       DSWL=ABS(SWLK-SWL)
       SWL=SWLK
       IF(SWL.GE.SMCC) THEN
        SWL=SMCC*0.8
        SMCC=SMCC*0.97
       ENDIF
       N=N+1 
      IF ( DSWL .GT. ERROR ) GOTO 10
       FRH2O_0 = AMAX1 ( (SMC - SWL*0.01 ), 0.02 )
      ENDIF      
   77 CONTINUE
      RETURN
      END     
