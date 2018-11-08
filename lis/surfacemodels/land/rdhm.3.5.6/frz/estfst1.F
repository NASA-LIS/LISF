c  VK 08/2012 scaling states option added
c  subroutine changes slightly: added nsoil_old (soil layer # from old parameters)
c  if nsoil_old .ge. 1 (scaling case) it estimates without regard missing values
c
      SUBROUTINE ESTFST1(NX,NSOIL,NUPL,NPIX,FGCO,STATE,TSOIL,ZSOIL,FGPM,
     +            TBOT,BRT,SMAX,SACST_PRV,SMC,SH2O,DTHR,DTFRZ)

#include "config.h"

      PARAMETER (T0=273.16)
C  DTHR - SAC-SMA TIME STEP IN HOURS 
C  DTSIM - SIMULATION TIME STEP OF FROZEN GROUND MODEL, IN SEC.
C  IT CAN BE LESS THAN SAC-SMA TIME STEP, BUT NOT GREATER.
      PARAMETER (DTSIM = 1800.)

C  FGPM ARRAY: 1-RSMAX, 2-CKSL, 3-ZBOT, 4-RTUP, 5-RTLW, 6-PSISAT, 7-SWLT
      REAL FGCO(*),TSOIL(*),FGPM(*),STATE(*),ZSOIL(*)
      REAL SMC(*),SH2O(*),SACST_PRV(*)
      INTEGER miss

C  DEFINE FROZEN GROUND SIMULATION TIME STEP AS A FRACTION OF SAC DTHR
C  WHICH IS IN HOURS
       RDT=DTHR*3600/DTSIM
       IDT=RDT+0.5
       DTFRZ=DTSIM
       IF(RDT .LT. 1.) THEN
        DTFRZ=DTHR*3600
       ELSE
        IF(MOD(DTHR*3600,DTFRZ) .NE. 0.) DTFRZ=DTHR*3600/IDT
       ENDIF
cc-       WRITE(*,*) ' SAC-SMA TIME STEP=',DTHR/24,'days '
cc-       WRITE(*,*) ' FROZEN GROUND TIME STEP=',DTFRZ/3600,'hr'
cc-       WRITE(*,*) ' NUMBER FROZEN GROUND LOOPS=',IDT      

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c VK 05/08 added explicitly soil moisture state grids: smc and sh20
c
c 3/10 surface water changes: smc(1) & sh2o(1) separated from others
C  DEFAULT MAX POROSITY OF THE RESIDUE LAYER IS 0.58
      if(smc(1) .lt. 0. ) then
       SMC(1)=FGPM(1)*0.15
      endif

      if( sh2o(1) .lt. 0.) then
       SH2O(1)=SMC(1)
      endif         

      miss=0

c VK 08/2012  scaling states option added  
      do i=2,nsoil
       if(smc(i) .lt. 0. .or. sh2o(i) .lt. 0.) miss=miss+1
      enddo

c VK 9/2012 If fill missing option, regenerate states even they available  
      if(nx == -9) miss=1
c
c VK 05/08 check if there are missing SM states
c          calculate them using SAC states and soil temperature
c VK 08/2012  scaling states option added 
      if(miss .gt. 0) then
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cvk 12/08        write(*,*) 'no frozen ground states: be estimated',miss      
C  DEFAULT MAX POROSITY OF THE RESIDUE LAYER IS 0.58
c 3/10       SMC(1)=FGPM(1)*0.15
c 3/10       SH2O(1)=SMC(1)
       
C  ESTIMATE UNFROZEN WATER STORAGES
       DZUP=ZSOIL(1)-ZSOIL(NUPL)
       DZLOW=ZSOIL(NUPL)-ZSOIL(NSOIL)  
       SUZ=STATE(1)+STATE(2)
       SLZ=STATE(3)+STATE(4)+STATE(5)
       SMCUZ=0.001*FGPM(4)*SUZ/DZUP+FGPM(7)
       SMCLZ=0.001*FGPM(5)*SLZ/DZLOW+FGPM(7)
       TBUP=TBND(TSOIL(1)+T0,TSOIL(2)+T0,ZSOIL,FGPM(3),1,NSOIL)
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
        IF(I .LE. NUPL) THEN
         SMC(I)=SMCUZ
         IF(TS .LE. T0) THEN
          SH2O(I)=FRH2O_356(TS,SMC(I),SMC(I),SMAX,BRT,FGPM(6),FGPM(2))
         ELSE
          SH2O(I)=SMC(I)
         ENDIF
         DSW=1000*(SH2O(I)-FGPM(7))*(ZSOIL(I-1)-ZSOIL(I))/FGPM(4)
         IF(DSW .GT. 0.) SUP=SUP+DSW        
        ELSE
         SMC(I)=SMCLZ
         IF(TS .LE. T0) THEN
          SH2O(I)=FRH2O_356(TS,SMC(I),SMC(I),SMAX,BRT,FGPM(6),FGPM(2))
         ELSE
          SH2O(I)=SMC(I)
         ENDIF
         DSW=1000*(SH2O(I)-FGPM(7))*(ZSOIL(I-1)-ZSOIL(I))/FGPM(5)
         IF(DSW .GT. 0.) SLW=SLW+DSW
        ENDIF 
       ENDDO
       l=0
       do i=1,nsoil
        if(tsoil(i) .lt. 0.) l=1
       enddo

c VK 9/2012 If fill missing option, regenerate states even they available
      if(nx == -9) then
       l=1
       nx=1
      endif

       IF(NX .GE. 1) THEN
        if(l .eq. 1) then
C  FROZEN STATES NOT AVAILABLE
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
        else
         do i=1,5
          fgco(i)=state(i)
         enddo
        endif   
       ENDIF   
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c VK 05/08 check if SAC-FRZ states available
      else
       if(nx .ge. 1) then
#ifdef CHPS
        write(MESSAGESTRING,*) 
     +             'ERROR: there are SM states but no SAC-FRZ states'
        write(MESSAGESTRING,*) (fgco(i),i=1,5)
        call logfromfortran( FATAL_LEVEL, MESSAGESTRING )
#else
        write(*,*) 'ERROR: there are SM states but no SAC-FRZ states'
        write(*,*) (fgco(i),i=1,5)
        stop
#endif
       endif
C       write(*,*) 'there are states',miss
      endif
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
C  SAVE PREVIOUS TOTAL SOIL MOISTURE STATES
      DO I=1,5
       SACST_PRV(I) = STATE(I)
      ENDDO

cc-       WRITE(*,*) '    INITIAL SOIL LAYER STATES: '
cc-       WRITE(*,'(A20,7F7.3)') ' LAYER DEPTH (M):   ',
cc-     +      ZSOIL(1),(ZSOIL(1)-ZSOIL(I),I=2,NSOIL)
cc-       WRITE(*,'(A20,7F7.2)') ' TOTAL MOISTURE:    ',
cc-     +                           (SMC(I),I=1,NSOIL)
cc-       WRITE(*,'(A20,7F7.2)') ' UNFROZEN MOISTURE: ',
cc-     +                           (SH2O(I),I=1,NSOIL)
cc-       WRITE(*,'(A20,5F7.2)') ' TEMPERATURE:       ',
cc-     +                          (TSOIL(I),I=1,NSOIL)       
      RETURN
      END

