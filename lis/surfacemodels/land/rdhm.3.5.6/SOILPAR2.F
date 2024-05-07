!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
C MEMBER SOILPAR1
C
      SUBROUTINE SOILPAR2(STXT,SAND,CLAY,SUPM,SLWM,SMAX,PSISAT,BRT,SWLT,
     +           QUARTZ,STYPE,NSOIL,NUPL,NSAC,ZSOIL,RTUP,RTLOW)

C#include "config.h"
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    PURPOSE:  CALCULATE SOIL PARAMETERS BASED ON SOIL TEXTURE (STXT),
C     AND UPPER (SUPM) AND LOWER (SLWM) ZONE STORAGES 
C     WRITTEN BY VICTOR KOREN - HRL   JANUARY 2000
C    
C     MODIFIED BY SHUGONG WANG NASA HSL  09/26/2013
C     ADD SAND AND CLAY AS INPUT 
C--------------------------------------------------------------------------
C        SMAX - SOIL POROSITY
C        PSISAT - SATURATION MATRIC POTENTIAL
C        BRT - SLOPE OF THE RETENTION CURVE
C        FRT = LOG(PSISAT)+BRT*LOG(SMAX)+2.
C        QUARTZ - FRACTION OF QUARTZ 
C        STYPE - SOIL TYPE (11 - IF COARSE SOIL, 12 - IF FINE)
C        NSOIL - TOTAL NUMBER OF SELECTED LAYERS
C        NUPL - NUMBER OF UPPER ZONE LAYERS + RESIDUE LAYER 
C        RTUP - =1, AND >1 ONLY WHEN UPPER ZONE DEEPPER THAN MAX SOIL DEPTH
C        RTLOW - THE SAME FOR LOWER ZONE 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      PARAMETER (NSOILMAX = 10)
      
C      REAL STXT, STYPE, SAND(12), CLAY(12), SDPTH(NSOILMAX)
      REAL STXT, STYPE, SAND, CLAY, SDPTH(NSOILMAX)
      REAL ZSOIL(*)
cc-      COMMON/FRZCNST/ FRST_FACT,ZBOT
      
C  DEFINITION OF DEFAULT SOIL LAYER THICKNESS
C  FIVE LAYERS ARE DIFINED. IF IT'S DESIRABLE TO USE MORE OR LESS 
C  LAYERS, SDPTH SHOULD BE CHANGED, HOWEVER, THE FIRST LAYER (VERY THIN 
C  RESIDUE LAYER) SHOULD NOT BE CHANGED 
cc      REAL SDPTH(10)/0.01, 0.11, 0.41, 0.81, 1.41, 0.,0.,0.,0.,0./
cc/was before_frz/      REAL SDPTH(10)/0.03, 0.14, 0.41, 0.80, 1.41, 0.,0.,0.,0.,0./
CC      REAL SDPTH/0.03, 0.14, 0.41, 0.83, 1.68, 0.,0.,0.,0.,0./
cc      REAL SDS(10)/0.03, 0.14, 0.41, 0.83, 1.69, 0.,0.,0.,0.,0./
      REAL SDS(10)/0.03, 0.14, 0.41, 0.81, 1.51, 0.,0.,0.,0.,0./
CC  this soil depth was only for Valdai to match 20, 40, 80 cm obs data
CC      REAL SDPTH(10)/0.03, 0.13, 0.43, 0.83, 1.43, 0.,0.,0.,0.,0./
      
C  SAND AND CLAY FRACTION TABLES DEPENDING ON SOIL TEXTURE CLASS (FROM 
C  CLASS 1 - SAND TO CLASS 12 - CLAY), FROM COSBY ET AL. [WWR, 20(6), 1984]
C    1 - SAND               5.2% (IN USA)
C    2 - LOAMY SAND         5.6%
C    3 - SANDY LOAM        18.1%
C    4 - SILTY LOAM        19.4%
C    5 - SILT        
C    6 - LOAM              12.7%
C    7 - SANDY CLAY LOAM    3.9%
C    8 - SILTY CLAY LOAM   12.0%
C    9 - CLAY LOAM          8.5%
C   10 - SANDY CLAY         0.7%
C   11 - SILTY CLAY         6.4%
C   12 - CLAY               7.5% 
C  FRACTION OF QUARTZ ASSUMED TO BE EQUAL FRACTION OF SAND AND ALL SOILS 
C  ASSUMED AS FINE, FROM PETERS-LIDARD ET AL. [J. ATM. SCI., 55, 1998]
C     
C  The following lines are commented by Shugong Wang 
C      DATA SAND/0.92,0.82,0.58,0.17,0.09,0.43,0.58,0.10,0.32,0.52,
C     +          0.06,0.22/
C      DATA CLAY/0.03,0.06,0.10,0.13,0.05,0.18,0.27,0.34,0.34,0.42,
C     +          0.47,0.58/
     
C  COSBY RELATIONSHIPS
C  	SMAX IS POROSITY
C	PSISAT IS 'SATURATION ' MATRIC POTENTIAL, IN M
C	BRT IS SLOPE OF RETENTION CURVE
C	SWLT IS WILTING POINT
C	QUARTZ IS QUARTZ FRACTION
C	STYPE IS SOIL TYPE: 11 - NATURAL COARSE, 12 - NATURAL FINE,
C			    21 - CRUSHED COARSE, 22 - CRUSHED FINE	
       
      ITXT=STXT+0.5
C      SMAX = -0.126*SAND(ITXT) + 0.489
      SMAX = -0.126*SAND + 0.489
cvk 12/05  Equation for PSISAT wrong: should be 7.74 coefficient not 0.759
cvk      PSISAT = 0.759*EXP(-3.02*SAND(ITXT))
C      PSISAT = 7.74*EXP(-3.02*SAND(ITXT))   
      PSISAT = 7.74*EXP(-3.02*SAND)   
C      BRT = 15.9*CLAY(ITXT) + 2.91
      BRT = 15.9*CLAY + 2.91
      SWLT = SMAX*(1500./PSISAT)**(-1./BRT)
C      QUARTZ = SAND(ITXT)
      QUARTZ = SAND
      STYPE = 12
      IQ=QUARTZ*100
cc-      if(ipr .gt. 0) then
cc-      WRITE(*,*) ' SOIL PROPERTIES:  NEW FROZEN GROUND VERSION'
cc-      WRITE(*,'(A,I2,A,F4.3,A,F4.3,A,I3,A)') '    TEXTURE=',ITXT,
cc-     +  ' POROSITY=',SMAX,' WILTING POINT=',SWLT,' QUARTZ=',IQ,'%'
cc-      WRITE(*,'(A,F4.2,A,F5.2)') '    SATURATED PREASURE=',
cc-     +   PSISAT,'m  BROOKS PARAMETER=',BRT
cc-      endif    

C  DEFINE SOIL LAYER THICKNESS & SOIL MOISTURE STATES OF THE LAYERS
      DO I=1,NSOILMAX
       SDPTH(I)=SDS(I)
       IF(SDPTH(I) .EQ. 0.) GOTO 20
      ENDDO
20    NSOIL=I-1
       
      ZUP=0.001*SUPM/(SMAX-SWLT)
      ZLOW=0.001*SLWM/(SMAX-SWLT)
      ZTOT=ZUP+ZLOW
      IF(ZUP .GT. (SDPTH(NSOIL)-SDPTH(1)) * 2 ) THEN  
       WRITE(*,*) 'WARNING: IN SOILPAR: TOO HIGH UPPER STORAGES:'
       WRITE(*,*) 'SUPM, ZUP:', SUPM,ZUP, 'smax swlt:',smax,swlt
C                 SDPTH(NSOIL),SDPTH(1),nsoil
C       WRITE(*,*) 'REDUCE UZTWM OR UZFWM, OR THEIR UPPER BOUNDS.'
cvk 12/08       STOP
      ENDIF  
      X=1.5*(SDPTH(NSOIL-1)-SDPTH(1))
      IF(ZUP .GT. X) SDPTH(NSOIL-1)=ZUP/1.5+SDPTH(1)

C   SPLIT THE UPPER ZONE INTO SOIL LAYERS 
      RTUP=1.
      RTLOW=1.
      ZSOIL(1)=-SDPTH(1)
cvk 1/2010      ZMAX=1.5*(SDPTH(NSOIL)-SDPTH(1))
      ZMAX=2.5*(SDPTH(NSOIL)-SDPTH(1))
      DO I=2,NSOIL
       xx1=1.5*(SDPTH(I)-SDPTH(1))
       IF(ZUP .LE. xx1 .OR. I .EQ. NSOIL) THEN
        NLOWL=I
        NUPL=I
        IF(ZUP .GT. ZMAX) THEN
         RTUP=ZMAX/ZUP
         ZUP=ZMAX
        ENDIF
        DO J=2,I
         ZSOIL(J)=ZSOIL(J-1)-ZUP*(SDPTH(J)-SDPTH(J-1))/
     +            (SDPTH(I)-SDPTH(1))
        ENDDO
        GOTO 1
       ENDIF
      ENDDO
       
    1 CONTINUE
C   SPLIT THE LOWER ZONE INTO SOIL LAYERS
      IF(NUPL .LT. NSOIL) THEN
       DO I=NUPL+1,NSOIL
cvk 10/2011        xx1=1.5*(SDPTH(I)-SDPTH(1))
cvk 10/2011        IF(ZTOT .LE. xx1 .OR. I .EQ. NSOIL) THEN
cvk 10/2011         IF(ZTOT .GT. ZMAX) THEN
cvk 10/2011          RTLOW=(ZMAX-ZUP)/(ZTOT-ZUP)
cvk 10/2011          ZTOT=ZMAX
cvk 10/2011         ENDIF
         DO J=NUPL+1,I
          ZSOIL(J)=ZSOIL(J-1)-(ZTOT-ZUP)*(SDPTH(J)-SDPTH(J-1))/
     +            (SDPTH(I)-SDPTH(NUPL))
         ENDDO
         NLOWL=I
cvk 10/2011         GOTO 2
cvk 10/2011        ENDIF 
       ENDDO
      ENDIF
    2 CONTINUE

      xx1=0.1*SDPTH(NSOIL)
      IF(NLOWL .LT. NSOIL .AND. -ZSOIL(NLOWL) .LT. xx1) THEN
       NSAC=NLOWL
       ZSOIL(NLOWL+1)=ZSOIL(NSAC)+ZSOIL(NSAC)-ZSOIL(NSAC-1)
       IF(-ZSOIL(NLOWL+1).LT.SDPTH(NSOIL)) ZSOIL(NLOWL+1)=-SDPTH(NSOIL)
       NSOIL=NLOWL+1
      ELSE 
       NSOIL=NLOWL
       NSAC=NSOIL
      ENDIF

cc-      if(ipr .gt. 0) then
cc-      write(*,*) ' nsoil=',nsoil,' nupl=',nupl,' nsac=',nsac
cc-      WRITE(*,*) 'UPPER ZONE MODEL LAYERS:'
cc-      WRITE(*,'(24X,A7,F4.2,A)') ' .00 - ',ZSOIL(1)-ZSOIL(2),'m'
cc-      IF(NUPL .GT. 2) WRITE(*,'(24X,F8.2,A3,F4.2,A)')  
cc-     +  (ZSOIL(1)-ZSOIL(I-1),' - ',ZSOIL(1)-ZSOIL(I),'m',I=3,NUPL)
cc-      WRITE(*,*) 'LOWER ZONE MODEL LAYERS:'
cc-      IF(NSAC .NE. NUPL) WRITE(*,'(24X,F4.2,A3,F4.2,A)') 
cc-     +  (ZSOIL(1)-ZSOIL(I-1),' - ',ZSOIL(1)-ZSOIL(I),'m',I=NUPL+1,NSAC)
cc-      IF(NSOIL .NE. NSAC) WRITE(*,'(A24,F4.2,A3,F4.2,A)') 
cc-     +  ' ADDITIONAL SOIL LAYER: ',
cc-     +  ZSOIL(1)-ZSOIL(NSAC),' - ',ZSOIL(1)-ZSOIL(NSOIL),'m' 
cc-      WRITE(*,*) ' UPPER ZONE RATIO = ',RTUP,
cc-     +    '  LOWER ZONE RATIO = ',RTLOW
cc-      endif	

      if(nupl .eq. nlowl) then
#if defined(CHPS)
       write(MESSAGESTRING,*) 'ERROR: NUPL>5',nupl,smax,nsac,supm,slowm,zup,
     +       zlow,(zsoil(i),i=1,nsoil)
        call logfromfortran( FATAL_LEVEL, MESSAGESTRING )
#else
       write(*,*) 'ERROR: NUPL>5',nupl,smax,nsac,supm,slowm,zup,
     +       zlow,(zsoil(i),i=1,nsoil)
       stop
#endif
      endif
        
      RETURN
      END     
