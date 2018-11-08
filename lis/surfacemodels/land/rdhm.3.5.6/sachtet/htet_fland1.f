C NEW VERSION OF SAC-SMA THAT ACCOUNTS FOR HEAT AND VEGETATION EFFECTS
C SAC-HTET VERSION 1 2010
C MEMBER FLAND1
C
      SUBROUTINE HTET_FLAND1(PXV,EDMND,TA,WE,AESC,SH,DT,SACST,FRZST,
     +             SACPAR,
     +             FRZPAR,NSOIL,NUPL,NSAC,IVERS,SURF,GRND,TET,SMC,
     +             SH2O,SACST_PRV,DTFRZ,IDTFRZ,FRZDUP,FRZDBT,FROST,
     +             TSINT, SWINT, SWHINT, DSINT, NDSINT, DSINTW, NDINTW,
     +             NORMALIZE, YEAR,MONTH,DAY,HOUR,
     +             fdown,q2,q2sat, solar,  nroot, pcinp, smax, sfld,
     +             fxexp, rtdis,
     +             dwsat, dksat, bexp, rdst, sfcprs, sfcspd, sfcref, 
     +             rsmin, rsmax, topt, rgl, hs, xlai, sfsl, 
     +             grn, grn_ind, ztmp, z0, 
     +             QUARTZ, bareadj, czil,  cm, ch, penpt,
     +             CHANLOSS, HRAPX, HRAPY, ERROR )

c  DT is in days here
C  DTFRZ IN SEC., IDTFRZ IS # FRZ_STEPS
C.......................................
C     THIS SUBROUTINE EXECUTES THE 'SAC-HTET ' OPERATION FOR ONE TIME
C         PERIOD.
C.......................................
C     SUBROUTINE INITIALLY WRITTEN AS SAC-SMA BY. . .
C            ERIC ANDERSON - HRL     APRIL 1979     VERSION 1
C     MODIFIED TO ACCOUNT FOR HEAT AND VEGETATION EFFECT (SAC-HTET) BY
C            VICTOR KOREN - HRL      OCTOBER 2010   VERSION 1 
C.......................................

CVK  FROZEN GROUND CHANGES
CVK  UZTWC,UZFWC,LZTWC,LZFSC,LZFPC ARE TOTAL WATER STORAGES
CVK  UZTWH,UZFWH,LZTWH,LZFSH,LZFPH ARE UNFROZEN WATER STORAGES

      PARAMETER (T0 = 273.16)
      REAL SACPAR(*),FRZPAR(*),SACST(*),FRZST(*),SACST_PRV(*)
      real smc(*),sh2o(*)
      integer year,month,day,hour, ERROR
      real smct(10), sh2ot(10), sh2ota(10)
      
c SACPAR() is array of original SAC parameters, and FRZPAR() is array
c of frozen ground parameters and calculated constants
c SACST() and FRZST() same for states
        
c  delited real FGCO(6),ZSOIL(6),TSOIL(8),FGPM(11)      
      REAL LZTWM,LZFSM,LZFPM,LZSK,LZPK,LZTWC,LZFSC,LZFPC
      REAL LZTWH,LZFSH,LZFPH

CVK----------------------------------------------------------------
CVK_02  NEW COMMON STATEMENT FOR DESIRED SOIL LAYERS
CVK     THIS VERSION HAS HARD CODED OUTPUT SOIL LAYERS
CVK     LATER ON IT SHOULD BE CHANGED TO MAKE THEM VARIABLE 
c      INTEGER NINT/5/,NINTW/5/
      INTEGER NDSINT,NDINTW, NORMALIZE
      REAL TSINT(*),SWINT(*),SWHINT(*)
ck      REAL DSINT(10)/0.075,0.15,0.35,0.75,1.5,0.0,0.0,0.,0.,0./
c      REAL DSINT(10)/0.10,0.40,0.6,0.75,1.5,0.0,0.0,0.,0.,0./      
ck      REAL DSINTW(10)/0.075,0.15,0.35,0.75,1.5,0.0,0.0,0.,0.,0./
      REAL DSINT(*), DSINTW(*)
c      REAL DSINTW(10)/0.10,0.40,0.6,0.75,1.5,0.0,0.0,0.,0.,0./
      REAL TSTMP(10),DSMOD(10),SWTMP(10),SWHTMP(10)
c      SAVE DSINT,DSINTW,NINT,NINTW
CVK----------------------------------------------------------------      

cc      DIMENSION EPDIST(24)
C     COMMON BLOCKS
cc      COMMON/FSMPM1/UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,REXP,LZTWM,
cc     1              LZFSM,LZFPM,LZSK,LZPK,PFREE,SIDE,SAVED,PAREA
CVK
CVK      COMMON/FPMFG1/FGPM(10)
c--      COMMON/FPMFG1/itta,FGPM(15),ivers,ifrze      

CVK_02  NEW COMMON BLOCK FOR INTERPOLATED SOIL TEMP AND SOIL PARAMETERS
c--      COMMON/TSLINT/TSINT(10),NINT,SWINT(10),SWHINT(10),NINTW
c--      COMMON/FRDSTFG/SMAX,PSISAT,BRT,SWLT,QUARTZ,STYPE,NUPL,NSAC,
c--     +               RTUZ,RTLZ,DZUP,DZLOW      
cc      COMMON/FRZCNST/ FRST_FACT,CKSOIL,ZBOT
c--      COMMON/FRZCNST/ FRST_FACT,ZBOT      
      
CVK
CVK  NEW FG VERSION PARAMETERS & SOIL LAYER DEFINITION:
CVK          FGPM(1) - SOIL TEXTURE CLASS
CVK          FGPM(2) - OPTIONAL, SOIL TEMPERATURE AT THE 3M DEPTH
CVK          FGPM(3) - OPTIONAL, POROSITY OF RESIDUE LAYER
CVK          PAR(18) [if no calb=FGPM(4)] - RUNOFF REDUCTION PARAMETER 1
CVK          PAR(19) [if no calb=FGPM(5)] - RUNOFF REDUCTION PARAMETER 2
CVK          PAR(20) [if no calb=FGPM(6)] - RUNOFF REDUCTION PARAMETER 3
CVK          FGPM(6) - RUNOFF REDUCTION PARAMETER 3 (FOR ERIC'S VERSION ONLY)
CVK          FGPM(7) - NUMBER OF SOIL LAYERS 
CVK          FGPM(8)-FGPM(15) - DEPTHS OF SOIL LAYERS (M), NEGATIVE.
CVK                             FIRST LAYER (RESIDUE) DEPTH=-0.03M

c--      COMMON/FSMCO1/UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC,FGCO(6),RSUM(7),
c--     1   PPE,PSC,PTA,PWE,PSH,TSOIL(8)     

c--      COMMON/FSUMS1/SROT,SIMPVT,SRODT,SROST,SINTFT,SGWFP,SGWFS,SRECHT,
c--     1              SETT,SE1,SE3,SE4,SE5

C
C    ================================= RCS keyword statements ==========
c--      CHARACTER*68     RCSKW1,RCSKW2
c--      DATA             RCSKW1,RCSKW2 /                                 '
c--     .$Source: /fs/hsmb5/hydro/CVS_root/gorms/sachtet/htet_fland1.f,v $
c--     . $',                                                             '
c--     .$Id: htet_fland1.f,v 1.8 2012-05-02 17:54:30 zcui Exp $
c--     . $' /
C    ===================================================================

cvk  ------ 3/2009  definition of new evapotranspiration variables  ------  
cvk  soil moisture related arrays here defined without 0.03 retentin layer
      real et(5),ets(5),smcs(5),sh2os(5),rtdis(5),zsoil(5)
      real rhstt(5),ai(5),bi(5),ci(5),sice(5),sh2ofg(5),sh2oa(5)
      real grn(12),grnday(1)
      integer nsoils,nroot,optpc,rdst
      integer iter/0/
      save iter

      real pcinp, smax, sfld, fxexp, dwsat, dksat, bexp, sfcprs
      real sfcspd, rsmin, rsmax, topt, rgl, hs, xlai, sfsl
      real solardt(144),tmax(5000),tmin(5000)
      real tair(100000),wspd(100000),pcp(100000)

      real ztmp, z0, cm, ch, pc, rc, rcs, rct, rcq, rcsoil, QUARTZ, t24
      real etp, rch, epsca, rr, flx2, bareadj, dslft, x1, x2, x3
      real HRAPX, HRAPY
      real czil, penpt, CHANLOSS
      
      logical snowng,frzgra
     
cvk next statement for tests only
      real lztb,lzhb,tbal(5),hbal(5)
      
cvk 08/17  to output basin average PED after Penman estimation
      real bsnavg_ped
      integer nbsnped,yearx,monthx,dayx,hourx
      save bsnavg_ped,nbsnped,yearx,monthx,dayx,hourx
      
C--------  start first part of output PED  ----------------
C  uncomment next commented lines if need PED output
C      if(iter .eq. 0) then
C       yearx=year
C       monthx=month
C       dayx=day
C       hourx=hour
C       bsnavg_ped=0.0
C       nbsnped=0
C       write(*,'(a)') 'RDHM OUTPUTS  PED  L    MM   1    '
C       write(*,'(a)') '01  1996 12   2002  1   F13.6'
C      endif
C---------------  end output PED  -------------------------
       ERROR = 0

cfews=========================================================================
cfews  define default values for canopy evaporation
cfews  in this version this effect is not accounted because no snow fall input
cfews  Later on it can be added by using snow17 sfall output
      cmc=0.0
      cmcmax=0.0
      cfactr=0.5
cfews=========================================================================      
      
c04 Define Noah parametric data by calling evap_init
cvk 08/09  this also includes reading daily tmax/tmin
C      if(iter .eq. 0) then
C       grnfile=''
C       grnfile=userfile(1:lnuserfl)
C        call EVAP_INIT(
C     +        
C     +       czil,ztmp,z0,grn,tair,wspd,
Ccfews     +        pcp,itair,iwspd,ipcp,diravg,grnfile,iprnt,bareadj,quartz,
C     +        pcp,itair,iwspd,ipcp,grnfile,iprnt,bareadj,quartz,
Ccfews     +        penpt,ptdst,isnmod,sfcref)
C     +        penpt,isnmod,sfcref)     
C
C       write(50,'(a,5f7.4)') '$Soil moisture:    ',(smc(i),i=2,nsoil)
C       write(50,'(a,5f7.4)') '$Liquid SM:        ',(sh2o(i),i=2,nsoil)
C       write(50,'(a,5f7.2)') '$Soil temperature: ',(frzst(i),i=2,nsoil)
C       write(50,'(a)') '$'
Cc       if(iprnt .eq. 1) close (50)
C
Ccvk 08/09  Estimate initial soil T at desired layers: only in 1st call 
      NMOD=NSOIL+1
      DSMOD(1)=-0.5*FRZPAR(10)
      TSTMP(1)=FRZST(1)
      DSMOD(NMOD)=-FRZPAR(5)
      TSTMP(NMOD)=FRZPAR(2)-T0
      DO i=2,nmod-1
       DSMOD(I)=-0.5*(FRZPAR(I+8)+FRZPAR(I+9))
       TSTMP(I)=FRZST(I)
      ENDDO 
      CALL SOIL_INT1(TSTMP,NMOD,DSMOD,DSINT,NDSINT,TSINT)
      
cvk 08/09  turn off/on Noah plant coefficient calculations 
cfews       if(pc .lt. 0.0) then
cfews        optpc=1
cfews       else
cfews        optpc=0
cfews       endif
       
cfews  make global constant 
cfews       ch = 1.E-4
cfews       cm = 1.E-4
cfews       windadjpc=1.0

cc       wndadjmax=0.0

cvk 10/09 -------------------------------  !!!!!  -------------------------
cvk 10/09  this is temporaly adjustment to SM for sites if swlt differs much 
cvk 10/09  from stxt estimate. It affects only swint output but not smc
       swlt=frzpar(9)
cfews       adjwlt=smcdry-swlt
       smcdry=swlt
cvk 10/09 -------------------------------  !!!!!  -------------------------       
       
       dthr=dt*24.
       dtm=dthr*60.
       dtsec=dtm*60.
       npix=1
       ntmx=0


     
C      endif
C      if( iter .EQ. 0 ) open( 50, file='dbg_output.txt' )
      iter = iter+1
      iprnt = 0

cfews  make cm the same as ch (Noah assumption)
C      cm = ch 
cvk 09/09 initialization of runoff excess from ssrt subroutine (entire column filled)
cvk       possible only in Noah redistribution
      rexces=0.0             

cfews  Excluded for fews version because no TS input
cfewscvk  if no grid, replace hourly temperature from TS
C      if(itair .eq. 0) TA=tair(iter)
cfews
cfewscvk  replace grid hourly precipitation from TS
cfewscvk  however, if TS missing value, use grid value
C      if(ipcp .eq. 0 .and. isnmod .eq. 0) pxv = pcp(iter)
cfews  end exclusion
      
cvk  if wind speed TS available replace default constant value
C      if(iwspd .eq. 1) sfcspd=wspd(iter)
c        write(*,'(i3,6f8.3)') iter,PXV,EDMND,TA,WE,AESC,SH
CVK  ADDITION FOR FROZEN DEPTH ESTIMATION
c--      SAVE IIPREV,FRSTPREV,FRZPREV
c--      IF(FGCO(1) .EQ. 0.) THEN
c--       IIPREV=0
c--       FRSTPREV=0.
c--       FRZPREV=0.
c--      ENDIF
CVK--------------------------------------

c define major parameters from the array
      UZTWM=SACPAR(1)
      UZFWM=SACPAR(2)
      ADIMP=SACPAR(5)
      LZTWM=SACPAR(9)
      LZFSM=SACPAR(10)
      LZFPM=SACPAR(11)
      PAREA=1.0-SACPAR(4)-ADIMP
      IF(IVERS .NE. 0) CKSL=FRZPAR(4)  

c define states from the array
      UZTWC=SACST(1)
      UZFWC=SACST(2)
      LZTWC=SACST(3)
      LZFSC=SACST(4)
      LZFPC=SACST(5)
      ADIMC=SACST(6)

      IF(IVERS .EQ. 0) THEN
CVK  OLD FROZEN GROUND VERSION: KEEP UNFROZEN WATER = TOTAL
C--       RUZICE=UZK
C--       RLZICE=LZSK
C--       RUZPERC=1.0
       UZTWH=UZTWC
       UZFWH=UZFWC
       LZTWH=LZTWC
       LZFSH=LZFSC
       LZFPH=LZFPC
      ELSE        
CVK  NEW FROZEN GROUND VERSION: USE ESTIMATED UNFROZEN WATER
CVK  REDEFINE UNFROZEN WATER VARIABLES IF THEY ARE NEGATIVE
C       
       UZTWH=FRZST(6)
       IF(UZTWH .LT. 0.) UZTWH=0.
       UZFWH=FRZST(7)
       IF(UZFWH .LT. 0.) UZFWH=0.
       LZTWH=FRZST(8)
       IF(LZTWH .LT. 0.) LZTWH=0.
       LZFSH=FRZST(9)
       IF(LZFSH .LT. 0.) LZFSH=0.
       LZFPH=FRZST(10)
       IF(LZFPH .LT. 0.) LZFPH=0.

CVK  RUZICE & RLZICE ARE REDUCTION OF FREE WATER MOVEMENT 
CVK  BASED ON KULIK'S THEORY: Kfrz = Kunfrz/(1+FGPM(4)*ICE)**2 
       ZUP=FRZPAR(9+NUPL)
       ZLW=FRZPAR(9+NSAC)
       RUZICE=0.001*FRZPAR(6)*(UZTWC-UZTWH+UZFWC-UZFWH)/
     +        (FRZPAR(10)-ZUP)
       RLZICE=0.001*FRZPAR(7)*(LZTWC-LZTWH+LZFSC-LZFSH)/
     +        (ZUP-ZLW)
       RUZPERC=1.0
       IF(RUZICE .EQ. 0.) THEN
        RUZICE = SACPAR(3)
       ELSE 
        RUZPERC=1.0/((1.+CKSL*RUZICE)**2)
        RUZICE=1.-EXP(LOG(1.-SACPAR(3))*RUZPERC)
       ENDIF
       IF(RLZICE .EQ. 0.) THEN
        RLZICE = SACPAR(12)
       ELSE  
        RLZICE=1.0/((1.+CKSL*RLZICE)**2)
        RLZICE=1.-EXP(LOG(1.-SACPAR(12))*RLZICE)
       ENDIF 
      ENDIF
c      if(uztwc .ne. uztwh) WRITE(*,*) 'ST1=',uztwc,uztwh
c      if(uzfwc .ne. uzfwh) WRITE(*,*) 'ST2=',uzfwc,uzfwh
c      if(lztwc .ne. lztwh) WRITE(*,*) 'ST3=',lztwc,lztwh
c      if(lzfsc .ne. lzfsh) WRITE(*,*) 'ST4=',lzfsc,lzfsh
c      if(lzfpc .ne. lzfph) WRITE(*,*) 'ST5=',lzfpc,lzfph

c test
c vk 6/2011       if(ipix .eq. 36) then
c vk 6/2011        ipix=iter-((iter-1)/78)*78
c vk 6/2011        write(*,'(2i5,3i3,2f5.0,11f11.4)') ipix,YEAR,MONTH,DAY,
c vk 6/2011     +      HOUR,hrapx,hrapy,frzpar(6),frzpar(7),
c vk 6/2011     +      pxv,ta,(frzst(i),i=1,nsoil),frzpar(2),frzpar(5)
c vk 6/2011        write(*,'(a3,5f13.7)') 'max',uztwm,uzfwm,lztwm,lzfsm,lzfpm
c vk 6/2011        write(*,'(a3,5f13.7)') 'tot',uztwc,uzfwc,lztwc,lzfsc,lzfpc
c vk 6/2011        write(*,'(a3,5f13.7)') 'liq',uztwh,uzfwh,lztwh,lzfsh,lzfph
c vk 6/2011        write(*,'(a3,5f12.7)') 'smc',(smc(i),i=1,nsoil)
c vk 6/2011        write(*,'(a3,5f12.7)') 'sh2',(sh2o(i),i=1,nsoil)
c vk 6/2011        write(*,'(a3,5f12.7)') 'sta',(frzst(i),i=1,nsoil)
c vk 6/2011       endif 
c vk 6/2011       if(iter .gt. 1000) stop

C.......................................
C     COMPUTE EVAPOTRANSPIRATION LOSS FOR THE TIME INTERVAL.
C        EDMND IS THE ET-DEMAND FOR THE TIME INTERVAL
cc      EDMND=EP*EPDIST(KINT)
cVK ADJUST EDMND FOR EFFECT OF SNOW & FOREST COVER.
c from EX1 OFS subroutine
cvk3/10      EDMND=(1.-(1.0-SACPAR(17))*AESC)*EDMND
    
cvk ********************   09/09   ************************************* 
cvk 09/09 added capability to run original SAC-HT evaporation component
c
cfews      if(sac_etdst .eq. -1) then
cfewsC
cfewsC     COMPUTE ET FROM UPPER ZONE.
cfewsCVK      E1=EDMND*(UZTWC/UZTWM)
cfewsCVK  ONLY UNFROZEN WATER CAN BE EVAPORATED
cfews      EDMND=(1.-(1.0-SACPAR(17))*AESC)*EDMND
cfews      E1=EDMND*(UZTWH/UZTWM)
cfews
cfews      RED=EDMND-E1
cfewsC     RED IS RESIDUAL EVAP DEMAND
cfewsCVK      UZTWC=UZTWC-E1
cfews      UZTWH=UZTWH-E1
cfews
cfews      E2=0.0
cfewsCV.K      IF(UZTWC.GE.0.) THEN
cfews      IF(UZTWH.GE.0.) THEN
cfewsCV.K    SUBTRACT ET FROM TOTAL WATER STORAGE
cfews       UZTWC=UZTWC-E1
cfews       GO TO 520
cfews      ENDIF 
cfews
cfewsC     E1 CAN NOT EXCEED UZTWC
cfewsCV.K      E1=E1+UZTWC
cfewsCV.K      UZTWC=0.0
cfews      E1=E1+UZTWH
cfews      UZTWH=0.0
cfewsCV.K   REDUCE TOTAL TENSION WATER BY ACTUAL E1
cfews      UZTWC=UZTWC-E1
cfews      IF(UZTWC .LT. 0.0) UZTWC=0.0
cfews            
cfews      RED=EDMND-E1
cfewsCV.K      IF(UZFWC.GE.RED) GO TO 221
cfews      IF(UZFWH.GE.RED) GO TO 521
cfews
cfewsC     E2 IS EVAP FROM UZFWC.
cfewsCV.K      E2=UZFWC
cfewsCV.K      UZFWC=0.0
cfews      E2=UZFWH
cfews      UZFWH=0.0
cfewsCV.K   REDUCE TOTAL FREE WATER BY ACTUAL E2
cfews      UZFWC=UZFWC-E2
cfews      IF(UZFWC .LT. 0.0) UZFWC=0.0
cfews            
cfews      RED=RED-E2
cfews      GO TO 525
cfews  521 E2=RED
cfewsCVK   SUBTRACT E2 FROM TOTAL & UNFROZEN FREE WATER STORAGES
cfews      UZFWC=UZFWC-E2
cfews      UZFWH=UZFWH-E2
cfews      RED=0.0
cfews  520 IF((UZTWC/UZTWM).GE.(UZFWC/UZFWM)) GO TO 525
cfewsC     UPPER ZONE FREE WATER RATIO EXCEEDS UPPER ZONE
cfewsC     TENSION WATER RATIO, THUS TRANSFER FREE WATER TO TENSION
cfews      UZRAT=(UZTWC+UZFWC)/(UZTWM+UZFWM)
cfews
cfewsCV.K  ACCOUNT FOR RATIO OF UNFROZEN WATER ONLY
cfewsCV.K  AND ADJUST FOUR SOIL STATES 
cfewsCV.K      UZTWC=UZTWM*UZRAT
cfewsCV.K      UZFWC=UZFWM*UZRAT
cfews      DUZTWC=UZTWM*UZRAT-UZTWC
cfews      IF(DUZTWC .GT. UZFWH) DUZTWC=UZFWH 
cfewsCV.K  TRANSFERED WATER CAN NOT EXCEED UNFROZEN FREE WATER
cfews      UZTWC=UZTWC+DUZTWC
cfews      UZTWH=UZTWH+DUZTWC
cfews      UZFWC=UZFWC-DUZTWC
cfews      UZFWH=UZFWH-DUZTWC
cfews            
cfewsCV.K  CHECK UNFROZEN WATER STORAGES TOO
cfews  525 IF (UZTWC.LT.0.00001) THEN
cfews       UZTWC=0.0
cfews       UZTWH=0.0
cfews      ENDIF 
cfews      IF (UZFWC.LT.0.00001) THEN
cfews       UZFWC=0.0
cfews       UZFWH=0.0
cfews      ENDIF 
cfewsC
cfewsC     COMPUTE ET FROM THE LOWER ZONE.
cfewsC     COMPUTE ET FROM LZTWC (E3)
cfewsCV.K      E3=RED*(LZTWC/(UZTWM+LZTWM))
cfewsCV.K      LZTWC=LZTWC-E3
cfewsCV.K      IF(LZTWC.GE.0.0) THEN
cfewsCV.K  ONLY UNFROZEN WATER CAN BE EVAPORATED
cfews      E3=RED*(LZTWH/(UZTWM+LZTWM))
cfews      LZTWH=LZTWH-E3
cfews      IF(LZTWH.GE.0.0) THEN
cfews       LZTWC=LZTWC-E3
cfews       GO TO 526
cfews      ENDIF
cfews       
cfewsC     E3 CAN NOT EXCEED LZTWC
cfewsCV.K      E3=E3+LZTWC
cfewsCV.K      LZTWC=0.0
cfews      E3=E3+LZTWH
cfews      LZTWH=0.0
cfewsCV.K   REDUCE TOTAL TENSION WATER BY E3
cfews       LZTWC=LZTWC-E3
cfews             
cfews  526 RATLZT=LZTWC/LZTWM
cfews      RATLZ=(LZTWC+LZFPC+LZFSC-SACPAR(16))/(LZTWM+LZFPM+LZFSM
cfews     +       -SACPAR(16))
cfews      IF(RATLZT.GE.RATLZ) GO TO 530
cfewsC     RESUPPLY LOWER ZONE TENSION WATER FROM LOWER
cfewsC     ZONE FREE WATER IF MORE WATER AVAILABLE THERE.
cfews      DEL=(RATLZ-RATLZT)*LZTWM
cfewsCV.K  ONLY UNFROZEN WATER CAN BE TRANSFERED
cfews      SFH=LZFSH+LZFPH
cfews      IF(DEL .GT. SFH) DEL=SFH
cfews      LZFSH=LZFSH-DEL
cfews      IF(LZFSH .GE. 0.0) THEN
cfewsC     TRANSFER FROM LZFSC TO LZTWC.      
cfews       LZFSC=LZFSC-DEL
cfews      ELSE
cfewsC     IF TRANSFER EXCEEDS LZFSC THEN REMAINDER COMES FROM LZFPC
cfews       LZFPC=LZFPC+LZFSH
cfews       LZFPH=LZFPH+LZFSH
cfews       xx=LZFSH+DEL
cfews       LZFSC=LZFSC-xx
cfews       LZFSH=0.0
cfews      ENDIF
cfews      LZTWC=LZTWC+DEL
cfews      LZTWH=LZTWH+DEL
cfews
cfewsCV.K      LZTWC=LZTWC+DEL
cfewsCV.K      LZFSC=LZFSC-DEL
cfewsCV.K      IF(LZFSC.GE.0.0) GO TO 230
cfewsCV.K      LZFPC=LZFPC+LZFSC
cfewsCV.K      LZFSC=0.0
cfews
cfewsCV.K  CHECK UNFROZEN WATER STORAGE
cfews  530 IF (LZTWC.LT.0.00001) THEN
cfews       LZTWC=0.0
cfews       LZTWH=0.0 
cfews      ENDIF 
cfewsC       
cfewscvk *********** end original SAC evaporation *******************************  
cfews
cfews      else
C
cvk ----- 3/2009 new evapotranspiration component  -----------------------
cvk array sh2os is the same as sh2o but reduced by first element, 
cvk retention layer, meaning that it starts from the actual first soil layer
       nsoils = nsoil-1
       do i=1,nsoils
        sh2os(i)=sh2o(i+1)
        smcs(i)=smc(i+1)
        sice(i)=smcs(i)-sh2os(i)
        zsoil(i)=frzpar(10+i)-frzpar(10)
       enddo
       
c  calculate daily (same hourly) greenness
       call get_daily_grn(year,month,day,npix,grn,grnday)
       shdfac=grnday(1)

c------------------ canopy resistance -----------------------------------
c 08/09  calculate canopy resistance if PC not constant 
cfews       if(optpc .eq. 1) then
cvk 08/09 Noah subroutines require Kelvin T
       tak=ta+T0
       tsl1=tsint(1)+T0
cvk 08/09 CALCULATE SLOPE OF SAT SPECIFIC HUMIDITY CURVE: DQSDT2
        dqsdt2 = dqsdt_ohd(tak,sfcprs)
cvk 08/09  virtual temperature at 1st soil and atmospheric layers
cvk 08/09  soil T from 1st interpolated layer (usually better to use 0.05m)
       thzv = tsl1*(1.0+0.61*q2)
       thav = (tak+0.0098*ztmp)*(1.0+0.61*q2)

       chref=ch
       cmref=cm
cvk 08/09  calculate surface layer exchange coefficients (cm not in use now)
       call sfcdif(ztmp,z0,thzv,thav,sfcspd,czil,cm,ch)

       if(pcinp .ne. -1.0) then
        pc=pcinp
       else 

cvk 08/09  select daily tmax/tmin and calculate time step solar radiation
cfews        if(mod(iter-1.0,24./dthr) .eq. 0.) then
cvk 08/09 ntmx controlls day loops
cfews         ntmx=ntmx+1
cfews         call soltmxmn(month,day,dtsec,tnoon,dtlocal,tmax(ntmx),
cfews     +                tmin(ntmx),rlat,rsmin,rsmax,rgl,xlai,sfsl,solardt)
cfews         nstep=0
cfews        endif

cvk 08/09 nstep controlls time loops per day
cfews        nstep=nstep+1
cfews        solar=solardt(nstep)

cvk 08/09  calculate actual q2 and saturated q2sat specific humidity, kg/kg 
cfews        call vapoprs(tak,sfcprs,q2sat,q2)
cfews   get q2 and q2sat from from 'do_sac'


cvk 08/09  calculate stress factors
        call canres(solar,ch,tak,q2,sfcprs,smcs,zsoil,nsoils,smcdry,
     &               sfld,rsmin,rc,pc,nroot,q2sat,dqsdt2,topt,rsmax, 
     &               rgl,hs,xlai,rcs,rct,rcq,rcsoil)

       endif
cvk  ----------------  end canopy resistance calculation   -----------

cvk  -----------------------------------------------------------------
cvk 4/2010  *********  test penman option  **********
cvk  t2v is a virtual air temperature
cvk  th2 is a potential air temperature at measurement hight
cvk  rlwdn is downwart longwave radiation from Noah 
       t2v=tak*(1.0+0.61*q2)
       th2=tak+(0.0098*ztmp)

cfews  move estimation of rlwdn and albedo to do_sac
cfews       emiss=1.0-0.261*EXP((-7.77E-4)*(273-tak)**2)
cfews       rlwdn=emiss*5.672E-8*tak**4
cc       albedo=soilalb+aesc*(snoalb-soilalb)
cfews       albedo=0.15+aesc*(0.7-0.15)
cfews       fdown=solar*(1.0-albedo)+rlwdn

cvk Calculate soil heat flux
       df1=CND_JOHNS(12.0,sh2os(1),SICE(1),SMAX,QUARTZ)
C ----------------------------------------------------------------------
C NEXT ADD SUBSURFACE HEAT FLUX REDUCTION EFFECT FROM THE 
C OVERLYING GREEN CANOPY, ADAPTED FROM SECTION 2.1.2 OF 
C PETERS-LIDARD ET AL. (1997, JGR, VOL 102(D4))
C ----------------------------------------------------------------------
       DF1 = DF1 * EXP(-2.0*SHDFAC)
       DSOIL = -0.5*(ZSOIL(1)+frzpar(10))
       SSOIL = DF1*(frzst(1)-frzst(2))/dsoil

cvk pxv is only liquid; so no frozen water
       snowng=.FALSE.
       frzgra=.FALSE.
              
cvk convert precip into kg/(m**2 sec)
       pxvs=pxv/dtsec

       call penman_ohd(tak,sfcprs,ch,t2v,th2,pxvs,fdown,t24,ssoil,q2,
     +              q2sat,etp,rch,epsca,rr,snowng,frzgra,dqsdt2,flx2) 
       etp=etp*dtsec
c         write(*,'(a,i5,5f10.5)') 'ch',iter,ch,sfcspd,etp,pc,rc
       if(etp .le. 0.0) then
c------------------------------Zero case evaporation  ---------------
cvk  If Penman generates negative etp (potential evaporation),
cvk  no direct or evapotranspiration occures; instead dew water
cvk  will be assumed and adds up to precipitation
        pxv=pxv-etp
        edir=0.0
        ec=0.0
        ett1=0.0
        do i=1,nsoils
         et(i)=0.0
        enddo
        etat=etp
        e1=0.0
        e2=0.0
        red=0.0
        edmnd=0.0
        edmndbare=0.0
       else       

c-----------------  Non-Zero case Evaporation  ------------------------
cvk  calculate all components of evapotranspiration
cvk all evaporation components are in mm/dt
cvk---------------------------------------------------------------
        
cvk  selection of ETP source 
cfews       if(penpt .gt. 0.0) then
cfews        edmnd=etp*penpt
cfews        edmndbare=edmnd
cfews       else
cfews        edmndbare=edmnd
cfews        EDMND=EDMND*rcs/(sfsl*dt)
cfews        if(ptdst .ne. 0.0) then
cfewscvk use non-uniform PET for direct evaporation
cfews         edmndbare=edmnd
cfews        endif

cvk 8/2011 addition
       if(penpt .ne. 0.0) then
c use Penman-based PET: if penpt = 999.0 or -999.0, estimate adjustment factors,
c                       if other values, use them as adjustment factors        
        if(abs(penpt) .eq. 999.0) then
         call penmanadj(grn_ind,year,month,day,penadj)
         edmnd=etp*penadj
        else 
         edmnd=etp*abs(penpt)
        endif 
        edmndbare=edmnd
       else 
        EDMND=EDMND*rcs/(sfsl*dt)
        edmndbare=edmnd

cvk  5/10 wind effect correction to PET
        ediff=abs(sfcref-sfcspd)
        if(sfcref .gt. 0.0 .and. ediff .gt. 0.01) then
         call sfcdif(ztmp,z0,thzv,thav,sfcref,czil,cmref,chref)
         call penman_ohd(tak,sfcprs,chref,t2v,th2,pxvs,fdown,t24,ssoil,q2,
     +              q2sat,etpref,rch,epsca,rr,snowng,frzgra,dqsdt2,flx2)
         etpref=etpref*dtsec
         if(etpref .le. 0.0) then
c vk 2011  changed 'windadjp' name to be more understandable
          windadj=windadjpc
         else 
          windadj=etp/etpref
          if(windadj .gt. 50.0) windadj=windadjpc 
         endif                  
c       if(windadj .gt. wndadjmax) wndadjmax=windadj
c       write(*,'(a,i7,3f6.3,3f6.2)') 'pmn',iter,etp,etpref,edmnd,
c     +                                 windadj,windadjpc  
         EDMND = EDMND*windadj
         windadjpc=windadj
        endif
         
       endif

cvk  08/17 ------  start  second part of output PED  ----------       
cvk  08/17  calculate basin average ped after penman estimation
c       if(yearx .ne. year .or. monthx .ne. month .or. 
c     +    dayx .ne. day .or. hourx .ne. hour) then
c        bsnavg_ped=bsnavg_ped/nbsnped
c        if(yearx .gt. 1999) then
c         ixy=yearx-2000
c        else
c         ixy=yearx-1900
c        endif  
c        write(*,4231) dayx,monthx,ixy,hourx+1,edmnd
c4231  format(10x,3i2,i4,F13.6)
c        yearx=year
c        monthx=month
c        dayx=day
c        hourx=hour
c        nbsnped=0
c        bsnavg_ped=0.0
c       endif
c       bsnavg_ped=bsnavg_ped+edmnd
c       nbsnped=nbsnped+1
c      write(*,'(i4,3f7.4,2f7.2)') nbsnped,edmnd,etp,etpref,windadj,ediff
cvk  08/17 ------------   end   --------------------------------                       

cvk 4/2010  bare soil only PET adjustment
cvk  PET will be reduced if greenness < bareadj
cvk assumed a priori relationship: PEadj=1.144*grn+0.202
cvk that leads bare soil effect factor 0.8
cvk

cfews   excluded bareadj option: no adjustment, just switch Ek-Chen
cfews       if(bareadj .gt. 0.0) then
cfews        if(shdfac .le. bareadj) then
cfews         edmndbare=edmndbare*(1.0-0.8*(1.0-shdfac))
cfews         if(edmndbare .lt. 0.0) edmndbare=0.0
cfews        endif
cfews       endif

cvk 4/2010  estimate upper zone weighted SM content to be used for direct evap.
cfews        if(diravg .eq. 0.0) then
         avgsmup=sh2os(1)
cfews        else 
cfews        avgsmup=0.0
cfews        do i=1,nupl-1
cfews         if(i .eq. 1) then
cfews          avgsmup=avgsmup-sh2os(i)*zsoil(i)
cfews         else 
cfews          avgsmup=avgsmup-sh2os(i)*(zsoil(i)-zsoil(i-1))
cfews         endif 
cfews        enddo
cfews        avgsmup=-avgsmup/zsoil(nupl-1)
cfews        endif

cvk  calculate evapotranspiration
cvk 6/10  make edmndbare = 0 if soil surface T<0
        if(frzst(1) .le. 0.0) then
c      write(*,'(a,i7,2f6.3,3f6.1,)') 'edmndbare=0:',
c     +      iter,edmndbare,smc(2)-sh2o(2),ta,frzst(1),frzst(2)
         edmndbare = 0.0
        endif 
         
        call EVAPO(etat,nsoils,cmc,EDMND,edmndbare,sh2os,smax,pc,
     +           swlt,sfld,shdfac,cmcmax,smcdry,cfactr,edir,ec,et,
     +           ett1,nroot,rtdis,fxexp,avgsmup,uztwh,uztwm,bareadj)

c        if(frzst(1) .le. 0.0) write(*,*) iter,edir
c         write(*,'(a,i6,f6.3,7f7.4)') 'evap',iter,pc,
c     +                  (et(ij),ij=1,nsoils),edir,ett1
c 3/10 adjustment to bare soil evaporation to be only from not snow covered area
c 3/10 however, evapotranspiration occures even from snow covered area
c 3/10 Noah model evaporate from snow cover but do not from vegetation if snow
        edir = (1-AESC)*edir 

cvk canopy water retention state definition
c actualy it's not incorporated yet because SAC does not know split
c between rain and melt. Later on need to be incorporated
        rhsct=shdfac*pxv-ec
        drip=0.
        excess=cmc+rhsct
        if(excess .gt. cmcmax) drip=excess-cmcmax
        px=pxv
        pxv=(1-shdfac)*pxv+drip 
       endif
cvk 4/2010  ---- end dew water case  ----------------             

c Cumulate evaporation for SAC zones
       ET1=0.0 
       ET2=0.0
       do i=1,nsoils
        sh2o(i+1)=sh2os(i)
        if(i .le. nupl-1) then
         ET1=ET1+ET(i)
        else
         ET2=ET2+ET(i)
        endif  
       enddo
       ET1=ET1+edir
       E3=ET2

cfewscvk -------------  Selection of SM redistribution  ----------------
cfewscvk       if(sac_etdst .eq. 1) then
cfewscvkcvk ---------  similar to original SAC option  --------------------------
cfewscvkc
cfewsc adjust canopy water storage
cfews        CMC = CMC + RHSCT
cfews        IF (CMC .LT. 1.E-20) CMC=0.0
cfews        CMC = MIN(CMC,CMCMAX)
cfews      
cfewscvk  3/2009 compute upper zone SM change
cfewscvk  3/2009 E1 and E2 evaporation from tension and free zone respectively 
cfews        E2=0.0
cfews        E1=ET1
cfews        RED=EDMND-E1
cfews        UZTWH=UZTWH-E1
cfews        IF(UZTWH.GE.0.) THEN
cfews         UZTWC=UZTWC-E1
cfews        ELSE
cfews         E1=E1+UZTWH
cfews         UZTWH=0.0
cfews         UZTWC=UZTWC-E1
cfewscc        IF(UZTWC .LT. 0.0) UZTWC=0.0
cfews         E2=ET1-E1
cfews         UZFWH=UZFWH-E2
cfews         IF(UZFWH .LT. 0.0)THEN
cfews          E2=UZFWH
cfews          UZFWH=0.0
cfews          RED=RED-E2
cfews         ELSE
cfews          E2=RED
cfews          UZFWH=UZFWH-E2
cfews          RED=0.0
cfews         ENDIF
cfews         UZFWC=UZFWC-E2 
cfews        ENDIF
cfews      
cfewscvk  3/2009  from original sac-ht
cfews  220    IF((UZTWC/UZTWM).GE.(UZFWC/UZFWM)) GO TO 225
cfewsC     UPPER ZONE FREE WATER RATIO EXCEEDS UPPER ZONE
cfewsC     TENSION WATER RATIO, THUS TRANSFER FREE WATER TO TENSION
cfews        UZRAT=(UZTWC+UZFWC)/(UZTWM+UZFWM)
cfews
cfewsCV.K  ACCOUNT FOR RATIO OF UNFROZEN WATER ONLY
cfewsCV.K  AND ADJUST FOUR SOIL STATES 
cfewsCV.K      UZTWC=UZTWM*UZRAT
cfewsCV.K      UZFWC=UZFWM*UZRAT
cfews        DUZTWC=UZTWM*UZRAT-UZTWC
cfews        IF(DUZTWC .GT. UZFWH) DUZTWC=UZFWH 
cfewsCV.K  TRANSFERED WATER CAN NOT EXCEED UNFROZEN FREE WATER
cfews        UZTWC=UZTWC+DUZTWC
cfews        UZTWH=UZTWH+DUZTWC
cfews        UZFWC=UZFWC-DUZTWC
cfews        UZFWH=UZFWH-DUZTWC
cfews
cfewsCV.K  CHECK UNFROZEN WATER STORAGES TOO
cfews  225   IF (UZTWC.LT.0.00001) THEN
cfews         UZTWC=0.0
cfews         UZTWH=0.0
cfews        ENDIF 
cfews        IF (UZFWC.LT.0.00001) THEN
cfews         UZFWC=0.0
cfews         UZFWH=0.0
cfews        ENDIF 
cfews
cfewscvk  3/2009     COMPUTE SM change in the lower zone
cfewsC     COMPUTE ET FROM tension water (E3)
cfewscvk  3/2009  original sac version evaporates only from tension
cfewscvk  3/2009  water. However, after subtraction evaporation from
cfewscvk  3/2009  LZTWH it checks ratio of tension-free water, and
cfewscvk  3/2009  may transfer part of free water to tension depending
cfewscvk  3/2009  on this ratio.
cfewscvk  3/2009  I will keep this assumption meaning that actual
cfewscvk  3/2009  evaporation may be less than defined by TRANS.
cfewscvk  3/2009  It is possible to generate TRANS evaporation based 
cfewscvk  3/2009  on tension water only. See if it'll be needed. 
cfews        E3=ET2
cfews        LZTWH=LZTWH-E3
cfews        IF(LZTWH.GE.0.0) THEN
cfews         LZTWC=LZTWC-E3
cfews        ELSE
cfews         E3=E3+LZTWH
cfews         LZTWH=0.0
cfews         LZTWC=LZTWC-E3
cfews        ENDIF
cfews
cfewscvk  3/2009  old version of redistribution water from free zone
cfews  226   RATLZT=LZTWC/LZTWM
cfews        RATLZ=(LZTWC+LZFPC+LZFSC-SACPAR(16))/(LZTWM+LZFPM+LZFSM
cfews     +       -SACPAR(16))
cfews        IF(RATLZT.GE.RATLZ) GO TO 230
cfewsC     RESUPPLY LOWER ZONE TENSION WATER FROM LOWER
cfewsC     ZONE FREE WATER IF MORE WATER AVAILABLE THERE.
cfews        DEL=(RATLZ-RATLZT)*LZTWM
cfewsCV.K  ONLY UNFROZEN WATER CAN BE TRANSFERED
cfews        SFH=LZFSH+LZFPH
cfews        IF(DEL .GT. SFH) DEL=SFH
cfews        LZFSH=LZFSH-DEL
cfews        IF(LZFSH .GE. 0.0) THEN
cfewsC     TRANSFER FROM LZFSC TO LZTWC.      
cfews         LZFSC=LZFSC-DEL
cfews        ELSE
cfewsC     IF TRANSFER EXCEEDS LZFSC THEN REMAINDER COMES FROM LZFPC
cfews         LZFPC=LZFPC+LZFSH
cfews         LZFPH=LZFPH+LZFSH
cfews         XX=LZFSH+DEL
cfews         LZFSC=LZFSC-XX
cfews         LZFSH=0.0
cfews        ENDIF
cfews        LZTWC=LZTWC+DEL
cfews        LZTWH=LZTWH+DEL
cfewsCV.K  CHECK UNFROZEN WATER STORAGE
cfews  230   IF (LZTWC.LT.0.00001) THEN
cfews         LZTWC=0.0
cfews         LZTWH=0.0 
cfews        ENDIF 
cfews
cfewscvk  3/2009 Actual evaporation from lower zone layers may be less
cfewscvk  3/2009 than estimated by TRANS because of evaporate only from tension water
cfewsC  adjust evaporation in case it was reduced because exceeded storages
cfews        if(ET1 .gt. 10E-10) then
cfews         aup=(E1+E2)/ET1
cfews        else
cfews         aup=1.0
cfews        endif  
cfews        if(ET2 .gt. 10E-10) then
cfews         alo=E3/ET2
cfews        else
cfews         alo=1.0
cfews        endif  
cfews
cfewscvk  subtract evaporation from each soil layer
cfews        do i=2,nsoil
cfews         ee=et(i-1)
cfews         if(i .le. nupl-1) then
cfews          if(i .eq. 2) ee=ee+edir
cfews          sh2o(i)=sh2o(i)-0.001*aup*ee/(frzpar(9+i-1)-frzpar(9+i))
cfews          smc(i)=smc(i)-0.001*aup*ee/(frzpar(9+i-1)-frzpar(9+i))                   
cfews         else
cfews          sh2o(i)=sh2o(i)-0.001*alo*ee/(frzpar(9+i-1)-frzpar(9+i))
cfews          smc(i)=smc(i)-0.001*alo*ee/(frzpar(9+i-1)-frzpar(9+i))
cfews         endif  
cfews        enddo 
cfewscvk 3/2009 -------------- end SAC redistribution  -------------------
      
cfewscvk       else

cvk---------  Noah SM redistribution option  -------------------------

cvk convert evaporation from mm/dt into m/sec
c        dtsec=dt*24.*3600.
        edirs=edir*0.001/dtsec
        do i=1,nsoils
         ets(i)=et(i)*0.001/dtsec
        enddo
C
C IF THE INFILTRATING EVAPORATION RATE IS NONTRIVIAL,
C   (WE CONSIDER NONTRIVIAL TO BE A EVAP TOTAL OVER THE TIME STEP 
C    EXCEEDING ONE ONE-THOUSANDTH OF THE WATER HOLDING CAPACITY OF 
C    THE FIRST SOIL LAYER)
C THEN CALL THE SRT/SSTEP SUBROUTINE PAIR TWICE IN THE MANNER OF 
C   TIME SCHEME "F" (IMPLICIT STATE, AVERAGED COEFFICIENT)
C   OF SECTION 2 OF KALNAY AND KANAMITSU (1988, MWR, VOL 116, 
C   PAGES 1945-1958)TO MINIMIZE 2-DELTA-T OSCILLATIONS IN THE 
C   SOIL MOISTURE VALUE OF THE TOP SOIL LAYER THAT CAN ARISE BECAUSE
C   OF THE EXTREME NONLINEAR DEPENDENCE OF THE SOIL HYDRAULIC 
C   DIFFUSIVITY COEFFICIENT AND THE HYDRAULIC CONDUCTIVITY ON THE
C   SOIL MOISTURE STATE
C OTHERWISE CALL THE SRT/SSTEP SUBROUTINE PAIR ONCE IN THE MANNER OF
C   TIME SCHEME "D" (IMPLICIT STATE, EXPLICIT COEFFICIENT) 
C   OF SECTION 2 OF KALNAY AND KANAMITSU
C PCPDRP IS UNITS OF KG/M**2/S OR MM/S, ZSOIL IS NEGATIVE DEPTH IN M 
C ----------------------------------------------------------------------

        nups=nupl-1
        etles=0.0

cvk 7/2010  new code: change Noah water redistribution  -----------------
c  Indroduced new arrays: smct(nsoils), sh2ot(nsoils), sh2ota(nsoils)
c
        frtup=0.001*uzfwc/(-zsoil(nups))
        frtlo=0.001*(lzfsc+lzfpc)/(zsoil(nups)-zsoil(nsoils))
cvk        frtup=0.001*uzfwc*FRZPAR(6)/(-zsoil(nups))
cvk        frtlo=0.001*(lzfsc+lzfpc)*FRZPAR(7)/(zsoil(nups)-zsoil(nsoils))
        do i=1,nsoils
c vk 6/2011 error         if(i .lt. nups) then
         if(i .le. nups) then
          smct(i)=smcs(i)-frtup
          sh2ot(i)=sh2os(i)-frtup
          if (smct(i) < 0.0) then
            smct(i) = 0.0
            !write(*,*) "SAC-HTET error: smct(", i , ") < 0.0"
            !write(*,*) "fixed it with 0.0", hrapx, hrapy
          endif

          if (sh2ot(i) < 0.0) then
            sh2ot(i) = 0.0
            !write(*,*) "SAC-HTET error: sh2ot(",i,") < 0.0" 
            !write(*,*) "fixed it with 0.0", hrapx, hrapy
          endif
         else
          smct(i)=smcs(i)-frtlo
          sh2ot(i)=sh2os(i)-frtlo
          
          if (smct(i) < 0.0) then
            smct(i) = 0.0
            !write(*,*) "SAC-HTET error: smct(",i,") < 0.0"
            !write(*,*) "fixed it with 0.0", hrapx, hrapy
          endif

          if (sh2ot(i) < 0.0) then
            sh2ot(i) = 0.0
            !write(*,*) "SAC-HTET error: sh2ot(",i,") < 0.0"
            !write(*,*) "fixed it with 0.0", hrapx, hrapy
          endif
         endif
        enddo           

c Calculate total moisture changes
cvk  7/2010  replaced old version   ------------------------
cvk  smcs - tens & free total water cont
cvk  sh2os - tens & free liquid water cont
cvk  smct - tens only total water cont
cvk  sh2ot - tens only liquid water cont
c
c vk 6/2011        if((et(1)+edir) .gt. (-zsoil(1))*smax) then
c vk 6/2011  Two iterrations will be performed if upper layer 
c vk 6/2011  evaporation exceeds 2% of maximum capacity of the layer
c
        eps=(et(1)+edir)/(-zsoil(1)*smax*10.0)
        if(eps .gt. 2.0) then
         call SRT(RHSTT,EDIRS,ETS,SH2OT,SH2OT,NSOILS,ZSOIL,
     &           DWSAT,DKSAT,SMAX,BEXP,SICE,AI,BI,CI,smct,rdst) 
         dummy=0.
         call SSTEP(SH2OFG,SH2OS,DUMMY,RHSTT,RHSCT,DTSEC,NSOILS,SMAX,
     &                CMCMAX,ZSOIL,SMCS,SICE,SMCDRY,AI,BI,CI,rexces,
     &                                                   nups,etles)

         DO K = 1,NSOIL
          SH2OA(K) = (SH2OS(K) + SH2OFG(K)) * 0.5
         END DO

        do i=1,nsoils
c vk 6/2011 error         if(i .lt. nups) then
         if(i .le. nups) then
          smct(i)=smcs(i)-frtup
          sh2ot(i)=sh2os(i)-frtup
          sh2ota(i)=sh2oa(i)-frtup
          if(sh2ot(i) .lt. 0.0) sh2ot(i) = 0.0
          if(sh2ota(i) .lt. 0.0) sh2ota(i) = 0.0
         else
          smct(i)=smcs(i)-frtlo
          sh2ot(i)=sh2os(i)-frtlo
          sh2ota(i)=sh2oa(i)-frtlo
          if(sh2ot(i) .lt. 0.0) sh2ot(i) = 0.0
          if(sh2ota(i) .lt. 0.0) sh2ota(i) = 0.0
         endif
        enddo    
         etles=0.0

         call SRT(RHSTT,EDIRS,ETS,SH2OT,SH2OtA,NSOILS,ZSOIL,
     &           DWSAT,DKSAT,SMAX,BEXP,SICE,AI,BI,CI,smct,rdst) 
     
         call SSTEP(SH2OS,SH2OS,CMC,RHSTT,RHSCT,DTSEC,NSOILS,SMAX,
     &             CMCMAX,ZSOIL,SMCS,SICE,SMCDRY,AI,BI,CI,rexces,
     &                                                nups,etles)

        else
         call SRT(RHSTT,EDIRS,ETS,SH2OT,SH2OT,NSOILS,ZSOIL,
     &           DWSAT,DKSAT,SMAX,BEXP,SICE,AI,BI,CI,smct,rdst) 

         call SSTEP(SH2OS,SH2OS,CMC,RHSTT,RHSCT,DTSEC,NSOILS,SMAX,
     &             CMCMAX,ZSOIL,SMCS,SICE,SMCDRY,AI,BI,CI,rexces,
     &                                                nups,etles)
        endif
cvk 4/2010 decrease et1 if evap below smcdry
        et1=et1+etles
        sdx=0.0
        do i=1,nsoils
         if(i .eq. 1) then
          xu=1000*(sh2os(i)-sh2o(i+1))*(-zsoil(1))
         else
          xu=1000*(sh2os(i)-sh2o(i+1))*(zsoil(i-1)-zsoil(i))
         endif
         sdx=sdx+xu
        enddo   

        dsup=0.
        dslo=0.
        do i=1,nsoils
         snew=sh2os(i)-smcdry
         if(snew .lt. 0.0) snew=0.0
         sold=sh2o(i+1)-smcdry
         if(i .eq. 1) then
          dz=-zsoil(1)*1000.
         else
          dz=(zsoil(i-1)-zsoil(i))*1000.
         endif 
         if(i .le. nupl-1) then
          dsup=dsup+(snew-sold)*dz
         else
          dslo=dslo+(snew-sold)*dz
         endif
         smc(i+1)=smcs(i)
         sh2o(i+1)=sh2os(i)
        enddo 

c Redistributed changes to the upper and lower zone storages
        xx=-1.
        call smc2sac(dsup,dslft,et1,uztwc,uzfwc,xx,uztwh,uzfwh,
     +                      xx,edmnd,e1,e2,red,uztwm,uzfwm,xx)

c test 6/08/2011
c         if(dslft .ne. 0.0) then
c          write(*,*) 'dslfup not 0:',YEAR,MONTH,DAY,HOUR,hrapx,hrapy
c          write(*,*) '  1 ',dsup,dslo,dslft,et1,uztwc,uzfwc,uztwh,
c     +               uzfwh,edmnd,e1,e2,red,uztwm,uzfwm,smax,smcdry
c          write(*,*) '  smcs',(smcs(i),i=1,nsoils)
c          write(*,*) ' sh2os',(sh2os(i),i=1,nsoils)
c         endif
         
        dslo=dslo+dslft
        call smc2sac(dslo,dslft,et2,lztwc,lzfsc,lzfpc,lztwh,lzfsh,
     +                    lzfph,edmnd,x1,x2,x3,lztwm,lzfsm,lzfpm)

cvk----------  end Noah redistribution option ---------------------

        SACST_PRV(1)=UZTWC
        SACST_PRV(2)=UZFWC
        SACST_PRV(3)=LZTWC
        SACST_PRV(4)=LZFSC
        SACST_PRV(5)=LZFPC
      
        frzst(6)=uztwh
        frzst(7)=uzfwh
        frzst(8)=lztwh
        frzst(9)=lzfsh
        frzst(10)=lzfph

cfews      endif
c  end evaporation component
c***************************************************************************      
   
cvk  3/2009  No chnages to original version
C     COMPUTE ET FROM ADIMP AREA.-E5
      E5=E1+(RED+E2)*((ADIMC-E1-UZTWC)/(UZTWM+LZTWM))

C      ADJUST ADIMC,ADDITIONAL IMPERVIOUS AREA STORAGE, FOR EVAPORATION.
      ADIMC=ADIMC-E5
      IF(ADIMC.GE.0.0) GO TO 231
C     E5 CAN NOT EXCEED ADIMC.
      E5=E5+ADIMC
      ADIMC=0.0
  231 E5=E5*ADIMP

C     E5 IS ET FROM THE AREA ADIMP.
C.......................................
C     COMPUTE PERCOLATION AND RUNOFF AMOUNTS.
      TWX=PXV+UZTWC-UZTWM       
C     TWX IS THE TIME INTERVAL AVAILABLE MOISTURE IN EXCESS
C     OF UZTW REQUIREMENTS.

cvk 2013  Channel loss addition changed directly to TWX
cvk 2013  In older version it added to PERC that can lead to overflow lower zone
      IF(CHANLOSS .LT. 0.0) CHANLOSS = 0.0
      closperdt=CHANLOSS*dtsec
      IF(TWX.GE.0.0) GO TO 232
C     ALL MOISTURE HELD IN UZTW--NO EXCESS.
      UZTWC=UZTWC+PXV
CV.K  ADJUST UNFROZEN TENSION WATER
      UZTWH=UZTWH+PXV      

cvk 2013      TWX=0.0
cvk 2013 add channel losses
      twx=closperdt
      GO TO 233
C      MOISTURE AVAILABLE IN EXCESS OF UZTWC STORAGE.
CV.K  232 UZTWC=UZTWM
  232 UZTWH=UZTWH+(UZTWM-UZTWC)
      UZTWC=UZTWM

cvk 2013 add channel losses
      twx=twx+closperdt
cvk 2013 add channel losses  233 ADIMC=ADIMC+PXV-TWX
cvk 2013 subtract added channel losses from additional impervious area
  233 ADIMC=ADIMC+PXV-(TWX-closperdt)
      
C
C     COMPUTE IMPERVIOUS AREA RUNOFF.
      ROIMP=PXV*SACPAR(4)
C      ROIMP IS RUNOFF FROM THE MINIMUM IMPERVIOUS AREA.
      SIMPVT=SIMPVT+ROIMP
C
C     INITIALIZE TIME INTERVAL SUMS.
      SBF=0.0
      SSUR=0.0
      SIF=0.0
      SPERC=0.0
      SDRO=0.0
      SPBF=0.0

cvk 09/09  addition to direct runoff if there was runoff excees in ssrt subr:
cvk        for Noah moisture redistribution option only needed
      SDRO=rexces
C
C     DETERMINE COMPUTATIONAL TIME INCREMENTS FOR THE BASIC TIME
C     INTERVAL
CV.K      NINC=1.0+0.2*(UZFWC+TWX)
CV.K  PERCOLATE UNFROZEN WATER ONLY
      NINC=1.0+0.2*(UZFWH+TWX)
C     NINC=NUMBER OF TIME INCREMENTS THAT THE TIME INTERVAL
C     IS DIVIDED INTO FOR FURTHER
C     SOIL-MOISTURE ACCOUNTING.  NO ONE INCREMENT
C     WILL EXCEED 5.0 MILLIMETERS OF UZFWC+PAV
      DINC=(1.0/NINC)*DT
C     DINC=LENGTH OF EACH INCREMENT IN DAYS.
      PINC=TWX/NINC
C     PINC=AMOUNT OF AVAILABLE MOISTURE FOR EACH INCREMENT.
C      COMPUTE FREE WATER DEPLETION FRACTIONS FOR
C     THE TIME INCREMENT BEING USED-BASIC DEPLETIONS
C      ARE FOR ONE DAY
CVK INTRODUCED REDUCTION (RUZICE & RLZICE) DUE FROZEN GROUND
CVK HOWEVER, PRIMARY RUNOFF IS UNCHANGED
CVK      DUZ=1.0-((1.0-UZK)**DINC)
CVK      DLZS=1.0-((1.0-LZSK)**DINC)
CVK  Linear transformation for frozen ground
cc      DUZ=1.0-((1.0-UZK*RUZICE)**DINC)
cc      DLZS=1.0-((1.0-LZSK*RLZICE)**DINC)
CVK  Non-linear (correct) transformation for frozen ground
      IF(IVERS .EQ. 0) THEN
       DUZ =1.0-((1.0-SACPAR(3))**DINC)
       DLZS=1.0-((1.0-SACPAR(12))**DINC)
      ELSE        
       DUZ=1.0-((1.0-RUZICE)**DINC)
       DLZS=1.0-((1.0-RLZICE)**DINC)
      ENDIF 
      DLZP=1.0-((1.0-SACPAR(13))**DINC)
C
C CVK  ADJUSTMENT TO DEPLETIONS DUE TO FROZEN WATER
         
C.......................................
C     START INCREMENTAL DO LOOP FOR THE TIME INTERVAL.
      DO 240 I=1,NINC
      ADSUR=0.0
C     COMPUTE DIRECT RUNOFF (FROM ADIMP AREA).
      RATIO=(ADIMC-UZTWC)/LZTWM
      IF (RATIO.LT.0.0) RATIO=0.0
      ADDRO=PINC*(RATIO**2)
C     ADDRO IS THE AMOUNT OF DIRECT RUNOFF FROM THE AREA ADIMP.
C
C     COMPUTE BASEFLOW AND KEEP TRACK OF TIME INTERVAL SUM.
CV.K      BF=LZFPC*DLZP
CV.K      LZFPC=LZFPC-BF
CV.K      IF (LZFPC.GT.0.0001) GO TO 234
CV.K      BF=BF+LZFPC
CV.K      LZFPC=0.0
CV.K  BASEFLOW FROM UNFROZEN WATER ONLY   
      BF=LZFPH*DLZP
      LZFPH=LZFPH-BF
      IF (LZFPH.GT.0.0001) THEN
       LZFPC=LZFPC-BF
       GO TO 234
      ENDIF
      BF=BF+LZFPH
      LZFPH=0.0
      LZFPC=LZFPC-BF
      IF(LZFPC .LE. 0.0001) LZFPC=0.0
CV.K-------------------------------------
C      
  234 SBF=SBF+BF
      SPBF=SPBF+BF
CV.K  SUPPLAMENTAL FLOW FROM UNFROZEN WATER ONLY (NOTE, DLZS
CV.K  NOTE, DLZS IS REDUCED DUE FROZEN GROUND
CV.K      BF=LZFSC*DLZS
CV.K      LZFSC=LZFSC-BF
CV.K      IF(LZFSC.GT.0.0001) GO TO 235
CV.K      BF=BF+LZFSC
CV.K      LZFSC=0.0
      BF=LZFSH*DLZS
      LZFSH=LZFSH-BF
      IF(LZFSH.GT.0.0001) THEN
cc?      IF(LZFSH.GT.0.0) THEN      
       LZFSC=LZFSC-BF
c         if(abs(lzfsc-lzfsh) .gt. 0.000001) then
c         if(abs(lzfsc-lzfsh) .gt. 0.000001) then
c         endif 
       GO TO 235
      ENDIF
      BF=BF+LZFSH
      LZFSH=0.0
      LZFSC=LZFSC-BF
      IF(LZFSC .LE. 0.0001) LZFSC=0.0   
CV.K--------------------------------------------
C       
  235 SBF=SBF+BF
C
C      COMPUTE PERCOLATION-IF NO WATER AVAILABLE THEN SKIP
ccvk      IF((PINC+UZFWC).GT.0.01) GO TO 251
      xx1=PINC+UZFWH

cVK-----------------------------------------------------------------
cVK 2013  NOTE       changed April 2013
cVK 2013 Original SAC-SMA has problem with 0.01 limit.
cVK 2013 If pinc < 0.01mm and uzfwc close to uzfwm, uzfwc can
cVK 2013 be more than uzfwm, and there is no checking of this.
cVK 2013 Therefore, I changed the limit to 0.0
c2013??      IF(xx1.GT.0.01) GO TO 251
cVK----------------------------------------------------------------
      IF(xx1.GT.0.0) GO TO 251

      UZFWC=UZFWC+PINC
CV.K  ADD TO UNFROZEN WATER ALSO
      UZFWH=UZFWH+PINC
      GO TO 249
  251 PERCM=LZFPM*DLZP+LZFSM*DLZS
CVK      PERC=PERCM*(UZFWC/UZFWM)
CV.K  USE ONLY UNFROZEN WATER RATIOS 
ccvk  new change: PERCOLATION REDUCED BY RUZPERC 
CC       PERC=PERCM*(UZFWH/UZFWM)*RUZICE
      PERC=PERCM*(UZFWH/UZFWM)
      IF(IVERS .NE. 0) PERC=PERC*RUZPERC
C--      PERC=PERCM*(UZFWH/UZFWM)*RUZPERC
CV.K      DEFR=1.0-((LZTWC+LZFPC+LZFSC)/(LZTWM+LZFPM+LZFSM))
cvk 6/22/00      DEFR=1.0-((LZTWH+LZFPH+LZFSH)/(LZTWM+LZFPM+LZFSM))
cvk  better to keep original definition of DEFR using total water
      DEFR=1.0-((LZTWC+LZFPC+LZFSC)/(LZTWM+LZFPM+LZFSM))
C     DEFR IS THE LOWER ZONE MOISTURE DEFICIENCY RATIO
c--      FR=1.0
C     FR IS THE CHANGE IN PERCOLATION WITHDRAWAL DUE TO FROZEN GROUND.
c--      FI=1.0
C     FI IS THE CHANGE IN INTERFLOW WITHDRAWAL DUE TO FROZEN GROUND.
c--      IF (IFRZE.EQ.0) GO TO 239
c--       UZDEFR=1.0-((UZTWC+UZFWC)/(UZTWM+UZFWM))
CVK
CVK     CALL FGFR1(DEFR,FR,UZDEFR,FI)
CVK      IF( IVERS .EQ. 1) THEN
CVK  IF IVERS=1, OLD VERSION; IF IVERS=2, NEW VERS. FROST INDEX,
CVK  BUT OLD VERS. OF PERCOLAT. AND INTERFLOW REDUCTION
c--      IF( IVERS .LE. 2) CALL FGFR1(DEFR,FR,UZDEFR,FI)
      
c--      IF(IVERS .EQ. 3 .AND. FGPM(5) .GT. 0.) THEN
CVK  OPTIONAL VERSION TO ACCOUNT FOR ADDITIONAL IMPERVIOUS
CVK  AREAS EFFECTS DUE FROZEN GROUND
c--       FR=1-SURFRZ1(FGCO(1),FGPM(6),FGPM(5))
c--       FI=FR
c--      ENDIF 
      
c--  239 PERC=PERC*(1.0+ZPERC*(DEFR**REXP))*FR
  239 PERC=PERC*(1.0+SACPAR(7)*(DEFR**SACPAR(8)))
C     NOTE...PERCOLATION OCCURS FROM UZFWC BEFORE PAV IS ADDED.
CV.K      IF(PERC.LT.UZFWC) GO TO 241
      IF(PERC.LT.UZFWH) GO TO 241
C      PERCOLATION RATE EXCEEDS UZFWH.
CV.K      PERC=UZFWC
      PERC=UZFWH
C     PERCOLATION RATE IS LESS THAN UZFWH.
  241 UZFWC=UZFWC-PERC
CV.K  ADJUST UNFROZEN STORAGE ALSO  
      UZFWH=UZFWH-PERC    
C     CHECK TO SEE IF PERCOLATION EXCEEDS LOWER ZONE DEFICIENCY.
      CHECK=LZTWC+LZFPC+LZFSC+PERC-LZTWM-LZFPM-LZFSM
      IF(CHECK.LE.0.0) GO TO 242
      PERC=PERC-CHECK
      UZFWC=UZFWC+CHECK
CV.K  ADJUST UNFROZEN STARAGE ALSO
      UZFWH=UZFWH+CHECK        

  242 SPERC=SPERC+PERC
C     SPERC IS THE TIME INTERVAL SUMMATION OF PERC
C
C     COMPUTE INTERFLOW AND KEEP TRACK OF TIME INTERVAL SUM.
C     NOTE...PINC HAS NOT YET BEEN ADDED
CV.K      DEL=UZFWC*DUZ*FI
CVK  INTERFLOW ALSO REDUCED DUE FROFEN GROUND (DUZ REDUCED BY RUZICE)
CVK  ADDITIONAL REDUCTION DUE IMPERVIOUS FROZEN AREAS (FI) IS OPTIONAL
CVK  IN THE NEW VERSION. BASIC OPTION IS FI=1
c--      DEL=UZFWH*DUZ*FI
      DEL=UZFWH*DUZ
      SIF=SIF+DEL
      UZFWC=UZFWC-DEL
CV.K  ADJUST UNFROZEN STORAGE ALSO
      UZFWH=UZFWH-DEL      

cvk 2013  NOTE: Channel losses addition performes earlier in the code now
cVK 7/22/2011 --------  NEW CHANGE, CHANNEL LOSSES  ------------------ 
cVK 7/22/2011 added channel losses to percolated water incementaly; --  
cVK 7/22/2011 assumes they go directly to the lower zone  ------------
cVK 7/22/2011 CHENLOSS is in mm/s therefore it' converted into mm/dt - 
c                                                                    -
cvk 2013      IF(CHANLOSS .LT. 0.0) CHANLOSS = 0.0
cvk 2013      PERC=PERC + CHANLOSS*dtsec/NINC
c                                                                    -
cVK 7/22/2011 --------------------------------------------------------       
      
C     DISTRIBE PERCOLATED WATER INTO THE LOWER ZONES
C     TENSION WATER MUST BE FILLED FIRST EXCEPT FOR THE PFREE AREA.
C     PERCT IS PERCOLATION TO TENSION WATER AND PERCF IS PERCOLATION
C         GOING TO FREE WATER.
      PERCT=PERC*(1.0-SACPAR(14))
      xx1=PERCT+LZTWC
      IF (xx1.GT.LZTWM) GO TO 243
      LZTWC=LZTWC+PERCT
CV.K  ADJUST UNFROZEN STORAGE ALSO
      LZTWH=LZTWH+PERCT      
      PERCF=0.0
      GO TO 244
  243 PERCF=PERCT+LZTWC-LZTWM
CV.K  CHANGE UNFROZEN WATER STORAGE
      LZTWH=LZTWH+LZTWM-LZTWC  
      LZTWC=LZTWM
C
C      DISTRIBUTE PERCOLATION IN EXCESS OF TENSION
C      REQUIREMENTS AMONG THE FREE WATER STORAGES.
  244 PERCF=PERCF+PERC*SACPAR(14)
      IF(PERCF.EQ.0.0) GO TO 245
      HPL=LZFPM/(LZFPM+LZFSM)
C     HPL IS THE RELATIVE SIZE OF THE PRIMARY STORAGE
C     AS COMPARED WITH TOTAL LOWER ZONE FREE WATER STORAGE.

c VK changed to account for ZERO MAX storage
      if(LZFPM .ne. 0.) then
       RATLP=LZFPC/LZFPM
      else
       RATLP = 1.
      endif
      if(LZFSM .ne. 0.) then
       RATLS=LZFSC/LZFSM
      else
       RATLS = 1.
      endif
        
C     RATLP AND RATLS ARE CONTENT TO CAPACITY RATIOS, OR
C     IN OTHER WORDS, THE RELATIVE FULLNESS OF EACH STORAGE
      FRACP=(HPL*2.0*(1.0-RATLP))/((1.0-RATLP)+(1.0-RATLS))
C     FRACP IS THE FRACTION GOING TO PRIMARY.
      IF (FRACP.GT.1.0) FRACP=1.0
      PERCP=PERCF*FRACP
      PERCS=PERCF-PERCP
C     PERCP AND PERCS ARE THE AMOUNT OF THE EXCESS
C     PERCOLATION GOING TO PRIMARY AND SUPPLEMENTAL
C      STORGES,RESPECTIVELY.
      LZFSC=LZFSC+PERCS
CV.K      IF(LZFSC.LE.LZFSM) GO TO 246
      IF(LZFSC.LE.LZFSM) THEN
       LZFSH=LZFSH+PERCS
       GO TO 246
      ENDIF
       
      PERCS=PERCS-LZFSC+LZFSM
CV.K  ADJUST UNFROZEN STORAGE ALSO
      LZFSH=LZFSH+PERCS
            
      LZFSC=LZFSM
  246 LZFPC=LZFPC+(PERCF-PERCS)
C     CHECK TO MAKE SURE LZFPC DOES NOT EXCEED LZFPM.
CV.K      IF (LZFPC.LE.LZFPM) GO TO 245
      IF (LZFPC.LE.LZFPM) THEN
       LZFPH=LZFPH+(PERCF-PERCS)
       GO TO 245
      ENDIF
       
      EXCESS=LZFPC-LZFPM
      LZTWC=LZTWC+EXCESS
CV.K  ADJUST UNFROZEN STORAGES ALSO
      LZTWH=LZTWH+EXCESS
      LZFPH=LZFPH+(PERCF-PERCS)-EXCESS
      LZFPC=LZFPM
C
C     DISTRIBUTE PINC BETWEEN UZFWC AND SURFACE RUNOFF.
  245 IF(PINC.EQ.0.0) GO TO 249
C     CHECK IF PINC EXCEEDS UZFWM
      xx1=PINC+UZFWC
      IF(xx1.GT.UZFWM) GO TO 248
C     NO SURFACE RUNOFF
      UZFWC=UZFWC+PINC
CV.K  ADJUST UNFROZEN STORAGE ALSO
      UZFWH=UZFWH+PINC
      GO TO 249
C
C     COMPUTE SURFACE RUNOFF (SUR) AND KEEP TRACK OF TIME INTERVAL SUM.
  248 SUR=PINC+UZFWC-UZFWM
      UZFWC=UZFWM
CV.K  ADJUST UNFROZEN STORAGE ALSO
      UZFWH=UZFWH+PINC-SUR
      SSUR=SSUR+SUR*PAREA
      ADSUR=SUR*(1.0-ADDRO/PINC)
C     ADSUR IS THE AMOUNT OF SURFACE RUNOFF WHICH COMES
C     FROM THAT PORTION OF ADIMP WHICH IS NOT
C     CURRENTLY GENERATING DIRECT RUNOFF.  ADDRO/PINC
C     IS THE FRACTION OF ADIMP CURRENTLY GENERATING
C     DIRECT RUNOFF.
      SSUR=SSUR+ADSUR*ADIMP
C
C     ADIMP AREA WATER BALANCE -- SDRO IS THE 6 HR SUM OF
C          DIRECT RUNOFF.
  249 ADIMC=ADIMC+PINC-ADDRO-ADSUR  
      xx1=UZTWM+LZTWM
      IF (ADIMC.LE.xx1) GO TO 247
      ADDRO=ADDRO+ADIMC-xx1
      ADIMC=xx1
  247 SDRO=SDRO+ADDRO*ADIMP
      IF (ADIMC.LT.0.00001) ADIMC=0.0
  240 CONTINUE

C.......................................
C     END OF INCREMENTAL DO LOOP.
C.......................................

C     COMPUTE SUMS AND ADJUST RUNOFF AMOUNTS BY THE AREA OVER
C     WHICH THEY ARE GENERATED.
cfews      if(sac_etdst .eq. -1) then
cfews       EUSED=E1+E2+E3
cfews      else
       EUSED=ET1+ET2
cfews      endif  
C     EUSED IS THE ET FROM PAREA WHICH IS 1.0-ADIMP-PCTIM
      SIF=SIF*PAREA
C
C     SEPARATE CHANNEL COMPONENT OF BASEFLOW
C     FROM THE NON-CHANNEL COMPONENT
      TBF=SBF*PAREA
C     TBF IS TOTAL BASEFLOW
      BFCC=TBF*(1.0/(1.0+SACPAR(15)))
C     BFCC IS BASEFLOW, CHANNEL COMPONENT
      BFP=SPBF*PAREA/(1.0+SACPAR(15))
      BFS=BFCC-BFP
      IF(BFS.LT.0.0)BFS=0.0
      BFNCC=TBF-BFCC
C     BFNCC IS BASEFLOW,NON-CHANNEL COMPONENT
C
C     ADD TO MONTHLY SUMS.
c--      SINTFT=SINTFT+SIF
c--      SGWFP=SGWFP+BFP
c--      SGWFS=SGWFS+BFS
c--      SRECHT=SRECHT+BFNCC
c--      SROST=SROST+SSUR
c--      SRODT=SRODT+SDRO
C
C     COMPUTE TOTAL CHANNEL INFLOW FOR THE TIME INTERVAL.
      TCI=ROIMP+SDRO+SSUR+SIF+BFCC
        GRND = SIF + BFCC   ! interflow is part of ground flow
CC	GRND = BFCC         ! interflow is part of surface flow
	SURF = TCI - GRND
C
C     COMPUTE E4-ET FROM RIPARIAN VEGETATION.
	E4=(EDMND-EUSED)*SACPAR(6)
C
C     SUBTRACT E4 FROM CHANNEL INFLOW
	TCI=TCI-E4
	IF(TCI.GE.0.0) GO TO 250
	E4=E4+TCI
	TCI=0.0
cc  250 SROT=SROT+TCI
250	CONTINUE
	GRND = GRND - E4
	IF (GRND .LT. 0.) THEN
	   SURF = SURF + GRND
	   GRND = 0.
	 IF (SURF .LT. 0.) SURF = 0.
	END IF
C
C     COMPUTE TOTAL EVAPOTRANSPIRATION-TET
      EUSED=EUSED*PAREA
      TET=EUSED+E5+E4
c--      SETT=SETT+TET
c--      SE1=SE1+E1*PAREA
c--      SE3=SE3+E3*PAREA
c--      SE4=SE4+E4
c--      SE5=SE5+E5
C     CHECK THAT ADIMC.GE.UZTWC
      IF (ADIMC.LT.UZTWC) ADIMC=UZTWC
C
c  Return back SAC states
      SACST(1)=UZTWC
      SACST(2)=UZFWC      
      SACST(3)=LZTWC
      SACST(4)=LZFSC
      SACST(5)=LZFPC
      SACST(6)=ADIMC

c new change: check negative states
      do i=1,6
       if(sacst(i) .lt. -1.0) then
        write(*,*) ' SAC state#',i,'<-1.',iter,sacst(i)
        write(*,*) '       :at HRAPX HRAPY=',HRAPX,HRAPY
        stop
        ERROR = 1
        return
       endif
       if(sacst(i) .lt. 0.0) sacst(i)=0.0
      enddo

      if(uztwh .lt. 0.0) uztwh=0.0
      if(uzfwh .lt. 0.0) uzfwh=0.0
      if(lztwh .lt. 0.0) lztwh=0.0
      if(lzfsh .lt. 0.0) lzfsh=0.0
      if(lzfph .lt. 0.0) lzfph=0.0
      if(sacst(1) .lt. uztwh) uztwh=sacst(1)
      if(sacst(2) .lt. uzfwh) uzfwh=sacst(2)
      if(sacst(3) .lt. lztwh) lztwh=sacst(3)
      if(sacst(4) .lt. lzfsh) lzfsh=sacst(4)
      if(sacst(5) .lt. lzfph) lzfph=sacst(5)

c new change        
       
CVK  NEW VERSION OF FROST INDEX  ------------------------------
       IF (IVERS .NE. 0) THEN
        IF(FRZST(6) .LT. 0.) THEN
         FRZST(6)=FRZST(6)+UZTWH
        ELSE
         FRZST(6)=UZTWH
        ENDIF
        IF(FRZST(7) .LT. 0.) THEN
         FRZST(7)=FRZST(7)+UZFWH
        ELSE
         FRZST(7)=UZFWH
        ENDIF
        IF(FRZST(8) .LT. 0.) THEN
         FRZST(8)=FRZST(8)+LZTWH
        ELSE
         FRZST(8)=LZTWH
        ENDIF
        IF(FRZST(9) .LT. 0.) THEN
         FRZST(9)=FRZST(9)+LZFSH
        ELSE
         FRZST(9)=LZFSH
        ENDIF
        IF(FRZST(10) .LT. 0.) THEN
         FRZST(10)=FRZST(10)+LZFPH
        ELSE
         FRZST(10)=LZFPH
        ENDIF

        CALL FROST2_1(PXV,TA,WE,AESC,SH,FRZPAR,SACPAR,FRZST,SACST,
     +        SACST_PRV,SMC,SH2O,DTFRZ,IDTFRZ,NSOIL,NUPL,NSAC,IVERS,
     +        FRZDUP,FRZDBT,FROST,SMAX,SMCDRY, HRAPX,HRAPY)

cvk 4/2009 adjusting sac liquid water variables after frost2_1 simulations 
       UZTWH=FRZST(6)
       IF(UZTWH .LT. 0.) UZTWH=0.
       UZFWH=FRZST(7)
       IF(UZFWH .LT. 0.) UZFWH=0.
       LZTWH=FRZST(8)
       IF(LZTWH .LT. 0.) LZTWH=0.
       LZFSH=FRZST(9)
       IF(LZFSH .LT. 0.) LZFSH=0.
       LZFPH=FRZST(10)
       IF(LZFPH .LT. 0.) LZFPH=0.

c check liquid states less or equal to total
       DO I = 1, 5
         IF( FRZST( I + 5 ) .GT. SACST(I) ) THEN
             FRZST( I + 5 ) = SACST( I )
         ENDIF
       ENDDO

CVK_02  NEW OPTION TO INTERPOLATE MODEL SOIL LAYER TEMP. INTO DESIRED LAYERS
      NMOD=NSOIL+1
cvk 08/09      DSMOD(1)=0.
      DSMOD(1)=-0.5*FRZPAR(10)
      TSTMP(1)=FRZST(1)
      SWTMP(1)=SMC(2)
      SWHTMP(1)=SH2O(2)
      DSMOD(NMOD)=-FRZPAR(5)
      TSTMP(NMOD)=FRZPAR(2)-T0
      SWTMP(NMOD)=SMAX
      SWHTMP(NMOD)=SMAX
      do i=2,nmod-1
       DSMOD(I)=-0.5*(FRZPAR(I+8)+FRZPAR(I+9))
       TSTMP(I)=FRZST(I)
       SWTMP(I)=SMC(I)
       SWHTMP(I)=SH2O(I)
      ENDDO 
cc-      do ii = 1, nmod
cc-        WRITE(*,*) SWTMP(ii), SWHTMP(ii), DSMOD(ii)
cc-      ENDDO

      CALL SOIL_INT1(TSTMP,NMOD,DSMOD,DSINT,NDSINT,TSINT)
      CALL SOIL_INT1(SWTMP,NMOD,DSMOD,DSINTW,NDINTW,SWINT)
      CALL SOIL_INT1(SWHTMP,NMOD,DSMOD,DSINTW,NDINTW,SWHINT)

cvk 10/09 SM adjustment to match site
      do jjj=1,NDINTW
cfews       swint(jjj)=swint(jjj)+adjwlt
       if(swint(jjj) .lt. 0.0) swint(jjj)=0.0
      enddo

cvk  1/2008 Option to generate normalized soil moisture content (SR)
      if(NORMALIZE .eq. 1) then
       DO I=1,NDINTW
        SWINT(I) =  (SWINT(I) -  smcdry)/(SMAX-smcdry)
        SWHINT(I) = (SWHINT(I) - smcdry)/(SMAX-smcdry)
       ENDDO
      endif 
cvk  1/2008  end soil moisture normalization
       
C--      DO I=1,NINTW
C--       IF(I .EQ. 1) THEN
C--        SWINT(I)=SWINT(I)*DSINTW(I)*1000.
C--        SWHINT(I)=SWHINT(I)*DSINTW(I)*1000
C--       ELSE	
C--        SWINT(I)=SWINT(I)*(DSINTW(I)-DSINTW(I-1))*1000.
C--        SWHINT(I)=SWHINT(I)*(DSINTW(I)-DSINTW(I-1))*1000
C--       ENDIF
C--      ENDDO 	

c print results
      if(iprnt .eq. 1) then
cfews       if(sac_etdst .ne. -1.) then
        flxup=dsup+et1      
        write(*,9393) year,month,day,hour+1,pxv,ta,uztwc,uzfwc,lztwc,
     +       lzfsc,lzfpc,sperc,et1,et2,edmndbare,edmnd,etp,tet,tci,
     +       flxup,edir,shdfac,solar,z0,ch,pc,rcs,rct,rcq,rcsoil
     +                       ,(smc(ii),sh2o(ii),ii=1,nsoil)
      endif 
9393  format(i5,3i4,f8.3,f6.1,5f7.2,8f8.4,f8.3,f8.4,f6.3,f7.1,2f7.
     +       4,5f6.3,10f7.4)
9394  format(a5,3a4,a8,  a6,  5a7,  8a8,  a8,  a8,  a6,  a7, 2a7, 5a6
     +                                                      ,10a7  )         
9395  format(i5,3i4,f8.3,f6.1,5f7.2,6f8.4)
9396  format(a5,3a4,a8,  a6,  5a7,  6a8)
c           if(iter .eq. 1) stop 
cc        WRITE (*,905) iter,(sacst(ii),ii=1,6),
cc     1  SPERC,ROIMP,SDRO,SSUR,SIF,BFS,BFP,TCI,
cc     2  EDMND,TET,PXV,(frzst(II),II=6,10),WE,SH,AESC,TA,frost,
cc     3  (frzst(II),II=1,nsoil),frzdup,frzdbt             
cc        write(*,977) smax,(swint(ii),swhint(ii),tstmp(ii),ii=1,nintw)
       ELSE
cc        WRITE (*,977) (sacst(ii),ii=1,6),
cc     1  SPERC,ROIMP,SDRO,SSUR,SIF,BFS,BFP,TCI,
cc     2  EDMND,TET,PXV,UZTWH,UZFWH,LZTWH,LZFSH,LZFPH,WE,SH,AESC,TA
       ENDIF
  905 FORMAT (i10,5x,F7.2,F7.3,F7.2,F7.3,2F7.2,7F7.3,2F8.3,F7.3,
     +        F9.4,5f7.2,2f6.1,F6.3,F9.4,f7.2,8f7.2,20F7.1)
  925 FORMAT (1H ,14x,13a7,2a8,a7,a9,5a7,3a6,a9,29a7) 
  977 FORMAT (1H ,14x,f7.2,6(3F7.2))
c             ,F7.3,F7.2,F7.3,2F7.2,7F7.3,2F8.3,F7.3,
c     +        F9.4,5F7.2,2f6.1,F6.3,F9.4,8f7.2,20F7.1)

C.......................................

      RETURN
      END
