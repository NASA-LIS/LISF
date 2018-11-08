C***********************************************************************
C****************     V. KOREN    03/2003      *************************             
C**  SUBROUTINE WAS CHANGED TO BE IMPLEMENTED INTO GRID BASED HL-RMS  **
C***********************************************************************
C
C MEMBER PACK19
C  (from old member FCPACK19)
C
ck      SUBROUTINE PACK19(KDA,KHR,NDT,TA,PX,PCTS,RSL,OWE,OSC,PGM,RM,TWE,
C                             LAST UPDATE: 06/22/95.14:05:09 BY $WC30EA
ck     1COVER,CWE,CAESC,IFUT,IDT,IBUG,IDN,IMN,IYR,IOUTYP,OPNAME)
C.......................................
C     THIS SUBROUTINE EXECUTES THE 'SNOW-17 ' OPERATIONAL FOR ONE
C        COMPUTATIONAL PERIOD.
C.......................................
C     SUBROUTINE INITIALLY WRITTEN BY...
C        ERIC ANDERSON - HRL   MAY 1980

C        UPDATED 4/15/00 BY V. KOREN TO ADD SNOW DEPTH CALCULATIONS
C.......................................
ck      DIMENSION PX(1),PCTS(1),RM(1),OPNAME(2)
C     COMMON BLOCKS
ck      INCLUDE 'common/ionum'
ck     INCLUDE 'common/fprog'
ck      COMMON/FDBUG/IODBUG,ITRACE,IDBALL,NDEBUG,IDEBUG(20)
ck      COMMON/FSNWUP/IUPWE,IUPSC
ck      COMMON/SNPM19/ALAT,SCF,MFMAX,MFMIN,NMF,UADJ,SI,MBASE,PXTEMP,
ck     1   PLWHC,TIPM,PA,ADC(11),LMFV,SMFV(12),LAEC,NPTAE,AE(2,14)
ck      COMMON/SNCO19/WE,NEGHS,LIQW,TINDEX,ACCMAX,SB,SBAESC,SBWS,
C
CVK  ADDED TWO MORE STATES: SNDPT & SNTMP
CVK     1   STORGE,AEADJ,NEXLAG,EXLAG(7)
ck     1   STORGE,AEADJ,NEXLAG,EXLAG(7),SNDPT,SNTMP
    
ck      COMMON/SUMS19/SPX,SSFALL,SRM,SMELT,SMELTR,SROBG,DSFALL,DRAIN,
ck     1 DQNET,DRSL,NDRSP
ck      COMMON/SNUP19/MFC,SFALLX,WINDC,SCTOL,WETOL,SNOF,UADJC
ck---------------------------------------------------------------------------
C
C    ================================= RCS keyword statements ==========
ck      CHARACTER*68     RCSKW1,RCSKW2
ck      DATA             RCSKW1,RCSKW2 /                                 '
ck     .$Source: /fs/hsmb5/hydro/CVS_root/gorms/snow17/pack19_rms.f,v $
ck     . $',                                                             '
ck     .$Id: pack19_rms.f,v 3.0 2009-12-30 16:49:07 zcui Exp $
ck     . $' /
C    ===================================================================
C
ck  KDA - day number for ICP print
ck  KHR - hour number for ICP print
ck  NDT - number of simulation (or precipitation) time intervals in air
ck        temperature time period, IDT  
ck  PGM - ground melt, mm per temperature time interval which is calculated 
ck       from daily ground melt parameter, and in PACK routine will be 
ck       converted into mm per precipitation time interval that=sim. time step
ck  IFUT - variable that defines the end time of MOD operation, and it involvs
ck         a few input (in this program from COMMON/SNUP19/) variables of MOD 
ck         operation: SFALLS, WINDC
ck  IDT - time interval for air temperature data. HL-RMS version will not use
ck        it here, but instead may use it in 'do_snow17' subroutine.    
ck  IBUG - index to print
ck  IDN, IMN & IYR - day, month & yera to print in ICP
ck  IOUTYP - index to print ICP information
ck  OPNAME - file name for ICP print
ck---------------------------------------------------------------------------
ck
ck  HL-RMS version  3/2003  VK
ck
ck  Variables that control ICP & Debug printing are excluded:
ck	KDA, KHR, IBUG, IYR, IOUTYP, OPNAME
ck  NDT & IDT are excluded from the parameter list assuming that 'do_snow17
ck          will make temperature and precip data over the same time interval
ck          that equals simulation time interval DT
ck  PGM will be in mm/day, and in subroutine will be converted into mm/DT
ck  IFUT will be 1 in this version assuming that there will be no MODs
ck  PX(1),PCTS(1),RM(1) arrays replaced by just simple variables
ck  All COMMON blocks are excluded, and needed variables are moved to 
ck             SUBROUTINE parameters
ck  Variables of cumulated water and heat balance in COMMON/SUMS19/ are not 
ck             included in Subroutine parameters, however they can be 
ck             outputed in the Subroutine
ck---------------------------------------------------------------------------
ck
      SUBROUTINE PACK19(IDD,IMN,DT,TA,PX,PCTS,RM,TWE,COVER,WE,NEGHS,
     +                     LIQW,TINDEX,ACCMAX,SNDPT,SNTMP,SB,SBAESC,
     +                      SBWS,STORGE,AEADJ,EXLAG,NEXLAG,ALAT,SCF,
     +                         MFMAX,MFMIN,NMF,UADJ,SI,MBASE,PXTEMP,
     +                        PLWHC,TIPM,GM,PA,LAEC,ADC,SNOF,DS,DTA,
     +                                                       SWITCH)

c
c----------------------------------------------------------------------------
ck  IDD & IMN - Day (in original SUBR. it was Julian day) and Month
ck   DT       - Simulation time step, in hours 
ck Inputs:
ck  TA & PX - air temperature and precipitation 
ck  PCTS - snow fraction of precipitation: if <0, no data, and 
ck         form of PX is defined based on 1) rain-snow elevation if LAEC 
ck         not =0, or 2) temperature if LAEC =0.
ck  RSL - for rain-snow elevation option only. If not useed, <0. 
ck  OWE & OSC - WE and % Snow cover for MOD function. (if MOD not used they
ck              non-defined,  <0)
c
ck Outputs: 
ck  RM - rain+melt (array replaced by a single value)
ck  TWE & COVER - final results of simulation (TWE total including liqw)
c
ck Snowmelt parameters from COMMON/SNPM19/ defined in EX19:
ck  ALAT - latitude
ck  SCF,MFMAX,MFMIN,NMF,UADJ,SI,MBASE,PXTEMP,PLWHC,TIPM - defined SNOW-17
ck         parameters
ck  GM - ground melt, mm/DT
ck  PA - Standard air presure; 
ck  ADC(11) - Array of snow cover fraction;
ck  LMFV - Melt factor variation type: =0 & ALAT<54, Default variation,
ck                                     =0 & ALAT>54, Alaska type
ck                                      NOT 0,       User specified
ck  SMFV(12) - NO in THIS VERSION: user specified melt factor variation 
ck  LAEC - Rain-snow split definition: =0, Based on temperature threshold,
ck                       not =0, Based on rain-snow elevation (array AE)
ck  NPTAE - NO in THIS VERSION: Number of elevation values in array(table) AE
ck  AE - NO in THIS VERSION: rain-snow elevation for form of precip. 
ck       In original version two dimensional array (2,14).      
c
ck States:
ck  WE     - Snow water equivalent (without liquid water), mm
ck  NEGHS  - Negative snow heat, mm
ck  LIQW   - Liquid water , mm
ck  TINDEX - Antecedent temperature index, Celsius
ck  ACCMAX - Cummulated snow water including liquid, mm
ck  SB, SBAESC, SBWS - Internal snow states during melt & new snow fall
ck  STORGE - Snow liquid water attenuation storage, mm
ck  AEADJ  - Adjusted areal snow cover fraction
ck  NEXLAG - Number of ordinates in lagged liquid water array (EXLAG)
ck  EXLAG(7) - Array of lagged liquid water values
ck  SNDPT  - Snow depth, cm
ck  SNTMP  - Average snow temperature, Celsius
c
ck MOD operation variables (should be defined in 'do_snow17' subroutine:
ck  MFC,SFALLX,WINDC,SCTOL,WETOL,SNOF,UADJC
ck  not available in this version although UADJC and states could be changed 
ck  from input deck 
c----------------------------------------------------------------------------
CVK  SWITCH VARIABLE CHANGE KIQUID WATER FREEZING VERSION
CVK  SWITCH=0, VICTOR'S VERSION
CVK  SWITCH=1, ERIC'S VERSION
c  
      INTEGER SWITCH
      REAL MFMAX,MFMIN,NMF,LIQW,NEGHS,MBASE,MELT,LIQWMX
      REAL ADC(11),EXLAG(7),SMFV(12)/12*-1.0/
      INTEGER IDANG(12)/285,316,345,10,40,71,101,132,163,193,224,254/

cvk 04/19/07 
      SAVE LMFV
C
cvk  test only. will be replaced by subroutine parameter
cvk      SWITCH = 0
      
CVK  CHANGES -------------------------------------------
cc      SAVE ISTRT,TPREV,DS
cc      INTEGER ISTRT/0/      

ck Convert date into number days from March 21
      IDN = IDANG(IMN) + IDD
ck  IFUT=1: no MOD option in this version
ck      IFUT = 1
ck  LMFV = 0: NO USER SPECIFIED MONTHLY MELT FACTOR. Array SMFV not defined
      LMFV = 0
            
CVK   CALCULATE AIR TEMPERATURE CHANGE, DTA
CVK   AND KEEP PREVIOUS TIME STEP TEMPERATURE
cc      IF(ISTRT .EQ. 0) THEN
cc       TPREV=TA
cc       IF(SNDPT .GT. 0.0) DS=0.1*WE/SNDPT
cc       ISTRT=1
cc      ENDIF 
cc      IF(SNDPT .LE. 0.0) THEN
cc       DTA=TA
cc      ELSE 
cc       DTA=TA-TPREV
cc       IF(TPREV .GT. 0. .AND. TA .GT. 0.) DTA=ABS(DTA)
cc       IF(TPREV .GT. 0. .AND. TA .LT. 0.) DTA=TA
cc      ENDIF 
cc      TPREV=TA
CVK-----------------------------------------------------
       
C     CONSTANTS
C     IF SNOWFALL EXCEEDS SNEW/HR--TINDEX=TPX
      SNEW=1.5
C     IF RAIN EXCEEDS RMIN/HR--USE RAIN-ON-SNOW MELT EQUATION
      RMIN=0.25
C     SBC=STEFAN/BOLTZMAN CONSTANT--MM/(((DEGK/100)**4)*HR)
      SBC=.0612

C     INITIAL VALUES
ck  Basic simulation time step is DT (hr)
ck  Recalculations into precip. time step excluded, all data are consistent  
ck  100 ITPX=IDT/NDT
ck      FITPX=ITPX
ck      FNDT=NDT
ck      GM=PGM/FNDT
ck      SFNEW=SNEW*FITPX
ck      RFMIN=RMIN*FITPX
ck      SBCI=SBC*FITPX

      SFNEW=SNEW*DT
      RFMIN=RMIN*DT
      SBCI=SBC*DT
      MC=0
ck      PSFALL=0.0
      PRAIN=0.0
      PQNET=0.0
      PSNWRO=0.0
      PROBG=0.0

CVK ----     V.I.K.  04/10/00  --------- 
ck      SXFALL=0.0
ck      SXMELT=0.0
ck      SXGSLOS=0.0
      SXRFRZ=0.0
C --------------------------------------     
C     CYCLE THROUGH THE COMPUTATIONAL PERIOD FOR EACH PRECIPITATION
C        TIME INTERVAL
ck      DO 200 I=1,NDT
ck      PXI=PX(I)
      PXI=PX      
      IF((PXI.EQ.0.0).AND.(WE.EQ.0.0)) GO TO 160
      SFALL=0.0
      CNHSPX=0.0
      RAIN=0.0
      RAINM=0.0
      IF(PXI.EQ.0.0) GO TO 110
C.......................................
C     DETERMINE FORM OF PRECIP. AND ACCUMULATE SNOW COVER IF SNOW.
ck      PCT=PCTS(I)
      PCT=PCTS
      IF(PCT.GT.1.0) PCT=1.0
      TPX=TA
      IF(PCT.LT.0.0) GO TO 102
C
C     FORM OF PRECIP. INPUT.
      FRACS=PCT
      FRACR=1.0-PCT
      GO TO 105
  102 IF (LAEC.EQ.0) GO TO 104
      GO TO 105
C
C     FORM OF PRECIP. BASED ON TEMPERATURE.
  104 IF(TPX.GT.PXTEMP) GO TO 103
C
C     SNOW
      FRACS=1.0
      FRACR=0.0
      GO TO 105
C
C     RAIN
  103 FRACS=0.0
      FRACR=1.0
  105 IF(FRACS.EQ.0.0) GO TO 109
C
C     ACCUMULATE SNOWFALL
      TS=TPX
      IF(TS.GT.0.0) TS=0.0
      SFALL=PXI*FRACS*SCF
ck      IF(IFUT.EQ.0) SFALL=SFALL*SFALLX
ck  Disabling cumulated variables to store (they do not used directly)
ck      SPX=SPX+SFALL
ck      SSFALL=SSFALL+SFALL
ck      DSFALL=DSFALL+SFALL
ck      PSFALL=PSFALL+SFALL
C
CEA  ***  Start of Eric changes  *******
CEA   REVISED CODE FOR ADJUSTING AESC STATES FOR NEW SNOWFALL
cc      xxs=we+liqw
      WELIQW = WE+LIQW
      IF (WELIQW .LT. SBWS) GO TO 106
C     WELIQW>SBWS - 100% COVER - SNOW ON BARE GROUND
ck  IF statements with expression do not work properly in Linux
ck  because register size does not match real type size
cck      IF ((WE+LIQW).LT.SBWS) GO TO 106
      SBWS=SBWS+0.75*SFALL
CEA      IF((SFALL.GE.SNOF).AND.(SB.GT.xxs)) SB=WE+LIQW
cck      IF((SFALL.GE.SNOF).AND.(SB.GT.WE+LIQW)) SB=WE+LIQW
      GO TO 107
CEA  106 IF(SFALL.GE.SNOF) SBWS=WE+LIQW+0.75*SFALL
  106 IF (WELIQW .GT. SB) GO TO 1061
C     ON DEPLETION CURVE OR .GE. AI
C     CHECK IF SUFFICIENT NEW SNOW TO LEAVE DEPLETION CURVE
      IF (SFALL .GE. SNOF) GO TO 1062
C     REMAIN ON DEPLETION CURVE
      SB=SB+SFALL
      SBWS=SB
C     SBAESC UPDATED IN AESC19
      GO TO 107
C     WELIQW>SB AND WELIQW LE SBWS - <100% COVER - SNOW ON BARE GROUND
 1061 IF (SFALL .GE. SNOF) GO TO 1062
      SBWS=SBWS+0.75*SFALL
      GO TO 107
 1062 SBWS=WELIQW+0.75*SFALL
CEA
CEA ***  End of Eric changes  ****
c
  107 WE=WE+SFALL  
C     IF WE+LIQW.GE.3*SB, ASSUME NEW ACCUMULATION PERIOD
      xxs=we+liqw
      xcst=3.0*SB
      IF(xxs .LT. xcst) GO TO 108
cck      IF(WE+LIQW.LT.3.0*SB) GO TO 108
      ACCMAX=WE+LIQW
      AEADJ=0.0
  108 CNHSPX=-TS*SFALL/160.0
      IF(SFALL.GT.SFNEW) TINDEX=TS
C
C     RAINFALL AND RAIN MELT.
  109 RAIN=PXI*FRACR
ck      SPX=SPX+RAIN
      PRAIN=PRAIN+RAIN
      IF(WE.EQ.0.0) GO TO 160
      DRAIN=DRAIN+RAIN
      TR=TPX
      IF(TR.LT.0.0) TR=0.0
      RAINM=0.0125*RAIN*TR
C.......................................
C     MELT AT GROUND-SNOW INTERFACE
  110 IF(WE.GT.GM) GO TO 111
      GMRO=WE+LIQW
      MELT=0.0
      ROBG=RAIN
      RAIN=0.0
ck      SROBG=SROBG+ROBG
      GO TO 150
  111 GMWLOS=(GM/WE)*LIQW
      GMSLOS=GM
C.......................................
C     COMPUTE SURFACE ENERGY EXCHANGE FOR THE COMPUTATIONAL PERIOD BASED
C        ON 100 PERCENT COVER AND NON-RAIN CONDITIONS -
      IF(MC.EQ.1) GO TO 115
      CALL MELT19(IDN,IMN,ALAT,TA,MELT,MFMAX,MFMIN,MBASE,TINDEX,TIPM,
     1   CNHS,NMF,LMFV,SMFV)
      MC=1
C.......................................
C     DETERMINE MELT FOR THE TIME INTERVAL - SURFACE ENERGY EXCHANGE
C        IS UNIFORM DURING THE COMPUTATIONAL PERIOD.
  115 CONTINUE
ck  115 CNHS=PCNHS/FNDT
      NR=1
      IF(RAIN.GT.RFMIN) GO TO 120
C
C     NON-RAIN OR LIGHT DIZZLE INTERVAL
ck      MELT=PMELT/FNDT
      MELT=MELT
ck  Excluding adjustment MOD
ck      MELT=MELT*MFC
      MELT=MELT+RAINM
      GO TO 130
C
C     RAIN INTERVAL.
  120 EA=2.7489E8*EXP(-4278.63/(TA+242.792))
C     ASSUME 90 PERCENT RELATIVE HUMIDITY DURING RAIN-ON-SNOW
      EA=0.90*EA
      TAK=(TA+273)*0.01
      TAK4=TAK*TAK*TAK*TAK
      QN=SBCI*(TAK4-55.55)
C
C     UADJC IS UADJ MOD MULTIPLIER added by mike smith 2/12/97
ck    UADJC and WINDC adjustment excluded in this version
      
      QE=8.5*(EA-6.11)*UADJ
ck      QE=8.5*(EA-6.11)*UADJ*UADJC
ck        IF(IFUT.EQ.0) QE=QE*WINDC
      QH=7.5*0.000646*PA*UADJ*TA
ck      QH=7.5*0.000646*PA*UADJ*TA*UADJC
ck        IF(IFUT.EQ.0) QH=QH*WINDC
      MELT=QN+QE+QH+RAINM
      IF(MELT.LT.0.0) MELT=0.0
      NR=0
C.......................................
C     COMPUTE AREAL EXTENT OF SNOW COVER BASED ON CONDITIONS AT THE
C        BEGINNING OF THE TIME INTERVAL ADJUSTED FOR NEW SNOWFALL.
C
  130 CALL AESC19(WE,LIQW,ACCMAX,SB,SBAESC,SBWS,SI,ADC,AEADJ,AESC,SNOF)
C
C     ADJUST VALUES FOR AESC.
      IF(AESC.EQ.1.0) GO TO 134
      MELT=MELT*AESC
      CNHS=CNHS*AESC
      GMWLOS=GMWLOS*AESC
      GMSLOS=GMSLOS*AESC
C.......................................
C     COMPUTE RAIN FALLING ON BARE GROUND
      ROBG=(1.0-AESC)*RAIN
      RAIN=RAIN-ROBG
      GO TO 135
  134 ROBG=0.0
C.......................................
C     COMPUTE SUM AND CHECK CNHS.
  135 CONTINUE
ck  135 SROBG=SROBG+ROBG
      xxs=cnhs+neghs 
      IF(xxs .LT. 0.0) CNHS=-1.0*NEGHS
cck      IF((CNHS+NEGHS).LT.0.0) CNHS=-1.0*NEGHS

CVK  CUMULATE FOR TIME PERIOD
ck      SXFALL=SXFALL+SFALL
ck      SXMELT=SXMELT+MELT
ck      SXGSLOS=SXGSLOS+GMSLOS
CVK -------------------------

C.......................................
C     ADJUST WE FOR SURFACE AND GROUND MELT
C     GROUND MELT
      WE=WE-GMSLOS
      LIQW=LIQW-GMWLOS
      GMRO=GMSLOS+GMWLOS
C
C     SURFACE MELT
      IF(MELT.LE.0.0) GO TO 137
      IF(MELT.LT.WE) GO TO 136
      MELT=WE+LIQW
      QNET=MELT
ck      DQNET=DQNET+QNET
ck      IF(NR.EQ.1) SMELT=SMELT+MELT
ck      IF(NR.EQ.0) SMELTR=SMELTR+MELT
      GO TO 150
  136 WE=WE-MELT
C     QNET=NET SURFACE ENERGY EXCHANGE IN MILLIMETERS WE.
  137 QNET=MELT-CNHS-CNHSPX
ck      DQNET=DQNET+QNET
ck      PQNET=PQNET+QNET
ck      IF(NR.EQ.1) SMELT=SMELT+MELT
ck      IF(NR.EQ.0) SMELTR=SMELTR+MELT
ck      IF(IBUG.EQ.2) CALL PRCO19
C.......................................
C     PERFORM HEAT AND WATER BALANCE FOR THE SNOW COVER.
      WATER=MELT+RAIN
      HEAT=CNHS+CNHSPX
      LIQWMX=PLWHC*WE
      NEGHS=NEGHS+HEAT
C     TEMPERATURE OF SNOW CAN NOT BE BELOW-52.8 DEGC
      IF(NEGHS.LT.0.0) NEGHS=0.0
      xcst=0.33*WE
      xxs1=WATER+LIQW
      xxs2=LIQWMX+NEGHS+PLWHC*NEGHS
      IF(NEGHS.GT.xcst) NEGHS=0.33*WE
      IF(xxs1 .LT. xxs2) GO TO 140
cck      IF(NEGHS.GT.0.33*WE) NEGHS=0.33*WE
cck      IF(WATER+LIQW.LT.(LIQWMX+NEGHS+PLWHC*NEGHS)) GO TO 140
C
C     EXCESS WATER EXISTS.
      EXCESS=WATER+LIQW-LIQWMX-NEGHS-PLWHC*NEGHS
      LIQW=LIQWMX+PLWHC*NEGHS
      WE=WE+NEGHS
C
CEA  CUMULATE REFROZEN WATER - PREVIOUSLY IGNORED AT THIS POINT
      SXRFRZ=SXRFRZ+NEGHS
      NEGHS=0.0
      GO TO 145
c
CVK  ***  Switch between 2 versions of freezing liquid water in snow:
CVK  - Eric NWSRFS version (switch=1), in input deck operation name snow171
CVK  - Victor version (switch=0), in input deck operation name snow17
CVK
  140 IF(SWITCH .EQ. 0) THEN
CVK  VICTOR'S VERSION
       XXS = WATER+LIQW
       IF(XXS .LT. NEGHS) GOTO 841
C  WATER+LIQW EXCEEDS NEGHS - LIQUID WATER CONTENT IS CHANGED
C  - LIQW INCREASES IF WATER > NEGHS
C  - LIQW DECREASES IF WATER < NEGHS AND NEGHS > 0.0               
       LIQW=LIQW+WATER-NEGHS
       WE=WE+NEGHS
C
CVK  CUMULATE REFROZEN WATER
       SXRFRZ=SXRFRZ+NEGHS
       NEGHS=0.0
       EXCESS=0.0
       GO TO 145
C
C   ALL WATER IS REFROZEN IN THE SNOW COVER.
  841  WE=WE+WATER+LIQW
       NEGHS=NEGHS-WATER-LIQW
CVK  CUMULATE REFROZEN WATER
       SXRFRZ=SXRFRZ+WATER+LIQW
       LIQW=0.0
       EXCESS=0.0
C
CEA END OF SECTION TESTING KOREN'S MODIFICATION
      ELSE
CVK ERIC'S VERSION
       IF(WATER.LT.NEGHS) GO TO 141
C     WATER EXCEEDS NEGHS - LIQUID WATER CONTENT IS INCREASED.
       LIQW=LIQW+WATER-NEGHS
       WE=WE+NEGHS
C
CVK  CUMULATE REFROZEN WATER
       SXRFRZ=SXRFRZ+NEGHS
       NEGHS=0.0
       EXCESS=0.0
       GO TO 145
C
C     ALL WATER IS REFROZEN IN THE SNOW COVER.
  141  WE=WE+WATER
       NEGHS=NEGHS-WATER
       EXCESS=0.0
C
CVK  CUMULATE REFROZEN WATER
       SXRFRZ=SXRFRZ+WATER
      ENDIF     
C
CVK********************************************************************
CVK  ***  OLD ONLY VICTOR'S CODE  ***
ck   V. Koren:  Correction of refreezing water condition. Instead of
ck             WATER < NEGHS IT USES WATER+LIQW < NEGHS because not
ck             just time interval WATER may be refreezed but cummulated
ck             before LIQW too     
ck  140 IF(WATER.LT.NEGHS) GO TO 141
CK  140 xxs=WATER+LIQW
CK      IF(xxs.LT.NEGHS) GO TO 141
cck  140 IF(WATER+LIQW.LT.NEGHS) GO TO 141
ck end correction
C     WATER EXCEEDS NEGHS - LIQUID WATER CONTENT IS INCREASED.
CK      LIQW=LIQW+WATER-NEGHS
CK      WE=WE+NEGHS
CVK  CUMULATE REFROZEN WATER
CK      SXRFRZ=SXRFRZ+NEGHS
CK      NEGHS=0.0
CK      EXCESS=0.0
CK      GO TO 145
C     ALL LIQUID WATER IS REFROZEN IN THE SNOW COVER.
ck   V. Koren: Correction of freezing liquid water (see comment above)
ck  141 WE=WE+WATER
ck      NEGHS=NEGHS-WATER
CK  141 WE=WE+WATER+LIQW
CK      NEGHS=NEGHS-WATER-LIQW
CK      SXRFRZ=SXRFRZ+WATER+LIQW
CK      LIQW=0.0
ck end correction
CK      EXCESS=0.0
CVK  CUMULATE REFROZEN WATER
CVK      SXRFRZ=SXRFRZ+WATER
C
CVK  ***   END ONLY VICTOR'S CODE  ***
CVK************************************************************************
C      
C     IF NO NEGATIVE HEAT - TINDEX MUST BE 0.0.
  145 IF(NEGHS.EQ.0.0) TINDEX=0.0
        
C.......................................
C     ROUTE EXCESS WATER THROUGH THE SNOW COVER.
      CALL ROUT19(DT,EXCESS,WE,AESC,STORGE,NEXLAG,EXLAG,PACKRO)
C.......................................
C     ADD GROUNDMELT RUNOFF TO SNOW COVER OUTFLOW.
      PACKRO=PACKRO+GMRO
      GO TO 190
C.......................................
C     SNOW GONE - SET ALL CARRYOVER TO NO SNOW CONDITIONS.
  150 TEX=0.0
      DO 151 N=1,NEXLAG
  151 TEX=TEX+EXLAG(N)
      PACKRO=GMRO+MELT+TEX+STORGE+RAIN
      CALL ZERO19(WE,NEGHS,LIQW,TINDEX,ACCMAX,SB,SBAESC,SBWS,
     +            STORGE,AEADJ,SNDPT,SNTMP,NEXLAG,EXLAG)  
      AESC=0.0
      GO TO 190
C.......................................
C     NO SNOW COVER - NO NEW SNOWFALL.
  160 ROBG=PXI
      PACKRO=0.0
      AESC=0.0
ck      SROBG=SROBG+ROBG
C.......................................
C     COMPUTE RAIN+MELT
ck  190 RM(I)=PACKRO+ROBG
ck      SRM=SRM+RM(I)
  190 RM=PACKRO+ROBG
ck      SRM=SRM+RM
ck      PSNWRO=PSNWRO+PACKRO
ck      PROBG=PROBG+ROBG
ck  200 CONTINUE
C     END OF COMPUTATIONAL PERIOD
C.......................................

C.......................................
C     SET SIMULATED AESC AND TOTAL WATER-EQUIVALENT.
      TEX=0.0
      DO 210 N=1,NEXLAG
  210 TEX=TEX+EXLAG(N)
CVK --  V.KOREN  04/05/00   ---------------------------------
CVK   CALL SNDEPTH SUBROUTINE TO CALCULATE SNOW DEPTH
      IF(WE .GT. 0.) THEN
       SLIQ=LIQW+TEX+STORGE
c07       SFALL=SFALL-MELT
c07       IF(SFALL .LT. 0.) SFALL=0.0
       CALL SNDEPTH(WE,SLIQ,SFALL,MELT,GMSLOS,SXRFRZ,TA,DTA,DT,
     +              SNDPT,DS,SNTMP)
      ELSE
       SNDPT=0.
       SNTMP=0.
       DS=0.1
      ENDIF

CVK ---------------------------------------------------------      
      TWE=WE+LIQW+TEX+STORGE
      IF (TWE.EQ.0.0) GO TO 215

C
C     COMPUTE AREAL EXTENT BASED ON CONDITIONS AT THE END OF THE PERIOD.
      CALL AESC19(WE,LIQW,ACCMAX,SB,SBAESC,SBWS,SI,ADC,AEADJ,AESC,SNOF)
  215 COVER=AESC

C
C     STORE VALUES SO COMPUTED VALUES WILL BE AVAILABLE FOR PRINTOUT
C        EVEN IF UPDATING OCCURS.
ck      CWE=TWE
ck      CAESC=COVER
C.......................................
C     UPDATING SECTION
ck      IF((OWE.LT.0.0).AND.(OSC.LT.0.0)) GO TO 280
ck      CALL UPDT19(OWE,OSC,TWE,COVER,IUPWE,IUPSC,WETOL,SCTOL,
ck     1  SI,ADC)
C.......................................
C     PASS VALUES TO GRAPHICS INTERFACE IF REQUESTED.
  280 CONTINUE
ck  280 IF (IOUTYP.EQ.0) GO TO 290
ck      LIQWMX=PLWHC*WE
ck      CALL ICP19(OPNAME,IMN,IYR,KDA,KHR,PRAIN,PSFALL,RSL,TA,
ck     1  PQNET,PSNWRO,PROBG,NEGHS,LIQW,LIQWMX,TWE,COVER,
ck     2  OWE,OSC,WE,TEX,STORGE,ACCMAX,SB,SBAESC,SBWS,
ck     3  AEADJ,TINDEX,MAINUM)
C.......................................

C     CHECK FOR DEBUG OUTPUT.
ck  290 IF(IBUG.EQ.0) RETURN
ck      WRITE(IODBUG,904) TWE,COVER,CWE,CAESC
ck  904 FORMAT(1H0,5X,17HOUTPUT DATA--TWE=,F6.1,2X,6HCOVER=,F4.2,2X,4HCWE=
ck     1,F6.1,2X,6HCAESC=,F4.2)
ck      WRITE(IODBUG,905) (RM(I),I=1,NDT)
ck  905 FORMAT(1H ,5X,3HRM=,24F5.1)
ck      IF ((TWE.NE.CWE).OR.(COVER.NE.CAESC)) CALL PRCO19
C.......................................
      RETURN
      END
