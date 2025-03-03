!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
      module clsmf25_model

!**** "Cleaned-up" version: July 2002 (rk)
!   $Id: catchment.F90,v 1.48.6.7.6.1 2011-04-01 19:37:11 smahanam Exp $      
!rr   ---------------------------
!rr
!rr   this version:!rr   catchment.f as from /land/koster/src on May 29, 2003
!rr   snow_may2003.f as from Stephen Dery by email on May 29, 2003
!rr   
!rr   several snow code bug fixes
!rr   snow water equivalent in [mm]=[kg/m2] *throughout*
!rr   order of arguments changed in snowrt
!rr   
!rr   additional modifications:
!rr   removed "mxnch" AND "mxchp"
!rr   commented out a bunch of write statements
!rr   - reichle, 29 May 03
!rr   ---------------------------

! koster+reichle, 13 Aug 2008: added optional output of tcX and qaX as 
!                               computed before AGCM adjustments
! reichle, 28 Oct 2010 - moved DZ, SHR, PHI, FSN, FWETL, FWETC to module catch_constants
!                      - renamed N_gndtmp -> N_gt
! reichle, 28 Oct 2010 - moved SURFLAY to GEOS_CatchGridComp, pass into catchment()
! reichle, 23 Nov 2010 - replaced PHIGT with POROS(N), ALHMGT with ALHM
! reichle, 30 Nov 2010 - zero-diff revisions and clean-up for off-line (land-only) MERRA
!                        replay capability
!                         - restored PHIGT, ALHMGT 
!                         - moved MIN_SNOW_MASS->MINSWE, DZ1MAX, and SATCAPFR to 
!                            catch_constants
!                         - moved "small" back from catch_constants() into snowrt()
        use LIS_coreMod
        use LIS_logMod
      use clsmf25_constants
      use clsmf25_sibalb_coeff

      IMPLICIT NONE

      private
      public clsmf25, partition, sibalb, rzequil, get_tf0d,                  &
                   calc_soil_moist, get_tf_nd

      real, parameter :: FSN = 1.e3*ALHMGT ! unit change J/kg -> J/m  
      
      CONTAINS

      SUBROUTINE CLSMF25 (                                                     &
                     NCH, DTSTEP, SFRAC,                                       &
                     cat_id,ITYP,DZSF,TRAINC,TRAINL, TSNOW, UM,                &
                     ETURB1, DEDQA1, DEDTC1, HSTURB1,DHSDQA1, DHSDTC1,         &
                     ETURB2, DEDQA2, DEDTC2, HSTURB2,DHSDQA2, DHSDTC2,         &
                     ETURB4, DEDQA4, DEDTC4, HSTURB4,DHSDQA4, DHSDTC4,         &
                     ETURBS, DEDQAS, DEDTCS, HSTURBS,DHSDQAS, DHSDTCS,         &
                     TM, QM, ra1, ra2, ra4, raS, SUNANG, PARDIR, PARDIF,       &
                     SWNETF,SWNETS,  HLWDWN, PSUR,  ZLAI,   GREEN,  Z2,        &
                     SQSCAT, RSOIL1, RSOIL2,   RDC,                            &
                     QSAT1, DQS1, ALW1, BLW1,  QSAT2, DQS2, ALW2, BLW2,        &
                     QSAT4, DQS4, ALW4, BLW4,  QSATS, DQSS, ALWS, BLWS,        &
                     BF1, BF2, BF3,VGWMAX,                                     &
                     CDCR1,CDCR2, psis, bee, poros, wpwet, cond, gnu,          &
                     ARS1,ARS2,ARS3,ARA1,ARA2,ARA3,ARA4,ARW1,ARW2,ARW3,ARW4,   &
                     tsa1,tsa2,tsb1,tsb2,atau,btau,BUG,                        &
                     TC1, TC2, TC4, QA1, QA2, QA4, CAPAC,                      &
                     CATDEF, RZEXC, srfexc, GHTCNT, TSURF,                     &
                     WESNN, HTSNNN, SNDZN,     EVAP, SHFLUX, RUNOFF,           &
                     EINT, ESOI, EVEG, ESNO,  BFLOW,RUNSRF,SMELT,              &
                     HLWUP,SWLAND,HLATN,QINFIL,AR1, AR2, RZEQ,                 &
                     GHFLUX, TPSN1, ASNOW0, TP1, TP2, TP3, TP4, TP5, TP6,      &
                     sfmc, rzmc, prmc, entot, wtot, WCHANGE, ECHANGE, HSNACC,  &
                     EVACC, SHACC,                                             &
                     SH_SNOW, AVET_SNOW, WAT_10CM, TOTWAT_SOIL, TOTICE_SOIL,   &
                     LHACC, TC1_0, TC2_0, TC4_0, QA1_0, QA2_0, QA4_0, fices_out,&
                     !ag(02Jan2021)
                     rivsto,fldsto,fldfrc)
        
!      use  clsmf25_diagn_routines

      IMPLICIT NONE

! -----------------------------------------------------------
!     INPUTS

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP, cat_id

      REAL, INTENT(IN) :: DTSTEP, SFRAC
      REAL, INTENT(IN), DIMENSION(NCH) :: DZSF, TRAINC, TRAINL, TSNOW,  UM,    &
                     ETURB1, DEDQA1, DEDTC1, HSTURB1,DHSDQA1, DHSDTC1,         &
                     ETURB2, DEDQA2, DEDTC2, HSTURB2,DHSDQA2, DHSDTC2,         &
                     ETURB4, DEDQA4, DEDTC4, HSTURB4,DHSDQA4, DHSDTC4,         &
                     ETURBS, DEDQAS, DEDTCS, HSTURBS,DHSDQAS, DHSDTCS,         &
                     TM, QM, ra1, ra2, ra4, raS, SUNANG, PARDIR, PARDIF,       &
                     SWNETF,SWNETS,  HLWDWN, PSUR,  ZLAI,   GREEN,  Z2,        &
                     SQSCAT, RSOIL1, RSOIL2,   RDC,                            &
                     QSAT1, DQS1, ALW1, BLW1,  QSAT2, DQS2, ALW2, BLW2,        &
                     QSAT4, DQS4, ALW4, BLW4,  QSATS, DQSS, ALWS, BLWS,        &
                     BF1, BF2, BF3,VGWMAX,                                     &
                     CDCR1,CDCR2, psis, bee, poros, wpwet, cond, gnu,          &
                     ARS1,ARS2,ARS3,ARA1,ARA2,ARA3,ARA4,ARW1,ARW2,ARW3,ARW4,   &
                     tsa1,tsa2,tsb1,tsb2,atau,btau

      LOGICAL, INTENT(IN) :: BUG
      

! -----------------------------------------------------------
!     PROGNOSTIC VARIABLES

      REAL, INTENT(INOUT), DIMENSION(NCH) ::                                   &
                     TC1, TC2, TC4, QA1, QA2, QA4, CAPAC,                      &
                     CATDEF, RZEXC, SRFEXC
 
      REAL, INTENT(INOUT), DIMENSION(N_GT, NCH) ::  GHTCNT
 
      REAL, INTENT(INOUT), DIMENSION(N_SNOW, NCH) :: WESNN, HTSNNN, SNDZN
 

! -----------------------------------------------------------
!     DIAGNOSTIC OUTPUT VARIABLES

      REAL, INTENT(OUT), DIMENSION(NCH) ::      EVAP, SHFLUX, RUNOFF,          &
                     EINT, ESOI, EVEG, ESNO, BFLOW,RUNSRF,SMELT,               &
                     HLWUP,SWLAND,HLATN,QINFIL,AR1, AR2, RZEQ,                 &
                     GHFLUX, TPSN1, ASNOW0, TP1, TP2, TP3, TP4, TP5, TP6,      &
                     sfmc, rzmc, prmc, entot, wtot, tsurf, WCHANGE, ECHANGE,   &
                     HSNACC, EVACC, SHACC

      REAL, INTENT(OUT), DIMENSION(NCH), OPTIONAL :: SH_SNOW, AVET_SNOW,       &
                     WAT_10CM, TOTWAT_SOIL, TOTICE_SOIL

      
      REAL, INTENT(OUT), DIMENSION(NCH), OPTIONAL :: LHACC

      REAL, INTENT(OUT), DIMENSION(NCH), OPTIONAL :: TC1_0,TC2_0,TC4_0
      REAL, INTENT(OUT), DIMENSION(NCH), OPTIONAL :: QA1_0,QA2_0,QA4_0 
      REAL, INTENT(OUT), DIMENSION(NCH,N_SNOW), OPTIONAL :: fices_out
! -----------------------------------------------------------
!     LOCAL VARIABLES

      INTEGER I,K,N,LAYER

      REAL, DIMENSION(NCH) ::  CSOIL, ASNOW, traincx, trainlx, EMAXRT,         &
            RC, SATCAP, SNWFRC, POTFRC,  ESNFRC, EVSNOW, SHFLUXS, HLWUPS,      &
            HFTDS1, HFTDS2, HFTDS4, DHFT1, DHFT2, DHFT4, TPSN2, TPSN3,         &
            QSATTC, DQSDTC, SWSRF1, SWSRF2, SWSRF4, AR4, RX11, RX21, RX12,     &
            RX14, RX24, RX22, EIRFRC, FCAN, THRU, RZEQOL, frice, srfmx,        &
            srfmn, RCST, EVAPFR, RCUN, PAR, PDIR, RDCX, EVAP1, EVAP2,          &
            EVAP4, SHFLUX1, SHFLUX2, SHFLUX4, HLWUP1, HLWUP2, HLWUP4,          &
            GHFLUX1, GHFLUX2, GHFLUX4, RZI, TC1SF, TC2SF, TC4SF, ar1old,       &
            ar2old, ar4old, GHFLUXS, DEDQA1X, DEDTC1X,                         &
            DHSDQA1X, DHSDTC1X, DEDQA2X, DEDTC2X, DHSDQA2X, DHSDTC2X,          &
            DEDQA4X, DEDTC4X, DHSDQA4X, DHSDTC4X, werror, sfmcun, rzmcun,      &
            prmcun,WTOT_ORIG,ENTOT_ORIG,HSNACC1,HSNACC2,HSNACC4,               &
              qa1_orig,qa2_orig,qa4_orig,tc1_orig,tc2_orig,tc4_orig,           &
              tcs_orig


      REAL, DIMENSION(N_GT) :: HT, TP, soilice

      REAL, DIMENSION(N_SNOW) :: TPSN, WESN, HTSNN, SNDZ, fices

      REAL, DIMENSION(N_SM) :: T1, AREA, tkgnd, fhgnd

      REAL :: EPFRC1, EPFRC2, EPFRC4, SUMEP, SUME, TC1SN, TC2SN, TC4SN,        &
              DTC1SN,DTC2SN,DTC4SN, RTBS1, RTBS2, RTBS4, ZBAR, THETAF,         &
              XFICE, FH21, FH21W, FH21I, FH21D, DFH21W, DFH21I, DFH21D,        &
              EVSN, SHFLS, HUPS, HCORR, SWNET0, HLWDWN0, TMPSNW, HLWTC,        &
              DHLWTC, HSTURB, DHSDEA, DHSDTC, ESATTC, ETURB, DEDEA, DEDTC,     &
              SNOWF, TS, fh31w, fh31i, fh31d, pr, ea, desdtc, areasc,          &
              pre, dummy1, dummy2, dummy3, areasc0, EDIF, EINTX,               &
              SCLAI, tsn1, tsn2, tsn3, hold, hnew, emaxrz, dedtc0,             &
              dhsdtc0, alhfsn, ADJ, raddn, zc1, tsnowsrf, dum, tsoil,          &
              QA1X, QA2X, QA4X, TC1X, TC2X, TC4X, TCSX,                        &
              EVAPX1,EVAPX2,EVAPX4,SHFLUXX1,SHFLUXX2,SHFLUXX4,EVEGFRC,         &
              EVAPXS,SHFLUXXS,phi


      LOGICAL :: ldum

      integer  n_out
!      integer  n_outs(7)
!      data     n_outs / 30380 , 64674 , 64676 , 64678 , 64681 , 64795 , 64802 /
      integer  n_outs(1)
      data     n_outs / 828/

      !ag (02Jan2021)
      real, intent(inout) :: rivsto(nch)               ! river storage
      real, intent(inout) :: fldsto(nch)               ! flood storage
      real, intent(inout) :: fldfrc(nch)               ! flooded fraction

!      do n=1,size(n_outs)
      do n=1,0
         n_out =  n_outs(n)
         write (*,*) 'INPUT catchment arguments for n_out: ',n_out

         write (*,*) 'NCH = ',NCH
         write (*,*) 'DTSTEP = ',DTSTEP
         write (*,*) 'SFRAC = ',SFRAC
         write (*,*) 'ITYP = ',ITYP(n_out)
         write (*,*) 'TRAINC = ',TRAINC(n_out)
         write (*,*) 'TRAINL = ',TRAINL(n_out)
         write (*,*) 'TSNOW = ',TSNOW(n_out)
         write (*,*) 'UM = ',UM(n_out)
         write (*,*) 'ETURB1 = ',ETURB1(n_out)
         write (*,*) 'DEDQA1 = ',DEDQA1(n_out)
         write (*,*) 'DEDTC1 = ',DEDTC1(n_out)
         write (*,*) 'HSTURB1 = ',HSTURB1(n_out)
         write (*,*) 'DHSDQA1 = ',DHSDQA1(n_out)
         write (*,*) 'DHSDTC1 = ',DHSDTC1(n_out)
         write (*,*) 'ETURB2 = ',ETURB2(n_out)
         write (*,*) 'DEDQA2 = ',DEDQA2(n_out)
         write (*,*) 'DECTC2 = ',DEDTC2(n_out)
         write (*,*) 'HSTURB2 = ',HSTURB2(n_out)
         write (*,*) 'DHSDQA2 = ',DHSDQA2(n_out)
         write (*,*) 'DHSDTC2 = ',DHSDTC2(n_out)
         write (*,*) 'ETURB4 = ',ETURB4(n_out)
         write (*,*) 'DEDQA4 = ',DEDQA4(n_out)
         write (*,*) 'DEDTC4 = ',DEDTC4(n_out)
         write (*,*) 'HSTURB4 = ',HSTURB4(n_out)
         write (*,*) 'HDSDQA4 = ',DHSDQA4(n_out)
         write (*,*) 'DHSDTC4 = ',DHSDTC4(n_out)
         write (*,*) 'ETURBS = ',ETURBS(n_out)
         write (*,*) 'DEDQAS = ',DEDQAS(n_out)
         write (*,*) 'DEDTCS = ',DEDTCS(n_out)
         write (*,*) 'HSTURBS = ',HSTURBS(n_out)
         write (*,*) 'DHSDQAS = ',DHSDQAS(n_out)
         write (*,*) 'DHSDTCS = ',DHSDTCS(n_out)
         write (*,*) 'TM = ',TM(n_out)
         write (*,*) 'QM = ',QM(n_out)
         write (*,*) 'ra1 = ',ra1(n_out)
         write (*,*) 'ra2 = ',ra2(n_out)
         write (*,*) 'ra4 = ',ra4(n_out)
         write (*,*) 'raS = ',raS(n_out)
         write (*,*) 'SUNANG = ',SUNANG(n_out)
         write (*,*) 'PARDIR = ',PARDIR(n_out)
         write (*,*) 'PARDIF = ',PARDIF(n_out)
         write (*,*) 'SWNETF = ',SWNETF(n_out)
         write (*,*) 'SWNETS = ',SWNETS(n_out)
         write (*,*) 'HLWDWN = ',HLWDWN(n_out)
         write (*,*) 'PSUR = ',PSUR(n_out)
         write (*,*) 'ZLAI = ',ZLAI(n_out)
         write (*,*) 'GREEN = ',GREEN(n_out)
         write (*,*) 'Z2 = ',Z2(n_out)
         write (*,*) 'SQSCAT = ',SQSCAT(n_out)
         write (*,*) 'RSOIL1 = ',RSOIL1(n_out)
         write (*,*) 'RSOIL2 = ',RSOIL2(n_out)
         write (*,*) 'RDC = ',RDC(n_out)
!         write (*,*) 'USFAC = ',U2FAC(n_out)
         write (*,*) 'QSAT1 = ',QSAT1(n_out)
         write (*,*) 'DQS1 = ',DQS1(n_out)
         write (*,*) 'ALW1 = ',ALW1(n_out)
         write (*,*) 'BLW1 = ',BLW1(n_out)
         write (*,*) 'QSAT2 = ',QSAT2(n_out)
         write (*,*) 'DQS2 = ',DQS2(n_out)
         write (*,*) 'ALW2 = ',ALW2(n_out)
         write (*,*) 'BLW2 = ',BLW2(n_out)
         write (*,*) 'QSAT4 = ',QSAT4(n_out)
         write (*,*) 'DQS4 = ',DQS4(n_out)
         write (*,*) 'ALW4 = ',ALW4(n_out)
         write (*,*) 'BLW4 = ',BLW4(n_out)
         write (*,*) 'QSATS = ',QSATS(n_out)
         write (*,*) 'DQSS = ',DQSS(n_out)
         write (*,*) 'ALWS = ',ALWS(n_out)
         write (*,*) 'BLWS = ',BLWS(n_out)
         write (*,*) 'BF1 = ',BF1(n_out)
         write (*,*) 'BF2 = ',BF2(n_out)
         write (*,*) 'BF3 = ',BF3(n_out)
         write (*,*) 'VGWMAX = ',VGWMAX(n_out)
         write (*,*) 'CDCR1 = ',CDCR1(n_out)
         write (*,*) 'CDCR2 = ',CDCR2(n_out)
         write (*,*) 'psis = ',psis(n_out)
         write (*,*) 'bee = ',bee(n_out)
         write (*,*) 'poros = ',poros(n_out)
         write (*,*) 'wpwet = ',wpwet(n_out)
         write (*,*) 'cond = ',cond(n_out)
         write (*,*) 'gnu = ',gnu(n_out)
         write (*,*) 'ars1 = ',ARS1(n_out)
         write (*,*) 'ars2 = ',ARS2(n_out)
         write (*,*) 'ars3 = ',ARS3(n_out)
         write (*,*) 'ara1 = ',ARA1(n_out)
         write (*,*) 'ara2 = ',ARA2(n_out)
         write (*,*) 'ara3 = ',ARA3(n_out)
         write (*,*) 'ara4 = ',ARA4(n_out)
         write (*,*) 'arw1 = ',ARW1(n_out)
         write (*,*) 'arw2 = ',ARW2(n_out)
         write (*,*) 'arw3 = ',ARW3(n_out)
         write (*,*) 'arw4 = ',ARW4(n_out)
         write (*,*) 'tsa1 = ',tsa1(n_out)
         write (*,*) 'tsa2 = ',tsa2(n_out)
         write (*,*) 'tsb1 = ',tsb1(n_out)
         write (*,*) 'tsb2 = ',tsb2(n_out)
         write (*,*) 'atau = ',atau(n_out)
         write (*,*) 'btau = ',btau(n_out)
         write (*,*) 'BUG = ',BUG
         write (*,*) 'TC1 = ',TC1(n_out)
         write (*,*) 'TC2 = ',TC2(n_out)
         write (*,*) 'TC4 = ',TC4(n_out)
         write (*,*) 'QA1 = ',QA1(n_out)
         write (*,*) 'QA2 = ',QA2(n_out)
         write (*,*) 'QA4 = ',QA4(n_out)
         write (*,*) 'CAPAC = ',CAPAC(n_out)
         write (*,*) 'CATDEF = ',CATDEF(n_out)
         write (*,*) 'RZEXC = ',RZEXC(n_out)
         write (*,*) 'SRFEXC = ',SRFEXC(n_out)
         write (*,*) 'GHTCNT = ',GHTCNT(:,n_out)
         write (*,*) 'TSURF = ',TSURF(n_out)
         write (*,*) 'WESNN = ',WESNN(:,n_out)
         write (*,*) 'HTSNNN = ',HTSNNN(:,n_out)
         write (*,*) 'SNDZN = ',SNDZN(:,n_out)
         write (*,*) 'EVAP = ',EVAP(n_out)
         write (*,*) 'SHFLUX = ',SHFLUX(n_out)
         write (*,*) 'RUNOFF = ',RUNOFF(n_out)
         write (*,*) 'EINT = ',EINT(n_out)
         write (*,*) 'ESOI = ',ESOI(n_out)
         write (*,*) 'EVEG = ',EVEG(n_out)
         write (*,*) 'ESNO = ',ESNO(n_out)
         write (*,*) 'BFLOW = ',BFLOW(n_out)
         write (*,*) 'RUNSRF = ',RUNSRF(n_out)
         write (*,*) 'SMELT = ',SMELT(n_out)
         write (*,*) 'HLWUP = ',HLWUP(n_out)
         write (*,*) 'HLATN = ',HLATN(n_out)
         write (*,*) 'QINFIL = ',QINFIL(n_out)
         write (*,*) 'AR1 = ',AR1(n_out)
         write (*,*) 'AR2 = ',AR2(n_out)
         write (*,*) 'RZEQ = ',RZEQ(n_out)
         write (*,*) 'GHFLUX = ',GHFLUX(n_out)
         write (*,*) 'TPSN1 = ',TPSN1(n_out)
         write (*,*) 'ASNOW0 = ',ASNOW0(n_out)
         write (*,*) 'TP1 = ',TP1(n_out)
         write (*,*) 'TP2 = ',TP2(n_out)
         write (*,*) 'TP3 = ',TP3(n_out)
         write (*,*) 'TP4 = ',TP4(n_out)
         write (*,*) 'TP5 = ',TP5(n_out)
         write (*,*) 'TP6 = ',TP6(n_out)
      enddo 
!rr ------------------------------------------------------------------      


!rr      DATA RCST /mxnch*1.E10/

      do n=1,nch
         rcst(n) = 1.E10
      end do
      

!**** ---------------------------------------------------
!**** PRE-PROCESS DATA AS NECESSARY:
!****

!rr      if (nch .gt. mxnch) then
!rr        write(*,*) 'catchment.f mxnch exceed: must be greater than',nch
!rr        stop
!rr        end if

      DO N=1,NCH
!       SATCAP(N) = 0.1 * ZLAI(N)
! change for convergence towards MOSAIC
!        SATCAP(N) = 0.2 * ZLAI(N) + 1.e-5
        SATCAP(N) = SATCAPFR * ZLAI(N) + 1.e-5
!!AMM        SATCAP(N) = 1.0 * ZLAI(N) + 1.e-5
        CSOIL(N)  = CSOIL_1
        if ( ityp(n) .ne. 1) CSOIL(N)  = CSOIL_2
        FCAN(N) = AMIN1( 1., AMAX1(0.,CAPAC(N)/SATCAP(N)) )
        SCLAI=amin1( 1., zlai(n)/2. )
        POTFRC(N)=FCAN(N)*SCLAI
        if(fcan(n) .lt. .1) POTFRC(N)=POTFRC(N)*(10.*fcan(n))

!       if(LIS_localPet.eq.165.and.n.eq.22527) then 
!          write(LIS_logunit,*) 'potf ',capac(n), satcap(n), &
!               fcan(n),sclai, potfrc(n)
!       endif
! Correction to RDC formulation -Randy Koster, 4/1/2011
!        RDCX(N)    = RDC(N)*SCLAI
        RDCX(N)    = RDC(N)

        DEDQA1X(N)  = AMAX1( DEDQA1(N), 500./ALHE )
        DEDTC1X(N)  = AMAX1( DEDTC1(N),   0. )
        DHSDQA1X(N) = AMAX1( DHSDQA1(N),   0. )
        DHSDTC1X(N) = AMAX1( DHSDTC1(N), -10. )

        DEDQA2X(N)  = AMAX1( DEDQA2(N), 500./ALHE )
        DEDTC2X(N)  = AMAX1( DEDTC2(N),   0. )
        DHSDQA2X(N) = AMAX1( DHSDQA2(N),   0. )
        DHSDTC2X(N) = AMAX1( DHSDTC2(N), -10. )

        DEDQA4X(N)  = AMAX1( DEDQA4(N), 500./ALHE )
        DEDTC4X(N)  = AMAX1( DEDTC4(N),   0. )
        DHSDQA4X(N) = AMAX1( DHSDQA4(N),   0. )
        DHSDTC4X(N) = AMAX1( DHSDTC4(N), -10. )


        qa1_orig(n)=qa1(n)
        qa2_orig(n)=qa2(n)
        qa4_orig(n)=qa4(n)
        tc1_orig(n)=tc1(n)
        tc2_orig(n)=tc2(n)
        tc4_orig(n)=tc4(n)


        if(ityp(n) .ge. 7) potfrc(n)=0.
!$$$        RA1(N)     = ONE / ( CD1(N) * max(UM(N),1.) )
!$$$        RA2(N)     = ONE / ( CD2(N) * max(UM(N),1.) )
!$$$        RA4(N)     = ONE / ( CD4(N) * max(UM(N),1.) )
!$$$        RAS(N)     = ONE / ( CDS(N) * max(UM(N),1.) ) 


!     HSNACC is an energy accounting term designed to account (among other,
!     lesser things) for the fact that snow is deposited at the snowpack
!     surface temperature while the atmosphere does not account for variations
!     in the heat content of deposited snow.
 
        HSNACC(N)=0.
        EVACC(N)=0.
        SHACC(N)=0.
        RUNSRF(N)=0.


!****   RESET LAND ICE VARIABLES, MAINTAINING TEMPS. AT EACH LAYER
        IF(ITYP(N) .EQ. 9) THEN
          if(sum(htsnnn(:,n)+wesnn(:,n))==0.) then
              TSN1=tc1(n)-273.16
              TSN2=tc1(n)-273.16
              TSN3=tc1(n)-273.16
            else
              TSN1=(HTSNNN(1,N)+WESNN(1,N)*ALHM)/(SCONST*WESNN(1,N)+1.e-5)
              TSN2=(HTSNNN(2,N)+WESNN(2,N)*ALHM)/(SCONST*WESNN(2,N)+1.e-5)
              TSN3=(HTSNNN(3,N)+WESNN(3,N)*ALHM)/(SCONST*WESNN(3,N)+1.e-5)
           endif
          WESNN(1,N)=.1
          WESNN(2,N)=.2
          WESNN(3,N)=.1
          HTSNNN(1,N)=-ALHM*WESNN(1,N)+TSN1*SCONST*WESNN(1,N)
          HTSNNN(2,N)=-ALHM*WESNN(2,N)+TSN1*SCONST*WESNN(2,N)
          HTSNNN(3,N)=-ALHM*WESNN(3,N)+TSN1*SCONST*WESNN(3,N)
          SNDZN(1,N)=WESNN(1,N)/.9
          SNDZN(2,N)=WESNN(2,N)/.9
          SNDZN(3,N)=WESNN(3,N)/.9
          POTFRC(N)=1.

       ENDIF

!****   RESET LAKE VARIABLES
        IF(ITYP(N) .EQ. 10) THEN
          CATDEF(N)=0.
          RZEXC(N)=0.
          SRFEXC(N)=0.
          SATCAP(N)=1000.
          CAPAC(N)=SATCAP(N)
          POTFRC(N)=1.
          ENDIF

        ENDDO

!**** ---------------------------------------------------
!**** DETERMINE INITIAL VALUE OF RZEQ:

      CALL RZEQUIL (                                                           &
                    NCH, ITYP, CATDEF, VGWMAX,CDCR1,CDCR2,WPWET,               &
                    ars1,ars2,ars3,ara1,ara2,ara3,ara4,                        &
                    arw1,arw2,arw3,arw4,                                       &
                    RZEQOL                                                     &
                   )

      IF (BUG) THEN
        WRITE(*,*) 'RZEQUIL OK'
        ENDIF


!rr   switched order of call to partition and do-loop for emaxrz & emaxrt
!rr   because srfmn was not initialized in the computation of emaxrt
!rr   reichle, Oct 22, 2003

!**** PARTITION CATCHMENT INTO EVAPORATION SUBREGIONS:
!****
      CALL PARTITION (                                                         &
                      NCH,DTSTEP,ITYP,DZSF,RZEXC,  RZEQOL,VGWMAX,CDCR1,CDCR2,  &
                      PSIS,BEE,poros,WPWET,                                    &
                      ars1,ars2,ars3,ara1,ara2,ara3,ara4,                      &
                      arw1,arw2,arw3,arw4,BUG,                                 &
                      SRFEXC,CATDEF,RUNSRF,                                    &
                      AR1, AR2, AR4,srfmx,srfmn,  SWSRF1,SWSRF2,SWSRF4,RZI     &
                     )


      DO N=1,NCH
         TSOIL=AR1(N)*TC1(N)+AR2(N)*TC2(N)+AR4(N)*TC4(N)
         WTOT_ORIG(N) = WESNN(1,N) + WESNN(2,N) + WESNN(3,N) +                 &
               (CDCR2(N)-CATDEF(N)) + RZEXC(N) +  SRFEXC(N) + CAPAC(N)
         ENTOT_ORIG(N) = HTSNNN(1,N) + HTSNNN(2,N) + HTSNNN(3,N) +             &
            TSOIL*CSOIL(N) + GHTCNT(1,N) + GHTCNT(2,N) +                       &
            GHTCNT(3,N) + GHTCNT(4,N) + GHTCNT(5,N) + GHTCNT(6,N)
         ENDDO


      DO N=1,NCH
         emaxrz=amax1(0.,RZEQOL(N)+RZEXC(N)-WPWET(N)*VGWMAX(N))
         EMAXRT(N)=(CAPAC(N)+emaxrz+(SRFEXC(N)-SRFMN(N)))/DTSTEP
         ENDDO
      
      do n=1,nch
        ar1old(n)=ar1(n)
        ar2old(n)=ar2(n)
        ar4old(n)=ar4(n)
        enddo

      IF (BUG) THEN
        WRITE(*,*) 'PARTITION OK'
        ENDIF

!**** ========================================================
!**** ENERGY BALANCES.

!**** COMPUTE "INITIAL ESTIMATE" OF HEAT FLUX TO DEEP SOIL (HFTDS)
!**** AND ITS DERIVATIVE WITH RESPECT TO TEMPERATURE (DHFTDS):

      DO N=1,NCH
        T1(1)=TC1(N)
        T1(2)=TC2(N)
        T1(3)=TC4(N)
        if (PHIGT<0.) then ! if statement for bkwd compatibility w/ off-line MERRA replay
           phi=POROS(N)
        else
           phi=PHIGT
        end if
        ZBAR=-SQRT(1.e-20+catdef(n)/bf1(n))-bf2(n)
        THETAF=.5
        DO LAYER=1,6
          HT(LAYER)=GHTCNT(LAYER,N)
          ENDDO

        CALL GNDTP0(                                                           &
                    T1,phi,ZBAR,THETAF,                                        &
                    HT,                                                        &
                    fh21w,fH21i,fh21d,dfh21w,dfh21i,dfh21D,tp                  &
                   )

        HFTDS1(N)=-FH21W
        HFTDS2(N)=-FH21I
        HFTDS4(N)=-FH21D
        DHFT1(N)=-DFH21W
        DHFT2(N)=-DFH21I
        DHFT4(N)=-DFH21D

        ENDDO 

      IF (BUG) THEN
        WRITE(*,*) 'HEAT FLUX INITIAL ESTIMATE OK'
        ENDIF

!**** -------------------------------------------------------------
!**** A. SNOW-FREE FRACTION.
!**** DETERMINE EVAPORATION, SENSIBLE HEAT FLUXES; UPDATE TEMPS:

      DO N=1,NCH
        PAR(N)    = PARDIR(N) + PARDIF(N) + 1.E-20
        PDIR(N)   = PARDIR(N) / PAR(N)
        TC1SF(N)  = TC1(N)
        TC2SF(N)  = TC2(N)
        TC4SF(N)  = TC4(N)
        ENDDO

      CALL RCUNST (                                                            &
                   NCH, ITYP, SUNANG, SQSCAT, PDIR, PAR, ZLAI, GREEN, BUG,     &
                   RCUN                                                        &
                  )

      IF (BUG) THEN
         WRITE(*,*) 'RCUNST OK'
      ENDIF

!**** 1. SATURATED FRACTION

      CALL ENERGY1 (                                                           &
                   NCH, DTSTEP, ITYP, UM, RCUN,                                &
                   ETURB1, DEDQA1X, DEDTC1X, HSTURB1, DHSDQA1X, DHSDTC1X,      &
                   QM,     RA1,   SWNETF,  HLWDWN, PSUR,                       &
                   RDCX,    HFTDS1, DHFT1,  QSAT1, DQS1, ALW1, BLW1,           &
                   EMAXRT,CSOIL,SWSRF1,POTFRC,.false.,                         &
                   TC1SF, QA1,                                                 &
                   EVAP1, SHFLUX1, HLWUP1, RX11, RX21, GHFLUX1, HSNACC1)

      IF (BUG) THEN
        WRITE(*,*) 'ENERGY1 OK'
        ENDIF

!**** 2. SUBSATURATED BUT UNSTRESSED FRACTION

!       if(LIS_localPet.eq.165) then 
!          write(LIS_logunit,*) 'potf1 ',potfrc(22527)
!       endif

!CC    print*,'energy2'
      CALL ENERGY2 (                                                           &
                   NCH, DTSTEP, ITYP, UM, RCUN,                                &
                   ETURB2, DEDQA2X, DEDTC2X, HSTURB2, DHSDQA2X, DHSDTC2X,      &
                   QM,     RA2,   SWNETF,  HLWDWN, PSUR,                       &
                   RDCX,    HFTDS2, DHFT2, QSAT2, DQS2, ALW2, BLW2,            &
                   EMAXRT,CSOIL,SWSRF2,POTFRC,.false., RZI, WPWET,             &
                   TC2SF, QA2,                                                 &
                   EVAP2, SHFLUX2, HLWUP2, RX12, RX22, GHFLUX2, HSNACC2)

      IF (BUG) THEN
        WRITE(*,*) 'ENERGY2 OK'
        ENDIF

!**** 3. WILTING FRACTION
!CC    print*,'energy4'
      CALL ENERGY4 (                                                           &
                   NCH, DTSTEP, ITYP, UM, RCST,                                &
                   ETURB4, DEDQA4X, DEDTC4X, HSTURB4, DHSDQA4X, DHSDTC4X,      &
                   QM,     RA4,   SWNETF,  HLWDWN, PSUR,                       &
                   RDCX,   HFTDS4, DHFT4, QSAT4, DQS4, ALW4, BLW4,             &
                   EMAXRT,CSOIL,SWSRF4,POTFRC,.false., WPWET,                  &
                   TC4SF, QA4,                                                 &
                   EVAP4, SHFLUX4, HLWUP4, RX14, RX24, GHFLUX4, HSNACC4)

      IF (BUG) THEN
        WRITE(*,*) 'ENERGY4 OK'
        ENDIF

!**** COMPUTE EIRFRC
      DO N=1,NCH
         
        RTBS1=RX11(N)*RX21(N)/(RX11(N)+RX21(N)+1.E-20)
        EPFRC1=POTFRC(N) * ( RA1(N) + RTBS1 ) / ( RA1(N) + POTFRC(N)*RTBS1 )
         
        RTBS2=RX12(N)*RX22(N)/(RX12(N)+RX22(N)+1.E-20)
        EPFRC2=POTFRC(N) * ( RA2(N) + RTBS2 ) / ( RA2(N) + POTFRC(N)*RTBS2 )
         
        RTBS4=RX14(N)*RX24(N)/(RX14(N)+RX24(N)+1.E-20)
        EPFRC4=POTFRC(N) * ( RA4(N) + RTBS4 ) / ( RA4(N) + POTFRC(N)*RTBS4 )
         
        SUMEP=EPFRC1*EVAP1(N)*AR1(N)+EPFRC2*EVAP2(N)*AR2(N)+                   &
              EPFRC4*EVAP4(N)*AR4(N)
        SUME=EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N)
        
        !   "quick fix" gone wrong in global AMSR-E assimilation
        !   trying to correct while staying as close as possible to past fix
        !   30 July 2007, reichle: 
        !   

!           EIRFRC(N)=SUMEP/(SUME+1.E-20)

        if (SUME/=-1.e-20) then
            EIRFRC(N)=SUMEP/(SUME+1.E-20)
          else
            EIRFRC(N)=SUMEP/(SUME+2.E-20)
        end if

        ENDDO

      IF (BUG) THEN
        WRITE(*,*) 'EIRFRC OK'
        ENDIF

!**** --------------------------------------------------------
!**** B. SNOW-COVERED FRACTION.

      DO N=1,NCH

        TS     = TM(N) 
        T1(1)  = TC1(N)-273.16 
        T1(2)  = TC2(N)-273.16 
        T1(3)  = TC4(N)-273.16 
        AREA(1)= AR1(N) 
        AREA(2)= AR2(N) 
        AREA(3)= AR4(N) 
        pr     = trainc(n)+trainl(n)+tsnow(n) 
        snowf  = tsnow(n) 
        dedea  = dedqas(n)*epsilon/psur(n) 
        dhsdea = dhsdqas(n)*epsilon/psur(n) 
        ea     = qm(n)*psur(n)/epsilon 
        esattc = qsats(n)*psur(n)/epsilon 
        desdtc = dqss(n)*psur(n)/epsilon 
        dedtc0  = dedtcs(n) + dedea*desdtc 
        dhsdtc0 = dhsdtcs(n) + dhsdea*desdtc 
        hsturb=hsturbs(n) 
        tkgnd(1)=1.8      !STEPH  
        tkgnd(2)=1.8 
        tkgnd(3)=1.8 
        raddn=hlwdwn(n)+swnets(n) 
        zc1=-.025 
        hups=0.0 
 
!**** 1. RUN SNOW MODEL: 
 
        do i=1,3 
          wesn(i)=wesnn(i,n)
          htsnn(i)=htsnnn(i,n) 
          sndz(i)=sndzn(i,n) 
          tpsn(i)=0.0 
          enddo 

!     TPSN1 is used as input here, contradicts "declaration" as output only.
!     EnKF has been using tpsn1 as part of state vector, which should fix this.
!     reichle, 18 Nov 02

!     Removed tpsn1 from state vector, now compute tpsn1 from prognostics
!     in process
!     reichle, 29 May 03

        call get_tf0d(htsnn(1),wesn(1),tsnowsrf,dum,ldum,ldum)
        tcs_orig(n)=tsnowsrf+tf
        if(wesn(1)+wesn(2)+wesn(3) .eq. 0.) tcs_orig(n)=                       &
                  amin1( tf, tc1_orig(n)*ar1(n)+tc2_orig(n)*ar2(n)+            &
                  tc4_orig(n)*(1.-ar1(n)-ar2(n)) )

        hlwtc=ALWS(N) + BLWS(N)*(TSNOWSRF+TF) 
        dhlwtc=BLWS(N)
        hcorr=0.

        CALL SNOWRT(n,                                                         &
                   t1,area,tkgnd,pr,snowf,ts,DTSTEP,                           &
                   eturbs(n),dedtc0,hsturb,dhsdtc0,hlwtc,dhlwtc,               &
                   desdtc,hups,raddn,zc1,                                      &
                   wesn,htsnn,sndz,                                            &
                   fices,tpsn,areasc,areasc0,pre,fhgnd,                        &
                   EVSN,SHFLS,alhfsn,hcorr                                     &
                  ) 
!        print *, n
!        if (N .eq.) then 
!           print *, 'Checking the parameters to SNOWRT at n=',n 
!           print *, 't1=', t1
!           print *, 'area=', area
!           print *, 'tkgnd=', tkgnd
!           print *, 'pr=', pr
!           print *, 'snowf', snowf
!           print *, 'ts=', ts
!           print *, 'DTSTEP=', DTSTEP
!           print *, 'etubrbs(n)=', eturbs(n)
!           print *, 'dedtc0=', dedtc0
!           print *, 'hsturb=', hsturb
!           print *, 'dhsdtc0=', dhsdtc0
!           print *, 'hlwtc=', hlwtc
!           print *, 'dhlwtc=', dhlwtc
!           print *, 'desdtc=', desdtc
!           print *, 'hups=', hups
!           print *, 'raddn=', raddn
!           print *, 'zc1=', zc1
!           print *, 'wesn=', wesn
!           print *, 'htsnn=', htsnn
!           print *, 'sndz=', sndz
!           print *, 'fices=', fices
!           print *, 'tpsn=', tpsn
!           print *, 'areasc=', areasc
!           print *, 'areasc0=', areasc0
!           print *, 'pre=', pre
!           print *, 'fhgnd=', fhgnd
!           print *, 'EVSN=', EVSN
!           print *, 'SHFLS=', SHFLS
!           print *, 'alhfsn=', alhfsn
!           print *, 'hcorr=', hcorr
!        endif

 
        TPSN1(N) = TPSN(1)+TF 
        TPSN2(N) = TPSN(2)+TF 
        TPSN3(N) = TPSN(3)+TF 

        SMELT(N) = PRE            
        fh31w=fhgnd(1) 
        fh31i=fhgnd(2) 
        fh31d=fhgnd(3) 
        asnow(n) = areasc 
        asnow0(n)= areasc0 
        HSNACC(N) = HSNACC(N) + (1.-ASNOW(N))*                                 &
             (HSNACC1(N)*AR1(N)+HSNACC2(N)*AR2(N)+HSNACC4(N)*AR4(N))           &
             + hcorr
 
!STEPH 
        do i=1,3 
          wesnn(i,n)=wesn(i) 
          htsnnn(i,n)=htsnn(i) 
          sndzn(i,n)=sndz(i) 
          enddo 
 
        traincx(n)= trainc(n)*(1.-areasc) 
        trainlx(n)= trainl(n)*(1.-areasc)

!**** 2. UPDATE SURFACE TEMPERATURE

        DTC1SN=((-(FH31W/(area(1)+1.e-20))-HFTDS1(N))*DTSTEP/CSOIL(N))/        &
                  (1.+DHFT1(N)*DTSTEP/CSOIL(N))
        DTC2SN=((-(FH31I/(area(2)+1.e-20))-HFTDS2(N))*DTSTEP/CSOIL(N))/        &
                  (1.+DHFT2(N)*DTSTEP/CSOIL(N))
        DTC4SN=((-(FH31D/(area(3)+1.e-20))-HFTDS4(N))*DTSTEP/CSOIL(N))/        &
                  (1.+DHFT4(N)*DTSTEP/CSOIL(N))
!        TC1SN=TC1(N)+DTC1SN
!        TC2SN=TC2(N)+DTC2SN
!        TC4SN=TC4(N)+DTC4SN

        TC1SN=TC1(N)+DTC1SN
        IF((TC1SN-TPSN3(N))*(TC1(N)-TPSN3(N)) .LT. 0.) THEN
          HSNACC(N)=HSNACC(N)+AREASC*AREA(1)*                                 &
                 (TC1SN-TPSN3(N))*CSOIL(N)/DTSTEP
          TC1SN=TPSN3(N)
          ENDIF

        TC2SN=TC2(N)+DTC2SN
        IF((TC2SN-TPSN3(N))*(TC2(N)-TPSN3(N)) .LT. 0.) THEN
          HSNACC(N)=HSNACC(N)+AREASC*AREA(2)*                                 &
                 (TC2SN-TPSN3(N))*CSOIL(N)/DTSTEP
          TC2SN=TPSN3(N)
          ENDIF

        TC4SN=TC4(N)+DTC4SN
        IF((TC4SN-TPSN3(N))*(TC4(N)-TPSN3(N)) .LT. 0.) THEN
          HSNACC(N)=HSNACC(N)+AREASC*AREA(3)*                                 &
                 (TC4SN-TPSN3(N))*CSOIL(N)/DTSTEP
          TC4SN=TPSN3(N)
          ENDIF


!          if(LIS_localPet.eq.129.and.n.eq.436) then 
!             write(LIS_logunit,*) 'mod ',areasc,TC4SF(N),TC4SN
!          endif

!rr        if ( TC1SN<240. .or. TC1SN>350.          ) then
!rr           write (*,*) 'bullshit2: ', N, TC1SN
!rr        end if
        
        TC1(N)=TC1SF(N)*(1-AREASC)+TC1SN*AREASC
        TC2(N)=TC2SF(N)*(1-AREASC)+TC2SN*AREASC
        TC4(N)=TC4SF(N)*(1-AREASC)+TC4SN*AREASC
        
!        if(LIS_localPet.eq.165.and.n.eq.22527) then 
!           write(LIS_logunit,*) 'tptr_a ',tc2(n),areasc, tc2sf(n),tc2sn
!        endif
        EVSNOW(N)=EVSN
        esno(n)=evsnow(n)*asnow(n)*DTSTEP ! to have esno in mm/20min (03-17-99)
        SHFLUXS(N)=SHFLS 
        HLWUPS(N) =HUPS 
        GHFLUXS(N)=AREA(1)*(HFTDS1(N)+DHFT1(N)*DTC1SN) +                       &
                   AREA(2)*(HFTDS2(N)+DHFT2(N)*DTC2SN) +                       &
                   AREA(3)*(HFTDS4(N)+DHFT4(N)*DTC4SN)
        ENDDO 
 


      IF (BUG) THEN 
        WRITE(*,*) 'SNOW FRACTION OK' 
        ENDIF 
 
      DO N=1,NCH 
        HLATN(N)=(1.-ASNOW(N))*                                                &
              (EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N))*ALHE           &
              +ASNOW(N)*EVSNOW(N)*ALHS 
        EVAP(N)=(1.-ASNOW(N))*                                                 &
              (EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N))                &
              +ASNOW(N)*EVSNOW(N) 
        EVAPFR(N)=(1.-ASNOW(N))*                                               &
              (EVAP1(N)*AR1(N)+EVAP2(N)*AR2(N)+EVAP4(N)*AR4(N))
        SHFLUX(N)=(1.-ASNOW(N))*                                               &
              (SHFLUX1(N)*AR1(N)+SHFLUX2(N)*AR2(N)+SHFLUX4(N)*AR4(N))          &
              +ASNOW(N)*SHFLUXS(N) 
        HLWUP(N)=(1.-ASNOW(N))*                                                &
              (HLWUP1(N)*AR1(N)+HLWUP2(N)*AR2(N)+HLWUP4(N)*AR4(N))             &
              +ASNOW(N)*HLWUPS(N) 
        SWLAND(N)=(1.-ASNOW(N))*SWNETF(N) + ASNOW(N)*SWNETS(N) 
        GHFLUX(N)=(1.-ASNOW(N))*                                               &
              (GHFLUX1(N)*AR1(N)+GHFLUX2(N)*AR2(N)+GHFLUX4(N)*AR4(N))          &
              +ASNOW(N)*GHFLUXS(N) 
        ENDDO 


      IF (BUG) THEN
        WRITE(*,*) 'ENERGY FLUXES OK'
        ENDIF


!****
!**** NOW ALLOW DEEPER SOIL TEMPERATURES TO BE UPDATED:

      DO N=1,NCH
        if (PHIGT<0.) then ! if statement for bkwd compatibility w/ off-line MERRA replay
           phi=POROS(N)
        else
           phi=PHIGT
        end if
        ZBAR=-SQRT(1.e-20+catdef(n)/bf1(n))-bf2(n)
        THETAF=.5
        DO LAYER=1,6
          HT(LAYER)=GHTCNT(LAYER,N)
          ENDDO
        FH21=-GHFLUX(N)

        CALL GNDTMP(                                                           &
              dtstep,phi,zbar,thetaf,fh21,                                     &
              ht,                                                              &
              xfice,tp, soilice)

        DO LAYER=1,6
          GHTCNT(LAYER,N)=HT(LAYER)
          ENDDO
        tp1(n)=tp(1)
        tp2(n)=tp(2)
        tp3(n)=tp(3)
        tp4(n)=tp(4)
        tp5(n)=tp(5)
        tp6(n)=tp(6)
        frice(n)=xfice
         
        ENDDO


      IF (BUG) THEN
        WRITE(*,*) 'DEEPER SOIL TEMPERATURES UPDATE OK'
        ENDIF

!**** ========================================================

!**** REMOVE EVAPORATED WATER FROM SURFACE RESERVOIRS:
!****
!**** (FIRST CORRECT FOR EXCESSIVE INTERCEPTION LOSS)

      DO N=1,NCH
        EINTX=EIRFRC(N)*EVAPFR(N)*DTSTEP
        IF(EINTX .GT. CAPAC(N)) THEN
          EDIF=(EINTX-CAPAC(N))/DTSTEP
!          EVACC(N)=EVACC(N)-EDIF
          EVAPFR(N)=EVAPFR(N)-EDIF
          EVAP(N)=EVAP(N)-EDIF
          HLATN(N)=HLATN(N)-EDIF*ALHE
          SHFLUX(N)=SHFLUX(N)+EDIF*ALHE
!          SHACC(N)=SHACC(N)+EDIF*ALHE
!          HSNACC(N)=HSNACC(N)+EDIF*ALHE
          EIRFRC(N)=CAPAC(N)/((EVAPFR(N)+1.E-20)*DTSTEP)
          ENDIF
        ENDDO


      CALL WUPDAT (                                                            &
                     NCH,   ITYP, DTSTEP, EVAPFR, SATCAP, TC1, RA1, RC,        &
                     RX11,RX21,RX12,RX22,RX14,RX24,                            &
                     AR1,AR2,AR4,CDCR1,EIRFRC,RZEQOL,srfmn,WPWET,VGWMAX,       &
                     CAPAC, RZEXC, CATDEF, SRFEXC,                             &
                     EINT, ESOI, EVEG                                          &
                    )

! ---------------------------------------------------------------------

      IF (BUG) THEN
        WRITE(*,*) 'WUPDAT OK'
        ENDIF

!**** REDISTRIBUTE MOISTURE BETWEEN RESERVOIRS:

      CALL RZDRAIN (                                                           &
                    NCH,DTSTEP,VGWMAX,SATCAP,RZEQOL,AR1,WPWET,                 &
                    tsa1,tsa2,tsb1,tsb2,atau,btau,CDCR2,poros,BUG,             &
                    CAPAC,RZEXC,SRFEXC,CATDEF,RUNSRF                           &
                    )

! ---------------------------------------------------------------------

      IF (BUG) THEN
        WRITE(*,*) 'RZDRAIN OK'
        ENDIF

!**** COMPUTE BASEFLOW FROM TOPMODEL EQUATIONS

      CALL BASE (                                                              &
                 NCH, DTSTEP,BF1, BF2, BF3, CDCR1, FRICE, COND, GNU, CDCR2,    &
                 CATDEF,                                                       &
                 BFLOW                                                         &
                )

! ---------------------------------------------------------------------

      IF (BUG) THEN
        WRITE(*,*) 'BASE OK'
        ENDIF

!**** UPDATE CANOPY INTERCEPTION; DETERMINE THROUGHFALL RATES.

      CALL INTERC (                                                            &
                   NCH, ITYP, DTSTEP, TRAINLX, TRAINCX, SMELT,                 &
                   SATCAP, CSOIL, SFRAC,BUG,                                   &
                   CAPAC,                                                      &
                   THRU                                                        &
                  )

      IF (BUG) THEN
        WRITE(*,*) 'INTERC OK'
        ENDIF

! ---------------------------------------------------------------------
!!ag (02Jan2021)
!!**** TWO-WAY COUPLING WITH ROUTING SCHEME: UPDATE THROUGHFALL RATES WITH SURFACE WATER STORAGE.
!n=640
!print*,1,'thru',n,thru(n),count(fldfrc==1),count(fldfrc==0)
!      do n=1,nch
!        !if flooded fraction flag is 1, i.e., if flooded fraction is above threshold, add river and flood storages to QINSUR
!!print*,n,fldfrc(n),rivsto(n),fldsto(n)
!        if(fldfrc(n)==1)then
!          thru(n) = thru(n) + (rivsto(n) + fldsto(n))*1000. !surface water storage units are in m/s (See HYMAP2_routing_run.F90 and clsmf25_getsws_hymap2.F90)
!        endif
!      enddo
!n=640 
!print*,1,'runsrf',n,runsrf(n)
! ---------------------------------------------------------------------

!**** DETERMINE SURFACE RUNOFF AND INFILTRATION RATES:

      CALL SRUNOFF (                                                           &
                    NCH,DTSTEP,AR1,ar2,ar4,THRU,frice,tp1,srfmx,BUG,           &
                    SRFEXC,RUNSRF,                                             &
                    QINFIL                                                     &
                   )

      IF (BUG) THEN
        WRITE(*,*) 'SRUNOFF'
        ENDIF
!n=670
!print*,2,'thru',n,thru(n)
!print*,2,'runsrf',n,runsrf(n)
!print*,2,'QINFIL',n,QINFIL(n)
!print*,2,'BFLOW',n,BFLOW(n)

!**** (ADD CHECK TO ENSURE RZEXC KEPT WITHIN BOUNDS IN SRUNOFF)
      
!**** RECOMPUTE RZEXC:

      CALL RZEQUIL (                                                           &
                    NCH, ITYP, CATDEF, VGWMAX,CDCR1,CDCR2,WPWET,               &
                    ars1,ars2,ars3,ara1,ara2,ara3,ara4,arw1,arw2,arw3,arw4,    &
                    RZEQ                                                       &
                   )

      IF (BUG) THEN
        WRITE(*,*) 'RZEQUIL'
        ENDIF

      DO N=1,NCH
        ADJ=0.5*(RZEQOL(N)-RZEQ(N))
        RZEXC(N)=RZEXC(N)+ADJ
        CATDEF(N)=CATDEF(N)+ADJ
        ! make sure catdef does not become negative
        ! reichle, Aug 16, 2002
        IF(CATDEF(N) .LT. 0.) THEN
           RUNSRF(N)=RUNSRF(N)-CATDEF(N)/DTSTEP
           CATDEF(N)=0.
           ENDIF
         ENDDO

!**** Correct energy imbalance due to changing areas:
  
      CALL CALC_SOIL_MOIST (                                                   &
          nch,ityp,dzsf,vgwmax,cdcr1,cdcr2,wpwet,poros,                        &
          psis,bee,ars1,ars2,ars3,ara1,ara2,ara3,ara4,arw1,arw2,arw3,arw4,     &
          srfexc,rzexc,catdef,                                                 &
          sfmc, rzmc, prmc,                                                    &
          werror, sfmcun, rzmcun, prmcun  )

! Add differences due to adjustments to land moisture prognostics
      do n=1,nch
         if(werror(n) .le. 0.) runsrf(n)=runsrf(n)-werror(n)/dtstep
         if(werror(n) .gt. 0.) then
           edif=werror(n)/dtstep
           EVAP(N)=EVAP(N)-EDIF
           HLATN(N)=HLATN(N)-EDIF*ALHE
           EVEGFRC=EVEG(N)/(EVEG(N)+ESOI(N)+1.E-20)
           EVEG(N)=EVEG(N)-EDIF*EVEGFRC*DTSTEP
           ESOI(N)=ESOI(N)-EDIF*(1.-EVEGFRC)*DTSTEP
           SHFLUX(N)=SHFLUX(N)+EDIF*ALHE
!           EVACC(N)=EVACC(n)-EDIF
!           SHACC(N)=SHACC(N)+EDIF*ALHE
           endif
         enddo


      CALL PARTITION (                                                         &
                      NCH,DTSTEP,ITYP,DZSF,RZEXC,  RZEQOL,VGWMAX,CDCR1,CDCR2,  &
                      PSIS,BEE,poros,WPWET,                                    &
                      ars1,ars2,ars3,ara1,ara2,ara3,ara4,                      &
                      arw1,arw2,arw3,arw4,BUG,                                 &
                      SRFEXC,CATDEF,RUNSRF,                                    &
                      AR1, AR2, AR4,srfmx,srfmn,  SWSRF1,SWSRF2,SWSRF4,RZI     &
                     )

      do n=1,nch
        hold=csoil(n)*(ar1old(n)*tc1(n)+ar2old(n)*tc2(n)+ar4old(n)*tc4(n))
        hnew=csoil(n)*(ar1(n)*tc1(n)+ar2(n)*tc2(n)+ar4(n)*tc4(n))
        shflux(n)=shflux(n)-(hnew-hold)/dtstep
!        SHACC(N)=SHACC(N)-(hnew-hold)/dtstep
        enddo
 
! ---------------------------------------------------------------------
!!ag (02Jan2021)
!!**** TWO-WAY COUPLING WITH ROUTING SCHEME: UPDATE THROUGHFALL RATES WITH SURFACE WATER STORAGE.
!n=670
!print*,3,'runsrf',n,runsrf(n)
!print*
!      do n=1,nch
!        !if flooded fraction flag is 0, i.e., if flooded fraction is below threshold, add river and flood storages to RUNSRF after vertical water balance
!        if(fldfrc(n)==0)then
!          runsrf(n) = runsrf(n) + (rivsto(n) + fldsto(n))*1000. !surface water storage units are in m/s (See HYMAP2_routing_run.F90 and clsmf25_getsws_hymap2.F90)
!        endif
!      enddo

!**** ---------------------------------------------------
!**** PROCESS DATA AS NECESSARY PRIOR TO RETURN:
!****
!**** ---------------------------------------------------


      DO N=1,NCH

        RUNOFF(N) = RUNSRF(N)+BFLOW(N)
        IF(CAPAC(N).LT.1.E-10) THEN
           RUNOFF(N) = RUNOFF(N)+CAPAC(N)/DTSTEP
           CAPAC(N) = 0.0
           endif

        EINT(N) = EINT(N) * ALHE / DTSTEP
        ESOI(N) = ESOI(N) * ALHE / DTSTEP
        EVEG(N) = EVEG(N) * ALHE / DTSTEP
        ESNO(N) = ESNO(N) * ALHS / DTSTEP
         
        TSOIL=AR1(N)*TC1(N)+AR2(N)*TC2(N)+AR4(N)*TC4(N)
        TSURF(N)=(1.-ASNOW(N))*TSOIL+ASNOW(N)*TPSN1(N)

        if(asnow(n) .eq. 0) then
          tpsn1(n)=TSURF(N)
          tpsn2(n)=TSURF(N)
          tpsn3(n)=TSURF(N)
       endif

        WTOT(N) = WESNN(1,N) + WESNN(2,N) + WESNN(3,N) +                       &
               (CDCR2(N)-CATDEF(N)) + RZEXC(N) +  SRFEXC(N) + CAPAC(N)

        ENTOT(N) = HTSNNN(1,N) + HTSNNN(2,N) + HTSNNN(3,N) +                   &
            TSOIL*CSOIL(N) + GHTCNT(1,N) + GHTCNT(2,N) +                       &
            GHTCNT(3,N) + GHTCNT(4,N) + GHTCNT(5,N) + GHTCNT(6,N)

        WCHANGE(N) = (WTOT(N)-WTOT_ORIG(N))/DTSTEP
        ECHANGE(N) = (ENTOT(N)-ENTOT_ORIG(N))/DTSTEP


! Perform check on sum of AR1 and AR2, to avoid calculation of negative 
! wilting fraction due to roundoff, outside of catchment:
        IF(AR1(N)+AR2(N) .GT. 1.) THEN
          IF(AR1(N) .GE. AR2(N)) AR1(N)=1.-AR2(N)
          IF(AR1(N) .LT. AR2(N)) AR2(N)=1.-AR1(N)
          ENDIF


!       Revise values of qa and tc prior to return to AGCM, to ensure that 
!       AGCM uses evaporation and sensible heat flux rates consistent
!       with those of land model (in situations of imposed water limits, 
!       etc.).

        QA1X=QA1_ORIG(N)
        QA2X=QA2_ORIG(N)
        QA4X=QA4_ORIG(N)
        TC1X=TC1_ORIG(N)
        TC2X=TC2_ORIG(N)
        TC4X=TC4_ORIG(N)
        TCSX=TCS_ORIG(N)

!       Make uncorrected values available for optional output 
!       (needed for off-line system)

        IF (PRESENT(TC1_0))  TC1_0(N)=TC1(N)
        IF (PRESENT(TC2_0))  TC2_0(N)=TC2(N)
        IF (PRESENT(TC4_0))  TC4_0(N)=TC4(N)
        IF (PRESENT(QA1_0))  QA1_0(N)=QA1(N)
        IF (PRESENT(QA2_0))  QA2_0(N)=QA2(N)
        IF (PRESENT(QA4_0))  QA4_0(N)=QA4(N)    
        
        

        IF(DEDQA1(N) .NE. 0.) QA1X = QA1_ORIG(N)+(EVAP1(N)-ETURB1(N))/DEDQA1(N)
        IF(DEDQA2(N) .NE. 0.) QA2X = QA2_ORIG(N)+(EVAP2(N)-ETURB2(N))/DEDQA2(N)
        IF(DEDQA4(N) .NE. 0.) QA4X = QA4_ORIG(N)+(EVAP4(N)-ETURB4(N))/DEDQA4(N)

        IF(DHSDTC1(N) .NE. 0.) TC1X = TC1_ORIG(N)+(SHFLUX1(N)-HSTURB1(N))/     &
                           DHSDTC1(N)
        IF(DHSDTC2(N) .NE. 0.) TC2X = TC2_ORIG(N)+(SHFLUX2(N)-HSTURB2(N))/     &
                           DHSDTC2(N)
        IF(DHSDTC4(N) .NE. 0.) TC4X = TC4_ORIG(N)+(SHFLUX4(N)-HSTURB4(N))/     &
                           DHSDTC4(N)

!       Ensure that modifications made to QA and TC are not too large:

        IF(ABS(QA1X-QA1(N)) .LE. 0.5*QA1(N)) THEN
            QA1(N)=QA1X
          ELSE
            QA1(N)=QA1(N)+SIGN(0.5*QA1(N),QA1X-QA1(N))
          ENDIF

         IF(ABS(QA2X-QA2(N)) .LE. 0.5*QA2(N)) THEN
            QA2(N)=QA2X
          ELSE
            QA2(N)=QA2(N)+SIGN(0.5*QA2(N),QA2X-QA2(N))
          ENDIF

        IF(ABS(QA4X-QA4(N)) .LE. 0.5*QA4(N)) THEN
            QA4(N)=QA4X
          ELSE
            QA4(N)=QA4(N)+SIGN(0.5*QA4(N),QA4X-QA4(N))
          ENDIF

        IF(ABS(TC1X-TC1(N)) .LE. 10.) THEN
            TC1(N)=TC1X
          ELSE
            TC1(N)=TC1(N)+SIGN(10.,TC1X-TC1(N))
          ENDIF

        IF(ABS(TC2X-TC2(N)) .LE. 10.) THEN
            TC2(N)=TC2X
          ELSE
            TC2(N)=TC2(N)+SIGN(10.,TC2X-TC2(N))
          ENDIF

        IF(ABS(TC4X-TC4(N)) .LE. 10.) THEN
            TC4(N)=TC4X
          ELSE
            TC4(N)=TC4(N)+SIGN(10.,TC4X-TC4(N))
          ENDIF


! EVACC and SHACC are the differences ("errors") between what the land surface
! model computes for the evaporation and sensible heat flux and what the AGCM
! will compute, since the latter is forced to compute these fluxes
! based on changes in near-surface humidity and temperature only -- and 
! because the fractions (AR1, AR2, AR4, and ASNOW0) provided back to the
! AGCM are not the same as those that went into computing the land model's
! fluxes, since the land model has to update those areas (based on the fluxes)
! as a matter of course.

        EVAPX1=ETURB1(N)+DEDQA1(N)*(QA1(N)-QA1_ORIG(N))
        EVAPX2=ETURB2(N)+DEDQA2(N)*(QA2(N)-QA2_ORIG(N))
        EVAPX4=ETURB4(N)+DEDQA4(N)*(QA4(N)-QA4_ORIG(N))
        EVAPXS=ETURBS(N)+DEDQAS(N)*DQSS(N)*(TPSN1(N)-TCS_ORIG(N))
        EVACC(N)=        (1.-ASNOW0(N))*                                       &
                        ( AR1(N)*EVAPX1+                                       &
                          AR2(N)*EVAPX2+                                       &
                          AR4(N)*EVAPX4 )                                      &
                      +  ASNOW0(N)*EVAPXS
        EVACC(N)=EVAP(N)-EVACC(N)


        ! added term for latent heat flux correction, reichle+qliu, 9 Oct 2008

        if(present(lhacc)) then
           LHACC(N)=  ALHE*(1.-ASNOW0(N))*                           &
                ( AR1(N)*EVAPX1+                                     &
                  AR2(N)*EVAPX2+                                     &
                  AR4(N)*EVAPX4 )                                    &
                + ALHS*ASNOW0(N)*EVAPXS
           LHACC(N)=HLATN(N)-LHACC(N)
           end if

        SHFLUXX1=HSTURB1(N)+DHSDTC1(N)*(TC1(N)-TC1_ORIG(N))
        SHFLUXX2=HSTURB2(N)+DHSDTC2(N)*(TC2(N)-TC2_ORIG(N))
        SHFLUXX4=HSTURB4(N)+DHSDTC4(N)*(TC4(N)-TC4_ORIG(N))
        SHFLUXXS=HSTURBS(N)+DHSDTCS(N)*(TPSN1(N)-TCS_ORIG(N))
        SHACC(N)=         (1.-ASNOW0(N))*                                      &
                        ( AR1(N)*SHFLUXX1+                                     &
                          AR2(N)*SHFLUXX2+                                     &
                          AR4(N)*SHFLUXX4 )                                    &
                      + ASNOW0(N)*SHFLUXXS
        SHACC(N)=SHFLUX(N)-SHACC(N)





! **** SPECIAL DIAGNOSTICS FOR AR5 DECADAL RUNS

        if (present(SH_SNOW))  SH_SNOW(N)=-ASNOW(N)*SHFLUXS(N)
        
        if (present(AVET_SNOW)) then
           CALL GET_TF_ND(N_SNOW, HTSNNN(:,N), WESNN(:,N), TPSN, FICES)
           AVET_SNOW(N)=(TPSN(1)+TF)*WESNN(1,N) + (TPSN(2)+TF)*WESNN(2,N) +       &
                (TPSN(3)+TF)*WESNN(3,N)
        end if

        if (present(WAT_10CM))  WAT_10CM(N)=0.1*(RZEQ(N)+RZEXC(N))+SRFEXC(N)

        if (present(TOTWAT_SOIL))  &
             TOTWAT_SOIL(N)=(CDCR2(N)/(1.-WPWET(N))-CATDEF(N)+RZEXC(N)+SRFEXC(N))

        if (present(TOTICE_SOIL)) then
           if (present(TOTWAT_SOIL)) then
              
              TOTICE_SOIL(N)=TOTWAT_SOIL(N)*FRICE(N)

           else
        
              TOTICE_SOIL(N)=1.e14    ! replace with MAPL no-data-value later
              
           end if
        end if

        do i = 1,N_SNOW
           fices_out(N,i) = fices(i)
        enddo

        ENDDO     

        
      do n=1,0 !size(n_outs)
         n_out =  n_outs(n)

         write (*,*) 'OUTPUT catchment arguments for n_out: ',n_out

         write (*,*) 'NCH = ',NCH
         write (*,*) 'DTSTEP = ',DTSTEP
         write (*,*) 'SFRAC = ',SFRAC
         write (*,*) 'ITYP = ',ITYP(n_out)
         write (*,*) 'TRAINC = ',TRAINC(n_out)
         write (*,*) 'TRAINL = ',TRAINL(n_out)
         write (*,*) 'TSNOW = ',TSNOW(n_out)
         write (*,*) 'UM = ',UM(n_out)
         write (*,*) 'ETURB1 = ',ETURB1(n_out)
         write (*,*) 'DEDQA1 = ',DEDQA1(n_out)
         write (*,*) 'DEDTC1 = ',DEDTC1(n_out)
         write (*,*) 'HSTURB1 = ',HSTURB1(n_out)
         write (*,*) 'DHSDQA1 = ',DHSDQA1(n_out)
         write (*,*) 'DHSDTC1 = ',DHSDTC1(n_out)
         write (*,*) 'ETURB2 = ',ETURB2(n_out)
         write (*,*) 'DEDQA2 = ',DEDQA2(n_out)
         write (*,*) 'DECTC2 = ',DEDTC2(n_out)
         write (*,*) 'HSTURB2 = ',HSTURB2(n_out)
         write (*,*) 'DHSDQA2 = ',DHSDQA2(n_out)
         write (*,*) 'DHSDTC2 = ',DHSDTC2(n_out)
         write (*,*) 'ETURB4 = ',ETURB4(n_out)
         write (*,*) 'DEDQA4 = ',DEDQA4(n_out)
         write (*,*) 'DEDTC4 = ',DEDTC4(n_out)
         write (*,*) 'HSTURB4 = ',HSTURB4(n_out)
         write (*,*) 'HDSDQA4 = ',DHSDQA4(n_out)
         write (*,*) 'DHSDTC4 = ',DHSDTC4(n_out)
         write (*,*) 'ETURBS = ',ETURBS(n_out)
         write (*,*) 'DEDQAS = ',DEDQAS(n_out)
         write (*,*) 'DEDTCS = ',DEDTCS(n_out)
         write (*,*) 'HSTURBS = ',HSTURBS(n_out)
         write (*,*) 'DHSDQAS = ',DHSDQAS(n_out)
         write (*,*) 'DHSDTCS = ',DHSDTCS(n_out)
         write (*,*) 'TM = ',TM(n_out)
         write (*,*) 'QM = ',QM(n_out)
         write (*,*) 'ra1 = ',ra1(n_out)
         write (*,*) 'ra2 = ',ra2(n_out)
         write (*,*) 'ra4 = ',ra4(n_out)
         write (*,*) 'raS = ',raS(n_out)
         write (*,*) 'SUNANG = ',SUNANG(n_out)
         write (*,*) 'PARDIR = ',PARDIR(n_out)
         write (*,*) 'PARDIF = ',PARDIF(n_out)
         write (*,*) 'SWNETF = ',SWNETF(n_out)
         write (*,*) 'SWNETS = ',SWNETS(n_out)
         write (*,*) 'HLWDWN = ',HLWDWN(n_out)
         write (*,*) 'PSUR = ',PSUR(n_out)
         write (*,*) 'ZLAI = ',ZLAI(n_out)
         write (*,*) 'GREEN = ',GREEN(n_out)
         write (*,*) 'Z2 = ',Z2(n_out)
         write (*,*) 'SQSCAT = ',SQSCAT(n_out)
         write (*,*) 'RSOIL1 = ',RSOIL1(n_out)
         write (*,*) 'RSOIL2 = ',RSOIL2(n_out)
         write (*,*) 'RDC = ',RDC(n_out)
!         write (*,*) 'USFAC = ',U2FAC(n_out)
         write (*,*) 'QSAT1 = ',QSAT1(n_out)
         write (*,*) 'DQS1 = ',DQS1(n_out)
         write (*,*) 'ALW1 = ',ALW1(n_out)
         write (*,*) 'BLW1 = ',BLW1(n_out)
         write (*,*) 'QSAT2 = ',QSAT2(n_out)
         write (*,*) 'DQS2 = ',DQS2(n_out)
         write (*,*) 'ALW2 = ',ALW2(n_out)
         write (*,*) 'BLW2 = ',BLW2(n_out)
         write (*,*) 'QSAT4 = ',QSAT4(n_out)
         write (*,*) 'DQS4 = ',DQS4(n_out)
         write (*,*) 'ALW4 = ',ALW4(n_out)
         write (*,*) 'BLW4 = ',BLW4(n_out)
         write (*,*) 'QSATS = ',QSATS(n_out)
         write (*,*) 'DQSS = ',DQSS(n_out)
         write (*,*) 'ALWS = ',ALWS(n_out)
         write (*,*) 'BLWS = ',BLWS(n_out)
         write (*,*) 'BF1 = ',BF1(n_out)
         write (*,*) 'BF2 = ',BF2(n_out)
         write (*,*) 'BF3 = ',BF3(n_out)
         write (*,*) 'VGWMAX = ',VGWMAX(n_out)
         write (*,*) 'CDCR1 = ',CDCR1(n_out)
         write (*,*) 'CDCR2 = ',CDCR2(n_out)
         write (*,*) 'psis = ',psis(n_out)
         write (*,*) 'bee = ',bee(n_out)
         write (*,*) 'poros = ',poros(n_out)
         write (*,*) 'wpwet = ',wpwet(n_out)
         write (*,*) 'cond = ',cond(n_out)
         write (*,*) 'gnu = ',gnu(n_out)
         write (*,*) 'ars1 = ',ARS1(n_out)
         write (*,*) 'ars2 = ',ARS2(n_out)
         write (*,*) 'ars3 = ',ARS3(n_out)
         write (*,*) 'ara1 = ',ARA1(n_out)
         write (*,*) 'ara2 = ',ARA2(n_out)
         write (*,*) 'ara3 = ',ARA3(n_out)
         write (*,*) 'ara4 = ',ARA4(n_out)
         write (*,*) 'arw1 = ',ARW1(n_out)
         write (*,*) 'arw2 = ',ARW2(n_out)
         write (*,*) 'arw3 = ',ARW3(n_out)
         write (*,*) 'arw4 = ',ARW4(n_out)
         write (*,*) 'tsa1 = ',tsa1(n_out)
         write (*,*) 'tsa2 = ',tsa2(n_out)
         write (*,*) 'tsb1 = ',tsb1(n_out)
         write (*,*) 'tsb2 = ',tsb2(n_out)
         write (*,*) 'atau = ',atau(n_out)
         write (*,*) 'btau = ',btau(n_out)
         write (*,*) 'BUG = ',BUG
         write (*,*) 'TC1 = ',TC1(n_out)
         write (*,*) 'TC2 = ',TC2(n_out)
         write (*,*) 'TC4 = ',TC4(n_out)
         write (*,*) 'QA1 = ',QA1(n_out)
         write (*,*) 'QA2 = ',QA2(n_out)
         write (*,*) 'QA4 = ',QA4(n_out)
         write (*,*) 'CAPAC = ',CAPAC(n_out)
         write (*,*) 'CATDEF = ',CATDEF(n_out)
         write (*,*) 'RZEXC = ',RZEXC(n_out)
         write (*,*) 'SRFEXC = ',SRFEXC(n_out)
         write (*,*) 'GHTCNT = ',GHTCNT(:,n_out)
         write (*,*) 'TSURF = ',TSURF(n_out)
         write (*,*) 'WESNN = ',WESNN(:,n_out)
         write (*,*) 'HTSNNN = ',HTSNNN(:,n_out)
         write (*,*) 'SNDZN = ',SNDZN(:,n_out)
         write (*,*) 'EVAP = ',EVAP(n_out)
         write (*,*) 'SHFLUX = ',SHFLUX(n_out)
         write (*,*) 'RUNOFF = ',RUNOFF(n_out)
         write (*,*) 'EINT = ',EINT(n_out)
         write (*,*) 'ESOI = ',ESOI(n_out)
         write (*,*) 'EVEG = ',EVEG(n_out)
         write (*,*) 'ESNO = ',ESNO(n_out)
         write (*,*) 'BFLOW = ',BFLOW(n_out)
         write (*,*) 'RUNSRF = ',RUNSRF(n_out)
         write (*,*) 'SMELT = ',SMELT(n_out)
         write (*,*) 'HLWUP = ',HLWUP(n_out)
         write (*,*) 'HLATN = ',HLATN(n_out)
         write (*,*) 'QINFIL = ',QINFIL(n_out)
         write (*,*) 'AR1 = ',AR1(n_out)
         write (*,*) 'AR2 = ',AR2(n_out)
         write (*,*) 'RZEQ = ',RZEQ(n_out)
         write (*,*) 'GHFLUX = ',GHFLUX(n_out)
         write (*,*) 'TPSN1 = ',TPSN1(n_out)
         write (*,*) 'ASNOW0 = ',ASNOW0(n_out)
         write (*,*) 'TP1 = ',TP1(n_out)
         write (*,*) 'TP2 = ',TP2(n_out)
         write (*,*) 'TP3 = ',TP3(n_out)
         write (*,*) 'TP4 = ',TP4(n_out)
         write (*,*) 'TP5 = ',TP5(n_out)
         write (*,*) 'TP6 = ',TP6(n_out)
      enddo
!rr ------------------------------------------------------------------      



      RETURN
      END SUBROUTINE clsmf25

!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN RCUNST ]
!****
      SUBROUTINE RCUNST (                                                      &
                         NCH, ITYP, SUNANG, SQSCAT, PDIR,                      &
                         PAR, ZLAI, GREEN,BUG,                                 &
                         RCUN                                                  &
                        )
!****
!****     This subroutine calculates the unstressed canopy resistance.
!**** (p. 1353, Sellers 1985.)  Extinction coefficients are computed first.
!****
      IMPLICIT NONE
!****
      LOGICAL, INTENT(IN) :: BUG
      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN), DIMENSION(NCH) :: SUNANG, PDIR, PAR, ZLAI,             &
            SQSCAT, GREEN

      REAL, INTENT(OUT), DIMENSION(NCH) :: RCUN

      REAL, DIMENSION(NTYPS) :: VGCHIL, VGZMEW, VGRST1, VGRST2, VGRST3


      INTEGER CHNO
      REAL  RHO4, EXTK1, EXTK2, RCINV, GAMMA, EKAT, DUM1, DUM2, DUM3,          &
              AA, BB, ZK, CC


      DATA VGCHIL /        0.1,        0.25,        0.01,        -0.3,         &
                          0.01,        0.20,         0.0,         0.0,         &
                           0.0,         0.0 /

      DATA VGZMEW/      0.9809,      0.9638,      0.9980,      1.0773,         &
                        0.9980,      0.9676,       1.000,       1.000,         &
                         1.000,       1.000 /

      DATA VGRST1 /     2335.9,      9802.2,      2869.7,      2582.0,         &
                       93989.4,      9802.2,         0.0,         0.0,         &
                           0.0,         0.0 /

      DATA VGRST2 /        0.0,        10.6,         3.7,         1.1,         &
                          0.01,        10.6,         0.0,         0.0,         &
                           0.0,         0.0 /

      DATA VGRST3 /      153.5,       180.0,       233.0,       110.0,         &
                         855.0,       180.0,         1.0,         1.0,         &
                           1.0,         1.0 /



      DO 100 ChNo = 1, NCH

!**** First compute optical parameters.
!**** (Note: CHIL is constrained to be >= 0.01, as in SiB calcs.)

      AA = 0.5 - (0.633 + 0.330*VGCHIL(ITYP(ChNo)))*VGCHIL(ITYP(ChNo))
      BB = 0.877 * ( ONE - 2.*AA )
      CC =  ( AA + BB*SUNANG(ChNo) ) / SUNANG(ChNo)

      EXTK1 =  CC * SQSCAT(ChNo)
      EXTK2 = (ONE / VGZMEW(ITYP(ChNo))) * SQSCAT(ChNo)

      DUM1 =      PDIR(ChNo)  *   CC
      DUM2 = (ONE-PDIR(ChNo)) * ( BB*(ONE/3.+PIE/4.) + AA*1.5 )

!**** Bound extinction coefficient by 50./ZLAI:

      ZK =     PDIR(ChNo) *AMIN1( EXTK1, 50./ZLAI(ChNo) ) +                    &
          (ONE-PDIR(ChNo))*AMIN1( EXTK2, 50./ZLAI(ChNo) )

!**** Now compute unstressed canopy resistance:

      GAMMA = VGRST1(ITYP(ChNo)) / VGRST3(ITYP(ChNo)) + VGRST2(ITYP(ChNo))

      EKAT = EXP( ZK*ZLAI(ChNo) )
      RHO4 = GAMMA / (PAR(ChNo) * (DUM1 + DUM2))

      DUM1 = (VGRST2(ITYP(ChNo)) - GAMMA) / (GAMMA + 1.E-20)
      DUM2 = (RHO4 * EKAT + ONE) / (RHO4 + ONE)
      DUM3 = ZK * VGRST3(ITYP(ChNo))

      RCINV = ( DUM1*ALOG(DUM2) + ZK*ZLAI(ChNo) ) / DUM3         
      rcinv = amax1(rcinv,0.)

      RCUN(ChNo) = ONE / (RCINV * GREEN(ChNo) + 1.E-10)

 100  CONTINUE


      RETURN
      END SUBROUTINE RCUNST
!****
!**** [ END RCUNST ]

!**** ===================================================
!**** ///////////////////////////////////////////////////
!**** ===================================================

      SUBROUTINE SRUNOFF (                                                     &
                          NCH,DTSTEP,AR1,ar2,ar4, THRU,frice,tp1,srfmx,BUG,    &
                          SRFEXC,RUNSRF,                                       &
                          QINFIL                                               &
                         )

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: AR1, ar2, ar4, THRU, frice, tp1,     &
             srfmx
      LOGICAL, INTENT(IN) :: BUG

      REAL, INTENT(INOUT), DIMENSION(NCH) ::  SRFEXC ,RUNSRF

      REAL, INTENT(OUT), DIMENSION(NCH) :: QINFIL

      INTEGER N
      REAL PTOTAL,srun0,frun,qin

!**** - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO N=1,NCH

        PTOTAL=THRU(N)
        frun=AR1(N)
        if(srfexc(n) .gt. 0.) then
!SVK: uncommented on 3/10/14
           frun=frun+ar2(n)*(srfexc(n)/(srfmx(n)+1.e-20))**2
           frun=frun+ar4(n)*(srfexc(n)/(srfmx(n)+1.e-20))**4
          endif
!        frun=frun+(1-frun)*frice(n)
        srun0=PTOTAL*frun

!**** Comment out this line in order to allow moisture
!**** to infiltrate soil:
!       if(tp1(n) .lt. 0.) srun0=ptotal

        if(ptotal-srun0 .gt. srfmx(n)-srfexc(n))                               &
                      srun0=ptotal-(srfmx(n)-srfexc(n)) 

        if (srun0 .gt. ptotal) then
!rr          write(*,*) 'srun0 > ptotal: N=',N
!rr          write(*,*) ' frice=',frice(n),' ar1=',ar1(n),' ptotal=',
!rr     &           ptotal,' tp1=',tp1(n)
!rr          write(*,*) ' ar2=',ar2(n),' ar4=',ar4(n),' srfexc=',
!rr     &           srfexc(n),' srfmx=',srfmx(n),' thru=',thru(n),
!rr     &           ' rzexc=',rzexc(n)
!rr          write(*,*) '=====> CORRECTION'
          srun0=ptotal
          endif

        RUNSRF(N)=RUNSRF(N)+srun0
        QIN=PTOTAL-srun0

        SRFEXC(N)=SRFEXC(N)+QIN
        RUNSRF(N)=RUNSRF(N)/DTSTEP
        QINFIL(N)=QIN/DTSTEP
         
        ENDDO
      
      RETURN
      END SUBROUTINE SRUNOFF

!**** ===================================================
!**** ///////////////////////////////////////////////////
!**** ===================================================

      SUBROUTINE RZDRAIN (                                                     &
                          NCH,DTSTEP,VGWMAX,SATCAP,RZEQ,AR1,WPWET,             &
                          tsa1,tsa2,tsb1,tsb2,atau,btau,CDCR2,poros,BUG,       &
                          CAPAC,RZEXC,SRFEXC,CATDEF,RUNSRF                     &
                         )

!-----------------------------------------------------------------
!        defines drainage timescales:
!             - tsc0, between srfex and rzex
!             - tsc2, between rzex and catdef
!        then defines correponding drainages
!        and updates the water contents
!-----------------------------------------------------------------

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN) ::  DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: VGWMAX, SATCAP, RZEQ, AR1, wpwet,    &
              tsa1, tsa2, tsb1, tsb2, atau, btau, CDCR2, poros
      LOGICAL, INTENT(IN) :: BUG

      REAL, INTENT(INOUT), DIMENSION(NCH) :: RZEXC, SRFEXC, CATDEF, CAPAC,     &
              RUNSRF


      INTEGER N
      REAL srflw,rzflw,FLOW,EXCESS,TSC0,tsc2,rzave,rz0,wanom,rztot,            &
            rzx,btaux,ax,bx,rzdif


!**** - - - - - - - - - - - - - - - - - - - - - - - - - 

      DO 100 N=1,NCH

!****   Compute equivalent of root zone excess in non-saturated area:
        rztot=rzeq(n)+rzexc(n)
        if(ar1(n).ne.1.) then
            rzave=(rztot-ar1(n)*vgwmax(n))/(1.-ar1(n))
            rzave=rzave*poros(n)/vgwmax(n)
          else
            rzave=poros(n)
          endif
  
! jeff 
        if (rzave .eq. 0.) then
          rzave=1.e-4
          print*,'problem: rzave = 0 in catchment',n
          end if
! jeff

        btaux=btau(n)
        if (srfexc(n) .lt. 0.) btaux=btau(n)*(poros(n)/rzave)
        rz0=amax1(0.001,rzave-srfexc(n)/(1000.*(-btaux)))
        tsc0=atau(n)/(rz0**3.)

        tsc0=tsc0*3600.
        if(tsc0.lt.dtstep) tsc0=dtstep

! ---------------------------------------------------------------------

        SRFLW=SRFEXC(N)*DTSTEP/TSC0
!rr   following inserted by koster Sep 22, 2003
        rzdif=rzave/poros(n)-wpwet(n)
!**** No moisture transport up if rz at wilting; employ ramping.
        if(rzdif.le.0. .and. srflw.lt.0.)  srflw=0.
        if(rzdif.gt.0. .and. rzdif.lt.0.01                                     &
                   .and. srflw.lt.0.) srflw=srflw*(rzdif/0.01)
        RZEXC(N)=RZEXC(N)+SRFLW
        SRFEXC(N)=SRFEXC(N)-SRFLW

!**** Topography-dependent tsc2, between rzex and catdef

        rzx=rzexc(n)/vgwmax(n)

        if(rzx .gt. .01) then
            ax=tsa1(n)
            bx=tsb1(n)
          elseif(rzx .lt. -.01) then
            ax=tsa2(n)
            bx=tsb2(n)
          else
            ax=tsa2(n)+(rzx+.01)*(tsa1(n)-tsa2(n))/.02
            bx=tsb2(n)+(rzx+.01)*(tsb1(n)-tsb2(n))/.02
          endif

        tsc2=exp(ax+bx*catdef(n))
        rzflw=rzexc(n)*tsc2*dtstep/3600.

        IF (CATDEF(N)-RZFLW .GT. CDCR2(N)) then
          RZFLW=CATDEF(N)-CDCR2(N)
          end if

        CATDEF(N)=CATDEF(N)-RZFLW
        RZEXC(N)=RZEXC(N)-RZFLW

!****   REMOVE ANY EXCESS FROM MOISTURE RESERVOIRS:

        IF(CAPAC(N) .GT. SATCAP(N)) THEN
          RZEXC(N)=RZEXC(N)+CAPAC(N)-SATCAP(N)
          CAPAC(N)=SATCAP(N)
          ENDIF

        IF(RZEQ(N) + RZEXC(N) .GT. VGWMAX(N)) THEN
          EXCESS=RZEQ(N)+RZEXC(N)-VGWMAX(N)
          RZEXC(N)=VGWMAX(N)-RZEQ(N)
          CATDEF(N)=CATDEF(N)-EXCESS
          ENDIF

        IF(CATDEF(N) .LT. 0.) THEN
          RUNSRF(N)=RUNSRF(N)-CATDEF(N)
          CATDEF(N)=0.
          ENDIF

  100 ENDDO

      RETURN
      END SUBROUTINE RZDRAIN
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN WUPDAT ]
!****
      SUBROUTINE WUPDAT (                                                      &
                           NCH,   ITYP, DTSTEP,  EVAP, SATCAP, TC, RA, RC,     &
                           RX11,RX21,RX12,RX22,RX14,RX24, AR1,AR2,AR4,CDCR1,   &
                           EIRFRC,RZEQ,srfmn,WPWET,VGWMAX,                     &
                           CAPAC, RZEXC, CATDEF, SRFEXC,                       &
                           EINT, ESOI, EVEG                                    &
                          )
!****
!**** THIS SUBROUTINE ALLOWS EVAPOTRANSPIRATION TO ADJUST THE WATER
!**** CONTENTS OF THE INTERCEPTION RESERVOIR AND THE SOIL LAYERS.
!****
      IMPLICIT NONE
!****
      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP
      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: EVAP, SATCAP, TC, RA, RC, RX11,      &
             RX21, RX12, RX22, RX14, RX24, AR1, AR2, AR4, CDCR1, EIRFRC,       &
             RZEQ, srfmn, WPWET, VGWMAX

      REAL, INTENT(INOUT), DIMENSION(NCH) :: CAPAC, CATDEF, RZEXC, SRFEXC

      REAL, INTENT(OUT), DIMENSION(NCH) :: EINT, ESOI, EVEG


      INTEGER CHNO
      REAL EGRO, CNDSAT, CNDUNS, ESATFR, cndv, cnds, WILT, egromx,rzemax

!****
!**** -----------------------------------------------------------------
      DO 100 CHNO = 1, NCH

!**** COMPUTE EFFECTIVE SURFACE CONDUCTANCES IN SATURATED AND UNSATURATED
!**** AREAS:

      CNDSAT=(AR1(CHNO)/RX11(CHNO)) + (AR1(CHNO)/RX21(CHNO))
      CNDUNS=(AR2(CHNO)/RX12(CHNO)) + (AR2(CHNO)/RX22(CHNO)) +                 &
             (AR4(CHNO)/RX14(CHNO)) + (AR4(CHNO)/RX24(CHNO))
     
      ESATFR=CNDSAT/(CNDSAT+CNDUNS)

!****
!**** PARTITION EVAP BETWEEN INTERCEPTION AND GROUND RESERVOIRS.
!****

      EINT(CHNO)=EIRFRC(CHNO)*EVAP(CHNO)*DTSTEP
      EGRO = EVAP(CHNO)*DTSTEP - EINT(CHNO)

!**** ENSURE THAT INDIVIDUAL CAPACITIES ARE NOT EXCEEDED.

      IF(EINT(CHNO) .GT. CAPAC(CHNO)) THEN
        EGRO=EGRO+EINT(CHNO)-CAPAC(CHNO)
        EINT(CHNO)=CAPAC(CHNO)
        ENDIF

! RK 09/16/03
      WILT=WPWET(CHNO)*VGWMAX(CHNO)
      rzemax=amax1(0.,RZEXC(CHNO)+RZEQ(CHNO)-WILT )
      egromx= rzemax + (srfexc(chno)-srfmn(chno))
      IF(EGRO .GT. egromx) THEN
! 06.02.98: the minimum is designed to prevent truncation errors 
        EINT(CHNO)=AMIN1(CAPAC(CHNO),EINT(CHNO)+EGRO-egromx)
        EGRO=egromx
        ENDIF
! RK 09/16/03
! RK: the above test ensures in particular that rzexc+rzeq never < 0.



      ESOI(CHNO)=AMIN1(SRFEXC(CHNO)-SRFMN(CHNO),                               &
              EGRO*(AR1(CHNO)*RX11(CHNO)/(RX11(CHNO)+RX21(CHNO)+1.E-20)        &
              +AR2(CHNO)*RX12(CHNO)/(RX12(CHNO)+RX22(CHNO)+1.E-20)             &
              +AR4(CHNO)*RX14(CHNO)/(RX14(CHNO)+RX24(CHNO)+1.E-20)))
      EVEG(CHNO)=EGRO-ESOI(CHNO)

!rdk   following inserted by koster Oct. 16, 2007
!**** special case if soil is at wilting point and rx14=rx24
!      if(rzemax .eq. 0.) then
!        esoi(chno)=AMIN1(SRFEXC(CHNO)-SRFMN(CHNO),EGRO)
!        EVEG(CHNO)=EGRO-ESOI(CHNO)
!        endif

      if(esoi(chno) .gt. SRFEXC(CHNO)-SRFMN(CHNO)) then
        esoi(chno)=SRFEXC(CHNO)-SRFMN(CHNO)
        EVEG(CHNO)=EGRO-ESOI(CHNO)
        endif

      if(eveg(chno) .gt. rzemax) then
        EVEG(CHNO)=rzemax
        esoi(chno)=egro-eveg(chno)
        endif


!****
!**** SPECIAL CASE FOR CONDENSATION:
      IF(EVAP(CHNO) .LT. 0.) THEN
        EINT(CHNO)=EVAP(CHNO)*DTSTEP
! 05.20.98: to prevent negative throughfall due to truncation errors 
        EINT(CHNO)=AMIN1(0.,EINT(CHNO))
        ESOI(CHNO)=0.
        EVEG(CHNO)=0.
        ENDIF

!****
!**** REMOVE MOISTURE FROM RESERVOIRS:
!****

      IF (CATDEF(CHNO) .LT. CDCR1(CHNO)) THEN
          CAPAC(CHNO) = AMAX1(0., CAPAC(CHNO) - EINT(CHNO))
          RZEXC(CHNO) = RZEXC(CHNO) - EVEG(CHNO)*(1.-ESATFR)
          SRFEXC(CHNO) = SRFEXC(CHNO) - ESOI(CHNO)*(1.-ESATFR)
          CATDEF(CHNO) = CATDEF(CHNO) + (ESOI(CHNO) + EVEG(CHNO))*ESATFR
! 05.12.98: first attempt to include bedrock
        ELSE
          CAPAC(CHNO) = AMAX1(0., CAPAC(CHNO) - EINT(CHNO))
          RZEXC(CHNO) = RZEXC(CHNO) -  EVEG(CHNO)
          SRFEXC(CHNO) = SRFEXC(CHNO) - ESOI(CHNO)
        ENDIF

!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE WUPDAT
!****
!**** [ END WUPDAT ]
!****
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN INTERC ]
!**** 
      SUBROUTINE INTERC (                                                      &
                         NCH, ITYP, DTSTEP, TRAINL, TRAINC,SMELT,              &
                         SATCAP, CSOIL, SFRAC,BUG,                             &
                         CAPAC,                                                &
                         THRU                                                  &
                        )
!****
!**** THIS ROUTINE USES THE PRECIPITATION FORCING TO DETERMINE 
!**** CHANGES IN INTERCEPTION AND SOIL MOISTURE STORAGE.
!**** Changes in snowcover are not treated here anymore.
!****
      IMPLICIT NONE

!****
      INTEGER, INTENT(IN) ::  NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP
      REAL, INTENT(IN) :: DTSTEP, SFRAC
      REAL, INTENT(IN), DIMENSION(NCH) :: TRAINL, TRAINC, SMELT, SATCAP,       &
            CSOIL 
      LOGICAL, INTENT(IN) :: BUG

      REAL, INTENT(INOUT), DIMENSION(NCH) :: CAPAC

      REAL, INTENT(OUT), DIMENSION(NCH) :: THRU


      INTEGER CHNO
      REAL WETINT, WATADD, CAVAIL, THRUC, TIMFRL, TIMFRC,                      &
           THRU1, THRU2, THRUL, XTCORR,SMPERS

      !DATA FWETL /0.02/, FWETC /0.02/
      DATA TIMFRL/1.00/
      DATA TIMFRC/0.333/
! value for GSWP
!      TIMFRC/0.125/  

!****
!**** ------------------------------------------------------------------
!**** LOOP OVER CHIPS:
      DO 100 CHNO = 1, NCH

!**** =======================================================
!****
!**** LOAD INTERCEPTION RESERVOIR.  STEP 1: LARGE SCALE CONDENSATION.
!****
!**** DETERMINE XTCORR, THE FRACTION OF A STORM THAT FALLS ON A PREVIOUSLY
!**** WET SURFACE DUE TO THE TIME CORRELATION OF PRECIPITATION POSITION.
!**** (TIME SCALE TIMFRL FOR LARGE SCALE STORMS SET TO ONE FOR FWETL=1
!**** TO REFLECT THE EFFECTIVE LOSS OF "POSITION MEMORY" WHEN STORM 
!**** COVERS ENTIRE GRID SQUARE.)

      XTCORR= (1.-TIMFRL) *                                                    &
            AMIN1( 1.,(CAPAC(CHNO)/SATCAP(CHNO))/(FWETL*SFRAC) )   

!****
!**** FILL INTERCEPTION RESERVOIR WITH PRECIPITATION.
!**** THRU1 IS FIRST CALCULATED AS THE AMOUNT FALLING THROUGH THE 
!****    CANOPY UNDER THE ASSUMPTION THAT ALL RAIN FALLS RANDOMLY.  
!****    ONLY A FRACTION 1-XTCORR FALLS RANDOMLY, THOUGH, SO THE RESULT 
!****    IS MULTIPLIED BY 1-XTCORR.
!****
      WATADD = TRAINL(CHNO)*DTSTEP + SMELT(CHNO)*DTSTEP
      CAVAIL = ( SATCAP(CHNO) - CAPAC(CHNO) ) * (FWETL*SFRAC)
      WETINT = CAPAC(CHNO)/SATCAP(CHNO)
      IF( WATADD*(1.-WETINT) .LT. CAVAIL ) THEN
          THRU1 = WATADD*WETINT
        ELSE
          THRU1 = (WATADD - CAVAIL)
        ENDIF
      THRU1=THRU1*(1.-XTCORR)

!**** THRU2 IS THE AMOUNT THAT FALLS IMMEDIATELY THROUGH THE CANOPY DUE
!**** TO 'POSITION MEMORY'.

      THRU2=XTCORR*WATADD

      THRUL=THRU1+THRU2

      CAPAC(CHNO)=CAPAC(CHNO)+WATADD-THRU1-THRU2

!****
!**** ---------------------------------------------------
!****
!**** STEP 2: MOIST CONVECTIVE PRECIPITATION.
!****
!**** DETERMINE XTCORR, THE FRACTION OF A STORM THAT FALLS ON A PREVIOUSLY
!**** WET SURFACE DUE TO THE TIME CORRELATION OF PRECIPITATION POSITION.
    
      XTCORR= (1.-TIMFRC) *                                                    &
           AMIN1( 1.,(CAPAC(CHNO)/SATCAP(CHNO))/(FWETC*SFRAC) )

!****
!**** FILL INTERCEPTION RESERVOIR WITH PRECIPITATION.
!**** THRU1 IS FIRST CALCULATED AS THE AMOUNT FALLING THROUGH THE 
!****    CANOPY UNDER THE ASSUMPTION THAT ALL RAIN FALLS RANDOMLY.  
!****    ONLY A FRACTION 1-XTCORR FALLS RANDOMLY, THOUGH, SO THE RESULT 
!****    IS MULTIPLIED BY 1-XTCORR.
!****
      WATADD = TRAINC(CHNO)*DTSTEP
      CAVAIL = ( SATCAP(CHNO) - CAPAC(CHNO) ) * (FWETC*SFRAC)
      WETINT = CAPAC(CHNO)/SATCAP(CHNO)
      IF( WATADD*(1.-WETINT) .LT. CAVAIL ) THEN
          THRU1 = WATADD*WETINT
        ELSE
          THRU1 = (WATADD - CAVAIL)
        ENDIF
      THRU1=THRU1*(1.-XTCORR)

!**** THRU2 IS THE AMOUNT THAT FALLS IMMEDIATELY THROUGH THE CANOPY DUE
!**** TO 'POSITION MEMORY'.

      THRU2=XTCORR*WATADD

      THRUC=THRU1+THRU2
      CAPAC(CHNO)=CAPAC(CHNO)+WATADD-THRU1-THRU2
!****
      IF (THRUL+THRUC .LT. -1.e-8) WRITE(*,*) 'THRU= ',                        &
          THRUL, THRUC, TRAINC(CHNO), TRAINL(CHNO), SMELT(CHNO)
      THRU(CHNO)=AMAX1(0., THRUL+THRUC)

 100  CONTINUE
!****
      RETURN
      END SUBROUTINE INTERC
!****
!**** [ END INTERC ]
!****
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
      SUBROUTINE BASE (                                                        &
                       NCH,DTSTEP,BF1,BF2,BF3,CDCR1,FRICE,COND,GNU,CDCR2,      &
                       CATDEF,                                                 &
                       BFLOW                                                   &
                      )

      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: BF1, BF2, BF3, CDCR1, FRICE, COND,   &
          GNU, CDCR2

      REAL, INTENT(INOUT), DIMENSION(NCH) :: CATDEF

      REAL, INTENT(OUT), DIMENSION(NCH) :: BFLOW


      INTEGER N
      REAL ZBAR, ashift

      data ashift/0./


      DO N=1,NCH
        ZBAR=SQRT(1.e-20+catdef(n)/bf1(n))-bf2(n)
! DMM: Fixing inconsistency in units between baseflow and cond. 09/01/2016
        BFLOW(N)=(1.-FRICE(N))*                                                &
              cond(n)*exp(-(bf3(n)-ashift)-gnu(n)*zbar)/gnu(n)
        IF (CATDEF(N) .GE. CDCR1(N)) BFLOW(N)=0.
        bflow(n)=amin1(cond(n),bflow(n))
! *1000 is to convert from m/s to mm/s
        bflow(n)=bflow(n)*1000.0
        CATDEF(N)=CATDEF(N)+BFLOW(N)*dtstep
        ENDDO


      RETURN
      END SUBROUTINE BASE


!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****

      SUBROUTINE PARTITION (                                                   &
                            NCH,DTSTEP,ITYP,DZSF,RZEXC,RZEQ,VGWMAX,CDCR1,CDCR2,&
                            PSIS,BEE,poros,WPWET,                              &
                            ars1,ars2,ars3,ara1,ara2,ara3,ara4,                &
                            arw1,arw2,arw3,arw4,BUG,                           &
                            srfexc,catdef,runsrf,                              &
                            AR1, AR2, AR4, srfmx, srfmn,                       &
                            SWSRF1,SWSRF2,SWSRF4,RZI                           &
                           )

      IMPLICIT NONE

! -------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: DZSF,RZEXC,RZEQ,VGWMAX,CDCR1,CDCR2,  &
                                          PSIS,BEE,poros,WPWET,                &
                                          ars1,ars2,ars3,ara1,ara2,ara3,ara4,  &
                                          arw1,arw2,arw3,arw4

      LOGICAL, INTENT(IN) :: BUG
! -------------------------------------------------------------------
      REAL, INTENT(INOUT), DIMENSION(NCH) :: srfexc,catdef,runsrf
! -------------------------------------------------------------------
      REAL, INTENT(OUT), DIMENSION(NCH) :: AR1, AR2, AR4, srfmx, srfmn,        &
                                           SWSRF1, SWSRF2, SWSRF4, RZI
! -------------------------------------------------------------------
      INTEGER :: N

      REAL :: cor, A150, W150, WMIN, AX, WMNEW, WRZ, TERM1, TERM2, TERM3,      &
              AREA0, AREA1, AREA2, AREA3, AREA4, ASCALE, WILT, D1, D2, CDI,    &
              DELTA1, DELTA2, DELTA4, MULTAR, CATDEFX, RZEQX, RZEQW, FACTOR,   &
              X0, RZEQY, CATDEFW, AR1W, ASUM, RZEQYI, RZEQWI, RZEQXI, AR20,    &
              ARG1, EXPARG1, ARG2, EXPARG2, ARG3, EXPARG3  !, surflay

      LOGICAL :: LSTRESS


      DATA LSTRESS/.FALSE./    !,surflay/20./

!****
!**** --------------------------------------------------       

!rr   next line for debugging, sep 23, 2003, reichle
!rr
!rr   write (*,*) 'entering partition()'

      DO N=1,NCH

        WILT=WPWET(N)
        WRZ=RZEXC(N)/VGWMAX(N)
        CATDEFX=AMIN1( CATDEF(N) , CDCR1(N) )

! CDI DEFINES IF THE SHAPE PARAMETER IS ADJUSTED IN ONE OR TWO SEGMENTS
        if (ara1(n) .ne. ara3(n)) then
            cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
          else
            cdi=0.
          endif

          if((1.+ars2(n)*CATDEFX+ars3(n)*CATDEFX*CATDEFX).eq.0) then 
             print*, 'divide by zero issue: ',n, catdefx, cdcr1(n), &
                  ars2(n), ars3(n)
          endif
        AR1(N)= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFX)                         &
                 /(1.+ars2(n)*CATDEFX+ars3(n)*CATDEFX*CATDEFX))) 

        if (CATDEFX .ge. cdi) then
            ax=ara3(n)*CATDEFX+ara4(n)
          else
            ax=ara1(n)*CATDEFX+ara2(n)
          endif

        WMIN=AMIN1(1.,AMAX1(0.,arw4(n)+(1.-arw4(n))*                           &
                 (1.+arw1(n)*CATDEFX)                                          &
                 /(1.+arw2(n)*CATDEFX+arw3(n)*CATDEFX*CATDEFX)))

!**** CRITICAL VALUE 1: AVERAGE MOISTURE IN ROOT ZONE AT WMIN
!**** ASSOCIATED WITH CATDEF.
        ARG1=AMAX1(-40., AMIN1(40., -AX*(1.-WMIN)))
        EXPARG1=EXP(ARG1)
        RZEQX=(WMIN-1.-(2./AX))*EXPARG1 + WMIN + (2./AX)
        RZEQXI=AX*EXPARG1 *                                                    &
          ( -1. -(2./AX) - (2./(AX*AX)) + WMIN + (WMIN/AX) )                   &
          + WMIN + 2./AX
        AR20=1.+(-AX-1.+AX*WMIN)*EXPARG1
!        if(ar20.eq.0) then 
!           print*, LIS_localPet, n,catdefx
!           print*, 'ar20 check ',LIS_localPet, n,catdefx
!           stop
!        endif
        RZEQXI=RZEQXI/(AR20+1.E-20)
!**** CRITICAL VALUE 2: AVERAGE MOISTURE IN ROOT ZONE WHEN WMIN
!**** IS EXACTLY AT WILTING POINT.
        ARG2=AMAX1(-40., AMIN1(40., -AX*(1.-WILT)))
        EXPARG2=EXP(ARG2)
        RZEQW=(WILT-1.-(2./AX))*EXPARG2 + WILT + (2./AX)
        RZEQWI=AX*EXPARG2 *                                                    &
         ( -1. -(2./AX) - (2./(AX*AX)) + WILT + (WILT/AX) )                    &
         + WILT + 2./AX
        AR20=1.+(-AX-1.+AX*WILT)*EXPARG2
        RZEQWI=RZEQWI/(AR20+1.E-20)

!**** SITUATION 1: CATDEF LE CDCR1
        IF(CATDEF(N) .LE. CDCR1(N)) THEN
          RZEQY=RZEQX+WRZ
          RZEQYI=RZEQXI+WRZ
          WMNEW=WMIN+WRZ
          ARG3=AMAX1(-40., AMIN1(40., -AX*(1.-WMNEW)))
          EXPARG3=EXP(ARG3)
          AREA1=(1.+AX-AX*WMIN)*EXPARG1
          AREA2=(1.+AX-AX*WMNEW)*EXPARG3
          IF(WMNEW .GE. WILT) THEN
            AR1(N)=AR1(N)+AREA2-AREA1
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF
          IF(WMNEW .LT. WILT) THEN
            AREA3=(1.+AX-AX*WILT)*EXPARG2
            AR1(N)=AR1(N)+AREA3-AREA1
            AR2(N)=1.-AR1(N)
            FACTOR=(RZEQX+WRZ-WILT)/(RZEQW-WILT)
            AR1(N)=AR1(N)*FACTOR
            AR2(N)=AR2(N)*FACTOR
            AR4(N)=1.-FACTOR
            ENDIF
          ENDIF

!**** SITUATION 2: CATDEF GT CDCR1
        IF(CATDEF(N) .GT. CDCR1(N)) THEN
          FACTOR=(CDCR2(N)-CATDEF(N))/(CDCR2(N)-CDCR1(N))
          RZEQY=WILT+(RZEQX-WILT)*FACTOR+WRZ
          RZEQYI=WILT+(RZEQXI-WILT)*FACTOR+WRZ

          IF(RZEQY .LT. WILT) THEN
            IF(RZEQY .LT. WILT-.001) THEN
!rr                WRITE(*,*) 'RZEXC WAY TOO LOW!  N=',N,' RZEQY=',RZEQY
!rr                WRITE(*,*) 'SRFEXC=',SRFEXC(N),'RZEXC=',RZEXC(N),
!rr     &                     'CATDEF=',CATDEF(N)
!             ELSE
!               WRITE(*,*) 'RZEXC TOO LOW  N=',N
              ENDIF
            RZEQY=WILT
            RZEQYI=WILT
            ENDIF

          IF(RZEQY .GE. RZEQX) THEN  ! RZEXC BRINGS MOISTURE ABOVE CDCR1 POINT
            WMNEW=WMIN+(RZEQY-RZEQX)
            ARG3=AMAX1(-40., AMIN1(40., -AX*(1.-WMNEW)))
            EXPARG3=EXP(ARG3)
            AREA1=(1.+AX-AX*WMIN)*EXPARG1
            AREA2=(1.+AX-AX*WMNEW)*EXPARG3
            AR1(N)=AR1(N)+(AREA2-AREA1)
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF

          IF(RZEQY .LT. RZEQX .AND. RZEQY .GE. RZEQW) THEN
            CATDEFW=CDCR2(N)+((RZEQW-WILT)/(RZEQX-WILT))*(CDCR1(N)-CDCR2(N))
            AR1W= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFW)                       &
                 /(1.+ars2(n)*CATDEFW+ars3(n)*CATDEFW*CATDEFW)))
            FACTOR=(RZEQY-RZEQW)/(RZEQX-RZEQW)
            AR1(N)=AR1W+FACTOR*(AR1(N)-AR1W)
            AR2(N)=1.-AR1(N)
            AR4(N)=0.
            ENDIF

          IF(RZEQY .LT. RZEQW) THEN
            CATDEFW=CDCR2(N)+((RZEQW-WILT)/(RZEQX-WILT))*(CDCR1(N)-CDCR2(N))
            AR1W= AMIN1(1.,AMAX1(0.,(1.+ars1(n)*CATDEFW)                       &
                 /(1.+ars2(n)*CATDEFW+ars3(n)*CATDEFW*CATDEFW)))
            AR1(N)=AR1W
            AR2(N)=1.-AR1(N)
            FACTOR=(RZEQY-WILT)/(RZEQW-WILT)
            AR1(N)=AR1(N)*FACTOR
            AR2(N)=AR2(N)*FACTOR
            AR4(N)=1.-FACTOR
            ENDIF

          ENDIF

        RZI(N)=RZEQYI

        SWSRF1(N)=1.
!mjs: changed .001 temporarily because of large bee.
        SWSRF2(N)=AMIN1(1., AMAX1(0.01, RZEQYI))
        SWSRF4(N)=AMIN1(1., AMAX1(0.01, WILT))

!**** EXTRAPOLATION OF THE SURFACE WETNESSES

! 1st step: surface wetness in the unstressed fraction without considering
!           the surface excess; we just assume an equilibrium profile from 
!           the middle of the root zone to the surface.

        SWSRF2(N)=((SWSRF2(N)**(-BEE(N))) - (.5/PSIS(N)))**(-1./BEE(N))
        SWSRF4(N)=((SWSRF4(N)**(-BEE(N))) - (.5/PSIS(N)))**(-1./BEE(N))

! srfmx is the maximum amount of water that can be added to the surface layer
! The choice of defining SWSRF4 like SWSRF2 needs to be better examined.
        srfmx(n)=ar2(n)*(1.-swsrf2(n))*(dzsf(n)*poros(n))
        srfmx(n)=srfmx(n)+ar4(n)*(1.-swsrf4(n))*(dzsf(n)*poros(n))
!**** For calculation of srfmn, assume surface moisture associated with
!**** AR1 is constantly replenished by water table.
        srfmn(n)=-(ar2(n)*swsrf2(n)+ar4(n)*swsrf4(n))*(dzsf(n)*poros(n))

        if(srfexc(n).gt.srfmx(n)) then
            cor=srfexc(n)-srfmx(n)     !  The correction is here
            srfexc(n)=srfmx(n)
            catdef(n)=catdef(n)-cor
            if(catdef(n).lt.0.) then
              runsrf(n)=runsrf(n)-catdef(n)/dtstep
              catdef(n)=0.
              endif
          else if(srfexc(n).lt.srfmn(n)) then
            cor=srfexc(n)-srfmn(n)
            catdef(n)=catdef(n)-cor
            srfexc(n)=srfmn(n)
          else
            cor=0.
          endif
          
        SWSRF2(N)=SWSRF2(N)+SRFEXC(N)/(dzsf(n)*poros(n)*(1.-ar1(n))+1.e-20)
        SWSRF2(N)=AMIN1(1., AMAX1(1.E-5, SWSRF2(N)))
        swsrf4(n)=swsrf4(n)+srfexc(n)/(dzsf(n)*poros(n)*(1.-ar1(n))+1.e-20)
        SWSRF4(N)=AMIN1(1., AMAX1(1.E-5, SWSRF4(N)))

        IF (AR1(N) .ge. 1.-1.E-5) then
          AR1(N)=1.
          AR2(N)=0.
          AR4(N)=0.
          SWSRF2(N)=1.
          SWSRF4(N)=wilt
          ENDIF

        IF (AR1(N) .LT. 0.) then
!rr          IF(AR1(N) .LT. -1.E-3) WRITE(*,*) 'AR1 TOO LOW: AR1=',AR1(N)
          AR1(N)=0.
          ENDIF
        ar1(n)=amax1(0., amin1(1., ar1(n)))
        ar2(n)=amax1(0., amin1(1., ar2(n)))
        ar4(n)=amax1(0., amin1(1., ar4(n)))
        asum=ar1(n)+ar2(n)+ar4(n)
        if(asum .lt. .9999 .or. asum .gt. 1.0001) then
          write(*,*) 'Areas do not add to 1: sum=',asum,'N=',n
       endif


        ENDDO

         
      RETURN
      END SUBROUTINE PARTITION

!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****

      SUBROUTINE energy1 (                                                     &
                       NCH, DTSTEP, ITYP, UM, RCIN,                            &
                       ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,         &
                       QM,     RA,   SWNET,  HLWDWN, PSUR,                     &
                       RDC,    HFTDS, DHFTDS,                                  &
                       QSATTC, DQSDTC, ALWRAD, BLWRAD,                         &
                       EMAXRT,CSOIL,SWSRF,POTFRC,BUG,                          &
                       TC, QA,                                                 &
                       EVAP, SHFLUX, HLWUP, RX1, RX2, GHFLUX, HSNACC)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: UM, RCIN, ETURB, HSTURB, QM, RA,     &
                    SWNET, HLWDWN, PSUR, RDC, HFTDS, DHFTDS, QSATTC, DQSDTC,   &
                    ALWRAD, BLWRAD, EMAXRT, CSOIL, SWSRF, POTFRC, DEDQA,       &
                    DEDTC, DHSDQA, DHSDTC
      LOGICAL, INTENT(IN) :: BUG


      REAL, INTENT(INOUT), DIMENSION(NCH) :: TC, QA

      REAL, INTENT(OUT), DIMENSION(NCH) :: EVAP, SHFLUX, HLWUP, RX1, RX2,      &
                    GHFLUX, HSNACC


      INTEGER ChNo, N
      REAL, DIMENSION(NCH) :: VPDSTR, ESATTX, VPDSTX, FTEMP, RC, EAX, TX,      &
                    RCX, DRCDTC, DUMMY,  FTEMPX, DRCDEA, DEDEA, DHSDEA, EM,    &
                    ESATTC, DESDTC, EA
      REAL  DELTC, DELEA

!**** - - - - - - - -
!****
      DATA DELTC /0.01/, DELEA /0.001/
!****
!**** --------------------------------------------------------------------- 
!****
!**** Expand data as specified by ITYP into arrays of size NCH.
!****

!****
!**** Pre-process input arrays as necessary:

      DO 100 ChNo = 1, NCH

!      DEDQA(CHNO)  = AMAX1(  DEDQA(CHNO), 500./ALHE )
!      DEDTC(CHNO)  = AMAX1(  DEDTC(CHNO),   0. )
!      DHSDQA(CHNO) = AMAX1( DHSDQA(CHNO),   0. )
!      DHSDTC(CHNO) = AMAX1( DHSDTC(CHNO), -10. )

      EM(CHNO)     = QM(CHNO) * PSUR(CHNO) / EPSILON
      EA(CHNO)     = QA(CHNO) * PSUR(CHNO) / EPSILON
      ESATTC(CHNO) = QSATTC(CHNO) * PSUR(CHNO) / EPSILON
      DESDTC(CHNO) = DQSDTC(CHNO) * PSUR(CHNO) / EPSILON
      DEDEA(CHNO)  = DEDQA(CHNO) * EPSILON / PSUR(CHNO)
      DHSDEA(CHNO) = DHSDQA(CHNO) * EPSILON / PSUR(CHNO)

 100  CONTINUE

!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
!****          ( VPDFAC  computes the vapor pressure deficit stress term,
!****           TMPFAC  computes the temperature stress term,
!****           RSURFP  computes rc given a parallel resist. from the
!****                   surface,
!****           RCANOP  computes rc corrected for snow and interception.)
!****

      CALL VPDFAC (                                                            &
                   NCH,  ITYP,  ESATTC, EA,                                    &
                   VPDSTR                                                      &
                  )


      CALL TMPFAC (                                                            &
                   NCH,  ITYP, TC,                                             &
                   FTEMP                                                       &
                  )


      DO N=1,NCH
        RC(N)=RCIN(N)/(VPDSTR(N)*FTEMP(N)+1.E-20)
        ENDDO

      CALL RSURFP1 (                                                           &
                   NCH, UM, RDC, SWSRF,ESATTC, EA,                             &
                   RC,                                                         &
                   RX1, RX2                                                    &
                  )


      CALL RCANOP (                                                            &
                   NCH, RA, ETURB, POTFRC,                                     &
                   RC                                                          &
                  )

!****
!**** -    -    -    -    -    -    -    -    -    -    -    -    -    -
!****
!**** Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
!****

      DO ChNo = 1, NCH
        TX(ChNo) = TC(ChNo) + DELTC
        ESATTX(ChNo) = ESATTC(ChNo) + DESDTC(CHNO) * DELTC
        EAX(ChNo) = EA(ChNo) + DELEA
        ENDDO

!****
!**** temperature:
      CALL VPDFAC (NCH, ITYP, ESATTX, EA, VPDSTX)
      CALL TMPFAC (NCH, ITYP, TX, FTEMPX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMPX(N)+1.E-20)
        ENDDO

      CALL RSURFP1 (NCH, UM, RDC, SWSRF, ESATTX, EA,                           &
                   RCX,                                                        &
                   DUMMY,DUMMY                                                 &
                  )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
!****
      DO  ChNo = 1, NCH
        DRCDTC(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELTC
        ENDDO
!****

!**** vapor pressure:
      CALL VPDFAC (NCH, ITYP, ESATTC, EAX, VPDSTX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMP(N)+1.E-20)
        ENDDO

      CALL RSURFP1 (NCH, UM, RDC, SWSRF, ESATTC, EAX,                          &
                   RCX,                                                        &
                   DUMMY,DUMMY                                                 &
                  )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
!****
      DO ChNo = 1, NCH
         DRCDEA(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELEA
         ENDDO
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 2: Solve the energy balance at the surface.
!****
      CALL FLUXES (                                                            &
                      NCH,   ITYP, DTSTEP, ESATTC, DESDTC,                     &
                    ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,             &
                       RC, DRCDEA, DRCDTC, SWNET, HLWDWN, ALWRAD, BLWRAD,      &
                       EM,  CSOIL,   PSUR, EMAXRT,  HFTDS, DHFTDS,             &
                       TC,     EA,                                             &
                     EVAP, SHFLUX,  HLWUP, GHFLUX,  HSNACC                     &
                  )

!****
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****
!****
!**** Process data for return to GCM:

      DO 2000 ChNo = 1, NCH
      QA(CHNO) = EA(CHNO) * EPSILON / PSUR(CHNO)
 2000 CONTINUE

!****
      RETURN
      END SUBROUTINE energy1

!****
!**** [ END CHIP ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
      SUBROUTINE energy2 (                                                     &
                       NCH, DTSTEP, ITYP, UM, RCIN,                            &
                       ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,         &
                       QM,     RA,   SWNET,  HLWDWN, PSUR,                     &
                       RDC,    HFTDS, DHFTDS,                                  &
                       QSATTC, DQSDTC, ALWRAD, BLWRAD,                         &
                       EMAXRT,CSOIL,SWSRF,POTFRC,BUG,RZI, WPWET,               &
                       TC, QA,                                                 &
                       EVAP, SHFLUX, HLWUP, RX1, RX2, GHFLUX, HSNACC)


      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) ::  ITYP

      REAL, INTENT(IN) ::  DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) ::  UM, RCIN, ETURB, HSTURB, QM, RA,    &
                    SWNET, HLWDWN, PSUR, RDC, HFTDS, DHFTDS, QSATTC, DQSDTC,   &
                    ALWRAD, BLWRAD, EMAXRT, CSOIL, SWSRF, POTFRC, RZI, WPWET,  &
                    DEDQA, DEDTC, DHSDQA, DHSDTC

      LOGICAL, INTENT(IN) ::   BUG


      REAL, INTENT(INOUT), DIMENSION(NCH) :: TC, QA

      REAL, INTENT(OUT), DIMENSION(NCH) :: EVAP, SHFLUX, HLWUP, RX1, RX2,      &
                    GHFLUX, HSNACC



      INTEGER ChNo, N
      REAL, DIMENSION(NCH) :: VPDSTR, ESATTX, VPDSTX, FTEMP, RC, EAX, TX,      &
                    RCX, DRCDTC, DUMMY, FTEMPX, DRCDEA, DEDEA, DHSDEA, EM,     &
                    ESATTC, DESDTC, EA, RSTFAC
      REAL DELTC, DELEA, STEXP, ATRANS, ASTRFR

!      PARAMETER (ASTRFR=.333)  ! STRESS TRANSITION POINT
!      PARAMETER (STEXP=1.)  ! STRESS RAMPING
! Changed SVK: on 3/21 based on Randy and Sarith's input
!      PARAMETER (ASTRFR=1.)  ! STRESS TRANSITION POINT
!      PARAMETER (STEXP=2.)  ! STRESS RAMPING
      DATA DELTC /0.01/, DELEA /0.001/


      
      ASTRFR = 0.333
      STEXP=2.0
!         ASTRFR = 0.5
!         STEXP=2.0
!         ASTRFR = 0.75
!         STEXP=2.0
!         ASTRFR = 0.667
!         STEXP=1.0
!****
!**** --------------------------------------------------------------------- 
!****
!**** Expand data as specified by ITYP into arrays of size NCH.
!****

!****
!**** Pre-process input arrays as necessary:

      DO 100 ChNo = 1, NCH

!      DEDQA(CHNO)  = AMAX1(  DEDQA(CHNO), 500./ALHE )
!      DEDTC(CHNO)  = AMAX1(  DEDTC(CHNO),   0. )
!      DHSDQA(CHNO) = AMAX1( DHSDQA(CHNO),   0. )
!      DHSDTC(CHNO) = AMAX1( DHSDTC(CHNO), -10. )

      EM(CHNO)     = QM(CHNO) * PSUR(CHNO) / EPSILON
      EA(CHNO)     = QA(CHNO) * PSUR(CHNO) / EPSILON
      ESATTC(CHNO) = QSATTC(CHNO) * PSUR(CHNO) / EPSILON
      DESDTC(CHNO) = DQSDTC(CHNO) * PSUR(CHNO) / EPSILON
      DEDEA(CHNO)  = DEDQA(CHNO) * EPSILON / PSUR(CHNO)
      DHSDEA(CHNO) = DHSDQA(CHNO) * EPSILON / PSUR(CHNO)

 100  CONTINUE

!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
!****          ( VPDFAC  computes the vapor pressure deficit stress term,
!****           TMPFAC  computes the temperature stress term,
!****           RSURFP  computes rc given a parallel resist. from the
!****                   surface,
!****           RCANOP  computes rc corrected for snow and interception.)
!****


!**** Compute water stress effect (RDK 03/30/06)

      CALL VPDFAC (                                                            &
                   NCH,  ITYP,  ESATTC, EA,                                    &
                   VPDSTR                                                      &
                  )


      CALL TMPFAC (                                                            &
                   NCH,  ITYP, TC,                                             &
                   FTEMP                                                       &
                  )

      DO N=1,NCH
        RC(N)=RCIN(N)/(VPDSTR(N)*FTEMP(N)+1.E-20)
     ENDDO
!       if(LIS_localPet.eq.165) then 
!          write(LIS_logunit,*) 'rc1 ',tc(22527), rc(22527)
!       endif


      DO CHNO = 1, NCH
        ATRANS = WPWET(CHNO)+ASTRFR*(1.-WPWET(CHNO))
        RSTFAC(CHNO)=AMAX1( 1.E-3, (RZI(CHNO)-WPWET(CHNO))/                    &
                                         (ATRANS-WPWET(CHNO)) )
        RSTFAC(CHNO)=AMIN1( 1., RSTFAC(CHNO))
        RC(ChNo) = RC(CHNO) / RSTFAC(CHNO)**STEXP
        RC(CHNO) = AMIN1 (RC(CHNO) , 1.E10)
        ENDDO
!       if(LIS_localPet.eq.165) then 
!          write(LIS_logunit,*) 'rc2 ',tc(22527), rc(22527)
!       endif

      CALL RSURFP2 (                                                           &
                   NCH, UM, RDC, SWSRF, ESATTC, EA, WPWET,                     &
                   RC,                                                         &
                   RX1, RX2                                                    &
                  )

!       if(LIS_localPet.eq.165) then 
!          write(LIS_logunit,*) 'rc3 ',tc(22527), rc(22527),potfrc(22527)
!       endif
      CALL RCANOP (                                                            &
                   NCH, RA, ETURB, POTFRC,                                     &
                   RC                                                          &
                  )
!       if(LIS_localPet.eq.165) then 
!          write(LIS_logunit,*) 'rc4 ',tc(22527), rc(22527)
!       endif

!****
!**** -    -    -    -    -    -    -    -    -    -    -    -    -    -
!****
!**** Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
!****

      DO ChNo = 1, NCH
        TX(ChNo) = TC(ChNo) + DELTC
        ESATTX(ChNo) = ESATTC(ChNo) + DESDTC(CHNO) * DELTC
        EAX(ChNo) = EA(ChNo) + DELEA
        ENDDO

!****
!**** temperature:
      CALL VPDFAC (NCH, ITYP, ESATTX, EA, VPDSTX)
      CALL TMPFAC (NCH, ITYP, TX, FTEMPX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMPX(N)+1.E-20)
        RCX(N) = RCX(N) / RSTFAC(N)**STEXP
        RCX(N) = AMIN1 (RCX(N) , 1.E10)
        ENDDO

      CALL RSURFP2 (NCH, UM, RDC, SWSRF,ESATTX, EA, WPWET,                     &
                   RCX,                                                        &
                   DUMMY,DUMMY                                                 &
                  )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
!****
      DO  ChNo = 1, NCH
        DRCDTC(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELTC
        ENDDO
!****

!**** vapor pressure:
      CALL VPDFAC (NCH, ITYP, ESATTC, EAX, VPDSTX)

      DO N=1,NCH
        RCX(N)=RCIN(N)/(VPDSTX(N)*FTEMP(N)+1.E-20)
        RCX(N) = RCX(N) / RSTFAC(N)**STEXP
        RCX(N) = AMIN1 (RCX(N) , 1.E10)
        ENDDO

      CALL RSURFP2 (NCH, UM, RDC, SWSRF, ESATTC, EAX, WPWET,                   &
                   RCX,                                                        &
                   DUMMY,DUMMY                                                 &
                  )
      CALL RCANOP (NCH, RA, ETURB, POTFRC, RCX)
!****
      DO ChNo = 1, NCH
         DRCDEA(ChNo) = (RCX(ChNo) - RC(ChNo)) / DELEA
         ENDDO
!****
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 2: Solve the energy balance at the surface.
!****
      CALL FLUXES (                                                            &
                      NCH,   ITYP, DTSTEP, ESATTC, DESDTC,                     &
                    ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,             &
                       RC, DRCDEA, DRCDTC,SWNET, HLWDWN, ALWRAD, BLWRAD,       &
                       EM,  CSOIL,   PSUR, EMAXRT, HFTDS, DHFTDS,              &
                       TC,     EA,                                             &
                     EVAP, SHFLUX,  HLWUP, GHFLUX, HSNACC                      &
                  )
!       if(LIS_localPet.eq.165) then 
!          write(LIS_logunit,*) 'e2 ',tc(22527)
!       endif

!****
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****
!****
!**** Process data for return to GCM:

      DO 2000 ChNo = 1, NCH
      QA(CHNO) = EA(CHNO) * EPSILON / PSUR(CHNO)
 2000 CONTINUE
!****

      RETURN
      END SUBROUTINE energy2

!****
!**** [ END CHIP ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
      SUBROUTINE energy4 (                                                     &
                       NCH, DTSTEP, ITYP, UM, RCIN,                            &
                       ETURB,  DEDQA,  DEDTC,  HSTURB, DHSDQA, DHSDTC,         &
                       QM,     RA,   SWNET,  HLWDWN, PSUR,                     &
                       RDC,    HFTDS, DHFTDS,                                  &
                       QSATTC, DQSDTC, ALWRAD, BLWRAD,                         &
                       EMAXRT,CSOIL,SWSRF,POTFRC,BUG,WPWET,                    &
                       TC, QA,                                                 &
                       EVAP, SHFLUX, HLWUP, RX1, RX2, GHFLUX, HSNACC           &
                       )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: UM, RCIN, ETURB, HSTURB, QM, RA,     &
                SWNET, HLWDWN, PSUR, RDC, HFTDS, DHFTDS, QSATTC, DQSDTC,       &
                ALWRAD, BLWRAD, EMAXRT, CSOIL, SWSRF, POTFRC, WPWET, DEDQA,    &
                DEDTC, DHSDQA, DHSDTC
      LOGICAL, INTENT(IN) ::  BUG

      REAL, INTENT(INOUT), DIMENSION(NCH) :: TC, QA

      REAL, INTENT(OUT), DIMENSION(NCH) :: EVAP, SHFLUX, HLWUP, RX1, RX2,      &
                GHFLUX, HSNACC


      INTEGER ChNo, N
      REAL, DIMENSION(NCH) :: DEDEA, DHSDEA, EM, ESATTC, DESDTC, EA, RC,       &
                DRCDTC, DRCDEA
      REAL  DELTC, DELEA

!****
      DATA DELTC /0.01/, DELEA /0.001/
!****
!**** --------------------------------------------------------------------- 
!****
!**** Expand data as specified by ITYP into arrays of size NCH.
!****

!****
!**** Pre-process input arrays as necessary:

      DO 100 ChNo = 1, NCH

!      DEDQA(CHNO)  = AMAX1(  DEDQA(CHNO), 500./ALHE )
!      DEDTC(CHNO)  = AMAX1(  DEDTC(CHNO),   0. )
!      DHSDQA(CHNO) = AMAX1( DHSDQA(CHNO),   0. )
!      DHSDTC(CHNO) = AMAX1( DHSDTC(CHNO), -10. )

      EM(CHNO)     = QM(CHNO) * PSUR(CHNO) / EPSILON
      EA(CHNO)     = QA(CHNO) * PSUR(CHNO) / EPSILON
      ESATTC(CHNO) = QSATTC(CHNO) * PSUR(CHNO) / EPSILON
      DESDTC(CHNO) = DQSDTC(CHNO) * PSUR(CHNO) / EPSILON
      DEDEA(CHNO)  = DEDQA(CHNO) * EPSILON / PSUR(CHNO)
      DHSDEA(CHNO) = DHSDQA(CHNO) * EPSILON / PSUR(CHNO)

 100  CONTINUE

!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 1: COMPUTE EFFECTIVE RESISTANCE RC FOR ENERGY BALANCE.
!****          ( VPDFAC  computes the vapor pressure deficit stress term,
!****           TMPFAC  computes the temperature stress term,
!****           RSURFP  computes rc given a parallel resist. from the
!****                   surface,
!****           RCANOP  computes rc corrected for snow and interception.)
!****

      DO N=1,NCH
        RC(N)=RCIN(N)
        ENDDO

      CALL RSURFP2 (                                                           &
                   NCH, UM, RDC, SWSRF, ESATTC, EA, WPWET,                     &
                   RC,                                                         &
                   RX1, RX2                                                    &
                  )

      CALL RCANOP (                                                            &
                   NCH, RA, ETURB, POTFRC,                                     &
                   RC                                                          &
                  )

!****
!**** -    -    -    -    -    -    -    -    -    -    -    -    -    -
!****
!**** Compute DRC/DT and DRC/DEA using temperature, v.p. perturbations:
!****

      DO ChNo = 1, NCH
        DRCDTC(CHNO)=0.
        DRCDEA(CHNO)=0.
        ENDDO

!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****       STEP 2: Solve the energy balance at the surface.
!****
      CALL FLUXES (                                                            &
                      NCH,   ITYP, DTSTEP, ESATTC, DESDTC,                     &
                    ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,             &
                       RC, DRCDEA, DRCDTC,SWNET, HLWDWN, ALWRAD, BLWRAD,       &
                       EM,  CSOIL,   PSUR, EMAXRT, HFTDS, DHFTDS,              &
                       TC,     EA,                                             &
                     EVAP, SHFLUX,  HLWUP, GHFLUX, HSNACC                      &
                   )

!****
!****
!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****
!****
!**** Process data for return to GCM:

      DO 2000 ChNo = 1, NCH
      QA(CHNO) = EA(CHNO) * EPSILON / PSUR(CHNO)
 2000 CONTINUE
!****

      RETURN
      END SUBROUTINE energy4

!****
!**** [ END CHIP ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!**** [ BEGIN FLUXES ]
!****
      SUBROUTINE FLUXES (                                                      &
                           NCH,   ITYP, DTSTEP, ESATTC, DESDTC,                &
                          ETURB,  DEDEA,  DEDTC, HSTURB, DHSDEA, DHSDTC,       &
                          RC, DRCDEA, DRCDTC,SWNET, HLWDWN, ALWRAD, BLWRAD,    &
                             EM,  CSOIL,   PSUR, EMAXRT,                       &
                          HFTDS, DHFTDS,                                       &
                             TC,     EA,                                       &
                           EVAP, SHFLUX,  HLWUP, GHFLUX, HSNACC                &
                         )
!****
!**** This subroutine computes the fluxes of latent and sensible heat
!**** from the surface through an energy balance calculation.
!****
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP

      REAL, INTENT(IN) :: DTSTEP
      REAL, INTENT(IN), DIMENSION(NCH) :: ESATTC, DESDTC, ETURB,  DEDEA,       &
                 DEDTC, HSTURB, DHSDEA, DHSDTC, RC, DRCDEA, DRCDTC, SWNET,     &
                 HLWDWN, ALWRAD, BLWRAD, EM,  CSOIL, PSUR, EMAXRT, HFTDS,      &
                 DHFTDS

      REAL, INTENT(INOUT), DIMENSION(NCH) :: TC, EA

      REAL, INTENT(OUT), DIMENSION(NCH) :: EVAP, SHFLUX, HLWUP, GHFLUX,        &
           HSNACC


      INTEGER ChNo
      REAL HLWTC,  CDEEPS, Q0,  RHOAIR,   CONST,  DHLWTC, EPLANT, A11, A12,    &
               A21, A22, F0, DEA, DTC, EANEW,  ESATNW,  EHARMN, DETERM, DENOM, &
               EDIF
      LOGICAL DEBUG, CHOKE


      DATA DEBUG /.FALSE./

!****
!**** -------------------------------------------------------------------

      HSNACC = 0.0  ! Initialize INTENT OUT variable before use -- LLT 

      DO 200 ChNo = 1, NCH
!****
      HLWTC = ALWRAD(CHNO) + BLWRAD(CHNO) * TC(CHNO)
      RHOAIR = PSUR(ChNo) * 100. / (RGAS * TC(ChNo))
      CONST = RHOAIR * EPSILON / PSUR(ChNo)
      DHLWTC = BLWRAD(CHNO)
!****
!**** Compute matrix elements A11, A22, AND Q0 (energy balance equation).
!****
      A11 = CSOIL(ChNo)/DTSTEP +                                               &
              DHLWTC +                                                         &
              DHSDTC(ChNo) +                                                   &
              ALHE*DEDTC(ChNo) +                                               &
              DHFTDS(CHNO)
      A12 = DHSDEA(ChNo) + ALHE * DEDEA(ChNo)
      Q0 =  SWNET(ChNo) +                                                      &
              HLWDWN(ChNo) -                                                   &
              HLWTC -                                                          &
              HSTURB(ChNo) -                                                   &
              ALHE * ETURB(ChNo) -                                             &
              HFTDS(CHNO)
!****
!**** Compute matrix elements A21, A22, and F0 (canopy water budget  
!**** equation) and solve for fluxes.  Three cases are considered:
!****
!**** 1. Standard case: RC>0.
!**** 2. RC = 0.  Can only occur if CIR is full or ETURB is negative.
!****
      CHOKE = .TRUE.
!       if(LIS_localPet.eq.165.and.chno.eq.22527) then 
!          write(LIS_logunit,*) 'fx ',tc(chno), rc(chno)
!       endif
      IF( RC(CHNO) .GT. 0.) THEN 
          EPLANT = CONST * (ESATTC(ChNo) - EA(ChNo)) / RC(ChNo)
          IF(EPLANT*ETURB(ChNo).GT.0.) THEN
              EHARMN = 2.*EPLANT*ETURB(CHNO) / (EPLANT + ETURB(ChNo))
            ELSE
              EHARMN=0.
            ENDIF
!****
!****            Some limitations to A21 and A22 are applied:
!****            we assume that the increase in plant evaporation
!****            due to an increase in either TC or EA balances 
!****            or outweighs any decrease due to RC changes.
!****

          A21 =  -DEDTC(ChNo)*RC(ChNo) +                                       &
            amax1(0., CONST*DESDTC(ChNo) - EHARMN*DRCDTC(ChNo) )
          A22 = -( RC(ChNo)*DEDEA(ChNo) +                                      &
                     amax1( 0., CONST + EHARMN*DRCDEA(ChNo) )   )

          F0 = RC(ChNo) * (ETURB(ChNo) - EPLANT)
          DETERM = AMIN1( A12*A21/(A11*A22) - 1., -0.1 )
          DEA = ( Q0*A21 - A11*F0 ) / ( DETERM * A11*A22 )
          DTC = ( Q0 - A12*DEA ) / A11
!       if(LIS_localPet.eq.165.and.ChNo.eq.22527) then 
!          write(LIS_logunit,*) 'f0 ',Q0, A12, DEA, A11, dtc
!       endif
          EVAP(ChNo) = ETURB(ChNo) + DEDEA(ChNo)*DEA + DEDTC(ChNo)*DTC
          SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC
          DENOM = DETERM * A11*A22
        ELSE
          CHOKE = .FALSE.
          A21 = -DESDTC(ChNo)
          A22 = 1.
          F0 = ESATTC(ChNo) - EA(ChNo)
          DEA = ( Q0*A21 - A11*F0 ) / ( A12*A21 - A11*A22 )
          DTC = ( Q0 - A12*DEA ) / A11
!       if(LIS_localPet.eq.165.and.ChNo.eq.22527) then 
!          write(LIS_logunit,*) 'f1 ',Q0, A12, DEA, A11, dtc
!       endif

          EVAP(ChNo) = ETURB(ChNo) + DEDEA(ChNo)*DEA + DEDTC(ChNo)*DTC
          SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC
          DENOM = A12 * A21 - A11*A22
        ENDIF

!**** - - - - - - - - - - - - - - - - - - - - - - -
!**** Adjustments

!**** 1. Adjust deltas and fluxes if all available water evaporates
!****    during time step:
!**** NOTE: SOME CALCS BELOW ASSUME CROSS DERIVATIVE TERMS ARE ZERO
!****
      IF( EVAP(CHNO) .GT. EMAXRT(CHNO) ) THEN
        CHOKE = .FALSE.
        DEA=SIGN(ETURB(CHNO),EA(CHNO))
        IF(DEDEA(CHNO) .NE. 0.) DEA = (EMAXRT(CHNO)-ETURB(CHNO))/DEDEA(CHNO) 
        DEA = EM(CHNO) - EA(CHNO)
        DTC =  (Q0 + ALHE*(ETURB(ChNo)-EMAXRT(CHNO)) - DHSDEA(CHNO)*DEA)       &
                /  ( A11 - ALHE*DEDTC(ChNo) )
!       if(LIS_localPet.eq.165.and.ChNo.eq.22527) then 
!          write(LIS_logunit,*) 'f2 ',dtc
!       endif
        EVAP(CHNO) = EMAXRT(CHNO)
        SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC
        ENDIF


!**** 2. Pathological cases. 

!**** CHECK THAT EVAP AND ETURB ARE STILL CONSISTENT
!**** NOTE: SOME CALCS BELOW ASSUME CROSS DERIVATIVE TERMS ARE ZERO

      IF( EVAP(CHNO)*(ETURB(CHNO)+DEDEA(CHNO)*DEA) .LT. 0. )  THEN 
        CHOKE = .FALSE.
        DEA=SIGN(ETURB(CHNO),EA(CHNO))
        IF(DEDEA(CHNO) .NE. 0.) DEA = -ETURB(CHNO)/DEDEA(CHNO) 
        DTC = ( Q0 + ALHE*ETURB(ChNo) - DHSDEA(CHNO)*DEA ) /                   &
                  ( A11 - ALHE*DEDTC(ChNo) )
!       if(LIS_localPet.eq.165.and.ChNo.eq.22527) then 
!          write(LIS_logunit,*) 'f3 ',tc(ChNo), dtc
!       endif
        EVAP(CHNO) = 0.
        SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC
        ENDIF



!**** 3. Excessive dea change: apply "choke".
! -- (03.09.98) : changed to correct the conservation.
!        if(LIS_localPet.eq.165.and.ChNo.eq.22527) then 
!           write(LIS_logunit,*) 'dea ',abs(dea), 0.5*EA(CHNO)
!        endif
      IF( CHOKE .AND. ABS(DEA) .GT. 0.5*EA(CHNO) ) THEN
        DEA = SIGN(.5*EA(CHNO),DEA)
        DTC = ( Q0 - A12*DEA ) / A11
!       if(LIS_localPet.eq.165.and.ChNo.eq.22527) then 
!          write(LIS_logunit,*) 'f4 ',tc(ChNo), dtc
!       endif
        EVAP(ChNo)   = ETURB(ChNo) + DEDEA(ChNo)*DEA + DEDTC(ChNo)*DTC
        SHFLUX(ChNo) = HSTURB(ChNo) + DHSDEA(ChNo)*DEA + DHSDTC(ChNo)*DTC

        IF(EVAP(CHNO) .GT. EMAXRT(CHNO)) THEN
          EDIF=EVAP(CHNO)-EMAXRT(CHNO)
          EVAP(CHNO) = EMAXRT(CHNO)
          SHFLUX(ChNo) = SHFLUX(CHNO)+EDIF*ALHE
!          HSNACC(ChNo) = HSNACC(CHNO)+EDIF*ALHE
          ENDIF

       ENDIF

!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!       if(LIS_localPet.eq.165.and.ChNo.eq.22527) then 
!          write(LIS_logunit,*) 'f5 ',tc(ChNo), dtc, ( Q0 - A12*DEA ) / A11
!       endif
      TC(ChNo) = TC(ChNo) + DTC
      EA(ChNo) = EA(ChNo) + DEA
      HLWUP(CHNO) = HLWTC + DHLWTC*DTC

      HLWTC = ALWRAD(CHNO) + BLWRAD(CHNO) * TC(CHNO)
! warning: this ghflux is the real ground heat flux, and does not include
! the temperature variation
      GHFLUX(CHNO)=HFTDS(CHNO)+DHFTDS(CHNO)*DTC

!**** Make sure EA remains positive

      EA(CHNO) = AMAX1(EA(CHNO), 0.0)

  200 CONTINUE

      RETURN
      END SUBROUTINE FLUXES
!****
!**** [ END FLUXES ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN VPDFAC ]
!****
      SUBROUTINE VPDFAC (                                                      &
                         NCH, ITYP, ESATTC, EA,                                &
                         VPDSTR                                                &
                        )
!****
!**** This subroutine computes the vapor pressure deficit stress.
!****
      IMPLICIT NONE

      INTEGER, INTENT(IN):: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP
      REAL, INTENT(IN), DIMENSION(NCH) :: ESATTC, EA
      REAL, INTENT(OUT), DIMENSION(NCH) :: VPDSTR

      INTEGER :: ChNo
      REAL, DIMENSION(NTYPS) :: VGDFAC
!****
      DATA VGDFAC /   .0273,    .0357,    .0310,    .0238,                     &
                      .0275,    .0275,       0.,       0.,                     &
                         0.,       0. /
!****
!**** -----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
!****
!      VPDSTR(ChNo) = 1. - (ESATTC(ChNo)-EA(ChNo)) * VGDFAC(ITYP(ChNo))
!      VPDSTR (ChNo) = AMIN1( 1., AMAX1( VPDSTR(ChNo), 1.E-10 ) )
      VPDSTR(CHNO) = 1.
!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE VPDFAC
!****
!**** [ END VPDFAC ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN TMPFAC ]
!****
      SUBROUTINE TMPFAC (                                                      &
                         NCH,  ITYP, TC,                                       &
                         FTEMP                                                 &
                        )
!****
!**** Compute temperature stress factor.
!****
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP
      REAL, INTENT(IN), DIMENSION(NCH) :: TC
      REAL, INTENT(OUT), DIMENSION(NCH) :: FTEMP

      INTEGER ChNo, TypPtr
      INTEGER, PARAMETER :: MEMFAC = 5
      REAL, DIMENSION(MEMFAC*NTYPS) :: VGTLL, VGTU, VGTCF1, VGTCF2, VGTCF3

      DATA VGTLL /MemFac*273., MemFac*273., MemFac*268., MemFac*283.,          &
                  MemFac*283., MemFac*273., MemFac*  0., MemFac*  0.,          &
                  MemFac*  0., MemFac*  0. /
      DATA VGTU /MemFac*318., MemFac*318., MemFac*313., MemFac*328.,           &
                 MemFac*323., MemFac*323., MemFac*  0., MemFac*  0.,           &
                 MemFac*  0., MemFac*  0. /
      DATA VGTCF1 / MemFac*-1.43549E-06,  MemFac*-6.83584E-07,                 &
                    MemFac* 1.67699E-07,  MemFac*-1.43465E-06,                 &
                    MemFac*-2.76097E-06,  MemFac*-1.58094E-07,                 &
                             MemFac* 0.,           MemFac* 0.,                 &
                             MemFac* 0.,           MemFac* 0. /
      DATA VGTCF2 / MemFac* 7.95859E-04,  MemFac* 3.72064E-04,                 &
                    MemFac*-7.65944E-05,  MemFac* 8.24060E-04,                 &
                    MemFac* 1.57617E-03,  MemFac* 8.44847E-05,                 &
                             MemFac* 0.,           MemFac* 0.,                 &
                             MemFac* 0. ,          MemFac* 0./
      DATA VGTCF3 / MemFac*-1.11575E-01,  MemFac*-5.21533E-02,                 &
                    MemFac* 6.14960E-03,  MemFac*-1.19602E-01,                 &
                    MemFac*-2.26109E-01,  MemFac*-1.27272E-02,                 &
                             MemFac* 0.,           MemFac* 0.,                 &
                             MemFac* 0.,           MemFac* 0. /
!****
!**** ----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
!****
      TypPtr = MOD(ChNo,MemFac) + (ITYP(ChNo)-1)*MemFac + 1
      FTEMP(ChNo) = (TC(ChNo) - VGTLL(TypPtr)) * (TC(ChNo) - VGTU(TypPtr)) *   &
                          ( VGTCF1(TypPtr)*TC(ChNo)*TC(ChNo) +                 &
                            VGTCF2(TypPtr)*TC(ChNo) + VGTCF3(TypPtr) )
      IF ( TC(ChNo) .LE. VGTLL(TypPtr) .OR. TC(ChNo) .GE. VGTU(TypPtr) )       &
            FTEMP (ChNo) = 1.E-10
      FTEMP(CHNO) = AMIN1( 1., AMAX1( FTEMP(ChNo), 1.E-10 ) )
!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE TMPFAC
!****
!**** [ END TMPFAC ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN RSURFP ]
!****
      SUBROUTINE RSURFP1 (                                                     &
                         NCH, UM, RDC, WET, ESATTC, EA,                        &
                         RC,                                                   &
                         RX1, RX2                                              &
                         )
!****
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN), DIMENSION(NCH) :: UM, RDC, WET, ESATTC, EA
      REAL, INTENT(INOUT), DIMENSION(NCH) :: RC
      REAL, INTENT(OUT), DIMENSION(NCH) :: RX1, RX2
! Changed SVK: on 3/21 based on Randy and Sarith's input
      real, parameter          :: RSSAT = 300.
      INTEGER ChNo
      REAL  U2, RSURF, HESAT
      
!****
!**** -----------------------------------------------------------------

      DO 100 ChNo = 1, NCH
!****
      U2 = UM(ChNo)

      RSURF = RDC(ChNo) / U2

      RX1(CHNO)=RC(CHNO)
      RX2(CHNO)=RSURF

      RC(ChNo) = RC(CHNO) * RSURF / ( RC(ChNo) + RSURF )
!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE RSURFP1
!****
!**** [ END RSURFP ]
!****
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN RSURFP ]
!****
      SUBROUTINE RSURFP2 (                                                     &
                         NCH, UM, RDC, WET, ESATTC, EA, WPWET,                 &
                         RC,                                                   &
                         RX1, RX2                                              &
                        )
!****
      IMPLICIT NONE


      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN), DIMENSION(NCH) :: UM, RDC, WET, ESATTC, EA, WPWET
      REAL, INTENT(INOUT), DIMENSION(NCH) :: RC
      REAL, INTENT(OUT), DIMENSION(NCH) :: RX1, RX2

      INTEGER ChNo
      REAL U2, RSURF, HESAT, RSWILT, RSSAT, ATERM, BTERM

! RDK 04/04/06
!  VALUES OF BARE SOIL SURFACE RESISTANCE AT WILTING POINT, SATURATION
! Changed SVK: on 3/21 based on Randy and Sarith's input
!      PARAMETER (RSWILT=500., RSSAT=25.) 
!      PARAMETER (RSWILT=2000., RSSAT=300.)
!****
!**** -----------------------------------------------------------------
      RSWILT=500.0
      RSSAT = 25.0
      
      DO 100 ChNo = 1, NCH
!****
      U2 = UM(ChNo)
      BTERM=(RSWILT-RSSAT) / (1./WPWET(CHNO)**2 -1.)
      ATERM=RSSAT-BTERM
      RSURF = RDC(ChNo) / U2 + ATERM + BTERM / (1.E-10 + WET(ChNo))**2

!**** Account for subsaturated humidity at soil surface:
!****
      HESAT = ESATTC(CHNO) * MIN( 1., WET(CHNO)*2. )
      IF( EA(CHNO) .LT. HESAT ) THEN
          RSURF=RSURF*( 1. + (ESATTC(CHNO)-HESAT)/(HESAT-EA(CHNO)) )
        ELSE
          RSURF=1.E10
        ENDIF


      RX1(CHNO)=RC(CHNO)
      RX2(CHNO)=RSURF

      RC(ChNo) = RC(CHNO) * RSURF / ( RC(ChNo) + RSURF )
!****
 100  CONTINUE
!****
      RETURN
      END SUBROUTINE RSURFP2
!****
!**** [ END RSURFP ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------
!****
!**** [ BEGIN RCANOP ]
!****
      SUBROUTINE RCANOP (                                                      &
                         NCH, RA, ETURB, POTFRC,                               &
                         RC                                                    &
                        )
!****
!**** The effective latent heat resistance RC depends on the quantity 
!**** of interception reservoir water and the snow cover.  POTFRC
!**** is the fraction of the tile from which potential evaporation
!**** occurs.
!****
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NCH
      REAL, INTENT(IN), DIMENSION(NCH) :: RA, ETURB, POTFRC
      REAL, INTENT(INOUT), DIMENSION(NCH) :: RC

      INTEGER N
!      REAL ETCRIT,RAMPFC

!**** (Note: ETCRIT arbitrarily set to ~-5 W/m2, or -2.e-6 mm/sec.)
!      DATA ETCRIT/ -2.E-6 /
!****
!**** -----------------------------------------------------------------

      DO N = 1, NCH

        RC(N)=RC(N)*(1.-POTFRC(N)) / ( 1.+POTFRC(N)*RC(N)/RA(N) )

!**** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!****   Assume RC=0 for condensation (dew).
!****   RAMPFC is used to ensure continuity in RC.

!_RDK   Remove zeroing of resistance for dew, due to stability problems
!        RAMPFC=ETURB(N)/ETCRIT
!        IF ( RAMPFC .GE. 0. ) RC(N) = RC(N)*(1.-RAMPFC)
!        IF ( RAMPFC .GT. 1. ) RC(N) = 0.
!****
        ENDDO

      RETURN
      END SUBROUTINE RCANOP
!****
!**** [ END RCANOP ]
!****
!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------

      SUBROUTINE RZEQUIL (                                                     &
                          NCH,ITYP,CATDEF,VGWMAX,CDCR1,CDCR2,WPWET,            &
                          ars1,ars2,ars3,ara1,ara2,ara3,ara4,                  &
                          arw1,arw2,arw3,arw4,                                 &
                          RZEQ                                                 &
                         )

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) :: ITYP
      REAL, INTENT(IN), DIMENSION(NCH) :: CATDEF, VGWMAX, CDCR1, CDCR2,        &
                   WPWET, ars1, ars2, ars3, ara1, ara2, ara3, ara4, arw1,      &
                   arw2, arw3, arw4

      REAL, INTENT(OUT), DIMENSION(NCH) :: RZEQ

      INTEGER N
      REAL AX,WMIN,ASCALE,cdi,wilt,catdefx,factor,ARG1,EXPARG1

! ----------------------------------------------------------------------

      DO N=1,NCH

        WILT=WPWET(N)
        CATDEFX=AMIN1( CATDEF(N) , CDCR1(N) )

! CDI DEFINES IF THE SHAPE PARAMETER IS ADJUSTED IN ONE OR TWO SEGMENTS
        if (ara1(n) .ne. ara3(n)) then
            cdi=(ara4(n)-ara2(n))/(ara1(n)-ara3(n))
          else
            cdi=0.
          endif

        if (CATDEFX .ge. cdi) then
            ax=ara3(n)*CATDEFX+ara4(n)
          else
            ax=ara1(n)*CATDEFX+ara2(n)
          endif

        WMIN=AMIN1(1.,AMAX1(0.,arw4(n)+(1.-arw4(n))*(1.+arw1(n)*CATDEFX)       &
                 /(1.+arw2(n)*CATDEFX+arw3(n)*CATDEFX*CATDEFX)))

        ARG1=AMAX1(-40., AMIN1(40., -AX*(1.-WMIN)))
        EXPARG1=EXP(ARG1)
        RZEQ(N)=(WMIN-1.-(2./AX))*EXPARG1 + WMIN + (2./AX)

        IF(CATDEF(N) .GT. CDCR1(N)) THEN
          FACTOR=(CDCR2(N)-CATDEF(N))/(CDCR2(N)-CDCR1(N))
          RZEQ(N)=WILT+(RZEQ(N)-WILT)*FACTOR
          ENDIF

! scaling:    
        RZEQ(N)=AMIN1(1.,AMAX1(0.,RZEQ(N)))
        RZEQ(N)=RZEQ(N)*VGWMAX(N)

      ENDDO

      RETURN
      END SUBROUTINE RZEQUIL

!**** -----------------------------------------------------------------
!**** /////////////////////////////////////////////////////////////////
!**** -----------------------------------------------------------------


      subroutine gndtp0(t1,phi,zbar,thetaf,ht,fh21w,fh21i,fh21d,                   &
                      dfh21w,dfh21i,dfh21d,tp)

! using a diffusion equation this code generates ground temperatures
! with depth given t1
!            *****************************************
!        input
!        dts     timestep in seconds
!        t1      terrestrial (layer 1) temperature in deg C
!        phi     porosity
!        zbar    mean depth to the water table. 
!        thetaf  mean vadose zone soil moisture factor (0-1) 
!        output,
!        ht      heat content in layers 2-7
!        tp      ground temperatures in layers 2-7 
!        tdeep   the temperature of the "deep"
!        f21     heat flux between layer 2 and the terrestrial layer (1)
!        df21    derivative of f21 with respect to temperature
!        xfice   a total soil column ice factor (0-1)
!             ***********************************

      REAL, INTENT(IN) :: phi, ZBAR, THETAF
      REAL, INTENT(IN), DIMENSION(*) :: HT
      REAL, INTENT(IN), DIMENSION(N_SM) :: T1

      REAL, INTENT(OUT) :: FH21W, FH21I, FH21D, DFH21W, DFH21I, DFH21D
      REAL, INTENT(OUT), DIMENSION(*) :: TP

      INTEGER L, K
      REAL, DIMENSION(N_GT) :: FICE, SHC, ZC, XKLH
      REAL, DIMENSION(N_GT+1) :: FH, ZB
      REAL SHW0, SHI0, SHR0, WS, XW, A1, TK1, A2, TK2, TK3, TKSAT,         &
           XWI, XD1, XD2, DENOM, XKLHW, TKDRY

      !data dz/0.0988,0.1952,0.3859,0.7626,1.5071,10.0/
      !DATA PHI/0.45/, FSN/3.34e+8/, SHR/2.4E6/

!     initialize parameters
      shw0=SHW*1000. ! PER M RATHER THAN PER KG 
      shi0=SHI*1000. ! PER M RATHER THAN PER KG 
      shr0=SHR*1000. ! PER M RATHER THAN PER KG 

! calculate the boundaries, based on the layer thicknesses(DZGT)

      zb(1)=-0.05
      zb(2)=zb(1)-DZGT(1)
      shc(1)=shr0*(1.-phi)*DZGT(1)
      zc(1)=0.5*(zb(1)+zb(2))

! evaluates the temperatures in the soil layers based on the heat values.  
!             ***********************************
!             input:
!             xw - water in soil layers, m
!             ht - heat in soil layers
!             fsn - heat of fusion of water   J/m
!             shc - specific heat capacity of soil
!             shi - specific heat capacity of ice
!             shw - specific heat capcity of water
!             snowd - snow depth, equivalent water m
!             output:
!             tp - temperature of layers, c
!             fice - fraction of ice of layers
!             pre - extra precipitation, i.e. snowmelt, m s-1
!             snowd - snow depth after melting, equivalent water m.
!             ***********************************
! determine fraction of ice and temp of soil layers based on layer
! heat and water content

      ws=phi*DZGT(1)  ! PORE SPACE IN LAYER 2
      xw=0.5*ws     ! ASSUME FOR THESE CALCULATIONS THAT THE PORE SPACE
                    ! IS ALWAYS HALF FILLED WITH WATER.  XW IS THE  
                    ! AMOUNT OF WATER IN THE LAYER.

      tp(1)=0.
      FICE(1) = AMAX1( 0., AMIN1( 1., -ht(1)/(fsn*xw) ) )

      IF(FICE(1) .EQ. 1.) THEN
          tp(1)=(ht(1)+xw*fsn)/(shc(1)+xw*shi0)
        ELSEIF(FICE(1) .EQ. 0.) THEN
          tp(1)=ht(1)/(shc(1)+xw*shw0)
        ELSE
          TP(1)=0.
        ENDIF
           
! evaluates:  layer thermal conductivities
! *****************************************
!             from farouki(cold regions sci and tech, 5, 1981,
!             67-75) the values for tk1,tk2,tk3 are as follows:
!             tk2=2.2**(phi(l,ibv)-xw), tk3=.57**xw, and
!             tk1=3**(1-phi(l,ibv) for %sand<.5, and
!             for computation purposes i have fit these eqs.
!             to 2nd order polynomials.
!             ***********************************
!             input:
!             sklh - soil heat conductivities of layers
!             zb - soil layer boundaries, m
!             zc - soil layer centers, m
!             DZGT - layer thickness, m
!             w - soil water content, m
!             phi - soil porosity, dimensionless
!             q - % sand, silt, clay, peat
!             fice - fraction of ice in layers
!             output:
!             xklh - thermal conductivity, w m-2 k-1
!             ***********************************
! lets get the thermal conductivity for the layers

      a1=1-phi
      tk1=1.01692+a1*(0.89865+1.06211*a1)
      xw=phi*(1.-fice(1))
      a2=phi-xw
      tk2=1.00543+a2*(0.723371+.464342*a2)
      tk3=0.998899+xw*(-0.548043+0.120291*xw)
      tksat=tk1*tk2*tk3

      xwi=1.0
      if (zbar .le. zb(2))then
            xwi=thetaf
         elseif (zbar .ge. zb(2) .and. zbar .le. zb(1))then
            xd1=zb(1)-zbar
            xd2=zbar-zb(2)
            xwi=((xd1*thetaf)+xd2)/(xd1+xd2)
         endif 

      xwi=min(xwi,1.)
      tkdry=0.226 ! = .039*0.45^(-2.2), from Farouki, p. 71
      xklh(1)=(tksat-tkdry)*xwi + tkdry
      xklhw=tksat

      denom=-.025-zc(1)
      fh21w=-xklhw*(t1(1)-273.16-tp(1))/denom
      fh21i=-xklh(1)*(t1(2)-273.16-tp(1))/denom
      fh21d=-xklh(1)*(t1(3)-273.16-tp(1))/denom
      dfh21w=-xklhw/denom
      dfh21i=-xklh(1)/denom
      dfh21d=dfh21i


      return
      end subroutine gndtp0





      subroutine gndtmp(dts,phi,zbar,thetaf,fh21,ht,xfice,tp, FICE)
! using a diffusion equation this code generates ground temperatures
! with depth given t1
!            *****************************************
!        input
!        dts     timestep in seconds
!        phi     porosity
!        t1      terrestrial (layer 1) temperature in deg C
!        zbar    mean depth to the water table. 
!        thetaf  mean vadose zone soil moisture factor (0-1) 
!        output,
!        ht      heat content in layers 2-7
!        tp      ground temperatures in layers 2-7 
!        tdeep   the temperature of the "deep"
!        f21     heat flux between layer 2 and the terrestrial layer (1)
!        df21    derivative of f21 with respect to temperature
!        xfice   a total soil column ice factor (0-1)
!             ***********************************

      REAL, INTENT(IN) :: DTS, phi, ZBAR, THETAF, FH21

      REAL, INTENT(INOUT), DIMENSION(*) :: HT

      REAL, INTENT(OUT) :: XFICE
      REAL, INTENT(OUT), DIMENSION(*) :: TP, FICE

      INTEGER L, LSTART, K
      REAL, DIMENSION(N_GT) :: ZC, SHC, XKLH
      REAL, DIMENSION(N_GT+1) :: FH, ZB
      REAL SHW0, SHI0, SHR0, WS, XW, A1, TK1, A2, TK2, TK3, TKSAT,      &
           XWI, XD1, XD2, XKTH, TKDRY

      !data dz/0.0988,0.1952,0.3859,0.7626,1.5071,10.0/
      !DATA PHI/0.45/, FSN/3.34e+8/, SHR/2.4E6/

! initialize parameters

      shw0=SHW*1000. ! PER M RATHER THAN PER KG 
      shi0=SHI*1000. ! PER M RATHER THAN PER KG 
      shr0=SHR*1000. ! PER M RATHER THAN PER KG 

!----------------------------------
! initialize fice in ALL components
! reichle, July 8, 2002

      do l=1,N_GT
        fice(l)=0.0
        enddo
!----------------------------------

! calculate the boundaries, based on the layer thicknesses(DZGT)

      zb(1)=-0.05    ! Bottom of surface layer, which is handled outside
                    ! this routine.
      do l=1,N_GT
        zb(l+1)=zb(l)-DZGT(l)
        shc(l)=shr0*(1-phi)*DZGT(l)
        enddo
      do l=1,N_GT
        zc(l)=0.5*(zb(l)+zb(l+1))
        enddo

 
! evaluates the temperatures in the soil layers based on the heat values.  
!             ***********************************
!             input:
!             xw - water in soil layers, m
!             ht - heat in soil layers
!             fsn - heat of fusion of water   J/m
!             shc - specific heat capacity of soil
!             shi - specific heat capacity of ice
!             shw - specific heat capcity of water
!             snowd - snow depth, equivalent water m
!             output:
!             tp - temperature of layers, c
!             fice - fraction of ice of layers
!             pre - extra precipitation, i.e. snowmelt, m s-1
!             snowd - snow depth after melting, equivalent water m.
!             ***********************************
! determine fraction of ice and temp of soil layers based on layer
! heat and water content

      do 10 k=1,N_GT
        ws=phi*DZGT(k)  ! PORE SPACE IN LAYER
        xw=0.5*ws   ! ASSUME FOR THESE CALCULATIONS THAT THE PORE SPACE
                    ! IS ALWAYS HALF FILLED WITH WATER.  XW IS THE  
                    ! AMOUNT OF WATER IN THE LAYER.
        tp(k)=0.
        fice(k)= AMAX1( 0., AMIN1( 1., -ht(k)/(fsn*xw) ) )

        IF(FICE(K) .EQ. 1.) THEN
            tp(k)=(ht(k)+xw*fsn)/(shc(k)+xw*shi0)
          ELSEIF(FICE(K) .EQ. 0.) THEN
            tp(k)=ht(k)/(shc(k)+xw*shw0)
          ELSE
            TP(K)=0.
          ENDIF

    10  continue
        
! evaluates:  layer thermal conductivities
! *****************************************
!             from farouki(cold regions sci and tech, 5, 1981,
!             67-75) the values for tk1,tk2,tk3 are as follows:
!             tk2=2.2**(phi(l,ibv)-xw), tk3=.57**xw, and
!             tk1=3**(1-phi(l,ibv) for %sand<.5, and
!             for computation purposes i have fit these eqs.
!             to 2nd order polynomials.
!             ***********************************
!             input:
!             sklh - soil heat conductivities of layers
!             zb - soil layer boundaries, m
!             zc - soil layer centers, m
!             DZGT - layer thickness, m
!             w - soil water content, m
!             phi - soil porosity, dimensionless
!             q - % sand, silt, clay, peat
!             fice - fraction of ice in layers
!             output:
!             xklh - thermal conductivity, w m-2 k-1
!             ***********************************



! lets get the thermal conductivity for the layers

      do k=1,N_GT

         a1=1-phi              ! ROCK FRACTION
         tk1=1.01692+a1*(0.89865+1.06211*a1)
         xw=phi*(1.-fice(k))   ! FOR SATURATED SOIL, XW HERE IS
                               !   THE LIQUID WATER FRACTION
         a2=phi-xw             ! FOR SATURATED SOIL, THE ICE FRACTION

         tk2=1.00543+a2*(0.723371+.464342*a2)
         tk3=0.998899+xw*(-0.548043+0.120291*xw)
         tksat=tk1*tk2*tk3

         xwi=1.0
         if (zbar .le. zb(k+1))then
             xwi=thetaf
           elseif (zbar .ge. zb(k+1) .and. zbar .le. zb(k))then
             xd1=zb(k)-zbar
             xd2=zbar-zb(k+1)
             xwi=((xd1*thetaf)+xd2)/(xd1+xd2)
           endif 

         xwi=min(xwi,1.)
         tkdry=0.226 ! = .039*0.45^(-2.2), from Farouki, p. 71
         xklh(k)=(tksat-tkdry)*xwi + tkdry

         enddo                



! evaluates heat flux between layers due to heat diffussion
!             ***********************************
!             input:
!             zb - soil layer boundaries, m
!             zc - soil layer centers, m
!             DZGT - layer thickness, m
!             fice - fraction of ice in layers
!             tp - temperature of layers, c
!             shw - specific heat of water
!             shi - specific heat of ice
!             fsn - heat of fusion    J/m
!             output:
!             fh - heat flux between layers
!             ***********************************


! total heat flux is via diffusion along the temperature gradient
      fh(N_GT+1)=0.
      fh(1)=fh21
      do k=2,N_GT
! THIS xkth is NEW (ie., Agnes corrected) - it should be fixed in all
! codes I'm using      
         xkth=((zb(k)-zc(k-1))*xklh(k-1)+(zc(k)-zb(k))*xklh(k))                &
               /(zc(k)-zc(k-1))     
         fh(k)=-xkth*(tp(k-1)-tp(k))/(zc(k-1)-zc(k))
         enddo

! update the heat contents in the model layers; ht(l)
! IF THERE'S SNOW THIS WILL HAVE TO BE MODIFIED L=1,N

      do k=1,N_GT
         ht(k)=ht(k)+(fh(k+1)-fh(k))*dts
         enddo
 


! evaluates the temperatures in the soil  layers based on the heat
! values.  
!             ***********************************
!             input:
!             xw - water in soil layers, m
!             ht - heat in soil layers
!             fsn - heat of fusion of water    J/m
!             shc - specific heat capacity of soil
!             shi - specific heat capacity of ice
!             shw - specific heat capcity of water
!             snowd - snow depth, equivalent water m
!             output:
!             tp - temperature of layers, c
!             fice - fraction of ice of layers
!             pre - extra precipitation, i.e. snowmelt, m s-1
!             snowd - snow depth after melting, equivalent water m.
!             ***********************************
! determine fraction of ice and temp of soil layers based on layer
! heat and water content
      do 1000 k=1,N_GT
         ws=phi*DZGT(k)          ! saturated water content
!            xl=l        
!            xw=(1/(7-xl))*ws
         xw=0.5*ws             ! For calculations here, assume soil
                               ! is half full of water.
         fice(k)=AMAX1( 0., AMIN1( 1., -ht(k)/(fsn*xw) ) )

         IF(FICE(K) .EQ. 1.) THEN
               tp(k)=(ht(k)+xw*fsn)/(shc(k)+xw*shi0)
            ELSEIF(FICE(K) .EQ. 0.) THEN
               tp(k)=ht(k)/(shc(k)+xw*shw0)
            ELSE
               TP(K)=0.
            ENDIF

 1000    continue
  
! determine the value of xfice
      xfice=0.0

      lstart=N_GT
      DO L=N_GT,1,-1
         IF(ZBAR .GE. ZB(L+1))THEN
            LSTART=L
            ENDIF
         ENDDO   

      do l=lstart,N_GT
         xfice=xfice+fice(l)
         enddo
      xfice=xfice/((N_GT+1)-lstart)      

      Return
      end subroutine gndtmp           




!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snowrt(n,                                                     &
           t1,area,tkgnd,precip,snowf,ts,dts, eturb,dedtc,hsturb,dhsdtc,       &
           hlwtc,dhlwtc,desdtc,hlwout,raddn,zc1,                               &
           wesn,htsnn,sndz,fices,tpsn,                                         &
           areasc,areasc0,pre,fhgnd,evap,shflux,lhflux,hcorr                   &
                       )

!*********************************************************************
! AUTHORS:  M. Stieglitz, M. Suarez, R. Koster & S. Dery.
! VERSION:  2003b - This version last updated:  05/30/03.
!*********
! INPUTS:
!*********
!  t1     : Temperature of catchment zones  [C]
!  ts     : Air temperature [K]
!  area   : Fraction of snow-free area in each catchment zone [0-1]
!  precip : Precipitation (Rain+snowfall) [kg/m^2/s == mm/s]
!  snowf  : Snowfall per unit area of catchment [kg/m^2/s == mm/s]
!  dts    : Time step  [s]
!  eturb  : Evaporation per unit area of snow [kg/m^2/s == mm/s]
!  dedtc  : d(eturb)/d(ts) [kg/m^2/s/K]
!  hsturb : Sensible heat flux per unit area of snow  [W/m^2]
!  dhsdtc : d(hsturb)/d(ts)  [W/m^2/K]
!  hlwtc  : Emitted IR per unit area of snow  [W/m^2]
!  dhlwtc : d(hlwtc)/d(ts)  [W/m^2/K]
!  raddn  : Net solar + incident terrestrial per unit area of snow [W/m^2]
!  tkgnd  : Thermal diffusivity of soil in catchment zones [W/m/K]
!  zc1    : Half-thickness of top soil layer [m]
!*********
! UPDATES:
!*********
!  wesn   : Layer water contents per unit area of catchment [kg/m^2]
!  htsnn  : Layer heat contents relative to liquid water at 0 C [J/m^2]
!  sndz   : Layer depths [m]
!*********
! OUTPUTS: 
!*********
!  tpsn   : Layer temperatures [C]
!  fices  : Layer frozen fraction [0-1]
!  areasc : Areal snow coverage at beginning of step [0-1]
!  areasc0: Areal snow coverage at end of step [0-1]
!  pre    : Liquid water flow from snow base [kg/m^2/s]
!  fhgnd  : Heat flux at snow base at catchment zones  [W/m^2]
!  hlwout : Final emitted IR flux per unit area of snow [W/m^2]
!  lhflux : Final latent heat flux per unit area of snow [W/m^2]
!  shflux : Final sensible heat flux per unit area of snow   [W/m^2]
!  evap   : Final evaporation per unit area of snow   [kg/m^2/s]
!*********************************************************************
! NOTA:  By convention, wesn is representative for a catchment area
! equal to 1 whereas sndz is relative to the area covered by snow only.
!*********************************************************************


      implicit none

!      real, parameter :: lhv    = 2.4548E6 !  2.5008e6   !  @ 0 C [J/kg]
!      real, parameter :: lhs    = 2.8368E6 !  2.8434e6 !  @ 0 C [J/kg]
!      real, parameter :: lhf    = (lhs-lhv)  !  @ 0 C [J/kg]
      real, parameter :: cpw    = 2065.22    !  @ 0 C [J/kg/K]
      
!rr      real, parameter :: cpw_liquid = 4185. ! [J/kg/K]

!      real, parameter :: tfrz   = 273.16     !  @ 0 C [K]
!      real, parameter :: rhofs  = 150.       !  [kg/m^3]
!      real, parameter :: rhoma  = 500.       !  [kg/m^3]
!      real, parameter :: rhow   = 1000.      !  [kg/m^3]
!      real, parameter :: wemin  = 13.        !  [kg/m^2]
      real, parameter :: snfr   = 0.01       !  holding capacity
      real, parameter :: small  = 1.e-6      !  small number
!      integer, parameter :: nlay = 3         !  number of layers
!      integer, parameter :: N_sm   = 3         !  number of zones
!      real   , parameter :: MIN_SNOW_MASS = .013 ! kg/M**2 equiv to 0.1% area

 
      integer, intent(in)  :: n 
      real,    intent(in ) :: t1(N_sm),area(N_sm),tkgnd(N_sm)
      real,    intent(in ) :: ts,precip,snowf,dts,dedtc,raddn,hlwtc
      real,    intent(in ) :: dhsdtc,desdtc,dhlwtc,eturb,hsturb,zc1
      real,    intent(inout)::wesn(N_snow),htsnn(N_snow),sndz(N_snow)
      real,    intent(out) :: tpsn(N_snow),fices(N_snow),fhgnd(N_sm)
      real,    intent(out) :: hlwout,lhflux,shflux,areasc0,evap,areasc

!Locals
      real :: tsx, mass,snowd,rainf,denom,alhv,lhturb,dlhdtc,hcorr,            &
           enew,eold,tdum,fnew,tnew,icedens,densfac,hnew,scale,t1ave,          &
           flxnet,fdum,dw,pre,waterin,waterout,snowin,snowout,                 &
           waterbal,precision,flow,term,dz,w(0:N_snow),HTSPRIME
      real, dimension(size(wesn)  ) :: tksn,cmpc,excs,dens
      real, dimension(size(wesn)  ) :: dtc,q,cl,cd,cr
      real, dimension(size(wesn)+1) :: fhsn,df
      real, dimension(size(wesn)  ) :: htest,ttest,ftest

      logical, dimension(size(wesn)  ) ::  ice1,tzero, ice10,tzero0

      integer :: i,izone
      logical :: logdum,kflag

       snowd = sum(wesn)
       snowin = snowd

!rr   correction for "cold" snow
       tsx   = min(ts-tf,0.)*cpw
       
!rr   correction for heat content of rain
!rr       tsx_rain = max(ts-tf,0.)*cpw_liquid
       
       df     = 0.
       dtc    = 0.
       tpsn   = 0.
       fices  = 0.
       areasc = 0.
       areasc0= 0.
       pre    = 0.
       fhgnd  = 0.
       hlwout = 0.
       shflux = 0.
       lhflux = 0.
       evap   = 0.
       excs   = 0.
       hcorr  = 0.
       dens   = rhofs
       rainf  = precip - snowf   ! [kg/m^2/s]

       
       if(snowd <= MINSWE) then ! no snow
!         Assume initial (very small) snow water melts; new area is based
!         on new snowfall

          areasc = min(snowd/wemin,1.)
          areasc0 = 0.
          pre = snowd/dts + areasc*rainf
          wesn  = 0.
          hcorr = hcorr + raddn*areasc + sum(htsnn)/dts
          htsnn = 0.
          sndz  = 0.

          if(snowf > 0.) then   ! only initialize with non-liquid part of precip
                                ! liquid part runs off (above)

             wesn    = snowf*dts/float(N_snow)  
             htsnn   = (tsx-alhm)*wesn
             areasc0 = min((snowf*dts)/wemin,1.)
             sndz = wesn/rhofs
!             hcorr  = hcorr - (tsx-alhm)*snowf    ! randy
             hcorr  = hcorr - tsx*snowf           ! randy
             call relayer(htsnn, wesn, sndz)
             call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)
             
          endif

        return ! if there was no snow at start of time step

     endif

!**** Determine the fractional snow coverage

       areasc = min(snowd/wemin,1.)

!**** Set the mean density & diffusivity of the layers

       do i=1,N_snow
         if(sndz(i) > 0) dens(i) = max(wesn(i)/(areasc*sndz(i)),rhofs)
       enddo
       tksn = 3.2217e-06*dens**2

!**** Determine temperature & frozen fraction of snow layers

         call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)

!**** Calculate the ground-snow energy flux at 3 zones

       denom = 1./(sndz(N_snow)*0.5-zc1)
       fhgnd = -sqrt(tkgnd*tksn(N_snow))*area*denom*(tpsn(N_snow)-t1)
       fhsn(N_snow+1) = sum(fhgnd)
       do i=1,N_sm
        df(N_snow+1)=df(N_snow+1)-sqrt(tkgnd(i)*tksn(N_snow))*area(i)*denom
       enddo

!**** Ensure against excessive heat flux between ground and snow:
!**** if heat flux found to cause the lowest snow layer temperature
!**** to "overshoot" (e.g. to become higher than the ground temperature
!**** when it had been lower), reduce the heat flux.  If the lowest 
!**** snow layer starts off at zero and the new temperature is greater
!**** than zero, reduce the heat flux to melt only half of the lowest
!**** layer snow.
!**** 
      t1ave=sum(t1*area)/sum(area)
      htest=htsnn
      htest(N_snow)=htest(N_snow)+fhsn(N_snow+1)*dts*areasc
      
      call get_tf_nd(N_snow, htest, wesn, ttest, ftest)

      scale=1.
      if((t1ave-tpsn(N_snow))*(t1ave-ttest(N_snow)) .lt. 0.) then
         scale=0.5*(tpsn(N_snow)-t1ave)/(tpsn(N_snow)-ttest(N_snow))
         endif
      if(tpsn(N_snow) .eq. 0. .and. ttest(N_snow) .gt. 0. .and.                &
                abs(fhsn(N_snow+1)) .gt. 1.e-10) then
         scale=(-0.5*htsnn(N_snow)/(dts*areasc))/fhsn(N_snow+1)
         endif

      fhsn(N_snow+1)=fhsn(N_snow+1)*scale
         df(N_snow+1)=df(N_snow+1)*scale
         fhgnd=fhgnd*scale


!**** Calculate heat fluxes between snow layers.
       do i=2,N_snow
          df(i) =  -sqrt(tksn(i-1)*tksn(i))/((sndz(i-1)+sndz(i))*0.5)
          fhsn(i)= df(i)*(tpsn(i-1)-tpsn(i))
      enddo
 
!**** Effective heat of vaporization includes bringing snow to 0 C

        alhv   = alhe + alhm                            !randy
!        alhv   = alhe + fices(1)*alhm + tpsn(1)*cpw    !randy

!**** Initial estimate of latent heat flux change with Tc

        lhturb = alhv*eturb
        dlhdtc = alhv*dedtc

!**** Initial estimate of net surface flux & its change with Tc

        fhsn(1) = lhturb + hsturb + hlwtc - raddn
        df(1)   = -(dlhdtc + dhsdtc + dhlwtc)

!**** Prepare array elements for solution & coefficient matrices.
!**** Terms are as follows:  left (cl), central (cd) & right (cr)
!**** diagonal terms in coefficient matrix & solution (q) terms.

        do i=1,N_snow

         call get_tf0d(htsnn(i),wesn(i),tdum,fdum, ice1(i),tzero(i))

         if(ice1(i)) then
           cl(i) = df(i)
           cd(i) = cpw*wesn(i)/dts - df(i) - df(i+1)
           cr(i) = df(i+1)
           q(i)  = fhsn(i+1)-fhsn(i)
         else
           cl(i) = 0.
           cd(i) = 1.
           cr(i) = 0.
           q(i)  = 0.
         endif

        enddo
        cl(1)    = 0.
        cr(N_snow) = 0.

        do i=1,N_snow-1
          if(.not.ice1(i)) cl(i+1) = 0.
        enddo

        do i=2,N_snow
          if(.not.ice1(i)) cr(i-1) = 0.
        enddo

!**** Solve the tri-diagonal matrix for implicit change in Tc.

        call TRID(dtc,cl,cd,cr,q,N_snow)

!**** Check temperature changes for passages across critical points,i.e.
!**** If implicit change has taken layer past melting/freezing, correct.

       do i=1,N_snow
          if(tpsn(i)+dtc(i) > 0. .or. htsnn(i)+wesn(i)*cpw*dtc(i) > 0.) then
             dtc(i)=-tpsn(i)
             endif
          if(.not.ice1(i)) dtc(i)=0.
          enddo
 
!**** Further adjustments; compute new values of h associated with
!**** all adjustments.

       eold=sum(htsnn)

       do i=1,N_snow

!**** Quick check for "impossible" condition:

          if(.not.tzero(i) .and. .not.ice1(i)) then
             write(*,*) 'bad snow condition: fice,tpsn =',fices(i),tpsn(i)
             stop
             endif

!****  Condition 1: layer starts fully frozen (temp < 0.)

          if(.not.tzero(i)) then
             tnew=tpsn(i)+dtc(i)
             fnew=1.

             endif

!****  Condition 2: layer starts with temp = 0, fices < 1.
!      Corrections for flxnet calculation: Koster, March 18, 2003.

          if(.not.ice1(i)) then
             tnew=0.
             if(i==1) flxnet= fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))              &
                   -fhsn(i)-df(i)*dtc(i)
             if(i > 1 .and. i < N_snow) flxnet=                                &
                    fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))                        &
                   -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
             if(i==N_snow) flxnet=fhsn(i+1)+df(i+1)*dtc(i)                     &
                   -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
             HTSPRIME=HTSNN(I)+AREASC*FLXNET*DTS
             call get_tf0d(HTSPRIME,wesn(i),  tdum,fnew,logdum,logdum)
             fnew=amax1(0.,  amin1(1.,  fnew))

             endif

!****  Condition 3: layer starts with temp = 0, fices = 1.
!      Corrections for flxnet calculation: Koster, March 18, 2003.

          if(ice1(i) .and. tzero(i)) then
             if(dtc(i) < 0.) then
                tnew=tpsn(i)+dtc(i)
                fnew=1.
                endif
             if(dtc(i) >= 0.) then
                tnew=0.
                if(i==1) flxnet=fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))            &
                       -fhsn(i)-df(i)*dtc(i)
                if(i > 1 .and. i < N_snow) flxnet=                             &
                    fhsn(i+1)+df(i+1)*(dtc(i)-dtc(i+1))                        &
                   -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))
                if(i==N_snow) flxnet=fhsn(i+1)+df(i+1)*dtc(i)                  &
                   -fhsn(i)-df(i)*(dtc(i-1)-dtc(i))

                HTSPRIME=HTSNN(I)+AREASC*FLXNET*DTS
                call get_tf0d(HTSPRIME,wesn(i),   tdum,fnew,logdum,logdum)
                fnew=amax1(0.,  amin1(1.,  fnew))
                endif
             endif

!**** Now update heat fluxes & compute sublimation or deposition.

         if(i == 1) then
             lhflux = lhturb + dlhdtc*dtc(1)
             shflux = hsturb + dhsdtc*dtc(1)
             hlwout = hlwtc  + dhlwtc*dtc(1)
             evap = lhflux/alhv
             dw = -evap*dts*areasc
             if(-dw > wesn(1) ) then
                dw = -wesn(1)
                evap = -dw/(dts*areasc)
!                shflux=shflux+(lhflux-evap*alhv)
                hcorr=hcorr+(lhflux-evap*alhv)*areasc
                lhflux=evap*alhv
                endif
             wesn(1)  = wesn(1) + dw
             denom = 1./dens(1)
             if(dw > 0.) denom = 1./rhoma
             sndz(1) = sndz(1) + dw*denom
             endif

         if(i == N_snow) then
             do izone=1,N_sm
                fhgnd(izone)=fhgnd(izone)+area(izone)*df(N_snow+1)*dtc(N_snow)
                enddo
             endif

!**** Now update thermodynamic quantities.

          htsnn(i)=(cpw*tnew-fnew*alhm)*wesn(i)
          tpsn(i) = tnew    
          fices(i)= fnew
        enddo

!**** Store excess heat in hcorr.

       enew=sum(htsnn)
       hcorr=hcorr-((enew-eold)/dts+areasc*(lhflux+shflux+hlwout-raddn)        &
                      -areasc*(fhsn(N_snow+1)+df(N_snow+1)*dtc(N_snow))        &
                      )

       call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)

!rr!**** Add rainwater and snow at ts., bal. budget with shflux.
!rr   (tried and failed 19 Jun 2003, reichle)
!rr
!rr       wesn (1) = wesn (1) + (rainf*areasc+snowf)*dts
!rr       htsnn(1) = htsnn(1) + (tsx -alhm)*(snowf*dts) + tsx_rain*rainf*dts
!rr       sndz (1) = sndz (1) + (snowf/rhofs)*dts
!rr       !  shflux   = shflux   + tsx*snowf                   ! randy
!rr       hcorr   = hcorr   - (tsx-alhm)*snowf - tsx_rain*rainf ! randy


!**** Add rainwater at 0 C, snow at ts., bal. budget with shflux.

       wesn (1) = wesn (1) + (rainf*areasc+snowf)*dts
       htsnn(1) = htsnn(1) + (tsx -alhm)*(snowf*dts)
       sndz (1) = sndz (1) + (snowf/rhofs)*dts
!       shflux  = shflux   + tsx*snowf          ! randy
!       hcorr   = hcorr   - (tsx-alhm)*snowf     ! randy
        hcorr   = hcorr   -  tsx*snowf          ! randy
       
       snowd=sum(wesn)

       call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)

!**** Move meltwater through the pack.
!**** Updated by Koster, August 27, 2002.

       pre = 0.
       flow = 0.

       do i=1,N_snow

        if(flow > 0.) then
         wesn (i) =  wesn(i) + flow
         call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)  
      endif

        pre  = (max((1.-fices(i))*wesn(i),0.))
        flow = 0.

        if(snowd > wemin) then

          icedens=wesn(i)*fices(i)/(sndz(i)+1.e-20)
          densfac=amax1(0., amin1(1., icedens/rhofs))
          term=densfac*snfr*(sndz(i)*rhow-wesn(i)*fices(i))
      
          if(pre > term) then
            pre = min(pre - term, wesn(i))
            wesn(i) = wesn(i) - pre
            flow = pre
          endif
        else
          wesn(i) = wesn(i) - pre
          flow = pre
       endif

!**** Adjust top layer snow depth to get proper density values
!**** But limit this change for large throughflow (STEPH 06/19/03)

        if(i==1)then
          dz=min(flow/dens(i),0.5*sndz(i))
          sndz(i)=sndz(i)-dz
        endif
       enddo

       pre = flow/dts
       snowd=sum(wesn)

!**** Update snow density by compaction (Pitman et al. 1991)

       excs = 0.
       mass = 0.
       w    = 0.

       if(snowd > wemin) then ! Compaction only after full coverage.

          do i=1,N_snow
             dens(i) = rhofs
             if(sndz(i)>0.) dens(i) = max(wesn(i)/(sndz(i)),rhofs)
          enddo
          
          cmpc    = exp(14.643 - (4000./min(tpsn+tf,tf))-.02*dens)

          do i=1,N_snow
             w(i) = wesn(i)
             mass = mass + 0.5*(w(i)+w(i-1))
             dens(i) = dens(i)*(1. + (dts*0.5e-7*9.81)*mass*cmpc(i))
             
!**** Clip densities below maximum value, adjust quantities accordingly
!**** while conserving heat & mass (STEPH 06/21/03).

             if(dens(i) > rhoma) then
                excs(i) = (dens(i)-rhoma)*sndz(i)
                wesn(i) = wesn(i) - excs(i)
                hnew = (cpw*tpsn(i)-fices(i)*alhm)*wesn(i)
                hcorr= hcorr+(htsnn(i)-hnew)/dts
                htsnn(i)= hnew
                dens(i) = rhoma
             endif
          enddo
       endif

       pre  = pre + sum(excs)/dts
       sndz = wesn/dens
       
!**** Restore layers to sigma values.
       
       do i=1,N_snow
         call get_tf0d(htsnn(i),wesn(i),tdum,fdum,ice10(i),tzero0(i))
         enddo

       call relayer(htsnn, wesn, sndz)

       call get_tf_nd(N_snow, htsnn, wesn, tpsn, fices)
       
!**** Check that (ice10,tzero) conditions are conserved through
!**** relayering process (or at least that (fices,tpsn) conditions don't 
!**** go through the (1,0) point); excess goes to hcorr.

       do i=1,N_snow
          kflag=.false.
          if(ice10(i).and.tzero0(i) .and.                                      &
             (fices(i) .ne. 1. .or. tpsn(i) .ne. 0.) ) kflag=.true.
          if(.not.ice10(i).and.tzero0(i) .and.                                 &
             (fices(i) .eq. 1. .and. tpsn(i) .lt. 0.) ) kflag=.true.
          if(ice10(i).and. .not.tzero0(i) .and.                                &
             (fices(i) .ne. 1. .and. tpsn(i) .eq. 0.) ) kflag=.true.

          if(kflag) then
             hnew=-alhm*wesn(i)
             hcorr=hcorr+(htsnn(i)-hnew)/dts
             htsnn(i)=hnew
             tpsn(i)=0.
             fices(i)=1.
             endif

          enddo


!**** Reset fractional area coverage.

       areasc0 = min(sum(wesn)/wemin,1.)
 
!**** Final check for water balance.

       waterin   = (rainf*areasc+snowf)*dts + max(dw,0.)
       waterout  = pre*dts - min(dw,0.)
       snowout   = sum(wesn)
       waterbal  = snowin + waterin - waterout - snowout
       precision = snowout*small
       if((waterbal > precision).and.(waterbal > small)) then
        write(*,*) 'Warning: Imbalance in snow water budget!', waterbal
        write(*,*) snowin,snowout,waterin, waterout
        write(*,*) precision, small
       endif

      return  !  end snow

      end subroutine snowrt

! **********************************************************************

      subroutine relayer(htsnn, wesn, sndz)
      
      implicit none
!      integer, parameter ::  N_snow = 3
      real, intent(inout) :: htsnn(N_snow),wesn(N_snow),sndz(N_snow)
      
      real, dimension(size(sndz),2)   :: ds
      real, dimension(size(sndz)+1)   :: sdold, sdnew
      real, dimension(size(sndz)+1,2) :: h, s
      
      integer :: i, j
      
!      real, parameter :: dz1max = 0.08   ! [m]
!      real, parameter :: wemin  = 13.0   ! [kg/m2]
      real, parameter :: small  = 1.e-20 
      real :: areasc,dz
      
!**** Initialize some variables.
      
      h  = 0.
      s  = 0.
      ds = 0.
      dz = 0.
      
      areasc = min(sum(wesn)/wemin,1.)
      
!**** Compute specific heat & water contents of old layers.

      do i=1,N_snow
         if (sndz(i) > 0.) then
            h(i,1) = htsnn(i)/sndz(i)
            h(i,2) =  wesn(i)/sndz(i)
         endif
      enddo
      
!**** Obtain old & new layer thicknesses & boundaries.
      
      sdold = 0.
      sdnew = 0.
      
      do i=N_snow,1,-1
         sdold(i) = sdold(i+1) + sndz(i)
      enddo
      
      sndz = sdold(1)/float(N_snow)
      if(sndz(1) > dz1max) then
         sndz(2:) = (sdold(1)-dz1max)/float(N_snow-1)
         sndz(1)  = dz1max
      endif
      
      do i=N_snow,1,-1
         sdnew(i) = sdnew(i+1) + sndz(i)
      enddo
      
!**** Since the snow boundary has moved, redistribute heat  
!     contents & water equivalents of old to new snow layers.
                                                              
      do i=1,N_snow
         
         j = i
         dz=sdnew(i+1)-sdold(i+1)
         if(dz < 0.) j = i + 1
         s(i+1,:) = h(j,:)*dz
         ds(i,:)  = s(i,:) - s(i+1,:)
      enddo

      htsnn = htsnn + ds(:,1)
      wesn  = wesn  + ds(:,2) 
      
      if(sum(wesn) < wemin) sndz = sndz /(areasc + small)
      return
      
      end subroutine relayer

! **********************************************************************

      subroutine get_tf0d(h,w,t,f,ice1,tzero)

      implicit none
      
      real, parameter :: cpw    = 2065.22   !  @ 0 C [J/kg/K]
!      real, parameter :: lhv    = 2.4548E6 !  2.5008e6   !  @ 0 C [J/kg]
!      real, parameter :: lhs    = 2.8368E6 !  2.8434e6 !  @ 0 C [J/kg]
!rr   real, parameter :: lhv    = 2.5008e6  !  @ 0 C [J/kg]
!rr   real, parameter :: lhs    = 2.8434e6  !  @ 0 C [J/kg]
!      real, parameter :: lhf    = (lhs-lhv) !  @ 0 C [J/kg]
      real, parameter :: tfac=1./cpw
      real, parameter :: ffac=1./alhm
      
      real, intent(in )   :: w, h
      real, intent(out)   :: t, f

      logical, intent(out) :: ice1,tzero
      
      real :: hbw
      
      hbw=0.
      if(w > 0.) hbw = h/w
      
      if(hbw < -1.00001*alhm) then
            t = (hbw+alhm)*tfac
            f = 1.
            ice1=.true.
            tzero=.false.
         elseif(hbw > -0.99999*alhm) then
            t = 0.
            f =-hbw*ffac
            ice1=.false.
            tzero=.true.
         else
            t = 0.
            f = 1.
            ice1=.true.
            tzero=.true.
         endif
      
      if(f < 0.) then
         t = hbw*tfac
         f = 0.
         endif
      
      if(w == 0.) then
         t = 0.
         f = 0.
         endif

      return

      end subroutine get_tf0d
      
! **********************************************************************
      
      subroutine get_tf_nd(N,h,w,t,f)
      
        use LIS_coreMod
        use LIS_logMod
!     n-dimensional version of get_tf
!     
!     avoid slow "where" statements
!     
!     can be called for any number of layers or catchments, for example

!     1.) call get_tf_nd( ncatm, htsnn1(1:ncatm), wesn1(1:ncatm),
!                         tpsn(1:ncatm),f(1:ncatm) )
!     
!     2.) call get_tf_nd(N_snow, h, w, t, f)

!     reichle, 22 Aug 2002
!     reichle, 29 Apr 2003 (updated parameter values)

      integer, intent(in) :: N
      
      real, dimension(n), intent(in)    :: h, w
      real, dimension(n), intent(out)   :: t, f
      
!     local variables
   
      real, parameter :: cpw    = 2065.22   !  @ 0 C [J/kg/K]
!      real, parameter :: lhv    = 2.4548E6 !  2.5008e6   !  @ 0 C [J/kg]
!      real, parameter :: lhs    = 2.8368E6 !  2.8434e6 !  @ 0 C [J/kg]
!rr   real, parameter :: lhv    = 2.5008e6  !  @ 0 C [J/kg]
!rr   real, parameter :: lhs    = 2.8434e6 !  @ 0 C [J/kg]
!      real, parameter :: lhf    = (lhs-lhv) !  @ 0 C [J/kg]
      real, parameter :: tfac=1./cpw
      real, parameter :: ffac=1./alhm
      
      integer :: i      
      
      real :: hbw
            
      do i=1,N
         
         if(w(i) .gt. 0.0) then
            hbw = h(i)/w(i)
         else
            hbw = 0.
         endif
         if(hbw .lt. -alhm) then
            t(i) = (hbw+alhm)*tfac
            f(i) = 1.
         elseif(hbw .gt. -alhm) then
            t(i) =  0.
            f(i) = -hbw*ffac
         else
            t(i) = 0.
            f(i) = 1.
         endif
         
         if(f(i) .lt. 0.) then
            t(i) = hbw*tfac
            f(i) = 0.
         endif
         
         if(w(i) .eq. 0.) then
            t(i) = 0.
            f(i) = 0.
         endif
         
      end do
      
      return
      
    end subroutine get_tf_nd
      
! **********************************************************************

      SUBROUTINE TRID(X,DD,D,RD,B,N)
      IMPLICIT NONE

      INTEGER,INTENT(IN) :: N
      REAL*4, INTENT(IN), DIMENSION(N) :: DD, RD
      REAL*4, INTENT(INOUT), DIMENSION(N) :: D, B
      REAL*4, INTENT(OUT),DIMENSION(N) :: X

      integer I,J
      real*4  RSF
      RSF=0.
      DO 10 I=2,N
         J=N+1-I
         if(D(J+1).ne.0.) RSF=RD(J)/D(J+1)
         D(J)=D(J)-DD(J+1)*RSF
   10 B(J)=B(J)- B(J+1)*RSF
      if(D(1).ne.0.) X(1)=B(1)/D(1)
      DO 20 J=2,N
   20 if(D(J).ne.0.) X(J)=(B(J)-DD(J)*X(J-1))/D(J)
      RETURN
      END SUBROUTINE TRID






      SUBROUTINE SIBALB (                                                      &
                         NCH, ITYP, VLAI, VGRN, ZTH,                           &
                         SCALVDR,SCALVDF,SCALIDR,SCALIDF,                      &
                         WESN,SNDZ,                                            &
      			 AVISDR, ANIRDR, AVISDF, ANIRDF,                       &
                         ASNVDR, ASNNDR, ASNVDF, ASNNDF                        &
              		)

      IMPLICIT NONE

! OUTPUTS:
! AVISDR:   visible, direct albedo.
! ANIRDR:   near infra-red, direct albedo.
! AVISDF:   visible, diffuse albedo.
! ANIRDF:   near infra-red, diffuse albedo.

! INPUTS:
! SCALVDR:  MODIS scale factor for visible, direct.
! SCALVDF:  MODIS scale factor for visible, diffuse.
! SCALIDR:  MODIS scale factor for NIR, direct.
! SCALIDF:  MODIS scale factor for NIR, diffuse.
! VLAI:     the leaf area index.
!VGRN:     the greenness index.
! ZTH:      The cosine of the solar zenith angle.
! SNW:      Snow cover in meters water equivalent.


      INTEGER, INTENT(IN) :: NCH
      INTEGER, INTENT(IN), DIMENSION(NCH) ::  ITYP
      REAL, INTENT(IN), DIMENSION(NCH) :: SCALVDR,SCALVDF,SCALIDR,SCALIDF,     &
                                             VLAI, VGRN,  ZTH
      REAL, INTENT(IN), DIMENSION(3,NCH) :: WESN, SNDZ
      REAL, INTENT(OUT), DIMENSION(NCH) :: AVISDR, ANIRDR, AVISDF,             &
                        ANIRDF, ASNVDR, ASNNDR, ASNVDF, ASNNDF


      REAL, PARAMETER :: ALVDRS = 0.100
      REAL, PARAMETER :: ALIDRS = 0.200
      REAL, PARAMETER :: ALVDRD = 0.300
      REAL, PARAMETER :: ALIDRD = 0.350
      REAL, PARAMETER :: ALVDRI = 0.700
      REAL, PARAMETER :: ALIDRI = 0.700


!      REAL, PARAMETER :: WEMIN  = 13.0   ! [KG/M2]

! ALVDRS:  Albedo of soil for visible   direct  solar radiation.
! ALIDRS:  Albedo of soil for infra-red direct  solar radiation.
! ALVDFS:  Albedo of soil for visible   diffuse solar radiation.
! ALIDFS:  Albedo of soil for infra-red diffuse solar radiation.

      INTEGER, PARAMETER :: NLAI = 14

      REAL, PARAMETER :: EPSLN = 1.E-6
      REAL, PARAMETER :: BLAI = 0.5
      REAL, PARAMETER :: DLAI = 0.5

      REAL, PARAMETER :: ALATRM = BLAI + (NLAI - 1) * DLAI - EPSLN

      INTEGER, PARAMETER :: NTYPS_SIB=9

      REAL :: SWE, TOTDEP, AREASC, DENSITY, DENS_EXC, FRACV, SNWMASK,          &
                  AMASK, ASNVDR_VEG, ASNNDR_VEG, ASNVDF_VEG, ASNNDF_VEG





! ITYP: Vegetation type as follows:
!                  1:  BROADLEAF EVERGREEN TREES
!                  2:  BROADLEAF DECIDUOUS TREES
!                  3:  NEEDLELEAF TREES
!                  4:  GROUND COVER
!                  5:  BROADLEAF SHRUBS
!                  6:  DWARF TREES (TUNDRA)
!                  7:  BARE SOIL
!                  8:  DESERT
!                  9:  ICE
!  NCH: Chip index
!

	INTEGER I, LAI
	REAL FAC, GAMMA, BETA, ALPHA, DX, DY, ALA, FVEG
        REAL, DIMENSION(2) :: GRN
        REAL, DIMENSION(4,NTYPS_SIB) :: SNWALB (4, NTYPS_SIB)
        REAL, DIMENSION(NTYPS_SIB) :: SNWMSK


      DATA GRN /0.33, 0.67/
      REAL, PARAMETER :: SNWALB_VISMAX = 0.70 
      REAL, PARAMETER :: SNWALB_VISMIN = 0.5
      REAL, PARAMETER :: SNWALB_NIRMAX = 0.5
      REAL, PARAMETER :: SNWALB_NIRMIN = 0.3
!      REAL, PARAMETER :: RHOFS = 150.                ! DENSITY OF FRESH SNOW
      REAL, DIMENSION(NTYPS_SIB) :: SNWMID

!       DATA SNWALB/.85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50,                                     &
!                   .85, .50, .85, .50                                      &
!      		  /

!      DATA SNWMSK/25., 5., 10., 0.2, 0.5, 0.2, 0.1, 0.1, 0.1/

!***  grassland and tundra values arbitrarily increased.
      DATA SNWMID /50.,30.,45.,20.,30.,20.,2.,2.,2./




! [ Definition of Functions: ]
!
!	REAL clsmf25_COEFFSIB

! --------------------------------------------------



!   Constants used in albedo calculations:

      REAL ALVDR (NLAI, 2, NTYPS_SIB)
      REAL BTVDR (NLAI, 2, NTYPS_SIB)
      REAL GMVDR (NLAI, 2, NTYPS_SIB)
      REAL ALIDR (NLAI, 2, NTYPS_SIB)
      REAL BTIDR (NLAI, 2, NTYPS_SIB)
      REAL GMIDR (NLAI, 2, NTYPS_SIB)
      
!  (Data statements for ALVDR described in full; data statements for
!   other constants follow same framework.)


!    BROADLEAF EVERGREEN (ITYP=4); GREEN=0.33; LAI: .5-7
	DATA (ALVDR (I, 1, 1), I = 1, 14)                                      &
      	  /0.0808, 0.0796, 0.0792, 0.0790, 10*0.0789/

!    BROADLEAF EVERGREEN (ITYP=4); GREEN=0.67; LAI: .5-7
	DATA (ALVDR (I, 2, 1), I = 1, 14)                                      &
      	  /0.0788, 0.0775, 0.0771, 0.0769, 10*0.0768/

!    BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.33; LAI: .5-7
	DATA (ALVDR (I, 1, 2), I = 1, 14)                                      &
      	  /0.0803, 0.0790, 0.0785, 0.0784, 3*0.0783, 7*0.0782/

!    BROADLEAF DECIDUOUS (ITYP=1); GREEN=0.67; LAI: .5-7
	DATA (ALVDR (I, 2, 2), I = 1, 14)                                      &
      	  /0.0782, 0.0770, 0.0765, 0.0763, 10*0.0762/

!    NEEDLELEAF (ITYP=3); GREEN=0.33; LAI=.5-7
	DATA (ALVDR (I, 1, 3), I = 1, 14)                                      &
      	  /0.0758, 0.0746, 0.0742, 0.0740, 10*0.0739/

!    NEEDLELEAF (ITYP=3); GREEN=0.67; LAI=.5-7
	DATA (ALVDR (I, 2, 3), I = 1, 14)                                      &
      	  /0.0683, 0.0672, 0.0667, 2*0.0665, 9*0.0664/

!    GROUNDCOVER (ITYP=2); GREEN=0.33; LAI=.5-7    
	DATA (ALVDR (I, 1, 4), I = 1, 14)                                      &
      	  /0.2436, 0.2470, 0.2486, 0.2494, 0.2498, 0.2500, 2*0.2501,           &
      		6*0.2502 /

!    GROUNDCOVER (ITYP=2); GREEN=0.67; LAI=.5-7
	DATA (ALVDR (I, 2, 4), I = 1, 14) /14*0.1637/

!    BROADLEAF SHRUBS (ITYP=5); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 5), I = 1, 14)                                      &
          /0.0807, 0.0798, 0.0794, 0.0792, 0.0792, 9*0.0791/

!    BROADLEAF SHRUBS (ITYP=5); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 5), I = 1, 14)                                      &
          /0.0787, 0.0777, 0.0772, 0.0771, 10*0.0770/

!    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.33,LAI=.5-7
        DATA (ALVDR (I, 1, 6), I = 1, 14)                                      &
          /0.0802, 0.0791, 0.0787, 0.0786, 10*0.0785/

!    DWARF TREES, OR TUNDRA (ITYP=6); GREEN=0.67,LAI=.5-7
        DATA (ALVDR (I, 2, 6), I = 1, 14)                                      &
          /0.0781, 0.0771, 0.0767, 0.0765, 0.0765, 9*0.0764/


!    BARE SOIL
	DATA (ALVDR (I, 1, 7), I = 1, 14) /14*ALVDRS/
	DATA (ALVDR (I, 2, 7), I = 1, 14) /14*ALVDRS/

!    DESERT
	DATA (ALVDR (I, 1, 8), I = 1, 14) /14*ALVDRD/
	DATA (ALVDR (I, 2, 8), I = 1, 14) /14*ALVDRD/

!    ICE
	DATA (ALVDR (I, 1, 9), I = 1, 14) /14*ALVDRI/
	DATA (ALVDR (I, 2, 9), I = 1, 14) /14*ALVDRI/
!****
!**** -------------------------------------------------
	DATA (BTVDR (I, 1, 1), I = 1, 14)                                      &
        /0.0153, 0.0372, 0.0506, 0.0587, 0.0630, 0.0652, 0.0663,               &
      	0.0668, 0.0671, 0.0672, 4*0.0673 /
	DATA (BTVDR (I, 2, 1), I = 1, 14)                                      &
     	  /0.0135, 0.0354, 0.0487, 0.0568, 0.0611, 0.0633, 0.0644,             &
     	0.0650, 0.0652, 0.0654, 0.0654, 3*0.0655 /
	DATA (BTVDR (I, 1, 2), I = 1, 14)                                      &
      	  /0.0148, 0.0357, 0.0462, 0.0524, 0.0554, 0.0569, 0.0576,             &
      	0.0579, 0.0580, 0.0581, 0.0581, 3*0.0582 /
	DATA (BTVDR (I, 2, 2), I = 1, 14)                                      &
      	  /0.0131, 0.0342, 0.0446, 0.0508, 0.0539, 0.0554, 0.0560,             &
      	0.0564, 0.0565, 5*0.0566 /
	DATA (BTVDR (I, 1, 3), I = 1, 14)                                      &
      	  /0.0108, 0.0334, 0.0478, 0.0571, 0.0624, 0.0652, 0.0666,             &
      	0.0673, 0.0677, 0.0679, 4*0.0680 /
	DATA (BTVDR (I, 2, 3), I = 1, 14)                                      &
      	  /0.0034, 0.0272, 0.0408, 0.0501, 0.0554, 0.0582, 0.0597,             &
      		0.0604, 0.0608, 0.0610, 4*0.0611 /
	DATA (BTVDR (I, 1, 4), I = 1, 14)                                      &
      	  /0.2050, 0.2524, 0.2799, 0.2947, 0.3022, 0.3059, 0.3076,             &
      		0.3085, 0.3088, 0.3090, 4*0.3091 /
	DATA (BTVDR (I, 2, 4), I = 1, 14)                                      &
      	  /0.1084, 0.1404, 0.1617, 0.1754, 0.1837, 0.1887, 0.1915,             &
      		0.1931, 0.1940, 0.1946, 0.1948, 0.1950, 2*0.1951  /
        DATA (BTVDR (I, 1, 5), I = 1, 14)                                      &
          /0.0203, 0.0406, 0.0548, 0.0632, 0.0679, 0.0703, 0.0716,             &
           0.0722, 0.0726, 0.0727, 0.0728, 0.0728, 0.0728, 0.0729 /
        DATA (BTVDR (I, 2, 5), I = 1, 14)                                      &
          /0.0184, 0.0385, 0.0526, 0.0611,  0.0658, 0.0683, 0.0696,            &
           0.0702, 0.0705, 0.0707, 4*0.0708 /
        DATA (BTVDR (I, 1, 6), I = 1, 14)                                      &
          /0.0199, 0.0388, 0.0494,  0.0554, 0.0584, 0.0599, 0.0606,            &
           0.0609, 0.0611, 5*0.0612  /
        DATA (BTVDR (I, 2, 6), I = 1, 14)                                      &
          /0.0181, 0.0371, 0.0476, 0.0537,  0.0568, 0.0583, 0.0590,            &
           0.0593, 0.0595, 0.0595, 4*0.0596 /
	DATA (BTVDR (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 7), I = 1, 14) /14*0./
	DATA (BTVDR (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 8), I = 1, 14) /14*0./
	DATA (BTVDR (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTVDR (I, 2, 9), I = 1, 14) /14*0./

!****
!**** -----------------------------------------------------------
	DATA (GMVDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.0814, 0.1361, 0.2078, 0.2650, 0.2986, 0.3169,  0.3265,            &
        	   0.3313, 0.3337, 0.3348, 0.3354, 0.3357, 2*0.3358 /
	DATA (GMVDR (I, 2, 1), I = 1, 14)                                      &
       	  /0.0760, 0.1336, 0.2034, 0.2622, 0.2969, 0.3159,  0.3259,            &
       	   0.3309, 0.3333, 0.3346, 0.3352, 0.3354, 2*0.3356 /
	DATA (GMVDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.0834, 0.1252, 0.1558, 0.1927, 0.2131,   0.2237, 0.2290,           &
       	   0.2315, 0.2327, 0.2332, 0.2335, 2*0.2336, 0.2337 /
	DATA (GMVDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.0789, 0.1235, 0.1531, 0.1912, 0.2122, 0.2232,  0.2286,            &
      	   0.2312, 0.2324, 0.2330, 0.2333, 0.2334, 2*0.2335 /
	DATA (GMVDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.0647, 0.1342, 0.2215, 0.2968, 0.3432, 0.3696, 0.3838,             &
       	   0.3912, 0.3950, 0.3968, 0.3978, 0.3982, 0.3984, 0.3985 /
	DATA (GMVDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.0258, 0.1227, 0.1999, 0.2825, 0.3339, 0.3634, 0.3794,             &
       	   0.3877, 0.3919, 0.3940, 0.3950, 0.3956, 0.3958, 0.3959 /
	DATA (GMVDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.3371, 0.5762, 0.7159, 0.7927, 0.8324, 0.8526,  0.8624,            &
       	   0.8671, 0.8693, 0.8704, 0.8709, 0.8710, 2*0.8712 /
	DATA (GMVDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.2634, 0.4375, 0.5532, 0.6291, 0.6763, 0.7048, 0.7213,             &
       	   0.7310, 0.7363, 0.7395, 0.7411, 0.7420, 0.7426, 0.7428 /
        DATA (GMVDR (I, 1, 5), I = 1, 14)                                      &
           /0.0971, 0.1544, 0.2511, 0.3157, 0.3548, 0.3768, 0.3886,            &
            0.3948, 0.3978, 0.3994, 0.4001, 0.4006, 0.4007, 0.4008 /
        DATA (GMVDR (I, 2, 5), I = 1, 14)                                      &
           /0.0924, 0.1470, 0.2458, 0.3123, 0.3527, 0.3756, 0.3877,            &
            0.3942, 0.3974, 0.3990, 0.3998, 0.4002, 0.4004, 0.4005 /
        DATA (GMVDR (I, 1, 6), I = 1, 14)                                      &
           /0.0970, 0.1355, 0.1841, 0.2230, 0.2447,  0.2561, 0.2617,           &
            0.2645, 0.2658, 0.2664, 0.2667, 3*0.2669 /
        DATA (GMVDR (I, 2, 6), I = 1, 14)                                      &
           /0.0934, 0.1337, 0.1812, 0.2213, 0.2437, 0.2554, 0.2613,            &
            0.2642, 0.2656, 0.2662, 0.2665, 0.2667, 0.2667, 0.2668 /
	DATA (GMVDR (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 7), I = 1, 14) /14*1./
	DATA (GMVDR (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 8), I = 1, 14) /14*1./
	DATA (GMVDR (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMVDR (I, 2, 9), I = 1, 14) /14*1./

!****
!****  -----------------------------------------------------------

	DATA (ALIDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.2867,  0.2840, 0.2828, 0.2822, 0.2819, 0.2818, 2*0.2817,          &
       	   6*0.2816 /
	DATA (ALIDR (I, 2, 1), I = 1, 14)                                      &
        	  /0.3564, 0.3573, 0.3577, 0.3580, 2*0.3581, 8*0.3582 /
	DATA (ALIDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.2848, 0.2819, 0.2804, 0.2798, 0.2795, 2*0.2793, 7*0.2792 /
	DATA (ALIDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.3544, 0.3550, 0.3553, 2*0.3555, 9*0.3556 /
	DATA (ALIDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.2350, 0.2311, 0.2293, 0.2285, 0.2281, 0.2280, 8*0.2279 /
	DATA (ALIDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.2474, 0.2436, 0.2418, 0.2410, 0.2406, 0.2405, 3*0.2404,           &
       	   5*0.2403 /
	DATA (ALIDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.5816, 0.6157, 0.6391, 0.6556, 0.6673, 0.6758, 0.6820,             &
       	   0.6866, 0.6899, 0.6924, 0.6943, 0.6956, 0.6966, 0.6974 /
	DATA (ALIDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.5489, 0.5770, 0.5955, 0.6079, 0.6163, 0.6221, 0.6261,             &
       	   0.6288, 0.6308, 0.6321, 0.6330, 0.6337, 0.6341, 0.6344 /
        DATA (ALIDR (I, 1, 5), I = 1, 14)                                      &
           /0.2845, 0.2837, 0.2832, 0.2831, 0.2830, 9*0.2829 /
        DATA (ALIDR (I, 2, 5), I = 1, 14)                                      &
           /0.3532, 0.3562, 0.3578,  0.3586, 0.3590, 0.3592, 0.3594,           &
            0.3594, 0.3594, 5*0.3595 /
        DATA (ALIDR (I, 1, 6), I = 1, 14)                                      &
           /0.2825, 0.2812, 0.2806, 0.2803, 0.2802, 9*0.2801 /
        DATA (ALIDR (I, 2, 6), I = 1, 14)                                      &
           /0.3512, 0.3538,  0.3552, 0.3559, 0.3562, 0.3564, 0.3565,           &
            0.3565, 6*0.3566 /
	DATA (ALIDR (I, 1, 7), I = 1, 14) /14*ALIDRS/
	DATA (ALIDR (I, 2, 7), I = 1, 14) /14*ALIDRS/
	DATA (ALIDR (I, 1, 8), I = 1, 14) /14*ALIDRD/
	DATA (ALIDR (I, 2, 8), I = 1, 14) /14*ALIDRD/
	DATA (ALIDR (I, 1, 9), I = 1, 14) /14*ALIDRI/
	DATA (ALIDR (I, 2, 9), I = 1, 14) /14*ALIDRI/

!****
!**** -----------------------------------------------------------
	DATA (BTIDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.1291, 0.1707, 0.1969, 0.2125, 0.2216,   0.2267, 0.2295,           &
       	   0.2311, 0.2319, 0.2323, 0.2326, 2*0.2327, 0.2328 /
	DATA (BTIDR (I, 2, 1), I = 1, 14)                                      &
       	  /0.1939, 0.2357, 0.2598, 0.2735, 0.2810,  0.2851, 0.2874,            &
       	   0.2885, 0.2892, 0.2895, 0.2897, 3*0.2898 /
	DATA (BTIDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.1217, 0.1522, 0.1713, 0.1820,   0.1879,  0.1910, 0.1926,          &
      	   0.1935, 0.1939, 0.1942, 2*0.1943, 2*0.1944 /
	DATA (BTIDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.1781, 0.2067, 0.2221, 0.2301,   0.2342,  0.2363, 0.2374,          &
       	   0.2379, 0.2382, 0.2383, 2*0.2384, 2*0.2385 /
	DATA (BTIDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.0846, 0.1299, 0.1614, 0.1814, 0.1935,   0.2004, 0.2043,           &
           0.2064, 0.2076, 0.2082, 0.2085, 2*0.2087, 0.2088 /
	DATA (BTIDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.0950, 0.1410, 0.1722, 0.1921, 0.2042, 0.2111,  0.2151,            &
       	   0.2172, 0.2184, 0.2191, 0.2194, 0.2196, 2*0.2197 /
	DATA (BTIDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.5256, 0.7444, 0.9908, 1.2700, 1.5680, 1.8505, 2.0767,             &
       	   2.2211, 2.2808, 2.2774, 2.2362, 2.1779, 2.1160, 2.0564 /
	DATA (BTIDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.4843, 0.6714, 0.8577, 1.0335, 1.1812, 1.2858, 1.3458,             &
       	   1.3688, 1.3685, 1.3546, 1.3360, 1.3168, 1.2989, 1.2838 /
	DATA (BTIDR (I, 1, 5), I = 1, 14)                                      &
           /0.1498, 0.1930, 0.2201, 0.2364, 0.2460, 0.2514, 0.2544,            &
            0.2560, 0.2569, 0.2574, 0.2577, 0.2578, 0.2579, 0.2579 /
        DATA (BTIDR (I, 2, 5), I = 1, 14)                                      &
           /0.2184, 0.2656, 0.2927, 0.3078, 0.3159,  0.3202, 0.3224,           &
            0.3235, 0.3241, 0.3244, 0.3245, 3*0.3246 /
        DATA (BTIDR (I, 1, 6), I = 1, 14)                                      &
           /0.1369, 0.1681, 0.1860, 0.1958, 0.2010,  0.2038, 0.2053,           &
            0.2060, 0.2064, 0.2066, 0.2067, 3*0.2068 /
        DATA (BTIDR (I, 2, 6), I = 1, 14)                                      &
           /0.1969, 0.2268, 0.2416,  0.2488, 0.2521, 0.2537, 0.2544,           &
            0.2547, 0.2548, 5*0.2549 / 
	DATA (BTIDR (I, 1, 7), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 7), I = 1, 14) /14*0./
	DATA (BTIDR (I, 1, 8), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 8), I = 1, 14) /14*0./
	DATA (BTIDR (I, 1, 9), I = 1, 14) /14*0./
	DATA (BTIDR (I, 2, 9), I = 1, 14) /14*0./

!****
!**** --------------------------------------------------------------
	DATA (GMIDR (I, 1, 1), I = 1, 14)                                      &
       	  /0.1582, 0.2581, 0.3227, 0.3635, 0.3882, 0.4026, 0.4108,             &
       	   0.4154, 0.4179, 0.4193, 0.4200, 0.4204, 0.4206, 0.4207 /
	DATA (GMIDR (I, 2, 1), I = 1, 14)                                      &
       	  /0.1934, 0.3141, 0.3818, 0.4200, 0.4415, 0.4533, 0.4598,             &
       	   0.4633, 0.4651, 0.4662, 0.4667, 0.4671, 2*0.4672 /
	DATA (GMIDR (I, 1, 2), I = 1, 14)                                      &
       	  /0.1347, 0.1871, 0.2277, 0.2515, 0.2651, 0.2727, 0.2768,             &
       	   0.2790, 0.2801, 0.2808, 0.2811, 0.2812, 0.2813, 0.2814 /
	DATA (GMIDR (I, 2, 2), I = 1, 14)                                      &
       	  /0.1440, 0.2217, 0.2629, 0.2839, 0.2947, 0.3003, 0.3031,             &
       	   0.3046, 0.3054, 0.3058, 0.3060, 2*0.3061, 0.3062 /
	DATA (GMIDR (I, 1, 3), I = 1, 14)                                      &
       	  /0.1372, 0.2368, 0.3235, 0.3839, 0.4229, 0.4465, 0.4602,             &
       	   0.4679, 0.4722, 0.4745, 0.4758, 0.4764, 0.4768, 0.4770 /
	DATA (GMIDR (I, 2, 3), I = 1, 14)                                      &
       	  /0.1435, 0.2524, 0.3370, 0.3955, 0.4332, 0.4563, 0.4697,             &
       	   0.4773, 0.4815, 0.4839, 0.4851, 0.4858, 0.4861, 0.4863 /
	DATA (GMIDR (I, 1, 4), I = 1, 14)                                      &
       	  /0.4298, 0.9651, 1.6189, 2.4084, 3.2992, 4.1928, 4.9611,             &
       	   5.5095, 5.8085, 5.9069, 5.8726, 5.7674, 5.6346, 5.4944 /
	DATA (GMIDR (I, 2, 4), I = 1, 14)                                      &
       	  /0.4167, 0.8974, 1.4160, 1.9414, 2.4147, 2.7803, 3.0202,             &
      	   3.1468, 3.1954, 3.1932, 3.1676, 3.1328, 3.0958, 3.0625 /
        DATA (GMIDR (I, 1, 5), I = 1, 14)                                      &
           /0.1959, 0.3203, 0.3985, 0.4472, 0.4766, 0.4937, 0.5034,            &
            0.5088, 0.5117, 0.5134, 0.5143, 0.5147, 0.5150, 0.5152 /
        DATA (GMIDR (I, 2, 5), I = 1, 14)                                      &
           /0.2328, 0.3859, 0.4734, 0.5227, 0.5498, 0.5644, 0.5720,            &
            0.5761, 0.5781, 0.5792, 0.5797, 0.5800, 0.5802, 0.5802 /
        DATA (GMIDR (I, 1, 6), I = 1, 14)                                      &
           /0.1447, 0.2244, 0.2698, 0.2953, 0.3094, 0.3170, 0.3211,            &
            0.3233, 0.3244, 0.3250, 0.3253, 0.3255, 0.3256, 0.3256 /
        DATA (GMIDR (I, 2, 6), I = 1, 14)                                      &
           /0.1643, 0.2624, 0.3110, 0.3347, 0.3461, 0.3517, 0.3543,            &
            0.3556, 0.3562, 0.3564, 0.3565, 0.3566, 0.3566, 0.3566 /
	DATA (GMIDR (I, 1, 7), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 7), I = 1, 14) /14*1./
	DATA (GMIDR (I, 1, 8), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 8), I = 1, 14) /14*1./
	DATA (GMIDR (I, 1, 9), I = 1, 14) /14*1./
	DATA (GMIDR (I, 2, 9), I = 1, 14) /14*1./

!**** -----------------------------------------------------------


!FPP$ EXPAND (clsmf25_COEFFSIB)

      DO I=1,NCH
        ALA = AMIN1 (AMAX1 (ZERO, VLAI(I)), ALATRM)
        LAI = 1 + MAX(0, INT((ALA-BLAI)/DLAI) )
        DX = (ALA - (BLAI+(LAI-1)*DLAI)) * (ONE/DLAI)
        DY = (VGRN(I)- GRN(1)) * (ONE/(GRN(2) - GRN(1)))

        ALPHA = clsmf25_COEFFSIB (ALVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        BETA  = clsmf25_COEFFSIB (BTVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        GAMMA = clsmf25_COEFFSIB (GMVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

        GAMMA = MAX(GAMMA,0.01)

        AVISDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
        AVISDF(I) = ALPHA-BETA                                                 &
                 + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))

        ALPHA = clsmf25_COEFFSIB (ALIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        BETA  = clsmf25_COEFFSIB (BTIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
        GAMMA = clsmf25_COEFFSIB (GMIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

        GAMMA = MAX(GAMMA,0.01)

        ANIRDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
        ANIRDF(I) = ALPHA-BETA                                                 &
                 + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))

! SCALE TO MODIS VALUES (SNOW-FREE)

	  AVISDR(I) = AVISDR(I) * SCALVDR(I)
          ANIRDR(I) = ANIRDR(I) * SCALIDR(I)
	  AVISDF(I) = AVISDF(I) * SCALVDF(I)
          ANIRDF(I) = ANIRDF(I) * SCALIDF(I)

! PROTECT AGAINST BAD SCALING

	  AVISDR(I) = AMIN1( 1., AMAX1( 0., AVISDR(I) ) )
          ANIRDR(I) = AMIN1( 1., AMAX1( 0., ANIRDR(I) ) )
	  AVISDF(I) = AMIN1( 1., AMAX1( 0., AVISDF(I) ) )
          ANIRDF(I) = AMIN1( 1., AMAX1( 0., ANIRDF(I) ) )

! SNOW ALBEDOES

        SWE=WESN(1,I)+WESN(2,I)+WESN(3,I)
        TOTDEP=SNDZ(1,I)+SNDZ(2,I)+SNDZ(3,I)
        AREASC = MIN(SWE/WEMIN,1.)
        DENSITY=(SWE/(AREASC+1.e-20)) / (TOTDEP+1.e-20)
        DENS_EXC=MAX(0., DENSITY-RHOFS)

        ASNVDR(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX - DENS_EXC*.0006)
        ASNNDR(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX - DENS_EXC*.0006)
        ASNVDF(I) = MAX(SNWALB_VISMIN, SNWALB_VISMAX - DENS_EXC*.0006)
        ASNNDF(I) = MAX(SNWALB_NIRMIN, SNWALB_NIRMAX - DENS_EXC*.0006)

! ACCOUNT FOR VEGETATION MASKING, FOR EACH COMPONENT

! A) FIRST DO MASKING IN VEGETATED FRACTION:
        FAC = SWE / (SWE + SNWMID(ITYP(I)))
        ASNVDR_VEG=AVISDR(I) + (ASNVDR(I)-AVISDR(I))*FAC
        ASNNDR_VEG=ANIRDR(I) + (ASNNDR(I)-ANIRDR(I))*FAC
        ASNVDF_VEG=AVISDF(I) + (ASNVDF(I)-AVISDF(I))*FAC
        ASNNDF_VEG=ANIRDF(I) + (ASNNDF(I)-ANIRDF(I))*FAC

! B) NOW ACCOUNT FOR SUBGRID VEGETATION FRACTION
        FVEG=AMIN1( 1., VLAI(I)/2. )
        ASNVDR(I)=ASNVDR(I)*(1.-FVEG)+ASNVDR_VEG*FVEG
        ASNNDR(I)=ASNNDR(I)*(1.-FVEG)+ASNNDR_VEG*FVEG
        ASNVDF(I)=ASNVDF(I)*(1.-FVEG)+ASNVDF_VEG*FVEG
        ASNNDF(I)=ASNNDF(I)*(1.-FVEG)+ASNNDF_VEG*FVEG


!        AMASK=FVEG*EXP(-TOTDEP/SNWMSK(ITYP(I)))

!
!	  AVISDR(I) = AVISDR(I) + (SNWALB(1,ITYP(I)) - AVISDR(I)) * FAC
!          ANIRDR(I) = ANIRDR(I) + (SNWALB(2,ITYP(I)) - ANIRDR(I)) * FAC
!	  AVISDF(I) = AVISDF(I) + (SNWALB(3,ITYP(I)) - AVISDF(I)) * FAC
!          ANIRDF(I) = ANIRDF(I) + (SNWALB(4,ITYP(I)) - ANIRDF(I)) * FAC
!	  ENDIF

!  (ORIGINAL NOTES FROM STIEGLITZ:)

!        ALBSNW = 0.913 - .0006*DENSITY
!        ALBSNW = AMIN1(1., AMAX1(ALBSNW,0.5))

!  a) formulation
!  While 3 basic land covers are allowed (bare soil, snow, and vegetation),
!  the grid cell is divided into 2 fractions: a vegetated fraction (A_v),
!  and a non-vegetated fraction (A_u).  In the absence of snow cover
!  A_v = m A_v0
!  A_u = 1 - A_v
!  where A_v0 is defined to be the grid cell specified value for the
!  vegetated fraction, A_u is the resultant bare soil fraction, and m, the
!  snow masking fraction, is unity.  In the absence of snow A_u is simply
!  the unvegetated portion of the grid cell.  In the presence of snow, the
!  green vegetation masked by snow is
!  m = exp( - d_s / d_m)
!  where d_s is the snow depth(not water equivalent) and d_m is the
!  vegetation specific masking depth.  A_u now becomes that portion of the
!  grid cell where vegetation is not visible
!  b) masking depths
!  DATA SMK/0.1,0.2,0.2,0.5,2.0,5.0,10.0,25.0/
!              SMK array types - desert, tundra, grass, shrub, woodland,
!  deciduous, evergreen, rain forest


!        FRACV=0.5
!        AMASK=FRACV*EXP(-TOTDEP/SNWMSK)
!        ALBSNW=ALBSNW*(1.-AMASK)+ALBAVE*AMASK
        
!        IF (SNW (I) .GT. ZERO) THEN
!	  FAC = SNW(I) / (SNW(I) + SNWMID(ITYP(I)))
!
!	  AVISDR(I) = AVISDR(I) + (SNWALB(1,ITYP(I)) - AVISDR(I)) * FAC
!          ANIRDR(I) = ANIRDR(I) + (SNWALB(2,ITYP(I)) - ANIRDR(I)) * FAC
!	  AVISDF(I) = AVISDF(I) + (SNWALB(3,ITYP(I)) - AVISDF(I)) * FAC
!          ANIRDF(I) = ANIRDF(I) + (SNWALB(4,ITYP(I)) - ANIRDF(I)) * FAC
!	  ENDIF

        ENDDO

      RETURN
      END SUBROUTINE SIBALB



  subroutine calc_soil_moist( &
       ncat,vegcls,dzsf,vgwmax,cdcr1,cdcr2,wpwet,poros, &
       psis,bee,ars1,ars2,ars3,ara1,ara2, &
       ara3,ara4,arw1,arw2,arw3,arw4, &
       srfexc,rzexc,catdef, &
       sfmc, rzmc, prmc,  &
       werror, sfmcun, rzmcun, prmcun )
    
    ! Calculate diagnostic soil moisture content from prognostic
    ! excess/deficit variables.
    !
    ! On input, also check validity of prognostic excess/deficit variables
    ! and modify if necessary.  Perturbed or updated excess/deficit variables 
    ! in data assimilation integrations may be unphysical.  
    ! Optional output "werror" contains excess or missing water related
    ! to inconsistency.
    !
    ! Optional outputs "smfcun", "rzmcun", "prmcun" are surface,
    ! root zone, and profile moisture content for unsaturated areas only,
    ! ie. excluding the saturated area of the catchment.    
    !
    ! NOTE: When calling with optional output arguments, use keywords
    !       unless arguments are in proper order!
    !       
    !       Example: 
    !       (don't want "werror" as output, but want "*mcun" output)
    !       
    !       call calc_soil_moist(         & 
    !            ncat, ...                &
    !            sfmc, rzmc, prmc,        &
    !            sfmcun=sfmc_unsat,   &
    !            rzmcun=rzmc_unsat,   & 
    !            prmcun=prmc_unsat )
    !
    ! replaces moisture_sep_22_2003.f (and older moisture.f)
    !
    ! koster+reichle, Feb 5, 2004
    !
    ! revised - koster+reichle, Mar 19, 2004
    !
    ! added optional *un output - koster+reichle, Apr 6, 2004
    !
    ! ----------------------------------------------------------------

    
    implicit none
    
    integer, parameter                   :: KSNGL=4
    integer,                  intent(in) :: ncat
    integer, dimension(ncat), intent(in) :: vegcls
    
    real(KIND=KSNGL), dimension(ncat), intent(in) :: dzsf,vgwmax,cdcr1,cdcr2
    real(KIND=KSNGL), dimension(ncat), intent(in) :: wpwet,poros,psis
    real(KIND=KSNGL), dimension(ncat), intent(in) :: bee,ars1
    real(KIND=KSNGL), dimension(ncat), intent(in) :: ars2,ars3,ara1,ara2,ara3
    real(KIND=KSNGL), dimension(ncat), intent(in) :: ara4,arw1,arw2,arw3,arw4
    
    real(KIND=KSNGL), dimension(ncat), intent(inout) :: srfexc, rzexc, catdef
    
    real(KIND=KSNGL), dimension(ncat), intent(out) :: sfmc, rzmc, prmc
    
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: werror
    
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: sfmcun
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: rzmcun
    real(KIND=KSNGL), dimension(ncat), intent(out), optional :: prmcun
    
    ! ----------------------------
    !    
    ! local variables
    
    integer :: n
    
    real(KIND=KSNGL), parameter :: dtstep_dummy = -9999.
    
    real(KIND=KSNGL), dimension(ncat) :: rzeq, runsrf_dummy, catdef_dummy
    real(KIND=KSNGL), dimension(ncat) :: ar1, ar2, ar4, prmc_orig
    real(KIND=KSNGL), dimension(ncat) :: srfmn, srfmx, swsrf1, swsrf2, swsrf4, rzi
    

    ! --------------------------------------------------------------------
    !
    ! compute soil water storage upon input [mm]
    
    do n=1,ncat
       prmc_orig(n) =                                                 &
            (cdcr2(n)/(1.-wpwet(n))-catdef(n)+rzexc(n)+srfexc(n))
    enddo
       
    ! -----------------------------------
    !
    ! check limits of catchment deficit
    !
    ! increased minimum catchment deficit from 0.01 to 1. to make the 
    ! check work with perturbed parameters and initial condition
    ! reichle, 16 May 01
    !
    ! IT REALLY SHOULD WORK WITH catdef > 0 (rather than >1.) ????
    ! reichle, 5 Feb 2004 
    
    do n=1,ncat     
       catdef(n)=max(1.,min(cdcr2(n),catdef(n)))
    end do
    
    ! ------------------------------------------------------------------
    !     
    ! check limits of root zone excess
    !     
    ! calculate root zone equilibrium moisture for given catchment deficit
    
    call rzequil( &
         ncat, vegcls, catdef, vgwmax, &
         cdcr1, cdcr2, wpwet, &
         ars1, ars2, ars3, ara1, ara2, ara3, ara4, &
         arw1, arw2, arw3, arw4, &
         rzeq)
    
    ! assume srfexc=0 and constrain rzexc appropriately
    ! (iteration would be needed to contrain srfexc and rzexc simultaneously)
    
    do n=1,ncat
       rzexc(n)=max(wpwet(n)*vgwmax(n)-rzeq(n),min(vgwmax(n)-rzeq(n),rzexc(n)))
    end do
    
    ! this translates into:
    !
    ! wilting level < rzmc < porosity
    ! 
    ! or more precisely:  wpwet*vgwmax < rzeq+rzexc < vgwmax
    ! 
    ! NOTE: root zone moisture is not allowed to drop below wilting level 
    
    ! -----------------------------------------------------------------
    !
    ! Call partition() for computation of surface moisture content.
    !
    ! Call to partition() also checks limits of surface excess.
    !
    ! Call partition with dtstep_dummy:
    !  In partition, dtstep is only used for a correction that
    !  puts water into runsrf (for which runsrf_dummy is used here).
    !  Also use catdef_dummy because partition() updates catdef
    !  whenever srfexc exceeds physical bounds, but this is not desired here.
    
    runsrf_dummy = 0.
    catdef_dummy = catdef          
    
    call partition( &
         ncat,dtstep_dummy,vegcls,dzsf,rzexc, &
         rzeq,vgwmax,cdcr1,cdcr2, &
         psis,bee,poros,wpwet, &
         ars1,ars2,ars3, &
         ara1,ara2,ara3,ara4, &
         arw1,arw2,arw3,arw4,.false., &
         srfexc,catdef_dummy,runsrf_dummy, &
         ar1, ar2, ar4,srfmx,srfmn, & 
         swsrf1,swsrf2,swsrf4,rzi &
         )
    
    ! compute surface, root zone, and profile soil moisture
    
    do n=1,ncat

       sfmc(n) = poros(n) *                                           &
            (swsrf1(n)*ar1(n) + swsrf2(n)*ar2(n) + swsrf4(n)*ar4(n))
       
       rzmc(n) = (rzeq(n)+rzexc(n)+srfexc(n))*poros(n)/vgwmax(n)
       
       ! compute revised soil water storage [mm]
       
       prmc(n) =                                                               &
            (cdcr2(n)/(1.-wpwet(n))-catdef(n)+rzexc(n)+srfexc(n))

       ! compute error in soil water storage [mm] (if argument is present)
       
       if (present(werror))  werror(n)=(prmc(n)-prmc_orig(n))
       
       ! convert to volumetric soil moisture
       ! note: dzpr = (cdcr2/(1-wpwet)) / poros 
       
       prmc(n) = prmc(n)*poros(n) / (cdcr2(n)/(1.-wpwet(n)))

       
       ! check for negative soil moisture 
       
       if ( (sfmc(n)<.0) .or. (rzmc(n)<.0) .or. (prmc(n)<.0) ) then
          
          write (*,*) 'FOUND NEGATIVE SOIL MOISTURE CONTENT.... stopping'
          write (*,*) n, sfmc(n), rzmc(n), prmc(n)
          stop
       end if

       ! compute moisture content in unsaturated areas [m3/m3] (if arg present)

       if (ar1(n)<1.) then       

          if (present(prmcun))  prmcun(n)=(prmc(n)-poros(n)*ar1(n))/(1.-ar1(n))
          if (present(rzmcun))  rzmcun(n)=(rzmc(n)-poros(n)*ar1(n))/(1.-ar1(n))
          if (present(sfmcun))  sfmcun(n)=(sfmc(n)-poros(n)*ar1(n))/(1.-ar1(n))

       else          
          
          if (present(prmcun))  prmcun(n)=poros(n)
          if (present(rzmcun))  rzmcun(n)=poros(n)
          if (present(sfmcun))  sfmcun(n)=poros(n)
          
       end if
       
    enddo

  return
    
  end subroutine calc_soil_moist

  ! ********************************************************************
    
      END MODULE clsmf25_model
