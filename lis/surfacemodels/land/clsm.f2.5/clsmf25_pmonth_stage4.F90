!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: clsmf25_pmonth
!
! !DESCRIPTION:
!
!
! !REVISION HISTORY:
! 23 Nov 2012: David Mocko, Added Catchment Fortuna-2.5
!
! !INTERFACE:
SUBROUTINE clsmf25_PMONTH (&
       NCH,    ITYP,   JDAY,   ALAT,&
       GREEN, ZLT, Z0, D, BUG,&
       SQSCAT, Z2, DZM,&
       RSOIL1, RSOIL2, SATCAP, RDC,  U2FAC &
       &                  )
!EOP
    ! reichle+qliu,  8 Oct 2008 - updated to match GEOS5-MERRA parameters
    ! reichle       10 Feb 2011 - moved here from file pmonth_stage4.F90
    
    IMPLICIT NONE
    
    !         This subroutine sets seasonally varying
    !      variables.  NOTE: LAI is stored as LAI (SIB) / VCOVER (SIB).
    !    
    INTEGER MONTHs, ntyps
    PARAMETER (MONTHS = 13)
    parameter (NTYPS = 10)
    !    
    LOGICAL BUG
    INTEGER NCH
    INTEGER ITYP(NCH), VEG, JDAY
    !    
    !     IntegerTYPe:  Vegetation types (biomes).
    !         1:  BrdEvr:  Broadleaf evergreen.
    !         2:  BrdDcd:  Broadleaf deciduous.
    !         3:  Needle:  Needleleaf trees.
    !         4:  GndCvr:  Groundcover (grass).
    !    
    REAL GREEN(NCH),  ZLT(NCH),    SQSCAT(NCH), Z0(NCH),&
         Z2(NCH),     D(NCH),      DZM(NCH),    RSOIL1(NCH),&
         RSOIL2(NCH), SATCAP(NCH), RDC(NCH),    U2FAC(NCH),&
         ALAT(NCH)
    !    
    INTEGER ChNo, K, KDAY, KLOW
    !    
    REAL VGRDC(MONTHS,NTYPS),  VGROTL(MONTHS,NTYPS),&
         VGZ0(MONTHS,NTYPS),   VGDD(MONTHS,NTYPS),&
         VGRF11(NTYPS),         VGRF12(NTYPS),&
         VGTR11(NTYPS),         VGTR12(NTYPS),&
         VGROCA(NTYPS),         VGROTD(NTYPS),&
         VGRDRS(NTYPS),         VGZ2(NTYPS),  &
         VGRT(NTYPS)
    
    !    
    REAL AKDAY0(MONTHS), ALPHAF, CUNI, DUM1, DUM2, ROOTL, SCAT,&
         VKC, VROOT, WHI, WLO, RLENGT, DPARAM
    !    
    !     --------------------------------------------------------------------
    !    
    !----- Sarith 4/19/06 -- from GEOScatchGridComp.F90
    !
    REAL, DIMENSION (NTYPS-2) :: VGRDA, VGRDB
    
    ! Correction to RDC formulation -Randy Koster, 4/1/2011
    !
    !data VGRDA / 285.9, 294.9, 652.9,  25.8,  100.7,   &
    !     22.9,  23.8, 23.8/
    !data VGRDB / 5.1 ,  7.2, 10.8,  4.8,  1.8,  5.1,  .000, .000/
    !
    data VGRDA / 285.9, 355.18, 660.24,  30.06,  100.7,  24.36,  23.8, 23.8/
    data VGRDB / 5.1 ,  7.2, 10.5,  4.8,  1.8,  5.1,  .000, .000/

    DATA VGRDC/285.87,285.87,285.87,285.87,285.87,285.87,285.87,&
         285.87, 285.87, 285.87, 285.87, 285.87, 285.87, &
         211.32, 211.32, 218.78, 243.40, 294.87, 345.90, 355.18,&
         341.84, 307.22, 244.84, 218.78, 211.32, 211.32, &
         565.41, 587.05, 623.46, 638.13, 652.86, 675.04, 660.24, &
         645.49, 638.13, 623.46, 587.05, 565.41, 565.41, &
         24.43,  24.63,  24.80,  24.96,  25.72,  27.74,  30.06,&
         28.86,  25.90,  25.11,  24.80,  24.63,  24.43, &
         103.60, 103.60, 102.35, 100.72, 100.72, 100.72, 100.72,&
         105.30, 107.94, 106.59, 104.49, 103.60, 103.60, &
         22.86,  22.86,  22.86,  22.86,  22.86,  23.01,  24.36,&
         24.69,  24.04,  22.86,  22.86,  22.86,  22.86, &
         23.76,  23.76,  23.76,  23.76,  23.76,  23.76,  23.76,&
         23.76,  23.76,  23.76,  23.76,  23.76,  23.76, &
         23.76,  23.76,  23.76,  23.76,  23.76,  23.76,  23.76,&
         23.76,  23.76,  23.76,  23.76,  23.76,  23.76, &
         23.76,  23.76,  23.76,  23.76,  23.76,  23.76,  23.76,&
         23.76,  23.76,  23.76,  23.76,  23.76,  23.76, &
         23.76,  23.76,  23.76,  23.76,  23.76,  23.76,  23.76,&
         23.76,  23.76,  23.76,  23.76,  23.76,  23.76/
    
    DATA VGROTL/19737.8,19737.8,19737.8,19737.8,19737.8,19737.8,&
         19737.8, 19737.8, 19737.8, 19737.8, 19737.8, 19737.8, 19737.8,&
         5010.0,  5010.0,  5270.0,  6200.0,  8000.0,  9700.0,  9500.0,&
         8400.0,  6250.0,  5270.0,  5010.0,  5010.0,  5010.0,&
         9000.0,  9200.0,  9533.3,  9666.7,  9800.0,  9866.7,  9733.3, &
         9666.7,  9533.3,  9200.0,  9000.0,  9000.0,  9000.0,&
         5500.0,  5625.0,  5750.0,  5875.0,  6625.0,  8750.0,  9375.0,&
         6875.0,  6000.0,  5750.0,  5625.0,  5500.0,  5500.0,&
         6500.0,  6000.0,  5500.0,  5500.0,  5500.0,  5500.0,  5500.0,&
         7500.0,  8500.0,  7000.0,  6500.0,  6500.0,  6500.0,&
         10625.0, 10625.0, 10625.0, 10625.0, 10625.0, 11250.0, 18750.0,&
         17500.0, 10625.0, 10625.0, 10625.0, 10625.0, 10625.0,&
         1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,&
         1.0,     1.0,     1.0,     1.0,     1.0,     1.0,&
         1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,&
         1.0,     1.0,     1.0,     1.0,     1.0,     1.0,&
         1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,&
         1.0,     1.0,     1.0,     1.0,     1.0,     1.0,&
         1.0,     1.0,     1.0,     1.0,     1.0,     1.0,     1.0,&
         1.0,     1.0,     1.0,     1.0,     1.0,     1.0/
      
    DATA VGZ0/2.6530, 2.6530, 2.6530, 2.6530, 2.6530, 2.6530, 2.6530, &
         2.6530, 2.6530, 2.6530, 2.6530, 2.6530, 2.6530,  &
         0.5200, 0.5200, 0.6660, 0.9100, 1.0310, 1.0440, 1.0420,&
         1.0370, 1.0360, 0.9170, 0.6660, 0.5200, 0.5200, &
         1.1120, 1.1030, 1.0880, 1.0820, 1.0760, 1.0680, 1.0730,&
         1.0790, 1.0820, 1.0880, 1.1030, 1.1120, 1.1120, &
         0.0777, 0.0778, 0.0778, 0.0779, 0.0778, 0.0771, 0.0759,&
         0.0766, 0.0778, 0.0779, 0.0778, 0.0778, 0.0777, &
         0.2450, 0.2450, 0.2270, 0.2000, 0.2000, 0.2000, 0.2000,&
         0.267,  0.292,  0.280,  0.258,  0.2450, 0.2450,  &
         0.0752, 0.0752, 0.0752, 0.0752, 0.0752, 0.0757, 0.0777,&
         0.0778, 0.0774, 0.0752, 0.0752, 0.0752, 0.0752, &
         0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, &
         0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, &
         0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, &
         0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, &
         0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, &
         0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112,&
         0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112, &
         0.0112, 0.0112, 0.0112, 0.0112, 0.0112, 0.0112/
    
    DATA VGDD/27.37, 27.37, 27.37, 27.37, 27.37, 27.37, 27.37,&
         27.37,   27.37,   27.37,   27.37,   27.37,   27.37, &
         13.66,   13.66,   14.62,   15.70,   16.33,   16.62,   16.66,&
         16.60,   16.41,   15.73,   14.62,   13.66,   13.66, &
         13.76,   13.80,   13.86,   13.88,   13.90,   13.93,   13.91,&
         13.89,   13.88,   13.86,   13.80,   13.76,   13.76, &
         0.218,   0.227,   0.233,   0.239,   0.260,   0.299,   0.325,&
         0.313,   0.265,   0.244,   0.233,   0.227,   0.218,  &
         2.813,   2.813,   2.662,   2.391,   2.391,   2.391,   2.391,&
         2.975,   3.138,   3.062,   2.907,   2.813,   2.813,  &
         0.10629, 0.10629, 0.10629, 0.10629, 0.10629, 0.12299, 0.21521,&
         0.22897, 0.19961, 0.10629, 0.10629, 0.10629, 0.10629,  &
         0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,&
         0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001, &
         0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,&
         0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001, &
         0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,&
         0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,&
         0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001,&
         0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001/
    
    DATA VGRF11 /0.10, 0.10, 0.07, 0.105 & 
         , 0.10, 0.10, .001, .001, .001, .001/
    
    DATA VGRF12 /0.16, 0.16, 0.16, 0.360&
         , 0.16, 0.16, .001, .001, .001, .001/
    
    DATA VGTR11 /0.05, 0.05, 0.05, 0.070&
         , 0.05, 0.05, .001, .001, .001, .001/
    
    DATA VGTR12 /.001, .001, .001,  .220&
         , .001, .001, .001, .001, .001, .001/
    
    DATA VGROCA / &
         0.384E-6, 0.384E-6, 0.384E-6, 0.384E-6, 0.384E-6, 0.384E-6,&
         .1E-6, .1E-6, .1E-6, .1E-6  / 
    
    DATA VGROTD /1.00, 1.00, 0.50, 0.50, 0.50, 0.20, 0.10, 0.10,&
         0.10, 0.10/ 
    
    DATA VGRDRS  /&
         0.75E13, 0.75E13, 0.75E13, 0.40E13, 0.75E13, 0.75E13,&
         0.10E13, 0.10E13, 0.10E13, 0.10E13  /
    
    ! use GEOS5/MERRA values for VGZ2, reichle+qliu, 8 Oct 2008
    ! DATA VGZ2 /35.0, 20.0, 17.0, 0.6, 5.0, 0.6, 0.1, 0.1, 0.1, 0.1/
    DATA VGZ2 /35.0, 20.0, 17.0, 0.6, 0.5, 0.6, 0.01, 0.01, 0.01, 0.01/ ! Dorman and Sellers (1989)
    
    DATA VGRT  / 19700., 7000., 9400., 7000., 7000., 14000., 1., 1., 1., 1./
    
    DATA AKDAY0 / 16.,  45.,  75., 105., 136., 166., 197., 228., 258.,&
         289., 319., 350., 381./
    
    DATA VKC /0.41/
    
    ! the following parameters must be consistent with: subroutine clsmf25()
    real,    parameter :: MIN_VEG_HEIGHT = 0.01           
    real,    parameter :: Z0_BY_ZVEG  = 0.13 
    
    DO ChNo = 1, NCH

       VEG = ITYP(ChNo)

       if(VEG.ne.9) then 
          !    
          !     DETERMINE WEIGHTING FACTOR:
          IF (ALAT(CHNO) .LT. 0.) THEN
             KDAY = MOD (182 + JDAY, 365)
          ELSE
             KDAY = MOD (365 + JDAY, 365)
          ENDIF
          
          !    
          KLOW = 12
          
          DO K = 1, 12
             IF (KDAY .GE. AKDAY0 (K)) KLOW = K
          END DO
          
          IF (KLOW .EQ. 12 .AND. KDAY .LT. AKDAY0 (12)) THEN
             WLO = (16. - KDAY) / 31.
          ELSE
             WLO = (AKDAY0 (KLOW + 1) - KDAY) /&
                  (AKDAY0 (KLOW + 1) - AKDAY0 (KLOW))
          ENDIF
          
          WHI = 1 - WLO
          
          !    
          !     COMPUTE SOME PARAMETERS DIRECTLY:
          !  LAI and type dependent parameters; RDC formulation can use veg frac in next version. see the last line in the do loop - Sarith
          !         RDC (ChNo) = WLO * VGRDC (KLOW, ITYP (ChNo)) +&
          !              WHI * VGRDC (KLOW + 1, ITYP (ChNo))
          
          ! use GEOS5/MERRA formulation for ROOTL, reichle+qliu, 8 Oct 2008
          !ROOTL = WLO * VGROTL (KLOW, ITYP (ChNo)) +&
          !     WHI * VGROTL (KLOW + 1, ITYP (ChNo))
          ROOTL = VGRT(VEG)
          RLENGT=WLO*VGZ0(KLOW,VEG) +&
               WHI*VGZ0(KLOW+1,VEG)
          DPARAM=WLO*VGDD(KLOW,VEG) +&
               WHI*VGDD(KLOW+1,VEG)
          
          VROOT = ROOTL * VGROCA (VEG)
          DUM1 = ALOG (VROOT / (1 - VROOT))
          DUM2 = 1. / (8. * 3.14159 * ROOTL)
          ALPHAF = DUM2 * (VROOT - 3. -2. * DUM1)
          RSOIL1 (ChNo) =&
               VGRDRS (VEG) / (ROOTL * VGROTD (VEG))
          RSOIL2 (ChNo) = ALPHAF / VGROTD (VEG)
          !    
          SCAT = GREEN (ChNo) *&
               (VGTR11 (VEG) + VGRF11 (VEG)) +&
               (1. - GREEN (ChNo)) * &
               (VGTR12 (VEG) + VGRF12 (VEG))
          SQSCAT (ChNo) = SQRT (1. - SCAT)
          
          IF ( (BUG .EQV. .TRUE.)  .AND. (ChNo .eq. 1)) THEN
             WRITE(*,*) 'PMON: GREEN(1)=', GREEN(ChNo)
             WRITE(*,*) 'PMON: SCAT=',SCAT
             WRITE(*,*) 'PMON: SQSCAT=',SQSCAT(ChNo)
          ENDIF
          
          !    
          DZM (ChNo) = 600.
          
          !      
          CUNI = (ALOG (0.025 * DZM (ChNo) / RLENGT) / VKC) + 8.4
          
          ! U2FAC not used in LDASsa for now, reichle+qliu, 8 Oct 2008
          !
          !U2FAC (ChNo) =&
          !     ALOG ((VGZ2 (ITYP (ChNo)) - DPARAM) / RLENGT) / &
          !     (CUNI * VKC)
          !
          U2FAC(ChNo) = -9999.
          
          Z2 (ChNo) = VGZ2 (VEG)
          !    
          SATCAP (ChNo) = 0.2 * ZLT (ChNo)
          ! Time varying roughness length as in the AGCM -Sarith 4/18/06
          Z0 (ChNo) =  Z0_BY_ZVEG*(Z2 (ChNo)- (Z2 (ChNo)-MIN_VEG_HEIGHT)*exp(-ZLT (ChNo)))


          ! Correction to RDC formulation -Randy Koster, 4/1/2011
          ! bad:       RDC = max(VGRDA(VEG)*min(VGRDB(VEG),LAI           ),0.001)
          ! correct:   RDC = max(VGRDA(VEG)*min(1.,        LAI/VGRDB(VEG)),0.001)

          !RDC(ChNo)   = max(VGRDA(ITYP (ChNo))*min(VGRDB(ITYP (ChNo)),       &
          !     ZLT (ChNo)),0.001)
          
          RDC(ChNo) = max( VGRDA(VEG)*min( 1., ZLT(ChNo)/VGRDB(VEG) ), 0.001)
          
       endif
    ENDDO
    
    !    
    RETURN
  END SUBROUTINE clsmf25_PMONTH
