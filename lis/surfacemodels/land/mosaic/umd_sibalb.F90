!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! umd_sibalb.f90
!
! DESCRIPTION:
!  To use the MOSAIC-SiB albedo calculation, 7 SiB vegetation types were
!  mapped to 13 UMD vegetation types and 4 coefficient arrays were defined.
!
!  Some old comments from calc_albedo.f:
   !  This is a subroutine for the LDAS driver
   !  i've commented out the program part of the this file.
   !  because the driver only needs the subroutine.
   !  -Jared
   !
   !
   !   	This program computes 4 albedos given :
   !    LAI, Greenness fraction, Cosine of the Zenith Angle,
    !    Snow depth (set to 0 now since snow-related albedo 
   !    calculations have been disabled in this subroutine and 
   !    are taken care of in the main code of Mosaic), Veg type,
   !    IRUN (set to 1) and Canopy Temp. (set to 0 b/c this 
   !    is only needed for snow related calculations that have
   !    been disabled in this subroutine 
   !
   !    Current version of mosaic uses avisdr,anirdr,avisdf,anirdf 
   !      that are all set to the albedo taken from EDAS data
   !    Using this subroutine, these four parameters can be
   !      individually computed and passed back to the main model.
   !	This subroutine is a modified version of the one written
   !	  by Randy Koster / GSFC---it was modified 
   !      by Brian Cosgrove / GSFC
!
!===========================================================================
! REVISION HISTORY:
!  04 Apr 2001: Urszula Jambor; Initial code, using old calc_albedo.f scheme
!  11 Feb 2002: Jon Gottschalck; Added check for low total LAI values for
!               because it caused GAMMA below to become very low and crash run
!               This check is needed because of the usage of satellite LAI data
!               Low LAI values less than 0.5 were not a part of the original 
!               Mosaic tabular values
!  07 August 2006: Jon Gottschalck; Modified check for low total LAI values
!               from 0.15 to 0.50. It was discovered as runs increased to 1 km
!               resolution. It was decided to limit the calculation of
!               the albedo parameters by using LAI values not less than 0.5.
!               The assumption used is that the albedo parameters calculated 
!               will not change significantly from LAI values of 0 to 0.5. 
!               This correction ONLY applies in the albedo parameters calculation, 
!               the original LAI value is used throughout the rest of the model.
!               TEMPLAI replaces VLAI in this subroutine to accomplish this.
!===========================================================================

subroutine umd_sibalb ( AVISDR, ANIRDR, AVISDF, ANIRDF,        &
                        VLAI, VGRN, ZTH, SNW, ITYP, IRUN, TC )

!***************************************************************************
!* OUTPUTS:*****************************************************************
!* AVISDR:   visible, direct albedo.                                       *
!* ANIRDR:   near infra-red, direct albedo.                                *
!* AVISDF:   visible, diffuse albedo.                                      *
!* ANIRDF:   near infra-red, diffuse albedo.                               *
!***************************************************************************
!* INPUTS:******************************************************************
!* VLAI:     the leaf area index.                                          *
!* VGRN:     the greenness index.                                          *
!* ZTH:      The cosine of the solar zenith angle.                         *
!* SNW:      Snow cover in meters water equivalent.                        *
!* ITYP:     Vegetation class                                              *
!* IRUN:     Set to 1                                                      *
!* TC:       Temperature of Canopy (confusion as to if it's used)          *
!* sib:      SiB albedo coefficients, mapped to UMD veg types              *
!***************************************************************************

  use sibalb_module    ! SiB-based coefficients for albedo calculation.
  IMPLICIT NONE
!  type (sibalbdec)     sib

!===  Begin definition and assignment of variables: =======================

  integer :: IRUN
  integer :: ITYP(IRUN)
  integer :: I, J, LAI
  integer, parameter :: NLAI = 14
  integer, parameter :: NTYPS = 9
  integer, parameter :: UMDNTYPS = 13     

  real, parameter :: EPSLN = 1.E-6
  real, parameter :: BLAI = 0.5
  real, parameter :: DLAI = 0.5
!======== ALATRM = (BLAI + (NLAI - 1) * DLAI - EPSLN)
  real, parameter :: ALATRM = ( 0.5 + ( 14 - 1 ) * 0.5 - 1.E-6 )

  real, parameter :: ZERO = 0.0
  real, parameter :: ONE  = 1.0
  
  real, dimension(IRUN) :: AVISDR, ANIRDR, AVISDF, ANIRDF
  real, dimension(IRUN) :: VLAI, VGRN, ZTH, SNW, TC, TEMPLAI
  
  real :: alb_coeff
  real :: FAC, GAMMA, BETA, ALPHA
  real :: DX, DY, ALA, WRMFAC
  
  real, dimension(2) :: GRN = (/0.33, 0.67/)
  real, parameter :: TICE = 273.16
  
  real, dimension(4,NTYPS) :: SNWALBold = &
       reshape((/.85, .50, .85, .50,      &
                 .85, .50, .85, .50,      &  
                 .85, .50, .85, .50,      &  
                 .85, .50, .85, .50,      & 
                 .85, .50, .85, .50,      &  
                 .85, .50, .85, .50,      &
                 .85, .50, .85, .50,      &
                 .85, .50, .85, .50,      &
                 .85, .50, .85, .50        /), (/ 4, NTYPS /))
  real, dimension(NTYPS) :: SNWMIDold = (/50.,50.,50.,2.,50.,2.,2.,2.,2./)
  real :: SNALB(4, UMDNTYPS)


!=== End of variable definitions and assignments =========================

!=== Begin calculating albedos ===========================================
!FPP$ EXPAND (COEFF)

  do I=1,IRUN

     TEMPLAI(I) = VLAI(I)

!=== Added check to make sure total LAI is above 0.15 to avoid model run termination     
!     IF (ITYP(I) .NE. 12 .AND. VLAI(I) .LT. 0.50) THEN
      IF (ITYP(I) .NE. 12 .AND. TEMPLAI(I) .LT. 0.50) THEN
         TEMPLAI(I) = 0.50
     ENDIF

!     ALA = AMIN1 (AMAX1 (ZERO, VLAI(I)), ALATRM)
     ALA = AMIN1 (AMAX1 (ZERO, TEMPLAI(I)), ALATRM)
     LAI = 1 + MAX(0, INT((ALA-BLAI)/DLAI) )
     DX = (ALA - (BLAI+(LAI-1)*DLAI)) * (ONE/DLAI)
     DY = (VGRN(I)- GRN(1)) * (ONE/(GRN(2) - GRN(1)))

     ALPHA = ALB_COEFF (sib%ALVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
     BETA  = ALB_COEFF (sib%BTVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
     GAMMA = ALB_COEFF (sib%GMVDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

     AVISDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
     AVISDF(I) = ALPHA-BETA                                &
          + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))

     ALPHA = ALB_COEFF (sib%ALIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
     BETA  = ALB_COEFF (sib%BTIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)
     GAMMA = ALB_COEFF (sib%GMIDR (1, 1, ITYP (I)), NLAI, LAI ,DX, DY)

     ANIRDR(I) = ALPHA - ZTH(I)*BETA / (GAMMA+ZTH(I))
     ANIRDF(I) = ALPHA-BETA                                &
          + 2.*BETA*GAMMA*(1.-GAMMA*ALOG((1.+GAMMA)/GAMMA))


!----------------------------------------------------
!	This section commented out b/c snow albedo already accounted
!	for in main Mosaic code.
!
!	If this part of the code will be used, then the SNWALBold
!	and SNWMIDold need to be mapped to UMD values, and the 
!	section of code below needs to be un-commented out.

!      EXP039: ALLOW REDUCTION IN ALBEDO WHEN SNOW IS CLOSE TO MELTING.
!      (FROM SIB). WE USE 1K (INSTEAD OF .1) BECAUSE RADIATION
!      IS CALLED EVERY 3 HOURS.

!	  IF (SNW (I) .GT. ZERO) THEN
!	   FAC = SNW(I) / (SNW(I) + SNWMID(ITYP(I)))

!           WRMFAC=1.0
!           IF(TC(I) .GT. TICE-1.0) WRMFAC=0.6

!	   AVISDR(I) = AVISDR(I) +
!     &  (SNWALB(1,ITYP(I))*WRMFAC - AVISDR(I)) * FAC
!	   ANIRDR(I) = ANIRDR(I) +
!     &  (SNWALB(2,ITYP(I))*WRMFAC - ANIRDR(I)) * FAC
!	   AVISDF(I) = AVISDF(I) +
!     &  (SNWALB(3,ITYP(I))*WRMFAC - AVISDF(I)) * FAC
!	   ANIRDF(I) = ANIRDF(I) +
!     &  (SNWALB(4,ITYP(I))*WRMFAC - ANIRDF(I)) * FAC
!	  ENDIF

  end do !i up to irun

end subroutine umd_sibalb


!=========================================================================
! FUNKSHUN!!!
!=========================================================================
REAL FUNCTION ALB_COEFF(TABLE, NTABL, LAI ,DX, DY)

  INTEGER NTABL, LAI
  REAL TABLE (NTABL, 2), DX, DY

  ALB_COEFF = (TABLE(LAI,  1)                               &
  + (TABLE(LAI  ,2) - TABLE(LAI  ,1)) * DY ) * (1.0-DX) &
  + (TABLE(LAI+1,1)                                     &
  + (TABLE(LAI+1,2) - TABLE(LAI+1,1)) * DY ) * DX

  RETURN
END FUNCTION ALB_COEFF

!=========================================================================
