!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module sibalb_module 
!BOP
! !MODULE: sibalb_module
!
! !DESCRIPTION:
!  In order to use the MOSAIC-SiB albedo calculation, SiB vegetation types
!  are mapped to UMD vegetation types.
!  For both SiB and UMD, 4 coefficient arrays are defined for calculating 
!  the albedos of soil in the subroutine sibalb, representing the 
!  following bands: \newline
!  visible, direct solar radiation, \newline
!  infra-red, direct solar radiation, \newline
!  visible, diffuse solar radiation, and \newline
!  infra-red, diffuse solar radiation. \newline
!
! MOSAIC/SiB ITYP: Vegetation type as follows: \newline
!                  1:  BROADLEAF EVERGREEN TREES \newline
!                  2:  BROADLEAF DECIDUOUS TREES \newline
!                  3:  NEEDLELEAF TREES \newline
!                  4:  GROUND COVER \newline
!                  5:  BROADLEAF SHRUBS \newline
!                  6:  DWARF TREES (TUNDRA) \newline
!                  7:  BARE SOIL \newline
!
! UMD ITYP: Vegetation type as follows: \newline
!                  1.  Evergreen Needleleaf Forest \newline
!                  2.  Evergreen Broadleaf Forest \newline
!                  3.  Deciduous Needleleaf Forest \newline
!                  4.  Deciduous Broadleaf Forest \newline
!                  5.  Mixed Cover \newline
!                  6.  Woodland \newline
!                  7.  Wooded Grassland \newline
!                  8.  Closed Shrubland \newline
!                  9.  Open Shrubland \newline
!                 10.  Grassland \newline
!                 11.  Cropland \newline
!                 12.  Bare Ground \newline
!                 13.  Urban and Built-Up \newline
!
! !REVISION HISTORY:
!  04 Apr 2001: Urszula Jambor; Initial code, using old calc_albedo.f scheme
!EOP
  IMPLICIT NONE
  public sibalbdec
  
  type sibalbdec 
         
!=== Constants used in albedo calculations: =========

!=== Original SiB albedo coefficients and SiB mapped to UMD coefficients
!=== dimension (NLAI=14, 2, NTYPS=9)
!=== dimension (NLAI=14, 2, UMDNTYPS=13)

     real, dimension(14, 2, 9) :: ALVDRold, BTVDRold, GMVDRold
     real, dimension(14, 2, 9) :: ALIDRold, BTIDRold, GMIDRold
     
     real, dimension(14, 2, 13) :: ALVDR, BTVDR, GMVDR
     real, dimension(14, 2, 13) :: ALIDR, BTIDR, GMIDR
     
  end type sibalbdec !sibalbdec
  
  type (sibalbdec)   sib
  
end module sibalb_module


