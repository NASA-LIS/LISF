!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_FORC_AttributesMod
!BOP
!
!
!
!  !REVISION HISTORY: 
!  14 Nov 2002    Sujay Kumar  Initial Specification
!  26 May 2011    Soni Yatheendradas    Potential ET variables added For FEWSNET 
!  14 Mar 2014    David Mocko: Added CAPE, CH, and CM to forcing variables
! 
  use ESMF

  implicit none 

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------

  public :: LIS_FORC_Tair
  public :: LIS_FORC_Qair
  public :: LIS_FORC_SWdown
  public :: LIS_FORC_SWdirect
  public :: LIS_FORC_SWdiffuse
  public :: LIS_FORC_LWdown
  public :: LIS_FORC_Wind_E
  public :: LIS_FORC_Wind_N
  public :: LIS_FORC_Psurf
  public :: LIS_FORC_Rainf
  public :: LIS_FORC_Snowf
  public :: LIS_FORC_CRainf

  public :: LIS_FORC_Forc_Hgt
  public :: LIS_FORC_Ch
  public :: LIS_FORC_Cm
  public :: LIS_FORC_Emiss
  public :: LIS_FORC_Q2sat
  public :: LIS_FORC_Cosz
  public :: LIS_FORC_Alb
  public :: LIS_FORC_XICE
  public :: LIS_FORC_QSFC

  public :: LIS_FORC_CHS2
  public :: LIS_FORC_CQS2
  public :: LIS_FORC_T2
  public :: LIS_FORC_Q2
  public :: LIS_FORC_TH2
  public :: LIS_FORC_TMN

  public :: LIS_FORC_LPressure
  public :: LIS_FORC_O3    !absorber
  
  public :: LIS_FORC_PET ! SY for FEWSNET
  public :: LIS_FORC_RefET ! SY for FEWSNET

  public :: LIS_FORC_CMFORC ! dmm
  public :: LIS_FORC_CHFORC ! dmm for NLDAS-2
  public :: LIS_FORC_CAPE ! dmm for NLDAS-2

  public :: LIS_FORC_PARDR
  public :: LIS_FORC_PARDF
  public :: LIS_FORC_SWNET
!<for vic>
  public :: LIS_FORC_SNOWFLAG
  public :: LIS_FORC_DENSITY
  public :: LIS_FORC_VAPORPRESS
  public :: LIS_FORC_VAPORPRESSDEFICIT
  public :: LIS_FORC_WIND
!</for vic>

  public :: LIS_FORC_Z0
  public :: LIS_FORC_GVF
!ccc - for CABLE
  public :: LIS_FORC_CO2
  
!EOP

  type, public :: forc_attrib_type
     character*100, allocatable :: varname(:)
     integer                :: selectOpt
     integer                :: vlevels
     character*10           :: units
  end type forc_attrib_type

  type(forc_attrib_type)   :: LIS_FORC_Tair
  type(forc_attrib_type)   :: LIS_FORC_Qair
  type(forc_attrib_type)   :: LIS_FORC_SWdown
  type(forc_attrib_type)   :: LIS_FORC_SWdirect
  type(forc_attrib_type)   :: LIS_FORC_SWdiffuse
  type(forc_attrib_type)   :: LIS_FORC_LWdown
  type(forc_attrib_type)   :: LIS_FORC_Wind_E
  type(forc_attrib_type)   :: LIS_FORC_Wind_N
  type(forc_attrib_type)   :: LIS_FORC_Psurf
  type(forc_attrib_type)   :: LIS_FORC_Rainf
  type(forc_attrib_type)   :: LIS_FORC_Snowf
  type(forc_attrib_type)   :: LIS_FORC_CRainf

  type(forc_attrib_type)   :: LIS_FORC_Forc_Hgt
  type(forc_attrib_type)   :: LIS_FORC_Ch
  type(forc_attrib_type)   :: LIS_FORC_Cm
  type(forc_attrib_type)   :: LIS_FORC_Emiss
  type(forc_attrib_type)   :: LIS_FORC_Q2sat
  type(forc_attrib_type)   :: LIS_FORC_Cosz
  type(forc_attrib_type)   :: LIS_FORC_Alb
  type(forc_attrib_type)   :: LIS_FORC_XICE

  type(forc_attrib_type)   :: LIS_FORC_QSFC
  type(forc_attrib_type)   :: LIS_FORC_CHS2
  type(forc_attrib_type)   :: LIS_FORC_CQS2
  type(forc_attrib_type)   :: LIS_FORC_T2
  type(forc_attrib_type)   :: LIS_FORC_Q2
  type(forc_attrib_type)   :: LIS_FORC_TH2
  type(forc_attrib_type)   :: LIS_FORC_TMN

  type(forc_attrib_type)   :: LIS_FORC_LPressure
  type(forc_attrib_type)   :: LIS_FORC_O3

  type(forc_attrib_type)   :: LIS_FORC_PET ! SY for FEWSNET
  type(forc_attrib_type)   :: LIS_FORC_RefET ! SY for FEWSNET

  type(forc_attrib_type)   :: LIS_FORC_CMFORC ! dmm
  type(forc_attrib_type)   :: LIS_FORC_CHFORC ! dmm for NLDAS-2
  type(forc_attrib_type)   :: LIS_FORC_CAPE ! dmm for NLDAS-2

  type(forc_attrib_type)   :: LIS_FORC_PARDR
  type(forc_attrib_type)   :: LIS_FORC_PARDF
  type(forc_attrib_type)   :: LIS_FORC_SWNET
!<for vic>
  type(forc_attrib_type)   :: LIS_FORC_SNOWFLAG
  type(forc_attrib_type)   :: LIS_FORC_DENSITY
  type(forc_attrib_type)   :: LIS_FORC_VAPORPRESS
  type(forc_attrib_type)   :: LIS_FORC_VAPORPRESSDEFICIT
  type(forc_attrib_type)   :: LIS_FORC_WIND
!</for vic>

  type(forc_attrib_type)   :: LIS_FORC_GVF
  type(forc_attrib_type)   :: LIS_FORC_Z0
!ccc - for CABLE
  type(forc_attrib_type)   :: LIS_FORC_CO2

end module LIS_FORC_AttributesMod

