!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_FORC_AttributesMod
!BOP
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

  public :: LDT_FORC_Tair
  public :: LDT_FORC_Qair
  public :: LDT_FORC_SWdown
  public :: LDT_FORC_SWdirect
  public :: LDT_FORC_SWdiffuse
  public :: LDT_FORC_LWdown
  public :: LDT_FORC_Wind_E
  public :: LDT_FORC_Wind_N
  public :: LDT_FORC_Psurf
  public :: LDT_FORC_Rainf
  public :: LDT_FORC_Snowf
  public :: LDT_FORC_CRainf
  public :: LDT_FORC_TotalPrecip

  public :: LDT_FORC_Forc_Hgt
  public :: LDT_FORC_Ch
  public :: LDT_FORC_Cm
  public :: LDT_FORC_Emiss
  public :: LDT_FORC_Q2sat
  public :: LDT_FORC_Cosz
  public :: LDT_FORC_Alb
  public :: LDT_FORC_XICE
  public :: LDT_FORC_QSFC

  public :: LDT_FORC_CHS2
  public :: LDT_FORC_CQS2
  public :: LDT_FORC_T2
  public :: LDT_FORC_Q2
  public :: LDT_FORC_TH2
  public :: LDT_FORC_TMN

  public :: LDT_FORC_LPressure
  public :: LDT_FORC_O3     ! absorber
  
  public :: LDT_FORC_PET    ! SY for FEWSNET
  public :: LDT_FORC_RefET  ! SY for FEWSNET

  public :: LDT_FORC_CMFORC ! dmm
  public :: LDT_FORC_CHFORC ! dmm for NLDAS-2
  public :: LDT_FORC_CAPE   ! dmm for NLDAS-2

  public :: LDT_FORC_PARDR
  public :: LDT_FORC_PARDF
  public :: LDT_FORC_SWNET
!<for vic>
  public :: LDT_FORC_SNOWFLAG
  public :: LDT_FORC_DENSITY
  public :: LDT_FORC_VAPORPRESS
  public :: LDT_FORC_VAPORPRESSDEFICIT
  public :: LDT_FORC_WIND
!</for vic>

!ccc - for CABLE
  public :: LDT_FORC_CO2
  
!EOP

  type, public :: LDT_forcDataEntry
     character*100, allocatable :: varname(:)
     character*100              :: standard_name
     character*100              :: short_name
     integer                :: selectOpt
     integer                :: vlevels
!     character*10           :: units
     character*20           :: units
     integer                :: nunits
     character(len=20)      :: dir
     integer                :: ndirs
     real                   :: valid_min
     real                   :: valid_max
     integer                :: selectProc
     integer                :: timeAvgOpt
     integer                :: varId_def
     integer                :: varId_opt
     integer, allocatable   :: count(:,:)
     real, allocatable      :: modelOutput(:,:,:) 
     character(len=20), allocatable :: unittypes(:)
     character(len=20), allocatable :: dirtypes(:)

  end type LDT_forcDataEntry

  type(LDT_forcDataEntry)   :: LDT_FORC_Tair
  type(LDT_forcDataEntry)   :: LDT_FORC_Qair
  type(LDT_forcDataEntry)   :: LDT_FORC_SWdown
  type(LDT_forcDataEntry)   :: LDT_FORC_SWdirect
  type(LDT_forcDataEntry)   :: LDT_FORC_SWdiffuse
  type(LDT_forcDataEntry)   :: LDT_FORC_LWdown
  type(LDT_forcDataEntry)   :: LDT_FORC_Wind_E
  type(LDT_forcDataEntry)   :: LDT_FORC_Wind_N
  type(LDT_forcDataEntry)   :: LDT_FORC_Psurf
  type(LDT_forcDataEntry)   :: LDT_FORC_Rainf
  type(LDT_forcDataEntry)   :: LDT_FORC_Snowf
  type(LDT_forcDataEntry)   :: LDT_FORC_CRainf
  type(LDT_forcDataEntry)   :: LDT_FORC_TotalPrecip

  type(LDT_forcDataEntry)   :: LDT_FORC_Forc_Hgt
  type(LDT_forcDataEntry)   :: LDT_FORC_Ch
  type(LDT_forcDataEntry)   :: LDT_FORC_Cm
  type(LDT_forcDataEntry)   :: LDT_FORC_Emiss
  type(LDT_forcDataEntry)   :: LDT_FORC_Q2sat
  type(LDT_forcDataEntry)   :: LDT_FORC_Cosz
  type(LDT_forcDataEntry)   :: LDT_FORC_Alb
  type(LDT_forcDataEntry)   :: LDT_FORC_XICE

  type(LDT_forcDataEntry)   :: LDT_FORC_QSFC
  type(LDT_forcDataEntry)   :: LDT_FORC_CHS2
  type(LDT_forcDataEntry)   :: LDT_FORC_CQS2
  type(LDT_forcDataEntry)   :: LDT_FORC_T2
  type(LDT_forcDataEntry)   :: LDT_FORC_Q2
  type(LDT_forcDataEntry)   :: LDT_FORC_TH2
  type(LDT_forcDataEntry)   :: LDT_FORC_TMN

  type(LDT_forcDataEntry)   :: LDT_FORC_LPressure
  type(LDT_forcDataEntry)   :: LDT_FORC_O3

  type(LDT_forcDataEntry)   :: LDT_FORC_PET ! SY for FEWSNET
  type(LDT_forcDataEntry)   :: LDT_FORC_RefET ! SY for FEWSNET

  type(LDT_forcDataEntry)   :: LDT_FORC_CMFORC ! dmm
  type(LDT_forcDataEntry)   :: LDT_FORC_CHFORC ! dmm for NLDAS-2
  type(LDT_forcDataEntry)   :: LDT_FORC_CAPE ! dmm for NLDAS-2

  type(LDT_forcDataEntry)   :: LDT_FORC_PARDR
  type(LDT_forcDataEntry)   :: LDT_FORC_PARDF
  type(LDT_forcDataEntry)   :: LDT_FORC_SWNET
!<for vic>
  type(LDT_forcDataEntry)   :: LDT_FORC_SNOWFLAG
  type(LDT_forcDataEntry)   :: LDT_FORC_DENSITY
  type(LDT_forcDataEntry)   :: LDT_FORC_VAPORPRESS
  type(LDT_forcDataEntry)   :: LDT_FORC_VAPORPRESSDEFICIT
  type(LDT_forcDataEntry)   :: LDT_FORC_WIND
!</for vic>

!ccc - for CABLE
  type(LDT_forcDataEntry)   :: LDT_FORC_CO2

end module LDT_FORC_AttributesMod

