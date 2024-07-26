!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: ESACCIsm_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! soil moisture Essential Climate Variable (ESACCI) product. 
!
! Soil moisture is recognized as an Essential Climate Variable (ECV) by NASA, 
! European Space Agency (ESA), and other agencies and institutions.  

! The Water Cycle Multi-mission Observation Strategy (WACMOS) and Climate Change 
! Initiative (CCI) Soil Moisture projects have released a 32-year harmonized
! product, from 1978-2010, and referred to as the ESACCI soil moisture dataset. 
! The product includes both active and passive microwave sensor retrievals.  The 
! active dataset provided by the University of Vienna (TU Wien) uses ERS-1, ERS-2 and
! METOP-A C-band scatterometers.  The passive microwave dataset is developed by
! VU University of Amsterdam and NASA, and includes NASA-based measurements from 
! Nimbus 7 SSMR, DMSP SSM/I, TRMM TMI and Aqua AMSR-E sensors.
!
! The merged dataset uses a fixed ranking system based on error derivations
! from triple collocation.  So if AMSR has an higher error than ASCAT, then ASCAT is used
! (and visa versa, see for example Liu et al.,RSE 2012 and Liu et al., HESS 2011).
! In the near future, a dynamic weighting function is planned to be used within the
! merging routine (considered research in progress).
!
! Further information about the product can be obtained at:
! http://www.esa-soilmoisture-cci.org/
!
! REFERENCES: 
! Liu, Y. Y., W. A. Dorigo, et al. (2012). "Trend-preserving blending of passive and 
! active microwave soil moisture retrievals." Remote Sensing of Environment 123: 280-297.
!
! Wagner, W., Dorigo, W., de Jeu, R., Fernandez, D., Benveniste, J., Haas, E., 
! and Ertl, M.: Fusion of active and passive microwave observations to create 
! an Essential Climate Variable data record on soil moisture. ISPRS Ann. 
! Photogramm. Remote Sens. Spatial Inf. Sci., I-7, 315-321, 
! doi:10.5194/isprsannals-I-7-315-2012, 2012.
!
!   
! !REVISION HISTORY: 
!  01 Oct 2012: Sujay Kumar, Initial Specification
!
module ESACCIsm_obsMod
! !USES: 
  use ESMF
  use map_utils
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE 
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ESACCIsm_obsinit
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: ESACCIsmobs
!EOP
  type, public :: esaccismobsdec

     character(len=LDT_CONST_PATH_LEN)          :: odir
     character(len=8)       :: sensor
     integer                :: mo
     real                   :: version
     real,    allocatable   :: smobs(:,:)
     integer                :: esaccinc, esaccinr
     type(proj_info)        :: esacciproj
     integer, allocatable   :: n11(:)
  end type esaccismobsdec

  type(esaccismobsdec), allocatable:: ESACCIsmobs(:)

contains
  
!BOP
! 
! !ROUTINE: ESACCIsm_obsInit
! \label{ESACCIsm_obsInit}
! 
! !INTERFACE: 
  subroutine ESACCIsm_obsinit()
! !USES: 
    use LDT_coreMod,    only : LDT_rc, LDT_config
    use LDT_DAobsDataMod, only : LDT_DAobsData, LDT_initializeDAobsEntry
    use LDT_timeMgrMod, only : LDT_clock, LDT_calendar
    use LDT_logMod,     only : LDT_verify, LDT_logunit

    implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading ESACCI soil moisture data. 
! 
!EOP
    integer            :: npts
    type(ESMF_TimeInterval) :: alarmInterval
    type(ESMF_Time)         :: alarmTime
    integer                 :: status, rc
    real                    :: gridDesci(20)
    integer                 :: n 

    allocate(ESACCIsmobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'ESA CCI soil moisture observation directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, ESACCIsmobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'ESA CCI soil moisture observation directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'ESA CCI soil moisture version of data:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, ESACCIsmobs(n)%version, &
            rc=status)
       call LDT_verify(status, &
            'ESA CCI soil moisture version of data: not defined')
    enddo

    call ESMF_ConfigFindLabel(LDT_config, &
         'ESA CCI soil moisture sensor type:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, ESACCIsmobs(n)%sensor, &
            rc=status)
       call LDT_verify(status, &
            'ESA CCI soil moisture sensor type: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(ESACCIsmobs(n)%smobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       ESACCIsmobs(n)%smobs = -9999.0
       
       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%soilmoist_obs, &
            "m3/m3",1,1)
       LDT_DAobsData(n)%soilmoist_obs%selectStats = 1
    
       ESACCIsmobs(n)%esaccinc = 1440
       ESACCIsmobs(n)%esaccinr = 720

       call map_set(PROJ_LATLON, -89.875,-179.875,&
            0.0, 0.25,0.25, 0.0,&
            ESACCIsmobs(n)%esaccinc,ESACCIsmobs(n)%esaccinr,&
            ESACCIsmobs(n)%esacciproj)
       
       gridDesci = 0 
       gridDesci(1) = 0 
       gridDesci(2) = 1440
       gridDesci(3) = 720
       gridDesci(4) = -89.875
       gridDesci(5) = -179.875
       gridDesci(6) = 128
       gridDesci(7) = 89.875
       gridDesci(8) = 179.875
       gridDesci(9) = 0.25
       gridDesci(10) = 0.25
       gridDesci(20) = 64
           
       allocate(ESACCIsmobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))       
       call neighbor_interp_input (n, gridDesci,&
            ESACCIsmobs(n)%n11)
    enddo
  end subroutine ESACCIsm_obsinit
     
end module ESACCIsm_obsMod
