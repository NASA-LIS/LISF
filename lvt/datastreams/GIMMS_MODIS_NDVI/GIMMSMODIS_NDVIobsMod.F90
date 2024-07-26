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
! 
! !MODULE: GIMMSMODIS_NDVIobsMod
! \label(GIMMSMODIS_NDVIobsMod)
!
! !INTERFACE:
module GIMMSMODIS_NDVIobsMod
! 
! !USES:   
  use ESMF

  implicit none

  PRIVATE 
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! !REVISION HISTORY: 
!  30 June 2016: Sujay Kumar & Kristi Arsenault, Initial Specification
! 
!EOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GIMMSMODIS_NDVIobsinit !Initializes structures for reading GIMMS MODIS NDVI data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GIMMSMODISNDVIobs !Object to hold GIMMS MODIS NDVI observation attributes
!EOP

  type, public :: gimmsmodisndvidec
     character*200           :: odir
     integer                 :: nc, nr
     real,    allocatable    :: rlat(:)
     real,    allocatable    :: rlon(:)
     integer, allocatable    :: n11(:)
     real                    :: gridDesc(50)
     logical                 :: startFlag
     real                    :: datares
  end type gimmsmodisndvidec
     
  type(gimmsmodisndvidec), save :: GIMMSMODISNDVIObs(2)

contains
  
!BOP
! 
! !ROUTINE: GIMMSMODIS_NDVIobsInit
! \label{GIMMSMODIS_NDVIobsInit}
!
! !INTERFACE: 
  subroutine GIMMSMODIS_NDVIobsinit(i)
! 
! !USES: 
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_logMod
    use LVT_timeMgrMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer,   intent(IN) :: i 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine initializes and sets up the data structures required
!   for reading the GIMMS MODIS NDVI data, including the computation of spatial 
!   interpolation weights. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer               :: status

    call ESMF_ConfigGetAttribute(LVT_Config, GIMMSMODISNDVIobs(i)%odir, &
         label='GIMMS MODIS NDVI data directory:',rc=status)
    call LVT_verify(status, 'GIMMS MODIS NDVI data directory: not defined')

    call LVT_update_timestep(LVT_rc, 86400)
    
    GIMMSMODISNDVIobs(i)%startFlag = .true. 

  end subroutine GIMMSMODIS_NDVIobsinit


end module GIMMSMODIS_NDVIobsMod
