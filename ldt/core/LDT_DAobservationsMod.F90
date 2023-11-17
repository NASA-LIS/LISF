!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LDT_DAobservationsMod
!BOP
!
! !MODULE: LDT_DAobservationsMod
! 
! !DESCRIPTION: 
!  The code in this file contains the basic datastructures and 
!  control routines for handling observational data
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
!  2 Dec 2021:   Mahdi Navari; modified to compute CDF for precipitation
! 
! !USES:       

  use LDT_DAobsDataMod

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_DAobsDataInit
  public :: LDT_readDAobsData
  
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
contains

!BOP
! 
! !ROUTINE: LDT_DAobsDataInit
! \label{LDT_DAobsDataInit}
! 
! !INTERFACE:   
  subroutine LDT_DAobsDataInit
! !USES: 
    use LDT_coreMod, only : LDT_rc
    use LDT_DAobsDataMod, only : LDT_DAobsEntryInit, LDT_DAobsData, &
         LDT_DAobsDataPtr
! 
! !DESCRIPTION: 
! 
!EOP
    implicit none

    integer :: i

    
    allocate(LDT_DAobsData(LDT_rc%nobs))
    allocate(LDT_DAobsDataPtr(LDT_rc%nobs,LDT_DA_MOC_COUNT))
! Initialize all options to zero by default and then let each observation
! plugin override them. 
    do i=1,LDT_rc%nobs
       call default_init_obsEntry(LDT_DAobsData(i)%swe_obs, "SWE")
       call default_init_obsEntry(LDT_DAobsData(i)%snowdepth_obs, "SnowDepth")
       call default_init_obsEntry(LDT_DAobsData(i)%soilmoist_obs, "SoilMoist")
       call default_init_obsEntry(LDT_DAobsData(i)%teff_obs, "SoilTeff")   !Y.Kwon
       call default_init_obsEntry(LDT_DAobsData(i)%tws_obs, "TWS")
       call default_init_obsEntry(LDT_DAobsData(i)%vod_obs, "VOD")
       call default_init_obsEntry(LDT_DAobsData(i)%lai_obs, "LAI")
       call default_init_obsEntry(LDT_DAobsData(i)%gvf_obs, "GVF")   !Y.Kwon
       call default_init_obsEntry(LDT_DAobsData(i)%totalprecip_obs, "TotalPrecip")
    enddo

    call daobservationsetup(trim(LDT_rc%obs_src)//char(0))
! it is assumed that 'observations' are always going to be in a grid space. 
! and processing is done only for one nest. 
    call LDT_DAobsEntryInit(1,LDT_rc%ngrid(1))

  end subroutine LDT_DAobsDataInit

  subroutine default_init_obsEntry(obsEntry, name)

    implicit none

    type(LDT_DAmetadataEntry) :: obsEntry
    character(len=*)        :: name
    
    obsEntry%selectOpt = 0 
    obsEntry%timeAvgOpt = 0 
    obsEntry%vlevels = 0 
    obsEntry%standard_name = name

  end subroutine default_init_obsEntry

!BOP
! 
! !ROUTINE: LDT_readDAobsData
! \label{LDT_readDAobsData}
! 
! !INTERFACE: 
  subroutine LDT_readDAobsData(n)
! !USES: 
    use LDT_coreMod,   only : LDT_rc

    implicit none
    
    integer              :: n 
! 
! !DESCRIPTION: 
! 
!  This subroutine reads the observations and processes them. 
! 
!EOP

    call readDAobservationSource(trim(LDT_rc%obs_src)//char(0),n)

  end subroutine LDT_readDAobsData

end module LDT_DAobservationsMod
