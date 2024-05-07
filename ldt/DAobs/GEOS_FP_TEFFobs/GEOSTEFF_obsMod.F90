!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !MODULE: GEOSTEFF_obsMod
! 
! !DESCRIPTION: 
! This module handles the observation plugin for the 
! GEOS FP soil temperature
!
! !REVISION HISTORY: 
!  29 Nov 2021: Yonghwan Kwon, Initial Specification
!
module GEOSTEFF_obsMod
! !USES:
  use ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: GEOSTeffobsinit !Initializes structures for reading GEOS soil Temp data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: GEOSTeffobs !Object to hold GEOSTeff observation attributes
!EOP
  type, public :: geoTeffdec

     character*100        :: odir
     integer              :: nc, nr
     real                 :: gridDesci(50)
     real,    allocatable :: teffobs(:,:)
     integer, allocatable :: n11(:)
     integer, allocatable :: n12(:)
     integer, allocatable :: n21(:)
     integer, allocatable :: n22(:)
     real,    allocatable :: w11(:)
     real,    allocatable :: w12(:)
     real,    allocatable :: w21(:)
     real,    allocatable :: w22(:)

  end type geoTeffdec

  type(geoTeffdec),allocatable :: GEOSTeffobs(:)

contains

!BOP
! 
! !ROUTINE: GEOSTeffobsinit
! \label{GEOSTeffobsinit}
!
! !INTERFACE: 
  subroutine GEOSTeffobsinit()
! 
! !USES:
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_timeMgrMod
    use LDT_logMod

    implicit none
! !ARGUMENTS:

! 
! !DESCRIPTION: 
!  This subroutine initializes and sets up the data structures required
!  for reading the GEOS Teff data.
! 
!EOP

    integer                 :: status, rc
    integer                 :: n

    allocate(GEOSTeffobs(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config, &
         'GEOS soil temperature data directory:', rc=status)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_Config, GEOSTeffobs(n)%odir, &
            rc=status)
       call LDT_verify(status, &
            'GEOS soil temperature data directory: not defined')
    enddo

    do n=1,LDT_rc%nnest

       allocate(GEOSTeffobs(n)%teffobs(LDT_rc%lnc(n),LDT_rc%lnr(n)))

       GEOSTeffobs(n)%teffobs = -9999.0

       call LDT_initializeDAobsEntry(LDT_DAobsData(n)%teff_obs, "K",1,1)
       LDT_DAobsData(1)%teff_obs%selectStats = 1

       GEOSTeffobs(n)%nc = 1152
       GEOSTeffobs(n)%nr = 721

       GEOSTeffobs(n)%gridDesci(1) = 0
       GEOSTeffobs(n)%gridDesci(2) = GEOSTeffobs(n)%nc
       GEOSTeffobs(n)%gridDesci(3) = GEOSTeffobs(n)%nr
       GEOSTeffobs(n)%gridDesci(4) = -90.0
       GEOSTeffobs(n)%gridDesci(5) = -180.0
       GEOSTeffobs(n)%gridDesci(6) = 128
       GEOSTeffobs(n)%gridDesci(7) = 90.0
       GEOSTeffobs(n)%gridDesci(8) = 179.6875
       GEOSTeffobs(n)%gridDesci(9) = 0.3125   !dlon
       GEOSTeffobs(n)%gridDesci(10) = 0.25    !dlat
       GEOSTeffobs(n)%gridDesci(20) = 64

       if(LDT_isLDTatAfinerResolution(n,0.3125)) then

          allocate(GEOSTeffobs(n)%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GEOSTeffobs(n)%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GEOSTeffobs(n)%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GEOSTeffobs(n)%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GEOSTeffobs(n)%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GEOSTeffobs(n)%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GEOSTeffobs(n)%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(GEOSTeffobs(n)%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))

          call bilinear_interp_input (n, &
               GEOSTeffobs(n)%gridDesci,&
               GEOSTeffobs(n)%n11,&
               GEOSTeffobs(n)%n12,&
               GEOSTeffobs(n)%n21,&
               GEOSTeffobs(n)%n22,&
               GEOSTeffobs(n)%w11,&
               GEOSTeffobs(n)%w12,&
               GEOSTeffobs(n)%w21,&
               GEOSTeffobs(n)%w22)

       else

          allocate(GEOSTeffobs(n)%n11(GEOSTeffobs(n)%nc*&
               GEOSTeffobs(n)%nr))

          call upscaleByAveraging_input (&
               GEOSTeffobs(n)%gridDesci,&
               LDT_rc%gridDesc(n,:),&
               GEOSTeffobs(n)%nc*&
               GEOSTeffobs(n)%nr,&
               LDT_rc%lnc(n)*LDT_rc%lnr(n),&
               GEOSTeffobs(n)%n11)

       endif
    enddo
  end subroutine GEOSTeffobsinit

end module GEOSTEFF_obsMod
