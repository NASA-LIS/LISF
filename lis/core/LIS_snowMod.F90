!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! Macros for tracing - Requires ESMF 7_1_0+
#ifdef ESMF_TRACE
#define TRACE_ENTER(region) call ESMF_TraceRegionEnter(region)
#define TRACE_EXIT(region) call ESMF_TraceRegionExit(region)
#else
#define TRACE_ENTER(region)
#define TRACE_EXIT(region)
#endif
module LIS_snowMod
!BOP
!
! !MODULE: LIS_snowMod
!
! !DESCRIPTION:
!  The code in this file implementts routines to read various sources
!  of snow depth data. Currently this module is used to ingest the 
!  SNODEP products from AFWA, which is obtained at 12z daily.  
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!
! !USES: 
  implicit none
  
  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_snow_setup       !allocates memory for required variables
  public :: LIS_snow_finalize    !cleanup allocated structures

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_snow_struc      !data structure containing snow variables
!EOP  
 
  type, public :: snow_type_dec
     real, allocatable :: snowdepth(:)
     real, allocatable :: sneqv(:)
  end type snow_type_dec

  type(snow_type_dec), allocatable :: LIS_snow_struc(:)

contains

!BOP
! 
! !ROUTINE: LIS_snow_setup
! \label{LIS_snow_setup}
! 
! !DESCRIPTION:
!
! Allocates memory for data structures used for reading 
! snow datasets. The snowdepth field is updated by the external
! files. The snow water equivalent fields are expected to be set
! by the model. 
! 
! !INTERFACE:
  subroutine LIS_snow_setup
! !USES:
    use ESMF
    use LIS_coreMod,    only : LIS_rc
!EOP
    integer :: n
    
    TRACE_ENTER("snow_setup")
    allocate(LIS_snow_struc(LIS_rc%nnest))

    do n=1,LIS_rc%nnest          
!       allocate(LIS_snow_struc(n)%snowdepth(LIS_rc%ntiles(n)))
       allocate(LIS_snow_struc(n)%snowdepth(LIS_rc%ngrid(n)))
       allocate(LIS_snow_struc(n)%sneqv(LIS_rc%ntiles(n)))
       LIS_snow_struc(n)%sneqv = 0.0
       LIS_snow_struc(n)%snowdepth = 0.0
       LIS_rc%snowsrc(n) = 1
    enddo
    TRACE_EXIT("snow_setup")
  end subroutine LIS_snow_setup

!BOP
! 
! !ROUTINE: LIS_snow_finalize
! \label{LIS_snow_finalize}
! 
! !DESCRIPTION:
!
! This routine cleans up snow-related structures 
! 
! !INTERFACE:
  subroutine LIS_snow_finalize
! !USES: 
    use LIS_coreMod, only : LIS_rc
!EOP
    implicit none

    deallocate(LIS_snow_struc)
  end subroutine LIS_snow_finalize

end module LIS_snowMod
