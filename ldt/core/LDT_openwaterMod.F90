!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
module LDT_openwaterMod
!BOP
!
! !MODULE: LDT_openwaterMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read openwater fraction
!  and type data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  openwater depth data and generates openwater surface fraction. 
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  15 Apr 2013: Kristi Arsenault: Modified for openwater datasets.
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_SurfaceTypeMod, only: LDT_assign_openwatersfctype
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_openwater_init    !allocates memory for required structures
!  public :: LDT_openwater_writeHeader
!  public :: LDT_openwater_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------

contains

!BOP
! 
! !ROUTINE: LDT_openwater_init
! \label{LDT_openwater_init}
! 
! !INTERFACE:
  subroutine LDT_openwater_init

! !USES:
  use LDT_domainMod, only: isSurfaceTypeSelected

! !DESCRIPTION:
!  Handles assigning open-water/openwater points, if
!  needing to account for these surface types when
!  running LIS-7.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[openwatersetup](\ref{openwatersetup}) \newline
!    calls the registry to invoke the openwater setup methods. 
!  \end{description}
!
!EOP
    implicit none
    integer   :: n
    integer   :: k, c ,r
    integer   :: rc
    real, allocatable :: openwater_fgrd(:,:,:)
! ____________________________________________

  do k = 1, LDT_rc%nsf_model_types
    if( LDT_rc%sf_model_type_name_select(k) == "Lake" ) then
   ! Since the openwater surface type is accounted for when lake
   !  surface is selected, including openwater type here is ignored ...
     return
    endif
  enddo

  do k = 1, LDT_rc%nsf_model_types
    if( LDT_rc%sf_model_type_name_select(k) == "Openwater" ) then

   !- Assign Open water surface types:
      do n = 1,LDT_rc%nnest
         allocate(openwater_fgrd( LDT_rc%lnc(n),LDT_rc%lnr(n),1)) !&
!                            LDT_LSMparam_struc(n)%openwater%num_bins) )

      !- Update surface type information:
         call LDT_assign_openwatersfctype( n, &
                  LDT_LSMparam_struc(n)%sfctype%num_bins,   &
                  LDT_LSMparam_struc(n)%sfctype%value,      &
                  LDT_LSMparam_struc(n)%dommask%value,     &
                  LDT_LSMparam_struc(n)%landcover%value,    &
                  openwater_fgrd )
!                  LDT_LSMparam_struc(n)%openwaterfrac%value    )

         write(LDT_logunit,*) "Finished assigning openwater surface types"
      enddo

    endif
  enddo

  end subroutine LDT_openwater_init


#if 0
  subroutine LDT_openwater_writeHeader(n,ftn,dimID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer   :: n
    integer   :: ftn
    integer   :: dimID(3)
    integer   :: tdimID(3)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
             LDT_LSMparam_struc(n)%inlandwatertype)

#endif
  end subroutine LDT_openwater_writeHeader

  subroutine LDT_openwater_writeData(n,ftn)

    integer     :: n 
    integer     :: ftn

    if( LDT_LSMparam_struc(n)%inlandwatertype%selectOpt.gt.0 ) &
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%inlandwatertype)

  end subroutine LDT_openwater_writeData

#endif

end module LDT_openwaterMod
