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
module Mosaic_parmsMod
!BOP
!
! !MODULE: Mosaic_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read Mosaic parameter
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  Mosaic parameter file data.
!
! !REVISION HISTORY:
!
!  04 Nov 2014: K. Arsenault: Added layers for Mosaic model
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: Mosaicparms_init    !allocates memory for required structures
  public :: Mosaicparms_writeHeader
  public :: Mosaicparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: Mosaic_struc

  type, public :: mos_type_dec

     real      :: mos_undef

! -  Mosaic-specific:
     type(LDT_paramEntry) :: mos  ! Mosaic parameters (collective)

  end type mos_type_dec

  type(mos_type_dec), allocatable :: Mosaic_struc(:)

contains

!BOP
! 
! !ROUTINE: Mosaicparms_init
! \label{Mosaicparms_init}
! 
! !INTERFACE:
  subroutine Mosaicparms_init(flag)

! !USES:
   use LDT_logMod,  only : LDT_verify, LDT_endrun, &
             LDT_getNextUnitNumber, LDT_releaseUnitNumber
!
! !DESCRIPTION:
!
! NASA's Mosaic land surface model parameters.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[mosParmssetup](\ref{mosParmssetup}) \newline
!    calls the registry to invoke the mosParms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer  :: flag
   integer  :: n
   integer  :: c,r,m,k
   integer  :: rc
   integer  :: file_status
   logical  :: file_exists
   logical  :: mos_select 

 ! _____________________________________________________________________

   allocate(Mosaic_struc(LDT_rc%nnest))

   mos_select = .false.
   do n=1,LDT_rc%nnest
      ! - Mosaic parameters:
      call set_param_attribs(Mosaic_struc(n)%mos,"Mosaic")
      if( Mosaic_struc(n)%mos%selectOpt == 1 ) then
         mos_select = .true.
      endif
   enddo

   if( mos_select ) then
     write(LDT_logunit,*)" - - - - - - - - - - Mosaic LSM Parameters - - - - - - - - - - - - -"

     write(LDT_logunit,*)"[INFO] NO Mosaic Model Parameters available yet, other than"
     write(LDT_logunit,*)"[INFO]  porosity and LAI."
      
   endif


 end subroutine Mosaicparms_init


 subroutine Mosaicparms_writeHeader(n,ftn,dimID,monthID)
   
   integer   :: n 
   integer   :: ftn
   integer   :: dimID(3)
   integer   :: monthID
   
   integer   :: t_dimID(3)
   
   if( Mosaic_struc(n)%mos%selectOpt == 1 ) then

   endif
   
   
 end subroutine Mosaicparms_writeHeader
 
  subroutine Mosaicparms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    if( Mosaic_struc(n)%mos%selectOpt == 1 ) then
!      call Mosaicparms_writeData(n,ftn)
!      call Mosaicparms_writeData(n,ftn)
    endif

  end subroutine Mosaicparms_writeData

!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name)

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name

! ____________________________________________________
    
   
   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = 1
   paramEntry%selectOpt = 1
   paramEntry%source = "Mosaic"
   paramEntry%units ="none"
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(short_name)

  end subroutine set_param_attribs

end module Mosaic_parmsMod

