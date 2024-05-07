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
module VIC_parmsMod
!BOP
!
! !MODULE: VIC_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read VIC parameter
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  VIC parameter file data.
!
! !REVISION HISTORY:
!
!  04 Nov 2013: K. Arsenault: Added layers for VIC model
!  18 Oct 2017: H. Beaudoing: Added variables read in config file
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
  public :: VICparms_init    !allocates memory for required structures
  public :: VICparms_writeHeader
  public :: VICparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: VIC_struc

  type, public :: vic_type_dec

     real      :: vic_undef

! -  VIC-specific:
     integer :: Nlayer             !soil layer default=3
     integer :: Nnode              !default=12
     integer :: Ndist              !default=12
     integer :: Nbands             !snow band default=25
     integer :: Nveg               !always 1??

     type(LDT_paramEntry) :: vic  ! VIC parameters (collective)

     real, allocatable  :: state_chunk(:,:)  ! 05/30/2014

  end type vic_type_dec

  type(vic_type_dec), allocatable :: VIC_struc(:)

contains

!BOP
! 
! !ROUTINE: VICparms_init
! \label{VICparms_init}
! 
! !INTERFACE:
  subroutine VICparms_init(flag)

! !USES:
   use LDT_logMod,    only : LDT_verify, LDT_endrun, &
             LDT_getNextUnitNumber, LDT_releaseUnitNumber
!
! !DESCRIPTION:
!
! Variable Infiltration Capacity (VIC) model parameters.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[vicParmssetup](\ref{vicParmssetup}) \newline
!    calls the registry to invoke the vicParms setup methods. 
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
   logical  :: vic_select 

   integer :: veg_tiling_scheme  ! 0=VIC, 1=LIS
   character(len=8) :: DIST_PRCP  !FALSE or TRUE

   ! _____________________________________________________________________

   allocate(VIC_struc(LDT_rc%nnest))

   vic_select = .false.
   do n=1,LDT_rc%nnest
      ! - VIC parameters:
      call set_param_attribs(VIC_struc(n)%vic,"VIC")
      if( VIC_struc(n)%vic%selectOpt == 1 ) then
         vic_select = .true.
      endif
   enddo

!hkb read in VIC layers from config
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,veg_tiling_scheme,&
           label="VIC412 veg tiling scheme:",default=1,rc=rc)
       if ( veg_tiling_scheme .eq. 1 .or. veg_tiling_scheme .eq. 0 ) then 
             VIC_struc(n)%Nveg = 1
       endif
      call ESMF_ConfigGetAttribute(LDT_config,VIC_struc(n)%Nlayer,&
           label="VIC412_NLAYER:",default=3,rc=rc)
      call ESMF_ConfigGetAttribute(LDT_config,VIC_struc(n)%Nnode, &
           label="VIC412_NODES:",default=12,rc=rc)
      call ESMF_ConfigGetAttribute(LDT_config,DIST_PRCP, &
           label="VIC412_DIST_PRCP:",default="FALSE",rc=rc)
       if ( DIST_PRCP .eq. "FALSE" ) then
            VIC_struc(n)%Ndist = 1
       else
            VIC_struc(n)%Ndist = 2
       endif
      call ESMF_ConfigGetAttribute(LDT_config,VIC_struc(n)%Nbands, &
           label="VIC412_SNOW_BAND:",default=25,rc=rc)
   enddo

   if( vic_select ) then
     write(LDT_logunit,*)" - - - - - - - - - - VIC LSM Parameters - - - - - - - - - - - - -"

     write(LDT_logunit,*)"[INFO] VIC Model Parameters are available for"
     write(LDT_logunit,*)"       specific landcover and landmask,"
     write(LDT_logunit,*)"       veg_tiling_scheme, NLAYER, NODES,"
     write(LDT_logunit,*)"       DIST_PRCP, and SNOW_BAND"

   endif


 end subroutine VICparms_init


 subroutine VICparms_writeHeader(n,ftn,dimID,monthID)
   
   integer   :: n 
   integer   :: ftn
   integer   :: dimID(3)
   integer   :: monthID
   
   integer   :: t_dimID(3)
   
   if( VIC_struc(n)%vic%selectOpt == 1 ) then

   endif
   
   
 end subroutine VICparms_writeHeader
 
  subroutine VICparms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    if( VIC_struc(n)%vic%selectOpt == 1 ) then
!      call VICparms_writeData(n,ftn)
!      call VICparms_writeData(n,ftn)
    endif

  end subroutine VICparms_writeData

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
   paramEntry%source = "VIC"
   paramEntry%units ="none"
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(short_name)

  end subroutine set_param_attribs

end module VIC_parmsMod

