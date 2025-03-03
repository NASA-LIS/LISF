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
module CLM2_parmsMod
!BOP
!
! !MODULE: CLM2_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read CLM2 parameter
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  CLM2 parameter file data.
!
! !REVISION HISTORY:
!
!  04 Nov 2014: K. Arsenault: Added layers for CLM2 model
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
  public :: CLM2parms_init    !allocates memory for required structures
  public :: CLM2parms_writeHeader
  public :: CLM2parms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: CLM2_struc

  type, public :: clm2_type_dec

     real      :: clm2_undef

! -  CLM2-specific:
     type(LDT_paramEntry) :: clm2  ! CLM2 parameters (collective)

  end type clm2_type_dec

  type(clm2_type_dec), allocatable :: CLM2_struc(:)

contains

!BOP
! 
! !ROUTINE: CLM2parms_init
! \label{CLM2parms_init}
! 
! !INTERFACE:
  subroutine CLM2parms_init(flag)

! !USES:
   use LDT_logMod,  only : LDT_verify, LDT_endrun, &
             LDT_getNextUnitNumber, LDT_releaseUnitNumber
!
! !DESCRIPTION:
!
! Community Land Model (CLM2) model parameters.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[clm2Parmssetup](\ref{clm2Parmssetup}) \newline
!    calls the registry to invoke the clm2Parms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer  :: n
   integer  :: flag
   integer  :: c,r,m,k
   integer  :: rc
   integer  :: file_status
   logical  :: file_exists
   logical  :: clm2_select 

 ! _____________________________________________________________________

   allocate(CLM2_struc(LDT_rc%nnest))

   clm2_select = .false.
   do n=1,LDT_rc%nnest
      ! - CLM2 parameters:
      call set_param_attribs(CLM2_struc(n)%clm2,"CLM2")
      if( CLM2_struc(n)%clm2%selectOpt == 1 ) then
         clm2_select = .true.
      endif
   enddo

   if( clm2_select ) then
     write(LDT_logunit,*)" - - - - - - - - - - CLM2 LSM Parameters - - - - - - - - - - - - -"

     write(LDT_logunit,*)" [INFO]  NO CLM2 Model Parameters available yet, other than"
     write(LDT_logunit,*)" [INFO]  soil color."
      
   endif


 end subroutine CLM2parms_init


 subroutine CLM2parms_writeHeader(n,ftn,dimID,monthID)
   
   integer   :: n 
   integer   :: ftn
   integer   :: dimID(3)
   integer   :: monthID
   
   integer   :: t_dimID(3)
   
   if( CLM2_struc(n)%clm2%selectOpt == 1 ) then

   endif
   
   
 end subroutine CLM2parms_writeHeader
 
  subroutine CLM2parms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    if( CLM2_struc(n)%clm2%selectOpt == 1 ) then
!      call CLM2parms_writeData(n,ftn)
!      call CLM2parms_writeData(n,ftn)
    endif

  end subroutine CLM2parms_writeData

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
   paramEntry%source = "CLM2"
   paramEntry%units ="none"
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(short_name)

  end subroutine set_param_attribs

end module CLM2_parmsMod

