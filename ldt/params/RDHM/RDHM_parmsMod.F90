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
module RDHM_parmsMod
!BOP
!
! !MODULE: RDHM_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read RDHM parameter
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  RDHM parameter file data.
!
! !REVISION HISTORY:
!
!  04 Nov 2013: K. Arsenault: Added layers for RDHM model
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
  use LDT_xmrg_reader
  use SACHTET_parmsMod
  use Snow17_parmsMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: RDHMparms_init    !allocates memory for required structures
  public :: RDHMparms_writeHeader
  public :: RDHMparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: RDHM_struc

  type, public :: rdhm_type_dec

! - RDHM 3.5.6
     character*140 :: rdhmconsts_table
     real          :: rdhm_undef

! -  RDHM-specific:
     type(LDT_paramEntry) :: rdhm356  ! RDHM v.3.5.6 parameters (collective)

  end type rdhm_type_dec

  type(rdhm_type_dec), allocatable :: RDHM_struc(:)

contains

!BOP
! 
! !ROUTINE: RDHMparms_init
! \label{RDHMparms_init}
! 
! !INTERFACE:
  subroutine RDHMparms_init(flag)

! !USES:
   use LDT_logMod,    only : LDT_verify, LDT_endrun, &
             LDT_getNextUnitNumber, LDT_releaseUnitNumber
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! Research Distributed Hydrological Model (RDHM), including
! the Sacramento Soil Moisture Accounting (SAC-SMA) 
! Heat transfer and evapotranspiration (HT-ET) and SNOW-17 
! parameter datasets.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[rdhmParmssetup](\ref{rdhmParmssetup}) \newline
!    calls the registry to invoke the rdhmParms setup methods. 
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
   logical  :: rdhm_select 
!   type(ESMF_Config)  :: rdhmconsts_table
!   type(LDT_fillopts) :: rdhm

   ! _____________________________________________________________________

   allocate(RDHM_struc(LDT_rc%nnest))

   rdhm_select = .false.
   do n=1,LDT_rc%nnest
      ! - RDHM (v3.5.6) parameters:
      call set_param_attribs(RDHM_struc(n)%rdhm356,"RDHM356")
      if( RDHM_struc(n)%rdhm356%selectOpt == 1 ) then
         rdhm_select = .true.
      endif
   enddo

   if( rdhm_select ) then
     write(LDT_logunit,*)" - - - - - - - - - - RDHM LSM Parameters - - - - - - - - - - - - -"

     call SACHTETParms_init( flag )
     call Snow17Parms_init( flag )

#if 0
   !- Load RDHM CONSTANTS input table file (filepath read-in from ldt.config file)
      call ESMF_ConfigFindLabel(LDT_config,"RDHM356 constants table:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,RDHM_struc(n)%rdhmconsts_table,rc=rc)
         call LDT_verify(rc,'RDHM356 constants table: not specified')

         inquire(file=trim(RDHM_struc(n)%rdhmconsts_table), exist=file_exists)
         if( .not. file_exists ) then
            write(LDT_logunit,*) "[ERR] RDHM Parameter Constants Table ",&
                 trim(RDHM_struc(n)%rdhmconsts_table)," does not exist."
            call LDT_endrun
         endif
         write(LDT_logunit,*) "Reading in RDHM Parameter Constants Table Entries"
         rdhmconsts_table = ESMF_ConfigCreate(rc=rc)
         call ESMF_ConfigLoadFile(rdhmconsts_table, &
              trim(RDHM_struc(n)%rdhmconsts_table), rc=rc)
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"RDHM356 universal undefined value:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,RDHM_struc(n)%rdhm_undef,rc=rc)
         call LDT_verify(rc,'RDHM356 universal undefined value: not specified')
      enddo
#endif

      ! -----
      
   endif

 end subroutine RDHMparms_init


 subroutine RDHMparms_writeHeader(n,ftn,dimID,monthID)
   
   integer   :: n 
   integer   :: ftn
   integer   :: dimID(3)
   integer   :: monthID
   
   integer   :: t_dimID(3)
   
   if( RDHM_struc(n)%rdhm356%selectOpt == 1 ) then

     call SACHTETparms_writeHeader(n,ftn,dimID,monthID)

     call Snow17Parms_writeHeader(n,ftn,dimID)

   endif
   
   
 end subroutine RDHMparms_writeHeader
 
  subroutine RDHMparms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    if( RDHM_struc(n)%rdhm356%selectOpt == 1 ) then
      call SACHTETparms_writeData(n,ftn)
      call Snow17parms_writeData(n,ftn)
    endif

  end subroutine RDHMparms_writeData

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
   paramEntry%source = "RDHM.3.5.6"
   paramEntry%units ="none"
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(short_name)

  end subroutine set_param_attribs

end module RDHM_parmsMod

