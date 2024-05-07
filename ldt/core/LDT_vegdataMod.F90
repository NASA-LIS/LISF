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
module LDT_vegdataMod
!BOP
!
! !MODULE: LDT_vegdataMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read leaf area index (LAI)
!  and stem area index (SAI) data. 
!
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  LAI/SAI climatology data and allows the users to 
!  specify the frequency of climatology (in months). 
!  The climatological data is temporally interpolated  
!  between months to the current simulation date. 
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  12 Feb 2013: KR Arsenault: Modified for CLSM LAI data.
!
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_vegdata_readParamSpecs
  public :: LDT_vegdata_init         ! allocates memory for required structures
  public :: LDT_vegdata_writeHeader
  public :: LDT_vegdata_writeData
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_vegdata_struc

  type, public :: vegdata_type_dec

     character(len=LDT_CONST_PATH_LEN)        :: drootfile
     type(LDT_paramEntry) :: rootdepth   ! Land cover-based root-depth

  end type vegdata_type_dec

  type(vegdata_type_dec), allocatable :: LDT_vegdata_struc(:)

contains

  subroutine LDT_vegdata_readParamSpecs
    
    character*100     :: source
    integer           :: rc
    integer           :: n
    
    allocate(LDT_vegdata_struc(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config,"Rootdepth data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!       call LDT_warning(rc,"Rootdepth data source: not defined")
       call LDT_set_param_attribs(rc,LDT_vegdata_struc(n)%rootdepth,&
            "ROOTDEPTH",source)
    enddo

  end subroutine LDT_vegdata_readParamSpecs
!BOP
! 
! !ROUTINE: LDT_vegdata_init
! \label{LDT_vegdata_init}
! 
! !INTERFACE:
  subroutine LDT_vegdata_init

! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: LDT_gridOptChecks
 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the LAI/SAI datasets.
! 
!  The routines invoked are: 
!  \begin{description}
!  \end{description}
!
!EOP
   implicit none
   integer :: n 
   integer :: rc
   logical :: rootdep_select 

   rootdep_select = .false.
   do n = 1, LDT_rc%nnest
      if(LDT_vegdata_struc(n)%rootdepth%selectOpt.eq.1) then
         rootdep_select = .true.
      endif
   enddo
   
!-- Root depth:
    if( rootdep_select ) then
       call ESMF_ConfigFindLabel(LDT_config,"Root depth file:", rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_vegdata_struc(n)%drootfile,rc=rc)
       enddo
    endif

    if( rootdep_select ) then
       allocate(LDT_vegdata_struc(n)%rootdepth%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            LDT_vegdata_struc(n)%rootdepth%vlevels))
    endif
   !- Root depth:
    if( rootdep_select ) then
       call readrootdepth(&
            trim(LDT_vegdata_struc(n)%rootdepth%source)//char(0),&
            n,LDT_vegdata_struc(n)%rootdepth%value)
    endif
    
  end subroutine LDT_vegdata_init
  
  subroutine LDT_vegdata_writeHeader(n,ftn,dimID,monthID)

    integer            :: n 
    integer            :: ftn
    integer            :: dimID(3)
    integer            :: monthID

    integer            :: tdimID(3)
    
    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)

    if( LDT_vegdata_struc(n)%rootdepth%selectOpt > 0 ) then 
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_vegdata_struc(n)%rootdepth)
    end if

   end subroutine LDT_vegdata_writeHeader

  subroutine LDT_vegdata_writeData(n,ftn)


    integer          :: n 
    integer          :: ftn

    if(LDT_vegdata_struc(n)%rootdepth%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_vegdata_struc(n)%rootdepth)
    endif

  end subroutine LDT_vegdata_writeData

end module LDT_vegdataMod
