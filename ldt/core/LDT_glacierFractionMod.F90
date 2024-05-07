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
module LDT_glacierFractionMod
!BOP
!
! !MODULE: LDT_glacierFractionMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read glacier mask data and 
! compute glacier fraction. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  glacier fraction. 
!
! !REVISION HISTORY:
!
!  23 Jun 2020: Mahdi Navari: This subroutine is based on LDT_glacierMod.F90 
!               by Sujay Kumar. It is modified for computing glacier fraction.
!
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
  public :: LDT_glacierfrac_readParamSpecs
  public :: LDT_glacierfrac_init    !allocates memory for required structures
  public :: LDT_glacierfrac_writeHeader
  public :: LDT_glacierfrac_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_glacierfrac_struc

  type, public :: glacierfrac_type_dec

     character(len=LDT_CONST_PATH_LEN)     :: glacierfracfile
     real, allocatable :: glacierfrac_gridDesc(:)
     character*50      :: glacierfrac_proj

     character*50      :: glacierfrac_gridtransform

     ! Glacier fraction parameters
     type(LDT_paramEntry) :: glacierfrac   ! Glacier gridcell fraction

  end type glacierfrac_type_dec

  type(glacierfrac_type_dec), allocatable :: LDT_glacierfrac_struc(:)

contains

  subroutine LDT_glacierfrac_readParamSpecs
    
    character*100     :: source
    integer           :: rc
    integer           :: n
    
    allocate(LDT_glacierfrac_struc(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config,"Glacier fraction data source:",rc=rc)
    do n=1,LDT_rc%nnest
       allocate(LDT_glacierfrac_struc(n)%glacierfrac_gridDesc(20))
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_glacierfrac_struc(n)%glacierfrac,&
            "GLACIERFRAC",source)
    enddo

  end subroutine LDT_glacierfrac_readParamSpecs


!BOP
! 
! !ROUTINE: LDT_glacierfrac_init
! \label{LDT_glacierfrac_init}
! 
! !INTERFACE:
  subroutine LDT_glacierfrac_init

! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_paramOptCheckMod, only: LDT_gridOptChecks !,&

!    use LDT_logMod,    only : LDT_verify

! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the glacier fraction datasets 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[glacierfractionsetup](\ref{glacierfractionsetup}) \newline
!    calls the registry to invoke the glacierfraction setup methods. 
!  \end{description}
!
!EOP
    implicit none
    integer  :: n
    integer  :: k
    integer  :: rc
    type(LDT_fillopts)  :: glacierfrac

    logical  :: glacierfrac_select

! _____________________________________________________________

   glacierfrac_select = .false.
   do n = 1, LDT_rc%nnest
      if( LDT_glacierfrac_struc(n)%glacierfrac%selectOpt.gt.0 ) then
        glacierfrac_select = .true.
      endif
   enddo

   if( glacierfrac_select ) then
     write(LDT_logunit,*)" - - - - - - - - - - Glacier Fraction Parameters - - - - - - - - - - - - -"
   endif

   do n=1,LDT_rc%nnest
   
      if( glacierfrac_select ) then 

         call set_glacierfraction_attribs( n, LDT_glacierfrac_struc(n)%glacierfrac%source )

         LDT_glacierfrac_struc(n)%glacierfrac%vlevels = LDT_glacierfrac_struc(n)%glacierfrac%num_bins
         allocate(LDT_glacierfrac_struc(n)%glacierfrac%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              LDT_glacierfrac_struc(n)%glacierfrac%vlevels))
      endif
   enddo

  ! Glacier fraction ldt.config entries:
    if( glacierfrac_select ) then

       call ESMF_ConfigFindLabel(LDT_config,"Glacier fraction map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_glacierfrac_struc(n)%glacierfracfile,rc=rc)
          call LDT_verify(rc,'Glacier fraction map: not specified')
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Glacier fraction spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_glacierfrac_struc(n)%glacierfrac_gridtransform,&
               rc=rc)
          call LDT_verify(rc,'Glacier fraction spatial transform: option not specified in the config file')
       enddo
     ! Set units and full names:
       do n=1,LDT_rc%nnest
          LDT_glacierfrac_struc(n)%glacierfrac%units="-"
          call setGlacierFracParmsFullnames( n, "glacierfrac", &
                  LDT_glacierfrac_struc(n)%glacierfrac%source )
       enddo

       ! Read in lat/lon extents and res for sources that are binary:
       LDT_glacierfrac_struc(:)%glacierfrac_proj = "none"
       do n=1, LDT_rc%nnest
        !if( index(LDT_glacierfrac_struc(n)%glacierfrac%source,"UserDerived").eq.1 ) then

          call ESMF_ConfigGetAttribute(LDT_config,LDT_glacierfrac_struc(n)%glacierfrac_proj,&
               label="Glacier fraction map projection:",rc=rc)
          call LDT_verify(rc,'Glacier fraction map projection: option not specified in the config file')

          call LDT_readDomainConfigSpecs("Glacier fraction", &
                    LDT_glacierfrac_struc(n)%glacierfrac_proj, &
                    LDT_glacierfrac_struc(n)%glacierfrac_gridDesc)

          if( LDT_glacierfrac_struc(n)%glacierfrac_proj == "latlon" ) then

            call LDT_gridOptChecks( n,"Glacier fraction", &
                     LDT_glacierfrac_struc(n)%glacierfrac_gridtransform, &
                     LDT_glacierfrac_struc(n)%glacierfrac_proj, &
                     LDT_glacierfrac_struc(n)%glacierfrac_gridDesc(9) )
          else
             !EMK...Handle unsupported map projection.
             write(LDT_logunit,*) &
                  '[ERR] Glacier fraction map projection only supports latlont'
             call LDT_endrun()
          endif
        !endif
       enddo
    end if

!-- Read Glacier maps:
    do n = 1,LDT_rc%nnest

 
    !- Glacier fraction:
       if( LDT_glacierfrac_struc(n)%glacierfrac%selectOpt == 1 ) then
          write(LDT_logunit,*) "Reading glacier fraction map: "//&
               trim(LDT_glacierfrac_struc(n)%glacierfracfile)
          call readglacierfrac(&
               trim(LDT_glacierfrac_struc(n)%glacierfrac%source)//char(0),&
               n,LDT_glacierfrac_struc(n)%glacierfrac%value(:,:,1))

!       print *, "checking mask values for:", LDT_glacierfrac_struc(n)%glacierfrac%short_name
!       param_waterval = 0.   !<-- CONUS map with lots of 0's ... can't work properly because of 0
!       fill_option = "average"
!       fill_value = 10.
!       fill_rad = 2.
!       call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
!                LDT_rc%glacier_gridtransform(n),                          &
!                LDT_glacierfrac_struc(n)%glacierfrac%num_bins,           &
!                LDT_glacierfrac_struc(n)%glacierfrac%value, param_waterval,  &
!                LDT_glacierfrac_struc(n)%landmask2%value, fill_option, fill_value, fill_rad )
       end if

    enddo

  end subroutine LDT_glacierfrac_init

  subroutine LDT_glacierfrac_writeHeader(n,ftn,dimID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    integer    :: n
    integer    :: ftn
    integer    :: dimID(3)
    integer    :: tdimID(3)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)

    if( LDT_glacierfrac_struc(n)%glacierfrac%selectOpt.gt.0 ) then
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            LDT_glacierfrac_struc(n)%glacierfrac)
    endif
         
#endif
  end subroutine LDT_glacierfrac_writeHeader

  subroutine LDT_glacierfrac_writeData(n,ftn)

    integer     :: n 
    integer     :: ftn


    if( LDT_glacierfrac_struc(n)%glacierfrac%selectOpt.gt.0 ) then
       call LDT_writeNETCDFdata(n,ftn,LDT_glacierfrac_struc(n)%glacierfrac)
    endif

  end subroutine LDT_glacierfrac_writeData

end module LDT_glacierFractionMod
