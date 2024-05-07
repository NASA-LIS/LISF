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
module LDT_irrigationMod
!BOP
!
! !MODULE: LDT_irrigationMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read irrigation fraction
!  and type data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  irrigation fraction and irrigation type data. 
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  10 Nov 2012: Kristi Arsenault: Modified for irrigation datasets.
!  12 Apr 2021: Wanshu Nie: Added support for reading in groundwater irrigation ratio  
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
  public :: LDT_irrig_readParamSpecs
  public :: LDT_irrigation_init    !allocates memory for required structures
  public :: LDT_irrigation_writeHeader
  public :: LDT_irrigation_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_irrig_struc

  type, public :: irrig_type_dec

     character(len=LDT_CONST_PATH_LEN)    :: irrigtypefile
     character(len=LDT_CONST_PATH_LEN)    :: irrigfracfile
     character(len=LDT_CONST_PATH_LEN)    :: irriggwratiofile
     real, allocatable :: irrig_gridDesc(:)
     character*50     :: irrig_proj

     character*50     :: irrigtype_gridtransform
     character*50     :: irrigfrac_gridtransform
     character*50     :: irriggwratio_gridtransform
     ! Irrigation parameters
     type(LDT_paramEntry) :: irrigtype   ! Type of irrigation
     type(LDT_paramEntry) :: irrigfrac   ! Irrigation gridcell fraction
     type(LDT_paramEntry) :: irriggwratio

  end type irrig_type_dec

  type(irrig_type_dec), allocatable :: LDT_irrig_struc(:)

contains

  subroutine LDT_irrig_readParamSpecs
    
    character*100     :: source
    integer           :: rc
    integer           :: n
    
    allocate(LDT_irrig_struc(LDT_rc%nnest))

    call ESMF_ConfigFindLabel(LDT_config,"Irrigation type data source:",rc=rc)
    do n=1,LDT_rc%nnest
       allocate(LDT_irrig_struc(n)%irrig_gridDesc(20))
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_irrig_struc(n)%irrigtype,&
            "IRRIGTYPE",source)
       if( LDT_rc%assimcropinfo(n) .and. rc /= 0 ) then
         call LDT_warning(rc,"WARNING: Irrigation type data source: not defined")
       endif
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Irrigation fraction data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!       call LDT_warning(rc,"Irrigation fraction data source: not defined")
       call LDT_set_param_attribs(rc,LDT_irrig_struc(n)%irrigfrac,&
            "IRRIGFRAC",source)
       if( LDT_rc%assimcropinfo(n) .and. rc /= 0 ) then
         call LDT_warning(rc,"WARNING: Irrigation fraction data source: not defined")
       endif
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Groundwater irrigation ratio data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_irrig_struc(n)%irriggwratio,&
            "irriggwratio",source)
       if( LDT_rc%assimcropinfo(n) .and. rc /= 0 ) then
         call LDT_warning(rc,"WARNING: groundwater irrigation ratio data source: not defined")
       endif
    enddo

  end subroutine LDT_irrig_readParamSpecs


!BOP
! 
! !ROUTINE: LDT_irrigation_init
! \label{LDT_irrigation_init}
! 
! !INTERFACE:
  subroutine LDT_irrigation_init

! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_paramOptCheckMod, only: LDT_gridOptChecks !,&
!           LDT_irrigationOptChecks
!    use LDT_logMod,    only : LDT_verify

! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the irrigation fraction datasets 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[irrigationsetup](\ref{irrigationsetup}) \newline
!    calls the registry to invoke the irrigation setup methods. 
!  \end{description}
!
!EOP
    implicit none
    integer  :: n
    integer  :: k
    integer  :: rc
    type(LDT_fillopts)  :: irrigtype
    type(LDT_fillopts)  :: irrigfrac
    type(LDT_fillopts)  :: irriggwratio

    logical  :: irrigfrac_select
    logical  :: irrigtype_select
    logical  :: irriggwratio_select

! _____________________________________________________________

   irrigfrac_select = .false.
   irrigtype_select = .false.
   irriggwratio_select = .false.
   do n = 1, LDT_rc%nnest
      if( LDT_irrig_struc(n)%irrigtype%selectOpt.gt.0 ) then
        irrigtype_select = .true.
      endif
      if( LDT_irrig_struc(n)%irrigfrac%selectOpt.gt.0 ) then
        irrigfrac_select = .true.
      endif
      if( LDT_irrig_struc(n)%irriggwratio%selectOpt.gt.0 ) then
        irriggwratio_select = .true.
      endif
   enddo

   if( irrigfrac_select .or. irrigtype_select .or. irriggwratio_select) then
     write(LDT_logunit,*)" - - - - - - - - - - Irrigation Parameters - - - - - - - - - - - - -"
   endif

   do n=1,LDT_rc%nnest
   
      if( irrigtype_select ) then 

         call set_irrigation_attribs( n, LDT_irrig_struc(n)%irrigtype%source )

         LDT_irrig_struc(n)%irrigtype%vlevels = LDT_irrig_struc(n)%irrigtype%num_bins
         allocate(LDT_irrig_struc(n)%irrigtype%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              LDT_irrig_struc(n)%irrigtype%vlevels))       
      endif
      if( irrigfrac_select ) then
         LDT_irrig_struc(n)%irrigfrac%vlevels = LDT_irrig_struc(n)%irrigfrac%num_bins
         allocate(LDT_irrig_struc(n)%irrigfrac%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              LDT_irrig_struc(n)%irrigfrac%vlevels))
      endif

      if( irriggwratio_select ) then

         call set_irriggwratio_attribs( n, LDT_irrig_struc(n)%irriggwratio%source )
         LDT_irrig_struc(n)%irriggwratio%vlevels = LDT_irrig_struc(n)%irriggwratio%num_bins
         allocate(LDT_irrig_struc(n)%irriggwratio%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              LDT_irrig_struc(n)%irriggwratio%vlevels))
      endif
   enddo

  ! Irrigation type ldt.config entries:
    if( irrigtype_select ) then 
       call ESMF_ConfigFindLabel(LDT_config,"Irrigation type map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_irrig_struc(n)%irrigtypefile,rc=rc)
          call LDT_verify(rc,'Irrigation type map: not specified')
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Irrigation type spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_irrig_struc(n)%irrigtype_gridtransform,&
               rc=rc)
          call LDT_verify(rc,'Irrigation type spatial transform: option not specified in the config file')
       enddo
     ! Set units and full names:
       do n=1,LDT_rc%nnest
          LDT_irrig_struc(n)%irrigtype%units="-"
          call setIrrigParmsFullnames( n, "irrigtype", &
                  LDT_irrig_struc(n)%irrigtype%source )
       enddo
    end if
  ! Irrigation fraction ldt.config entries:
    if( irrigfrac_select ) then

       call ESMF_ConfigFindLabel(LDT_config,"Irrigation fraction map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_irrig_struc(n)%irrigfracfile,rc=rc)
          call LDT_verify(rc,'Irrigation fraction map: not specified')
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Irrigation fraction spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_irrig_struc(n)%irrigfrac_gridtransform,&
               rc=rc)
          call LDT_verify(rc,'Irrigation fraction spatial transform: option not specified in the config file')
       enddo
     ! Set units and full names:
       do n=1,LDT_rc%nnest
          LDT_irrig_struc(n)%irrigfrac%units="-"
          call setIrrigParmsFullnames( n, "irrigfrac", &
                  LDT_irrig_struc(n)%irrigfrac%source )
       enddo

! Groundwater irrigation ratio ldt.config entries:
    if( irriggwratio_select ) then
       call ESMF_ConfigFindLabel(LDT_config,"Groundwater irrigation ratio map:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_irrig_struc(n)%irriggwratiofile,rc=rc)
          call LDT_verify(rc,'Groundwater irrigation ratio map: not specified')
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Groundwater irrigation ratio spatial transform:",rc=rc)
       do n=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_irrig_struc(n)%irriggwratio_gridtransform,&
               rc=rc)
          call LDT_verify(rc,'Groundwater irrigation ratio spatial transform: option not specified in the config file')
       enddo
     ! Set units and full names:
       do n=1,LDT_rc%nnest
          LDT_irrig_struc(n)%irriggwratio%units="-"
          call setIrrigParmsFullnames( n, "irriggwratio", &
                  LDT_irrig_struc(n)%irriggwratio%source )
       enddo

    end if

       ! Read in lat/lon extents and res for sources that are binary:
       LDT_irrig_struc(:)%irrig_proj = "none"
       do n=1, LDT_rc%nnest
        if( index(LDT_irrig_struc(n)%irrigfrac%source,"UserDerived").eq.1 ) then

          call ESMF_ConfigGetAttribute(LDT_config,LDT_irrig_struc(n)%irrig_proj,&
               label="Irrigation fraction map projection:",rc=rc)
          call LDT_verify(rc,'Irrigation fraction map projection: option not specified in the config file')

          call LDT_readDomainConfigSpecs("Irrigation fraction", &
                    LDT_irrig_struc(n)%irrig_proj, &
                    LDT_irrig_struc(n)%irrig_gridDesc)

          if( LDT_irrig_struc(n)%irrig_proj == "latlon" ) then

            call LDT_gridOptChecks( n,"Irrigation fraction", &
                     LDT_irrig_struc(n)%irrigfrac_gridtransform, &
                     LDT_irrig_struc(n)%irrig_proj, &
                     LDT_irrig_struc(n)%irrig_gridDesc(9) )
          else
             !EMK...Handle unsupported map projection.
             write(LDT_logunit,*) &
                  '[ERR] Irrigation fraction map projection only supports latlont'
             call LDT_endrun()
          endif
        endif
       enddo
    end if

!-- Read Irrigation maps:
    do n = 1,LDT_rc%nnest

    !- Irrigation type:
       if( LDT_irrig_struc(n)%irrigtype%selectOpt == 1 ) then
          write(LDT_logunit,*) "Reading irrigation type: "//&
               trim(LDT_irrig_struc(n)%irrigtypefile)
          call readirrigtype( trim(LDT_irrig_struc(n)%irrigtype%source)//char(0),&
                n, LDT_irrig_struc(n)%irrigtype%value,   &
                LDT_irrig_struc(n)%irrigtype%num_bins )
       endif
 
    !- Irrigation fraction:
       if( LDT_irrig_struc(n)%irrigfrac%selectOpt == 1 ) then
          write(LDT_logunit,*) "Reading irrigation map: "//&
               trim(LDT_irrig_struc(n)%irrigfracfile)
          call readirrigfrac(&
               trim(LDT_irrig_struc(n)%irrigfrac%source)//char(0),&
               n,LDT_irrig_struc(n)%irrigfrac%value(:,:,1))

!       print *, "checking mask values for:", LDT_irrig_struc(n)%irrigfrac%short_name
!       param_waterval = 0.   !<-- CONUS map with lots of 0's ... can't work properly because of 0
!       fill_option = "average"
!       fill_value = 10.
!       fill_rad = 2.
!       call LDT_contIndivParam_Fill(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
!                LDT_rc%irrig_gridtransform(n),                          &
!                LDT_irrig_struc(n)%irrigfrac%num_bins,           &
!                LDT_irrig_struc(n)%irrigfrac%value, param_waterval,  &
!                LDT_irrig_struc(n)%landmask2%value, fill_option, fill_value, fill_rad )
       end if

    !- Groundwater irrigation ratio:
       if( LDT_irrig_struc(n)%irriggwratio%selectOpt == 1 ) then
          write(LDT_logunit,*) "Reading groundwater irrigation ratio: "//&
               trim(LDT_irrig_struc(n)%irriggwratiofile)
          write(LDT_logunit,*) "Groundwater irrigation ratio data source: "//&
               trim(LDT_irrig_struc(n)%irriggwratio%source)
          call readirriggwratio( trim(LDT_irrig_struc(n)%irriggwratio%source)//char(0),&
                n, LDT_irrig_struc(n)%irriggwratio%value,   &
                LDT_irrig_struc(n)%irriggwratio%num_bins )
       endif

    enddo

  end subroutine LDT_irrigation_init

  subroutine LDT_irrigation_writeHeader(n,ftn,dimID)

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

    if( LDT_irrig_struc(n)%irrigtype%selectOpt.gt.0 ) then
       call LDT_verify(nf90_def_dim(ftn,'irrigtypes',&
            LDT_irrig_struc(n)%irrigtype%vlevels,tdimID(3)))

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
                LDT_irrig_struc(n)%irrigtype)
    endif
          
    if( LDT_irrig_struc(n)%irrigfrac%selectOpt.gt.0 ) then
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            LDT_irrig_struc(n)%irrigfrac)
    endif
         
    if( LDT_irrig_struc(n)%irriggwratio%selectOpt.gt.0 ) then
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            LDT_irrig_struc(n)%irriggwratio)
    endif

#endif
  end subroutine LDT_irrigation_writeHeader

  subroutine LDT_irrigation_writeData(n,ftn)

    integer     :: n 
    integer     :: ftn

    if( LDT_irrig_struc(n)%irrigtype%selectOpt.gt.0 ) then
       call LDT_writeNETCDFdata(n,ftn,LDT_irrig_struc(n)%irrigtype)
    endif

    if( LDT_irrig_struc(n)%irrigfrac%selectOpt.gt.0 ) then
       call LDT_writeNETCDFdata(n,ftn,LDT_irrig_struc(n)%irrigfrac)
    endif

    if( LDT_irrig_struc(n)%irriggwratio%selectOpt.gt.0 ) then
       call LDT_writeNETCDFdata(n,ftn,LDT_irrig_struc(n)%irriggwratio)
    endif

  end subroutine LDT_irrigation_writeData

end module LDT_irrigationMod
