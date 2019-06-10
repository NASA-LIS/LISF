!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v1.0
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
!
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
  public :: LDT_irrig_readParamSpecs
  public :: LDT_irrigation_init    !allocates memory for required structures
  public :: LDT_irrigation_writeHeader
  public :: LDT_irrigation_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_irrig_struc

  type, public :: irrig_type_dec

     character*140    :: irrigtypefile
     character*140    :: irrigfracfile
     real, allocatable :: irrig_gridDesc(:)
     character*50     :: irrig_proj

     character*50     :: irrigtype_gridtransform
     character*50     :: irrigfrac_gridtransform

     character*50     :: irrigfrac_typeopt

     ! Irrigation parameters
     type(LDT_paramEntry) :: irrigtype   ! Type of irrigation
     type(LDT_paramEntry) :: irrigfrac   ! Irrigation gridcell fraction

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
!    use LDT_logMod,    only : LDT_verify
!    use LDT_paramOptCheckMod, only: LDT_irrigationOptChecks, &
!                       LDT_gridOptChecks

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

    logical  :: irrigfrac_select
    logical  :: irrigtype_select

! _____________________________________________________________

   irrigfrac_select = .false.
   irrigtype_select = .false.
   do n = 1, LDT_rc%nnest
      if( LDT_irrig_struc(n)%irrigtype%selectOpt.gt.0 ) then
        irrigtype_select = .true.
      endif
      if( LDT_irrig_struc(n)%irrigfrac%selectOpt.gt.0 ) then
        irrigfrac_select = .true.
      endif
   enddo

   if( irrigfrac_select .or. irrigtype_select ) then
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
       ! Option to select irrigation type, if GRIPC selected: 
       if( LDT_irrig_struc(1)%irrigfrac%source == "GRIPC" )then
         call ESMF_ConfigFindLabel(LDT_config,"Irrigation fraction type option:",rc=rc)
         do n=1,LDT_rc%nnest
            call ESMF_ConfigGetAttribute(LDT_config,&
                 LDT_irrig_struc(n)%irrigfrac_typeopt,&
                 default="sprinkler",&
                 rc=rc)
            call LDT_verify(rc,'Irrigation fraction type option: option not specified in the config file')
            select case ( LDT_irrig_struc(n)%irrigfrac_typeopt )
             case( "sprinkler", "paddy", "sprinkler+paddy" )
               write(LDT_logunit,*) "[INFO]  For GRIPC irrigation fraction, " 
               write(LDT_logunit,*) "  the type of irrigation selected is: ",&
                    trim(LDT_irrig_struc(n)%irrigfrac_typeopt)
             case default
               write(LDT_logunit,*) "[ERR]  Incorrect GRIPC irrigation fraction type selected."
               write(LDT_logunit,*) "  Check your ldt.config file and enter one of these options:"
               write(LDT_logunit,*) "    sprinkler | paddy | sprinkler+paddy "
               call LDT_endrun
            end select  
         enddo
       end if
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

  end subroutine LDT_irrigation_writeData

end module LDT_irrigationMod
