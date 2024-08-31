!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module GeoWRSI_parmsMod
!BOP
!
! !MODULE: GeoWRSI_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read WRSI 
!  parameter datasets.
!
!  \subsubsection{Overview}
!  These routines in this module provides routines to read the 
!  WRSI model data, which primarily rely on the *BIL file-formatted
!  module (located in the libs directory).
!
! !REVISION HISTORY:
!
!  08 Aug 2005: Sujay Kumar; Initial implementation
!  20 Oct 2013: K. Arsenault: Updated for WRSI LSM parameters.
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_gfracMod
  use LDT_paramDataMod
  use LDT_logMod

  implicit none

  PRIVATE

!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: GeoWRSIparms_init    !allocates memory for required structures
  public :: GeoWRSIparms_writeHeader
  public :: GeoWRSIparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: GEOWRSI_struc

  type, public :: geowrsi_type_dec
     real                   :: wrsiparms_gridDesc(20)
     character*50           :: wrsiparms_proj
     character*50           :: wrsiparms_gridtransform

!  - WRSI model parameters:
     character*140  :: wrsimask_file   ! WRSI regional mask file
     character*140  :: lgp_file        ! Length of growing period file
     character*140  :: whc_file        ! Water holding capacity file
     character*140  :: sosclim_file    ! SOS climatology file
     character*140  :: wrsiclim_file   ! WRSI climatology file
     character*140  :: sos_file        ! Start-of-season (SOS) file
     character*140  :: sosanom_file    ! SOS anomaly file


! -  WRSI model-specifc:
     type(LDT_paramEntry) :: lgp         ! Length of growing period (WRSI)
     type(LDT_paramEntry) :: whc         ! Water holding capacity (WRSI)
     type(LDT_paramEntry) :: sos         ! Start of season (SOS; WRSI)
     type(LDT_paramEntry) :: sosclim     ! SOS climatology (WRSI)
     type(LDT_paramEntry) :: sosanom     ! SOS anomaly (WRSI)
     type(LDT_paramEntry) :: wrsiclim    ! Water requirement-standard index climatology (WRSI)
     type(LDT_paramEntry) :: wrsimask    ! WRSI regional landmask

  end type geowrsi_type_dec

  type(geowrsi_type_dec), pointer :: GeoWRSI_struc(:)


contains

!BOP
! 
! !ROUTINE: GeoWRSIparms_init
! \label{GeoWRSIparms_init}
! 
! !INTERFACE:
  subroutine GeoWRSIparms_init(flag)
! !USES:
    use LDT_logMod,    only : LDT_verify
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the WRSI LSM Parameter datasets.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[wrsiParms](\ref{wrsiParms}) \newline
!    calls the registry to invoke the WRSI parms setup methods. 
!  \end{description}
!
!EOP
    implicit none
    integer   :: flag
    integer   :: n
    integer   :: rc
    character*50           :: wrsiparms_proj
! __________________________

    allocate(GeoWRSI_struc(LDT_rc%nnest))
    do n=1,LDT_rc%nnest

! - WRSI parameter:
       call set_param_attribs(GeoWRSI_struc(n)%lgp,"LGP",&
            units="dekad", &
            full_name="GeoWRSIv2 Length of growing period")
       call set_param_attribs(GeoWRSI_struc(n)%whc,"WHC",&
            units="mm", &
            full_name="GeoWRSIv2 Water holding capacity")
       call set_param_attribs(GeoWRSI_struc(n)%sosclim,"SOSCLIM",&
            units="dekad", &
            full_name="GeoWRSIv2 SOS climatology")
       call set_param_attribs(GeoWRSI_struc(n)%wrsiclim,"WRSICLIM",&
            units="-", &
            full_name="GeoWRSIv2 WRSI climatology")
       call set_param_attribs(GeoWRSI_struc(n)%wrsimask,"WRSIMASK",&
            units="-", &
            full_name="GeoWRSIv2 WRSI regional mask")
       call set_param_attribs(GeoWRSI_struc(n)%sos,"SOS",&
            units="dekad", &
            full_name="GeoWRSIv2 SOS (from run-time)")
       call set_param_attribs(GeoWRSI_struc(n)%sosanom,"SOSANOM",&
            units="dekad", &
            full_name="GeoWRSIv2 SOS anomaly")
! --- 
    enddo

    if( GeoWRSI_struc(1)%wrsimask%selectOpt == 1 ) then
       write(LDT_logunit,*)" - - - - - - - - - WRSI LSM Parameters - - - - - - - - - - - -"
    endif

!    allocate(LDT_rc%wrsiparms_gridDesc(LDT_rc%nnest,20))
    do n=1,LDT_rc%nnest
    
       allocate(GeoWRSI_struc(n)%lgp%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            GeoWRSI_struc(n)%lgp%vlevels))
       allocate(GeoWRSI_struc(n)%whc%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            GeoWRSI_struc(n)%whc%vlevels))
       
       allocate(GeoWRSI_struc(n)%sos%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            GeoWRSI_struc(n)%sos%vlevels))
       
       allocate(GeoWRSI_struc(n)%sosclim%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            GeoWRSI_struc(n)%sosclim%vlevels))

       allocate(GeoWRSI_struc(n)%sosanom%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            GeoWRSI_struc(n)%sosanom%vlevels))
       
       allocate(GeoWRSI_struc(n)%wrsiclim%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            GeoWRSI_struc(n)%wrsiclim%vlevels))
       
       allocate(GeoWRSI_struc(n)%wrsimask%value(&
            LDT_rc%lnc(n),LDT_rc%lnr(n),&
            GeoWRSI_struc(n)%wrsimask%vlevels))

    enddo


!-- WRSI Parameter Config Entries:

    call ESMF_ConfigFindLabel(LDT_config,"WRSI landmask file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,GeoWRSI_struc(n)%wrsimask_file,rc=rc)
       call LDT_verify(rc,'WRSI landmask file: not specified')
    enddo
#if 0
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%wrsiparms_proj,&
         label="WRSI map projection:",rc=rc)
    call LDT_verify(rc,'WRSI map projection: option not specified in the config file')
    
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%wrsiparms_gridtransform,&
         label="WRSI spatial transform:",rc=rc)
    call LDT_verify(rc,'WRSI transform: option not specified in the config file')
    
!      call LDT_readDomainConfigSpecs("WRSI",LDT_rc%wrsiparms_proj,LDT_rc%wrsiparms_gridDesc)
#endif
    call ESMF_ConfigFindLabel(LDT_config,"WRSI length of growing period file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,GeoWRSI_struc(n)%lgp_file,rc=rc)
       call LDT_verify(rc,'WRSI length of growing period file: not specified')
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"WRSI water holding capacity file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,GeoWRSI_struc(n)%whc_file,rc=rc)
       call LDT_verify(rc,'WRSI water holding capacity file: not specified')
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"WRSI SOS file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,GeoWRSI_struc(n)%sos_file,rc=rc)
       call LDT_verify(rc,'WRSI SOS file: not specified')
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"WRSI SOS anomaly file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,GeoWRSI_struc(n)%sosanom_file,rc=rc)
       call LDT_verify(rc,'WRSI SOS anomaly file: not specified')
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"WRSI SOS climatology file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,GeoWRSI_struc(n)%sosclim_file,rc=rc)
       call LDT_verify(rc,'WRSI SOS climatology file: not specified')
    enddo
    call ESMF_ConfigFindLabel(LDT_config,"WRSI WRSI climatology file:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,GeoWRSI_struc(n)%wrsiclim_file,rc=rc)
       call LDT_verify(rc,'WRSI WRSI climatology file: not specified')
    enddo

!-- Read in WRSI Parameter Datasets:
    do n=1,LDT_rc%nnest

   !- WRSI parameters:
       write(LDT_logunit,*) "Reading WRSI mask from "//trim(GeoWRSI_struc(n)%wrsimask_file)
       call read_GeoWRSI2_landmask(&
            n,GeoWRSI_struc(n)%wrsimask%value(:,:,1))
       write(LDT_logunit,*) "Done reading WRSI mask values."

       write(LDT_logunit,*) "Reading WRSI LGP from "//trim(GeoWRSI_struc(n)%lgp_file)
       call read_GeoWRSI2_lgp(&
            n,GeoWRSI_struc(n)%lgp%value(:,:,1))
       write(LDT_logunit,*) "Done reading WRSI LGP values."

       write(LDT_logunit,*) "Reading WRSI WHC from "//trim(GeoWRSI_struc(n)%whc_file)
       call read_GeoWRSI2_whc(&
            n,GeoWRSI_struc(n)%whc%value(:,:,1)) 
       write(LDT_logunit,*) "Done reading WRSI WHC values."
       
       write(LDT_logunit,*) "Reading WRSI SOS Clim from "//trim(GeoWRSI_struc(n)%sosclim_file)
       call read_GeoWRSI2_sosclim(&
            n,GeoWRSI_struc(n)%sosclim%value(:,:,1)) 
       write(LDT_logunit,*) "Done reading WRSI SOS Climatology values."

       write(LDT_logunit,*) "Reading WRSI (WRSI) Clim from "//trim(GeoWRSI_struc(n)%wrsiclim_file)
       call read_GeoWRSI2_wrsiclim(&
            n,GeoWRSI_struc(n)%wrsiclim%value(:,:,1)) 
       write(LDT_logunit,*) "Done reading WRSI (WRSI) Climatology values."

       if( trim(GeoWRSI_struc(n)%sos_file) .eq. "none" ) then
          write(*,*)  "No SOS file assigned ..."
          GeoWRSI_struc(n)%sos%value(:,:,1) = LDT_rc%udef
       else
          write(LDT_logunit,*) "Reading WRSI SOS values from "//trim(GeoWRSI_struc(n)%sos_file)
          call read_GeoWRSI2_sos(&
               n,GeoWRSI_struc(n)%sos%value(:,:,1)) 
          write(LDT_logunit,*) "Done reading WRSI SOS values."
       endif
       if( trim(GeoWRSI_struc(n)%sosanom_file) .eq. "none" ) then
          write(*,*)  "No SOS anomaly file assigned ..."
          GeoWRSI_struc(n)%sosanom%value(:,:,1) = LDT_rc%udef
       else
          write(LDT_logunit,*) "Reading WRSI SOS Anom from "//trim(GeoWRSI_struc(n)%sosanom_file)
          call read_GeoWRSI2_sosanom(&
                n,GeoWRSI_struc(n)%sosanom%value(:,:,1)) 
          write(LDT_logunit,*) "Done reading WRSI SOS anomalies values."
       endif
    enddo

  end subroutine GeoWRSIparms_init

  subroutine GeoWRSIparms_writeHeader(n,ftn,dimID)

    integer   :: n
    integer   :: ftn
    integer   :: dimID(3)
!    integer   :: monthID

    integer   :: t_dimID(3)

    if(LDT_gfrac_struc(n)%gfrac%selectOpt.gt.0) then
       if(LDT_gfrac_struc(n)%gfracInterval.eq."monthly") then !monthly
          t_dimID(1) = dimID(1)
          t_dimID(2) = dimID(2)
!          t_dimID(3) = monthID
       endif
    end if

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         GeoWRSI_struc(n)%lgp)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         GeoWRSI_struc(n)%whc)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         GeoWRSI_struc(n)%sos)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         GeoWRSI_struc(n)%sosclim)    
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         GeoWRSI_struc(n)%sosanom)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         GeoWRSI_struc(n)%wrsiclim)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         GeoWRSI_struc(n)%wrsimask)
    
  end subroutine GeoWRSIparms_writeHeader

  subroutine GeoWRSIparms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    call LDT_writeNETCDFdata(n,ftn,GeoWRSI_struc(n)%lgp)
    call LDT_writeNETCDFdata(n,ftn,GeoWRSI_struc(n)%whc)
    call LDT_writeNETCDFdata(n,ftn,GeoWRSI_struc(n)%sos)
    call LDT_writeNETCDFdata(n,ftn,GeoWRSI_struc(n)%sosclim)
    call LDT_writeNETCDFdata(n,ftn,GeoWRSI_struc(n)%sosanom)
    call LDT_writeNETCDFdata(n,ftn,GeoWRSI_struc(n)%wrsiclim)
    call LDT_writeNETCDFdata(n,ftn,GeoWRSI_struc(n)%wrsimask)

  end subroutine GeoWRSIparms_writeData

!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name, &
                               units, full_name )

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name
   character(len=*),     optional     :: units
   character(len=*),     optional     :: full_name

   character(20) :: unit_temp
   character(100):: name_temp
! ____________________________________________________

   if(present(units)) then
      unit_temp = units
   else
      unit_temp = "none"
   endif
   if(present(full_name)) then
      name_temp = full_name
   else
      name_temp = trim(short_name)
   endif

   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = 1
   paramEntry%selectOpt = 1
   paramEntry%source = "GeoWRSIv2.0"
   paramEntry%units = unit_temp
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module GeoWRSI_parmsMod
