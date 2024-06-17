!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module AquaCrop_parmsMod
!BOP
!
! !MODULE: AquaCrop_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read greenness fraction
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the AquaCrop
!  crop type from the AC_Crop.Inventory
!
! !REVISION HISTORY:
!
!  10 May 2024; Michel Becthold, Louise Busschaert, initial implementation
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
  public :: AquaCropParms_init    !allocates memory for required structures
  public :: AquaCropParms_writeHeader
  public :: AquaCropParms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: AquaCrop_struc

  type, public :: aquacrop_type_dec
      ! -  AquaCrop LSM-specific:
      type(LDT_paramEntry) :: cropt   ! crop type
      type(LDT_paramEntry) :: comp_size ! compartment size
      integer :: nlayers !  number of soil layers
      real :: lthickness(5) ! thickness of layers, max 5 layers for AC
      integer :: max_comp ! 12 by default
  end type aquacrop_type_dec

  type(aquacrop_type_dec), allocatable :: AquaCrop_struc(:)



contains

  subroutine AquaCropParms_init(flag)
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the AquaCrop fraction datasets 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[acParmssetup](\ref{acParmssetup}) \newline
!    calls the registry to invoke the noahParms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer  :: flag
   integer  :: n,i,c,r,m
   integer  :: rc

! _____________________________________________________________________

  allocate(AquaCrop_struc(LDT_rc%nnest))
    do n=1,LDT_rc%nnest
    AquaCrop_struc(n)%max_comp = 12
    ! Set crop type
      call set_param_attribs(Aquacrop_struc(n)%cropt, "AC_CROPT",&
          units="-", &
          full_name="Aquacrop crop type")
      call ESMF_ConfigFindLabel(LDT_config,"AquaCrop crop type data source:",rc=rc)
      call ESMF_ConfigGetAttribute(LDT_config,AquaCrop_struc(n)%cropt%source,rc=rc)
      call LDT_verify(rc,"AquaCrop crop type data source: not defined")
      allocate(Aquacrop_struc(n)%cropt%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              Aquacrop_struc(n)%cropt%vlevels))
      select case (AquaCrop_struc(n)%cropt%source)
        case( "CONSTANT" )
          call read_CONSTANT_AC_crop(&
                    n,AquaCrop_struc(n)%cropt%value(:,:,1))
        case default
          write(LDT_logunit,*) "[WARN] crop type data source not valid for AquaCrop."
          write(LDT_logunit,*) "  Please select: CONSTANT"
          write(LDT_logunit,*) "Program stopping ..."
          call LDT_endrun
      end select
    ! End crop type
  
    ! Define compartment size
    
      call set_param_attribs(Aquacrop_struc(n)%comp_size, "AC_comp_size",&
          units="m", &
          full_name="Aquacrop compartment size")
      Aquacrop_struc(n)%comp_size%vlevels = AquaCrop_struc(n)%max_comp
      Aquacrop_struc(n)%comp_size%num_bins = AquaCrop_struc(n)%max_comp
      allocate(Aquacrop_struc(n)%comp_size%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              Aquacrop_struc(n)%comp_size%vlevels))
      call define_AC_compartments(n, AquaCrop_struc(n)%comp_size%value(:,:,:))
    enddo ! End nest
    
  end subroutine AquaCropParms_init

 subroutine AquaCropParms_writeHeader(n,ftn,dimID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif

    integer   :: n,i 
    integer   :: ftn
    integer   :: dimID(3)

    integer   :: ndimID(3)  ! 3D, vlevel>1
    character(25) :: str

    ndimID(1) = dimID(1)
    ndimID(2) = dimID(2)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            Aquacrop_struc(n)%cropt)
    call LDT_verify(nf90_def_dim(ftn,'AC_max_compartments',&
          AquaCrop_struc(n)%max_comp,ndimID(3)))
    ndimID(1) = dimID(1)
    ndimID(2) = dimID(2)
    call LDT_writeNETCDFdataHeader(n,ftn,ndimID,&
            Aquacrop_struc(n)%comp_size)
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOIL_LAYERS", &
        AquaCrop_struc(n)%nlayers))
    do i=1,AquaCrop_struc(n)%nlayers
      write (str, '(i0)') i
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"THICKNESS_LAYER_"//trim(str),&
                      AquaCrop_struc(n)%lthickness(i)))
    enddo
#endif

  end subroutine AquaCropParms_writeHeader

  subroutine AquaCropParms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    call LDT_writeNETCDFdata(n,ftn,AquaCrop_struc(n)%cropt)
    call LDT_writeNETCDFdata(n,ftn,AquaCrop_struc(n)%comp_size)

  end subroutine AquaCropParms_writeData


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
   paramEntry%source = "AquaCrop"
   paramEntry%units = unit_temp
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module AquaCrop_parmsMod
