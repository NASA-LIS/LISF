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
module SnowModel_parmsMod
!BOP
!
! !MODULE: SnowModel_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read SnowModel parameter
!  data. 
!  \subsubsection{Overview}
!  This routines in this module provides routines to read the 
!  SnowModel parameter file data.
!
! !REVISION HISTORY:
!
!  16 July 2020: K. Arsenault: SnowModel model parameters
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
  public :: SnowModelparms_init    !allocates memory for required structures
  public :: SnowModelparms_writeHeader
  public :: SnowModelparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: SnowModel_struc

  type, public :: sm_type_dec

! -  SnowModel-specific:
     integer :: Nveg         

     real           :: topoveg_gridDesc(20)
     character(140) :: topoveg_file
     character(50)  :: topoveg_gridtransform
     character(50)  :: topoveg_proj

     type(LDT_paramEntry) :: sm  ! SnowModel parameters (collective)
     type(LDT_paramEntry) :: sm_topo      ! SnowModel-specific topo map
     type(LDT_paramEntry) :: sm_vege      ! SnowModel-specific veg map

  end type sm_type_dec

  type(sm_type_dec), allocatable :: SnowModel_struc(:)

contains

!BOP
! 
! !ROUTINE: SnowModelparms_init
! \label{SnowModelparms_init}
! 
! !INTERFACE:
  subroutine SnowModelparms_init

! !USES:
   use LDT_logMod,  only : LDT_verify, LDT_endrun, &
             LDT_getNextUnitNumber, LDT_releaseUnitNumber
!
! !DESCRIPTION:
!
! Glen Liston's SnowModel model parameters.
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[smParmssetup](\ref{smParmssetup}) \newline
!    calls the registry to invoke the SnowModel Parms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer  :: n
   integer  :: c,r,m,k
   integer  :: rc
   integer  :: file_status
   logical  :: file_exists
   logical  :: sm_select 
   logical  :: check_data
   character*50  :: topoveg_proj

   type(LDT_fillopts) :: topoveg

!   integer :: veg_tiling_scheme  ! 0=SnowModel, 1=LIS

   ! _____________________________________________________________________

   allocate(SnowModel_struc(LDT_rc%nnest))

   sm_select = .false.
   do n=1,LDT_rc%nnest
      ! - SnowModel parameters:
      call set_param_attribs(SnowModel_struc(n)%sm,"SnowModel")
      if( SnowModel_struc(n)%sm%selectOpt == 1 ) then
         sm_select = .true.

         call set_param_attribs(SnowModel_struc(n)%sm_topo,"SMTOPO",&
            units="m", &
            full_name="SnowModel topography")

         call set_param_attribs(SnowModel_struc(n)%sm_vege,"SMVEG",&
            units="-", &
            full_name="SnowModel landcover classes")
      endif
   enddo

   if( sm_select ) then
     write(LDT_logunit,*)" - - - - - - - - - - SnowModel LSM Parameters - - - - - - - - - - - - -"
     write(LDT_logunit,*)"[INFO] SnowModel Model Parameters are available for"
     write(LDT_logunit,*)"       specific landcover, landmask and topographic maps"
   endif

   ! SnowModel specific topo-vege input parameter file 
   !  (Typically a Grads binary file read in the preprocess.f routine)
   check_data = .false.

   call ESMF_ConfigFindLabel(LDT_config,"Snowmodel topo-veg data source:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,SnowModel_struc(n)%sm_topo%source,rc=rc)
      call LDT_verify(rc,"Snowmodel topo-veg data source: not defined")

      if( SnowModel_struc(n)%sm_topo%source.eq."none" ) then
         SnowModel_struc(n)%sm_topo%selectOpt = 0
         SnowModel_struc(n)%sm_vege%selectOpt = 0
      endif
      if( SnowModel_struc(n)%sm_topo%selectOpt.eq.1 ) then
         check_data = .true.
         allocate(SnowModel_struc(n)%sm_topo%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SnowModel_struc(n)%sm_topo%num_bins))
         allocate(SnowModel_struc(n)%sm_vege%value(&
              LDT_rc%lnc(n),LDT_rc%lnr(n),&
              SnowModel_struc(n)%sm_vege%num_bins))
      endif
   enddo

   if(check_data) then
      write(LDT_logunit,*)"[INFO] SnowModel Topo-Vege Fields "

      call ESMF_ConfigFindLabel(LDT_config,"Snowmodel topo-veg map:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SnowModel_struc(n)%topoveg_file,rc=rc)
      enddo
      call ESMF_ConfigGetAttribute(LDT_config,topoveg_proj,&
           label="Snowmodel topo-veg map projection:",rc=rc)
      call LDT_verify(rc,'Snowmodel topo-veg map projection: option not specified in the config file')
      SnowModel_struc(:)%topoveg_proj = topoveg_proj

      ! SnowModel parameter projection check:
      do n=1,LDT_rc%nnest
         if( LDT_rc%lis_map_proj(n) .ne. SnowModel_struc(n)%topoveg_proj ) then
            write(LDT_logunit,*)"[ERR] SnowModel parameter projection should be set "
            write(LDT_logunit,*)"  the same as the LIS domain projection, for now.  "
            write(LDT_logunit,*)"  Future options will be provided to process input "
            write(LDT_logunit,*)"  SnowModel parameters for LIS model runs."
            call LDT_endrun
         endif
      enddo

      call ESMF_ConfigFindLabel(LDT_config,"Snowmodel topo-veg spatial transform:",rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,SnowModel_struc(n)%topoveg_gridtransform,&
              rc=rc)
         call LDT_verify(rc,'Snowmodel topo-veg spatial transform: option not specified in the config file')
      enddo
      do n=1,LDT_rc%nnest
         if( SnowModel_struc(n)%topoveg_gridtransform .ne. "none" ) then
            write(LDT_logunit,*)"[ERR] Currently, spatial grid transform of SnowModel parameter"
            write(LDT_logunit,*)"  is not supported but will be in future versions."
            call LDT_endrun
         endif
      enddo

      ! Read in SnowModel topo-veg "fill" options:
      topoveg%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, topoveg%filltype, &
           label="Snowmodel topo-veg fill option:",rc=rc)
      call LDT_verify(rc,"Snowmodel topo-veg fill option: option not specified in the config file")

      if( topoveg%filltype == "neighbor" ) then
         call ESMF_ConfigGetAttribute(LDT_config, topoveg%fillradius, &
              label="Snowmodel topo-veg fill radius:",rc=rc)
         call LDT_verify(rc,"Snowmodel topo-veg fill radius: option not specified in the config file")

         call ESMF_ConfigGetAttribute(LDT_config, topoveg%fillvalue, &
              label="Snowmodel topo-veg fill value:",rc=rc)
         call LDT_verify(rc,"Snowmodel topo-veg fill value: option not specified in the config file")

      elseif( topoveg%filltype == "none" ) then
         write(LDT_logunit,*) " -- 'NONE' Parameter-Mask Agreement Option Selected for Snowmodel topo-veg"
      else
         write(LDT_logunit,*) "[ERR] Fill option for Snowmodel topo-veg is not valid: ",&
              trim(topoveg%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none or neighbor "
         call LDT_endrun
      end if

      ! Read in topo-veg map (used in original SnowModel):
      do n=1,LDT_rc%nnest
         select case ( SnowModel_struc(n)%sm_topo%source )

          case ( "Grads_binary" )
            call read_SMGradsBin_topoveg(n,&
                 SnowModel_struc(n)%sm_topo%value, &
                 SnowModel_struc(n)%sm_vege%value )

!          case ( "CONSTANT" )
!            call read_CONSTANT_topoveg(n,&
!                      SnowModel_struc(n)%topoveg%value)

          case default
            write(LDT_logunit,*) "[ERR] SnowModel Topo-veg data source has not been selected."
            write(LDT_logunit,*) "  Your SnowModel LSM will not run without this parameter set."
            write(LDT_logunit,*) "  Please select one of the following: "
            write(LDT_logunit,*) " -- Grads_binary | Native | none  "
            call LDT_endrun
         end select

      end do

   endif


 end subroutine SnowModelparms_init


 subroutine SnowModelparms_writeHeader(n,ftn,dimID)
   
   integer   :: n 
   integer   :: ftn
   integer   :: dimID(3)
   
   if( SnowModel_struc(n)%sm%selectOpt == 1 ) then

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             SnowModel_struc(n)%sm_topo)

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
             SnowModel_struc(n)%sm_vege)

   endif
   
   
 end subroutine SnowModelparms_writeHeader
 
  subroutine SnowModelparms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn

    if( SnowModel_struc(n)%sm%selectOpt == 1 ) then
      call LDT_writeNETCDFdata(n,ftn,SnowModel_struc(n)%sm_topo)
      call LDT_writeNETCDFdata(n,ftn,SnowModel_struc(n)%sm_vege)
    endif

  end subroutine SnowModelparms_writeData

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
   paramEntry%source = "SnowModel"
   paramEntry%units = unit_temp
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(name_temp)

  end subroutine set_param_attribs

end module SnowModel_parmsMod

