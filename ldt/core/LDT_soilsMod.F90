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
module LDT_soilsMod
!BOP
!
! !MODULE: LDT_soilsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read various sources of
!  soil parameter data. 
! 
!  \subsubsection{Overview}
!  This module provides routines for reading and manipulating various 
!  parameters related to soil properties. The following ldtt of 
!  soil parameters are currently supported. 
!  \begin{description}
!   \item[sand, silt, clay fractions]
!   \item[soil color data]
!   \item[soil texture data]
!   \item[soil porosity data]
!   \item[saturated hydraulic conductivity data]
!   \item[thermal conductivity data]
!   \item[b parameter data]
!   \item[hydraulic conductivity data]
!   \item[quartz data]
!   \item[wilting point wetness data]
!   \item[soil depth data]
!   \item[bedrock depth data]
!   \item[hydrological soils group data]
!  \end{description}
!
! !REVISION HISTORY:
!
!  21 Oct 2005: Sujay Kumar; Initial implementation
!  21 Nov 2012: K. Arsenault; Include additional soil parameters
!
  use ESMF
#if ( defined SPMD )
  use mpi
#endif
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_coreMod
  use LDT_logMod
  use LDT_paramMaskCheckMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LDT_soils_readParamSpecs
  public :: LDT_soils_init  ! initializes data structures and memory
  public :: LDT_soils_writeHeader
  public :: LDT_soils_writeData
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LDT_soils_struc

  type, public :: soils_type_dec
     character*50         :: soilclr_proj
     character*50         :: soilclr_gridtransform
     real                 :: soilclr_gridDesc(20)
   ! Soil hydraulic parameters
     type(LDT_paramEntry) :: color          ! Soil color: soil albedo parameter
     type(LDT_paramEntry) :: porosity       ! Soil porosity
     type(LDT_paramEntry) :: soildepth      ! Soil depth 
     type(LDT_paramEntry) :: hsg            ! Hydraulic soil groups (STATSGO v1)
     type(LDT_paramEntry) :: bulkdens       ! Bulk density (STATSGO v1)
     type(LDT_paramEntry) :: watercap       ! Water capacity (STATSGO v1)
     type(LDT_paramEntry) :: rockvol        ! Rock volume (STATSGO v1)
     type(LDT_paramEntry) :: rockfragcls    ! Rock fragment class (STATSGO v1)
     type(LDT_paramEntry) :: permrate       ! Permeability rate (STATSGO v1)
     type(LDT_paramEntry) :: texture_nlyrs  ! Soil texture - Number of layers (STATSGO v1)

  end type soils_type_dec

  type(soils_type_dec), allocatable :: LDT_soils_struc(:)

!EOP

!BOP 
! 
! !ROUTINE: LDT_soils_init
! \label{LDT_soils_init}
! 
! !INTERFACE:
  interface LDT_soils_init
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure soils_init_LIS
     module procedure soils_init_LISHydro
! 
! !DESCRIPTION:
! This interface provides the initialization routines for processing
! soil parameters both in the standard LIS preprocessing mode 
! as well as in the LISHydro(WRFHydro) preprocessing mode.
!EOP 
  end interface


!BOP 
! 
! !ROUTINE: LDT_soils_writeHeader 
! \label{LDT_soils_writeHeader}
! 
! !INTERFACE:
  interface LDT_soils_writeHeader
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure soils_writeHeader_LIS
     module procedure soils_writeHeader_LISHydro
! 
! !DESCRIPTION:
! This interface provides routines for writing NETCDF header for soils data
! both in the standard LIS preprocessing mode as well as in 
! the LISHydro(WRFHydro) preprocessing mode.
!EOP 
  end interface


!BOP 
! 
! !ROUTINE: LDT_soils_writeData
! \label{LDT_soils_writeData}
! 
! !INTERFACE:
  interface LDT_soils_writeData
! !PRIVATE MEMBER FUNCTIONS: 
     module procedure soils_writeData_LIS
     module procedure soils_writeData_LISHydro
! 
! !DESCRIPTION:
! This interface provides routines for writing soils data both 
! in the standard LIS preprocessing mode as well as in the LISHydro(WRFHydro) 
! preprocessing mode.
!EOP 
  end interface


 contains

!BOP
! 
! !ROUTINE: LDT_soils_readParamSpecs
! \label{LDT_soils_readParamSpecs}
! 
! !INTERFACE: 
   subroutine LDT_soils_readParamSpecs
! 
! !DESCRIPTION: 
!  This subroutine reads the soil parameter specifications 
!  from the config file.  
!EOP

     character*100    :: source
     integer          :: rc
     integer          :: n
     
     rc = 0

     allocate(LDT_soils_struc(LDT_rc%nnest))

  !- Soil textures at multiple soil depths:
  ! - Currently only works with SAC-HTET/RDHM:
     do n=1,LDT_rc%nnest
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%texture_nlyrs,&
             "TEXTURE_NLYRS",LDT_LSMparam_struc(n)%texture%source)
     end do

  !- Soil color data source check:
     call ESMF_ConfigFindLabel(LDT_config,"Soil color data source:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%color,&
             "SOILCOLOR",source)
     enddo
   ! LSM-required parameter check:
     if( index(LDT_rc%lsm,"CLM") == 1 ) then
       if( rc /= 0 ) then
          call LDT_warning(rc,"[WARN] Soil color data source: not defined")
       endif
     endif
     
  !- Porosity data source check:
     call ESMF_ConfigFindLabel(LDT_config,"Porosity data source:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!        call LDT_warning(rc,"Porosity data source: not defined")
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%porosity,&
             "POROSITY",source)
     enddo
   ! LSM-required parameter check:
     if( index(LDT_rc%lsm,"CLSM") == 1 ) then
       if( rc /= 0 ) then
          call LDT_warning(rc,"[WARN] Porosity data source: not defined")
       endif
     endif
     
   ! Hydrologic Soil Groups (STATSGO v1 only at this time):
     call ESMF_ConfigFindLabel(LDT_config,"Hydrologic soil group source:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!        call LDT_warning(rc,"Hydrologic soil group source: not defined")
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%hsg,&
             "HSG",source)
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Bulk density data source:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!        call LDT_warning(rc,"Bulk density data source: not defined")
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%bulkdens,&
             "BULKDENSITY",source)
     enddo

     call ESMF_ConfigFindLabel(LDT_config,"Water capacity data source:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!        call LDT_warning(rc,"Water capacity data source: not defined")
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%watercap,&
            "WATERCAPACITY",source)
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Rock volume data source:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!        call LDT_warning(rc,"Rock volume data source: not defined")
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%rockvol,&
             "ROCKVOL",source)
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Rock frag class data source:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!        call LDT_warning(rc,"Rock frag class data source: not defined")
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%rockfragcls,&
             "ROCKFRAGCLS",source)
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Permeability data source:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!        call LDT_warning(rc,"Permeability data source: not defined")
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%permrate,&
             "PERMRATE",source)
     enddo
     
     call ESMF_ConfigFindLabel(LDT_config,"Soildepth data source:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!        call LDT_warning(rc,"Soildepth data source: not defined")
        call LDT_set_param_attribs(rc,LDT_soils_struc(n)%soildepth,&
             "SOILDEPTH",source)
     enddo


   end subroutine LDT_soils_readParamSpecs

!BOP
! 
! !ROUTINE: soils_init_LIS
! \label{soils_init_LIS}
! 
! !INTERFACE:
  subroutine soils_init_LIS()

! !USES:
    use LDT_coreMod,   only : LDT_rc, LDT_config
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_paramOptCheckMod, only: LDT_soilsOptChecks, & 
                       LDT_gridOptChecks
! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! soils datasets and reads the soils data based on the 
! choice of options specified in the ldt configuration, 
! in the standard LIS preprocessing mode. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readsoiltexture](\ref{readsoiltexture}) \newline
!    invokes the generic method in the registry to read
!    the soil texture data
!   \item[readsoilfrac](\ref{readsoilfrac}) \newline
!    invokes the generic method in the registry to read
!    the sand, clay and silt fraction data
!   \item[readcolor](\ref{readcolor}) \newline
!    invokes the generic method in the registry to read
!    the soil color data
!   \item[readporosity](\ref{readporosity}) \newline
!    invokes the generic method in the registry to read
!    the porosity data
!   \item[readsoildepth](\ref{readsoildepth}) \newline
!    invokes the generic method in the registry to read
!    the soil depth parameter data
!   \item[readbdrckdpth](\ref{readsoildepth}) \newline
!    invokes the generic method in the registry to read
!    the depth to bedrock parameter data
!   \item[readwpwet](\ref{readwpwet}) \newline
!    invokes the generic method in the registry to read
!    the wilting poing wetness parameter data
!   \item[readquartz](\ref{readquartz}) \newline
!    invokes the generic method in the registry to read
!    the quartz data
!   \item[readhsg](\ref{readhsg}) \newline
!    invokes the generic method in the registry to read
!    the hydrological soils group data
!  \end{description}
!
!EOP
    implicit none
    integer  :: n, i, c, r
    integer  :: rc
    character*50        :: soilclr_proj
    real, allocatable   :: soilfrac_array(:,:,:)
    character(50)       :: tile_temp
    real, allocatable   :: soilclr_gridDesc(:,:)
    type(LDT_fillopts)  :: soiltext
    type(LDT_fillopts)  :: soilfrac
    type(LDT_fillopts)  :: soilcolor
    type(LDT_fillopts)  :: soildepth
    type(LDT_fillopts)  :: rootdepth
    type(LDT_fillopts)  :: porosity
    logical             :: soil_select
    logical             :: check_data
    
! _____________________________________________________________________________
    

  ! Plugin:  Soil Texture Maps - Set number of texture types (attributes):
    if(LDT_LSMparam_struc(1)%texture%selectOpt.eq.1) then 
       call setTextureattribs(trim(LDT_LSMparam_struc(1)%texture%source)//char(0))  
    endif

  ! Plugin:  HSG Maps - Set number of Hydrologic Soil Group types (attributes):
    if(LDT_soils_struc(1)%hsg%selectOpt.eq.1) then
       call setHSGattribs(trim(LDT_soils_struc(1)%hsg%source)//char(0))
    endif

    soil_select = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%texture%selectOpt == 1 .or. &
           LDT_LSMparam_struc(n)%sand%selectOpt   == 1 .or. &
           LDT_soils_struc(n)%porosity%selectOpt  == 1 ) then 
         soil_select = .true.
       end if
    enddo

    if(soil_select) &
         write(LDT_logunit,*)" - - - - - - - - - - Soils Parameters - - - - - - - - - - - - - -"

    
    allocate(LDT_rc%soil_gridDesc(LDT_rc%nnest,20))
    allocate(LDT_rc%soiltext_gridDesc(LDT_rc%nnest,20))
    allocate(LDT_rc%hsg_gridDesc(LDT_rc%nnest,20))

    allocate(soilclr_gridDesc(LDT_rc%nnest,20))

    do n=1, LDT_rc%nnest

       LDT_rc%soil_classification(n) = "Soil texture not selected"
    !- Soil texture only:
       if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then

          LDT_LSMparam_struc(n)%texture%vlevels = &
                 LDT_LSMparam_struc(n)%texture%num_bins
          allocate(LDT_LSMparam_struc(n)%texture%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%texture%num_bins))

!          if( (index(LDT_rc%lsm,"RDHM") == 1 .or. &
!              index(LDT_rc%lsm,"SACHTET") == 1 ) .and. &
!              LDT_rc%create_soilparms_option == "create" ) then
          allocate(LDT_soils_struc(n)%texture_nlyrs%value(&
                   LDT_rc%lnc(n),LDT_rc%lnr(n),&
                   LDT_soils_struc(n)%texture_nlyrs%vlevels))
!          endif

       endif

    !- Activated soil fractions (combined for sand/silt/clay):
       if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1) then

          call ESMF_ConfigFindLabel(LDT_config,"Soil fraction number of bands:",rc=rc)
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_LSMparam_struc(n)%sand%num_bins,rc=rc)
          call LDT_verify(rc,"Soil fraction number of bands: not defined")

          LDT_LSMparam_struc(n)%clay = LDT_LSMparam_struc(n)%sand
          LDT_LSMparam_struc(n)%silt = LDT_LSMparam_struc(n)%sand
          LDT_LSMparam_struc(n)%gravel = LDT_LSMparam_struc(n)%sand
          LDT_LSMparam_struc(n)%soilsfgrd = LDT_LSMparam_struc(n)%sand

          LDT_LSMparam_struc(n)%sand%vlevels = LDT_LSMparam_struc(n)%sand%num_bins
          LDT_LSMparam_struc(n)%clay%vlevels = LDT_LSMparam_struc(n)%clay%num_bins
          LDT_LSMparam_struc(n)%silt%vlevels = LDT_LSMparam_struc(n)%silt%num_bins
          LDT_LSMparam_struc(n)%gravel%vlevels = LDT_LSMparam_struc(n)%silt%num_bins
          LDT_LSMparam_struc(n)%soilsfgrd%vlevels = LDT_LSMparam_struc(n)%sand%num_bins

          allocate(LDT_LSMparam_struc(n)%soilsfgrd%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%sand%num_bins))

          allocate(LDT_LSMparam_struc(n)%sand%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%sand%num_bins))

          allocate(LDT_LSMparam_struc(n)%clay%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%clay%num_bins))

          allocate(LDT_LSMparam_struc(n)%silt%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%silt%num_bins))

          allocate(LDT_LSMparam_struc(n)%gravel%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%gravel%num_bins))
       endif

    !- Soil color:
       if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
          LDT_soils_struc(n)%color%vlevels = LDT_soils_struc(n)%color%num_bins
          allocate(LDT_soils_struc(n)%color%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_soils_struc(n)%color%num_bins))
       endif
    !- Soil hydraulic properties:
       if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
          if( LDT_soils_struc(n)%porosity%source == "FAO" ) then
             LDT_soils_struc(n)%porosity%num_bins = 3
          else
             LDT_soils_struc(n)%porosity%num_bins = 1
          endif
          LDT_soils_struc(n)%porosity%vlevels = &
              LDT_soils_struc(n)%porosity%num_bins
          allocate(LDT_soils_struc(n)%porosity%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_soils_struc(n)%porosity%num_bins))
       endif
       if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then
          allocate(LDT_soils_struc(n)%soildepth%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_soils_struc(n)%soildepth%num_bins))
       endif

!       if(LDT_soils_struc(n)%rootdepth%selectOpt.eq.1) then
!          allocate(LDT_soils_struc(n)%rootdepth%value(&
!               LDT_rc%lnc(n),LDT_rc%lnr(n),&
!               LDT_soils_struc(n)%rootdepth%num_bins))
!       endif

       if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then
          LDT_soils_struc(n)%hsg%vlevels = LDT_soils_struc(n)%hsg%num_bins
          allocate(LDT_soils_struc(n)%hsg%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_soils_struc(n)%hsg%num_bins))
       endif

    enddo
    
!-- Read in LDT config file entries for soil parameter files:

    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo

 !- Soil texture:
    if( check_data) then 

      do n=1,LDT_rc%nnest
         if( INDEX(LDT_LSMparam_struc(n)%texture%source,"STATSGO") > 0 ) then
            LDT_rc%soil_classification(n) = "STATSGO"
         elseif( INDEX(LDT_LSMparam_struc(n)%texture%source,"CONSTANT") > 0 ) then
            LDT_rc%soil_classification(n) = "STATSGO"
         elseif( INDEX(LDT_LSMparam_struc(n)%texture%source,"ISRIC") > 0 ) then
            LDT_rc%soil_classification(n) = "STATSGO"
         else
            LDT_rc%soil_classification(n) = LDT_LSMparam_struc(n)%texture%source
         endif
      enddo
      call ESMF_ConfigFindLabel(LDT_config,"Soil texture map:", rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%txtfile(n),rc=rc)
      enddo
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soiltext_proj,&
           label="Soil texture map projection:",rc=rc)
      call LDT_verify(rc,'Soil texture map projection: option not specified in the config file')

      call ESMF_ConfigFindLabel(LDT_config,"Soil texture spatial transform:",&
           rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soiltext_gridtransform(n),&
           rc=rc)
         call LDT_verify(rc,'Soil texture spatial transform: option not specified in the config file')
      enddo

      soiltext%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, soiltext%filltype, &
           label="Soil texture fill option:",rc=rc)
      call LDT_verify(rc,"Soil texture fill option: option not specified in the config file")

      if( soiltext%filltype == "neighbor" ) then
         call ESMF_ConfigGetAttribute(LDT_config, soiltext%fillradius, &
              label="Soil texture fill radius:",rc=rc)
         call LDT_verify(rc,"Soil texture fill radius: option not specified in the config file")

         call ESMF_ConfigGetAttribute(LDT_config, soiltext%fillvalue, &
             label="Soil texture fill value:",rc=rc)
         call LDT_verify(rc,"Soil texture fill value: option not specified in the config file")

         ! EMK...Allow option to use separate fill value south of 60S.
         ! If option not specified, just use the regular fill value.
         call ESMF_ConfigFindLabel(LDT_config, &
              label="Soil texture fill value for Antarctica:", &
              rc=rc)
         if (rc == ESMF_SUCCESS) then
            call ESMF_ConfigGetAttribute(LDT_config, &
                 soiltext%fillvalue_antarctica, &
                 label="Soil texture fill value for Antarctica:", &
                 rc=rc)
         else
            soiltext%fillvalue_antarctica = soiltext%fillvalue
         end if

         ! EMK...Allow option to force exclusion of water points when
         ! filling soil texture.  Necessary when using STATSGOFAO south of
         ! 60S (texture set to all water) and multiple surface model types
         ! are used.
         ! If option not specified, set to false.

         call ESMF_ConfigFindLabel(LDT_config, &
              label="Soil texture force exclusion of water points during fill:",&
              rc=rc)
         if (rc == ESMF_SUCCESS) then
            call ESMF_ConfigGetAttribute(LDT_config, soiltext%force_exclude_water, &
                 label="Soil texture force exclusion of water points during fill:",rc=rc)
         else
            soiltext%force_exclude_water = .false.
         end if

       elseif( soiltext%filltype == "none" ) then
         write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for Soil texture"
       else
         write(LDT_logunit,*) "[ERR] Fill option for Soil texture is not valid: ",trim(soiltext%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none or neighbor "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
      end if

    endif

 !- Soil fractions:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1.or.&
            LDT_soils_struc(n)%porosity%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo

    if(check_data) then 
       call ESMF_ConfigFindLabel(LDT_config,"Sand fraction map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%safile(i),rc=rc)
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Clay fraction map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%clfile(i),rc=rc)
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Silt fraction map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%sifile(i),rc=rc)
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Gravel fraction map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gravelfile(i),rc=rc)
       enddo

       soilfrac%filltype = "none"
       call ESMF_ConfigGetAttribute(LDT_config, soilfrac%filltype, &
            label="Soils fill option:",rc=rc)
       call LDT_verify(rc,"Soils fill option: option not specified in the config file")
       if( soilfrac%filltype == "neighbor" ) then
          call ESMF_ConfigGetAttribute(LDT_config, soilfrac%fillradius, &
              label="Soils fill radius:",rc=rc)
          call LDT_verify(rc,"Soils fill radius: option not specified in the config file")

          call ESMF_ConfigGetAttribute(LDT_config, soilfrac%fillvalue, &
               label="Soils fill value:",rc=rc)
          call LDT_verify(rc,"Soils fill value: option not specified in the config file")
       elseif( soilfrac%filltype == "none" ) then
          write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for Soil fractions"
       else
         write(LDT_logunit,*) "[ERR] Fill option for Soil fractions is not valid: ",trim(soilfrac%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
       end if
    end if
    
 !- Soil color:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_soils_struc(n)%color%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo
    if(check_data) then 
       call ESMF_ConfigFindLabel(LDT_config,"Soil color map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%iscfile(i),rc=rc)
       enddo
       call ESMF_ConfigGetAttribute(LDT_config,soilclr_proj,&
            label="Soil color map projection:",rc=rc)
       call LDT_verify(rc,'Soil color map projection: option not specified in the config file')

       call ESMF_ConfigFindLabel(LDT_config,"Soil color spatial transform:",&
            rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_soils_struc(i)%soilclr_gridtransform,&
               label="Soil color spatial transform:",rc=rc)
          call LDT_verify(rc,'Soil color spatial transform: option not specified in the config file')
       enddo
    end if
    
 !- Slope type:

 !- Soil hydraulic properties:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo
    if(check_data) then 
      call ESMF_ConfigFindLabel(LDT_config,"Porosity map:", rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%pofile(i),rc=rc)
      enddo
      if( soilfrac%filltype == "neighbor" ) then
         call ESMF_ConfigGetAttribute(LDT_config, porosity%fillvalue, &
              label="Porosity fill value:",rc=rc)
         call LDT_verify(rc,"Porosity fill value: option not specified in the config file")
      endif
    endif

    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo
    if(check_data) then 
      call ESMF_ConfigFindLabel(LDT_config,"Hydrologic soil group map:", rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%hsgfile(i),rc=rc)
      enddo
    endif
   
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo
    if(check_data) then 
       call ESMF_ConfigFindLabel(LDT_config,"Soil depth map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%dsoilfile(i),rc=rc)
       enddo
    end if

    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soils_proj,&
         label="Soils map projection:",rc=rc)
    call LDT_verify(rc,'Soils map projection: option not specified in the config file')

    call ESMF_ConfigFindLabel(LDT_config,"Soils spatial transform:",rc=rc)
    do i=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soils_gridtransform(i),&
            rc=rc)
       call LDT_verify(rc,'Soils spatial transform: option not specified in the config file')
    enddo


!- Check config options for each soil-related parameter:

   do n=1, LDT_rc%nnest

    if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1 .or. &
       LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
       if( index(LDT_LSMparam_struc(n)%sand%source,"Native").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%sand%source,"STATSGOv1").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%sand%source,"CONSTANT").eq.0 ) then
         call LDT_readDomainConfigSpecs("Soils", LDT_rc%soils_proj, LDT_rc%soil_gridDesc)
         call LDT_gridOptChecks( n, "Soils", LDT_rc%soils_gridtransform(n), &
                                 LDT_rc%soils_proj, LDT_rc%soil_gridDesc(n,9) )
       endif
       call LDT_soilsOptChecks(n, "Soils", LDT_LSMparam_struc(n)%sand%source, &
                      LDT_rc%soils_gridtransform(n) )
    endif

    if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then

     ! Don't need to read in "Native" grid extents/resolution, just for LIS inputs
       if( index(LDT_LSMparam_struc(n)%texture%source,"Native").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%texture%source,"CONSTANT").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%texture%source,"ISRIC").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%texture%source,"Special").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%texture%source,"STATSGOv1").eq.0 ) then
          call LDT_readDomainConfigSpecs("Soil texture", LDT_rc%soiltext_proj, &
                                         LDT_rc%soiltext_gridDesc)
         if( LDT_rc%soiltext_proj == "latlon" ) then
           call LDT_gridOptChecks( n, "Soil texture", LDT_rc%soiltext_gridtransform(n), &
                     LDT_rc%soiltext_proj, LDT_rc%soiltext_gridDesc(n,9) )
         endif
       endif

       call LDT_soilsOptChecks(n, "Soil texture", LDT_LSMparam_struc(n)%texture%source, &
                    LDT_rc%soiltext_gridtransform(n) )
    endif

    if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
       call LDT_readDomainConfigSpecs("Soil color", &
            soilclr_proj, soilclr_gridDesc)
       call LDT_gridOptChecks( n, "Soil color", &
            LDT_soils_struc(n)%soilclr_gridtransform, &
            soilclr_proj, soilclr_gridDesc(n,9) )
       call LDT_soilsOptChecks(n, "Soil color", &
            soilclr_proj, LDT_soils_struc(n)%soilclr_gridtransform )
    endif
    
    LDT_soils_struc(n)%soilclr_proj = soilclr_proj
    LDT_soils_struc(n)%soilclr_gridDesc(:) = soilclr_gridDesc(n,:)

    if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then
       call readsoildepth(&
            trim(LDT_soils_struc(n)%soildepth%source)//char(0),&
            n,LDT_soils_struc(n)%soildepth%value,&
            LDT_LSMparam_struc(n)%landmask%value)
       
    endif

  end do
   
!- Read in soil parameter value information: 
   do n=1,LDT_rc%nnest

   != Soil texture:
      if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then

         call readsoiltexture(&
              trim(LDT_LSMparam_struc(n)%texture%source)//char(0), &
              n, LDT_LSMparam_struc(n)%texture%num_bins,           &
              LDT_LSMparam_struc(n)%texture%value,                 &
              LDT_soils_struc(n)%texture_nlyrs%value )

       ! Fill where parameter values are missing compared to land/water mask:
         if( soiltext%filltype == "neighbor" ) then 
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%texture%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%texture%short_name)
            if( INDEX(LDT_LSMparam_struc(n)%texture%source,"STATSGO") > 0 .or. &
               trim(LDT_LSMparam_struc(n)%texture%source)=="CONSTANT" .or. &
               trim(LDT_LSMparam_struc(n)%texture%source)=="Special" ) then
               soiltext%watervalue = 14.
            elseif(INDEX(LDT_LSMparam_struc(n)%texture%source,"ZOBLER") >0) then 
               soiltext%watervalue = 1.0
            elseif(INDEX(LDT_LSMparam_struc(n)%texture%source,"ISRIC") >0) then 
               soiltext%watervalue = 14.
            else
               soiltext%watervalue = LDT_rc%udef
            endif
            ! EMK...Add optional settings for Antarctica and for forcibly ignoring
            ! water points.
            call LDT_discreteParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                     LDT_rc%soiltext_gridtransform(n),                       &
                     LDT_LSMparam_struc(n)%texture%num_bins,              &
                     LDT_LSMparam_struc(n)%texture%value, soiltext%watervalue, &
                     LDT_LSMparam_struc(n)%landmask2%value,               &
                     soiltext%filltype, soiltext%fillvalue, &
                     soiltext%fillradius, &
                     soiltext%fillvalue_antarctica, &
                     soiltext%force_exclude_water)

         endif
      endif
      
   != Soil fractions (combined):
      if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1) then
         LDT_LSMparam_struc(n)%sand%short_name = "SAND"
         LDT_LSMparam_struc(n)%sand%standard_name = "Sand fraction"
         LDT_LSMparam_struc(n)%clay%short_name = "CLAY"
         LDT_LSMparam_struc(n)%clay%standard_name = "Clay fraction"
         LDT_LSMparam_struc(n)%silt%short_name = "SILT"
         LDT_LSMparam_struc(n)%silt%standard_name = "Silt fraction"
         LDT_LSMparam_struc(n)%gravel%short_name = "GRAVEL"
         LDT_LSMparam_struc(n)%gravel%standard_name = "Gravel fraction"
         
         LDT_LSMparam_struc(n)%soilsfgrd%short_name = "SOILSFGRD"
         LDT_LSMparam_struc(n)%soilsfgrd%standard_name = "Soil area fraction"
         
         call readsoilfrac(&
              trim(LDT_LSMparam_struc(n)%sand%source)//char(0),n, &
              LDT_LSMparam_struc(n)%sand%num_bins,   &
              LDT_LSMparam_struc(n)%soilsfgrd%value, &
              LDT_LSMparam_struc(n)%sand%value,      &
              LDT_LSMparam_struc(n)%clay%value,      &
              LDT_LSMparam_struc(n)%silt%value)

       ! Fill where parameter values are missing compared to land/water mask:
         if( soilfrac%filltype == "neighbor" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%sand%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%sand%short_name)
            soilfrac%watervalue = LDT_rc%udef

         !- Tiled output data fields:
            if( LDT_rc%soils_gridtransform(n) == "average" ) then
               allocate( soilfrac_array(LDT_rc%lnc(n),LDT_rc%lnr(n),3) )
               soilfrac_array(:,:,1) = LDT_LSMparam_struc(n)%sand%value(:,:,1)
               soilfrac_array(:,:,2) = LDT_LSMparam_struc(n)%silt%value(:,:,1)
               soilfrac_array(:,:,3) = LDT_LSMparam_struc(n)%clay%value(:,:,1)
               tile_temp = "tile"

               if( LDT_LSMparam_struc(n)%sand%source == "CONSTANT" ) then
                  call LDT_contTileParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
                       tile_temp, 3, soilfrac_array, LDT_LSMparam_struc(n)%soilsfgrd%value, &
                       soilfrac%watervalue, LDT_LSMparam_struc(n)%landmask2%value, &
                       soilfrac%filltype, soilfrac%fillvalue, soilfrac%fillradius )
                  LDT_LSMparam_struc(n)%sand%value(:,:,1) = soilfrac_array(:,:,1)
                  LDT_LSMparam_struc(n)%silt%value(:,:,1) = soilfrac_array(:,:,1)
                  LDT_LSMparam_struc(n)%clay%value(:,:,1) = soilfrac_array(:,:,1)

               else
                  call LDT_contTileParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                       tile_temp, 3, soilfrac_array, soilfrac_array, soilfrac%watervalue,&
                       LDT_LSMparam_struc(n)%landmask2%value, soilfrac%filltype,       &
                       soilfrac%fillvalue, soilfrac%fillradius )
                  LDT_LSMparam_struc(n)%sand%value(:,:,1) = soilfrac_array(:,:,1)
                  LDT_LSMparam_struc(n)%silt%value(:,:,1) = soilfrac_array(:,:,2)
                  LDT_LSMparam_struc(n)%clay%value(:,:,1) = soilfrac_array(:,:,3)
               endif

!               do r = 1, LDT_rc%lnr(n) 
!                  do c = 1, LDT_rc%lnc(n) 
!                     if( sum(soilfrac_array(c,r,1:3)) > 0.97 ) then
!                       LDT_LSMparam_struc(n)%soilsfgrd%value(c,r,1) = 1.0
!                     endif 
!                  enddo
!               enddo
               deallocate( soilfrac_array )
            endif
         elseif( soilfrac%filltype == "none" ) then
!            Do nothing ...
         else
           write(LDT_logunit,*)"[INFO] 'neighbor' is the only option that currently "
           write(LDT_logunit,*)"  exists for soils fraction parameter-mask agreement. "
           write(LDT_logunit,*)"  Please select: 'neighbor' for fill option ... Stopping."
           call LDT_endrun
         endif
      endif
      
   != Soil Color (used mainly with CLM, BATS, etc. models):
      if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
         call readcolor(&
              trim(LDT_soils_struc(n)%color%source)//char(0),&
              n,LDT_soils_struc(n)%color%value)

#if 0
       ! Fill where parameter values are missing compared to land/water mask:
         if( soilcolor%filltype == "neighbor" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_soils_struc(n)%color%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_soils_struc(n)%color%short_name)
            soilcolor%watervalue = LDT_rc%udef
!            soilcolor%filltype = "neighbor"
!            soilcolor%fillvalue = 4.
!            soilcolor%fillradius = 2.
            call LDT_discreteParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 LDT_soils_struc(n)%soilclr_gridtransform,         &
                 LDT_soils_struc(n)%color%num_bins,                &
                 LDT_soils_struc(n)%color%value, soilcolor%watervalue, &
                 LDT_LSMparam_struc(n)%landmask2%value,               &
                 soilcolor%filltype, soilcolor%fillvalue, soilcolor%fillradius )
         endif
#endif
      endif

      if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
         call readporosity(&
              trim(LDT_soils_struc(n)%porosity%source)//char(0),&
              n,LDT_soils_struc(n)%porosity%value,&
              LDT_LSMparam_struc(n)%landmask%value)

!         if( porosity%filltype == "average" .or. porosity%filltype == "neighbor" ) then
         if( soilfrac%filltype == "average" .or. soilfrac%filltype == "neighbor" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                               trim(LDT_soils_struc(n)%porosity%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                               trim(LDT_soils_struc(n)%porosity%short_name)
            porosity%watervalue = LDT_rc%udef
!            porosity%filltype = "average"
!            porosity%fillvalue = 0.32
!            porosity%fillradius = 3.
            call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 LDT_rc%soils_gridtransform(n),                               &
                 LDT_soils_struc(n)%porosity%num_times,                 &
                 LDT_soils_struc(n)%porosity%value, porosity%watervalue, &
                 LDT_LSMparam_struc(n)%landmask2%value, soilfrac%filltype, &
                 porosity%fillvalue, soilfrac%fillradius )
          endif
      endif

   != Read in hydrological soils group (HSG):
      if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then
         call readhsg(&
              trim(LDT_soils_struc(n)%hsg%source)//char(0),&
              n,LDT_soils_struc(n)%hsg%value,&
              LDT_LSMparam_struc(n)%landmask%value)
      endif

!      if(LDT_LSMparam_struc(n)%rootdepth%selectOpt.eq.1) then
!         call readrootdepth(&
!              trim(LDT_LSMparam_struc(n)%rootdepth%source)//char(0),&
!              n,LDT_LSMparam_struc(n)%rootdepth%value)

!       rootdepth%filltype = "average"
!       rootdepth%watervalue = LDT_rc%udef
!       rootdepth%fillvalue = 1.0
!       rootdepth%fillradius = 3.
!       call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
!                LDT_rc%soils_gridtransform(n),                           &
!                LDT_LSMparam_struc(n)%rootdepth%num_bins,             &
!                LDT_LSMparam_struc(n)%rootdepth%value, param_waterval, &
!                LDT_LSMparam_struc(n)%landmask2%value, rootdepth%filltype,&
!                rootdepth%fillvalue, rootdepth%fillradius )
!      endif

   enddo
   
 end subroutine Soils_init_LIS

!BOP
! 
! !ROUTINE: soils_init_LISHydro
! \label{soils_init_LISHydro}
! 
! !INTERFACE:
  subroutine soils_init_LISHydro(flag)

! !USES:
    use LDT_coreMod,   only : LDT_rc, LDT_config
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_paramOptCheckMod, only: LDT_soilsOptChecks, & 
                       LDT_gridOptChecks
! 
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! soils datasets and reads the soils data based on the 
! choice of options specified in the ldt configuration, 
! in the LISHydro preprocessing mode 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[readsoiltexture](\ref{readsoiltexture}) \newline
!    invokes the generic method in the registry to read
!    the soil texture data
!   \item[readsoilfrac](\ref{readsoilfrac}) \newline
!    invokes the generic method in the registry to read
!    the sand, clay and silt fraction data
!   \item[readcolor](\ref{readcolor}) \newline
!    invokes the generic method in the registry to read
!    the soil color data
!   \item[readporosity](\ref{readporosity}) \newline
!    invokes the generic method in the registry to read
!    the porosity data
!   \item[readsoildepth](\ref{readsoildepth}) \newline
!    invokes the generic method in the registry to read
!    the soil depth parameter data
!   \item[readbdrckdpth](\ref{readsoildepth}) \newline
!    invokes the generic method in the registry to read
!    the depth to bedrock parameter data
!   \item[readwpwet](\ref{readwpwet}) \newline
!    invokes the generic method in the registry to read
!    the wilting poing wetness parameter data
!   \item[readquartz](\ref{readquartz}) \newline
!    invokes the generic method in the registry to read
!    the quartz data
!   \item[readhsg](\ref{readhsg}) \newline
!    invokes the generic method in the registry to read
!    the hydrological soils group data
!  \end{description}
!
!EOP
    implicit none
    integer  :: n, i, c, r, t
    integer  :: rc
    character*50        :: soilclr_proj
    real, allocatable   :: soilfrac_array(:,:,:)
    character(50)       :: tile_temp
    real, allocatable   :: soilclr_gridDesc(:,:)
    type(LDT_fillopts)  :: soiltext
    type(LDT_fillopts)  :: soilfrac
    type(LDT_fillopts)  :: soilcolor
    type(LDT_fillopts)  :: soildepth
    type(LDT_fillopts)  :: rootdepth
    type(LDT_fillopts)  :: porosity
    logical             :: soil_select
    logical             :: check_data
    integer             :: flag

    integer             :: nc,nr, ntxtypes
    real, allocatable   :: sctdom(:,:)
    real                :: maxv, domv
    real, allocatable   :: texture1(:,:,:)    
! _____________________________________________________________________________
    

  ! Plugin:  Soil Texture Maps - Set number of texture types (attributes):
    if(LDT_LSMparam_struc(1)%texture%selectOpt.eq.1) then 
       call setTextureattribs(trim(LDT_LSMparam_struc(1)%texture%source)//char(0))  
    endif

  ! Plugin:  HSG Maps - Set number of Hydrologic Soil Group types (attributes):
    if(LDT_soils_struc(1)%hsg%selectOpt.eq.1) then
       call setHSGattribs(trim(LDT_soils_struc(1)%hsg%source)//char(0))
    endif

    soil_select = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%texture%selectOpt == 1 .or. &
           LDT_LSMparam_struc(n)%sand%selectOpt   == 1 .or. &
           LDT_soils_struc(n)%porosity%selectOpt  == 1 ) then 
         soil_select = .true.
       end if
    enddo

    if(soil_select) &
         write(LDT_logunit,*)" - - - - - - - - - - Soils Parameters - - - - - - - - - - - - - -"

    
    allocate(LDT_rc%soil_gridDesc(LDT_rc%nnest,20))
    allocate(LDT_rc%soiltext_gridDesc(LDT_rc%nnest,20))
    allocate(LDT_rc%hsg_gridDesc(LDT_rc%nnest,20))

    allocate(soilclr_gridDesc(LDT_rc%nnest,20))

    do n=1, LDT_rc%nnest

       allocate(LDT_LSMparam_struc(n)%sctdom%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               1))
       !print *, "Shape %sctdom:", shape(LDT_LSMparam_struc(n)%sctdom%value)
       LDT_rc%soil_classification(n) = "Soil texture not selected"
    !- Soil texture only:
       if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then

          LDT_LSMparam_struc(n)%texture%vlevels = &
                 LDT_LSMparam_struc(n)%texture%num_bins
          allocate(LDT_LSMparam_struc(n)%texture%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%texture%num_bins))

!          if( (index(LDT_rc%lsm,"RDHM") == 1 .or. &
!              index(LDT_rc%lsm,"SACHTET") == 1 ) .and. &
!              LDT_rc%create_soilparms_option == "create" ) then
          allocate(LDT_soils_struc(n)%texture_nlyrs%value(&
                   LDT_rc%lnc(n),LDT_rc%lnr(n),&
                   LDT_soils_struc(n)%texture_nlyrs%vlevels))
!          endif

       endif

    !- Activated soil fractions (combined for sand/silt/clay):
       if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1) then

          call ESMF_ConfigFindLabel(LDT_config,"Soil fraction number of bands:",rc=rc)
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_LSMparam_struc(n)%sand%num_bins,rc=rc)
          call LDT_verify(rc,"Soil fraction number of bands: not defined")

          LDT_LSMparam_struc(n)%clay = LDT_LSMparam_struc(n)%sand
          LDT_LSMparam_struc(n)%silt = LDT_LSMparam_struc(n)%sand
          LDT_LSMparam_struc(n)%gravel = LDT_LSMparam_struc(n)%sand
          LDT_LSMparam_struc(n)%soilsfgrd = LDT_LSMparam_struc(n)%sand

          LDT_LSMparam_struc(n)%sand%vlevels = LDT_LSMparam_struc(n)%sand%num_bins
          LDT_LSMparam_struc(n)%clay%vlevels = LDT_LSMparam_struc(n)%clay%num_bins
          LDT_LSMparam_struc(n)%silt%vlevels = LDT_LSMparam_struc(n)%silt%num_bins
          LDT_LSMparam_struc(n)%gravel%vlevels = LDT_LSMparam_struc(n)%silt%num_bins
          LDT_LSMparam_struc(n)%soilsfgrd%vlevels = LDT_LSMparam_struc(n)%sand%num_bins

          allocate(LDT_LSMparam_struc(n)%soilsfgrd%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%sand%num_bins))

          allocate(LDT_LSMparam_struc(n)%sand%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%sand%num_bins))

          allocate(LDT_LSMparam_struc(n)%clay%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%clay%num_bins))

          allocate(LDT_LSMparam_struc(n)%silt%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%silt%num_bins))

          allocate(LDT_LSMparam_struc(n)%gravel%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_LSMparam_struc(n)%gravel%num_bins))
       endif

    !- Soil color:
       if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
          LDT_soils_struc(n)%color%vlevels = LDT_soils_struc(n)%color%num_bins
          allocate(LDT_soils_struc(n)%color%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_soils_struc(n)%color%num_bins))
       endif
    !- Soil hydraulic properties:
       if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
          if( LDT_soils_struc(n)%porosity%source == "FAO" ) then
             LDT_soils_struc(n)%porosity%num_bins = 3
          else
             LDT_soils_struc(n)%porosity%num_bins = 1
          endif
          LDT_soils_struc(n)%porosity%vlevels = &
              LDT_soils_struc(n)%porosity%num_bins
          allocate(LDT_soils_struc(n)%porosity%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_soils_struc(n)%porosity%num_bins))
       endif
       if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then
          allocate(LDT_soils_struc(n)%soildepth%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_soils_struc(n)%soildepth%num_bins))
       endif

!       if(LDT_soils_struc(n)%rootdepth%selectOpt.eq.1) then
!          allocate(LDT_soils_struc(n)%rootdepth%value(&
!               LDT_rc%lnc(n),LDT_rc%lnr(n),&
!               LDT_soils_struc(n)%rootdepth%num_bins))
!       endif

       if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then
          LDT_soils_struc(n)%hsg%vlevels = LDT_soils_struc(n)%hsg%num_bins
          allocate(LDT_soils_struc(n)%hsg%value(&
               LDT_rc%lnc(n),LDT_rc%lnr(n),&
               LDT_soils_struc(n)%hsg%num_bins))
       endif

    enddo
    
!-- Read in LDT config file entries for soil parameter files:

    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo

 !- Soil texture:
    if( check_data) then 

      do n=1,LDT_rc%nnest
         if( INDEX(LDT_LSMparam_struc(n)%texture%source,"STATSGO") > 0 ) then
            LDT_rc%soil_classification(n) = "STATSGO"
         elseif( INDEX(LDT_LSMparam_struc(n)%texture%source,"CONSTANT") > 0 ) then
            LDT_rc%soil_classification(n) = "STATSGO"
         elseif( INDEX(LDT_LSMparam_struc(n)%texture%source,"ISRIC") > 0 ) then
            LDT_rc%soil_classification(n) = "STATSGO"
         else
            LDT_rc%soil_classification(n) = LDT_LSMparam_struc(n)%texture%source
         endif
      enddo
      call ESMF_ConfigFindLabel(LDT_config,"Soil texture map:", rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%txtfile(n),rc=rc)
      enddo
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soiltext_proj,&
           label="Soil texture map projection:",rc=rc)
      call LDT_verify(rc,'Soil texture map projection: option not specified in the config file')

      call ESMF_ConfigFindLabel(LDT_config,"Soil texture spatial transform:",&
           rc=rc)
      do n=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soiltext_gridtransform(n),&
           rc=rc)
         call LDT_verify(rc,'Soil texture spatial transform: option not specified in the config file')
      enddo

      soiltext%filltype = "none"
      call ESMF_ConfigGetAttribute(LDT_config, soiltext%filltype, &
           label="Soil texture fill option:",rc=rc)
      call LDT_verify(rc,"Soil texture fill option: option not specified in the config file")

      if( soiltext%filltype == "neighbor" ) then
         call ESMF_ConfigGetAttribute(LDT_config, soiltext%fillradius, &
              label="Soil texture fill radius:",rc=rc)
         call LDT_verify(rc,"Soil texture fill radius: option not specified in the config file")

         call ESMF_ConfigGetAttribute(LDT_config, soiltext%fillvalue, &
             label="Soil texture fill value:",rc=rc)
         call LDT_verify(rc,"Soil texture fill value: option not specified in the config file")

         ! EMK...Allow option to use separate fill value south of 60S.
         ! If option not specified, just use the regular fill value.
         call ESMF_ConfigFindLabel(LDT_config, &
              label="Soil texture fill value for Antarctica:", &
              rc=rc)
         if (rc == ESMF_SUCCESS) then
            call ESMF_ConfigGetAttribute(LDT_config, &
                 soiltext%fillvalue_antarctica, &
                 label="Soil texture fill value for Antarctica:", &
                 rc=rc)
         else
            soiltext%fillvalue_antarctica = soiltext%fillvalue
         end if

         ! EMK...Allow option to force exclusion of water points when
         ! filling soil texture.  Necessary when using STATSGOFAO south of
         ! 60S (texture set to all water) and multiple surface model types
         ! are used.
         ! If option not specified, set to false.

         call ESMF_ConfigFindLabel(LDT_config, &
              label="Soil texture force exclusion of water points during fill:",&
              rc=rc)
         if (rc == ESMF_SUCCESS) then
            call ESMF_ConfigGetAttribute(LDT_config, soiltext%force_exclude_water, &
                 label="Soil texture force exclusion of water points during fill:",rc=rc)
         else
            soiltext%force_exclude_water = .false.
         end if

       elseif( soiltext%filltype == "none" ) then
         write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for Soil texture"
       else
         write(LDT_logunit,*) "[ERR] Fill option for Soil texture is not valid: ",trim(soiltext%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none or neighbor "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
      end if

    endif

 !- Soil fractions:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1.or.&
            LDT_soils_struc(n)%porosity%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo

    if(check_data) then 
       call ESMF_ConfigFindLabel(LDT_config,"Sand fraction map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%safile(i),rc=rc)
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Clay fraction map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%clfile(i),rc=rc)
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Silt fraction map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%sifile(i),rc=rc)
       enddo
       call ESMF_ConfigFindLabel(LDT_config,"Gravel fraction map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%gravelfile(i),rc=rc)
       enddo

       soilfrac%filltype = "none"
       call ESMF_ConfigGetAttribute(LDT_config, soilfrac%filltype, &
            label="Soils fill option:",rc=rc)
       call LDT_verify(rc,"Soils fill option: option not specified in the config file")
       if( soilfrac%filltype == "neighbor" ) then
          call ESMF_ConfigGetAttribute(LDT_config, soilfrac%fillradius, &
              label="Soils fill radius:",rc=rc)
          call LDT_verify(rc,"Soils fill radius: option not specified in the config file")

          call ESMF_ConfigGetAttribute(LDT_config, soilfrac%fillvalue, &
               label="Soils fill value:",rc=rc)
          call LDT_verify(rc,"Soils fill value: option not specified in the config file")
       elseif( soilfrac%filltype == "none" ) then
          write(LDT_logunit,*) "[INFO] 'NONE' Parameter-Mask Agreement Option Selected for Soil fractions"
       else
         write(LDT_logunit,*) "[ERR] Fill option for Soil fractions is not valid: ",trim(soilfrac%filltype)
         write(LDT_logunit,*) "  Please select one of these:  none, neighbor or average "
         write(LDT_logunit,*) "  Programming stopping ..."
         call LDT_endrun
       end if
    end if
    
 !- Soil color:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_soils_struc(n)%color%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo
    if(check_data) then 
       call ESMF_ConfigFindLabel(LDT_config,"Soil color map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%iscfile(i),rc=rc)
       enddo
       call ESMF_ConfigGetAttribute(LDT_config,soilclr_proj,&
            label="Soil color map projection:",rc=rc)
       call LDT_verify(rc,'Soil color map projection: option not specified in the config file')

       call ESMF_ConfigFindLabel(LDT_config,"Soil color spatial transform:",&
            rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,&
               LDT_soils_struc(i)%soilclr_gridtransform,&
               label="Soil color spatial transform:",rc=rc)
          call LDT_verify(rc,'Soil color spatial transform: option not specified in the config file')
       enddo
    end if
    
 !- Slope type:

 !- Soil hydraulic properties:
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo
    if(check_data) then 
      call ESMF_ConfigFindLabel(LDT_config,"Porosity map:", rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%pofile(i),rc=rc)
      enddo
      if( soilfrac%filltype == "neighbor" ) then
         call ESMF_ConfigGetAttribute(LDT_config, porosity%fillvalue, &
              label="Porosity fill value:",rc=rc)
         call LDT_verify(rc,"Porosity fill value: option not specified in the config file")
      endif
    endif

    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo
    if(check_data) then 
      call ESMF_ConfigFindLabel(LDT_config,"Hydrologic soil group map:", rc=rc)
      do i=1,LDT_rc%nnest
         call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%hsgfile(i),rc=rc)
      enddo
    endif
   
    check_data = .false. 
    do n=1,LDT_rc%nnest
       if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then 
          check_data = .true. 
       endif
    enddo
    if(check_data) then 
       call ESMF_ConfigFindLabel(LDT_config,"Soil depth map:", rc=rc)
       do i=1,LDT_rc%nnest
          call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%dsoilfile(i),rc=rc)
       enddo
    end if

    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soils_proj,&
         label="Soils map projection:",rc=rc)
    call LDT_verify(rc,'Soils map projection: option not specified in the config file')

    call ESMF_ConfigFindLabel(LDT_config,"Soils spatial transform:",rc=rc)
    do i=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%soils_gridtransform(i),&
            rc=rc)
       call LDT_verify(rc,'Soils spatial transform: option not specified in the config file')
    enddo


!- Check config options for each soil-related parameter:

   do n=1, LDT_rc%nnest

    if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1 .or. &
       LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
       if( index(LDT_LSMparam_struc(n)%sand%source,"Native").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%sand%source,"STATSGOv1").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%sand%source,"CONSTANT").eq.0 ) then
         call LDT_readDomainConfigSpecs("Soils", LDT_rc%soils_proj, LDT_rc%soil_gridDesc)
         call LDT_gridOptChecks( n, "Soils", LDT_rc%soils_gridtransform(n), &
                                 LDT_rc%soils_proj, LDT_rc%soil_gridDesc(n,9) )
       endif
       call LDT_soilsOptChecks(n, "Soils", LDT_LSMparam_struc(n)%sand%source, &
                      LDT_rc%soils_gridtransform(n) )
    endif

    if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then

     ! Don't need to read in "Native" grid extents/resolution, just for LIS inputs
       if( index(LDT_LSMparam_struc(n)%texture%source,"Native").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%texture%source,"CONSTANT").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%texture%source,"ISRIC").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%texture%source,"Special").eq.0 .and. &
           index(LDT_LSMparam_struc(n)%texture%source,"STATSGOv1").eq.0 ) then
          call LDT_readDomainConfigSpecs("Soil texture", LDT_rc%soiltext_proj, &
                                         LDT_rc%soiltext_gridDesc)
         if( LDT_rc%soiltext_proj == "latlon" ) then
           call LDT_gridOptChecks( n, "Soil texture", LDT_rc%soiltext_gridtransform(n), &
                     LDT_rc%soiltext_proj, LDT_rc%soiltext_gridDesc(n,9) )
         endif
       endif

       call LDT_soilsOptChecks(n, "Soil texture", LDT_LSMparam_struc(n)%texture%source, &
                    LDT_rc%soiltext_gridtransform(n) )
    endif

    if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
       call LDT_readDomainConfigSpecs("Soil color", &
            soilclr_proj, soilclr_gridDesc)
       call LDT_gridOptChecks( n, "Soil color", &
            LDT_soils_struc(n)%soilclr_gridtransform, &
            soilclr_proj, soilclr_gridDesc(n,9) )
       call LDT_soilsOptChecks(n, "Soil color", &
            soilclr_proj, LDT_soils_struc(n)%soilclr_gridtransform )
    endif
    
    LDT_soils_struc(n)%soilclr_proj = soilclr_proj
    LDT_soils_struc(n)%soilclr_gridDesc(:) = soilclr_gridDesc(n,:)

    if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then
       call readsoildepth(&
            trim(LDT_soils_struc(n)%soildepth%source)//char(0),&
            n,LDT_soils_struc(n)%soildepth%value,&
            LDT_LSMparam_struc(n)%landmask%value)
       
    endif

  end do
   
!- Read in soil parameter value information: 
   do n=1,LDT_rc%nnest

   != Soil texture:
      if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then

         call readsoiltexture(&
              trim(LDT_LSMparam_struc(n)%texture%source)//char(0), &
              n, LDT_LSMparam_struc(n)%texture%num_bins,           &
              LDT_LSMparam_struc(n)%texture%value,                 &
              LDT_soils_struc(n)%texture_nlyrs%value )

       ! Fill where parameter values are missing compared to land/water mask:
         if( soiltext%filltype == "neighbor" ) then 
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%texture%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%texture%short_name)
            if( INDEX(LDT_LSMparam_struc(n)%texture%source,"STATSGO") > 0 .or. &
               trim(LDT_LSMparam_struc(n)%texture%source)=="CONSTANT" .or. &
               trim(LDT_LSMparam_struc(n)%texture%source)=="Special" ) then
               soiltext%watervalue = 14.
            elseif(INDEX(LDT_LSMparam_struc(n)%texture%source,"ZOBLER") >0) then 
               soiltext%watervalue = 1.0
            elseif(INDEX(LDT_LSMparam_struc(n)%texture%source,"ISRIC") >0) then 
               soiltext%watervalue = 14.
            else
               soiltext%watervalue = LDT_rc%udef
            endif
            ! EMK...Add optional settings for Antarctica and for forcibly ignoring
            ! water points.
            call LDT_discreteParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                     LDT_rc%soiltext_gridtransform(n),                       &
                     LDT_LSMparam_struc(n)%texture%num_bins,              &
                     LDT_LSMparam_struc(n)%texture%value, soiltext%watervalue, &
                     LDT_LSMparam_struc(n)%landmask2%value,               &
                     soiltext%filltype, soiltext%fillvalue, &
                     soiltext%fillradius, &
                     soiltext%fillvalue_antarctica, &
                     soiltext%force_exclude_water)

         endif
      endif
      
   != Soil fractions (combined):
      if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1) then
         LDT_LSMparam_struc(n)%sand%short_name = "SAND"
         LDT_LSMparam_struc(n)%sand%standard_name = "Sand fraction"
         LDT_LSMparam_struc(n)%clay%short_name = "CLAY"
         LDT_LSMparam_struc(n)%clay%standard_name = "Clay fraction"
         LDT_LSMparam_struc(n)%silt%short_name = "SILT"
         LDT_LSMparam_struc(n)%silt%standard_name = "Silt fraction"
         LDT_LSMparam_struc(n)%gravel%short_name = "GRAVEL"
         LDT_LSMparam_struc(n)%gravel%standard_name = "Gravel fraction"
         
         LDT_LSMparam_struc(n)%soilsfgrd%short_name = "SOILSFGRD"
         LDT_LSMparam_struc(n)%soilsfgrd%standard_name = "Soil area fraction"
         
         call readsoilfrac(&
              trim(LDT_LSMparam_struc(n)%sand%source)//char(0),n, &
              LDT_LSMparam_struc(n)%sand%num_bins,   &
              LDT_LSMparam_struc(n)%soilsfgrd%value, &
              LDT_LSMparam_struc(n)%sand%value,      &
              LDT_LSMparam_struc(n)%clay%value,      &
              LDT_LSMparam_struc(n)%silt%value)

       ! Fill where parameter values are missing compared to land/water mask:
         if( soilfrac%filltype == "neighbor" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%sand%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_LSMparam_struc(n)%sand%short_name)
            soilfrac%watervalue = LDT_rc%udef

         !- Tiled output data fields:
            if( LDT_rc%soils_gridtransform(n) == "average" ) then
               allocate( soilfrac_array(LDT_rc%lnc(n),LDT_rc%lnr(n),3) )
               soilfrac_array(:,:,1) = LDT_LSMparam_struc(n)%sand%value(:,:,1)
               soilfrac_array(:,:,2) = LDT_LSMparam_struc(n)%silt%value(:,:,1)
               soilfrac_array(:,:,3) = LDT_LSMparam_struc(n)%clay%value(:,:,1)
               tile_temp = "tile"

               if( LDT_LSMparam_struc(n)%sand%source == "CONSTANT" ) then
                  call LDT_contTileParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n),    &
                       tile_temp, 3, soilfrac_array, LDT_LSMparam_struc(n)%soilsfgrd%value, &
                       soilfrac%watervalue, LDT_LSMparam_struc(n)%landmask2%value, &
                       soilfrac%filltype, soilfrac%fillvalue, soilfrac%fillradius )
                  LDT_LSMparam_struc(n)%sand%value(:,:,1) = soilfrac_array(:,:,1)
                  LDT_LSMparam_struc(n)%silt%value(:,:,1) = soilfrac_array(:,:,1)
                  LDT_LSMparam_struc(n)%clay%value(:,:,1) = soilfrac_array(:,:,1)

               else
                  call LDT_contTileParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n),     &
                       tile_temp, 3, soilfrac_array, soilfrac_array, soilfrac%watervalue,&
                       LDT_LSMparam_struc(n)%landmask2%value, soilfrac%filltype,       &
                       soilfrac%fillvalue, soilfrac%fillradius )
                  LDT_LSMparam_struc(n)%sand%value(:,:,1) = soilfrac_array(:,:,1)
                  LDT_LSMparam_struc(n)%silt%value(:,:,1) = soilfrac_array(:,:,2)
                  LDT_LSMparam_struc(n)%clay%value(:,:,1) = soilfrac_array(:,:,3)
               endif

!               do r = 1, LDT_rc%lnr(n) 
!                  do c = 1, LDT_rc%lnc(n) 
!                     if( sum(soilfrac_array(c,r,1:3)) > 0.97 ) then
!                       LDT_LSMparam_struc(n)%soilsfgrd%value(c,r,1) = 1.0
!                     endif 
!                  enddo
!               enddo
               deallocate( soilfrac_array )
            endif
         elseif( soilfrac%filltype == "none" ) then
!            Do nothing ...
         else
           write(LDT_logunit,*)"[INFO] 'neighbor' is the only option that currently "
           write(LDT_logunit,*)"  exists for soils fraction parameter-mask agreement. "
           write(LDT_logunit,*)"  Please select: 'neighbor' for fill option ... Stopping."
           call LDT_endrun
         endif
      endif
      
   != Soil Color (used mainly with CLM, BATS, etc. models):
      if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
         call readcolor(&
              trim(LDT_soils_struc(n)%color%source)//char(0),&
              n,LDT_soils_struc(n)%color%value)

#if 0
       ! Fill where parameter values are missing compared to land/water mask:
         if( soilcolor%filltype == "neighbor" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_soils_struc(n)%color%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                                  trim(LDT_soils_struc(n)%color%short_name)
            soilcolor%watervalue = LDT_rc%udef
!            soilcolor%filltype = "neighbor"
!            soilcolor%fillvalue = 4.
!            soilcolor%fillradius = 2.
            call LDT_discreteParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 LDT_soils_struc(n)%soilclr_gridtransform,         &
                 LDT_soils_struc(n)%color%num_bins,                &
                 LDT_soils_struc(n)%color%value, soilcolor%watervalue, &
                 LDT_LSMparam_struc(n)%landmask2%value,               &
                 soilcolor%filltype, soilcolor%fillvalue, soilcolor%fillradius )
         endif
#endif
      endif

      if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
         call readporosity(&
              trim(LDT_soils_struc(n)%porosity%source)//char(0),&
              n,LDT_soils_struc(n)%porosity%value,&
              LDT_LSMparam_struc(n)%landmask%value)

!         if( porosity%filltype == "average" .or. porosity%filltype == "neighbor" ) then
         if( soilfrac%filltype == "average" .or. soilfrac%filltype == "neighbor" ) then
            write(LDT_logunit,*) "Checking/filling mask values for: ", &
                               trim(LDT_soils_struc(n)%porosity%short_name)
            write(fill_logunit,*) "Checking/filling mask values for: ", &
                               trim(LDT_soils_struc(n)%porosity%short_name)
            porosity%watervalue = LDT_rc%udef
!            porosity%filltype = "average"
!            porosity%fillvalue = 0.32
!            porosity%fillradius = 3.
            call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
                 LDT_rc%soils_gridtransform(n),                               &
                 LDT_soils_struc(n)%porosity%num_times,                 &
                 LDT_soils_struc(n)%porosity%value, porosity%watervalue, &
                 LDT_LSMparam_struc(n)%landmask2%value, soilfrac%filltype, &
                 porosity%fillvalue, soilfrac%fillradius )
          endif
      endif

   != Read in hydrological soils group (HSG):
      if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then
         call readhsg(&
              trim(LDT_soils_struc(n)%hsg%source)//char(0),&
              n,LDT_soils_struc(n)%hsg%value,&
              LDT_LSMparam_struc(n)%landmask%value)
      endif

!      if(LDT_LSMparam_struc(n)%rootdepth%selectOpt.eq.1) then
!         call readrootdepth(&
!              trim(LDT_LSMparam_struc(n)%rootdepth%source)//char(0),&
!              n,LDT_LSMparam_struc(n)%rootdepth%value)

!       rootdepth%filltype = "average"
!       rootdepth%watervalue = LDT_rc%udef
!       rootdepth%fillvalue = 1.0
!       rootdepth%fillradius = 3.
!       call LDT_contIndivParam_Fill( n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
!                LDT_rc%soils_gridtransform(n),                           &
!                LDT_LSMparam_struc(n)%rootdepth%num_bins,             &
!                LDT_LSMparam_struc(n)%rootdepth%value, param_waterval, &
!                LDT_LSMparam_struc(n)%landmask2%value, rootdepth%filltype,&
!                rootdepth%fillvalue, rootdepth%fillradius )
!      endif

    ! _____________________SCT_DOM_____________________________

       !print*, LDT_LSMparam_struc(n)%texture%value
       !SCT_DOM = "TEXTURE"
       ntxtypes = LDT_LSMparam_struc(n)%texture%num_bins  !size(landcover,3)
       !print *,"ntxtypes:", ntxtypes
       allocate(sctdom(LDT_rc%lnc(n),LDT_rc%lnr(n)))  !LDT_rc%lnc(n),LDT_rc%lnr(n)
       !print *, "shape(sctdom):", shape(sctdom)
       allocate(texture1(LDT_rc%lnc(n),LDT_rc%lnr(n),ntxtypes))
       texture1 = LDT_LSMparam_struc(n)%texture%value(:,:,1:16) ! LDT_LSMparam_struc(n)%soilsfgrd%short_name
       !texture1 = LDT_LSMparam_struc(n)%soilsfgrd%value
       !print*,"texture:", texture1
       nr=size(texture1,1)
       nc=size(texture1,2)
       !print *,"sctdom:", nr, nc,LDT_rc%lnr(n),LDT_rc%lnc(n)
       !ntxtypes= LDT_LSMparam_struc(n)%texture%num_bins  !size(landcover,3)
       do r = 1,LDT_rc%lnr(n)    ! nr
           do c = 1, LDT_rc%lnc(n) ! nc
              maxv = -1
              domv = -1
               do  t = 1, ntxtypes
                  if  (texture1(c,r,t).gt.maxv) then
                    maxv = texture1(c,r,t)   ! + maxv
                    domv = t
                    !domv = texture1(c,r,t)+domv
                    !print *, domv 
                 end if
               enddo
                !if (domv.ne.14) then   !non water pixels 
               sctdom(c,r) = domv
                   !print *, sctdom(c,r)
                !else
               !    sctdom(c,r) = -9999
                   !print *, sctdom(c,r)
                !endif
                !print *,c, r, sctdom(c,r)
           enddo
       enddo


    !    print*,"print sctdom", sctdom
    !    print*,"print shape(sctdom):", shape(sctdom)
        LDT_LSMparam_struc(n)%sctdom%value(:,:,1) =sctdom
        LDT_LSMparam_struc(n)%sctdom%short_name    = "SCT_DOM"
        !print *, LDT_LSMparam_struc(n)%sctdom%value(:,:,1)

    !-----------------------------------------------------   


   enddo
   
 end subroutine soils_init_LISHydro

!BOP
! !ROUTINE: soils_writeHeader_LIS
! \label{soils_writeHeader_LIS}
! 
! !INTERFACE: 
 subroutine soils_writeHeader_LIS(n,ftn,dimID)
! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
! 
! !DESCRIPTION: 
!
!  This routine writes the NetCDF headers for soils data fields
!  in the standard preprocessing mode for LIS. 
!
!EOP
    integer    :: n 
    integer    :: ftn
    integer    :: dimID(3)
    integer    :: tdimID(3)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)

 !- Soil types:
    if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then

       call LDT_verify(nf90_def_dim(ftn,'soiltypes',&
            LDT_LSMparam_struc(n)%texture%vlevels,tdimID(3)))

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_LSMparam_struc(n)%texture)

     ! Add if creating soils for SAC-HTET:
       if( (index(LDT_rc%lsm,"RDHM") == 1 .or. &
            index(LDT_rc%lsm,"SACHTET") == 1 ) .and. & 
            LDT_rc%create_soilparms_option == "create" ) then

         call LDT_verify(nf90_def_dim(ftn,'texturelayers',&
              LDT_soils_struc(n)%texture_nlyrs%vlevels,tdimID(3)))

         call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
              LDT_soils_struc(n)%texture_nlyrs)
       endif
    endif

 !- Soil fractions (combined):
    if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1) then

      call LDT_verify(nf90_def_dim(ftn,'soilfracbins',&
           LDT_LSMparam_struc(n)%sand%num_bins,tdimID(3)))

      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%soilsfgrd)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%sand)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%clay)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%silt)
!      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
!           LDT_LSMparam_struc(n)%gravel)
    endif

 !- Hydrologic soil groups:
    if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then
       call LDT_verify(nf90_def_dim(ftn,'hsgtypes',&
            LDT_soils_struc(n)%hsg%vlevels,tdimID(3)))

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_soils_struc(n)%hsg)
    endif
    if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_soils_struc(n)%color)
    endif
    if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
       call LDT_verify(nf90_def_dim(ftn,'porosity_levels',&
            LDT_soils_struc(n)%porosity%vlevels,tdimID(3)))

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_soils_struc(n)%porosity)
    endif
    if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then
       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_soils_struc(n)%soildepth)
    endif

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOILTEXT_SCHEME", &
         LDT_rc%soil_classification(1)))

  ! Attributes serving Noah-MP only (at this time):
    if ((LDT_rc%lsm.eq."Noah-MP.3.6").or.(LDT_rc%lsm.eq."Noah-MP.4.0.1")) then
    ! Number of soil types:
      if( LDT_rc%soil_classification(1) == "STATSGO" ) then
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SOILTYPES", &
              19))
      elseif( LDT_rc%soil_classification(1) == "Special" ) then
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SOILTYPES", &
              14))
      elseif( LDT_rc%soil_classification(1) == "ISRIC" ) then
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SOILTYPES", &
              19))
      else
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SOILTYPES", &
              19))
      endif

    ! Number of slope types:
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SLOPETYPES", &
              9))
    endif

#endif
  end subroutine Soils_writeHeader_LIS

!BOP
! !ROUTINE: soils_writeHeader_LISHydro
! \label{soils_writeHeader_LISHydro}
!
! !INTERFACE: 
 subroutine soils_writeHeader_LISHydro(n,ftn,dimID,flag)
! !USES: 
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
! 
! !DESCRIPTION: 
!
!  This routine writes the NetCDF headers for soils data fields
!  in the preprocessing mode for LISHydro. 
!
!EOP
    integer    :: n 
    integer    :: ftn
    integer    :: dimID(4)
    integer    :: tdimID(4)
    integer    :: flag
    integer    :: sctdomId   
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)
    tdimID(4) = dimID(4)

 !-Dominent Texture (SCT_DOM)
       ! SCT_DOM field attributes: !
    call LDT_verify(nf90_def_var(ftn, &
            "SCT_DOM", &
            nf90_float, (/tdimID(1:2),tdimID(4)/),LDT_LSMparam_struc(n)%sctdomId), &
            'nf90_def_var failed for LDT_LSMparam_struc(n)%sctdom')
    !print *,"dimID of sctdom =", (/tdimID(1:2),tdimID(3)/)


    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%sctdomId, &
            "standard_name","Dominant category"),&
            'nf90_put_att failed for sctdom:standard_name')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%sctdomId, &
            "units","-"),&
            'nf90_put_att failed for sctdom:units')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%sctdomId, &
            "scale_factor",1.0),&
            'nf90_put_att failed for sctdom:scale_factor')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%sctdomId, &
            "add_offset",0.0),&
            'nf90_put_att failed for sctdom:add_offset')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%sctdomId, &
            "missing_value",-9999.),&
            'nf90_put_att failed for sctdom:missing_value')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%sctdomId, &
            "vmin",0.0),&
            'nf90_put_att failed for sctdom:vmin')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%sctdomId, &
            "vmax",0.0),&
            'nf90_put_att failed for sctdom:vmax')
    call LDT_verify(nf90_put_att(ftn,LDT_LSMparam_struc(n)%sctdomId, &
            "stagger","M"),&
            'nf90_put_att failed for sctdom:stagger')


 !- Soil types:
    if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then

       call LDT_verify(nf90_def_dim(ftn,'soil_cat',&
            LDT_LSMparam_struc(n)%texture%vlevels,tdimID(3)))  !'soiltypes' is replaced from 'soil_cat'
 
       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_LSMparam_struc(n)%texture,flag)

     ! Add if creating soils for SAC-HTET:
       if( (index(LDT_rc%lsm,"RDHM") == 1 .or. &
            index(LDT_rc%lsm,"SACHTET") == 1 ) .and. & 
            LDT_rc%create_soilparms_option == "create" ) then

         call LDT_verify(nf90_def_dim(ftn,'texturelayers',&
              LDT_soils_struc(n)%texture_nlyrs%vlevels,tdimID(3)))

         call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
              LDT_soils_struc(n)%texture_nlyrs,flag)
       endif
    endif

 !- Soil fractions (combined):
    if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1) then

      call LDT_verify(nf90_def_dim(ftn,'soilfracbins',&
           LDT_LSMparam_struc(n)%sand%num_bins,tdimID(3)))

      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%soilsfgrd,flag)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%sand,flag)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%clay,flag)
      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
           LDT_LSMparam_struc(n)%silt,flag)
!      call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
!           LDT_LSMparam_struc(n)%gravel)
    endif

 !- Hydrologic soil groups:
    if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then
       call LDT_verify(nf90_def_dim(ftn,'hsgtypes',&
            LDT_soils_struc(n)%hsg%vlevels,tdimID(3)))

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_soils_struc(n)%hsg,flag)
    endif
    if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_soils_struc(n)%color,flag)
    endif
    if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
       call LDT_verify(nf90_def_dim(ftn,'porosity_levels',&
            LDT_soils_struc(n)%porosity%vlevels,tdimID(3)))

       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_soils_struc(n)%porosity,flag)
    endif
    if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then
       call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
            LDT_soils_struc(n)%soildepth,flag)
    endif

    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOILTEXT_SCHEME", &
         LDT_rc%soil_classification(1)))
    call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ISOILWATER", &
         14))

  ! Attributes serving Noah-MP only (at this time):
    if( LDT_rc%lsm == "Noah-MP.3.6" ) then
    ! Number of soil types:
      if( LDT_rc%soil_classification(1) == "STATSGO" ) then
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SOILTYPES", &
              19))
      elseif( LDT_rc%soil_classification(1) == "Special" ) then
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SOILTYPES", &
              14))
      elseif( LDT_rc%soil_classification(1) == "ISRIC" ) then
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SOILTYPES", &
              19))
      else
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SOILTYPES", &
              19))
      endif

    ! Number of slope types:
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"NUMBER_SLOPETYPES", &
              9))
    endif

#endif
  end subroutine soils_writeHeader_LISHydro

!BOP
! !ROUTINE: soils_writeData_LIS
! \label{soils_writeData_LIS}
! 
! !INTERFACE: 
  subroutine soils_writeData_LIS(n,ftn)
! !USES: 
    use LDT_coreMod, only : LDT_rc
! 
! !DESCRIPTION: 
!  This routine writes the soils data fields to the NetCDF file
!  in the standard preprocessing mode for LIS.  
!
!EOP
    integer      :: n 
    integer      :: ftn
    integer      :: ierr

    if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%texture)

     ! Add if creating soils for SAC-HTET:
       if( (index(LDT_rc%lsm,"RDHM") == 1 .or. &
            index(LDT_rc%lsm,"SACHTET") == 1 ) .and. &
            LDT_rc%create_soilparms_option == "create" ) then
         call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%texture_nlyrs)
       endif
    endif
    if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1) then
!    print*, 'processed surface type ', LDT_localPet
#if ( defined SPMD )
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%soilsfgrd)
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%sand)
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%clay)
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%silt)
!       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%gravel)
    endif
    if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%color)
    endif
    if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%porosity)
    endif
    if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%soildepth)
    endif
    if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%hsg)
    endif
!    if(LDT_soils_struc(n)%rootdepth%selectOpt.eq.1) then
!       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%rootdepth)
!    endif

  end subroutine Soils_writeData_LIS

!BOP
! 
! !ROUTINE: soils_writeData_LISHydro
! \label{soils_writeData_LISHydro}
!
! !INTERFACE: 
  subroutine soils_writeData_LISHydro(n,ftn,flag)
! !USES: 
    use LDT_coreMod, only : LDT_rc
    use netcdf    
! !DESCRIPTION: 
!  This routine writes the soils data fields to the NetCDF file
!  in the standard preprocessing mode for LISHydro.  
!
!EOP
    integer      :: n 
    integer      :: ftn
    integer      :: ierr
    integer      :: flag
    integer      :: nc,nr    

    !write SCT_DOM
    nr=size(LDT_LSMparam_struc(n)%sctdom%value,1)
    nc=size(LDT_LSMparam_struc(n)%sctdom%value,2)
    !print *, "nr and nc SCT_DOM writeData", nr, nc
    call LDT_verify(nf90_put_var(ftn,&
            LDT_LSMparam_struc(n)%sctdomId, LDT_LSMparam_struc(n)%sctdom%value(:,:,1),&
            (/1,1/),(/nr, nc/)),&
           'nf90_put_att failed for sctdom')

    if(LDT_LSMparam_struc(n)%texture%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%texture)

     ! Add if creating soils for SAC-HTET:
       if( (index(LDT_rc%lsm,"RDHM") == 1 .or. &
            index(LDT_rc%lsm,"SACHTET") == 1 ) .and. &
            LDT_rc%create_soilparms_option == "create" ) then
         call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%texture_nlyrs)
       endif
    endif
    if(LDT_LSMparam_struc(n)%sand%selectOpt.eq.1) then
!    print*, 'processed surface type ', LDT_localPet
#if ( defined SPMD )
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%soilsfgrd)
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%sand)
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%clay)
       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%silt)
!       call LDT_writeNETCDFdata(n,ftn,LDT_LSMparam_struc(n)%gravel)
    endif
    if(LDT_soils_struc(n)%color%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%color)
    endif
    if(LDT_soils_struc(n)%porosity%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%porosity)
    endif
    if(LDT_soils_struc(n)%soildepth%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%soildepth)
    endif
    if(LDT_soils_struc(n)%hsg%selectOpt.eq.1) then
       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%hsg)
    endif
!    if(LDT_soils_struc(n)%rootdepth%selectOpt.eq.1) then
!       call LDT_writeNETCDFdata(n,ftn,LDT_soils_struc(n)%rootdepth)
!    endif

  end subroutine soils_writeData_LISHydro
  
end module LDT_soilsMod
