!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"
module LDT_paramProcMod
!BOP
!
! !MODULE: LDT_paramProcMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manages different running
!   domain implementations
!
! !REVISION HISTORY: 
!  02 Apr 2012:  Sujay Kumar;  Initial Specification
! 
  use ESMF
  use LDT_SurfaceTypeMod
  use LDT_LMLCMod
  use LDT_gfracMod
  use LDT_albedoMod
  use LDT_soilsMod
  use LDT_laisaiMod
  use LDT_vegdataMod
  use LDT_topoMod
  use LDT_glacierMod

  use LDT_LSMparamProcMod
  use LDT_routingParamProcMod
  use LDT_lakeParamProcMod
  use LDT_openwaterMod

  use LDT_metforcingParmsMod

  use LDT_climateParmsMod
  use LDT_irrigationMod
  use LDT_LSMCropModifier_Mod

!  use LDT_metforcingParmsMod

  use LDT_logMod
  use LDT_paramDataMod
  use LDT_OPTUEMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_paramProcConfig
  public :: LDT_paramProcInit   
  public :: LDT_readProcParamInit   
  public :: LDT_paramProcWriteHeader
  public :: LDT_paramProcWrite
  public :: LDT_paramProcWriteFinalize
!EOP

contains

  
!BOP
! !ROUTINE: LDT_paramProcConfig
! \label{LDT_paramProcConfig}
!
!  This subroutine handles reading in the path/name
!   of the file to be written out containing processed
!   parameters.  Also top allocate statements for 
!   LSM parameter and forcing structures are set.
!
! !INTERFACE: 
  subroutine LDT_paramProcConfig()

! !USES:
    use LDT_coreMod, only : LDT_rc, LDT_config
    use LDT_logMod,  only : LDT_logunit, LDT_verify
    use LDT_paramMaskCheckMod

    integer   :: n 
    integer   :: rc
! ____________________________________________

    allocate(LDT_LSMparam_struc(LDT_rc%nnest))
    allocate(LDT_rc%nmaskpts(LDT_rc%nnest))  ! Used mainly for CLSMF2.5
    if( LDT_rc%nmetforc > 0 ) then           ! Metforcing parameter arrays
       allocate(LDT_force_struc(LDT_rc%nnest,LDT_rc%nmetforc_parms))
    endif

    call ESMF_ConfigFindLabel(LDT_config,&
         "Processed LSM parameter filename:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,&
            LDT_LSMparam_struc(n)%param_filename,rc=rc)
       call LDT_verify(rc,"Processed LSM parameter filename: not defined")
    enddo

  end subroutine LDT_paramProcConfig


!BOP
! !ROUTINE: LDT_paramProcInit
! \label{LDT_paramProcInit}
!
! !INTERFACE: 
  subroutine LDT_paramProcInit()

! !USES:
    use LDT_coreMod, only : LDT_rc, LDT_config
    use LDT_logMod,  only : LDT_logunit, LDT_verify
    use LDT_paramMaskCheckMod

    integer   :: n 
! ____________________________________________

    write(LDT_logunit,*) "LSM User-selected:  ",trim(LDT_rc%lsm)

    write(LDT_logunit,*) " - - - - - - MODEL PARAMETERS (NOT Selected) - - - - - - - "

 !- Read Parameter Source Options:
    call LDT_readParamSpecs()

 !- Set/read-in Parameter-Mask Fill Log File Path:
    call LDT_ParamMaskFill_Log()

!-- Initialize and/or read in each parameter selected in ldt.config:
    call LDT_surfacetype_init
    call LDT_LMLC_init
    call LDT_lakeparams_init
    call LDT_openwater_init
    call LDT_glacier_init

    call LDT_soils_init
    call LDT_greenness_init
    call LDT_albedo_init
    call LDT_laisai_init
    call LDT_topo_init
    call LDT_LSMparams_init

 !- Meteorological forcing setup:
    call LDT_forcingparms_init
    call LDT_climateParms_init

    call LDT_routingParams_init
    call LDT_irrigation_init
    call LDT_LSMCropMod_init

  end subroutine LDT_paramProcInit

!BOP
! !ROUTINE: LDT_readProcParamInit
! \label{LDT_readProcParamProcInit}
!
! !INTERFACE: 
  subroutine LDT_readProcParamInit()

! !USES:
    use LDT_coreMod, only : LDT_rc, LDT_config
    use LDT_logMod,  only : LDT_logunit, LDT_verify
    use LDT_paramMaskCheckMod
!
! !DESCRIPTION:
!   This routine looks for which data sources to be
!    read in, setting the filepath/name of the 
!    parameter-mask data fill log, and sets the main
!    non-LSM type parameters required to run LDT.
!
!EOP
    integer   :: n
    integer   :: rc
! ____________________________________________

    write(LDT_logunit,*) "[INFO] - - - - - - MODEL PARAMETERS (NOT Selected) - - - - - - - "

 !- Read Parameter Source Options:
    call LDT_readParamSpecs()

 !- Set/read-in Parameter-Mask Fill Log File Path:
    call LDT_ParamMaskFill_Log()

!-- Initialize and/or read in each parameter selected in ldt.config:
    call LDT_surfacetype_init
    call LDT_LMLC_init
    call LDT_lakeparams_init
    call LDT_openwater_init

    call LDT_topo_init

  end subroutine LDT_readProcParamInit


!BOP
! !ROUTINE: LDT_readParamSpecs
! \label{LDT_readParamSpecs}
!
! !INTERFACE:
  subroutine LDT_readParamSpecs()

! !DESCRIPTION:
!   This routine specifies the names of the data 
!    sources for each parameter type.
!
! !USES:
    use LDT_coreMod,  only : LDT_rc, LDT_config
    use LDT_domainMod,only : isSurfaceTypeSelected

    integer             :: n 
    integer             :: rc
    integer             :: m
    logical             :: const_lc
    character*50        :: const_lctype
    character*100       :: source
! _________________________________________

! - Universal LSM-based parameters -

 !- Read in Landmask Data Source Option:
    call ESMF_ConfigFindLabel(LDT_config,"Landmask data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_verify(rc,"Landmask data source: not defined")
       call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%landmask,&
            "LANDMASK",source)
       call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%dommask,&
            "DOMAINMASK",source)
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Regional mask data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_warning(rc,"Regional mask data source: not defined")
       call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%regmask,&
            "REGIONMASK",source)
    enddo

 !- Read in Landcover Data Source option: 
    call ESMF_ConfigFindLabel(LDT_config,"Landcover data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       if( isSurfaceTypeSelected(1) ) then
          call LDT_verify(rc,"Landcover data source: not defined")
          call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%landcover,&
               "LANDCOVER",source)
          call setLandcoverCategories(n,source)
        ! Note: Landcover source/classification options will be merged 
        !        and new options to be specified in future LDT versions.
       endif
    enddo
  ! CONSTANT LC Option: Read in associated landcover classification type: 
    const_lc = .false.
    do n=1,LDT_rc%nnest
      if( LDT_LSMparam_struc(n)%landcover%source == "CONSTANT" ) then
        const_lc = .true.
      endif
    end do
    if( const_lc ) then
       call ESMF_ConfigFindLabel(LDT_config,"Landcover classification:",rc=rc)
       do n=1,LDT_rc%nnest
          if( LDT_LSMparam_struc(n)%landcover%source == "CONSTANT" ) then
            call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%lc_type(n),rc=rc)
            call LDT_verify(rc,'Landcover classification: not specified')
            select case( LDT_rc%lc_type(n) )
             case( "UMD" )
               LDT_LSMparam_struc(n)%landcover%num_bins = 14
             case( "IGBPNCEP" )
               LDT_LSMparam_struc(n)%landcover%num_bins = 20
             case( "USGS" )
               LDT_LSMparam_struc(n)%landcover%num_bins = 24
             case( "MOSAIC" )
               LDT_LSMparam_struc(n)%landcover%num_bins = 7
             case( "ISA" )
               LDT_LSMparam_struc(n)%landcover%num_bins = 13
             case( "CLM45" )
               LDT_LSMparam_struc(n)%landcover%num_bins = 36
             case( "Bondville" )
               LDT_LSMparam_struc(n)%landcover%num_bins = 20
             case default
               print *, "[ERR] CONSTANT Landcover classification not recognized."
               print *, "  Options:  UMD, IGBPNCEP, USGS, MOSAIC, ISA "
               print *, " Stopping ..."
               call LDT_endrun
            end select
          else
            call ESMF_ConfigGetAttribute(LDT_config,const_lctype,rc=rc)
            call LDT_verify(rc,'Landcover classification: not specified')
          endif
       enddo
    endif
    !- Set number parameter layers based on "vertical" levels:
    !- (will be eventually removed)
    do n=1,LDT_rc%nnest
       LDT_rc%nt = LDT_LSMparam_struc(n)%landcover%num_bins
    enddo

 !- Read in Lakecover data source option:
    call ESMF_ConfigFindLabel(LDT_config,"Lakecover data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
      if( isSurfaceTypeSelected(2) ) then
         call LDT_verify(rc,"Lakecover data source: not defined")
         call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%lakecover,&
              "LAKECOVER",source)
         call setLakecoverCategories(n,source)
      endif
    enddo

    !- glacier mask data source check:
    call ESMF_ConfigFindLabel(LDT_config,"Glacier mask data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       if(isSurfaceTypeSelected(3)) then 
          call LDT_verify(rc,"Glacier mask data source: not defined")
          call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%glaciermask,&
             "GLACIERMASK",source)
          call setGlacierMaskCategories(n,source)
       endif
     enddo

 !- Read in Soil texture data source option:
    call ESMF_ConfigFindLabel(LDT_config,"Soil texture data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!       call LDT_warning(rc,"Soil texture data source: not defined")
       call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%texture,&
            "TEXTURE",source)
       if(rc.eq.0) call setTextureCategories(n,source)
    enddo
  ! LSM-required parameter check:
    if( index(LDT_rc%lsm,"Noah") == 1 .or. &
        index(LDT_rc%lsm,"RDHM") == 1 .or. &
        index(LDT_rc%lsm,"SACHTET") == 1 ) then
      if( rc /= 0 ) then
         call LDT_warning(rc,"WARNING: Soil texture data source: not defined")
         print *, "WARNING: Soil texture data source: not defined"
      endif
    endif

 !- Read in Soil Fraction Data Source Option:
    call ESMF_ConfigFindLabel(LDT_config,"Soil fraction data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%sand,&
            "SAND",source)
    enddo
  ! LSM-required parameter check:
!    if( index(LDT_rc%lsm,"Noah") == 1 .or. &
!        index(LDT_rc%lsm,"CLM" )) then
      if( rc /= 0 .and. LDT_LSMparam_struc(1)%sand%selectOpt .ne. 1) then
         call LDT_warning(rc,"WARNING: Soil fraction data source: not defined")
      endif
!    endif

 !- Read in elevation, slope, aspect Data Source Options:
    call ESMF_ConfigFindLabel(LDT_config,"Elevation data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!       call LDT_warning(rc,"Elevation data source: not defined")
       call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%elevation,&
            "ELEVATION",source)
    enddo
  ! CHECK: If Metforcing lapse-rate correction is turned on:
    if( LDT_rc%nmetforc > 0 .and. rc /= 0 ) then
      do m=1,LDT_rc%nmetforc
        if( LDT_rc%met_ecor(m) == "lapse-rate" ) then
          write(LDT_logunit,*) "Lapse-rate adjustment turned on for Metforcing: ",&
                m,", ",trim(LDT_rc%metforc(m))
          call LDT_warning(rc,"WARNING: Elevation data source: not defined")
        endif
      enddo
    endif

    call ESMF_ConfigFindLabel(LDT_config,"Slope data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_warning(rc,"Slope data source: not defined")
       call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%slope,&
            "SLOPE",source)
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Aspect data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
       call LDT_warning(rc,"Aspect data source: not defined")
       call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%aspect,&
            "ASPECT",source)
    enddo

    call ESMF_ConfigFindLabel(LDT_config,"Curvature data source:",rc=rc)
    do n=1,LDT_rc%nnest
       call ESMF_ConfigGetAttribute(LDT_config,source,rc=rc)
!       call LDT_warning(rc,"Curvature data source: not defined")
       call LDT_set_param_attribs(rc,LDT_LSMparam_struc(n)%curvature,&
            "CURVATURE",source)
    enddo

    LDT_rc%create_soilparms_option = "readin"
    if( index(LDT_rc%lsm,"RDHM") == 1 .or. &
        index(LDT_rc%lsm,"SACHTET") == 1 ) then
      call ESMF_ConfigGetAttribute(LDT_config, LDT_rc%create_soilparms_option, &
                label="Create or readin soil parameters:",rc=rc)
      call LDT_verify(rc,"Create or readin soil parameters: option not specified in the config file")
    endif
    
    call LDT_soils_readParamSpecs()
    call LDT_gfrac_readParamSpecs()
    call LDT_albedo_readParamSpecs()
    call LDT_laisai_readParamSpecs()
    call LDT_LSMCrop_readParamSpecs()
    call LDT_irrig_readParamSpecs()
    call LDT_climate_readParamSpecs()
    call LDT_vegdata_readParamSpecs()
    call LDT_glacier_readParamSpecs()

  end subroutine LDT_readParamSpecs


!BOP
! !ROUTINE: LDT_paramProcWriteHeader
! \label{LDT_paramProcWriteHeader}
!
! !DESCRIPTION:
!  This routine writes the header information
!   for the output parameter file in netcdf format.
!
! !INTERFACE: 
  subroutine LDT_paramProcWriteHeader(n)

! !USES:
    use LDT_coreMod
    use LDT_logMod
    use map_utils
    use CLM45_parmsMod
    
    integer, intent(in)   :: n 

    integer               :: iret
    integer               :: m
    integer               :: dimID(3), monthID, qID
    integer               :: bdimID(3)
    integer,allocatable   :: met_dimID(:,:)
    integer               :: tdimID
    character(len=8)      :: date
    character(len=10)     :: time
    character(len=5)      :: zone
    integer, dimension(8) :: values
    integer               :: c,r

    integer               :: shuffle, deflate, deflate_level
    integer               :: bufsize
! ________________________________________________________

    bufsize = 4

!    LDT_rc%gnc_buf(n) = LDT_rc%gnc(n) + bufsize
!    LDT_rc%gnr_buf(n) = LDT_rc%gnr(n) + bufsize

!SVK-edit
    if(LDT_masterproc) then 
       shuffle = NETCDF_shuffle
       deflate = NETCDF_deflate
       deflate_level =NETCDF_deflate_level
    
#if(defined USE_NETCDF3 || defined USE_NETCDF4)

#if(defined USE_NETCDF3) 
       iret=nf90_create(path=trim(LDT_LSMparam_struc(n)%param_filename),&
            cmode=nf90_clobber, ncid=LDT_LSMparam_struc(n)%param_file_ftn)
       call LDT_verify(iret,'creating netcdf file failed in LDT_paramProcWrite')
#endif
#if(defined USE_NETCDF4) 
       iret=nf90_create(path=trim(LDT_LSMparam_struc(n)%param_filename),&
            cmode=nf90_netcdf4, ncid=LDT_LSMparam_struc(n)%param_file_ftn)
       call LDT_verify(iret,'creating netcdf file failed in LDT_paramProcWrite')
#endif    

!-- General header:
       call date_and_time(date,time,zone,values)

!-- Write out dimensions headers:
       allocate( met_dimID(LDT_rc%nmetforc_parms,3) )
       
       !- Grid-domain dimensions:
       call LDT_verify(nf90_def_dim(LDT_LSMparam_struc(n)%param_file_ftn,'east_west',LDT_rc%gnc(n),dimID(1)))
       call LDT_verify(nf90_def_dim(LDT_LSMparam_struc(n)%param_file_ftn,'north_south',LDT_rc%gnr(n),dimID(2)))
       
       call LDT_verify(nf90_def_dim(LDT_LSMparam_struc(n)%param_file_ftn,'east_west_b',LDT_rc%gnc_buf(n),bdimID(1)))
       call LDT_verify(nf90_def_dim(LDT_LSMparam_struc(n)%param_file_ftn,'north_south_b',LDT_rc%gnr_buf(n),bdimID(2)))
    
  ! Forcing domain:
       if( LDT_rc%nmetforc > 0 ) then
          do m = 1, LDT_rc%nmetforc_parms
             call LDT_verify(nf90_def_dim(LDT_LSMparam_struc(n)%param_file_ftn,&
                  'east_west_'//trim(LDT_rc%metforc_parmsrc(m)),&
                  LDT_rc%met_nc(m),met_dimID(m,1)))
             call LDT_verify(nf90_def_dim(LDT_LSMparam_struc(n)%param_file_ftn,&
                  'north_south_'//trim(LDT_rc%metforc_parmsrc(m)),&
                  LDT_rc%met_nr(m),met_dimID(m,2)))
          end do
       end if

 !- Time-based dimensions:
       if(LDT_rc%monthlyData(n)) then 
          call LDT_verify(nf90_def_dim(LDT_LSMparam_struc(n)%param_file_ftn,'month',12,monthID))
       endif
       if(LDT_rc%quarterlyData(n)) then 
          call LDT_verify(nf90_def_dim(LDT_LSMparam_struc(n)%param_file_ftn,'quarter',4,qID))
       endif
       
       call LDT_verify(nf90_def_dim(LDT_LSMparam_struc(n)%param_file_ftn,'time',1,tdimID))
       call LDT_verify(nf90_def_var(LDT_LSMparam_struc(n)%param_file_ftn,'time',nf90_float,dimids=tdimID,&
            varID=LDT_LSMparam_struc(n)%xtimeID))
       
 !- LIS-Domain Grid file attributes:
       select case ( LDT_rc%lis_map_proj )
          
       case ( "latlon" )
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "EQUIDISTANT CYLINDRICAL"))
!Hiroko: special handling for CLM-4.5 testcase to make output identical
! to Yudong's input_lis for 0.9x1.25 global resolution; latitude has to
! be as the below specified corner values for LIS.
    if (index(LDT_rc%lsm,"CLM.4.5") == 1 .and. &
        LDT_rc%gridDesc(n,9) == 1.25     .and. &
        LDT_rc%gridDesc(n,2) == 288      .and. &
        LDT_rc%gridDesc(n,3) == 192) then
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               -89.529))
!          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
!               0.625))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"NORTH_EAST_CORNER_LAT", &
               90.393))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"NORTH_EAST_CORNER_LON", &
               LDT_rc%gridDesc(n,8)))
!          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX", &
!               1.25))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,9)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY", &
               0.942))       
    else  ! all other resolutions and domains
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
!          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"NORTH_EAST_CORNER_LAT", &
!               LDT_rc%gridDesc(n,7)))
!          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"NORTH_EAST_CORNER_LON", &
!               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,9)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,10)))       
    endif
          
       case ( "mercator" )
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "MERCATOR"))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))
          
       case ( "lambert" )    ! Lambert conformal
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "LAMBERT CONFORMAL"))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT2", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))

       case ( "polar" )    ! Polar stereographic
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "POLAR STEREOGRAPHIC"))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"ORIENT", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))

       case ( "hrap" )        ! HRAP - Grid 240 
          !    -- HRAP Grid over the Contiguous United States and Puerto Rico (polar stereographic)
          !     - Grid 240 (based on http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html)
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "HRAP"))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"ORIENT", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))

       case ( "gaussian" )
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION","GAUSSIAN"))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT",LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON",LDT_rc%gridDesc(n,5)))
          !       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT1",LDT_rc%gridDesc(n,10)))
          !       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"ORIENT",LDT_rc%gridDesc(n,7)))
          !       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"STANDARD_LON",LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX",LDT_rc%gridDesc(n,9)))
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"NUMBER OF LAT CIRCLES",LDT_rc%gridDesc(n,10)))
       end select
       
 !- Forcing-Domain Grid file attributes:
       do m = 1, LDT_rc%nmetforc_parms
          select case ( LDT_rc%met_proj(m) )
             
          case ( "latlon" )
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION_"//LDT_rc%metforc_parmsrc(m), &
                  "EQUIDISTANT CYLINDRICAL"))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,4)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,5)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"NORTH_EAST_CORNER_LAT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,7)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"NORTH_EAST_CORNER_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,8)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,9)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,10)))
             
          case ( "mercator" )
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION_"//LDT_rc%metforc_parmsrc(m), &
                  "MERCATOR"))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,4)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,5)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT1_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,10)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"STANDARD_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,11)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,8)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,9)))
             
          case ( "lambert" )    ! Lambert conformal
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION_"//LDT_rc%metforc_parmsrc(m), &
                  "LAMBERT CONFORMAL"))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,4)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,5)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT1_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,10)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT2_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,7)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"STANDARD_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,11)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,8)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,9)))

          case ( "polar" )    ! Polar stereographic
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION_"//LDT_rc%metforc_parmsrc(m), &
                  "POLAR STEREOGRAPHIC"))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,4)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,5)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT1_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,10)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"ORIENT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,7)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"STANDARD_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,11)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,8)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,9)))

          case ( "hrap" )        ! HRAP - Grid 240 
             !    -- HRAP Grid over the Contiguous United States and Puerto Rico (polar stereographic)
             !     - Grid 240 (based on http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html)
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION_"//LDT_rc%metforc_parmsrc(m), &
                  "HRAP"))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,4)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,5)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"TRUELAT1_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,10)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"ORIENT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,7)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"STANDARD_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,11)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DX_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,8)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DY_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,9)))

          case ( "gaussian" )    ! Gaussian coordinate system
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"MAP_PROJECTION_"//LDT_rc%metforc_parmsrc(m), &
                  "GAUSSIAN"))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"NORTH_WEST_CORNER_LAT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,4)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"NORTH_WEST_CORNER_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,5)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LAT_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,7)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"SOUTH_EAST_CORNER_LON_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,8)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"DI_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,9)))
             call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"N_"//LDT_rc%metforc_parmsrc(m), &
                  LDT_rc%met_gridDesc(m,10)))

          end select
       end do  ! end forcing loop

       !- Include water points:
       if( LDT_rc%inc_water_pts) then
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"INC_WATER_PTS", &
               "true"))
       else
          call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"INC_WATER_PTS", &
               "false"))
       endif

       ! - Write Parameter Header Data:
       call writeParamHeaders(n, LDT_LSMparam_struc(n)%param_file_ftn, dimID, met_dimID, monthID, qID)

       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"title", &
            "Land Data Toolkit (LDT) parameter-processed output"))
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"institution", &
            "NASA GSFC Hydrological Sciences Laboratory"))

#ifndef LDT_SKIP_HISTORY
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"history", &
            "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
            date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)))
#endif
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"references", &
            "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"))
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,NF90_GLOBAL,"comment", &
            "website: http://lis.gsfc.nasa.gov/"))
       
! - write lat, lon:
       LDT_LSMparam_struc(n)%xlat%short_name    = "lat"
       LDT_LSMparam_struc(n)%xlat%standard_name = "latitude"
       LDT_LSMparam_struc(n)%xlat%units         = "degrees_north"
       
       LDT_LSMparam_struc(n)%xlon%short_name    = "lon"
       LDT_LSMparam_struc(n)%xlon%standard_name = "longitude"
       LDT_LSMparam_struc(n)%xlon%units         = "degrees_east"
       allocate(LDT_LSMparam_struc(n)%xlat%value(LDT_rc%gnc(n),LDT_rc%gnr(n),1))
       allocate(LDT_LSMparam_struc(n)%xlon%value(LDT_rc%gnc(n),LDT_rc%gnr(n),1))
       
!Hiroko: special handling for CLM-4.5 to make output identical to Yudong's input_lis
    if( index(LDT_rc%lsm,"CLM.4.5") == 1 ) then
       do r=1,LDT_rc%gnr(n)
          do c=1,LDT_rc%gnc(n)
             LDT_LSMparam_struc(n)%xlat%value(c,r,1) = &
                           real(CLM45_struc(n)%yc%dvalue(c,r,1,1)) 
             LDT_LSMparam_struc(n)%xlon%value(c,r,1) = &
                           real(CLM45_struc(n)%xc%dvalue(c,r,1,1)) 
             !correction for CLM-4.5 lon range 0-360
             if (LDT_LSMparam_struc(n)%xlon%value(c,r,1) .ge. 180.0 ) then
              LDT_LSMparam_struc(n)%xlon%value(c,r,1) = &
                LDT_LSMparam_struc(n)%xlon%value(c,r,1) - 360.0
             endif
          enddo
       enddo
    else
       do r=1,LDT_rc%gnr(n)
          do c=1,LDT_rc%gnc(n)
             call ij_to_latlon(LDT_domain(n)%ldtglbproj,&
                  real(c), real(r), LDT_LSMparam_struc(n)%xlat%value(c,r,1),&
                  LDT_LSMparam_struc(n)%xlon%value(c,r,1))
          enddo
       enddo
    endif
       
       ! - write lat_b, lon_b:
       LDT_LSMparam_struc(n)%xlat_b%short_name    = "lat_b"
       LDT_LSMparam_struc(n)%xlat_b%standard_name = "latitude_b"
       LDT_LSMparam_struc(n)%xlat_b%units         = "degrees_north"
       
       LDT_LSMparam_struc(n)%xlon_b%short_name    = "lon_b"
       LDT_LSMparam_struc(n)%xlon_b%standard_name = "longitude_b"
       LDT_LSMparam_struc(n)%xlon_b%units         = "degrees_east"
   
!    allocate(LDT_LSMparam_struc(n)%xlat_b%value(&
!         -1:LDT_rc%gnc(n)+2,-1:LDT_rc%gnr(n)+2,1))
!    allocate(LDT_LSMparam_struc(n)%xlon_b%value(&
!         -1:LDT_rc%gnc(n)+2,-1:LDT_rc%gnr(n)+2,1))
       allocate(LDT_LSMparam_struc(n)%xlat_b%value(&
            LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n),1))
       allocate(LDT_LSMparam_struc(n)%xlon_b%value(&
            LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n),1))
       
!Hiroko: special handling for CLM-4.5 to make output identical to Yudong's input_lis
    if( index(LDT_rc%lsm,"CLM.4.5") == 1 ) then
       do r=-1,LDT_rc%gnr(n)+2
          do c=-1,LDT_rc%gnc(n)+2
             call ij_to_latlon(LDT_domain(n)%ldtproj,&
                  real(c), real(r), LDT_LSMparam_struc(n)%xlat_b%value(c+2,r+2,1),&
                  LDT_LSMparam_struc(n)%xlon_b%value(c+2,r+2,1))
             !correction for CLM-4.5 lon range 0-360
             if (LDT_LSMparam_struc(n)%xlon_b%value(c+2,r+2,1) .ge. 180.0) then
              LDT_LSMparam_struc(n)%xlon_b%value(c+2,r+2,1) = &
                LDT_LSMparam_struc(n)%xlon_b%value(c+2,r+2,1) - 360.0
             endif
          enddo
       enddo
    else
       do r=-1,LDT_rc%gnr(n)+2
          do c=-1,LDT_rc%gnc(n)+2
             call ij_to_latlon(LDT_domain(n)%ldtproj,&
                  real(c), real(r), LDT_LSMparam_struc(n)%xlat_b%value(c+2,r+2,1),&
                  LDT_LSMparam_struc(n)%xlon_b%value(c+2,r+2,1))
          enddo
       enddo
    endif
       
!! Xlat field attributes: !!

#if(defined USE_NETCDF3 || defined USE_NETCDF4)    
       call LDT_verify(nf90_def_var(LDT_LSMparam_struc(n)%param_file_ftn,&
            trim(LDT_LSMparam_struc(n)%xlat%short_name),&
            nf90_float, dimids=dimID(1:2),varid=LDT_LSMparam_struc(n)%xlatid),&
            'nf90_def_var failed for LDT_LSMparam_struc(n)%xlat')

#if(defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_LSMparam_struc(n)%param_file_ftn,&
            LDT_LSMparam_struc(n)%xlatid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for LDT_LSMparam_struc(n)%xlat')
#endif
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatid, &
            "standard_name",trim(LDT_LSMparam_struc(n)%xlat%standard_name)),&
            'nf90_put_att failed for xlat:standard_name')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatid, &
            "units",trim(LDT_LSMparam_struc(n)%xlat%units)),&
            'nf90_put_att failed for xlat:units')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatid, &
            "scale_factor",1.0),&
            'nf90_put_att failed for xlat:scale_factor')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatid, &
            "add_offset",0.0),&
            'nf90_put_att failed for xlat:add_offset')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatid, &
            "missing_value",LDT_rc%udef),&
            'nf90_put_att failed for xlat:missing_value')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatid, &
            "vmin",0.0),&
            'nf90_put_att failed for xlat:vmin')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatid, &
            "vmax",0.0),&
            'nf90_put_att failed for xlat:vmax')

       !! Xlon field attributes: !!
       call LDT_verify(nf90_def_var(LDT_LSMparam_struc(n)%param_file_ftn,&
            trim(LDT_LSMparam_struc(n)%xlon%short_name),&
            nf90_float, dimids=dimID(1:2),varid=LDT_LSMparam_struc(n)%xlonid),&
            'nf90_def_var failed for LDT_LSMparam_struc(n)%xlon')

#if(defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_LSMparam_struc(n)%param_file_ftn,&
            LDT_LSMparam_struc(n)%xlonid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for LDT_LSMparam_struc(n)%xlon')
#endif
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonid, &
            "standard_name",trim(LDT_LSMparam_struc(n)%xlon%standard_name)),&
            'nf90_put_att failed for xlon:standard_name')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonid, &
            "units",trim(LDT_LSMparam_struc(n)%xlon%units)),&
            'nf90_put_att failed for xlon:units')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonid, &
            "scale_factor",1.0),&
            'nf90_put_att failed for xlon:scale_factor')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonid, &
            "add_offset",0.0),&
            'nf90_put_att failed for xlon:add_offset')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonid, &
            "missing_value",LDT_rc%udef),&
            'nf90_put_att failed for xlon:missing_value')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonid, &
            "vmin",0.0),&
            'nf90_put_att failed for xlon:vmin')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonid, &
            "vmax",0.0),&
            'nf90_put_att failed for xlon:vmax')

       !! Xlat_b field attributes: !!
       call LDT_verify(nf90_def_var(LDT_LSMparam_struc(n)%param_file_ftn,&
            trim(LDT_LSMparam_struc(n)%xlat_b%short_name),&
            nf90_float, dimids=bdimID(1:2),varid=LDT_LSMparam_struc(n)%xlatbid),&
            'nf90_def_var failed for LDT_LSMparam_struc(n)%xlat_b')

#if(defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_LSMparam_struc(n)%param_file_ftn,&
            LDT_LSMparam_struc(n)%xlatbid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for LDT_LSMparam_struc(n)%xlat_b')
#endif

       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatbid, &
            "standard_name",trim(LDT_LSMparam_struc(n)%xlat_b%standard_name)),&
            'nf90_put_att failed for xlat_b:standard_name')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatbid, &
            "units",trim(LDT_LSMparam_struc(n)%xlat_b%units)),&
            'nf90_put_att failed for xlat_b:units')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatbid, &
            "scale_factor",1.0),&
            'nf90_put_att failed for xlat_b:scale_factor')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatbid, &
            "add_offset",0.0),&
            'nf90_put_att failed for xlat_b:add_offset')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatbid, &
            "missing_value",LDT_rc%udef),&
            'nf90_put_att failed for xlat_b:missing_value')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatbid, &
            "vmin",0.0),&
            'nf90_put_att failed for xlat_b:vmin')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlatbid, &
            "vmax",0.0),&
            'nf90_put_att failed for xlat_b:vmax')

       !! Xlon_b field attributes: !!
       call LDT_verify(nf90_def_var(LDT_LSMparam_struc(n)%param_file_ftn,&
            trim(LDT_LSMparam_struc(n)%xlon_b%short_name),&
            nf90_float, dimids=bdimID(1:2),varid=LDT_LSMparam_struc(n)%xlonbid),&
            'nf90_def_var failed for LDT_LSMparam_struc(n)%xlon_b')

#if(defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_LSMparam_struc(n)%param_file_ftn,&
            LDT_LSMparam_struc(n)%xlonbid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for LDT_LSMparam_struc(n)%xlon_b')
#endif
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonbid, &
            "standard_name",trim(LDT_LSMparam_struc(n)%xlon_b%standard_name)),&
            'nf90_put_att failed for xlon_b:standard_name')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonbid, &
            "units",trim(LDT_LSMparam_struc(n)%xlon_b%units)),&
            'nf90_put_att failed for xlon_b:units')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonbid, &
            "scale_factor",1.0),&
            'nf90_put_att failed for xlon_b:scale_factor')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonbid, &
            "add_offset",0.0),&
            'nf90_put_att failed for xlon_b:add_offset')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonbid, &
            "missing_value",LDT_rc%udef),&
            'nf90_put_att failed for xlon_b:missing_value')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonbid, &
            "vmin",0.0),&
            'nf90_put_att failed for xlon_b:vmin')
       call LDT_verify(nf90_put_att(LDT_LSMparam_struc(n)%param_file_ftn,LDT_LSMparam_struc(n)%xlonbid, &
            "vmax",0.0),&
            'nf90_put_att failed for xlon_b:vmax')

#endif
!End netcdf USE directives
    endif
  end subroutine LDT_paramProcWriteHeader

!BOP
! !ROUTINE: LDT_paramProcWrite
! \label{LDT_paramProcWrite}
!
! !INTERFACE: 
  subroutine LDT_paramProcWrite(n)

! !USES:
    use LDT_coreMod
    use LDT_logMod
    use map_utils
    
    integer, intent(in)   :: n 

    integer               :: m
    integer               :: dimID(3), monthID, qID
    integer,allocatable   :: met_dimID(:,:)
    integer               :: tdimID, xtimeID

!SVK-edit    
    if(LDT_masterproc) then 
       call LDT_verify(nf90_enddef(LDT_LSMparam_struc(n)%param_file_ftn))
       
       call LDT_verify(nf90_put_var(LDT_LSMparam_struc(n)%param_file_ftn,&
            LDT_LSMparam_struc(n)%xlatid, LDT_LSMparam_struc(n)%xlat%value(:,:,1),&
            (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
            'nf90_put_att failed for xlat')
       
       call LDT_verify(nf90_put_var(LDT_LSMparam_struc(n)%param_file_ftn,&
            LDT_LSMparam_struc(n)%xlonid, LDT_LSMparam_struc(n)%xlon%value(:,:,1),&
            (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
            'nf90_put_att failed for xlon')

       call LDT_verify(nf90_put_var(LDT_LSMparam_struc(n)%param_file_ftn,&
            LDT_LSMparam_struc(n)%xlatbid, LDT_LSMparam_struc(n)%xlat_b%value(:,:,1),&
            (/1,1/),(/LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n)/)),&
            'nf90_put_att failed for xlat_b')
       
       call LDT_verify(nf90_put_var(LDT_LSMparam_struc(n)%param_file_ftn,&
            LDT_LSMparam_struc(n)%xlonbid, LDT_LSMparam_struc(n)%xlon_b%value(:,:,1),&
            (/1,1/),(/LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n)/)),&
            'nf90_put_att failed for xlon_b')
       
       call LDT_verify(nf90_put_var(LDT_LSMparam_struc(n)%param_file_ftn,&
            LDT_LSMparam_struc(n)%xtimeID,0.0))

    endif

! - Write Parameter Output Data:
    call writeParamData(n,LDT_LSMparam_struc(n)%param_file_ftn)

#endif    

  end subroutine LDT_paramProcWrite

!BOP
! !ROUTINE: LDT_paramProcWriteFinalize
! \label{LDT_paramProcWriteFinalize}
!
! !INTERFACE: 
  subroutine LDT_paramProcWriteFinalize(n)

! !USES:
    use LDT_coreMod
    use LDT_logMod
    use map_utils
#if ( defined SPMD )
    use mpi
#endif
! !DESCRIPTION:
!
!EOP
    
    integer, intent(in)   :: n 
    integer               :: ierr
! ________________________________________________________

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    if(LDT_masterproc) then 
! - Close file:
       call LDT_verify(nf90_close(LDT_LSMparam_struc(n)%param_file_ftn))
    endif
#endif    
!    print*, 'reached finalize ', LDT_localPet
#if ( defined SPMD )
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

    write(unit=LDT_logunit,fmt=*) "-- Finished writing parameters to netcdf output file -- "

  end subroutine LDT_paramProcWriteFinalize


  subroutine writeParamHeaders(n, ftn, dimID, met_dimID, monthID, qID)
    
! !USES:
    use LDT_coreMod, only : LDT_rc, LDT_config
    integer     :: n 
    integer     :: ftn
    integer     :: dimID(3)
    integer     :: met_dimID(LDT_rc%nmetforc_parms,3)
    integer     :: monthID
    integer     :: qID
    
    call LDT_LMLC_writeHeader(n,ftn,dimID)
    call LDT_surfacetype_writeHeader(n,ftn,dimID)
    call LDT_LSMCropMod_writeHeader(n,ftn,dimID)
    call LDT_lakeparams_writeHeader(n,ftn,dimID,monthID)

    call LDT_soils_writeHeader(n,ftn,dimID)
    call LDT_topo_writeHeader(n,ftn,dimID)
    call LDT_greenness_writeHeader(n,ftn,dimID,monthID)
    call LDT_albedo_writeHeader(n,ftn,dimID,monthID,qID)
    call LDT_laisai_writeHeader(n,ftn,dimID,monthID)

    call LDT_LSMparams_writeHeader(n,ftn,dimID,monthID)
    call LDT_routingparams_writeHeader(n,ftn,dimID,monthID)
    call LDT_irrigation_writeHeader(n,ftn,dimID)
    call LDT_climateParms_writeHeader(n,ftn,dimID,met_dimID,monthID)

! - Forcing-specific parameter headers
    call LDT_forcingParms_writeHeader(n,ftn,dimID,met_dimID)

!OPT/UE parameters
    if (LDT_rc%runmode.eq."OPTUE parameter processing") then
       call LDT_optue_writeHeader(n,ftn,dimID)
    endif

  end subroutine writeParamHeaders


  subroutine writeParamData(n, ftn)

! !USES:
    use LDT_coreMod

    integer  :: n 
    integer  :: ftn

    call LDT_LMLC_writeData(n,ftn)
    call LDT_surfacetype_writeData(n,ftn)

    call LDT_LSMCropMod_writeData(n,ftn)
    call LDT_lakeparams_writeData(n,ftn)

    call LDT_soils_writeData(n,ftn)
    call LDT_topo_writeData(n,ftn)
    call LDT_greenness_writeData(n,ftn)
    call LDT_albedo_writeData(n,ftn)
    call LDT_laisai_writeData(n,ftn)

    call LDT_LSMparams_writeData(n,ftn)
    call LDT_routingparams_writeData(n,ftn)
    call LDT_irrigation_writeData(n,ftn)
    call LDT_climateParms_writeData(n,ftn)

! - Forcing-specific data
    call LDT_forcingParms_writeData(n,ftn)

! OPT/UE parameters
    if (LDT_rc%runmode.eq."OPTUE parameter processing") then
       call LDT_optue_writeData(n,ftn)
    endif

  end subroutine writeParamData

  subroutine setLakecoverCategories(n,source)
    integer          :: n 
    character(len=*) :: source

    select case( source )
     case ( "GLDBv1", "GLDBv2" ) 
       LDT_LSMparam_struc(n)%lakecover%num_bins = 1
       LDT_LSMparam_struc(n)%lakecover%standard_name = &
            "(FLake) GLDB lake fraction map"
     case default
       print*, "The Lake Cover Type, ",trim(source),", is not known. "
       print*, "Please define setLakecoverCategories" 
       print*, "Stopping ..."
       stop
    end select

  end subroutine setLakecoverCategories

  subroutine setGlacierMaskCategories(n,source)
    integer          :: n 
    character(len=*) :: source

    select case( source )
     case ( "GLIMS", "MODIS_Native" ) 
       LDT_LSMparam_struc(n)%glaciermask%num_bins = 1
       LDT_LSMparam_struc(n)%glaciermask%standard_name = &
            "GLIMS glacier mask map"
     case default
       print*, "The glaciermask source ",trim(source),", is not known. "
       print*, "Please define setGlaciermaskCategories" 
       print*, "Stopping ..."
       stop
    end select

  end subroutine setGlaciermaskCategories

  subroutine setTextureCategories(n,source)
    integer          :: n 
    character(len=*) :: source

    if(source.ne."none") then 
       if(source.eq."STATSGOFAO_Native") then 
          LDT_LSMparam_struc(n)%texture%num_bins = 16
          LDT_LSMparam_struc(n)%texture%standard_name = &
            "(NCAR) STATSGO+FAO blended soil texture map"

       elseif(source.eq."STATSGOFAO_LIS") then 
          LDT_LSMparam_struc(n)%texture%num_bins = 16
          LDT_LSMparam_struc(n)%texture%standard_name = &
            "(LIS-modified) STATSGO+FAO blended soil texture map"

       elseif(source.eq."STATSGOv1") then 
          LDT_LSMparam_struc(n)%texture%num_bins = 16
          LDT_LSMparam_struc(n)%texture%standard_name = &
            "(USDA/PSU) STATSGOv1 CONUS-only soil texture map"

       elseif(source.eq."ZOBLER_GFS") then 
          LDT_LSMparam_struc(n)%texture%num_bins = 9
          LDT_LSMparam_struc(n)%texture%standard_name = &
            "ZOBLER GFS soil texture map"

       elseif(source.eq."CONSTANT") then
          LDT_LSMparam_struc(n)%texture%num_bins = 16
          LDT_LSMparam_struc(n)%texture%standard_name = &
            "CONSTANT soil texture map"

       elseif(source.eq."ISRIC") then
          LDT_LSMparam_struc(n)%texture%num_bins = 13
          LDT_LSMparam_struc(n)%texture%standard_name = &
            "ISRIC soil texture map"

       elseif(source.eq."Special") then  
          LDT_LSMparam_struc(n)%texture%num_bins = 14
          LDT_LSMparam_struc(n)%texture%standard_name = &
            "Special soil texture map"

       else
          print*, 'Please define setTextureCategories for '//trim(source)
          stop
       endif
    endif
  end subroutine setTextureCategories

  subroutine setLandcoverCategories(n,source)

    use LDT_coreMod, only : LDT_rc
    integer          :: n 
    character(len=*) :: source

 ! Note: landcover "source" options will be updated in future

    LDT_LSMparam_struc(:)%landcover%vlevels = 1

    select case( source )

      case( "AVHRR", "AVHRR_GFS", "SACHTET.3.5.6", "RDHM.3.5.6" )
        LDT_rc%lc_type(n) = "UMD"
        LDT_LSMparam_struc(n)%landcover%num_bins = 14
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "AVHRR UMD landcover map"

      case( "MODIS_Native", "MODIS_LIS" )
        LDT_rc%lc_type(n) = "IGBPNCEP"
        LDT_LSMparam_struc(n)%landcover%num_bins = 20
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "MODIS-IGBP (NCEP-modified) landcover map"

      case( "MODIS_Native_PFT" )
        LDT_rc%lc_type(n) = "JULES_PFT"
        LDT_LSMparam_struc(n)%landcover%num_bins = 9
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "JULES PFT landcover map based on MODIS-IGBP (NCEP-modified)"
      
      case( "UKMO_IGBP_Native_PFT" )
        LDT_rc%lc_type(n) = "JULES_PFT"
        LDT_LSMparam_struc(n)%landcover%num_bins = 9
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "JULES PFT landcover map based on UKMO IGBP "

      case( "UM_Native_Ancillary" )
        LDT_rc%lc_type(n) = "JULES_PFT"
        LDT_LSMparam_struc(n)%landcover%num_bins = 10
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "UM CAP generated land cover data "


      case( "USGS_Native", "USGS_LIS" ) 
        LDT_rc%lc_type(n) = "USGS"
        LDT_LSMparam_struc(n)%landcover%num_bins = 24
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "USGS landcover map"

      case( "MOSAIC", "CLSMF2.5" ) 
        LDT_rc%lc_type(n) = "MOSAIC"
        LDT_LSMparam_struc(n)%landcover%num_bins = 7
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "CLSMF2.5/MOSAIC landcover map"

      case( "ISA" ) 
        LDT_rc%lc_type(n) = "ISA"
        LDT_LSMparam_struc(n)%landcover%num_bins = 13
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "ISA landcover map"

      case( "ALMIPII" ) 
        LDT_rc%lc_type(n) = "ECOCLIMAP2"
        LDT_LSMparam_struc(n)%landcover%num_bins = 12
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "ECOCLIMAP2 landcover map"

      case( "VIC411", "VIC412" ) 
        LDT_rc%lc_type(n) = "UMD"
        LDT_LSMparam_struc(n)%landcover%num_bins = 14
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "VIC-based AVHRR/UMD landcover map"

      case( "CLM45" ) 
        LDT_rc%lc_type(n) = "CLM45"
        LDT_LSMparam_struc(n)%landcover%num_bins = 36
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "CLM-4.5 PFT and landunits landcover map"

      case( "CONSTANT" ) 
        LDT_LSMparam_struc(n)%landcover%num_bins = 13
        LDT_LSMparam_struc(n)%landcover%standard_name = &
            "'CONSTANT' landcover field"
     case default
       print*, 'Please define setLandcoverCategories '
       print*, ' Stopping ... '
       stop
    end select


  end subroutine setLandcoverCategories

end module LDT_paramProcMod

