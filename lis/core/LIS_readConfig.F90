!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LIS_misc.h"

!BOP
!
! !ROUTINE: LIS_readConfig
!  \label{LIS_readConfig}
!
! !REVISION HISTORY:
!
!  23 Feb 2001: Urszula Jambor; Added GEOS or GDAS forcing option
!  27 Mar 2001: Jon Gottschalck; Revision of subroutine by implementing namelists
!  15 Apr 2002: Urszula Jambor; Added ECMWF forcing options, also
!               adding 1 & 1/2 degree GLDAS domain options.
!  28 Apr 2002: Kristi Arsenault; Added NOAH LSM code
!  14 Nov 2003: Sujay Kumar; Modified card file that includes regional
!               modeling options
!  02 Feb 2006: Sujay Kumar; Switched to the Inpak format
!  19 Jan 2007: Chuck Alonge; Added read for Snow Depth (snowh)
!               and added parameter output option (wparm)
!  29 Dec 2007: Marv Freimund; Used trim on filenames
!  17 Jan 2011: David Mocko, added max/min greenness & slope type
!  15 May 2023  Sujay Kumar, added support for 1D lat/lon output for latlon and merc projections
!
! !INTERFACE:
subroutine LIS_readConfig()
! !USES:
  use ESMF 
  use LIS_coreMod,     only : LIS_rc, LIS_config, &
                              LIS_localPet, LIS_npes, LIS_masterproc
  use LIS_histDataMod, only : LIS_histData
  use LIS_timeMgrMod,  only : LIS_date2time, LIS_parseTimeString
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_logMod
  use LIS_mpiMod, only: LIS_mpi_comm ! EMK
!
! !DESCRIPTION:
!
!  The code in this file initializes the LIS configuration management utility.
!  The generic, model-independent runtime specifications of a LIS simulation
!  are read by this routine. The routine also initializes the LIS log buffers
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     convert the date format to a floating point time format
!  \end{description}
!EOP
  implicit none

  character*100  :: temp1
  integer        :: i,rc
  character*1    :: fproc(4)
  integer        :: doy
  real           :: gmt
  integer        :: ninsts
  integer        :: npert_forc, npert_state, npert_obs
  logical        :: exists
  character*10   :: time
  integer        :: ios
  integer        :: final_dirpos
  character(len=LIS_CONST_PATH_LEN) :: diag_fname
  character(len=LIS_CONST_PATH_LEN) :: diag_dir
  integer :: ierr ! EMK
  integer, external  :: LIS_create_subdirs
! ______________________________________________________________

  if ( LIS_masterproc ) then

     inquire(file=trim(LIS_rc%lis_config_file), exist=exists)
     if( .not. exists ) then    
        write(*,*) "[ERR] LIS config file, ",     &
                          trim(LIS_rc%lis_config_file), &
                          ", does not exist."
        write(*,*) " Also, if the LIS config file you wanted to run with" 
        write(*,*) "  is not 'lis.config', please put an '-f' or '--file'"
        write(*,*) "  in front of your intended LIS config file."
        write(*,*) " This LIS run is stopping here ..."
        call LIS_endrun
     endif
  endif
  LIS_config = ESMF_ConfigCreate(rc=rc)
  call ESMF_ConfigLoadFile(LIS_config,LIS_rc%lis_config_file,rc=rc)

!------------------------------------------------------------------------
! Open runtime diagnostics file
!------------------------------------------------------------------------
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%dfile,label="Diagnostic output file:",&
       rc=rc)
  call LIS_verify(rc,'Diagnostic output file: not defined')

 ! Make the diagnostic file directory names/path:
  diag_fname = LIS_rc%dfile
  final_dirpos = scan(diag_fname, "/", BACK = .TRUE.)
  if(final_dirpos.ne.0) then 
    diag_dir = diag_fname(1:final_dirpos)
    !EMK...Only a single invocation is needed.
    !ios = LIS_create_subdirs(len_trim(diag_dir),trim(diag_dir))
    if (LIS_masterproc) then
       ios = LIS_create_subdirs(len_trim(diag_dir),trim(diag_dir))
       if (ios .ne. 0) then
          write(LIS_logunit,*)'[ERR] Problem creating directory ', &
               trim(diag_dir)
          flush(LIS_logunit)
       end if
    end if
  endif

! EMK... Make sure diagnostic file directory has been created before 
! continuing.
#if (defined SPMD)
  call mpi_barrier(LIS_mpi_comm, ierr)
#endif

  write(unit=temp1,fmt='(i4.4)') LIS_localPet
  read(unit=temp1,fmt='(4a1)')fproc
  LIS_rc%dfile = trim(LIS_rc%dfile)//"."//fproc(1)//fproc(2)//fproc(3)//fproc(4)
  open(unit=LIS_logunit,file=trim(LIS_rc%dfile))


!------------------------------------------------------------------------
! Reading in parameters that need to be initialized
! to avoid any problems later.
!------------------------------------------------------------------------

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%runmode,&
       label="Running mode:",rc=rc)
  call LIS_verify(rc,'Running mode: option not specified in the config file')

!  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%lis_map_proj,&
!       label="Map projection of the LIS domain:",rc=rc)
!  call LIS_verify(rc,'Map projection of the LIS domain: option not specified in the config file')

  ! CM Grabs new optional lis.config entry for the number of dimensions of the lat/lon fields
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nlatlon_dimensions,&
       label="Number of dimensions in the lat/lon output fields:",rc=rc)
  call LIS_warning(rc, 'Number of dimensions in the lat/lon output fields: option not specified in the config file. Assigning value to "1D"')

    ! CM If the user did not specify the number of dimension for the lat/lon fields, use 1D. In LIS_domainMod, this will switch to 2D for all projections except latlon. 
    if ( rc /= 0 ) then
      LIS_rc%nlatlon_dimensions = "1D"
    endif   
  
    ! CM If the user specified an invalid dimension option, assign to "1D"
    if ( LIS_rc%nlatlon_dimensions /= "1D" .AND. LIS_rc%nlatlon_dimensions /= "2D" ) then
      call LIS_warning(1,'Invalid lis.config entry, "Number of dimensions in the lat/lon output fields:" Assigning value to "1D"')
      LIS_rc%nlatlon_dimensions = "1D"
    endif

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nnest,&
       label="Number of nests:",rc=rc)
  call LIS_verify(rc,'Number of nests: option not specified in the config file')
  
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nsf_model_types,&
       label="Number of surface model types:",rc=rc)
  call LIS_verify(rc,'Number of surface model types: option not specified in the config file')

  LIS_rc%max_model_types = 5

  LIS_rc%lsm_index       = 1
  LIS_rc%lake_index      = 2
  LIS_rc%glacier_index   = 3
  LIS_rc%wetland_index   = 4
  LIS_rc%openwater_index = 5

  LIS_rc%glaciermodel = "none"
  allocate(LIS_rc%iterationId(LIS_rc%nnest))
  LIS_rc%iterationId = 1
   
  allocate(LIS_rc%sf_model_type_name_select(LIS_rc%nsf_model_types))
  allocate(LIS_rc%sf_model_type_select(LIS_rc%nsf_model_types))

  allocate(LIS_rc%sf_model_type_name(LIS_rc%max_model_types))
  allocate(LIS_rc%sf_model_type(LIS_rc%max_model_types))
  allocate(LIS_rc%npatch(LIS_rc%nnest, LIS_rc%max_model_types))

  LIS_rc%sf_model_type_name(1) = "LSM"
  LIS_rc%sf_model_type(1)      = LIS_rc%lsm_index
  LIS_rc%sf_model_type_name(2) = "Lake"
  LIS_rc%sf_model_type(2)      = LIS_rc%lake_index
  LIS_rc%sf_model_type_name(3) = "Glacier"
  LIS_rc%sf_model_type(3)      = LIS_rc%glacier_index
  LIS_rc%sf_model_type_name(4) = "Wetland"
  LIS_rc%sf_model_type(4)      = LIS_rc%wetland_index
  LIS_rc%sf_model_type_name(5) = "Openwater"
  LIS_rc%sf_model_type(5)      = LIS_rc%openwater_index
  
  call ESMF_ConfigFindLabel(LIS_config,"Surface model types:",rc=rc)
  do i=1,LIS_rc%nsf_model_types
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%sf_model_type_name_select(i),rc=rc)
     call LIS_verify(rc,"Surface model types: not defined")
     
     call LIS_mapSurfaceModelType(LIS_rc%sf_model_type_name_select(i), &
          LIS_rc%sf_model_type_select(i))
  enddo

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nmetforc,&
       label="Number of met forcing sources:",rc=rc)
  call LIS_verify(rc,'Number of met forcing sources: not specified')

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%metforc_blend_alg,&
       label="Blending method for forcings:",rc=rc)
  call LIS_verify(rc,'Blending method for forcings: not specified')

  allocate(LIS_rc%metforc(LIS_rc%nmetforc))
  allocate(LIS_rc%met_nf(LIS_rc%nmetforc))
  allocate(LIS_rc%met_nensem(LIS_rc%nmetforc))
  allocate(LIS_rc%met_nperforc(LIS_rc%nmetforc))
  allocate(LIS_rc%met_ecor(LIS_rc%nmetforc))
  allocate(LIS_rc%pcp_downscale(LIS_rc%nmetforc))
  allocate(LIS_rc%met_interp(LIS_rc%nmetforc))
  allocate(LIS_rc%met_upscale(LIS_rc%nmetforc))
  allocate(LIS_rc%met_tinterp(LIS_rc%nmetforc))
  allocate(LIS_rc%met_proj(LIS_rc%nmetforc))
!  allocate(LIS_rc%metforc_ensmem(LIS_rc%nmetforc))

! Default initialization
  LIS_rc%met_nensem = 1
  LIS_rc%met_nperforc = 1
  LIS_rc%met_proj = "none"

  call ESMF_ConfigFindLabel(LIS_config,"Met forcing sources:",rc=rc)
  do i=1,LIS_rc%nmetforc
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%metforc(i),rc=rc)
     call LIS_verify(rc,"Met forcing sources: not defined")
  enddo

!  call ESMF_ConfigFindLabel(LIS_config,"Met forcing chosen ensemble member:",rc=rc)
!  do i=1,LIS_rc%nmetforc
!     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%metforc_ensmem(i),rc=rc)
!     call LIS_warning(rc,"Met forcing chosen ensemble member: not defined")
!  enddo

  call ESMF_ConfigFindLabel(LIS_config,&
       "Topographic correction method (met forcing):",rc=rc)
  call LIS_verify(rc,'Topographic correction method (met forcing): not defined')
  do i=1,LIS_rc%nmetforc
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%met_ecor(i),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,&
       "Enable spatial downscaling of precipitation:",rc=rc)
  call LIS_verify(rc,'Enable spatial downscaling of precipitation: not defined')
  do i=1,LIS_rc%nmetforc
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%pcp_downscale(i),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,&
       "Spatial interpolation method (met forcing):",rc=rc)
  call LIS_verify(rc,'Spatial interpolation method (met forcing): not defined')
  do i=1,LIS_rc%nmetforc
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%met_interp(i),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,&
       "Spatial upscaling method (met forcing):",rc=rc)
  call LIS_verify(rc,'Spatial upscaling method (met forcing): not defined')
  do i=1,LIS_rc%nmetforc
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%met_upscale(i),rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,&
       "Temporal interpolation method (met forcing):",rc=rc)
  call LIS_verify(rc,'Temporal interpolation method (met forcing): not defined')
  do i=1,LIS_rc%nmetforc
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%met_tinterp(i),rc=rc)
  enddo

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%zterp_correction, &
       label="Enable new zterp correction (met forcing):",default=.false.,rc=rc)
  

  allocate(LIS_rc%nts(LIS_rc%nnest))

  allocate(LIS_rc%paramfile(LIS_rc%nnest))
  allocate(LIS_rc%usemaskmap(LIS_rc%nnest))
  allocate(LIS_rc%uselcmap(LIS_rc%nnest))
  allocate(LIS_rc%usetexturemap(LIS_rc%nnest))
  allocate(LIS_rc%usesoilfractionmap(LIS_rc%nnest))
  allocate(LIS_rc%usesoilcolormap(LIS_rc%nnest))
  allocate(LIS_rc%useelevationmap(LIS_rc%nnest))
  allocate(LIS_rc%useslopemap(LIS_rc%nnest))
  allocate(LIS_rc%useaspectmap(LIS_rc%nnest))
  allocate(LIS_rc%usecurvaturemap(LIS_rc%nnest))
  allocate(LIS_rc%uselaimap(LIS_rc%nnest))
  allocate(LIS_rc%usesaimap(LIS_rc%nnest))
  allocate(LIS_rc%usealbedomap(LIS_rc%nnest))
  allocate(LIS_rc%usemxsnalbmap(LIS_rc%nnest))
  allocate(LIS_rc%usegreennessmap(LIS_rc%nnest))
  allocate(LIS_rc%useroughnessmap(LIS_rc%nnest))
  allocate(LIS_rc%useemissmap(LIS_rc%nnest))
  allocate(LIS_rc%useporositymap(LIS_rc%nnest))
  allocate(LIS_rc%usepsisatmap(LIS_rc%nnest))
  allocate(LIS_rc%useksatmap(LIS_rc%nnest))
  allocate(LIS_rc%usebexpmap(LIS_rc%nnest))
  allocate(LIS_rc%usequartzmap(LIS_rc%nnest))
  allocate(LIS_rc%usesnowmap(LIS_rc%nnest))
  allocate(LIS_rc%snowsrc(LIS_rc%nnest))

  LIS_rc%snowsrc = 0 
  LIS_rc%usemaskmap = "none" 
  LIS_rc%uselcmap = "none" 
  LIS_rc%usetexturemap = "none" 
  LIS_rc%usesoilfractionmap = "none" 
  LIS_rc%usesoilcolormap = "none" 
  LIS_rc%useelevationmap = "none" 
  LIS_rc%useslopemap = "none" 
  LIS_rc%useaspectmap = "none" 
  LIS_rc%usecurvaturemap = "none" 
  LIS_rc%uselaimap = "none" 
  LIS_rc%usesaimap = "none" 
  LIS_rc%usealbedomap = "none" 
  LIS_rc%usemxsnalbmap = "none" 
  LIS_rc%usegreennessmap = "none" 
  LIS_rc%useroughnessmap = "none" 
  LIS_rc%useemissmap = "none"
  LIS_rc%useporositymap = "none" 
  LIS_rc%usepsisatmap = "none" 
  LIS_rc%useksatmap = "none" 
  LIS_rc%usebexpmap = "none" 
  LIS_rc%usequartzmap = "none" 
  LIS_rc%usesnowmap = 0  

  call ESMF_ConfigFindLabel(LIS_config,"LIS domain and parameter data file:",&
       rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%paramfile(i),rc=rc)
     call LIS_verify(rc,"LIS domain and parameter data file: not defined")
  enddo 

  call ESMF_ConfigFindLabel(LIS_config,"Landmask data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usemaskmap(i),rc=rc)
     call LIS_verify(rc,"Landmask data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Landcover data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%uselcmap(i),rc=rc)
     call LIS_verify(rc,"Landcover data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Soil texture data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usetexturemap(i),rc=rc)
     call LIS_verify(rc,"Soil texture data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Soil fraction data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usesoilfractionmap(i),rc=rc)
     call LIS_verify(rc,"Soil fraction data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Soil color data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usesoilcolormap(i),rc=rc)
     call LIS_verify(rc,"Soil color data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Elevation data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%useelevationmap(i),rc=rc)
     call LIS_verify(rc,"Elevation data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Slope data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%useslopemap(i),rc=rc)
     call LIS_verify(rc,"Slope data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Aspect data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%useaspectmap(i),rc=rc)
     call LIS_verify(rc,"Aspect data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Curvature data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usecurvaturemap(i),rc=rc)
     call LIS_verify(rc,"Curvature data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"LAI data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%useLAImap(i),rc=rc)
     call LIS_verify(rc,"LAI data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"SAI data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%useSAImap(i),rc=rc)
     call LIS_verify(rc,"SAI data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Albedo data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usealbedomap(i),rc=rc)
     call LIS_verify(rc,"Albedo data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Max snow albedo data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usemxsnalbmap(i),rc=rc)
     call LIS_verify(rc,"Max snow albedo data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Greenness data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usegreennessmap(i),rc=rc)
     call LIS_verify(rc,"Greenness data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Roughness data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%useroughnessmap(i),rc=rc)
     call LIS_verify(rc,"Roughness data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Porosity data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%useporositymap(i),rc=rc)
     call LIS_verify(rc,"Porosity data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Ksat data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%useksatmap(i),rc=rc)
     call LIS_verify(rc,"Ksat data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"B parameter data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usebexpmap(i),rc=rc)
     call LIS_verify(rc,"B parameter data source: not defined")
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Quartz data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%usequartzmap(i),rc=rc)
     call LIS_verify(rc,"Quartz data source: not defined")
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"Emissivity data source:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%useemissmap(i),rc=rc)
     call LIS_verify(rc,"Emissivity data source: not defined")
  enddo

  LIS_rc%tbot_update_lag = 0
  call ESMF_ConfigFindLabel(LIS_config,"TBOT lag skin temperature update option:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%tbot_update_lag,rc=rc)

  LIS_rc%tbot_lagday = 0
  call ESMF_ConfigFindLabel(LIS_config,"TBOT skin temperature lag days:",rc=rc)
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%tbot_lagday,rc=rc)

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%forcvarlistfile,&
       label="Forcing variables list file:",rc=rc)
  call LIS_verify(rc,"Forcing variables list file: not defined")

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%wopt,&
       label="Output methodology:",&
       rc=rc)
  call LIS_verify(rc,'Output methodology: not defined')

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%wopt_rst,&
       label="Output model restart files:",&
       rc=rc)
  call LIS_verify(rc,'Output model restart files: not defined')

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%wout,&
       label="Output data format:",&
       rc=rc)
  call LIS_verify(rc,'Output data format: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%wstyle,&
       label="Output naming style:",&
       rc=rc)
  call LIS_verify(rc,'Output naming style: not defined')
  !EMK Extra info required
  if (LIS_rc%wstyle == "557WW streamflow convention" .or. &
       LIS_rc%wstyle == "557WW medium range forecast convention") then
     call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%security_class, &
          label="AGRMET security classification:", rc=rc)
     call LIS_verify(rc, 'AGRMET security classification: option not specified in the config file')
     call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%distribution_class, &
          label="AGRMET distribution classification:", rc=rc)
     call LIS_verify(rc, 'AGRMET distribution classification: option not specified in the config file')
     call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%data_category, &
          label="AGRMET data category:", rc=rc)
     call LIS_verify(rc, 'AGRMET data category: option not specified in the config file')
     call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%area_of_data, &
          label="AGRMET area of data:", rc=rc)
     call LIS_verify(rc, 'AGRMET area of data: option not specified in the config file')
  endif

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%sout,&
       label="Enable output statistics:",default=.false.,&
       rc=rc)

  if (LIS_rc%wout.eq."grib1" .or. LIS_rc%wout.eq."grib2") then
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%grib_table,&
        label="Output GRIB Table Version:",rc=rc)
     call LIS_verify(rc,'Output GRIB Table Version: not defined')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%grib_center_id,&
        label="Output GRIB Center Id:",rc=rc)
     call LIS_verify(rc,'Output GRIB Center Id: not defined')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%grib_subcenter_id,&
        label="Output GRIB Subcenter Id:",rc=rc)
     call LIS_verify(rc,'Output GRIB Subcenter Id: not defined')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%grib_grid_id,&
        label="Output GRIB Grid Id:",rc=rc)
     call LIS_verify(rc,'Output GRIB Grid Id: not defined')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%grib_process_id,&
        label="Output GRIB Process Id:",rc=rc)
     call LIS_verify(rc,'Output GRIB Process Id: not defined')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%grib_packing_type,&
        label="Output GRIB Packing Type:",rc=rc)
     call LIS_verify(rc,'Output GRIB Packing Type: not defined')
  endif

!  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%plevel,label="Logging level:",&
!       rc=rc)
!  call LIS_verify(rc,'Logging level: not defined')

  LIS_rc%plevel = 1
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%startcode,label="Start mode:",&
       rc=rc)
  call LIS_verify(rc,'Start mode: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%syr,label="Starting year:",&
       rc=rc)
  call LIS_verify(rc,'Starting year: not defined')

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%smo,label="Starting month:",&
       rc=rc)
  call LIS_verify(rc,'Starting month: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%sda,label="Starting day:",&
       rc=rc)
  call LIS_verify(rc,'Starting day: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%shr,label="Starting hour:",&
       rc=rc)
  call LIS_verify(rc,'Starting hour: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%smn,label="Starting minute:",&
       rc=rc)
  call LIS_verify(rc,'Starting minute: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%sss,label="Starting second:",&
       rc=rc)
  call LIS_verify(rc,'Starting second: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%eyr1,label="Ending year:",&
       rc=rc)
  call LIS_verify(rc,'Ending year: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%emo1,label="Ending month:",&
       rc=rc)
  call LIS_verify(rc,'Ending month: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%eda1,label="Ending day:",&
       rc=rc)
  call LIS_verify(rc,'Ending day: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%ehr1,label="Ending hour:",&
       rc=rc)
  call LIS_verify(rc,'Ending hour: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%emn1,label="Ending minute:",&
       rc=rc)
  call LIS_verify(rc,'Ending minute: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%ess1,label="Ending second:",&
       rc=rc)
  call LIS_verify(rc,'Ending second: not defined')

  LIS_rc%eyr = LIS_rc%eyr1
  LIS_rc%emo = LIS_rc%emo1
  LIS_rc%eda = LIS_rc%eda1
  LIS_rc%ehr = LIS_rc%ehr1
  LIS_rc%emn = LIS_rc%emn1
  LIS_rc%ess = LIS_rc%ess1


! Time step is initialized to a maximum value of 1 day (in seconds). The time manager will be reinitialized to 
! the minimum timestep among different model components. 
  LIS_rc%ts = 86400.0
  LIS_rc%nts = 86400.0

  call ESMF_ConfigGetAttribute(LIS_config,time,&
       label="LIS time window interval:",default="1mo",rc=rc)

  call LIS_parseTimeString(time,LIS_rc%twInterval)

  allocate(LIS_rc%nensem(LIS_rc%nnest))

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%udef,label="Undefined value:",&
       rc=rc)
  call LIS_verify(rc,'Undefined value: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%odir,label="Output directory:",&
       rc=rc)
  call LIS_verify(rc,'Output directory: not defined')

  call ESMF_ConfigFindLabel(LIS_config,"Number of ensembles per tile:",rc=rc)
  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nensem(i),rc=rc)
     call LIS_verify(rc,'Number of ensembles per tile: not defined')
  enddo

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%surface_maxt,&
       label="Maximum number of surface type tiles per grid:",rc=rc)
  call LIS_verify(rc,'Maximum number of surface type tiles per grid: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%surface_minp,&
       label="Minimum cutoff percentage (surface type tiles):",rc=rc)
  call LIS_verify(rc,&
       'Minimum cutoff percentage (surface type tiles): not defined')

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%soilt_maxt,&
       label="Maximum number of soil texture tiles per grid:",rc=rc)
  call LIS_verify(rc,'Maximum number of soil texture tiles per grid: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%soilt_minp,&
       label="Minimum cutoff percentage (soil texture tiles):",rc=rc)
  call LIS_verify(rc,&
       'Minimum cutoff percentage (soil texture tiles): not defined')

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%soilf_maxt,&
       label="Maximum number of soil fraction tiles per grid:",rc=rc)
  call LIS_verify(rc,'Maximum number of soil fraction tiles per grid: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%soilf_minp,&
       label="Minimum cutoff percentage (soil fraction tiles):",rc=rc)
  call LIS_verify(rc,&
       'Minimum cutoff percentage (soil fraction tiles): not defined')

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%elev_maxt,&
       label="Maximum number of elevation bands per grid:",rc=rc)
  call LIS_verify(rc,'Maximum number of elevation bands per grid: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%elev_minp,&
       label="Minimum cutoff percentage (elevation bands):",rc=rc)
  call LIS_verify(rc,&
       'Minimum cutoff percentage (elevation bands): not defined')


  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%slope_maxt,&
       label="Maximum number of slope bands per grid:",rc=rc)
  call LIS_verify(rc,'Maximum number of slope bands per grid: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%slope_minp,&
       label="Minimum cutoff percentage (slope bands):",rc=rc)
  call LIS_verify(rc,&
       'Minimum cutoff percentage (slope bands): not defined')


  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%aspect_maxt,&
       label="Maximum number of aspect bands per grid:",rc=rc)
  call LIS_verify(rc,'Maximum number of aspect bands per grid: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%aspect_minp,&
       label="Minimum cutoff percentage (aspect bands):",rc=rc)
  call LIS_verify(rc,&
       'Minimum cutoff percentage (aspect bands): not defined')


  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%decompose_by_processes,&
       label="Decompose by processes:",default=.false.,rc=rc)
  if ( .not. LIS_rc%decompose_by_processes ) then
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%npesx,&
          label="Number of processors along x:",rc=rc)
     call LIS_verify(rc,'Number of processors along x: not defined')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%npesy,&
          label="Number of processors along y:",rc=rc)
     call LIS_verify(rc,'Number of processors along y: not defined')
  else
     ! Decomposing by processes.  Ignoring both
     !    "Number of processors along x:"
     ! and
     !    "Number of processors along y:"
     LIS_rc%npesx = 0
     LIS_rc%npesy = 0
  endif
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%halox,&
       label="Halo size along x:",rc=rc)
  call LIS_verify(rc,'Halo size along x: not defined')
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%haloy,&
       label="Halo size along y:",rc=rc)
  call LIS_verify(rc,'Halo size along y: not defined')
  if(LIS_npes.eq.1) then
     LIS_rc%npesx = 1
     LIS_rc%npesy = 1
  endif

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%ndas,&
       label="Number of data assimilation instances:",rc=rc)  
  call LIS_verify(rc,"Number of data assimilation instances: not defined")

  ninsts = max(1,LIS_rc%ndas)
!  ninsts_state  = max(LIS_rc%ndas, LIS_rc%npert_state)
!  ninsts_obs    = max(LIS_rc%ndas, LIS_rc%npert_obs)
!  ninsts = max(ninsts_state, ninsts_obs)

  allocate(LIS_rc%daalg(LIS_rc%ndas))
  allocate(LIS_rc%biasalg(LIS_rc%ndas))
  allocate(LIS_rc%biasrst(LIS_rc%ndas))
  allocate(LIS_rc%biasrstInterval(LIS_rc%ndas))
  allocate(LIS_rc%daset(ninsts))
  allocate(LIS_rc%dascaloption(LIS_rc%ndas))
  allocate(LIS_rc%nstvars(ninsts))
  allocate(LIS_rc%incroption(LIS_rc%ndas))
  allocate(LIS_rc%nobtypes(ninsts))
  allocate(LIS_rc%daoutInterval(LIS_rc%ndas))
  allocate(LIS_rc%perturb_obs(ninsts))
  allocate(LIS_rc%perturb_state(ninsts))
  allocate(LIS_rc%pertobsInterval(ninsts))
  allocate(LIS_rc%pertstateInterval(ninsts))
  allocate(LIS_rc%progpertattribFile(ninsts))
  allocate(LIS_rc%progattribFile(ninsts))
  allocate(LIS_rc%obspertattribFile(ninsts))
  allocate(LIS_rc%obsattribFile(ninsts))
  allocate(LIS_rc%biasOptionsFile(LIS_rc%ndas))  
  allocate(LIS_rc%biasrstFile(LIS_rc%ndas))  
  allocate(LIS_rc%pertrestartFile(LIS_rc%nnest))

  allocate(LIS_rc%wensems(LIS_rc%ndas))
  allocate(LIS_rc%wobs(LIS_rc%ndas))
  allocate(LIS_rc%winnov(LIS_rc%ndas))

  npert_forc = 0 
  npert_state = 0 
  npert_obs = 0 


  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%perturb_forcing,&
       label="Forcing perturbation algorithm:",rc=rc)
  call LIS_verify(rc,'Forcing perturbation algorithm: not defined')

  if(LIS_rc%perturb_forcing.ne."none") then 
     npert_forc = 1
  endif
  call ESMF_ConfigFindLabel(LIS_config, "Observation perturbation algorithm:",&
       rc=rc)
  do i=1,ninsts
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%perturb_obs(i),rc=rc)     
     call LIS_verify(rc,'Observation perturbation algorithm: not defined')
     if(LIS_rc%perturb_obs(i).ne."none") then 
        npert_obs = npert_obs+ 1
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"State perturbation algorithm:",rc=rc)
  do i=1,ninsts
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%perturb_state(i),rc=rc)
     call LIS_verify(rc,'State perturbation algorithm: not defined')
     if(LIS_rc%perturb_state(i).ne."none") then 
        npert_state = npert_state + 1
     endif
  enddo

  call ESMF_ConfigGetAttribute(LIS_config,time, &
       label="Forcing perturbation frequency:",rc=rc)
  if(LIS_rc%perturb_forcing.ne."none") then 
     call LIS_verify(rc,'Forcing perturbation frequency: not defined')  
     call LIS_parseTimeString(time,LIS_rc%pertforcInterval)
  endif

  call ESMF_ConfigFindLabel(LIS_config,"Observation perturbation frequency:",&
       rc=rc)
  do i=1,ninsts
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     if(LIS_rc%perturb_obs(i).ne."none") then 
        call LIS_verify(rc,'Observation perturbation frequency: not defined')
        
        call LIS_parseTimeString(time,LIS_rc%pertobsInterval(i))
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"State perturbation frequency:",rc=rc)
  do i=1,ninsts
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     if(LIS_rc%perturb_state(i).ne."none") then 
        call LIS_verify(rc,'State perturbation frequency: not defined')
        call LIS_parseTimeString(time,LIS_rc%pertstateInterval(i))
     endif
  enddo

!  if(npert_forc.ne.0.or.npert_state.ne.0.or.npert_obs.ne.0) then 
  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%pertrestart,&
       label="Perturbations start mode:",rc=rc)
  call LIS_verify(rc,'Perturbations start mode: not specified')
  
  call ESMF_ConfigGetAttribute(LIS_config,time,&
       label="Perturbations restart output interval:",rc=rc)
  call LIS_verify(rc,'Perturbations restart output interval: not specified')
  
  call LIS_parseTimeString(time,LIS_rc%pertrestartInterval)

  LIS_rc%pert_bias_corr = 1
!  if(npert_forc.ne.0.or.npert_state.ne.0) then 
!     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%pert_bias_corr,&
!          label="Apply perturbation bias correction:",rc=rc)
!     call LIS_verify(rc,'Apply perturbation bias correction: not specified')
!  endif

!  if(npert_forc.ne.0.or.npert_state.ne.0.or.npert_obs.ne.0) then 
  call ESMF_ConfigFindLabel(LIS_config,"Perturbations restart filename:",rc=rc)
  do i=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%pertrestartfile(i), rc=rc)
     call LIS_verify(rc,'Perturbations restart filename: not specified')
  enddo
!  endif

  if(npert_forc.ne.0) then 
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%forcattribfile,&
          label="Forcing attributes file:",rc=rc)
     call LIS_verify(rc,'Forcing attributes file: not defined')
     
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%forcpertattribfile,&
          label="Forcing perturbation attributes file:",rc=rc)
     call LIS_verify(rc,'Forcing perturbation attributes file: not defined')
  endif

!  if(npert_state.ne.0) then 
  call ESMF_ConfigFindLabel(LIS_config,"State attributes file:", rc=rc)
  do i=1,ninsts
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%progattribFile(i),rc=rc)
     call LIS_verify(rc,'State attributes file: not defined')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"State perturbation attributes file:",&
       rc=rc)
  do i=1,ninsts
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%progpertattribFile(i),rc=rc)
     call LIS_verify(rc,'State perturbation attributes file: not defined')
  enddo
!  endif
!  if(npert_state.ne.0.or.LIS_rc%ndas.gt.0) then 
!     call ESMF_ConfigFindLabel(LIS_config,"Number of state variables:",rc=rc)
!     do i=1,ninsts
!        call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nstvars(i),rc=rc)
!        call LIS_verify(rc,'Data assimilation number of state variables: not defined')
!     enddo
!  
!  endif

!  if(npert_obs.ne.0) then 
  call ESMF_ConfigFindLabel(LIS_config,&
       "Observation attributes file:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%obsattribFile(i),rc=rc)
     call LIS_verify(rc,'Observation attributes file: not defined')
  enddo
  call ESMF_ConfigFindLabel(LIS_config,&
       "Observation perturbation attributes file:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%obspertattribFile(i),rc=rc)
     call LIS_verify(rc,'Observation perturbation attributes file: not defined')
  enddo
!  endif

!  if(npert_forc.ne.0) then 
!     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nforcepert, &
!          label="Number of forcing fields to be perturbed:",rc=rc)
!     call LIS_verify(rc,'Number of forcing fields to be perturbed: not defined')
!  endif

  LIS_rc%nperts = 0 

  if(LIS_rc%ndas.gt.0) then 
     LIS_rc%nperts = LIS_rc%ndas
  else
     do i=1,ninsts
        if(LIS_rc%perturb_state(i).ne."none".or.&
             LIS_rc%perturb_obs(i).ne."none") then 
           LIS_rc%nperts = 1
        endif
     enddo
  endif

  call ESMF_ConfigFindLabel(LIS_config, "Data assimilation algorithm:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%daalg(i),rc=rc)
     call LIS_verify(rc,'Data assimilation algorithm: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Data assimilation set:",rc=rc)
  do i=1,ninsts
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%daset(i),rc=rc)
     call LIS_verify(rc,'Data assimilation set: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Data assimilation exclude analysis increments:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%incroption(i),rc=rc)
     if(LIS_rc%daalg(i).ne."none") then 
        call LIS_verify(rc,'Data assimilation exclude analysis increments: not defined')
     endif
  enddo
! increment option : 1 - exclude (apply only bias) 0 - include (apply bias + analysis)
! increment option : 1 - apply both, 0 - apply bias only. 

  call ESMF_ConfigFindLabel(LIS_config,"Data assimilation number of observation types:",rc=rc)
  do i=1,ninsts
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%nobtypes(i),rc=rc)
     call LIS_verify(rc,'Data assimilation number of observation types: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Data assimilation output interval for diagnostics:",&
       rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     if(LIS_rc%daalg(i).ne."none") then 
        call LIS_verify(rc,'Data assimilation output interval for diagnostics: not defined')

        call LIS_parseTimeString(time,LIS_rc%daoutInterval(i))
     endif     
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Data assimilation output ensemble spread:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%wensems(i),rc=rc)
     if(LIS_rc%daalg(i).ne."none") then 
        call LIS_verify(rc,'Data assimilation output ensemble spread: not defined')
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Data assimilation output processed observations:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%wobs(i),rc=rc)
     if(LIS_rc%daalg(i).ne."none") then 
        call LIS_verify(rc,'Data assimilation output processed observations: not defined')
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Data assimilation output innovations:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%winnov(i),rc=rc)
     if(LIS_rc%daalg(i).ne."none") then 
        call LIS_verify(rc,'Data assimilation output innovations: not defined')
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Data assimilation scaling strategy:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%dascaloption(i),rc=rc)
     if(LIS_rc%daalg(i).ne."none") then 
        call LIS_verify(rc,'Data assimilation scaling strategy: not defined')
     endif
  enddo


  call ESMF_ConfigFindLabel(LIS_config,"Bias estimation algorithm:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%biasalg(i),rc=rc)
     call LIS_verify(rc,'Bias estimation algorithm: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Bias estimation attributes file:",&
       rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%biasOptionsFile(i),rc=rc)
     if(LIS_rc%biasalg(i).ne."none") then 
        call LIS_verify(rc,'Bias estimation attributes file: not defined')
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Bias estimation restart output frequency:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     if(LIS_rc%biasalg(i).ne."none") then 
        call LIS_verify(rc,'Bias estimation restart output frequency: not defined')

        call LIS_parseTimeString(time,LIS_rc%biasrstInterval(i))
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Bias estimation start mode:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%biasrst(i),rc=rc)
     if(LIS_rc%biasalg(i).ne."none") then 
        call LIS_verify(rc,'Bias estimation start mode: not defined')
     endif
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"Bias estimation restart file:",rc=rc)
  do i=1,LIS_rc%ndas
     call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%biasrstFile(i),rc=rc)
     if(LIS_rc%biasalg(i).ne."none") then 
        call LIS_verify(rc,'Bias restart file: not defined')
     endif
  enddo
  
  call ESMF_ConfigGetAttribute(LIS_config, LIS_rc%routingmodel, &
       label="Routing model:",default="none", rc=rc)

 LIS_rc%endcode = 1

88 format(a4,25x,a3,5x,16a)
89 format(20x,a49)
!------------------------------------------------------------------------
! Set Time
!------------------------------------------------------------------------
  LIS_rc%YR= LIS_rc%SYR
  LIS_rc%MO= LIS_rc%SMO
  LIS_rc%DA= LIS_rc%SDA
  LIS_rc%HR= LIS_rc%SHR
  LIS_rc%MN= LIS_rc%SMN
  LIS_rc%SS= LIS_rc%SSS


  call LIS_date2time(LIS_rc%time,LIS_rc%doy,LIS_rc%gmt, &
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn,LIS_rc%ss)

  write(unit=LIS_logunit,FMT=*)'[INFO] *** NASA Land Information System (LIS) ***'
  write(unit=LIS_logunit,FMT=*)'[INFO] starting time: ',LIS_rc%smo,'/',LIS_rc%sda,'/',LIS_rc%syr
  write(unit=LIS_logunit,FMT=*)'[INFO] ending time: ',LIS_rc%emo,'/',LIS_rc%eda,'/',LIS_rc%eyr
  write(unit=LIS_logunit,FMT=*)'  '


  allocate(LIS_histData(LIS_rc%nnest))

  do i=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%syear,&
          label="Output start year:",default=LIS_rc%syr,rc=rc)
     call LIS_verify(rc,'Output start year: not specified')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%smonth,&
          label="Output start month:",default=LIS_rc%smo,rc=rc)
     call LIS_verify(rc,'Output start month: not specified')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%sday,&
          label="Output start day:",default=LIS_rc%sda,rc=rc)
     call LIS_verify(rc,'Output start day: not specified')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%shour,&
          label="Output start hour:",default=LIS_rc%shr,rc=rc)
     call LIS_verify(rc,'Output start hour: not specified')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%smin,&
          label="Output start minutes:",default=LIS_rc%smn,rc=rc)
     call LIS_verify(rc,'Output start minutes: not specified')
     call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%ssec,&
          label="Output start seconds:",default=LIS_rc%sss,rc=rc)
     call LIS_verify(rc,'Output start seconds: not specified')
     
     call LIS_date2time(LIS_histData(i)%time, doy, gmt, &
          LIS_histData(i)%syear,           &
          LIS_histData(i)%smonth,          &
          LIS_histData(i)%sday,            &
          LIS_histData(i)%shour,           &
          LIS_histData(i)%smin,            &
          LIS_histData(i)%ssec)

     write(unit=LIS_logunit,FMT=*) '[INFO] Output start time:',  &
          LIS_histData(i)%time,  &
          LIS_histData(i)%syear,  &
          LIS_histData(i)%smonth, &
          LIS_histData(i)%sday,   &
          LIS_histData(i)%shour,  &
          LIS_histData(i)%smin,   &
          LIS_histData(i)%ssec
  enddo

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%output_at_specifictime,&
       label="Output at specific time only:",default=0,rc=rc)

  if(LIS_rc%output_at_specifictime.eq.1) then 
     do i=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%month,&
             label="Specific output writing time (month):",default=-1,rc=rc)
        call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%day,&
             label="Specific output writing time (day):",default=-1,rc=rc)
        call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%hour,&
             label="Specific output writing time (hour):",default=-1,rc=rc)
        call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%min,&
             label="Specific output writing time (minute):",default=-1,rc=rc)
        call ESMF_ConfigGetAttribute(LIS_config,LIS_histData(i)%sec,&
             label="Specific output writing time (second):",default=-1,rc=rc)
     enddo
  endif

!
!------------------------------------------------------------------------
!  Setting Satellite LAI variables
!------------------------------------------------------------------------
  LIS_rc%laitime = 0.0
  LIS_rc%saitime = 0.0

  allocate(LIS_rc%rstflag(LIS_rc%nnest))
  allocate(LIS_rc%gridchange(LIS_rc%nnest))
  allocate(LIS_rc%tscount(LIS_rc%nnest))


  allocate(LIS_rc%mfile(LIS_rc%nnest))
  allocate(LIS_rc%vfile(LIS_rc%nnest))
  allocate(LIS_rc%vfile_form(LIS_rc%nnest))
  allocate(LIS_rc%elevfile(LIS_rc%nnest))
  allocate(LIS_rc%slfile(LIS_rc%nnest))
  allocate(LIS_rc%aspfile(LIS_rc%nnest))
  allocate(LIS_rc%curvfile(LIS_rc%nnest))
  allocate(LIS_rc%txtfile(LIS_rc%nnest))
  allocate(LIS_rc%safile(LIS_rc%nnest))
  allocate(LIS_rc%clfile(LIS_rc%nnest))
  allocate(LIS_rc%sifile(LIS_rc%nnest))
  allocate(LIS_rc%iscfile(LIS_rc%nnest))
  allocate(LIS_rc%pofile(LIS_rc%nnest))
  allocate(LIS_rc%psisatfile(LIS_rc%nnest))
  allocate(LIS_rc%ksatfile(LIS_rc%nnest))
  allocate(LIS_rc%bexpfile(LIS_rc%nnest))
  allocate(LIS_rc%qzfile(LIS_rc%nnest))
  allocate(LIS_rc%albfile(LIS_rc%nnest))
  allocate(LIS_rc%mxsnal(LIS_rc%nnest))
  allocate(LIS_rc%tbotfile(LIS_rc%nnest))
  allocate(LIS_rc%shdmaxfile(LIS_rc%nnest))
  allocate(LIS_rc%shdminfile(LIS_rc%nnest))
  allocate(LIS_rc%slopetypefile(LIS_rc%nnest))
!  allocate(LIS_rc%laifile(LIS_rc%nnest))
!  allocate(LIS_rc%saifile(LIS_rc%nnest))
!  allocate(LIS_rc%tile_coord_file(LIS_rc%nnest))
!  allocate(LIS_rc%tile_veg_file(LIS_rc%nnest))
  allocate(LIS_rc%outputSpecFile(LIS_rc%nnest))

  LIS_rc%rstflag = 1
  LIS_rc%gridchange = 1
  
!------------------------------------------------------------------------
! Select which vegetation tile space and mask files
!------------------------------------------------------------------------
  call LIS_initialize_registries

  allocate(LIS_rc%DAincrMode(LIS_rc%nnest))
  
  LIS_rc%DAincrMode   = 1
  LIS_rc%forecastMode = 0 

end subroutine LIS_readConfig
