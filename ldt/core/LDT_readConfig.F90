!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: LDT_readConfig
!  \label{LDT_readConfig}
!
! !REVISION HISTORY:
!
!  02 Oct 2008: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine LDT_readConfig(configfile)

! !USES:
  use ESMF 
  use LDT_coreMod,         only : LDT_rc, LDT_config, LDT_localPet, LDT_npes
  use LDT_DAobsDataMod
  use LDT_DAmetricsDataMod,only : LDT_DAmetricsPtr
  use LDT_logMod,          only : LDT_logunit, LDT_verify, LDT_endrun, LDT_warning
  use LDT_runmode_pluginMod, only : LDT_runmode_plugin
  use LDT_param_pluginMod
  use LDT_ANNdata_pluginMod
  use LDT_metforcScale_pluginMod
  use LDT_pluginIndices

!
! !DESCRIPTION:
!
!  The code in this file initializes the LDT configuration management utility.
!  The runtime specifications of a LDT simulation
!  are read by this routine.
!EOP
  implicit none
  
  character(len=*), intent(in) :: configfile

  integer        :: rc
  integer        :: n, i, k
  character*100  :: temp1
  character*1    :: fproc(4)

  integer        :: ios
  integer        :: final_dirpos
  character(100) :: diag_fname
  character(100) :: diag_dir
  integer, external :: LDT_create_subdirs

!____________________________________________________________

  LDT_config = ESMF_ConfigCreate(rc=rc)
  call LDT_verify(rc,'problem in creating LDT_config object')

  call ESMF_ConfigLoadFile(LDT_config,trim(configfile),rc=rc)
  call LDT_verify(rc,'problem in loading ldt.config')

!== Read in common LDT config inputs: ==

  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%runmode,&
       label="LDT running mode:",rc=rc)
  call LDT_verify(rc,'LDT running mode: option not specified in the config file')

! Read in and Set LDT run-time log-diagnostic filename path:
  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%diagfile,&
       label="LDT diagnostic file:",rc=rc)
  call LDT_verify(rc,'LDT diagnostic file: not defined')

! Make the diagnostic file directory names/path:
  diag_fname = LDT_rc%diagfile
  final_dirpos = scan(diag_fname, "/", BACK = .TRUE.)
  if(final_dirpos.ne.0) then
     diag_dir = diag_fname(1:final_dirpos)
     ios = LDT_create_subdirs(len_trim(diag_dir),trim(diag_dir))
  endif

  write(unit=temp1,fmt='(i4.4)') LDT_localPet
  read (unit=temp1,fmt='(4a1)' ) fproc

  LDT_rc%diagfile = trim(LDT_rc%diagfile)//"."//fproc(1)//fproc(2)//fproc(3)//fproc(4)
  open(unit=LDT_logunit,file=trim(LDT_rc%diagfile))

! ---

! Read in Mask-Parameter Fill diagnostic filepath name:
  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%mpfillfile,&
       label="Mask-parameter fill diagnostic file:", &
       default="MaskParamFill.log",rc=rc)
  call LDT_verify(rc,'Mask-parameter fill diagnostic file: not defined')

! Make the diagnostic file directory names/path:
  diag_fname = LDT_rc%mpfillfile
  final_dirpos = scan(diag_fname, "/", BACK = .TRUE.)
  if(final_dirpos.ne.0) then
     diag_dir = diag_fname(1:final_dirpos)
     ios = LDT_create_subdirs(len_trim(diag_dir),trim(diag_dir))
  endif

! ---

  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%odir,&
       label="LDT output directory:",&
       rc=rc)
  call LDT_verify(rc,'LDT output directory: not defined')
!  write(unit=LDT_logunit,fmt=*) "LDT output directory: ",LDT_rc%odir

  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%nnest,&
       label="LIS number of nests:",rc=rc)
  call LDT_verify(rc,'LIS number of nests: option not specified in the config file')

!- Surface type model options:

  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%nsf_model_types,&
       label="Number of surface model types:",rc=rc)
  call LDT_verify(rc,'Number of surface model types: option not specified in the config file')

  allocate(LDT_rc%sf_model_type_name_select(LDT_rc%nsf_model_types))
  allocate(LDT_rc%sf_model_type_select(LDT_rc%nsf_model_types))

  call ESMF_ConfigFindLabel(LDT_config,"Surface model types:",rc=rc)
  do i=1,LDT_rc%nsf_model_types
     call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%sf_model_type_name_select(i),rc=rc)
     call LDT_verify(rc,"Surface model types: not defined")

     call LDT_mapSurfaceModelType(LDT_rc%sf_model_type_name_select(i), &
          LDT_rc%sf_model_type_select(i))
  enddo

  do k = 1, LDT_rc%nsf_model_types
    if( LDT_rc%sf_model_type_name_select(k) .eq. "LSM"  .or. &
        LDT_rc%sf_model_type_name_select(k) .eq. "Lake" .or. &
        LDT_rc%sf_model_type_name_select(k) .eq. "Glacier" .or. &
        LDT_rc%sf_model_type_name_select(k) .eq. "Openwater" ) then
    else
      write(*,*) "[ERR] Only 'LSM', 'Lake', 'Glacier' or 'Openwater' surface "
      write(*,*) "    model type options are currently available. Please"
      write(*,*) "    select one of these options for now."
      write(*,*) "Stopping LDT run ..."
      call LDT_endrun
    endif
  enddo

  LDT_rc%max_model_types = 5

  LDT_rc%lsm_index       = 1
  LDT_rc%lake_index      = 2
  LDT_rc%glacier_index   = 3
  LDT_rc%wetland_index   = 4
  LDT_rc%openwater_index = 5

  allocate(LDT_rc%sf_model_type_name(LDT_rc%max_model_types))
  allocate(LDT_rc%sf_model_type(LDT_rc%max_model_types))

  LDT_rc%sf_model_type_name(1) = "LSM"
  LDT_rc%sf_model_type(1)      = LDT_rc%lsm_index
  LDT_rc%sf_model_type_name(2) = "Lake"
  LDT_rc%sf_model_type(2)      = LDT_rc%lake_index
  LDT_rc%sf_model_type_name(3) = "Glacier"
  LDT_rc%sf_model_type(3)      = LDT_rc%glacier_index
  LDT_rc%sf_model_type_name(4) = "Wetland"
  LDT_rc%sf_model_type(4)      = LDT_rc%wetland_index
  LDT_rc%sf_model_type_name(5) = "Openwater"
  LDT_rc%sf_model_type(5)      = LDT_rc%openwater_index

!- LSM parameters:
   call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%lsm,&
        label="Land surface model:",default="none",rc=rc)
   call LDT_verify(rc,'Land surface model: option not specified in the config file')

!- Lake parameters:
   call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%lakemodel,&
        label="Lake model:",default="none",rc=rc)
   call LDT_verify(rc,'Lake model: option not specified in the config file')

   call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%routingmodel,&
        label="Routing model:",default="none",rc=rc)
   call LDT_verify(rc,'Routing model: option not specified in the config file')

 ! LSM/Crop-specific entries:
   allocate(LDT_rc%assimcropinfo(LDT_rc%nnest))
   LDT_rc%assimcropinfo = .false.
   if( trim(LDT_rc%runmode) == "LSM parameter processing" ) then
     call ESMF_ConfigFindLabel(LDT_config,"Incorporate crop information:",rc=rc)
     do n=1,LDT_rc%nnest
        call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%assimcropinfo(n),rc=rc)
        call LDT_verify(rc,'Incorporate crop information: not specified')
     enddo
   endif

!-- Meteorological Forcing inputs: --

   LDT_rc%nmetforc = 0
   call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%nmetforc,&
        label="Number of met forcing sources:",rc=rc)
   call LDT_verify(rc,'Number of met forcing sources: not specified')

 ! Determine sources of metforcing options selected:
   if( LDT_rc%nmetforc > 0 ) then

     allocate(LDT_rc%metforc(LDT_rc%nmetforc))
     LDT_rc%metforc = "none"
     LDT_rc%nmetforc_parms = 0

   ! Read in Met forcing sources:
     call ESMF_ConfigFindLabel(LDT_config,"Met forcing sources:",rc=rc)
     do i=1,LDT_rc%nmetforc
        call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%metforc(i),rc=rc)
        call LDT_verify(rc,"Met forcing sources: not defined")

        if( LDT_rc%metforc(i) == "GDAS" ) then
           LDT_rc%nmetforc_parms = LDT_rc%nmetforc_parms + 6
        elseif( LDT_rc%metforc(i) == "ECMWF" ) then
           LDT_rc%nmetforc_parms = LDT_rc%nmetforc_parms + 8
        else
           LDT_rc%nmetforc_parms = LDT_rc%nmetforc_parms + 1
        endif
     end do

   ! Forcing-array
     allocate(LDT_rc%met_gridtransform(LDT_rc%nmetforc))
     allocate(LDT_rc%met_ecor(LDT_rc%nmetforc))
     allocate(LDT_rc%met_nf(LDT_rc%nmetforc))
     allocate(LDT_rc%met_ts(LDT_rc%nmetforc))
     allocate(LDT_rc%met_tinterp(LDT_rc%nmetforc))
     allocate(LDT_rc%met_zterp(LDT_rc%nmetforc))
     allocate(LDT_rc%met_validhr(LDT_rc%nmetforc))
     allocate(LDT_rc%met_nensem(LDT_rc%nmetforc))

   ! Forcing Parameter-array
     allocate(LDT_rc%metforc_parms(LDT_rc%nmetforc_parms))
     allocate(LDT_rc%metforc_parmsrc(LDT_rc%nmetforc_parms))
     allocate(LDT_rc%met_gridtransform_parms(LDT_rc%nmetforc_parms))
     allocate(LDT_rc%met_ecor_parms(LDT_rc%nmetforc_parms))
     allocate(LDT_rc%met_gridDesc(LDT_rc%nmetforc_parms,20))
     allocate(LDT_rc%met_proj(LDT_rc%nmetforc_parms))
     allocate(LDT_rc%met_nc(LDT_rc%nmetforc_parms))
     allocate(LDT_rc%met_nr(LDT_rc%nmetforc_parms))

     ! Default initialization
     LDT_rc%met_nensem = 1

     LDT_rc%met_gridtransform = "none"
     LDT_rc%met_gridtransform_parms = "none"
     call ESMF_ConfigFindLabel(LDT_config,&
          "Met spatial transform methods:",rc=rc)
     call LDT_verify(rc,'Met spatial transform methods: not defined')
     do i=1,LDT_rc%nmetforc
        call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%met_gridtransform(i),rc=rc)
     enddo

     LDT_rc%met_ecor = "none"
     LDT_rc%met_ecor_parms = "none"
     call ESMF_ConfigFindLabel(LDT_config,&
          "Topographic correction method (met forcing):",rc=rc)
     call LDT_verify(rc,'Topographic correction method (met forcing): not defined')
     do i=1,LDT_rc%nmetforc
        call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%met_ecor(i),rc=rc)
     enddo

     k = 0
     do i = 1, LDT_rc%nmetforc 
     !- Account for GDAS' multi-grids:
        if( LDT_rc%metforc(i) == "GDAS" ) then
           k = k + 1
           LDT_rc%metforc_parms(k+0) = "GDAS_T126" 
           LDT_rc%metforc_parms(k+1) = "GDAS_T170" 
           LDT_rc%metforc_parms(k+2) = "GDAS_T254" 
           LDT_rc%metforc_parms(k+3) = "GDAS_T382" 
           LDT_rc%metforc_parms(k+4) = "GDAS_T574" 
           LDT_rc%metforc_parms(k+5) = "GDAS_T1534" 
           LDT_rc%met_gridtransform_parms(k:k+5) = LDT_rc%met_gridtransform(i)
           LDT_rc%met_ecor_parms(k:k+5) = LDT_rc%met_ecor(i)
           LDT_rc%metforc_parmsrc(k:k+5) = LDT_rc%metforc_parms(k:k+5)
           k = k + 5
     !- Account for ECMWF' multi-grids:
        elseif( LDT_rc%metforc(i) == "ECMWF" ) then
           k = k + 1
           LDT_rc%metforc_parms(k+0) = "ECMWF_S23R4"
           LDT_rc%metforc_parms(k+1) = "ECMWF_S25R1"
           LDT_rc%metforc_parms(k+2) = "ECMWF_S30R1"
           LDT_rc%metforc_parms(k+3) = "ECMWF_S33R1"
           LDT_rc%metforc_parms(k+4) = "ECMWF_S35R2"
           LDT_rc%metforc_parms(k+5) = "ECMWF_S35R3"
           LDT_rc%metforc_parms(k+6) = "ECMWF_S36R1"
           LDT_rc%metforc_parms(k+7) = "ECMWF_S37R2"
           LDT_rc%met_gridtransform_parms(k:k+7) = LDT_rc%met_gridtransform(i)
           LDT_rc%met_ecor_parms(k:k+7) = LDT_rc%met_ecor(i)
           LDT_rc%metforc_parmsrc(k:k+7) = LDT_rc%metforc_parms(k:k+7)
           k = k + 7
     !- All other forcing grids:
        else
           k = k + 1
           LDT_rc%metforc_parms(k) = LDT_rc%metforc(i)
           LDT_rc%met_ecor_parms(k) = LDT_rc%met_ecor(i)
           LDT_rc%met_gridtransform_parms(k) = LDT_rc%met_gridtransform(i)
           LDT_rc%metforc_parmsrc(k) = LDT_rc%metforc_parms(k)

         ! Reassign forcing dataset names with "()" for the output *netcdf file:
           select case( LDT_rc%metforc_parms(k) )
             case( "AGRMET radiation (polar stereographic)" )
               LDT_rc%metforc_parmsrc(k) = "AGRMET_radiation_ps"
             case( "AGRMET radiation (latlon)" )
               LDT_rc%metforc_parmsrc(k) = "AGRMET_radiation_latlon"
             case( "GEOS5 forecast" )
               LDT_rc%metforc_parmsrc(k) = "GEOS5_fcst"
             case( "GDAS(LSWG)" )
               LDT_rc%metforc_parmsrc(k) = "GDAS_LSWG"
             case( "TRMM 3B42V6" )
               LDT_rc%metforc_parmsrc(k) = "TRMM_3B42V6"
             case( "TRMM 3B42V7" )
               LDT_rc%metforc_parmsrc(k) = "TRMM_3B42V7"
             case( "TRMM 3B42RTV7" )
               LDT_rc%metforc_parmsrc(k) = "TRMM_3B42RTV7"
             case( "CPC CMORPH" )
               LDT_rc%metforc_parmsrc(k) = "CPC_CMORPH"
             case( "CPC STAGEII" )
               LDT_rc%metforc_parmsrc(k) = "CPC_STAGEII"
             case( "CPC STAGEIV" )
               LDT_rc%metforc_parmsrc(k) = "CPC_STAGEIV"
             case( "RFE2(daily)" )
               LDT_rc%metforc_parmsrc(k) = "RFE2_daily"
             case( "RFE2(gdas)" )
               LDT_rc%metforc_parmsrc(k) = "RFE2_gdas"
             case( "Noah Bondville" )
               LDT_rc%metforc_parmsrc(k) = "Noah_Nondville"
             case( "FASST test" )
               LDT_rc%metforc_parmsrc(k) = "FASST_test"
             case( "TRIGRS test" )
               LDT_rc%metforc_parmsrc(k) = "TRIGRS_test"
             case( "Rhone AGG" )
               LDT_rc%metforc_parmsrc(k) = "Rhone_AGG"
             case( "VIC processed forcing" )
               LDT_rc%metforc_parmsrc(k) = "VIC_forcing"
             case default
             ! Assign Metforcing sources from readin index options:
               LDT_rc%metforc_parmsrc(k) = LDT_rc%metforc_parms(k)
           end select

        endif
     end do
!     do i = 1, LDT_rc%nmetforc; print *, i, LDT_rc%metforc(i); enddo
!     do i = 1, LDT_rc%nmetforc_parms; print *, i, LDT_rc%metforc_parms(i); enddo

   ! Read in entries for metforcing processing only:
     LDT_rc%met_zterp = .false.

     call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%zterp_correction, &
          label="Enable new zterp correction (met forcing):",default=.false.,rc=rc)

     if( trim(LDT_rc%runmode) == "Metforce processing" .or. &
         trim(LDT_rc%runmode) == "Metforce temporal downscaling" .or. &
         trim(LDT_rc%runmode) == "Statistical downscaling of met forcing"  ) then

        LDT_rc%metforc_blend_alg = "none"
        call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%metforc_blend_alg,&
             label="Blending method for forcings:",rc=rc)
        call LDT_verify(rc,'Blending method for forcings: not specified')
 
        LDT_rc%met_tinterp = "none"
        call ESMF_ConfigFindLabel(LDT_config,&
            "Temporal interpolation method (met forcing):",rc=rc)
        call LDT_verify(rc,'Temporal interpolation method (met forcing): not defined')
        do i=1,LDT_rc%nmetforc
           call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%met_tinterp(i),rc=rc)
        enddo

     endif

  end if  ! End met forcing check

! ________________________________

  allocate(LDT_rc%lis_map_proj(LDT_rc%nnest))
  allocate(LDT_rc%lis_map_resfactor(LDT_rc%nnest))

  call ESMF_ConfigFindLabel(LDT_config, &
       "Map projection of the LIS domain:",rc=rc)
  do n=1,LDT_rc%nnest     
     call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%lis_map_proj(n),rc=rc)
     call LDT_verify(rc,'Map projection of the LIS domain: option not specified in the config file')
  enddo

  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%udef,label="Undefined value:",&
       rc=rc)
  call LDT_verify(rc,'Undefined value: not defined')

! Option to add buffer around parameter grid domain:
  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%add_buffer,&
       label="Add buffer to parameter grid domain:",&
       default=0,rc=rc)
  call LDT_verify(rc,'Add buffer to parameter grid domain: not defined')

  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%x_buffer,&
       label="Buffer count in x-direction:",&
       default=5,rc=rc)
  call LDT_verify(rc,'Buffer count in x-direction: not defined')

  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%y_buffer,&
       label="Buffer count in y-direction:",&
       default=5,rc=rc)
  call LDT_verify(rc,'Buffer count in y-direction: not defined')


! Set number of processors along x and y directions:
  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%npesx,&
       label="Number of processors along x:", &
       default=1,rc=rc)
  call LDT_verify(rc,'Number of processors along x: not defined')
  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%npesy,&
       label="Number of processors along y:", &
       default=1,rc=rc)
  call LDT_verify(rc,'Number of processors along y: not defined')

  LDT_rc%halox = 0
  LDT_rc%haloy = 0

  if(LDT_npes.eq.1) then
     LDT_rc%npesx = 1
     LDT_rc%npesy = 1
  endif

  if(LDT_rc%npesx*LDT_rc%npesy.ne.LDT_npes) then
     write(unit=LDT_logunit,fmt=*) "Layout does not match the number of processors ..."
     write(unit=LDT_logunit,fmt=*) "npex, npey, ",LDT_rc%npesx,"x",LDT_rc%npesy,"!=",LDT_npes
     write(unit=LDT_logunit,fmt=*) "Stopping program ..."
     call LDT_endrun()
  endif

! Output file writing formats and directory/file structures:
  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%wopt,&
       label="Output methodology:",&
       default="2d gridspace",&
       rc=rc)
  call LDT_verify(rc,'Output methodology: not defined')

  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%wout,&
       label="Output data format:",&
       default="netcdf",&
       rc=rc)
  call LDT_verify(rc,'Output data format: not defined')

  call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%wstyle,&
       label="Output naming style:",&
       default="3 level hierarchy",&
       rc=rc)
  call LDT_verify(rc,'Output naming style: not defined')

! --------------------------------
! Allocate LDT config inputs for number of nests:
! LIS run domain grid:
  allocate(LDT_rc%gridDesc(LDT_rc%nnest,20))

! Main parameter domain grids:
  allocate(LDT_rc%lc_gridtransform(LDT_rc%nnest))    
  allocate(LDT_rc%mask_gridtransform(LDT_rc%nnest))          
  allocate(LDT_rc%reg_gridtransform(LDT_rc%nnest))           
  allocate(LDT_rc%soils_gridtransform(LDT_rc%nnest))         
  allocate(LDT_rc%soiltext_gridtransform(LDT_rc%nnest))      
  allocate(LDT_rc%topo_gridtransform(LDT_rc%nnest))          

  allocate(LDT_rc%lc_type(LDT_rc%nnest))
  allocate(LDT_rc%mask_source(LDT_rc%nnest))
  allocate(LDT_rc%soil_classification(LDT_rc%nnest))
  allocate(LDT_rc%mask_type(LDT_rc%nnest))
  allocate(LDT_rc%gridcell_water_frac(LDT_rc%nnest))
  allocate(LDT_rc%gridcell_glacier_frac(LDT_rc%nnest))

  allocate(LDT_rc%cliplandmask(LDT_rc%nnest))

  allocate(LDT_rc%mfile(LDT_rc%nnest))
  allocate(LDT_rc%vfile(LDT_rc%nnest))
  allocate(LDT_rc%regfile(LDT_rc%nnest))
  allocate(LDT_rc%elevfile(LDT_rc%nnest))
  allocate(LDT_rc%slfile(LDT_rc%nnest))
  allocate(LDT_rc%aspfile(LDT_rc%nnest))
  allocate(LDT_rc%curvfile(LDT_rc%nnest))

  allocate(LDT_rc%txtfile(LDT_rc%nnest))
  allocate(LDT_rc%safile(LDT_rc%nnest))
  allocate(LDT_rc%clfile(LDT_rc%nnest))
  allocate(LDT_rc%sifile(LDT_rc%nnest))
  allocate(LDT_rc%gravelfile(LDT_rc%nnest))
  allocate(LDT_rc%iscfile(LDT_rc%nnest))
  allocate(LDT_rc%pofile(LDT_rc%nnest))
  allocate(LDT_rc%psisatfile(LDT_rc%nnest))
  allocate(LDT_rc%ksatfile(LDT_rc%nnest))
  allocate(LDT_rc%bexpfile(LDT_rc%nnest))
  allocate(LDT_rc%qzfile(LDT_rc%nnest))
  allocate(LDT_rc%dsoilfile(LDT_rc%nnest))
  allocate(LDT_rc%bdrckdepfile(LDT_rc%nnest))
  allocate(LDT_rc%hsgfile(LDT_rc%nnest))
  allocate(LDT_rc%bulkdensfile(LDT_rc%nnest))
  allocate(LDT_rc%domrockfile(LDT_rc%nnest))
  allocate(LDT_rc%rockvolfile(LDT_rc%nnest))
  allocate(LDT_rc%awcfile(LDT_rc%nnest))
  allocate(LDT_rc%permabfile(LDT_rc%nnest))

!  allocate(LDT_rc%tbotfile(LDT_rc%nnest))
!  allocate(LDT_rc%slopetypefile(LDT_rc%nnest))

  allocate(LDT_rc%gdasT126elevfile(LDT_rc%nnest))
  allocate(LDT_rc%gdasT170elevfile(LDT_rc%nnest))
  allocate(LDT_rc%gdasT254elevfile(LDT_rc%nnest))
  allocate(LDT_rc%gdasT382elevfile(LDT_rc%nnest))
  allocate(LDT_rc%gdasT574elevfile(LDT_rc%nnest))
  allocate(LDT_rc%gdasT1534elevfile(LDT_rc%nnest))

  allocate(LDT_rc%elevfileifs23r4(LDT_rc%nnest))
  allocate(LDT_rc%elevfileifs25r1(LDT_rc%nnest))
  allocate(LDT_rc%elevfileifs30r1(LDT_rc%nnest))
  allocate(LDT_rc%elevfileifs33r1(LDT_rc%nnest))
  allocate(LDT_rc%elevfileifs35r2(LDT_rc%nnest))
  allocate(LDT_rc%elevfileifs35r3(LDT_rc%nnest))
  allocate(LDT_rc%elevfileifs36r1(LDT_rc%nnest))
  allocate(LDT_rc%elevfileifs37r2(LDT_rc%nnest))

  allocate(LDT_rc%outputSpecFile(LDT_rc%nnest))

  allocate(LDT_rc%monthlyData(LDT_rc%nnest))
  allocate(LDT_rc%quarterlyData(LDT_rc%nnest))
  allocate(LDT_rc%nts(LDT_rc%nnest))
  allocate(LDT_rc%tscount(LDT_rc%nnest))
  allocate(LDT_rc%rstflag(LDT_rc%nnest))
  allocate(LDT_rc%gridchange(LDT_rc%nnest))

  allocate(LDT_rc%glaciermask(LDT_rc%nnest))

  LDT_rc%tscount = 0 
  LDT_rc%rstflag = 1
  LDT_rc%gridchange = 1

!== Read in LDT config model parameter inputs: ==

!- Call each run mode and parameter module plugin routine:
  call LDT_runmode_plugin
  call LDT_landcover_plugin
  call LDT_soils_plugin
  call LDT_topo_plugin    
  call LDT_laisai_plugin
  call LDT_gfrac_plugin
  call LDT_irrigation_plugin
  call LDT_alb_plugin

  call LDT_climate_plugin
  call LDT_forcingparams_plugin
  call LDT_LSMparam_plugin
  call LDT_routingparam_plugin
  call LDT_lakeparam_plugin

  call LDT_glacier_plugin

  call LDT_ANNinputdata_plugin
  call LDT_ANNoutputdata_plugin

  call LDT_timedscale_plugin

! by default set the number of passes through the time loop as 1

  allocate(LDT_DAmetricsPtr(LDT_DA_MOC_COUNT))

  LDT_rc%endcode = 1
  LDT_rc%pass = 1 
  LDT_rc%wout_form = 1
  LDT_rc%lis_wopt = 2
  LDT_rc%monthlyData = .false. 
  LDT_rc%quarterlyData = .false. 

end subroutine LDT_readConfig
