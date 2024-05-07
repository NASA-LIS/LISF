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
! !ROUTINE: LVT_readConfig
!  \label{LVT_readConfig}
!
! !INTERFACE:
subroutine LVT_readConfig(configfile)
! 
! !USES:
  use ESMF 
  use LVT_coreMod
  use LVT_logMod
  use LVT_statsDataMod
  use LVT_runmode_pluginMod
  use LVT_metric_pluginMod
  use LVT_timeMgrMod
  use LVT_pluginIndices

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  The code in this file initializes the LVT configuration management utility.
!  The runtime specifications of a LVT simulation are read by this routine.
! 
! !FILES USED:
!
! !REVISION HISTORY:
!  02 Oct 2008: Sujay Kumar; Initial Specification
! 
!EOP
  
  character(len=*) :: configfile

  character*100 :: temp1
  character*1   :: fproc(4)
  logical       :: exists
  integer       :: i,k,rc
  integer       :: twsmooth
  character*30  :: scInterval
  character*10  :: time
  integer        :: ios,final_dirpos
  character(100) :: diag_fname
  character(100) :: diag_dir
  integer, external  :: LVT_create_subdirs

  inquire(file=(configfile), exist=exists)
  if ( .not. exists ) then
     write(LVT_logunit,*) '[ERR]: LVT config file ',      &
          (configfile), &
          ' does not exist'
     call LVT_endrun
  endif

  LVT_config = ESMF_ConfigCreate(rc=rc)
  call LVT_verify(rc,'problem in creating LVT_config object')

  call ESMF_ConfigLoadFile(LVT_config,trim(configfile),rc=rc)
  call LVT_verify(rc,'problem in loading lvt.config')

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%diagfile,&
       label="LVT diagnostic file:",rc=rc)
  if( rc > 0 ) then
    write(*,*) "[ERR] An error was returned in reading in "
    write(*,*) "[ERR]  the LVT diagnostic file ... "
    write(*,*) "[ERR] Please check your lvt.config file:  "
    write(*,*) "      ",trim(configfile)
    write(*,*) " Stopping ..."
    stop
  endif
  call LVT_verify(rc,'LVT diagnostic file: not defined')

! Make the diagnostic file directory names/path:
  diag_fname = LVT_rc%diagfile
  final_dirpos = scan(diag_fname, "/", BACK = .TRUE.)
  if(final_dirpos.ne.0) then 
     diag_dir = diag_fname(1:final_dirpos)
     ios = LVT_create_subdirs(len_trim(diag_dir),trim(diag_dir))
  endif

  write(unit=temp1,fmt='(i4.4)') LVT_localPet
  read(unit=temp1,fmt='(4a1)')fproc
  
  LVT_rc%diagfile = trim(LVT_rc%diagfile)//"."//fproc(1)//fproc(2)//fproc(3)//fproc(4)
  open(unit=LVT_logunit,file=(LVT_rc%diagfile))

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%runmode,&
       label="LVT running mode:",rc=rc)
  call LVT_verify(rc,'LVT running mode: option not specified in the config file')

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%domain,&
       label="Map projection of the LVT analysis:",rc=rc)
  if(rc.ne.0) then 
     write(LVT_logunit,*) &
          '[INFO] Map projection of the LVT analysis: option not specified in the config file'
     write(LVT_logunit,*) "[INFO] Please note that the option 'Map projection of the LIS run:' is"
     write(LVT_logunit,*) "[INFO] now deprecated. It should be replaced with the"
     write(LVT_logunit,*) "[INFO] entry - 'Map projection of the LVT analysis:' in the lvt.config file"
     call LVT_endrun()
  endif

  LVT_rc%max_model_types = 5
  LVT_rc%lsm_index = 1
  LVT_rc%lake_index = 2
  LVT_rc%glacier_index = 3
  LVT_rc%wetland_index = 4
  LVT_rc%openwater_index = 5
    

  LVT_rc%obs_duplicate = .false. 
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%nDataStreams,&
       label="Number of analysis sources:",default=2,rc=rc)
  call LVT_verify(rc,'Number of analysis sources: option not specified in the config file')

  allocate(LVT_rc%obssource(LVT_rc%nDataStreams))
  allocate(LVT_rc%resetFlag(LVT_rc%nDataStreams))
  LVT_rc%resetFlag = .true. 

  call ESMF_ConfigFindLabel(LVT_config,"Analysis data sources:",rc=rc)
  do i=1,LVT_rc%nDataStreams
     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%obssource(i),rc=rc)
     if(rc.ne.0) then 
        write(LVT_logunit,*) "[INFO] Analysis data sources: not defined in lvt.config"
        write(LVT_logunit,*) "[INFO] Please note that the option 'Observation source:' is"
        write(LVT_logunit,*) "[INFO] now deprecated. It should be replaced with the"
        write(LVT_logunit,*) "[INFO] entry - 'Analysis data sources:' in the lvt.config file"
        call LVT_endrun()
     endif
  enddo

  LVT_rc%lis_output_obs = .false. 
  do i=1,LVT_rc%nDataStreams
     if(LVT_rc%obssource(i).eq."LIS output") then 
        LVT_rc%lis_output_obs = .true. 
     endif
  enddo
  
  if(LVT_rc%obssource(1).eq.LVT_rc%obssource(2)) then 
     LVT_rc%obs_duplicate = .true. 
  endif

  if(LVT_rc%lis_output_obs) then 
     if((LVT_rc%obssource(1).ne."LIS output").and.&
          (LVT_rc%obssource(2).eq."LIS output")) then 
        write(LVT_logunit,*) "[ERR] Please specify 'LIS output' as the first"
        write(LVT_logunit,*) "[ERR] analysis data source when LIS output is"
        write(LVT_logunit,*) "[ERR] compared to a non-LIS data"
        call LVT_endrun()
     endif
  endif
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%lvt_out_format,&
       label="LVT output format:",rc=rc)
  call LVT_verify(rc,'LVT output format: option not specified in the config file')

  if (LVT_rc%lvt_out_format.eq."grib1" .or. LVT_rc%lvt_out_format.eq."grib2") then
     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%grib_table,&
          label="Output GRIB Table Version:",rc=rc)
     call LVT_verify(rc,'Output GRIB Table Version: not defined')
     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%grib_center_id,&
        label="Output GRIB Center Id:",rc=rc)
     call LVT_verify(rc,'Output GRIB Center Id: not defined')
     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%grib_subcenter_id,&
        label="Output GRIB Subcenter Id:",rc=rc)
     call LVT_verify(rc,'Output GRIB Subcenter Id: not defined')
     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%grib_grid_id,&
        label="Output GRIB Grid Id:",rc=rc)
     call LVT_verify(rc,'Output GRIB Grid Id: not defined')
     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%grib_process_id,&
        label="Output GRIB Process Id:",rc=rc)
     call LVT_verify(rc,'Output GRIB Process Id: not defined')
     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%grib_packing_type,&
        label="Output GRIB Packing Type:",rc=rc)
     call LVT_verify(rc,'Output GRIB Packing Type: not defined')
  endif

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%lvt_wopt,&
       label="LVT output methodology:",&
       rc=rc)
  call LVT_verify(rc,'LVT output methodology: not defined')



  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%syr,label="Starting year:",&
       rc=rc)
  call LVT_verify(rc,'Starting year: not defined')

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%smo,label="Starting month:",&
       rc=rc)
  call LVT_verify(rc,'Starting month: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%sda,label="Starting day:",&
       rc=rc)
  call LVT_verify(rc,'Starting day: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%shr,label="Starting hour:",&
       rc=rc)
  call LVT_verify(rc,'Starting hour: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%smn,label="Starting minute:",&
       rc=rc)
  call LVT_verify(rc,'Starting minute: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%sss,label="Starting second:",&
       rc=rc)
  call LVT_verify(rc,'Starting second: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%eyr1,label="Ending year:",&
       rc=rc)
  call LVT_verify(rc,'Ending year: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%emo1,label="Ending month:",&
       rc=rc)
  call LVT_verify(rc,'Ending month: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%eda1,label="Ending day:",&
       rc=rc)
  call LVT_verify(rc,'Ending day: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%ehr1,label="Ending hour:",&
       rc=rc)
  call LVT_verify(rc,'Ending hour: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%emn1,label="Ending minute:",&
       rc=rc)
  call LVT_verify(rc,'Ending minute: not defined')
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%ess1,label="Ending second:",&
       rc=rc)
  call LVT_verify(rc,'Ending second: not defined')

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%startmode,label="Start mode:",&
       rc=rc)
  call LVT_verify(rc,'Start mode: not defined')

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%wrst,&
       label="LVT output restart files:",default=0, rc=rc)
  call LVT_verify(rc,'LVT output restart files: not defined')

  LVT_rc%restartInterval = 0 

  if(LVT_rc%wrst.ne.0) then 
     call ESMF_ConfigGetAttribute(LVT_config,time,&
          label="LVT restart output interval:",rc=rc)
     call LVT_verify(rc,'LVT restart output interval: not defined')
     
     call LVT_parseTimeString(time,LVT_rc%restartInterval)
 
     if(LVT_rc%startmode.eq."restart") then 
        call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%rstfile,&
             label="LVT restart filename:",rc=rc)
        call LVT_verify(rc,'LVT restart filename: not defined')
     endif
  endif

  LVT_rc%eyr = LVT_rc%eyr1
  LVT_rc%emo = LVT_rc%emo1
  LVT_rc%eda = LVT_rc%eda1
  LVT_rc%ehr = LVT_rc%ehr1
  LVT_rc%emn = LVT_rc%emn1
  LVT_rc%ess = LVT_rc%ess1

  LVT_rc%YR= LVT_rc%SYR
  LVT_rc%MO= LVT_rc%SMO
  LVT_rc%DA= LVT_rc%SDA
  LVT_rc%HR= LVT_rc%SHR
  LVT_rc%MN= LVT_rc%SMN
  LVT_rc%SS= LVT_rc%SSS

  call ESMF_ConfigFindLabel(LVT_config,label="LVT clock timestep:",rc=rc)
  call ESMF_ConfigGetAttribute(LVT_config,time,rc=rc)
  if(rc.ne.0) then 
     write(LVT_logunit,*) '[INFO] LVT clock timestep: not defined'
       
     write(LVT_logunit,*) "[INFO] Please note that the option 'LIS output timestep:' is"
     write(LVT_logunit,*) "[INFO] now deprecated. It should be replaced with the"
     write(LVT_logunit,*) "[INFO] entry - 'LVT clock timestep:' in the lvt.config file"
     write(LVT_logunit,*) "[INFO] If LIS output is one of the observational sources, "
     write(LVT_logunit,*) "[INFO] please set the LVT clock timestep to be the same"
     write(LVT_logunit,*) "[INFO] as that of the LIS output frequency"
     call LVT_endrun()
  endif

  if(time.eq."dekad") then 
     LVT_rc%ts = 864000  !10 days
     LVT_rc%tsconv = "dekad"
  else
     LVT_rc%tsconv = "regular"
     call LVT_parseTimeString(time,LVT_rc%ts)
  endif

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%udef,label="Undefined value:",&
       rc=rc)
  call LVT_verify(rc,'Undefined value: not defined')

!  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%npesx,&
!       label="Number of processors along x:",rc=rc)
!  call LVT_verify(rc,'Number of processors along x: not defined')
!  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%npesy,&
!       label="Number of processors along y:",rc=rc)
!  call LVT_verify(rc,'Number of processors along y: not defined')

  
  LVT_rc%npesx = 1
  LVT_rc%npesy = 1

  LVT_rc%halox = 0 
  LVT_rc%haloy = 0 

  if(LVT_npes.eq.1) then
     LVT_rc%npesx = 1
     LVT_rc%npesy = 1
  endif

  if(LVT_rc%npesx*LVT_rc%npesy.ne.LVT_npes) then 
     write(unit=LVT_logunit,fmt=*) '[ERR] Layout does not match the number of processors...'
     write(unit=LVT_logunit,fmt=*) '[ERR] npex, npey, ',LVT_rc%npesx,'x',LVT_rc%npesy,'!=',LVT_npes
     write(unit=LVT_logunit,fmt=*) '[ERR] Stopping program ..'
     call LVT_endrun()
  endif

  call ESMF_ConfigGetAttribute(LVT_config,time, &
       label="Metrics computation frequency:",&
       rc=rc)
  if(rc.ne.0) then 
     write(LVT_logunit,*) "[INFO] Please note that the option 'Temporal averaging interval:' is"
     write(LVT_logunit,*) "[INFO] now deprecated. It should be replaced with the"
     write(LVT_logunit,*) "[INFO] entry - 'Metrics computation frequency:' in the lvt.config file"
     call LVT_endrun()
  endif

  if(time.eq."dekad") then 
     LVT_rc%tavgInterval= 864000  !10 days
  else
     call LVT_parseTimeString(time,LVT_rc%tavgInterval)
  endif

  allocate(LVT_rc%tlag(LVT_rc%nDataStreams))

  call ESMF_ConfigFindLabel(LVT_config,"Temporal lag in metrics computations:",rc=rc)
  do i=1,LVT_rc%nDataStreams
     call ESMF_ConfigGetAttribute(LVT_config,time,default="0ss",rc=rc)
     if(time.eq."dekad") then 
        write(LVT_logunit,*) '[WARN] dekad option is not supported for temporal lag '
     else
        call LVT_parseTimeString(time,LVT_rc%tlag(i))
     endif

  enddo
     
  call ESMF_ConfigGetAttribute(LVT_config,time, &
       label="Metrics output frequency:",&
       rc=rc)
  if(rc.ne.0) then 
     write(LVT_logunit,*) "[INFO] Please note that the option 'Stats output interval:' is"
     write(LVT_logunit,*) "[INFO] now deprecated. It should be replaced with the"
     write(LVT_logunit,*) "[INFO] entry - 'Metrics output frequency:' in the lvt.config file"
     call LVT_endrun()
  endif

  if(time.eq."dekad") then 
     LVT_rc%statswriteint = 864000  !10 days 
  else
     call LVT_parseTimeString(time,LVT_rc%statswriteint)
  endif
  
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%smoothObs, &
       label="Apply temporal smoothing to obs:",&
       default = 0,&
       rc=rc)
  call LVT_verify(rc,'Apply temporal smoothing to obs: not defined')

  if(LVT_rc%smoothObs.gt.0) then 
     call ESMF_ConfigGetAttribute(LVT_config,time,&
          label="Obs temporal smoothing window half length:",&
          rc=rc)
     call LVT_verify(rc,'Obs temporal smoothing window half length: not defined')

     call LVT_parseTimeString(time, twsmooth)
     call ESMF_TimeIntervalSet(LVT_obsSmTwL, s=twsmooth, &
          rc=rc)
     call LVT_verify(rc, 'TimeIntervalSet failed in LVT_readConfig')

     call ESMF_ConfigGetAttribute(LVT_config,time,&
          label="Obs temporal smoothing window interval:",&
          rc=rc)
     call LVT_verify(rc,'Obs temporal smoothing window interval: not defined')

     call LVT_parseTimeString(time, twsmooth)
          call ESMF_TimeIntervalSet(LVT_obsSmTwI, s=twsmooth, &
          rc=rc)
     call LVT_verify(rc, 'TimeIntervalSet failed in LVT_readConfig')
  endif


  LVT_rc%use_shift_mo = 1
  if(LVT_rc%tavgInterval.eq.31536000.0.or.&
       LVT_rc%statswriteint.eq.31536000) then 
      call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%use_shift_mo, &
           label="Starting month if a shifted year definition is used in temporal averaging:",&
           default=1,&
           rc=rc)
      call LVT_verify(rc,'Starting month if a shifted year definition is used in temporal averaging: not defined')
   endif
  
   call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%dataMask,&
        label="Apply external mask:",&
        default = 0,&
        rc=rc)
  call LVT_verify(rc,'Apply external mask: not defined')

  if(LVT_rc%dataMask.ne.0) then 
     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%maskdir,&
          label="External mask directory:",&
          rc=rc)
     call LVT_verify(rc,'External mask directory: not defined')
  endif
  
!monthly seasonal mask
  if(LVT_rc%dataMask.eq.3) then 

     call ESMF_ConfigFindLabel(LVT_config,"Temporal (monthly) mask flags:",rc=rc)
     do i=1,12
        call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%monthly_mask(i),rc=rc)
        call LVT_verify(rc,'Temporal (monthly) mask flags: option not specified in the config file')
     enddo
  endif

  call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%data_based_strat, &
       label="External data-based stratification: ",default=0,rc=rc)
  call LVT_verify(rc,"External data-based stratification: not defined")

  call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%var_based_strat, &
       label="Variable-based stratification: ",default =0, rc=rc)
  call LVT_verify(rc,"Variable-based stratification: not defined")

  if(LVT_rc%var_based_strat.ge.1) then 
     call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%vname_strat, &
          label="Stratification variable:",rc=rc)
     call LVT_verify(rc,&
          "Stratification variable: not defined")
     
     call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%strat_var_threshold, &
          label="Stratification threshold: ",rc=rc)
     call LVT_verify(rc,&
          "Stratification threshold: not defined")
  endif

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%obsCountThreshold,&
       label="Observation count threshold:",&
       default=0,&
       rc=rc)
  call LVT_verify(rc,'Observation count threshold: not defined')

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%statsodir,&
       label="Metrics output directory:",&
       rc=rc)
  if(rc.ne.0) then
     write(LVT_logunit,*) "[INFO] Please note that the option 'Stats output directory:' is"
     write(LVT_logunit,*) "[INFO] now deprecated. It should be replaced with the"
     write(LVT_logunit,*) "[INFO] entry - 'Metrics output directory:' in the lvt.config file"
     call LVT_endrun()
  endif
  call LVT_verify(rc,'Metrics output directory: not defined')

  call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%nensem,&
       label="Number of ensembles in the LVT analysis:", default=1, rc=rc)
  
  call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%noebal_refet,&
       label="Calculate reference ET without energy balance:", default=0, rc=rc)

  if(LVT_rc%runmode.eq.LVT_DataCompId) then 

     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%computeEnsMetrics,&
          label="Compute ensemble metrics:",default=0,rc=rc)
     call LVT_verify(rc,"Compute ensemble metrics: not defined")
     
     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%computeICMetrics,&
          label="Compute information theory metrics:",default=0,rc=rc)
     call LVT_verify(rc,"Compute information theory metrics: not defined")
  

!----------------------------------------------------------------------------
! If information theory metrics are enabled, then other options are turned
! off
!----------------------------------------------------------------------------
     if(LVT_rc%computeICmetrics.eq.1) then 
        LVT_rc%computeEnsMetrics = 0
     endif
     
     call ESMF_ConfigGetAttribute(LVT_config, scInterval, &
          label="Seasonal cycle interval type:",&
          default="monthly",&
          rc=rc)
     call LVT_verify(rc,'Seasonal cycle interval type: not defined')
     if(scInterval.eq."monthly") then 
        LVT_rc%scInterval = 1
     elseif(scInterval.eq."3 monthly") then 
        LVT_rc%scInterval = 2
     elseif(scInterval.eq."6 monthly") then 
        LVT_rc%scInterval = 6
     elseif(scInterval.eq."yearly") then 
        LVT_rc%scInterval = 12
     elseif(scInterval.eq."3 monthly WY") then
        LVT_rc%scInterval = 21 !scInterval is 2 with a variation of the 3 monthly
     else
        write(LVT_logunit,*) '[ERR] Seasonal cycle interval type must be -- '
        write(LVT_logunit,*) '[ERR] monthly, 3 monthly, 3 monthly WY, 6 monthly or yearly'
        call LVT_endrun()
     endif
     
  ! 1= month, 2-3 months. 
     call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%scCountThreshold, &
          label="Seasonal cycle minimum count threshold:",&
          default=0,&
          rc=rc)
     call LVT_verify(rc,'Seasonal cycle minimum count threshold: not defined')

     call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%adcCountThreshold, &
          label="Average diurnal cycle minimum count threshold:",&
          default=0,&
          rc=rc)
     call LVT_verify(rc,'Average diurnal cycle minimum count threshold: not defined')
  
     call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%pval_CI,&
          label="Confidence interval (%): ",default=95.0,rc=rc)
     call LVT_verify(rc,"Confidence interval (%): not defined")     
     LVT_rc%pval_CI = (100-LVT_rc%pval_CI)/100.0

  endif
! The innovation distribution can be computed only in the DA mode.   
!  if(LVT_rc%runmode.eq.LVT_dastatId) then 
!     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%computeInnovDist, &
!          label="Compute DA innovation distribution:", rc=rc)
!     call LVT_verify(rc, 'Compute DA innovation distribution: not defined')
!     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%computeGain, &
!          label="Compute DA analysis gain:", rc=rc)
!     call LVT_verify(rc, 'Compute DA analysis gain: not defined')
!     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%computeSpread, &
!          label="Compute DA ensemble spread:", rc=rc)
!     call LVT_verify(rc, 'Compute DA ensemble spread: not defined')
!     call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%nstvars, &
!          label="Number of state variables in the DA update:", rc=rc)
!     call LVT_verify(rc, 'Number of state variables in the DA update: not defined')
!  endif


  ! The following lines are added by Shugong Wang for backward support of LIS 6 
  call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%lis_version, &
       label="LIS version:", default=7, rc=rc)  ! 6, 7, ...
  ! call LVT_verify(rc,"LIS version: not defined")
  ! No need to vary. If there is no "LIS version", default version is 7.      
  ! read in "experiment code"
  if(LVT_rc%lis_version == 6) then ! read in experiment code 
     call ESMF_ConfigGetAttribute(LVT_config, LVT_rc%expcode, &
          label="Experiment code:", rc=rc)  
     call LVT_verify(rc,"Experiment code: not defined. It is required for LIS 6 results!")
     
  endif

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%lsm,&
       label="Land surface model:",rc=rc)  
  if(LVT_rc%lis_version == 6) then ! read in experiment code 
     call LVT_verify(rc,'Land surface model: option not specified in the config file. It is required for LIS 6 results!')
  endif
  ! End of code by Shugong Wang

  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%paramfile,&
       label="Input domain and mask data file:",rc=rc)
  if(rc.ne.0) then 
     write(LVT_logunit,*) &
          '[INFO] Input domain and mask data file: option not specified in the config file'
     write(LVT_logunit,*) "[INFO] Please note that the option 'LIS parameter data file:' is"
     write(LVT_logunit,*) "[INFO] now deprecated. It should be replaced with the"
     write(LVT_logunit,*) "[INFO] entry - 'Input domain and mask data file:' in the lvt.config file"
     call LVT_endrun()
  endif
    
  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%spTransformAnlysDomain,&
       label="Spatial transform method for generating the LVT analysis domain:",&
       default="none", rc=rc)
 
  call LVT_runmode_plugin
  call LVT_metric_plugin

  if(LVT_rc%statswriteint.lt.LVT_rc%tavgInterval) then 
     write(LVT_logunit,*) '[ERR] Stats output interval should not be '
     write(LVT_logunit,*) '[ERR] less than time averaging interval... '
     call LVT_endrun()
  endif

  LVT_rc%endcode = 1

!  if((LVT_rc%anlys_data_class).eq."LSM") then 
!     LVT_rc%model_name = (LVT_rc%lsm)
!  elseif((LVT_rc%anlys_data_class).eq."Routing") then 
!     LVT_rc%model_name =(LVT_rc%routing_model)
!  elseif((LVT_rc%anlys_data_class).eq."RTM") then 
!     LVT_rc%model_name = (LVT_rc%rtm_model)
!  endif


!  call ESMF_ConfigFindLabel(LVT_config,"LVT datastream attributes file:",&
!       rc=rc)
!  if(rc.ne.0) then 
!     write(LVT_logunit,*) &
!          '[INFO] LVT datastream attributes file: option not specified in the config file'
!     write(LVT_logunit,*) "[INFO] Please note that the option 'LIS output attributes file:' is"
!     write(LVT_logunit,*) "[INFO] now deprecated. It should be replaced with the"
!     write(LVT_logunit,*) "[INFO] entry - 'LVT datastream attributes file:' in the lvt.config file"
!     call LVT_endrun()
!  endif
!  call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%outputSpecFile,rc=rc)
!  write(LVT_logunit,*) '[INFO] Opening Model Output Attributes File', LVT_rc%outputSpecFile

!The following are being hardcoded since AFWA is the only organization that 
!employs the following settings

  LVT_rc%security_class='U'
  LVT_rc%distribution_class = 'C'
  !LVT_rc%data_category = 'ANLYS'
  LVT_rc%data_category = 'C'
  LVT_rc%area_of_data = 'GLOBAL'

  allocate(LVT_metricsPtr(LVT_NMETRICS))
  LVT_rc%maskflag = 0

! to remove eventually
  LVT_rc%nnest = 1
  LVT_rc%odir = ""

  LVT_rc%curr_pass = 1
  LVT_rc%pass = 1
  if(LVT_rc%runmode.eq.LVT_benchMarkId) then 
     LVT_rc%pass = 2
  endif
end subroutine LVT_readConfig
