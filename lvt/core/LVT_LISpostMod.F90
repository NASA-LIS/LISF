!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
#include "LVT_NetCDF_inc.h"
!BOP
! 
! !MODULE: LVT_LISpostMod
! \label(LVT_domanMod)
!
! !INTERFACE:
module LVT_LISpostMod
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_fileIOMod
  use LVT_LISoutputHandlerMod
  use map_utils
  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
  PRIVATE
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! The code in this file provides interfaces to manages different running
! domain implementations
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initialize_LISpost
  public :: LVT_process_LISoutput

  public :: LVT_LISpost
  
!EOP

  type, public :: lispost_type_dec
     integer   :: npes
     integer   :: nfields
     integer   :: nlevels
     integer   :: nsoilmoistlayers, nsoiltemplayers
     real, pointer :: soilmoistlyrthk(:), soiltemplyrthk(:)
  end type lispost_type_dec
  
  type(lispost_type_dec) :: LVT_LISpost
  
contains
!BOP
! 
! !ROUTINE: LVT_initialize_LISpost
! \label{LVT_initialize_LISpost}
!
! !INTERFACE: 
  subroutine LVT_initialize_LISpost()
! 
! !USES:

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! 
!EOP

    integer           :: i, rc

    call ESMF_ConfigGetAttribute(LVT_config,LVT_LISpost%npes,&
         label="Number of processors used in the LIS output generation:",rc=rc)
    call LVT_verify(rc,'Number of processors used in the LIS output generation: option not specified in the config file')
    call ESMF_ConfigGetAttribute(LVT_config,LVT_LISpost%nfields,&
         label="Number of fields in the LIS output:",rc=rc)
    call LVT_verify(rc,'Number of fields in the LIS output: option not specified in the config file')
    call ESMF_ConfigGetAttribute(LVT_config,LVT_LISpost%nlevels,&
         label="Number of vertical levels in the LIS output:",rc=rc)
    call LVT_verify(rc,'Number of vertical levels in the LIS output: option not specified in the config file')
    call ESMF_ConfigGetAttribute(LVT_config, LVT_LISpost%nsoilmoistlayers, &
         label="LIS output number of soil moisture layers:", default=0, rc=rc)
    call ESMF_ConfigGetAttribute(LVT_config, LVT_LISpost%nsoiltemplayers, &
         label="LIS output number of soil temperature layers:", default=0, rc=rc)
    
    if (LVT_LISpost%nsoilmoistlayers .gt. 0) then
       allocate(LVT_LISpost%soilmoistlyrthk(LVT_LISpost%nsoilmoistlayers))

       call ESMF_ConfigFindLabel(LVT_config, "LIS output soil moisture layer thickness:", rc=rc)
       call LVT_verify(rc, 'LIS output soil moisture layer thickness: option not specified in the config file')

       do i=1, LVT_LISpost%nsoilmoistlayers
          call ESMF_ConfigGetAttribute(LVT_config,LVT_LISpost%soilmoistlyrthk(i), rc=rc)
          call LVT_verify(rc, 'The number of soil moisture layers indicated in the  &
               config file is greater than the number of thickness values provided')
       enddo
    endif
    
    if (LVT_LISpost%nsoiltemplayers .gt. 0) then
       allocate(LVT_LISpost%soiltemplyrthk(LVT_LISpost%nsoiltemplayers))

       call ESMF_ConfigFindLabel(LVT_config, "LIS output soil temperature layer thickness:", rc=rc)
       call LVT_verify(rc, 'LIS output soil temperature layer thickness: option not specified in the config file')

       do i=1, LVT_LISpost%nsoilmoistlayers
          call ESMF_ConfigGetAttribute(LVT_config,LVT_LISpost%soiltemplyrthk(i), rc=rc)
          call LVT_verify(rc, 'The number of soil temperature layers indicated in the & 
               config file is greater than the number of thickness values provided')
       enddo
    endif

  end subroutine LVT_initialize_LISpost

!BOP
! 
! !ROUTINE: LVT_process_LISoutput
! \label{LVT_process_LISoutput}
!
! !INTERFACE:   
  subroutine LVT_process_LISoutput

!
!DESCRIPTION:
! This routine reads the distributed binary LIS outputs and quilts
! them into a single NetCDF file.
!
    integer       :: i,l,m,t,c,r
    integer       :: source
    integer       :: ftn(LVT_LISpost%npes),ftn_nc
    integer       :: lnc(LVT_LISpost%npes),lnr(LVT_LISpost%npes)
    logical       :: file_exists
    character*200 :: froot,temp1,fname_out
    character*100  :: vname,units
    character*200  :: fname(LVT_LISpost%npes)
    real,          allocatable :: var(:,:,:,:)
    real           :: gvar(LVT_rc%gnc,LVT_rc%gnr,&
         LVT_LIS_rc(1)%nensem,LVT_LISpost%nlevels)
    real           :: lat(LVT_rc%gnr)
    real           :: lon(LVT_rc%gnc)
    character*1    :: fproc(4)
    integer        :: dimID(4),tdimID,xtimeID,varId
    integer        :: latid, lonid
    integer        :: iret
    integer        :: nlevels
    logical        :: read_flag
    integer        :: ews_ind(LVT_LISpost%npes), ewe_ind(LVT_LISpost%npes)
    integer        :: nss_ind(LVT_LISpost%npes), nse_ind(LVT_LISpost%npes)
    integer        :: fill_value
    integer :: shuffle, deflate, deflate_level    

    source = 1
    if(LVT_LIS_rc(1)%anlys_data_class.eq."LSM") then       
       call LVT_create_output_filename(LVT_LIS_rc(1)%nest, &
            source,froot, &
            'SURFACEMODEL',&
            writeint = LVT_LIS_rc(1)%ts, &
            wout = LVT_LIS_rc(1)%format,&
            style=LVT_LIS_rc(1)%style, &
            odir=LVT_LIS_rc(1)%odir)

       call LVT_create_output_filename(LVT_LIS_rc(1)%nest, &
            source,fname_out, &
            'SURFACEMODEL',&
            writeint = LVT_LIS_rc(1)%ts, &
            wout = "netcdf",&
            style=LVT_LIS_rc(1)%style, &
            odir=LVT_LIS_rc(1)%odir)
       
    elseif(LVT_LIS_rc(1)%anlys_data_class.eq."Routing") then 
       call LVT_create_output_filename(LVT_LIS_rc(1)%nest, &
            source, froot, &
            'ROUTING',&
            writeint = LVT_rc%ts, wout = LVT_LIS_rc(1)%format,&
            style=LVT_LIS_rc(1)%style,&
            odir=LVT_LIS_rc(1)%odir)
       call LVT_create_output_filename(LVT_LIS_rc(1)%nest, &
            source, fname_out, &
            'ROUTING',&
            writeint = LVT_rc%ts, &
            wout = "netcdf",&
            style=LVT_LIS_rc(1)%style,&
            odir=LVT_LIS_rc(1)%odir)
       
    elseif(LVT_LIS_rc(1)%anlys_data_class.eq."RTM") then 
       call LVT_create_output_filename(LVT_LIS_rc(1)%nest,&
            source, froot, &
            'RTM',&
            writeint = LVT_rc%ts, wout = LVT_LIS_rc(1)%format,&
            style=LVT_LIS_rc(1)%style, &
            odir=LVT_LIS_rc(1)%odir)
    elseif(LVT_LIS_rc(1)%anlys_data_class.eq."Irrigation") then 
       call LVT_create_output_filename(LVT_LIS_rc(1)%nest, &
            source, froot, &
            'IRRIGATION',&
            writeint = LVT_rc%ts, wout = LVT_LIS_rc(1)%format,&
            style=LVT_LIS_rc(1)%style,&
            odir=LVT_LIS_rc(1)%odir)
    endif

    if(LVT_LIS_rc(1)%format.eq."distributed binary") then

       read_flag = .true.
       
       do i=1,LVT_LISpost%npes
          write(temp1,'(i4.4)') i-1
          read(temp1,fmt='(4a1)') fproc
          
          fname(i) = trim(froot)//'.'//fproc(1)//fproc(2)//fproc(3)//fproc(4)
          
          inquire(file=trim(fname(i)),exist=file_exists)
       
          if(file_exists) then 
             
             ! write(LVT_logunit,*) '[INFO] Reading LIS output ',trim(fname(i))

             ftn(i) = LVT_getNextUnitNumber()
             open(ftn(i),file=trim(fname(i)),form='unformatted')
             read(ftn(i)) ews_ind(i)
             read(ftn(i)) ewe_ind(i)
             read(ftn(i)) nss_ind(i)
             read(ftn(i)) nse_ind(i)

             lnc(i) = (ewe_ind(i) - ews_ind(i))+1
             lnr(i) = (nse_ind(i) - nss_ind(i))+1
          else
             read_flag = .false. 
          endif
       enddo

       if(read_flag) then

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          iret = nf90_create(path=fname_out,cmode=nf90_hdf5,&
               ncid = ftn_nc)
          call LVT_verify(iret,'creating netcdf file failed in LVT_LISpostMod')
          
          write(LVT_logunit,*) '[INFO] writing '//trim(fname_out)
          if(LVT_LIS_rc(1)%wopt.eq."2d gridspace") then
             call LVT_verify(nf90_def_dim(ftn_nc,'east_west',LVT_rc%gnc,&
                  dimID(1)),&
                  'nf90_def_dim for east_west failed in LVT_LISpostMod')
             call LVT_verify(nf90_def_dim(ftn_nc,'north_south',LVT_rc%gnr,&
                  dimID(2)),&
                  'nf90_def_dim for north_south failed in LVT_LISpostMod')
             call LVT_verify(nf90_def_dim(ftn_nc,'SoilMoist_profiles',&
                  LVT_LISpost%nlevels,&
                  dimID(3)),&
                  'nf90_def_dim for SoilMoist_profiles failed in LVT_LISpostMod')        
          elseif(LVT_LIS_rc(1)%wopt.eq."2d ensemble gridspace") then 
             call LVT_verify(nf90_def_dim(ftn_nc,'east_west',LVT_rc%gnc,&
                  dimID(1)),&
                  'nf90_def_dim for east_west failed in LVT_LISpostMod')
             call LVT_verify(nf90_def_dim(ftn_nc,'north_south',LVT_rc%gnr,&
                  dimID(2)),&
                  'nf90_def_dim for north_south failed in LVT_LISpostMod')
             call LVT_verify(nf90_def_dim(ftn_nc,'ensemble',LVT_LIS_rc(1)%nensem,&
                  dimID(3)),&
                  'nf90_def_dim for ensemble failed in LVT_LISpostMod')
             call LVT_verify(nf90_def_dim(ftn_nc,'SoilMoist_profiles',&
                  LVT_LISpost%nlevels,&
                  dimID(4)),&
                  'nf90_def_dim for SoilMoist_profiles failed in LVT_LISpostMod')
             
          endif
          
          call LVT_verify(nf90_def_dim(ftn_nc,'time',1,tdimID),&
               'nf90_def_dim for time failed in LVT_LISpostMod')
          call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"missing_value", LVT_rc%udef),&
               'nf90_put_att for missing_value failed in LVT_LISpostMod')
          
          call LVT_verify(nf90_def_var(ftn_nc,'time',&
               nf90_float,dimids = tdimID, varID=xtimeID),&
               'nf90_def_var for time failed in LVT_LISpostMod')
          
          call LVT_verify(nf90_def_var(ftn_nc,'lat',&
               nf90_float,dimids=dimID(2),varid=latid),&
               'nf90_def_att for lat failed in LVT_LISpostMod')
          
          call LVT_verify(nf90_def_var(ftn_nc,'lon',&
               nf90_float,dimids=dimID(1),varid=lonid),&
               'nf90_def_att for lat failed in LVT_LISpostMod')    
          
    ! Global attributes
          if (LVT_LISpost%nsoilmoistlayers .gt. 0) then
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"NUM_SOIL_MOIST_LAYERS", &
                  LVT_LISpost%nsoilmoistlayers),&
                  'nf90_put_att for NUM_SOIL_MOIST_LAYERS failed in LVT_LISpostMod')
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOIL_MOIST_LAYER_THICKNESSES", &
                  LVT_LISpost%soilmoistlyrthk),&
                  'nf90_put_att for SOIL_MOIST_LAYER_THICKNESSES failed in LVT_LISpostMod')
          endif
          
          if (LVT_LISpost%nsoiltemplayers .gt. 0) then
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"NUM_SOIL_TEMP_LAYERS", &
                  LVT_LISpost%nsoiltemplayers),&
                  'nf90_put_att for NUM_SOIL_TEMP_LAYERS failed in LVT_LISpostMod')
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOIL_TEMP_LAYER_THICKNESSES", &
                  LVT_LISpost%soiltemplyrthk),&
                  'nf90_put_att for SOIL_TEMP_LAYER_THICKNESSES failed in LVT_LISpostMod')
          endif

          call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"title", &
               "LIS land surface model output"),&
               'nf90_put_att for title failed in LVT_LISpostMod')
          call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"institution", &
               "NASA GSFC"),&
               'nf90_put_att for institution failed in LVT_LISpostMod')
          call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"source",&
               trim(LVT_LIS_rc(1)%model_name)),&
               'nf90_put_att for source failed in LVT_LISpostMod')
!    call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"history", &
!         "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
!         date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)),&
!         'nf90_put_att for history failed in LVT_LISpostMod')
          call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"references", &
               "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"),&
               'nf90_put_att for references failed in LVT_LISpostMod')
          call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"conventions", &
               "CF-1.6"),'nf90_put_att for conventions failed in LVT_LISpostMod')
          call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"comment", &
               "website: http://lis.gsfc.nasa.gov/"),&
               'nf90_put_att for comment failed in LVT_LISpostMod')
          
          if(trim(LVT_rc%domain).eq."latlon") then !latlon
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"MAP_PROJECTION", &
                  "EQUIDISTANT CYLINDRICAL"))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(9)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(10)))       
             
          elseif(trim(LVT_rc%domain).eq."mercator") then
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"MAP_PROJECTION", &
                  "MERCATOR"))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))       
          elseif(trim(LVT_rc%domain).eq."lambert") then
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"MAP_PROJECTION", &
                  "LAMBERT CONFORMAL"))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
            LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"TRUELAT2", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
             
          elseif(trim(LVT_rc%domain).eq."polar") then
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"MAP_PROJECTION", &
                  "POLAR STEREOGRAPHIC"))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
                  LVT_rc%gridDesc(4)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
                  LVT_rc%gridDesc(5)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"TRUELAT1", &
                  LVT_rc%gridDesc(10)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"ORIENT", &
                  LVT_rc%gridDesc(7)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"STANDARD_LON", &
                  LVT_rc%gridDesc(11)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"DX", &
                  LVT_rc%gridDesc(8)))
             call LVT_verify(nf90_put_att(ftn_nc,NF90_GLOBAL,"DY", &
                  LVT_rc%gridDesc(9)))
             
          endif
          
          call LVT_verify(nf90_enddef(ftn_nc),'nf90_enddef failed in LVT_LISpostMod')
          
          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                call ij_to_latlon(LVT_domain%lvtproj, &
                     float(c), float(r), lat(r), lon(c))
             enddo
          enddo
          
          call LVT_verify(nf90_put_var(ftn_nc,xtimeID,0.0),&
               'nf90_put_var failed for time in LVT_LISpostMod')
          call LVT_verify(nf90_put_var(ftn_nc,latid,lat),&
               'nf90_put_var failed for lat in LVT_LISpostMod')    
          call LVT_verify(nf90_put_var(ftn_nc,lonid,lon),&
               'nf90_put_var failed for lat in LVT_LISpostMod')
          
#endif              
          
          do t=1,LVT_LISpost%nfields
             
             do i=1,LVT_LISPost%npes
                read(ftn(i)) nlevels
                read(ftn(i)) vname
                read(ftn(i)) units
             
                if(nlevels.eq.1) then 
                   allocate(var(lnc(i),lnr(i),LVT_LIS_rc(1)%nensem,1))
                   
                   do m=1,LVT_LIS_rc(1)%nensem
                      read(ftn(i)) var(:,:,m,1)
                   enddo
                   
                   gvar(ews_ind(i):ewe_ind(i),nss_ind(i):nse_ind(i),:,1) = &
                        var(:,:,:,1)
                   deallocate(var)
                else
                   allocate(var(lnc(i),lnr(i),LVT_LIS_rc(1)%nensem,nlevels))
                   
                   do l=1,nlevels                                            
                      do m=1,LVT_LIS_rc(1)%nensem
                         read(ftn(i)) var(:,:,m,l)
                      enddo
                   enddo
                   
                   gvar(ews_ind(i):ewe_ind(i),nss_ind(i):nse_ind(i),:,:) = &
                        var(:,:,:,:)
                   deallocate(var)

                endif

             enddo
             
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
             shuffle = NETCDF_shuffle
             deflate = NETCDF_deflate
             deflate_level =NETCDF_deflate_level
             
             call LVT_verify(nf90_redef(ftn_nc),&
                  'nf90_redef failed in LVT_LISpostMod')
             
             if(LVT_LIS_rc(1)%wopt.eq."2d gridspace") then

                if(nlevels.eq.1) then 
                   call LVT_verify(nf90_def_var(ftn_nc,&
                        trim(vname),&
                        nf90_float,&
                        dimids = dimID(1:2), varID=varId),&
                        'nf90_def_var for '//trim(vname)//&
                        'failed in LVT_LISpostMod')
                   
!                   call LVT_verify(nf90_def_var_fill(ftn_nc,&
!                        varId, &
!                        1,fill_value), 'nf90_def_var_fill failed for '//&
!                        vname)
                   
                   call LVT_verify(nf90_def_var_deflate(ftn_nc,&
                        varId,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(vname)//&
                        'failed in LVT_LISpostMod')
                else
                   call LVT_verify(nf90_def_var(ftn_nc,&
                        trim(vname),&
                        nf90_float,&
                        dimids = dimID(1:3), varID=varId),&
                        'nf90_def_var for '//trim(vname)//&
                        'failed in LVT_LISpostMod')
                   
!                   call LVT_verify(nf90_def_var_fill(ftn_nc,&
!                        varId, &
!                        1,fill_value), 'nf90_def_var_fill failed for '//&
!                        vname)
                   
                   call LVT_verify(nf90_def_var_deflate(ftn_nc,&
                        varId,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(vname)//&
                        'failed in LVT_LISpostMod')

                endif
                
             elseif(LVT_LIS_rc(1)%wopt.eq."2d ensemble gridspace") then

                if(nlevels.eq.1) then 
                   call LVT_verify(nf90_def_var(ftn_nc,&
                        trim(vname),&
                        nf90_float,&
                        dimids = dimID(1:3), varID=varId),&
                        'nf90_def_var for '//trim(vname)//&
                        'failed in LVT_LISpostMod')
                   
!                   call LVT_verify(nf90_def_var_fill(ftn_nc,&
!                        varId, &
!                        1,fill_value), 'nf90_def_var_fill failed for '//&
!                        vname)
                   
                   call LVT_verify(nf90_def_var_deflate(ftn_nc,&
                        varId,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(vname)//&
                        'failed in LVT_LISpostMod')
                else
                   call LVT_verify(nf90_def_var(ftn_nc,&
                        trim(vname),&
                        nf90_float,&
                        dimids = dimID(1:4), varID=varId),&
                        'nf90_def_var for '//trim(vname)//&
                        'failed in LVT_LISpostMod')
                   
!                   call LVT_verify(nf90_def_var_fill(ftn_nc,&
!                        varId, &
!                        1,fill_value), 'nf90_def_var_fill failed for '//&
!                        vname)
                   
                   call LVT_verify(nf90_def_var_deflate(ftn_nc,&
                        varId,&
                        shuffle, deflate, deflate_level),&
                        'nf90_def_var_deflate for '//trim(vname)//&
                        'failed in LVT_LISpostMod')
                   
                endif
             endif
             
             call LVT_verify(nf90_put_att(ftn_nc,varid,"units",&
                  trim(units)),&
                  'nf90_put_att for units failed in LVT_process_LISoutput')
             call LVT_verify(nf90_put_att(ftn_nc,varid,"_FillValue",&
                  LVT_rc%udef),&
                  'nf90_put_att for _FillValue failed in LVT_process_LISoutput')
             call LVT_verify(nf90_enddef(ftn_nc))             
             
             if(LVT_LIS_rc(1)%wopt.eq."2d gridspace") then
                if(nlevels.eq.1) then
                   
                   call LVT_verify(nf90_put_var(ftn_nc,varid,gvar(:,:,1,1)),&
                        'nf90_put_var failed in LVT_LISpostMod')
                else
                   call LVT_verify(nf90_put_var(ftn_nc,varid,gvar(:,:,1,:)),&
                        'nf90_put_var failed in LVT_LISpostMod')
                endif
             elseif(LVT_LIS_rc(1)%wopt.eq."2d ensemble gridspace") then
                if(nlevels.eq.1) then
                   call LVT_verify(nf90_put_var(ftn_nc,varid,gvar(:,:,:,1)),&
                        'nf90_put_var failed in LVT_LISpostMod')
                else
                   call LVT_verify(nf90_put_var(ftn_nc,varid,gvar(:,:,:,:)),&
                        'nf90_put_var failed in LVT_LISpostMod')
                endif
             endif
          end do
          call LVT_verify(nf90_close(ftn_nc),&
               'nf90_close failed in LVT_LISpostMod')
#endif

          do i=1,LVT_LISPost%npes
             call LVT_releaseUnitNumber(ftn(i))
          enddo
       end if

    end if

  end subroutine LVT_process_LISoutput

end module LVT_LISpostMod

