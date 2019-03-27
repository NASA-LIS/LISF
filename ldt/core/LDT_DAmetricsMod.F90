!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
#include "LDT_NetCDF_inc.h"
module LDT_DAmetricsMod
!BOP
! !MODULE: LDT_DAmetricsMod
!
! !DESCRIPTION:
!  The code in this file controls the flow of various statistics computations
!
! !REVISION HISTORY:
!  02 Oct 2008: Sujay Kumar; Initial version
!
  use ESMF
  use LDT_DAmetricsDataMod
  use LDT_DAobsDataMod
  use LDT_paramDataMod

  use LDT_logMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none 
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_DAmetricsInit
  public :: LDT_diagnoseDAobsMetrics
  public :: LDT_readDAdataMask
  public :: LDT_computeDAobsMetrics
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
!EOP
 
contains

!BOP
! 
! !ROUTINE: LDT_DAmetricsInit
! \label{LDT_DAmetricsInit}
!
! !INTERFACE:   
  subroutine LDT_DAmetricsInit
! !USES: 
    use ESMF
    use map_utils
    use LDT_coreMod,       only : LDT_rc, LDT_domain
    use LDT_timeMgrMod,    only : LDT_clock, LDT_calendar, LDT_seconds2time

    implicit none
! 
! !DESCRIPTION: 
!  This routine initializes data structures required for various statistics
!  computations. 
!EOP
    integer               :: nsize
    integer               :: n 
    character*100         :: fname_domain
    integer               :: dimID(4)
    integer               :: bdimID(3)
    character(len=8)      :: date
    character(len=10)     :: time
    character(len=5)      :: zone
    integer, dimension(8) :: values
    integer               :: c,r,iret
    integer               :: lmaskid
    integer               :: xlatid, xlonid
    integer               :: xlatbid, xlonbid
    integer               :: shuffle, deflate, deflate_level    
    real,  allocatable    :: xlat(:,:),xlon(:,:)
    real,  allocatable    :: xlat_b(:,:),xlon_b(:,:)

    n = 1
    nsize = LDT_rc%ngrid(n)

    call registerMetricsEntry(LDT_DA_MOC_SWE,nsize,&
         LDT_DAobsData(n)%swe_obs,LDT_DAmetrics%swe)
    call registerMetricsEntry(LDT_DA_MOC_SNOWDEPTH, nsize,&
         LDT_DAobsData(n)%snowdepth_obs,LDT_DAmetrics%snowdepth)
    call registerMetricsEntry(LDT_DA_MOC_SOILMOIST,nsize,&
         LDT_DAobsData(n)%soilmoist_obs,LDT_DAmetrics%soilmoist)
    call registerMetricsEntry(LDT_DA_MOC_TWS,nsize,&
         LDT_DAobsData(n)%tws_obs,LDT_DAmetrics%tws)
    call registerMetricsEntry(LDT_DA_MOC_VOD,nsize,&
         LDT_DAobsData(n)%vod_obs,LDT_DAmetrics%vod)
    call registerMetricsEntry(LDT_DA_MOC_LAI,nsize,&
         LDT_DAobsData(n)%lai_obs,LDT_DAmetrics%lai)
!------------------------------------------------------------------------
! the generation of the obsgrid only doesn't require a pass through the
! data
!------------------------------------------------------------------------
    if(LDT_rc%comp_obsgrid.eq.1) then 

       allocate(xlat(LDT_rc%gnc(n),LDT_rc%gnr(n)))
       allocate(xlon(LDT_rc%gnc(n),LDT_rc%gnr(n)))
       allocate(xlat_b(LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n)))
       allocate(xlon_b(LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n)))
        
       call system('mkdir -p '//(LDT_rc%odir))
       write(LDT_logunit,*) "Writing to LDT output directory: ",&
            trim(LDT_rc%odir)
       
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
       shuffle = NETCDF_shuffle
       deflate = NETCDF_deflate
       deflate_level =NETCDF_deflate_level
       
#if (defined USE_NETCDF3)
       fname_domain = trim(LDT_rc%odir)//'/'//&
            trim(LDT_rc%dapreprocfile)//'_domain.nc'
       write(LDT_logunit,*) 'Writing CDF domain file ',trim(fname_domain)
       iret=nf90_create(path=trim(fname_domain),cmode=nf90_clobber,&
            ncid=LDT_rc%ftn_DAobs_domain)
#endif
#if (defined USE_NETCDF4)
       fname_domain = trim(LDT_rc%odir)//'/'//&
            trim(LDT_rc%dapreprocfile)//'_domain.nc'
       write(LDT_logunit,*) 'Writing CDF domain file ',trim(fname_domain)
       iret=nf90_create(path=trim(fname_domain),cmode=nf90_netcdf4,&
            ncid=LDT_rc%ftn_DAobs_domain)
#endif
      
!domain file
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 

       call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,&
            'east_west',LDT_rc%gnc(n),dimID(1)))
       call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,&
            'north_south',LDT_rc%gnr(n),dimID(2)))
       call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,'ntimes',&
            LDT_rc%cdf_ntimes,dimID(3)))
       
       call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,'east_west_b',&
            LDT_rc%gnc_buf(n),bdimID(1)))
       call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,'north_south_b',&
            LDT_rc%gnr_buf(n),bdimID(2)))
       
       if(trim(LDT_rc%lis_map_proj).eq."latlon") then !latlon
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain, &
               NF90_GLOBAL, "MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL"))
          
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "DX", &
               LDT_rc%gridDesc(n,9)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "DY", &
               LDT_rc%gridDesc(n,10)))       
          
       elseif(trim(LDT_rc%lis_map_proj).eq."mercator") then 
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "MAP_PROJECTION", &
               "MERCATOR"))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "DY", &
               LDT_rc%gridDesc(n,9)))
          
       elseif(trim(LDT_rc%lis_map_proj).eq."lambert") then !lambert conformal
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "MAP_PROJECTION", &
               "LAMBERT CONFORMAL"))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "TRUELAT2", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "DY", &
               LDT_rc%gridDesc(n,9)))
          
       elseif(trim(LDT_rc%lis_map_proj).eq."polar") then ! polar stereographic
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "MAP_PROJECTION", &
               "POLAR STEREOGRAPHIC"))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "ORIENT", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "DY", &
               LDT_rc%gridDesc(n,9)))
       elseif(trim(LDT_rc%lis_map_proj).eq."ease V2") then ! ease V2
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "MAP_PROJECTION", &
               "EASE V2"))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          if(LDT_rc%gridDesc(n,10).eq.0.09) then 
             call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
                  "GRIDTYPE", &
                  "M09"))         
          elseif(LDT_rc%gridDesc(n,10).eq.0.36) then 
             call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
                  "GRIDTYPE", &
                  "M36"))         
          endif
       endif
       
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,&
            NF90_GLOBAL,"missing_value", -9999.0))       
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "temporal_resolution_CDF", LDT_rc%cdf_ntimes))             
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,&
            NF90_GLOBAL,"title", &
            "Land Data Toolkit (LDT) output"))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,&
            NF90_GLOBAL,"institution", &
            "NASA GSFC Hydrological Sciences Laboratory"))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,"history", &
            "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
            date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,"references", &
            "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,"comment", &
            "website: http://lis.gsfc.nasa.gov/"))
       

       !landmask field
       call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
            "LANDMASK",&
            nf90_float, dimids=dimID(1:2),varid=lmaskid),&
            'nf90_def_var failed for LANDMASK')
    
#if(defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
            lmaskid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for LANDMASK')
#endif
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
            "standard_name","LANDMASK"),&
            'nf90_put_att failed for LANDMASK:standard_name')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
            "units","-"),&
            'nf90_put_att failed for LANDMASK:units')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
            "scale_factor",1.0),&
            'nf90_put_att failed for LANDMASK:scale_factor')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
            "add_offset",0.0),&
            'nf90_put_att failed for LANDMASK:add_offset')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
            "missing_value",LDT_rc%udef),&
            'nf90_put_att failed for LANDMASK:missing_value')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
            "vmin",0.0),&
            'nf90_put_att failed for LANDMASK:vmin')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
            "vmax",0.0),&
            'nf90_put_att failed for LANDMASK:vmax')
       
       !lat field
       call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
            "lat",&
            nf90_float, dimids=dimID(1:2),varid=xlatid),&
            'nf90_def_var failed for xlat')
       
#if(defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
            xlatid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for xlat')
#endif
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
            "standard_name","latitude"),&
            'nf90_put_att failed for xlat:standard_name')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
            "units","degrees_north"),&
            'nf90_put_att failed for xlat:units')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
            "scale_factor",1.0),&
            'nf90_put_att failed for xlat:scale_factor')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
            "add_offset",0.0),&
            'nf90_put_att failed for xlat:add_offset')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
            "missing_value",LDT_rc%udef),&
            'nf90_put_att failed for xlat:missing_value')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
            "vmin",0.0),&
            'nf90_put_att failed for xlat:vmin')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
            "vmax",0.0),&
            'nf90_put_att failed for xlat:vmax')
       
       !! Xlon field attributes: !!
       call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
            "lon",&
            nf90_float, dimids=dimID(1:2),varid=xlonid),&
            'nf90_def_var failed for xlon')
       
#if(defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
            xlonid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for LDT_LSMparam_struc(n)%xlon')
#endif
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
            "standard_name","longitude"),&
            'nf90_put_att failed for xlon:standard_name')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
            "units","degrees_east"),&
            'nf90_put_att failed for xlon:units')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
            "scale_factor",1.0),&
            'nf90_put_att failed for xlon:scale_factor')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
            "add_offset",0.0),&
            'nf90_put_att failed for xlon:add_offset')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
            "missing_value",LDT_rc%udef),&
            'nf90_put_att failed for xlon:missing_value')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
            "vmin",0.0),&
            'nf90_put_att failed for xlon:vmin')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
            "vmax",0.0),&
            'nf90_put_att failed for xlon:vmax')
       
       !! Xlat_b field attributes: !!
       call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
            "lat_b",&
            nf90_float, dimids=bdimID(1:2),varid=xlatbid),&
            'nf90_def_var failed for xlat_b')
       
#if(defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
            xlatbid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for xlat_b')
#endif
       
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
            "standard_name","latitude_b"),&
            'nf90_put_att failed for xlat_b:standard_name')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
            "units","degrees_north"),&
            'nf90_put_att failed for xlat_b:units')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
            "scale_factor",1.0),&
            'nf90_put_att failed for xlat_b:scale_factor')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
            "add_offset",0.0),&
            'nf90_put_att failed for xlat_b:add_offset')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
            "missing_value",LDT_rc%udef),&
            'nf90_put_att failed for xlat_b:missing_value')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
            "vmin",0.0),&
            'nf90_put_att failed for xlat_b:vmin')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
            "vmax",0.0),&
            'nf90_put_att failed for xlat_b:vmax')
       
       !! Xlon_b field attributes: !!
       call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
            "lon_b",&
            nf90_float, dimids=bdimID(1:2),varid=xlonbid),&
            'nf90_def_var failed for xlon_b')
       
#if(defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
            xlonbid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed for xlon_b')
#endif
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
            "standard_name","longitude_b"),&
            'nf90_put_att failed for xlon_b:standard_name')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
            "units","degrees_east"),&
            'nf90_put_att failed for xlon_b:units')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
            "scale_factor",1.0),&
            'nf90_put_att failed for xlon_b:scale_factor')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
            "add_offset",0.0),&
            'nf90_put_att failed for xlon_b:add_offset')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
            "missing_value",LDT_rc%udef),&
            'nf90_put_att failed for xlon_b:missing_value')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
            "vmin",0.0),&
            'nf90_put_att failed for xlon_b:vmin')
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
            "vmax",0.0),&
            'nf90_put_att failed for xlon_b:vmax')
    
#endif
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
       call LDT_verify(nf90_enddef(LDT_rc%ftn_DAobs_domain))
#endif

       do r=1,LDT_rc%gnr(n)
          do c=1,LDT_rc%gnc(n)
             call ij_to_latlon(LDT_domain(n)%ldtglbproj,&
                  real(c), real(r), xlat(c,r),&
                  xlon(c,r))
       enddo
    enddo

    do r=-1,LDT_rc%gnr(n)+2
       do c=-1,LDT_rc%gnc(n)+2
          call ij_to_latlon(LDT_domain(n)%ldtproj,&
               real(c), real(r),xlat_b(c+2,r+2),&
               xlon_b(c+2,r+2))
       enddo
    enddo

#if (defined USE_NETCDF3 || defined USE_NETCDF4)        
    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         lmaskid, LDT_LSMparam_struc(n)%landmask%value(:,:,1),&
         (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_att failed for LANDMASK')

    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         xlatid, xlat,&
         (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_att failed for xlat')
    
    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         xlonid, xlon,&
         (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_att failed for xlon')
    
    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         xlatbid, xlat_b,&
         (/1,1/),(/LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n)/)),&
         'nf90_put_att failed for xlat_b')
    
    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         xlonbid, xlon_b,&
         (/1,1/),(/LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n)/)),&
         'nf90_put_att failed for xlon_b')
    
#endif
       iret=nf90_close(LDT_rc%ftn_DAobs_domain)
       write(LDT_logunit,*) 'Successfully wrote CDF file ',trim(fname_domain)
#endif
       
       deallocate(xlat)
       deallocate(xlon)
       deallocate(xlat_b)
       deallocate(xlon_b)
    endif
    
  end subroutine LDT_DAmetricsInit

!BOP
! !ROUTINE: registerMetricsEntry
! \label{registerMetricsEntry}
! 
! !INTERFACE: 
  subroutine registerMetricsEntry(ldt_moc_index,nsize, obs, metrics)
! !USES: 
    use LDT_coreMod,  only : LDT_rc

! !ARGUMENTS:     
    integer                 :: ldt_moc_index
    integer                 :: nsize
    type(LDT_DAmetaDataEntry) :: obs
    type(DAmetricsEntry), target :: metrics
! 
! !DESCRIPTION: 
!  This routine initializes the objects to hold different statistics 
!  computations
!
!   The arguments are: 
!   \begin{description}
!    \item[nsize] number of obsing points 
!    \item[obs] object to hold obs variable information
!    \item[metrics] object to hold statistics computations
!   \end{description}
!EOP
    
    LDT_DAmetricsPtr(ldt_moc_index)%dataEntryPtr => metrics
    call initMetricsEntry(nsize, obs, metrics)

  end subroutine registerMetricsEntry

!BOP
! !ROUTINE: initMetricsEntry
! \label{initMetricsEntry}
! 
! !INTERFACE: 
  subroutine initMetricsEntry(nsize, obs, metrics)
! !USES: 
    use LDT_coreMod,  only : LDT_rc

! !ARGUMENTS:     
    integer                 :: nsize
    type(LDT_DAmetaDataEntry) :: obs
    type(DAmetricsEntry) :: metrics
! 
! !DESCRIPTION: 
!  This routine initializes the objects to hold different statistics 
!  computations
!
!   The arguments are: 
!   \begin{description}
!    \item[nsize] number of obsing points 
!    \item[obs] object to hold obs variable information
!    \item[metrics] object to hold statistics computations
!   \end{description}
!EOP

    metrics%standard_name = obs%standard_name
    metrics%selectOpt = obs%selectStats

    if(metrics%selectOpt.eq.1) then 
       if(LDT_rc%comp_cdf.eq.1) then 
          allocate(metrics%mask(nsize,LDT_rc%cdf_ntimes, obs%vlevels))
          allocate(metrics%maxval(nsize, LDT_rc%cdf_ntimes,obs%vlevels))
          allocate(metrics%minval(nsize, LDT_rc%cdf_ntimes,obs%vlevels))
          metrics%mask = LDT_rc%udef
          metrics%maxval = -1000000.0
          metrics%minval =  1000000.0
          allocate(metrics%count_drange_total(nsize, LDT_rc%cdf_ntimes,&
               obs%vlevels))
          metrics%count_drange_total = 0

          allocate(metrics%cdf_bincounts(nsize,&
               LDT_rc%cdf_ntimes, obs%vlevels,LDT_rc%cdf_nbins))
          allocate(metrics%delta(nsize,LDT_rc%cdf_ntimes,obs%vlevels))
          allocate(metrics%xrange(nsize,&
               LDT_rc%cdf_ntimes, obs%vlevels,LDT_rc%cdf_nbins))
          allocate(metrics%cdf(nsize,&
               LDT_rc%cdf_ntimes, obs%vlevels,LDT_rc%cdf_nbins))
          
          metrics%cdf_bincounts = 0 
          metrics%delta = 0 
          metrics%xrange = 0 
          metrics%cdf = 0 

          allocate(metrics%sx_mu(nsize, LDT_rc%cdf_ntimes, obs%vlevels))
          allocate(metrics%mu(nsize, LDT_rc%cdf_ntimes, obs%vlevels))
          allocate(metrics%count_mu(nsize, LDT_rc%cdf_ntimes,obs%vlevels))

          metrics%sx_mu = 0 
          metrics%mu = 0 
          metrics%count_mu = 0 

          allocate(metrics%sx_sigma(nsize, LDT_rc%cdf_ntimes,obs%vlevels))
          allocate(metrics%sxx_sigma(nsize, LDT_rc%cdf_ntimes,obs%vlevels))
          allocate(metrics%sigma(nsize, LDT_rc%cdf_ntimes,obs%vlevels))
          allocate(metrics%count_sigma(nsize, LDT_rc%cdf_ntimes,obs%vlevels))

          metrics%sx_sigma = 0 
          metrics%sxx_sigma = 0 
          metrics%sigma = 0 
          metrics%count_sigma = 0           

       endif
    endif

  end subroutine initMetricsEntry

!BOP
! 
! !ROUTINE: LDT_diagnoseDAobsMetrics
! \label{LDT_diagnoseDAobsMetrics}
! 
! !INTERFACE:
  subroutine LDT_diagnoseDAobsMetrics(n, pass)
! !USES: 
    use LDT_coreMod,   only : LDT_rc
    use LDT_DrangeMod,    only : LDT_diagnoseDrange
    use LDT_MuMod,    only : LDT_diagnoseMu
    use LDT_SigmaMod,    only : LDT_diagnoseSigma
    use LDT_CDFMod,       only : LDT_diagnoseCDF
! 
! !DESCRIPTION: 
! This routine invokes the methods for computing the specified statistics
! 
!EOP

    implicit none

    integer       :: pass
    integer       :: n 

    if(LDT_rc%computeFlag) then 
       if(LDT_rc%comp_cdf.eq.1) then           
          if(pass.eq.1) then 
             call LDT_diagnoseDrange(n)
             call LDT_diagnoseMu(n)
             call LDT_diagnoseSigma(n)
          endif
          if(pass.eq.2) then 
             call LDT_diagnoseCDF(n)
          endif
       endif
       
    endif

  end subroutine LDT_diagnoseDAobsMetrics


!BOP
! 
! !ROUTINE: LDT_computeDAobsMetrics
! \label{LDT_computeDAobsMetrics}
! 
! !INTERFACE:
  subroutine LDT_computeDAobsMetrics(n, pass)
! !USES: 
    use LDT_coreMod,   only : LDT_rc
    use LDT_DrangeMod,    only : LDT_computeDrange
    use LDT_MuMod,    only : LDT_computeMu
    use LDT_SigmaMod,    only : LDT_computeSigma
    use LDT_CDFMod,       only : LDT_computeCDF

! 
! !DESCRIPTION: 
! This routine invokes the methods for computing the specified statistics
! 
!EOP

    implicit none

    character*100 :: fname_cdf
    character*100 :: fname_domain
    integer       :: pass
    integer       :: rc
    integer       :: n 
    integer       :: iret


    if(LDT_rc%computeFlag) then 
       if(LDT_rc%comp_cdf.eq.1) then           
          if(pass.eq.1) then 
             call LDT_computeDrange(n)          
             call LDT_computeMu(n)
             call LDT_computeSigma(n)
          endif
          if(pass.eq.2) then 
             call LDT_computeCDF(n)
          endif
       endif
       
       if(pass.eq.2.and.LDT_rc%comp_cdf.eq.1) then 
          
          if(LDT_rc%endtime.eq.1) then 
             call system('mkdir -p '//(LDT_rc%odir))
             write(LDT_logunit,*) "Writing to LDT output directory: ",&
                   trim(LDT_rc%odir)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
#if (defined USE_NETCDF3)
             if(LDT_rc%comp_cdf.eq.1) then 
                fname_cdf = trim(LDT_rc%odir)//'/'//&
                     trim(LDT_rc%dapreprocfile)//'.nc'
                write(LDT_logunit,*) 'Writing CDF file ',trim(fname_cdf)
                iret=nf90_create(path=trim(fname_cdf),cmode=nf90_clobber,&
                     ncid=LDT_rc%ftn_cdf)
                fname_domain = trim(LDT_rc%odir)//'/'//&
                     trim(LDT_rc%dapreprocfile)//'_domain.nc'
                write(LDT_logunit,*) 'Writing CDF domain file ',trim(fname_domain)
                iret=nf90_create(path=trim(fname_domain),cmode=nf90_clobber,&
                     ncid=LDT_rc%ftn_DAobs_domain)
             endif
#endif
#if (defined USE_NETCDF4)
             if(LDT_rc%comp_cdf.eq.1) then 
                fname_cdf = trim(LDT_rc%odir)//'/'//&
                     trim(LDT_rc%dapreprocfile)//'.nc'
                write(LDT_logunit,*) 'Writing CDF file ',trim(fname_cdf)
                iret=nf90_create(path=trim(fname_cdf),cmode=nf90_netcdf4,&
                     ncid=LDT_rc%ftn_cdf)
                fname_domain = trim(LDT_rc%odir)//'/'//&
                     trim(LDT_rc%dapreprocfile)//'_domain.nc'
                write(LDT_logunit,*) 'Writing CDF domain file ',trim(fname_domain)
                iret=nf90_create(path=trim(fname_domain),cmode=nf90_netcdf4,&
                     ncid=LDT_rc%ftn_DAobs_domain)
             endif
#endif
             
             call outputFinalMetrics(n,pass)
             
             if(LDT_rc%comp_cdf.eq.1) then 
                iret=nf90_close(LDT_rc%ftn_cdf)
                write(LDT_logunit,*) 'Successfully wrote CDF file ',trim(fname_cdf)
                iret=nf90_close(LDT_rc%ftn_DAobs_domain)
                write(LDT_logunit,*) 'Successfully wrote CDF file ',trim(fname_domain)
             endif
             
#endif
          endif
       endif

       LDT_rc%computeFlag = .false. 
       
       call LDT_resetDAobsData(n)

    endif

  end subroutine LDT_computeDAobsMetrics

!BOP
! !ROUTINE: outputFinalMetrics
! \label{outputFinalMetrics}
! 
! !INTERFACE: 
  subroutine outputFinalMetrics(n,pass)
    use LDT_coreMod
    use map_utils

    implicit none
! !ARGUMENTS: 
    integer          :: n
    integer          :: pass
! 
! !DESCRIPTION: 
!  This routine writes the final set of statistics at the end of the analysis. 
!
!EOP
    integer               :: index
    integer               :: c,r
    integer               :: dimID(4)
    integer               :: bdimID(3)
    character(len=8)      :: date
    character(len=10)     :: time
    character(len=5)      :: zone
    integer, dimension(8) :: values
    integer               :: lmaskid
    integer               :: xlatid, xlonid
    integer               :: xlatbid, xlonbid
    integer               :: shuffle, deflate, deflate_level    
    real                  :: xlat(LDT_rc%gnc(n),LDT_rc%gnr(n))
    real                  :: xlon(LDT_rc%gnc(n),LDT_rc%gnr(n))
    real                  :: xlat_b(LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n))
    real                  :: xlon_b(LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n))

!global headers
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    call date_and_time(date,time,zone,values)

    call LDT_verify(nf90_def_dim(LDT_rc%ftn_cdf,'ngrid',&
         LDT_rc%glbngrid(n),dimID(1)))
    call LDT_verify(nf90_def_dim(LDT_rc%ftn_cdf,'ntimes',&
         LDT_rc%cdf_ntimes,dimID(2)))
    call LDT_verify(nf90_def_dim(LDT_rc%ftn_cdf,'nbins',&
         LDT_rc%cdf_nbins,dimID(4)))

    call LDT_verify(nf90_put_att(LDT_rc%ftn_cdf,NF90_GLOBAL,&
         "missing_value", -9999.0))          
    call LDT_verify(nf90_put_att(LDT_rc%ftn_cdf,NF90_GLOBAL,&
         "temporal_resolution_CDF", LDT_rc%cdf_ntimes))          
    call LDT_verify(nf90_put_att(LDT_rc%ftn_cdf,NF90_GLOBAL,&
         "title", &
         "Land Data Toolkit (LDT) output"))
    call LDT_verify(nf90_put_att(LDT_rc%ftn_cdf,NF90_GLOBAL,&
         "institution", &
         "NASA GSFC Hydrological Sciences Laboratory"))
    call LDT_verify(nf90_put_att(LDT_rc%ftn_cdf,NF90_GLOBAL,&
         "history", &
         "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
         date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)))
    call LDT_verify(nf90_put_att(LDT_rc%ftn_cdf,NF90_GLOBAL,"references", &
         "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"))
    call LDT_verify(nf90_put_att(LDT_rc%ftn_cdf,NF90_GLOBAL,"comment", &
         "website: http://lis.gsfc.nasa.gov/"))
#endif
    do index=1,LDT_DA_MOC_COUNT
       call writeFinalSingleEntryHeader(LDT_DAmetricsPtr(index)%dataEntryPtr,&
            LDT_DAobsDataPtr(1,index)%dataEntryPtr, dimID)
    enddo

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    call LDT_verify(nf90_enddef(LDT_rc%ftn_cdf))
#endif

    do index=1,LDT_DA_MOC_COUNT
       call writeFinalSingleEntry(pass,LDT_DAmetricsPtr(index)%dataEntryPtr,&
            LDT_DAobsDataPtr(1,index)%dataEntryPtr)
    enddo

!domain file
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 

    call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,&
         'east_west',LDT_rc%gnc(n),dimID(1)))
    call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,&
         'north_south',LDT_rc%gnr(n),dimID(2)))
    call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,'ntimes',&
         LDT_rc%cdf_ntimes,dimID(3)))
    
    call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,'east_west_b',&
         LDT_rc%gnc_buf(n),bdimID(1)))
    call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,'north_south_b',&
         LDT_rc%gnr_buf(n),bdimID(2)))
    
    if(trim(LDT_rc%lis_map_proj).eq."latlon") then !latlon
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain, &
            NF90_GLOBAL, "MAP_PROJECTION", "EQUIDISTANT CYLINDRICAL"))

       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LAT", &
            LDT_rc%gridDesc(n,4)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LON", &
            LDT_rc%gridDesc(n,5)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "DX", &
            LDT_rc%gridDesc(n,9)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "DY", &
            LDT_rc%gridDesc(n,10)))       
       
    elseif(trim(LDT_rc%lis_map_proj).eq."mercator") then 
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "MAP_PROJECTION", &
            "MERCATOR"))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LAT", &
            LDT_rc%gridDesc(n,4)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LON", &
            LDT_rc%gridDesc(n,5)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "TRUELAT1", &
            LDT_rc%gridDesc(n,10)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "STANDARD_LON", &
            LDT_rc%gridDesc(n,11)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "DX", &
            LDT_rc%gridDesc(n,8)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "DY", &
            LDT_rc%gridDesc(n,9)))
       
    elseif(trim(LDT_rc%lis_map_proj).eq."lambert") then !lambert conformal
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "MAP_PROJECTION", &
            "LAMBERT CONFORMAL"))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LAT", &
            LDT_rc%gridDesc(n,4)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LON", &
            LDT_rc%gridDesc(n,5)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "TRUELAT1", &
            LDT_rc%gridDesc(n,10)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "TRUELAT2", &
            LDT_rc%gridDesc(n,7)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "STANDARD_LON", &
            LDT_rc%gridDesc(n,11)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "DX", &
            LDT_rc%gridDesc(n,8)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "DY", &
            LDT_rc%gridDesc(n,9)))
       
    elseif(trim(LDT_rc%lis_map_proj).eq."polar") then ! polar stereographic
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "MAP_PROJECTION", &
            "POLAR STEREOGRAPHIC"))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LAT", &
            LDT_rc%gridDesc(n,4)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LON", &
            LDT_rc%gridDesc(n,5)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "TRUELAT1", &
            LDT_rc%gridDesc(n,10)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "ORIENT", &
            LDT_rc%gridDesc(n,7)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "STANDARD_LON", &
            LDT_rc%gridDesc(n,11)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "DX", &
            LDT_rc%gridDesc(n,8)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "DY", &
            LDT_rc%gridDesc(n,9)))
    elseif(trim(LDT_rc%lis_map_proj).eq."ease V2") then ! ease V2
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "MAP_PROJECTION", &
            "EASE V2"))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LAT", &
            LDT_rc%gridDesc(n,4)))
       call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
            "SOUTH_WEST_CORNER_LON", &
            LDT_rc%gridDesc(n,5)))
       if(LDT_rc%gridDesc(n,10).eq.0.09) then 
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "GRIDTYPE", &
               "M09"))         
       elseif(LDT_rc%gridDesc(n,10).eq.0.36) then 
          call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
               "GRIDTYPE", &
               "M36"))         
       endif
    endif
 
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,&
         NF90_GLOBAL,"missing_value", -9999.0))       
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,&
         "temporal_resolution_CDF", LDT_rc%cdf_ntimes))             
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,&
         NF90_GLOBAL,"title", &
         "Land Data Toolkit (LDT) output"))
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,&
         NF90_GLOBAL,"institution", &
         "NASA GSFC Hydrological Sciences Laboratory"))
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,"history", &
         "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
         date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)))
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,"references", &
         "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"))
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,NF90_GLOBAL,"comment", &
         "website: http://lis.gsfc.nasa.gov/"))


!landmask field
    call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
         "LANDMASK",&
         nf90_float, dimids=dimID(1:2),varid=lmaskid),&
         'nf90_def_var failed for LANDMASK')
    
#if(defined USE_NETCDF4) 
    call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
         lmaskid, shuffle, deflate, deflate_level),&
         'nf90_def_var_deflate failed for LANDMASK')
#endif
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
         "standard_name","LANDMASK"),&
         'nf90_put_att failed for LANDMASK:standard_name')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
         "units","-"),&
         'nf90_put_att failed for LANDMASK:units')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
         "scale_factor",1.0),&
         'nf90_put_att failed for LANDMASK:scale_factor')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
         "add_offset",0.0),&
         'nf90_put_att failed for LANDMASK:add_offset')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
         "missing_value",LDT_rc%udef),&
         'nf90_put_att failed for LANDMASK:missing_value')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
         "vmin",0.0),&
         'nf90_put_att failed for LANDMASK:vmin')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,lmaskid, &
         "vmax",0.0),&
         'nf90_put_att failed for LANDMASK:vmax')

!lat field
    call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
         "lat",&
         nf90_float, dimids=dimID(1:2),varid=xlatid),&
         'nf90_def_var failed for xlat')
    
#if(defined USE_NETCDF4) 
    call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
         xlatid, shuffle, deflate, deflate_level),&
         'nf90_def_var_deflate failed for xlat')
#endif
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
         "standard_name","latitude"),&
         'nf90_put_att failed for xlat:standard_name')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
         "units","degrees_north"),&
         'nf90_put_att failed for xlat:units')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
         "scale_factor",1.0),&
         'nf90_put_att failed for xlat:scale_factor')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
         "add_offset",0.0),&
         'nf90_put_att failed for xlat:add_offset')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
         "missing_value",LDT_rc%udef),&
         'nf90_put_att failed for xlat:missing_value')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
         "vmin",0.0),&
         'nf90_put_att failed for xlat:vmin')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatid, &
         "vmax",0.0),&
         'nf90_put_att failed for xlat:vmax')
    
    !! Xlon field attributes: !!
    call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
         "lon",&
         nf90_float, dimids=dimID(1:2),varid=xlonid),&
         'nf90_def_var failed for xlon')
    
#if(defined USE_NETCDF4) 
    call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
         xlonid, shuffle, deflate, deflate_level),&
         'nf90_def_var_deflate failed for LDT_LSMparam_struc(n)%xlon')
#endif
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
         "standard_name","longitude"),&
         'nf90_put_att failed for xlon:standard_name')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
         "units","degrees_east"),&
         'nf90_put_att failed for xlon:units')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
         "scale_factor",1.0),&
         'nf90_put_att failed for xlon:scale_factor')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
         "add_offset",0.0),&
         'nf90_put_att failed for xlon:add_offset')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
         "missing_value",LDT_rc%udef),&
         'nf90_put_att failed for xlon:missing_value')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
         "vmin",0.0),&
         'nf90_put_att failed for xlon:vmin')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonid, &
         "vmax",0.0),&
         'nf90_put_att failed for xlon:vmax')
    
    !! Xlat_b field attributes: !!
    call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
         "lat_b",&
         nf90_float, dimids=bdimID(1:2),varid=xlatbid),&
         'nf90_def_var failed for xlat_b')
    
#if(defined USE_NETCDF4) 
    call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
         xlatbid, shuffle, deflate, deflate_level),&
         'nf90_def_var_deflate failed for xlat_b')
#endif
    
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
         "standard_name","latitude_b"),&
         'nf90_put_att failed for xlat_b:standard_name')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
         "units","degrees_north"),&
         'nf90_put_att failed for xlat_b:units')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
         "scale_factor",1.0),&
         'nf90_put_att failed for xlat_b:scale_factor')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
         "add_offset",0.0),&
         'nf90_put_att failed for xlat_b:add_offset')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
         "missing_value",LDT_rc%udef),&
         'nf90_put_att failed for xlat_b:missing_value')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
         "vmin",0.0),&
         'nf90_put_att failed for xlat_b:vmin')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlatbid, &
         "vmax",0.0),&
         'nf90_put_att failed for xlat_b:vmax')
    
    !! Xlon_b field attributes: !!
    call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
         "lon_b",&
         nf90_float, dimids=bdimID(1:2),varid=xlonbid),&
         'nf90_def_var failed for xlon_b')
    
#if(defined USE_NETCDF4) 
    call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
         xlonbid, shuffle, deflate, deflate_level),&
         'nf90_def_var_deflate failed for xlon_b')
#endif
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
         "standard_name","longitude_b"),&
         'nf90_put_att failed for xlon_b:standard_name')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
         "units","degrees_east"),&
         'nf90_put_att failed for xlon_b:units')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
         "scale_factor",1.0),&
         'nf90_put_att failed for xlon_b:scale_factor')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
         "add_offset",0.0),&
         'nf90_put_att failed for xlon_b:add_offset')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
         "missing_value",LDT_rc%udef),&
         'nf90_put_att failed for xlon_b:missing_value')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
         "vmin",0.0),&
         'nf90_put_att failed for xlon_b:vmin')
    call LDT_verify(nf90_put_att(LDT_rc%ftn_DAobs_domain,xlonbid, &
         "vmax",0.0),&
         'nf90_put_att failed for xlon_b:vmax')
    
#endif
    do index=1,LDT_DA_MOC_COUNT
       call writeFinalSingleDomainEntryHeader(&
            LDT_DAmetricsPtr(index)%dataEntryPtr,&
            LDT_DAobsDataPtr(1,index)%dataEntryPtr, dimID)
    enddo

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
    call LDT_verify(nf90_enddef(LDT_rc%ftn_DAobs_domain))
#endif

    do r=1,LDT_rc%gnr(n)
       do c=1,LDT_rc%gnc(n)
          call ij_to_latlon(LDT_domain(n)%ldtglbproj,&
               real(c), real(r), xlat(c,r),&
               xlon(c,r))
       enddo
    enddo

    do r=-1,LDT_rc%gnr(n)+2
       do c=-1,LDT_rc%gnc(n)+2
          call ij_to_latlon(LDT_domain(n)%ldtproj,&
               real(c), real(r),xlat_b(c+2,r+2),&
               xlon_b(c+2,r+2))
       enddo
    enddo

#if (defined USE_NETCDF3 || defined USE_NETCDF4)        
    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         lmaskid, LDT_LSMparam_struc(n)%landmask%value(:,:,1),&
         (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_att failed for LANDMASK')

    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         xlatid, xlat,&
         (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_att failed for xlat')
    
    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         xlonid, xlon,&
         (/1,1/),(/LDT_rc%gnc(n),LDT_rc%gnr(n)/)),&
         'nf90_put_att failed for xlon')
    
    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         xlatbid, xlat_b,&
         (/1,1/),(/LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n)/)),&
         'nf90_put_att failed for xlat_b')
    
    call LDT_verify(nf90_put_var(LDT_rc%ftn_DAobs_domain,&
         xlonbid, xlon_b,&
         (/1,1/),(/LDT_rc%gnc_buf(n),LDT_rc%gnr_buf(n)/)),&
         'nf90_put_att failed for xlon_b')
    
#endif
    do index=1,LDT_DA_MOC_COUNT
       call writeFinalSingleDomainEntry(pass,&
            LDT_DAmetricsPtr(index)%dataEntryPtr,&
            LDT_DAobsDataPtr(1,index)%dataEntryPtr)
    enddo

  end subroutine OutputFinalMetrics


!BOP
! !ROUTINE: writeFinalSingleEntryHeader
! \label{writeFinalSingleEntryHeader}
! 
! !INTERFACE:   
  subroutine writeFinalSingleEntryHeader(metrics, obs, dimID)
! !USES: 
    use LDT_coreMod,    only  : LDT_rc
    use LDT_historyMod, only  : LDT_writevar_gridded

    implicit none
! !ARGUMENTS: 
    type(DAmetricsEntry)      :: metrics
    type(LDT_DAmetadataEntry) :: obs
    integer                 :: dimID(4)
!
! !DESCRIPTION: 
!  This routine writes the specified set of statistics for a 
!  single variable at the end of the analysis. 
!EOP

    integer    :: k 
    integer    :: varid1, varid2
    integer    :: i,c,r
    integer    :: n 
    character*100 :: vname
    integer :: shuffle, deflate, deflate_level

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    n = 1
    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
       vname = trim(obs%standard_name)//'_levels'
       call LDT_verify(nf90_def_dim(LDT_rc%ftn_cdf,trim(vname),&
            obs%vlevels,dimID(3)))
       

       call LDT_verify(nf90_def_var(LDT_rc%ftn_cdf,&
            trim(obs%standard_name)//'_xrange',&
            nf90_float, dimids = dimID, varid=obs%varID(1)))

#if (defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_cdf,&
            obs%varID(1), shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed in LDT_DAmetricsMod')
#endif

       call LDT_verify(nf90_def_var(LDT_rc%ftn_cdf,&
            trim(obs%standard_name)//'_mu',&
            nf90_float, dimids = dimID(1:3), varid=obs%varID(2)))

#if (defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_cdf,&
            obs%varID(2), shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed in LDT_DAmetricsMod')
#endif
       call LDT_verify(nf90_def_var(LDT_rc%ftn_cdf,&
            trim(obs%standard_name)//'_sigma',&
            nf90_float, dimids = dimID(1:3), varid=obs%varID(3)))

#if (defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_cdf,&
            obs%varID(3), shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed in LDT_DAmetricsMod')
#endif
       call LDT_verify(nf90_def_var(LDT_rc%ftn_cdf,&
            trim(obs%standard_name)//'_CDF',&
            nf90_float, dimids = dimID, varid=obs%varID(4)))

#if (defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_cdf,&
            obs%varID(4), shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed in LDT_DAmetricsMod')
#endif

!       call LDT_verify(nf90_def_var(LDT_rc%ftn_cdf,&
!            trim(obs%standard_name)//'_COUNTS',&
!            nf90_float, dimids = dimID, varid=obs%varID(5)))

!#if (defined USE_NETCDF4) 
!       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_cdf,&
!            obs%varID(5), shuffle, deflate, deflate_level),&
!            'nf90_def_var_deflate failed in LDT_DAmetricsMod')
!#endif

    endif
#endif
    
  end subroutine writeFinalSingleEntryHeader


!BOP
! !ROUTINE: writeFinalSingleEntry
! \label{writeFinalSingleEntry}
! 
! !INTERFACE:   
  subroutine writeFinalSingleEntry(pass,metrics, obs)
! !USES: 
    use LDT_coreMod,    only  : LDT_rc
    use LDT_historyMod, only  : LDT_writevar_gridded

    implicit none
! !ARGUMENTS: 
    integer                 :: pass
    type(DAmetricsEntry)      :: metrics
    type(LDT_DAmetadataEntry) :: obs
!
! !DESCRIPTION: 
!  This routine writes the specified set of statistics for a 
!  single variable at the end of the analysis. 
!EOP

    integer    :: k 
    integer    :: varid1, varid2
    integer    :: i,c,r
    integer    :: n 

    n = 1
    if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
       
       if(pass.eq.2.and.LDT_rc%comp_cdf.eq.1) then 
          call LDT_writevar_gridded(n,LDT_rc%ftn_cdf,&
               metrics%xrange, LDT_rc%cdf_ntimes, &
               obs%vlevels, LDT_rc%cdf_nbins,&
               obs%varID(1))
          call LDT_writevar_gridded(n,LDT_rc%ftn_cdf,&
               metrics%mu, LDT_rc%cdf_ntimes, obs%vlevels, &
               obs%varID(2))
          call LDT_writevar_gridded(n,LDT_rc%ftn_cdf,&
               metrics%sigma, LDT_rc%cdf_ntimes,obs%vlevels,  &
               obs%varID(3))
          call LDT_writevar_gridded(n,LDT_rc%ftn_cdf,&
               metrics%cdf, LDT_rc%cdf_ntimes, &
               obs%vlevels, LDT_rc%cdf_nbins,&
               obs%varID(4))
!          call LDT_writevar_gridded(n,LDT_rc%ftn_cdf,&
!               float(metrics%cdf_bincounts), LDT_rc%cdf_ntimes, &
!               obs%vlevels, LDT_rc%cdf_nbins,&
!               obs%varID(5))
       endif
    endif

  end subroutine writeFinalSingleEntry

!BOP
! !ROUTINE: writeFinalSingleDomainEntryHeader
! \label{writeFinalSingleDomainEntryHeader}
! 
! !INTERFACE:   
  subroutine writeFinalSingleDomainEntryHeader(metrics, obs, dimID)
! !USES: 
    use LDT_coreMod,    only  : LDT_rc
    use LDT_historyMod, only  : LDT_writevar_gridded

    implicit none
! !ARGUMENTS: 
    type(DAmetricsEntry)      :: metrics
    type(LDT_DAmetadataEntry) :: obs
    integer                 :: dimID(4)
!
! !DESCRIPTION: 
!  This routine writes the specified set of statistics for a 
!  single variable at the end of the analysis. 
!EOP

    integer    :: k 
    integer    :: varid1, varid2
    integer    :: i,c,r
    integer    :: n 
    character*100 :: vname
    integer :: shuffle, deflate, deflate_level

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    n = 1

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
       vname = trim(obs%standard_name)//'_levels'
       call LDT_verify(nf90_def_dim(LDT_rc%ftn_DAobs_domain,trim(vname),&
            obs%vlevels,dimID(4)))

       call LDT_verify(nf90_def_var(LDT_rc%ftn_DAobs_domain,&
            trim(obs%standard_name)//'_domain',&
            nf90_float, dimids = dimID, varid=obs%varID(1)))

#if (defined USE_NETCDF4) 
       call LDT_verify(nf90_def_var_deflate(LDT_rc%ftn_DAobs_domain,&
            obs%varID(1), shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate failed in LDT_DAmetricsMod')
#endif
    endif
#endif
    
  end subroutine writeFinalSingleDomainEntryHeader


!BOP
! !ROUTINE: writeFinalSingleDomainEntry
! \label{writeFinalSingleDomainEntry}
! 
! !INTERFACE:   
  subroutine writeFinalSingleDomainEntry(pass,metrics, obs)
! !USES: 
    use LDT_coreMod,    only  : LDT_rc
    use LDT_historyMod, only  : LDT_writevar_gridded

    implicit none
! !ARGUMENTS: 
    integer                 :: pass
    type(DAmetricsEntry)      :: metrics
    type(LDT_DAmetadataEntry) :: obs
!
! !DESCRIPTION: 
!  This routine writes the specified set of statistics for a 
!  single variable at the end of the analysis. 
!EOP

    integer    :: k 
    integer    :: varid1, varid2
    integer    :: i,j,c,r
    integer    :: n 

    n = 1
    if(obs%selectOpt.eq.1.and.metrics%selectOpt.eq.1) then 
       
       if(pass.eq.2.and.LDT_rc%comp_cdf.eq.1) then 
          do k=1,obs%vlevels
             do j=1,LDT_rc%cdf_ntimes
                call LDT_writevar_gridded(n,LDT_rc%ftn_DAobs_domain,&
                     metrics%mask(:,j,k), obs%varID(1), dim1=k,&
                     dim2=j,wopt = "2d gridspace" )
             enddo
          enddo
       endif
    endif

  end subroutine writeFinalSingleDomainEntry

!BOP
! !ROUTINE: LDT_readDAdataMask
! \label{LDT_readDAdataMask}
! 
! !INTERFACE: 
  subroutine LDT_readDAdataMask(n)

! !USES:     
    use LDT_coreMod,    only : LDT_rc
    use LDT_historyMod, only : LDT_readvar_gridded

    implicit none
! 
! !DESCRIPTION: 
!   This routine reads the external data mask, which is used to 
!   screen grid points, both spatially and temporally. 
!EOP
    integer       :: n 
    character*100 :: maskfile
    logical       :: file_exists
    real          :: datamask(LDT_rc%lnc(n), LDT_rc%lnr(n))
    integer       :: ftn
    integer       :: c,r


    if(LDT_rc%applyMask.eq.1) then  !read the external mask, timestamped file
       call create_datamask_filename(maskfile)
       inquire(file=trim(maskfile), exist=file_exists)
       if(file_exists) then 
          write(LDT_logunit,*) 'Reading mask ',trim(maskfile)
          ftn = LDT_getNextUnitNumber()
          open(ftn,file=trim(maskfile), form='unformatted')

          call LDT_readvar_gridded(n,ftn,datamask)

          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(datamask(c,r).ne.LDT_rc%udef) then 
                   LDT_rc%datamask(c,r) = 1
                else
                   LDT_rc%datamask(c,r) = 0
                endif
             enddo
          enddo
          call LDT_releaseUnitNumber(ftn)          
       else
          LDT_rc%datamask = 0
       endif

    elseif(LDT_rc%applyMask.eq.2) then  !read the static mask file
       inquire(file=trim(LDT_rc%maskdir), exist=file_exists)

       if(file_exists) then 
          write(LDT_logunit,*) 'Reading mask ',trim(LDT_rc%maskdir)
          ftn = LDT_getNextUnitNumber()
          open(ftn,file=trim(LDT_rc%maskdir), form='unformatted')
          
          call LDT_readvar_gridded(n,ftn,datamask)
          
          do r=1,LDT_rc%lnr(n)
             do c=1,LDT_rc%lnc(n)
                if(datamask(c,r).ne.LDT_rc%udef) then 
                   LDT_rc%datamask(c,r) = 1
                else
                   LDT_rc%datamask(c,r) = 0
                endif
             enddo
          enddo
          call LDT_releaseUnitNumber(ftn)          
       else
          LDT_rc%datamask = 0 
       endif
    else
       LDT_rc%datamask = 1
    endif
    
  end subroutine LDT_readDAdataMask


!BOP
! !ROUTINE: create_datamask_filename
! \label{create_datamask_filename}
! 
! !INTERFACE: 
  subroutine create_datamask_filename(maskfile)
! !USES: 
    use LDT_coreMod, only : LDT_rc

    implicit none
! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped filename for the external 
!  mask data
!EOP
    character(len=*)    :: maskfile

    character*12        :: fdate

    write(unit=fdate, fmt='(i4.4,i2.2,i2.2,i2.2,i2.2)') LDT_rc%yr, &
         LDT_rc%mo, LDT_rc%da, LDT_rc%hr, LDT_rc%mn
!    maskfile = trim(LDT_rc%maskdir)//&
!         '/'//trim(fdate)//'.1gs4r'
    maskfile = trim(LDT_rc%maskdir)//&
         '/'//trim(fdate)//'.d01.gs4r'

  end subroutine create_datamask_filename

  
end module LDT_DAmetricsMod
