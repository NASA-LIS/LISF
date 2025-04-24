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
#include "LIS_NetCDF_inc.h"

!BOP
! !ROUTINE: RAPID_routing_writerst
! \label{RAPID_routing_writerst}
!
! !REVISION HISTORY:
! 13 Jul 2021: Yeosang Yoon;  Initial implementation
! 31 Aug 2022: Yeosang Yoon;  fix code to product rst file on time
! 20 Mar 2024: Yeosang Yoon; Support to run with ensemble mode

subroutine RAPID_routing_writerst(n)

!
! !DESCRIPTION:
!  This routine writes restart files for RAPID. The restart files
!  are in NetCDF format.
!

  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_timeMgrMod
  use RAPID_routingMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
  
  integer, intent(in)   :: n 
  
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer               :: ftn
  integer               :: status
  logical               :: alarmCheck
  integer               :: dimid_time, dimid_riv_bas, dim_Qout(2)
  integer               :: dimid_ens, dimid_Qout_ens(3)            !2d ensemble gridspace
  integer               :: varid_Qout
  integer               :: shuffle, deflate, deflate_level

  integer, dimension(8) :: values
  character(len=8)      :: date
  character(len=10)     :: time
  character(len=5)      :: zone

  integer :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/
  integer :: mo, mn, ss

  ! check leap year
  if((mod(LIS_rc%yr,4) .eq. 0 .and. mod(LIS_rc%yr, 100).ne.0) &
       .or.(mod(LIS_rc%yr,400) .eq.0)) then
     days(2) = 29
  else
     days(2) = 28
  endif

  alarmCheck = .false.
  ! "30-day" month interval based alarms:
  if(RAPID_routing_struc(n)%rstInterval == 2592000.0) then
     if (LIS_rc%da == days(LIS_rc%mo) .and. (LIS_rc%hr == 23) &
         .and. (LIS_rc%mn == 45)) then
        alarmCheck = .true.
     endif
  else ! daily interval
     alarmCheck = LIS_isAlarmRinging(LIS_rc,&
         "RAPID router restart alarm")
  endif

  ! for NetCDF deflate compression option
  shuffle = NETCDF_shuffle
  deflate = NETCDF_deflate
  deflate_level =NETCDF_deflate_level

  if(alarmCheck .or. LIS_rc%endtime==1) then
     call LIS_create_output_directory('ROUTING')
     call LIS_create_restart_filename(n,filename,&
          'ROUTING','RAPID_router',wformat="netcdf")
     write(LIS_logunit,*) '[INFO] Writing RAPID routing restart ',trim(filename)

     if(LIS_masterproc) then
#if (defined USE_NETCDF4)
        status = nf90_create(path=filename, cmode=nf90_netcdf4, ncid = ftn)
        call LIS_verify(status,"Error in nf90_open in RAPID_routing_writerst")
#endif
#if (defined USE_NETCDF3)
        status = nf90_create(Path = filename, cmode = nf90_clobber, ncid = ftn)
        call LIS_verify(status, "Error in nf90_open in RAPID_routing_writerst")
#endif
        ! Define the dimensions.
        status = nf90_def_dim(ftn,'time',1,dimid_time)
        call LIS_verify(status, "Error in nf90_def_dim for time in RAPID_routing_writerst")
        status = nf90_def_dim(ftn,'rivid',RAPID_routing_struc(n)%n_riv_bas,dimid_riv_bas)
        call LIS_verify(status, "Error in nf90_def_dim for rivid in RAPID_routing_writerst")
        if(RAPID_routing_struc(n)%useens==2) then         ! ensemble mode
           status = nf90_def_dim(ftn,'ensemble',LIS_rc%nensem(n),dimid_ens)
           call LIS_verify(status, "Error in nf90_def_dim for ensemble in RAPID_routing_writerst")
           dimid_Qout_ens = (/dimid_riv_bas, dimid_ens, dimid_time/)
        else
           dim_Qout = (/dimid_riv_bas, dimid_time/)
        endif
        
        ! Define variables
        if(RAPID_routing_struc(n)%useens==2) then         ! ensemble mode
           status = nf90_def_var(ftn,"Qout",NF90_REAL,dimid_Qout_ens,varid_Qout)
           call LIS_verify(status, "Error in nf90_def_var for Qout in RAPID_ruting_writerst")
        else
           status = nf90_def_var(ftn,"Qout",NF90_REAL,dim_Qout,varid_Qout)
           call LIS_verify(status, "Error in nf90_def_var for Qout in RAPID_routing_writerst")
        endif

        !Define compression parameters
        call LIS_verify(nf90_def_var_deflate(ftn,varid_Qout,shuffle, deflate, deflate_level), &
                        'Error in nf90_def_var_deflate for Qout in RAPID_routing_writerst')

        ! Define variable attributes
        call LIS_verify(nf90_put_att(ftn,varid_Qout,'long_name','average river water discharge ' &
                        // 'downstream of each river reach'), 'nf90_put_att failed for long_name')
        call LIS_verify(nf90_put_att(ftn,varid_Qout,'unit','m3 s-1'), &
                        'nf90_put_att failed for unit')
        call LIS_verify(nf90_put_att(ftn,varid_Qout,'coordinates','lon lat'), &
                        'nf90_put_att failed for coordinates')
        call LIS_verify(nf90_put_att(ftn,varid_Qout,'grid_mapping','crs'), &
                        'nf90_put_att failed for grid_mapping')
        call LIS_verify(nf90_put_att(ftn,varid_Qout,'cell_methods','time: mean'), &
                        'nf90_put_att failed for cell_methods')

        ! Define global attributes
        call date_and_time(date,time,zone,values)
        call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value",LIS_rc%udef), &
             'nf90_put_att failed for missing_value')
        call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"title","LIS RAPID model restart"), &
             'nf90_put_att failed for title')
        call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution",trim(LIS_rc%institution)), &
             'nf90_put_att failed for institution')
        call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"source",'RAPID'), &
             'nf90_put_att failed for source')
        call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
             "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
             date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)), &
             'nf90_put_att failed for history')
        call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
             "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"), &
             'nf90_put_att failed for references')
        call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
             "CF-1.6"),'nf90_put_att failed for conventions') !CF version 1.6
        call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
             "website: http://lis.gsfc.nasa.gov/"),&
             'nf90_put_att failed for comment')

        ! End define mode. 
        call LIS_verify(nf90_enddef(ftn),'Error in ncf90_enddef in RAPID_routing_writerst')

        !Write data
        if(RAPID_routing_struc(n)%useens==2) then         ! ensemble mode
           status = nf90_put_var(ftn,varid_Qout,RAPID_routing_struc(n)%Qout_ens,&
                    (/1,1,1/), (/RAPID_routing_struc(n)%n_riv_bas,LIS_rc%nensem(n),1/))
           call LIS_verify(status,'Error in nf90_put_var in RAPID_routing_writerst')
        else
           status = nf90_put_var(ftn,varid_Qout,RAPID_routing_struc(n)%Qout,&
                    (/1,1/), (/RAPID_routing_struc(n)%n_riv_bas,1/))
           call LIS_verify(status,'Error in nf90_put_var in RAPID_routing_writerst')
        endif

        ! Close the file.
        call LIS_verify(nf90_close(ftn), "Error in nf90_close in RAPID_routing_writerst")
     endif !if(LIS_masterproc) then
  endif

end subroutine RAPID_routing_writerst

