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
! !ROUTINE: RAPID_routing_output
! \label{RAPID_routing_output}
!
! !REVISION HISTORY:
! 31 Jan 2022: Yeosang Yoon;  Initial implementation
! 08 Aug 2022: Yeosang Yoon;  Add code to check 'Output methodology' option
! 31 Aug 2022: Yeosang Yoon;  Fix to broadcast a message (output file)
! 27 Apr 2023: Eric Kemp; Updated length of output file.
! 20 Mar 2024: Yeosang Yoon; Support to run with ensemble mode

subroutine RAPID_routing_output(n)
  
  use ESMF
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_histDataMod
  use LIS_historyMod
  use LIS_fileIOMod
  use RAPID_routingMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
  
  integer, intent(in)   :: n 
  
  character(len=LIS_CONST_PATH_LEN) :: filename ! EMK
  integer               :: ftn
  integer               :: status
  integer               :: mo, da
  logical               :: alarmCheck
  integer               :: dimid_time, dimid_riv_bas
  integer               :: dimid_rivid, dimid_nv, dimid_nerr, dimid_Qout(2)
  integer               :: dimid_ens, dimid_Qout_ens(3), varid_ens, i          ! ensemble mode
  integer, allocatable  :: ensval(:)
  integer               :: varid_Qout, varid_Qout_err,varid_rivid,varid_time
  integer               :: varid_time_bnds, varid_lon, varid_lat,varid_crs
  integer               :: shuffle, deflate, deflate_level
 
  integer, dimension(8) :: values
  character(len=8)      :: xtime_begin_date
  character(len=6)      :: xtime_begin_time
  character(len=50)     :: xtime_units, xtime_timeInc
  character(len=8)      :: date
  character(len=10)     :: time
  character(len=5)      :: zone

  shuffle = NETCDF_shuffle
  deflate = NETCDF_deflate
  deflate_level =NETCDF_deflate_level

  alarmCheck = .false.
  if(trim(LIS_rc%wopt).ne."none") then
     if(LIS_rc%time >= LIS_histData(n)%time) then
        if(LIS_rc%output_at_specifictime.eq.1) then
           if(LIS_histData(n)%month.eq.-1) then
              mo = LIS_rc%mo
           else
              mo = LIS_histData(n)%month
           endif
           if(LIS_histData(n)%day.eq.-1) then
              da = LIS_rc%da
           else
              da = LIS_histData(n)%day
           endif

           if(LIS_rc%mo.eq.mo.and.&
              LIS_rc%da.eq.da.and.&
              LIS_rc%hr.eq.LIS_histData(n)%hour.and.&
              LIS_rc%mn.eq.LIS_histData(n)%min.and.&
              LIS_rc%ss.eq.LIS_histData(n)%sec) then
              alarmCheck = .true.
           endif
        else
           alarmCheck = LIS_isAlarmRinging(LIS_rc,&
                "RAPID router output alarm")
        endif
     
        if(alarmCheck) then
           ! Create output file name
           call LIS_create_output_directory('ROUTING')
           call LIS_create_output_filename(n,filename,&
                   model_name='ROUTING',&
                   writeint=RAPID_routing_struc(n)%outInterval)
           write(LIS_logunit,*) '[INFO] Writing routing model output to:  ',trim(filename)
 
           if(LIS_masterproc) then
              ! Create file
#if (defined USE_NETCDF4)
              status = nf90_create(path=filename, cmode=nf90_netcdf4, ncid = ftn)
              call LIS_verify(status,"Error in nf90_open in RAPID_routing_output")
#endif
#if (defined USE_NETCDF3)
              status = nf90_create(Path = filename, cmode = nf90_clobber, ncid = ftn)
              call LIS_verify(status, "Error in nf90_open in RAPID_routing_output")
#endif
              status = nf90_def_dim(ftn,'rivid',RAPID_routing_struc(n)%n_riv_bas,dimid_rivid)
              call LIS_verify(status, "nf90_put_att failed for rivid in RAPID_routing_output")
              status = nf90_def_dim(ftn,'time',1,dimid_time)
              call LIS_verify(status, "nf90_put_att failed for time in RAPID_routing_output")
              !status = nf90_def_dim(ftn,'nv',2,dimid_nv)
              !call LIS_verify(status, "nf90_put_att failed for nv in RAPID_routing_output")
              !status = nf90_def_dim(ftn,'nerr',3,dimid_nerr)
              !call LIS_verify(status, "nf90_put_att failed for nerr in RAPID_routing_output")

              if(RAPID_routing_struc(n)%useens==2) then         ! ensemble mode
                 status = nf90_def_dim(ftn,'ensemble',LIS_rc%nensem(n),dimid_ens)
                 call LIS_verify(status, "Error in nf90_def_dim for ensemble in RAPID_routing_output")
                 dimid_Qout_ens = (/dimid_rivid, dimid_ens, dimid_time/)
              else
                 !dimid_Qout = (/dimid_rivid, dimid_time/)
                 dimid_Qout = (/dimid_time, dimid_rivid/)
              endif

              !Define variables
              if(RAPID_routing_struc(n)%useens==2) then         ! ensemble mode
                 status = nf90_def_var(ftn,"ensemble",NF90_INT,dimid_ens,varid_ens)
                 call LIS_verify(status, "nf90_def_var failed for crs in RAPID_routing_output")
                 status = nf90_def_var(ftn,"Qout",NF90_REAL,dimid_Qout_ens,varid_Qout)
                 call LIS_verify(status, "Error in nf90_def_var for Qout in RAPID_ruting_output")
              else
                 status = nf90_def_var(ftn,"Qout",NF90_REAL,dimid_Qout,varid_Qout)
                 call LIS_verify(status, "nf90_def_var failed for Qout in RAPID_routing_output")
              endif
              !status = nf90_def_var(ftn,"Qout_err",NF90_REAL,&
              !         (/dimid_rivid,dimid_nerr/),varid_Qout_err)
              !call LIS_verify(status, "nf90_def_var failed for Qout_err in RAPID_routing_output")
              status = nf90_def_var(ftn,"rivid",NF90_INT,dimid_rivid,varid_rivid)
              call LIS_verify(status, "nf90_def_var failed for rivid in RAPID_routing_output")
              status = nf90_def_var(ftn,"time",NF90_FLOAT,dimid_time,varid_time)
              call LIS_verify(status, "nf90_def_var failed for time in RAPID_routing_output")
              !status = nf90_def_var(ftn,"time_bnds",NF90_INT,&
              !         (/dimid_nv,dimid_time/),varid_time_bnds)
              !call LIS_verify(status, "nf90_def_var failed for time_bnds in RAPID_routing_output")
              status = nf90_def_var(ftn,"lon",NF90_FLOAT,dimid_rivid,varid_lon)
              call LIS_verify(status, "nf90_def_var failed for lon in RAPID_routing_output")
              status = nf90_def_var(ftn,"lat",NF90_FLOAT,dimid_rivid,varid_lat)
              call LIS_verify(status, "nf90_def_var failed for lat in RAPID_routing_output")
              status = nf90_def_var(ftn,"crs",NF90_INT,varid_crs)
              call LIS_verify(status, "nf90_def_var failed for crs in RAPID_routing_output")
          
              !Define compression parameters
              call LIS_verify(nf90_def_var_deflate(ftn,varid_Qout,shuffle,deflate,deflate_level),&
                           'nf90_def_var_deflate failed for Qout in RAPID_routing_output')  
              !call LIS_verify(nf90_def_var_deflate(ftn,varid_Qout_err,shuffle,deflate,deflate_level),&
              !                'nf90_def_var_deflate failed for Qout_err in RAPID_routing_output')
              call LIS_verify(nf90_def_var_deflate(ftn,varid_rivid,shuffle,deflate,deflate_level),&
                           'nf90_def_var_deflate failed for rivid in RAPID_routing_output')
              call LIS_verify(nf90_def_var_deflate(ftn,varid_time,shuffle,deflate,deflate_level),&
                           'nf90_def_var_deflate failed for time in RAPID_routing_output')
              !call LIS_verify(nf90_def_var_deflate(ftn,varid_time_bnds,shuffle,deflate,deflate_level),&
              !                'nf90_def_var_deflate failed for time_bnds in RAPID_routing_output')
              call LIS_verify(nf90_def_var_deflate(ftn,varid_lon,shuffle,deflate,deflate_level),&
                           'nf90_def_var_deflate failed for lon in RAPID_routing_output')
              call LIS_verify(nf90_def_var_deflate(ftn,varid_lat,shuffle,deflate,deflate_level),&
                           'nf90_def_var_deflate failed for lat in RAPID_routing_output')
              !call LIS_verify(nf90_def_var_deflate(ftn,varid_crs,shuffle,deflate,deflate_level),&
              !             'nf90_def_var_deflate failed for crs in RAPID_routing_output')
              if(RAPID_routing_struc(n)%useens==2) then         ! ensemble mode
                 call LIS_verify(nf90_def_var_deflate(ftn,varid_ens,shuffle,deflate,deflate_level),&
                              'nf90_def_var_deflate failed for ens in RAPID_routing_output')
              endif

              ! Define variable attributes
              call LIS_verify(nf90_put_att(ftn,varid_Qout,'long_name','average river water ' & 
                                        // 'discharge downstream of each river reach'),   &
                                           'nf90_put_att failed for Qout')
              call LIS_verify(nf90_put_att(ftn,varid_Qout,'units','m3 s-1'), &
                           'nf90_put_att failed for Qout')
              call LIS_verify(nf90_put_att(ftn,varid_Qout,'coordinates','lon lat'), &
                           'nf90_put_att failed for Qout')
              call LIS_verify(nf90_put_att(ftn,varid_Qout,'grid_mapping','crs'), &
                           'nf90_put_att failed for Qout')
              call LIS_verify(nf90_put_att(ftn,varid_Qout,'cell_methods','time: mean'), &
                           'nf90_put_att failed for Qout')

              !call LIS_verify(nf90_put_att(ftn,varid_Qout_err,'long_name','average river water '     &
              !                             //'discharge uncertainty downstream of each river reach'),&
              !                'nf90_put_att failed for Qout_err')
              !call LIS_verify(nf90_put_att(ftn,varid_Qout_err,'units','m3 s-1'), &
              !                'nf90_put_att failed for Qout_err')
              !call LIS_verify(nf90_put_att(ftn,varid_Qout_err,'coordinates','lon lat'), &
              !                'nf90_put_att failed for Qout_err')
              !call LIS_verify(nf90_put_att(ftn,varid_Qout_err,'grid_mapping','crs'), &
              !                'nf90_put_att failed for Qout_err')
              !call LIS_verify(nf90_put_att(ftn,varid_Qout_err,'cell_methods','time: mean'), &
              !                'nf90_put_att failed for Qout_err')

              call LIS_verify(nf90_put_att(ftn,varid_rivid,'long_name','unique identifier for each '& 
                                        //'river reach'),                                        &
                           'nf90_put_att failed for rivid')
              call LIS_verify(nf90_put_att(ftn,varid_rivid,'units','1'), &
                           'nf90_put_att failed for rivid')
              call LIS_verify(nf90_put_att(ftn,varid_rivid,'cf_role','timeseries_id'), &
                           'nf90_put_att failed for rivid')

              ! LIS output is always writing output for a single time
              write(xtime_units,200) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, &
                   LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
200           format('minutes since ',I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':', &
                   I2.2,':',I2.2)
              write(xtime_begin_date, fmt='(I4.4,I2.2,I2.2)') &
                   LIS_rc%yr, LIS_rc%mo, LIS_rc%da
              write(xtime_begin_time, fmt='(I2.2,I2.2,I2.2)') &
                   LIS_rc%hr, LIS_rc%mn, LIS_rc%ss
              write(xtime_timeInc, fmt='(I20)') nint(RAPID_routing_struc(n)%outInterval)

              call LIS_verify(nf90_put_att(ftn,varid_time,'unit',trim(xtime_units)), &
                           'nf90_put_att failed for time')
              call LIS_verify(nf90_put_att(ftn,varid_time,'long_name','time'), &
                           'nf90_put_att failed for time')
              call LIS_verify(nf90_put_att(ftn,varid_time,'time_increment',trim(adjustl(xtime_timeInc))), &
                           'nf90_put_att failed for time')
              call LIS_verify(nf90_put_att(ftn,varid_time,'begin_date',xtime_begin_date), &
                           'nf90_put_att failed for time')
              call LIS_verify(nf90_put_att(ftn,varid_time,'begin_time',xtime_begin_time), &
                           'nf90_put_att failed for time')
              !call LIS_verify(nf90_put_att(ftn,varid_time,'bounds','time_bnds'), &
              !                'nf90_put_att failed for time')

              call LIS_verify(nf90_put_att(ftn,varid_lon,'standard_name','longitude'), &
                           'nf90_put_att failed for lon')
              call LIS_verify(nf90_put_att(ftn,varid_lon,'long_name','longitude of a point '&
                                        // 'related to each river reach'),                &
                           'nf90_put_att failed for lon')
              call LIS_verify(nf90_put_att(ftn,varid_lon,'units','degrees_east'), &
                           'nf90_put_att failed for lon')
              call LIS_verify(nf90_put_att(ftn,varid_lon,'axis','X'), &
                           'nf90_put_att failed for lon')

              call LIS_verify(nf90_put_att(ftn,varid_lat,'standard_name','latitude'), &
                           'nf90_put_att failed for lat')
              call LIS_verify(nf90_put_att(ftn,varid_lat,'long_name','latitude of a point '&
                                        // 'related to each river reach'),                &
                           'nf90_put_att failed for lat')
              call LIS_verify(nf90_put_att(ftn,varid_lat,'units','degrees_north'), &
                           'nf90_put_att failed for lat')
              call LIS_verify(nf90_put_att(ftn,varid_lat,'axis','Y'), &
                           'nf90_put_att failed for lat')

              call LIS_verify(nf90_put_att(ftn,varid_crs,'grid_mapping_name','latitude_longitude'), &
                           'nf90_put_att failed for crs')
              call LIS_verify(nf90_put_att(ftn,varid_crs,'semi_major_axis','6378137.'), &
                           'nf90_put_att failed for crs')
              call LIS_verify(nf90_put_att(ftn,varid_crs,'inverse_flattening','298.257223563'), &
                           'nf90_put_att failed for crs')

              if(RAPID_routing_struc(n)%useens==2) then         ! ensemble mode
                 call LIS_verify(nf90_put_att(ftn,varid_ens,"units","ensemble number"), &
                              'nf90_put_att failed for ensemble units')
                 call LIS_verify(nf90_put_att(ftn,varid_ens, "long_name","Ensemble numbers"),&
                              'nf90_put_att failed for ensemble long_name failed')
                 allocate(ensval(LIS_rc%nensem(n)))
                 do i = 1, LIS_rc%nensem(n)
                    ensval(i) = i
                 enddo
              endif

              ! Define global attributes
              call date_and_time(date,time,zone,values)
              call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value",LIS_rc%udef), &
                   'nf90_put_att failed for missing_value')
              call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"title","LIS RAPID model output"), &
                   'nf90_put_att failed for title')
              call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution",trim(LIS_rc%institution)), &
                   'nf90_put_att failed for institution')
              call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"source",'RAPID v1.8.0'), &
                   'nf90_put_att failed for source')
              call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
                   "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
                   date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)), &
                   'nf90_put_att failed for history')
              call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"references",       &
                           'https://github.com/c-h-david/ra'//              &
                           'pid/, http://dx.doi.org/10.1175/2011JHM1345.1'),&
                           'nf90_put_att failed for references')
              call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
                   "CF-1.6"),'nf90_put_att failed for conventions') !CF version 1.6
              call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
                   'website: http://lis.gsfc.nasa.gov/'),&
                   'nf90_put_att failed for comment')

              ! End define mode.
              call LIS_verify(nf90_enddef(ftn),'Error in ncf90_enddef in RAPID_routing_output')

              ! Write data
              status = nf90_put_var(ftn,varid_rivid,RAPID_routing_struc(n)%riv_bas_id)
              call LIS_verify(status,'Error in nf90_put_var for rivid in RAPID_routing_output')
              status = nf90_put_var(ftn,varid_lon,RAPID_routing_struc(n)%riv_tot_lon)
              call LIS_verify(status,'Error in nf90_put_var for lon in RAPID_routing_output')
              status = nf90_put_var(ftn,varid_lat,RAPID_routing_struc(n)%riv_tot_lat)
              call LIS_verify(status,'Error in nf90_put_var for lat in RAPID_routing_output')
              status = nf90_put_var(ftn,varid_time,0.0)
              call LIS_verify(status,'Error in nf90_put_var for time in RAPID_routing_output')
              status = nf90_put_var(ftn,varid_crs,0)
              call LIS_verify(status,'Error in nf90_put_var for crs in RAPID_routing_output')

              if(RAPID_routing_struc(n)%useens==2) then         ! ensemble mode
                 status = nf90_put_var(ftn,varid_ens,ensval,(/1/),(/LIS_rc%nensem(n)/))
                 call LIS_verify(status,'Error in nf90_put_var for ens in RAPID_routing_output')
                 status = nf90_put_var(ftn,varid_Qout,RAPID_routing_struc(n)%Qout_ens,&
                          (/1,1,1/), (/RAPID_routing_struc(n)%n_riv_bas,LIS_rc%nensem(n),1/))
                 call LIS_verify(status,'Error in nf90_put_var in RAPID_routing_output')
              else
                 status = nf90_put_var(ftn,varid_Qout,&
                          reshape(RAPID_routing_struc(n)%Qout,(/1,RAPID_routing_struc(n)%n_riv_bas/)))
                 call LIS_verify(status,'Error in nf90_put_var for Qout in RAPID_routing_output')
              endif

              ! Close the file.
              call LIS_verify(nf90_close(ftn), "Error in nf90_close in RAPID_routing_output")
           endif
        endif
     endif
  endif
end subroutine RAPID_routing_output
