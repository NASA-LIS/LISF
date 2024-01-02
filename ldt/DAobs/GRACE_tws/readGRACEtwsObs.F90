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
!BOP
! 
! !ROUTINE: readGRACEtwsObs
! \label{readGRACEtwsObs}
! 
! !REVISION HISTORY: 
!  14 Mar 2013: Sujay Kumar, Initial Specification
!
!  21 DEc 2017: Augusto Getirana, Account for surface water storage (SWS) by subtracting the SWS signal from GRACE TWS. 
!  This process results in the land water storage (LWS), which is consistent with the LSM water storage. 
!  Getirana, A., et al., 2017b. Rivers and floodplains as key components of global terrestrial water storage variability. 
!  Geophysical Research Letters, 44. DOI: 10.1002/2017GL074684.
! 
! !INTERFACE: 
subroutine readGRACEtwsObs(n)
! !USES:   
  use LDT_coreMod
  use LDT_logMod
  use LDT_historyMod
  use LDT_DAobsDataMod
  use LDT_timeMgrMod
  use LDT_constantsMod
  use GRACEtws_obsMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
! 
! !DESCRIPTION: 
! 
! This subroutine reads the raw GRACE data, computes the anomalies and 
! generates a new set of GRACE observations by incorporating the anomalies
! into the model generated terrestrial water storage observations. 
!
! TWS outputs from LIS is expected to be in units of mm.
!
!EOP
  character(len=LDT_CONST_PATH_LEN) :: fname,filename
  integer               :: c,r,c1,r1,k,t,iret
  integer               :: ftn
  integer               :: yr,mo,da,hr
  integer               :: timeId, tId, tbId, lweId, scaleId,errId,err1Id
  real                  :: tws_data_ip(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                  :: tws_data(GRACEtwsobs%nc,GRACEtwsobs%nr)
  real                  :: output_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real,    allocatable  :: output_basin_data(:)
  integer, allocatable  :: output_basin_data_count(:)
  real,    allocatable  :: var_merr(:)
  real,    allocatable  :: var_lerr(:)
  real,    allocatable  :: basin_merr(:)
  real,    allocatable  :: basin_lerr(:)
  integer, allocatable  :: nerr(:)
  real                  :: expdbm,expdbl
  real                  :: lat1,lon1,lat2,lon2,dist,betam, betal
  real                  :: output_err_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                  :: output_merr_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                  :: output_lerr_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1             :: output_bitmap(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1             :: output_err_bitmap(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                  :: input_data(GRACEtwsobs%gracenc*GRACEtwsobs%gracenr)
  logical*1             :: input_bitmap(GRACEtwsobs%gracenc*GRACEtwsobs%gracenr)
  integer               :: npts(GRACEtwsobs%gracenc,GRACEtwsobs%gracenr) 
  logical               :: file_exists
  real                  :: dt
  integer               :: cat_val, twsid
  integer               :: currTime,yyyy,mm,dd,hh
  character*10          :: ftime
  logical               :: valid_data
  integer               :: md_nc

  !ag (21 Dec 2017)
  real                  :: sws_data_ip(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                  :: sws_data(GRACEtwsobs%nc,GRACEtwsobs%nr)
  real                  :: sws_time(GRACEtwsobs%nc,GRACEtwsobs%nr)
  real                  :: avg
  integer               :: swsid,date,iloc(1)

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
!
! At initial time, read the data into memory. At all other times, 
! simply index into the correct temporal location of the data. 
!
  if( GRACEtwsobs%startMode ) then 
     GRACEtwsobs%startMode = .false. 

     ! Read 'raw' GRACE obs file:
     inquire(file=trim(GRACEtwsobs%gracefile),exist=file_exists) 
     if(file_exists) then 
        write(LDT_logunit,*) '[INFO] Reading GRACE data '//&
             trim(GRACEtwsobs%gracefile)
        call LDT_verify(nf90_open(path=trim(GRACEtwsobs%gracefile),&
             mode=nf90_nowrite, ncid = ftn), &
             'nf90_open failed in readGRACEtwsObs ')
        call LDT_verify(nf90_inq_dimid(ftn,'time',timeId),&
             'nf90_inq_dimid failed in readGRACEtwsObs')
        call LDT_verify(nf90_inquire_dimension(ftn,timeId,&
             len=GRACEtwsobs%tdims),&
             'nf90_inq_dimension failed in readGRACEtwsObs')

        allocate(GRACEtwsobs%tvals(GRACEtwsobs%tdims))
        allocate(GRACEtwsobs%time_bounds(2,GRACEtwsobs%tdims))
        allocate(GRACEtwsobs%lwe_thickness(GRACEtwsobs%gracenc, &
             GRACEtwsobs%gracenr,GRACEtwsobs%tdims))
        allocate(GRACEtwsobs%twsavg(GRACEtwsobs%gracenc, &
             GRACEtwsobs%gracenr))

        call LDT_verify(nf90_inq_varid(ftn,'lwe_thickness',lweId),&
             'nf90_inq_varid failed for lwe_thickness in readGRACEtwsObs')
        call LDT_verify(nf90_get_var(ftn,lweId,GRACEtwsobs%lwe_thickness),&
             'nf90_get_var failed for lwe_thickness in readGRACEtwsObs')
        
        call LDT_verify(nf90_inq_varid(ftn,'time_bounds',tbId),&
             'nf90_inq_varid failed for time_bounds in readGRACEtwsObs')
        call LDT_verify(nf90_get_var(ftn,tbId,GRACEtwsobs%time_bounds),&
             'nf90_get_var failed for time_bounds in readGRACEtwsObs')
        
        call LDT_verify(nf90_inq_varid(ftn,'time',tId),&
             'nf90_inq_varid failed for time in readGRACEtwsObs')
        call LDT_verify(nf90_get_var(ftn,tId,GRACEtwsobs%tvals),&
             'nf90_get_var failed for time in readGRACEtwsObs')
        
        call LDT_verify(nf90_close(ftn))
        
        ! Optional GRACE scale factor file to scale data:
        if( trim(GRACEtwsobs%gracescalefile) .eq. "none" ) then
           write(LDT_logunit,*) '[INFO] Using unscaled GRACE data ...'

        else  ! Read in file
           write(LDT_logunit,*) '[INFO] Reading GRACE scalefactor file: ' &
                //trim(GRACEtwsobs%gracescalefile)
           call LDT_verify(nf90_open(path=trim(GRACEtwsobs%gracescalefile),&
                mode=nf90_nowrite, ncid = ftn), &
                'nf90_open failed in readGRACEtwsObs scale factor')
           allocate(GRACEtwsobs%scalefactor(GRACEtwsobs%gracenc, &
                GRACEtwsobs%gracenr))
           
           call LDT_verify(nf90_inq_varid(ftn,'SCALE_FACTOR',scaleId), &
                'nf90_inq_varid failed for scale_factor in readGRACEtwsObs')
           call LDT_verify(nf90_get_var(ftn,scaleId,GRACEtwsobs%scalefactor),&
                'nf90_get_var failed for scale_factor in readGRACEtwsObs') 
           call LDT_verify(nf90_close(ftn))
        end if
        
        ! Read GRACE obs measurement and 'leakage' error file:
        if( GRACEtwsobs%graceerrfile .ne. "none" ) then
           write(LDT_logunit,*) '[INFO] Reading GRACE error file: ' &
                //trim(GRACEtwsobs%graceerrfile)
           call LDT_verify(nf90_open(path=trim(GRACEtwsobs%graceerrfile),&
                mode=nf90_nowrite, ncid = ftn), &
                'nf90_open failed in readGRACEtwsObs error file')

           allocate(GRACEtwsobs%merror(GRACEtwsobs%gracenc, &
                GRACEtwsobs%gracenr))  
           allocate(GRACEtwsobs%lerror(GRACEtwsobs%gracenc, &
                GRACEtwsobs%gracenr))  
           call LDT_verify(nf90_inq_varid(ftn,'MEASUREMENT_ERROR',errId),&
                'nf90_inq_varid failed for MEASUREMENT_ERROR in readGRACEtwsObs')
           call LDT_verify(nf90_get_var(ftn,errId,GRACEtwsobs%merror),&
                'nf90_get_var failed for MEASUREMENT_ERROR in readGRACEtwsObs') 

           call LDT_verify(nf90_inq_varid(ftn,'LEAKAGE_ERROR',err1Id),&
                'nf90_inq_varid failed for LEAKAGE_ERROR in readGRACEtwsObs')
           call LDT_verify(nf90_get_var(ftn,errId,GRACEtwsobs%lerror),&
                'nf90_get_var failed for LEAKAGE_ERROR in readGRACEtwsObs')  
           call LDT_verify(nf90_close(ftn))
        end if
        
        npts = 0 
        GRACEtwsobs%twsavg = 0 
        
        ! Loop over GRACE time dimensions:
        do k=1,GRACEtwsobs%tdims
           currtime = GRACEtwsobs%tvals(k)*24.0 + GRACEtwsobs%refTime
           call LDT_julhr_date(currtime, yr, mo, da, hr)

           ! Within baseline GRACE data period for calculating TWS climatology:
           if(yr.ge.GRACEtwsobs%b_syr.and.&
                yr.le.GRACEtwsobs%b_eyr) then

              ! If no scaling applied:
              if(GRACEtwsobs%gracescalefile == "none") then  ! BZ
                 do r=1,GRACEtwsobs%gracenr
                    do c=1,GRACEtwsobs%gracenc
                       if(GRACEtwsobs%lwe_thickness(c,r,k).ne.32767.0) then 
                          GRACEtwsobs%twsavg(c,r) = GRACEtwsobs%twsavg(c,r) + & 
                               GRACEtwsobs%lwe_thickness(c,r,k)
                          npts(c,r) = npts(c,r) + 1
                       endif
                    enddo
                 enddo
                 
              ! If scaling applied:
              else 
                 do r=1,GRACEtwsobs%gracenr
                    do c=1,GRACEtwsobs%gracenc
                       if(GRACEtwsobs%lwe_thickness(c,r,k).ne.32767.0) then 
                          GRACEtwsobs%twsavg(c,r) = GRACEtwsobs%twsavg(c,r) + & 
                               GRACEtwsobs%lwe_thickness(c,r,k) * &
                               GRACEtwsobs%scalefactor(c,r)
                          npts(c,r) = npts(c,r) + 1
                       endif
                    enddo
                 enddo
              endif
              
           endif   ! Conditional for baseline averaging period
        enddo

        ! Calculate average and convert to mm.
        do r=1,GRACEtwsobs%gracenr
           do c=1,GRACEtwsobs%gracenc
              if( npts(c,r).gt.0 ) then   ! Convert to mm. 
                 GRACEtwsobs%twsavg(c,r) = GRACEtwsobs%twsavg(c,r)&
                                         /npts(c,r)
              else
                 GRACEtwsobs%twsavg(c,r) = LDT_rc%udef
              endif
           enddo
        enddo
        write(LDT_logunit,*) '[INFO] Finished reading GRACE data '//&
             trim(GRACEtwsobs%gracefile)

        !ag (21 Dec 2017)
        allocate(GRACEtwsobs%swsavg(GRACEtwsobs%gracenc, &
             GRACEtwsobs%gracenr,GRACEtwsobs%tdims))
        allocate(GRACEtwsobs%nswsavg(GRACEtwsobs%gracenc, &
             GRACEtwsobs%gracenr,GRACEtwsobs%tdims))
        allocate(GRACEtwsobs%swsdate(GRACEtwsobs%tdims))
        GRACEtwsobs%swsavg = 0 
        GRACEtwsobs%nswsavg = 0 
        
        !get GRACE dates to match with SWS simulation
        do k=1,GRACEtwsobs%tdims
          currTime = GRACEtwsobs%tvals(k)*24+GRACEtwsobs%reftime
          call LDT_julhr_date(currTime, yyyy,mm,dd,hh)
          if(k>1)then
            !if 2 grace observations fall in the same month (this can be improved)
            if(GRACEtwsobs%swsdate(k-1)==yyyy*100+mm)then
              if(mm==12)then
                mm=1
              else
                mm=mm+1
              endif
            endif
          endif
          GRACEtwsobs%swsdate(k)=yyyy*100+mm
        enddo


     else
        write(LDT_logunit,*) '[ERR] GRACE raw obs file '&
             //trim(GRACEtwsobs%gracefile)//& 
             'does not exist ...'
        call LDT_endrun()
     endif
  endif

  ! Account for surface water storage (SWS) by subtracting the SWS signal from GRACE TWS. This process results
  ! in the land water storage (LWS), which is consistent with the LSM water storage. 
  ! Getirana, A., et al., 2017b. Rivers and floodplains as key components of global terrestrial water storage variability. GRL, 44. DOI: 10.1002/2017GL074684.
  if(GRACEtwsobs%swsflag.eq.1)then
  
    !Read during the first pass for averaging:
    if(LDT_rc%pass_id.eq.1) then 
      date=LDT_rc%yr*100+LDT_rc%mo
      
      !Read data for dates with GRACE obs
      if(minval(abs(date-GRACEtwsobs%swsdate))==0)then
              
        call create_lsm_twsoutput_filename(GRACEtwsobs%nest, &
             GRACEtwsobs%format,&
             fname,GRACEtwsobs%odir, GRACEtwsobs%wstyle, &
             GRACEtwsobs%wopt,'ROUTING')

        inquire(file=trim(fname),exist=file_exists)
        if(file_exists) then 
           write(LDT_logunit,*) '[INFO] reading routing scheme output ',trim(fname)
        
           if(GRACEtwsobs%format.eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
              
              iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
              call LDT_verify(iret, 'Error opening file '//trim(fname))
              
              iret = nf90_inq_varid(ftn, "SWS_tavg", swsid)
              call LDT_verify(iret, 'Error in nf90_inq_varid: SWS_tavg')
              
              iret = nf90_get_var(ftn, swsid, sws_data,start = (/1,1/),&
                   count = (/GRACEtwsobs%nc,GRACEtwsobs%nr/))
              call LDT_verify(iret,'Error in nf90_get_var: SWS_tavg')

              iret = nf90_close(ftn)
              call LDT_verify(iret,'Error in nf90_close')
#endif
              call transformLISoutToGRACEgrid(n,sws_data,sws_data_ip)

              iloc(:)=minloc(abs(date-GRACEtwsobs%swsdate))
              k=iloc(1)
              do r=1,LDT_rc%lnr(n)
                 do c=1,LDT_rc%lnc(n)
                    if(sws_data_ip(c,r).ne.LDT_rc%udef) then 
                       GRACEtwsobs%swsavg(c,r,k) = &
                            GRACEtwsobs%swsavg(c,r,k) + sws_data_ip(c,r)
                       GRACEtwsobs%nswsavg(c,r,k) = &
                            GRACEtwsobs%nswsavg(c,r,k) + 1
                    endif
                 enddo
              enddo

           endif
        else
           write(LDT_logunit,*) '[ERR] LIS file '//trim(fname)
           write(LDT_logunit,*) '[ERR] not found. Program stopping ...'
           call LDT_endrun()
        endif
      endif

      ! At the end of the first cycle. 
      if(LDT_rc%endtime.eq.1) then 
        ! Loop over GRACE time dimensions:
        do r=1,LDT_rc%lnr(n)
          do c=1,LDT_rc%lnc(n)
            do k=1,GRACEtwsobs%tdims
              if(GRACEtwsobs%nswsavg(c,r,k)>0.and.GRACEtwsobs%lwe_thickness(c,r,k)/=32767) then 
                 GRACEtwsobs%swsavg(c,r,k) = &
                      GRACEtwsobs%swsavg(c,r,k)/&
                      GRACEtwsobs%nswsavg(c,r,k) 
              else
                 GRACEtwsobs%swsavg(c,r,k) = 0
              endif
            enddo
            if(count(GRACEtwsobs%nswsavg(c,r,:)>0)>0)then
              avg=sum(GRACEtwsobs%swsavg(c,r,:),GRACEtwsobs%nswsavg(c,r,:)>0)/count(GRACEtwsobs%nswsavg(c,r,:)>0)
              where(GRACEtwsobs%nswsavg(c,r,:)>0)GRACEtwsobs%swsavg(c,r,:)=GRACEtwsobs%swsavg(c,r,:)-avg
            endif
          enddo
        enddo
      endif
      
    endif
  endif

  
  ! Read during the first pass for averaging:
  if(LDT_rc%pass_id.eq.1) then 
     
     call create_lsm_twsoutput_filename(GRACEtwsobs%nest, &
          GRACEtwsobs%format,&
          fname,GRACEtwsobs%odir, GRACEtwsobs%wstyle, &
          GRACEtwsobs%wopt,'SURFACEMODEL')
     
     ! Average only between baseline years:
     if(LDT_rc%yr.ge.GRACEtwsobs%b_syr.and.&
          LDT_rc%yr.le.GRACEtwsobs%b_eyr) then

        inquire(file=trim(fname),exist=file_exists)
        if(file_exists) then 
           write(LDT_logunit,*) '[INFO] reading LSM output ',trim(fname)
        
           if(GRACEtwsobs%format.eq."netcdf") then 
#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
              
              iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
              call LDT_verify(iret, 'Error opening file '//trim(fname))
              
              iret = nf90_inq_varid(ftn, "TWS_tavg", twsid)
              call LDT_verify(iret, 'Error in nf90_inq_varid: TWS_tavg')
              
              iret = nf90_get_var(ftn, twsid, tws_data,start = (/1,1/),&
                   count = (/GRACEtwsobs%nc,GRACEtwsobs%nr/))
              call LDT_verify(iret,'Error in nf90_get_var: TWS_tavg')

              iret = nf90_close(ftn)
              call LDT_verify(iret,'Error in nf90_close')
#endif
              call transformLISoutToGRACEgrid(n,tws_data,tws_data_ip)

              do r=1,LDT_rc%lnr(n)
                 do c=1,LDT_rc%lnc(n)
                    if(tws_data_ip(c,r).ne.LDT_rc%udef) then 
                       GRACEtwsobs%lisavg(c,r) = &
                            GRACEtwsobs%lisavg(c,r) + tws_data_ip(c,r)
                       GRACEtwsobs%nlisavg(c,r) = &
                            GRACEtwsobs%nlisavg(c,r) + 1
                    endif
                 enddo
              enddo

           endif
        else
           write(LDT_logunit,*) '[ERR] LIS file '//trim(fname)
           write(LDT_logunit,*) '[ERR] not found. Program stopping ...'
           call LDT_endrun()
        endif
     endif
     ! At the end of the first cycle. 

     if(LDT_rc%pass_id.eq.1.and.LDT_rc%endtime.eq.1) then 
        do r=1,LDT_rc%lnr(n)
           do c=1,LDT_rc%lnc(n)
              if(GRACEtwsobs%nlisavg(c,r).gt.0) then 
                 GRACEtwsobs%lisavg(c,r) = &
                      GRACEtwsobs%lisavg(c,r)/&
                      GRACEtwsobs%nlisavg(c,r) 
              else
                 GRACEtwsobs%lisavg(c,r) = LDT_rc%udef
              endif
           enddo
        enddo
     endif

  ! Second pass through files:
  elseif(LDT_rc%pass_id.eq.2) then 
     
     call LDT_get_julhr(LDT_rc%yr,LDT_rc%mo,LDT_rc%da,&
          LDT_rc%hr,LDT_rc%mn,LDT_rc%ss,currTime)
     
     dt = float((currTime-GRACEtwsobs%refTime))/24.0

    if(GRACEtwsobs%datasource.eq."GRACE TWS Mascon 0.25 deg") then 
       md_nc=720
    elseif(GRACEtwsobs%datasource.eq."GRACE TWS Mascon 0.5 deg") then 
       md_nc=360
    elseif(GRACEtwsobs%datasource.eq."GRACE TWS Original 1 deg") then 
       md_nc=180
    endif

     do k=1,GRACEtwsobs%tdims
        if(floor(GRACEtwsobs%tvals(k)).eq.dt) then 
           !  Interpolate the anomaly map to the LIS output grid
           ! 
           input_data = LDT_rc%udef
           
           ! NO scaling performed:
           if(GRACEtwsobs%gracescalefile .eq. "none") then       
            do r=1,GRACEtwsobs%gracenr
               do c=1,GRACEtwsobs%gracenc               
                  if(GRACEtwsobs%twsavg(c,r).ne.LDT_rc%udef.or.&
                       GRACEtwsobs%lwe_thickness(c,r,k).ne.32767) then 
                     ! Convert GRACE 0 to 360dg domain to -180 to 180E:
!                     if(c.le.180) then 
                     if(c.le.md_nc) then 
!                        input_data(c+(r-1)*GRACEtwsobs%gracenc+180) = & 
                        input_data(c+(r-1)*GRACEtwsobs%gracenc+md_nc) = & 
                             (GRACEtwsobs%lwe_thickness(c,r,k)-&
                             GRACEtwsobs%twsavg(c,r))
                     else
!                        input_data(c+(r-1)*GRACEtwsobs%gracenc-180) = & 
                        input_data(c+(r-1)*GRACEtwsobs%gracenc-md_nc) = & 
                             (GRACEtwsobs%lwe_thickness(c,r,k)-&
                             GRACEtwsobs%twsavg(c,r))
                     endif
                  endif
               enddo
            enddo
            
           ! Scaling performed:
           else

            do r=1,GRACEtwsobs%gracenr
               do c=1,GRACEtwsobs%gracenc               
                  if(GRACEtwsobs%twsavg(c,r).ne.LDT_rc%udef.or.&
                       GRACEtwsobs%lwe_thickness(c,r,k).ne.32767) then 
                     ! Convert GRACE 0 to 360dg domain to -180 to 180E:
!                     if(c.le.180) then 
                     if(c.le.md_nc) then 
!                        input_data(c+(r-1)*GRACEtwsobs%gracenc+180) = & 
                        input_data(c+(r-1)*GRACEtwsobs%gracenc+md_nc) = & 
                             (GRACEtwsobs%lwe_thickness(c,r,k) &
                             *GRACEtwsobs%scalefactor(c,r)-&
                             GRACEtwsobs%twsavg(c,r))
                     else
!                        input_data(c+(r-1)*GRACEtwsobs%gracenc-180) = & 
                        input_data(c+(r-1)*GRACEtwsobs%gracenc-md_nc) = & 
                             (GRACEtwsobs%lwe_thickness(c,r,k) &
                             *GRACEtwsobs%scalefactor(c,r)-&
                             GRACEtwsobs%twsavg(c,r))
                     endif
                  endif
               enddo
            enddo
            
         endif  ! End of Scaling conditional  

         input_bitmap = .false. 
         do t=1,GRACEtwsobs%gracenc*GRACEtwsobs%gracenr
            if(input_data(t).ne.LDT_rc%udef) then 
               input_bitmap(t) = .true.
            endif
         enddo
         
         call neighbor_interp(LDT_rc%gridDesc(n,:),&
              input_bitmap, input_data, &
              output_bitmap, output_data, & 
              GRACEtwsobs%gracenc*GRACEtwsobs%gracenr, & 
              LDT_rc%lnc(n)*LDT_rc%lnr(n),&
              LDT_domain(n)%lat, LDT_domain(n)%lon,&
              GRACEtwsobs%n111,&
              LDT_rc%udef, iret)

         ! Account for the measurement ("merror") and  "leakage" errors:
         if(GRACEtwsobs%graceerrfile .ne. "none") then            
            
            do r=1,GRACEtwsobs%gracenr
               do c=1,GRACEtwsobs%gracenc               
                  if(GRACEtwsobs%merror(c,r).ne.LDT_rc%udef.or.&
                       GRACEtwsobs%lwe_thickness(c,r,k).ne.32767) then 
                     ! Convert GRACE 0 to 360dg domain to -180 to 180E:
!                     if(c.le.180) then 
!                        input_data(c+(r-1)*GRACEtwsobs%gracenc+180) = & 
                     if(c.le.md_nc) then 
                        input_data(c+(r-1)*GRACEtwsobs%gracenc+md_nc) = & 
                             GRACEtwsobs%merror(c,r)
                     else
!                        input_data(c+(r-1)*GRACEtwsobs%gracenc-180) = & 
                        input_data(c+(r-1)*GRACEtwsobs%gracenc-md_nc) = & 
                             GRACEtwsobs%merror(c,r)
                     endif
                  endif
               enddo
            enddo
            
            call neighbor_interp(LDT_rc%gridDesc(n,:),&
                 input_bitmap, input_data, &
                 output_err_bitmap, output_merr_data, & 
                 GRACEtwsobs%gracenc*GRACEtwsobs%gracenr, & 
                 LDT_rc%lnc(n)*LDT_rc%lnr(n),&
                 LDT_domain(n)%lat, LDT_domain(n)%lon,&
                 GRACEtwsobs%n111,&
                 LDT_rc%udef, iret)
            
            do r=1,GRACEtwsobs%gracenr
               do c=1,GRACEtwsobs%gracenc               
                  if(GRACEtwsobs%lerror(c,r).ne.LDT_rc%udef.or.&
                       GRACEtwsobs%lwe_thickness(c,r,k).ne.32767) then 
                     ! Convert GRACE 0 to 360dg domain to -180 to 180E:
!                     if(c.le.180) then 
!                        input_data(c+(r-1)*GRACEtwsobs%gracenc+180) = & 
                     if(c.le.md_nc) then 
                        input_data(c+(r-1)*GRACEtwsobs%gracenc+md_nc) = & 
                             GRACEtwsobs%lerror(c,r)
                     else
!                        input_data(c+(r-1)*GRACEtwsobs%gracenc-180) = & 
                        input_data(c+(r-1)*GRACEtwsobs%gracenc-md_nc) = & 
                             GRACEtwsobs%lerror(c,r)
                     endif
                  endif
               enddo
            enddo
            
            call neighbor_interp(LDT_rc%gridDesc(n,:),&
                          input_bitmap, input_data, &
                          output_err_bitmap, output_lerr_data, & 
                          GRACEtwsobs%gracenc*GRACEtwsobs%gracenr, & 
                          LDT_rc%lnc(n)*LDT_rc%lnr(n),&
                          LDT_domain(n)%lat, LDT_domain(n)%lon,&
                          GRACEtwsobs%n111,&
                          LDT_rc%udef, iret)
            
            
            ! Process and average GRACE to the basin-scale:
            if(GRACEtwsobs%process_basin_scale.eq.1) then 
               allocate(var_merr(GRACEtwsobs%basin_cat_max))
               allocate(var_lerr(GRACEtwsobs%basin_cat_max))
               allocate(nerr(GRACEtwsobs%basin_cat_max))
               allocate(basin_merr(GRACEtwsobs%basin_cat_max))
               allocate(basin_lerr(GRACEtwsobs%basin_cat_max))

               var_merr = 0
               var_lerr = 0. 
               nerr  = 0 
               betam = 300. ! km ~ measurement error decorrelation length
               betal = 100. ! km ~ leakage error decorrelation length

               ! All locations with a given category:
               do t=1,GRACEtwsobs%basin_cat_max
                  do r=1,LDT_rc%lnr(n)
                     do c=1,LDT_rc%lnc(n)
                        if(GRACEtwsobs%basin_cat(c,r).eq.t) then 
                           
                           lat1 = LDT_domain(n)%lat(c+(r-1)*LDT_rc%lnc(n))
                           lon1 = LDT_domain(n)%lon(c+(r-1)*LDT_rc%lnc(n))
                           nerr(t) = nerr(t)+1

                           do r1=1,LDT_rc%lnr(n)
                              do c1=1,LDT_rc%lnc(n)
                                 if(GRACEtwsobs%basin_cat(c1,r1).eq.t) then 
                                    
                                    lat2 = LDT_domain(n)%lat(c1+(r1-1)*LDT_rc%lnc(n))
                                    lon2 = LDT_domain(n)%lon(c1+(r1-1)*LDT_rc%lnc(n))
                                    
                                    ! lon, lat in degs, dist in km
                                    dist = sqrt(((lon1-lon2)*cos(lat1))**2.+&
                                         (lat1-lat2)**2.) * (LDT_CONST_PI/180.0) * 6371. 
                                    
                                    expdbm = exp(-(dist**2.)/(2.*betam**2.))
                                    expdbl = exp(-(dist**2.)/(2.*betal**2.))

                                    var_merr(t) = var_merr(t) + &
                                         output_merr_data(c+(r-1)*LDT_rc%lnc(n))* &
                                         output_merr_data(c1+(r1-1)*LDT_rc%lnc(n)) *&
                                         expdbm
                                    var_lerr(t) = var_lerr(t) + &
                                         output_lerr_data(c+(r-1)*LDT_rc%lnc(n)) * &
                                         output_lerr_data(c1+(r1-1)*LDT_rc%lnc(n)) * &
                                         expdbl
                                 endif
                              enddo
                           enddo

                        endif
                     enddo
                  enddo
                  
                  basin_merr(t) = sqrt(var_merr(t))/float(nerr(t))
                  basin_lerr(t) = sqrt(var_lerr(t))/float(nerr(t))
               enddo
               
               do r=1,LDT_rc%lnr(n)
                  do c=1,LDT_rc%lnc(n)
                     if(output_err_bitmap(c+(r-1)*LDT_rc%lnc(n))) then 
                        t = nint(GRACEtwsobs%basin_cat(c,r))
                        if(t.gt.0) then 
                           output_err_data(c+(r-1)*LDT_rc%lnc(n)) = sqrt(basin_merr(t)**2+&
                                basin_lerr(t)**2)
                        else
                           output_err_data(c+(r-1)*LDT_rc%lnc(n)) = LDT_rc%udef
                        endif
                     else
                        output_err_data(c+(r-1)*LDT_rc%lnc(n)) = LDT_rc%udef
                     endif
                  enddo
               enddo
               deallocate(var_merr)
               deallocate(var_lerr)
               deallocate(nerr)
               deallocate(basin_merr)
               deallocate(basin_lerr)

            ! Work on gridded GRACE obs:
            else
               
               do t=1,LDT_rc%lnc(n)*LDT_rc%lnr(n)
                  if(output_err_bitmap(t)) then 
                     output_err_data(t) = sqrt(output_merr_data(t)**2+&
                          output_lerr_data(t)**2)
                  else
                     output_err_data(t) = LDT_rc%udef
                  endif
               enddo
            endif
         ! When not accounting for error in the data:
         else
            output_err_data = LDT_rc%udef
         endif
         !
         ! Generate processed observations by adding the anomalies to the 
         ! LIS output data. 
         ! 
         if(GRACEtwsobs%process_basin_scale.eq.1) then 
            
            allocate(output_basin_data(GRACEtwsobs%basin_cat_max))
            allocate(output_basin_data_count(GRACEtwsobs%basin_cat_max))
            
            output_basin_data = 0.0
            output_basin_data_count = 0
            
            do r=1,LDT_rc%lnr(n)
               do c=1,LDT_rc%lnc(n)
                  if(GRACEtwsobs%basin_cat(c,r).ne.LDT_rc%udef) then 
                     cat_val = nint(GRACEtwsobs%basin_cat(c,r))
                     output_basin_data(cat_val) = &
                          output_basin_data(cat_val) + & 
                          output_data(c+(r-1)*LDT_rc%lnc(n))
                     output_basin_data_count(cat_val) = &
                          output_basin_data_count(cat_val) + 1 
                  endif
               enddo
            enddo
            
            do t=1,GRACEtwsobs%basin_cat_max
               output_basin_data(t) = &
                    output_basin_data(t)/&
                    float(output_basin_data_count(t))
            enddo
            
            do r=1,LDT_rc%lnr(n)
               do c=1,LDT_rc%lnc(n)
                  if(output_data(c+(r-1)*LDT_rc%lnc(n)).ne.LDT_rc%udef) then
                     if(GRACEtwsobs%basin_cat(c,r).ne.LDT_rc%udef) then
                        cat_val = nint(GRACEtwsobs%basin_cat(c,r))
                        output_data(c+(r-1)*LDT_rc%lnc(n)) = &
                             output_basin_data(cat_val)
                     endif
                  endif
               enddo
            enddo
            
            deallocate(output_basin_data)
            deallocate(output_basin_data_count)
         endif
            
         do r=1,LDT_rc%lnr(n)
            do c=1,LDT_rc%lnc(n)
               if(GRACEtwsobs%lisavg(c,r).gt.0.0 .and. & !ag
                    output_data(c+(r-1)*LDT_rc%lnc(n)).ne.LDT_rc%udef) then 
                  output_data(c+(r-1)*LDT_rc%lnc(n)) = &
                       GRACEtwsobs%lisavg(c,r) - GRACEtwsobs%swsavg(c,r,k) + & !ag
                       output_data(c+(r-1)*LDT_rc%lnc(n))*10.0 !to mm. 
               else
                  output_data(c+(r-1)*LDT_rc%lnc(n)) = &
                       LDT_rc%udef
               endif
            enddo
         enddo
         
         ! Write processed observation at the GRACE observation timestamp
         valid_data = .false. 
         do r=1,LDT_rc%lnr(n)
            do c=1,LDT_rc%lnc(n)
               if(output_data(c+(r-1)*LDT_rc%lnc(n)).ne.LDT_rc%udef) then
                  valid_data = .true.
               endif
            enddo
         enddo
         
         ! Write out processed GRACE output data:
         if(valid_data) then
            ftn = LDT_getNextUnitNumber()
            currTime = GRACEtwsobs%tvals(k)*24+GRACEtwsobs%reftime
            call LDT_julhr_date(currTime, yyyy,mm,dd,hh)
            
            write(unit=ftime,fmt='(i4.4,i2.2)') yyyy,mm
            filename = trim(LDT_rc%odir)//'/GRACE_obs_'//trim(ftime)//'.bin'
            
            write(LDT_logunit,*) "[INFO] writing processed GRACE obs "//&
                 trim(filename)
            open(ftn,file=trim(filename), form='unformatted')
            write(ftn) output_data
            write(ftn) output_err_data
            call LDT_releaseUnitNumber(ftn)
            
         endif

      endif
   enddo
endif
#endif

end subroutine readGRACEtwsObs

!BOP
!
! !ROUTINE: create_lsm_twsoutput_filename
! \label{create_lsm_twsoutput_filename}
!
! !INTERFACE:
subroutine create_lsm_twsoutput_filename(n, form, fname, odir, wstyle, wopt,mname)
! !USES:
   use LDT_coreMod,  only : LDT_rc
   use LDT_logMod
   use LDT_constantsMod, only : LDT_CONST_PATH_LEN

   implicit none 
! !ARGUMENTS:
   integer,   intent(IN)        :: n
   character(len=*)             :: fname
   character(len=*)             :: form
   character(len=*)             :: odir
   character(len=*)             :: wstyle
   character(len=*)             :: wopt
   character(len=*)             :: mname !ag (21Dec2017)
! 
! !DESCRIPTION:  
!  Create the file name for the output data files. It creates both the GSWP
!  style of output filenames and the standard LIS style. The convention used
!  in LIS creates a filename in a hierarchical style (output directory, 
!  model name, date, file extention)
!
!  2 level hierarchy
!  \begin{verbatim}
!   <output directory>/<model name>/LIS_HIST_<yyyymmddhhmnss>.<extension>
!  \end{verbatim}
!  3 level hierarchy
!  \begin{verbatim}
!   <output directory>/<model name>/<yyyymm>/LIS_HIST_<yyyymmddhhmnss>.<extension>
!  \end{verbatim}
!  4 level hierarchy
!  \begin{verbatim}
!   <output directory>/<model name>/<yyyy>/<yyyymm>/LIS_HIST_<yyyymmddhhmnss>.<extension>
!  \end{verbatim}
!  WMO convention
!  \begin{verbatim}
!   <output directory>/<AFWA Weather product style>
!  \end{verbatim}
!   A filename in the convention of weather products (such as): \newline
!   {\small
!   PS.AFWA\_SC.U\_DI.C\_DC.ANLYS\_GP.LIS\_GR.C0P25DEG\_AR.GLOBAL\_PA.03-HR-SUM\_DD.YYYYMMDD\_DT.HH00\_DF.GR1 \newline
!   }
!   where                             \newline
!    PS = Product source              \newline
!    SC = security classification     \newline
!    DI = distribution classification \newline
!    DC = data category               \newline
!    GP = generating process          \newline
!    GR = grid                        \newline
!    AR = area of data                \newline
!    PA = parameter                   \newline
!    DD = date                        \newline
!    DT = data time                   \newline
!    DF = data format                 \newline
! 
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest
!   \item [fname]
!     the created file name. 
!   \item [model\_name]
!    string describing the name of the model 
!   \item [writeint]
!    output writing interval  of the model
!   \item [style]
!    style option as described above
!  \end{description}
!EOP
   character(len=8)        :: date
   character(len=10)       :: time
   character(len=5)        :: zone
   integer, dimension(8)   :: values
   !character(len=20)       :: mname !ag (21Dec2017)
   character(len=10)       :: cdate
   character(len=14)       :: cdate1
   character(len=2)        :: fint
   character(len=10)       :: fres
   character(len=10)       :: fres2
   character(len=10)       :: fres3
   character*1             :: fres1(10)
   character(len=1)        :: fproj
   integer                 :: curr_mo = 0
   character(len=LDT_CONST_PATH_LEN) :: dname
   character(len=200), save :: out_fname
   integer                  :: i, c

   !ag (21Dec2017)
   !mname = 'SURFACEMODEL' 

   if(wstyle.eq."4 level hierarchy") then 
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4)') LDT_rc%yr
      dname = trim(dname)//trim(cdate)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2, i2.2)') LDT_rc%yr, LDT_rc%mo, LDT_rc%da
      dname = trim(dname)//trim(cdate)
      
      out_fname = trim(dname)//'/LIS_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( form )
      case ( "binary" )
         if(wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_lsm_twsoutput_filename -- '// &
              'Unrecognized output format')
         call LDT_endrun 
      endselect

   elseif(wstyle.eq."3 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      write(unit=cdate, fmt='(i4.4, i2.2)') LDT_rc%yr, LDT_rc%mo
      dname = trim(dname)//trim(cdate)//'/'

      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( form )
       case ("binary")
         if(wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
       case ("grib1")
         out_fname = trim(out_fname)//'.grb'
       case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
       case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
       case default
         call ldt_log_msg('ERR: create_lsm_twsoutput_filename -- '// &
              'Unrecognized form value')
         call LDT_endrun 
      endselect

   elseif(wstyle.eq."2 level hierarchy") then
      write(unit=cdate1, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
           LDT_rc%yr, LDT_rc%mo, &
           LDT_rc%da, LDT_rc%hr,LDT_rc%mn
      
      dname = trim(odir)//'/'
      dname = trim(dname)//trim(mname)//'/'
      
      out_fname = trim(dname)//'LIS_HIST_'//trim(cdate1)
      
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
      out_fname = trim(out_fname)//trim(cdate)
      
      select case ( form )
      case ("binary")
         if(wopt.eq."1d tilespace") then 
            out_fname = trim(out_fname)//'.ts4r'
         elseif(wopt.eq."2d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         elseif(wopt.eq."1d gridspace") then 
            out_fname = trim(out_fname)//'.gs4r'
         endif
      case ("grib1")
         out_fname = trim(out_fname)//'.grb'
      case ("netcdf")
         out_fname = trim(out_fname)//'.nc'
      case ("grib2")
         out_fname = trim(out_fname)//'.gr2'
      case default
         call ldt_log_msg('ERR: create_lsm_twsoutput_filename -- '// &
              'Unrecognized form value')
         call LDT_endrun 
      endselect

   elseif(wstyle.eq."WMO convention") then 
      write(LDT_logunit,*) '[WARN] WMO convention style not currently supported '
   endif

   fname = out_fname

 end subroutine create_lsm_twsoutput_filename


 subroutine transformLISoutToGRACEgrid(n,tws_inp,tws_out)

  use LDT_coreMod
  use GRACEtws_obsMod

  integer,    intent(in) :: n
  real                   :: tws_inp(GRACEtwsobs%nc*GRACEtwsobs%nr)
  real                   :: tws_out(LDT_rc%lnc(n),LDT_rc%lnr(n))
!
! !DESCRIPTION: 
!  This routine interpolates or upscales the input data to 
!  the LDT grid. If the input data is finer than the LDT
!  grid, the input data is upscaled. If the input data is
!  coarser, then it is interpolated to the LDT grid. 
!
!EOP
  integer         :: ios
  integer         :: c,r
  logical*1       :: tws_data_b(GRACEtwsobs%nc*GRACEtwsobs%nr)
  real            :: twsobs_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1       :: twsobs_b_ip(GRACEtwsobs%nc*GRACEtwsobs%nr)
  
   do r=1,GRACEtwsobs%nr
      do c=1, GRACEtwsobs%nc
         if(tws_inp(c+(r-1)*GRACEtwsobs%nc).ne.LDT_rc%udef) then 
            tws_data_b(c+(r-1)*GRACEtwsobs%nc) = .true. 
         else
            tws_data_b(c+(r-1)*GRACEtwsobs%nc) = .false.
         endif
      enddo
   enddo

   if(LDT_isLDTatAfinerResolution(n,GRACEtwsobs%datares)) then 

!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
      call bilinear_interp(LDT_rc%gridDesc(n,:),&
           tws_data_b, tws_inp, twsobs_b_ip, twsobs_ip, &
           GRACEtwsobs%nc*GRACEtwsobs%nr, &
           LDT_rc%lnc(n)*LDT_rc%lnr(n), &
           LDT_domain(n)%lat, LDT_domain(n)%lon,&
           GRACEtwsobs%w11, GRACEtwsobs%w12, &
           GRACEtwsobs%w21, GRACEtwsobs%w22, &
           GRACEtwsobs%n11, GRACEtwsobs%n12, &
           GRACEtwsobs%n21, GRACEtwsobs%n22, &
           LDT_rc%udef, ios)
   else
      call upscaleByAveraging(&
           GRACEtwsobs%nc*GRACEtwsobs%nr,&
           LDT_rc%lnc(n)*LDT_rc%lnr(n),LDT_rc%udef, &
           GRACEtwsobs%n11,tws_data_b, tws_inp, twsobs_b_ip,twsobs_ip)
      
   endif
   
   do r=1,LDT_rc%lnr(n)
      do c=1,LDT_rc%lnc(n)
         if(twsobs_b_ip(c+(r-1)*LDT_rc%lnc(n))) then 
            tws_out(c,r) = twsobs_ip(c+(r-1)*LDT_rc%lnc(n))
         else
            tws_out(c,r) = LDT_rc%udef
         endif
      enddo
   enddo
  
 end subroutine transformLISoutToGRACEgrid
