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
! !ROUTINE: readsimGRACEJPLObs
! \label{readsimGRACEJPLObs}
! 
! !REVISION HISTORY: 
!  24 Feb 2015: Sujay Kumar, Initial Specification
!  12 Dec 2019: Eric Kemp, Added interface block for calling
!               create_lsm_output_filename.
! !INTERFACE: 
subroutine readsimGRACEJPLObs(n)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_historyMod
  use LDT_DAobsDataMod
  use LDT_timeMgrMod
  use simGRACEJPL_obsMod

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
!  TWS outputs from LIS is expected to be in mm
!EOP 
  character(len=LDT_CONST_PATH_LEN)         :: fname,filename
  integer               :: c,r,k,kk,t,iret
  integer               :: ftn
  integer               :: refyr,refmo,refda,refhr,refmn,refss
  integer               :: yr,mo,da,hr
  integer               :: timeId, tId, tbId, lweId, scaleId,errId,err1Id
  real                  :: tws_data(LDT_rc%lnc(n),LDT_rc%lnr(n))
  real                  :: output_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                  :: output_err_data(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1             :: output_bitmap(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical*1             :: output_err_bitmap(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                  :: input_data(simGRACEJPLobs%gracenc*simGRACEJPLobs%gracenr)
  logical*1             :: input_bitmap(simGRACEJPLobs%gracenc*simGRACEJPLobs%gracenr)
  integer               :: npts(simGRACEJPLobs%gracenc,simGRACEJPLobs%gracenr)
  logical               :: file_exists
  real                  :: dt
  integer               :: currTime,yyyy,mm,dd,hh
  character*10          :: ftime
  integer               :: status,sec
  logical               :: valid_data
  integer                 :: dt1
  type(ESMF_TimeInterval) :: dayInt,deltat
  type(ESMF_Time)         :: ctime,reftime

  !EMK...Add interface block for create_lsm_output_filename.  This is required
  !per the Fortran 90 standard since the subroutine has optional arguments but
  !is not defined in a module.  TODO:  Put the subroutine in a module.
  interface
     subroutine create_lsm_output_filename(n, form, fname, odir, wstyle, wopt,&
          run_dd, map_proj, security_class,   &
          distribution_class, data_category,  &
          area_of_data, write_interval)
        integer,   intent(IN)        :: n
        character(len=*)             :: fname
        character(len=*)             :: form
        character(len=*)             :: odir
        character(len=*)             :: wstyle
        character(len=*)             :: wopt
        real, dimension(8), optional :: run_dd
        character(len=*),   optional :: map_proj
        character(len=*),   optional :: security_class
        character(len=*),   optional :: distribution_class
        character(len=*),   optional :: data_category
        character(len=*),   optional :: area_of_data
        character(len=*),   optional :: write_interval
     end subroutine create_lsm_output_filename
  end interface


!
! At initial time, read the data into memory. At all other times, 
! simply index into the correct temporal location of the data. 
!
  output_err_data = LDT_rc%udef

  if(simGRACEJPLobs%startMode) then 
     simGRACEJPLobs%startMode = .false. 

     allocate(simGRACEJPLobs%lwe_thickness(simGRACEJPLobs%gracenc, &
          simGRACEJPLobs%gracenr,simGRACEJPLobs%tdims))
     allocate(simGRACEJPLobs%twsavg(simGRACEJPLobs%gracenc, &
          simGRACEJPLobs%gracenr))
     
     refyr = 2006
     refmo = 1
     refda = 1
     refhr = 0 
     refmn = 0 
     refss = 0 

     call ESMF_TimeSet(reftime,yy=refyr,mm=refmo,&
          dd=refda,h=refhr,m=refmn,s=refss,rc=status)
     call ESMF_TimeIntervalSet(dayInt,d=1,rc=status)
     
     cTime = refTime

     kk = 1
     do while(kk.le.simGRACEJPLobs%tdims)
        
        deltat =  (cTime-refTime)        
        call ESMF_TimeIntervalGet(deltat,d=dt1,s=sec, rc=status)

        if(floor(simGRACEJPLobs%tvals(kk)).eq.dt1) then 

           call ESMF_TimeGet(cTime,yy=refyr,mm=refmo,&
                dd=refda,h=refhr,m=refmn,s=refss,rc=status)
           
           call createSimGRACEJPLfilename(simGRACEJPLobs%gracedir,&
                refyr, refmo, refda, simGRACEJPLobs%config, filename)
            
           inquire(file=trim(filename),exist=file_exists) 

           if(file_exists) then 
              write(LDT_logunit,*) 'Reading simulated GRACE data '//&
                   trim(filename)
              ftn = LDT_getNextUnitNumber()
              open(ftn,file=filename,form='unformatted')
              read(ftn) simGRACEJPLobs%lwe_thickness(:,:,kk)
              call LDT_releaseUnitNumber(ftn)
!           else
!              write(LDT_logunit,*) 'simulated GRACE file '&
!                   //trim(filename)//& 
!                   'does not exist....'
!              call LDT_endrun()
           endif
           kk= kk + 1
        endif
        cTime = cTime + dayInt
     enddo

     npts = 0 
     simGRACEJPLobs%twsavg = 0 
     
     do k=1,simGRACEJPLobs%tdims
        do r=1,simGRACEJPLobs%gracenr

           do c=1,simGRACEJPLobs%gracenc
              if(simGRACEJPLobs%lwe_thickness(c,r,k).ne.32767.0) then 
                 simGRACEJPLobs%twsavg(c,r) = simGRACEJPLobs%twsavg(c,r) + & 
                      simGRACEJPLobs%lwe_thickness(c,r,k)
                 npts(c,r) = npts(c,r) + 1
              endif
           enddo
        enddo                  
     enddo
     do r=1,simGRACEJPLobs%gracenr
        do c=1,simGRACEJPLobs%gracenc
           if(npts(c,r).gt.0) then !convert to mm. 
              simGRACEJPLobs%twsavg(c,r) = simGRACEJPLobs%twsavg(c,r)&
                   /npts(c,r)
           else
              simGRACEJPLobs%twsavg(c,r) = LDT_rc%udef
           endif
        enddo
     enddo
     write(LDT_logunit,*) 'Finished reading GRACE data '//&
          trim(filename)
  endif

  
  if(LDT_rc%pass_id.eq.1) then ! read during the first pass for averaging
     
     call create_lsm_output_filename(simGRACEJPLobs%nest, &
          simGRACEJPLobs%format,&
          fname, simGRACEJPLobs%odir, simGRACEJPLobs%wstyle, &
          simGRACEJPLobs%wopt)
     
     !average only between baseline years. 
     if(LDT_rc%yr.ge.simGRACEJPLobs%b_syr.and.&
          LDT_rc%yr.le.simGRACEJPLobs%b_eyr) then
        inquire(file=trim(fname),exist=file_exists)
        
        if(file_exists) then 
           write(LDT_logunit,*) 'reading LSM output ',trim(fname)
        
           if(simGRACEJPLobs%format.eq."netcdf") then 
              
              iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
              call LDT_verify(iret, 'Error opening file '//trim(fname))
              
              call LDT_readLISSingleNetcdfVar(n,ftn, "TWS_tavg",&
                   1,tws_data)
              iret = nf90_close(ftn)
              call LDT_verify(iret,'Error in nf90_close')

              do r=1,LDT_rc%lnr(n)
                 do c=1,LDT_rc%lnc(n)
                    if(tws_data(c,r).ne.LDT_rc%udef) then 
                       simGRACEJPLobs%lisavg(c,r) = &
                            simGRACEJPLobs%lisavg(c,r) + tws_data(c,r)
                       simGRACEJPLobs%nlisavg(c,r) = &
                            simGRACEJPLobs%nlisavg(c,r) + 1
                    endif
                 enddo
              enddo
           endif
        else
           write(LDT_logunit,*) 'LIS file '//trim(fname)
           write(LDT_logunit,*) 'not found. Program stopping....'
           call LDT_endrun()
        endif
     endif
     ! At the end of the first cycle. 
!
     if(LDT_rc%pass_id.eq.1.and.LDT_rc%endtime.eq.1) then 
        do r=1,LDT_rc%lnr(n)
           do c=1,LDT_rc%lnc(n)
              if(simGRACEJPLobs%nlisavg(c,r).gt.0) then 
                 simGRACEJPLobs%lisavg(c,r) = &
                      simGRACEJPLobs%lisavg(c,r)/&
                      simGRACEJPLobs%nlisavg(c,r) 
              else
                 simGRACEJPLobs%lisavg(c,r) = LDT_rc%udef
              endif
           enddo
        enddo
     endif
  elseif(LDT_rc%pass_id.eq.2) then 
     
     call LDT_get_julhr(LDT_rc%yr,LDT_rc%mo,LDT_rc%da,&
          LDT_rc%hr,LDT_rc%mn,LDT_rc%ss,currTime)
     
     dt =  float((currTime-simGRACEJPLobs%refTime))/24.0
     
     do k=1,simGRACEJPLobs%tdims
        if(floor(simGRACEJPLobs%tvals(k)).eq.dt) then 
           !  Interpolate the anomaly map to the LIS output grid
           ! 
           input_data = LDT_rc%udef
	   
           do r=1,simGRACEJPLobs%gracenr
              do c=1,simGRACEJPLobs%gracenc               
                 if(simGRACEJPLobs%twsavg(c,r).ne.LDT_rc%udef.or.&
                      simGRACEJPLobs%lwe_thickness(c,r,k).ne.32767) then 
                    input_data(c+(r-1)*simGRACEJPLobs%gracenc) = & 
                         (simGRACEJPLobs%lwe_thickness(c,r,k) - & 
                         simGRACEJPLobs%twsavg(c,r))
                 endif
              enddo
           enddo
            
           input_bitmap = .false. 
           do t=1,simGRACEJPLobs%gracenc*simGRACEJPLobs%gracenr
              if(input_data(t).ne.LDT_rc%udef) then 
                 input_bitmap(t) = .true.
              endif
           enddo
         
           call bilinear_interp(LDT_rc%gridDesc(n,:),&
                input_bitmap, input_data, &
                output_bitmap, output_data, & 
                simGRACEJPLobs%gracenc*simGRACEJPLobs%gracenr, & 
                LDT_rc%lnc(n)*LDT_rc%lnr(n),&
                LDT_domain(n)%lat, LDT_domain(n)%lon,&
                simGRACEJPLobs%w11, simGRACEJPLobs%w12,&
                simGRACEJPLobs%w21,simGRACEJPLobs%w22,&
                simGRACEJPLobs%n11,simGRACEJPLobs%n12,&
                simGRACEJPLobs%n21,simGRACEJPLobs%n22,&
                LDT_rc%udef, iret)                       
!
! Generate processed observations by adding the anomalies to the 
! LIS output data. 
! 
            do r=1,LDT_rc%lnr(n)
               do c=1,LDT_rc%lnc(n)
                  if(simGRACEJPLobs%lisavg(c,r).ne.LDT_rc%udef.and.&
                       output_data(c+(r-1)*LDT_rc%lnc(n)).ne.LDT_rc%udef) then 
                     output_data(c+(r-1)*LDT_rc%lnc(n)) = &
                          simGRACEJPLobs%lisavg(c,r) + &
                          output_data(c+(r-1)*LDT_rc%lnc(n))*10.0 !to mm. 
                  else
                     output_data(c+(r-1)*LDT_rc%lnc(n)) = &
                          LDT_rc%udef
                  endif
               enddo
            enddo
            
!Write processed observation at the GRACE observation timestamp
            valid_data = .false. 
            do r=1,LDT_rc%lnr(n)
               do c=1,LDT_rc%lnc(n)
                  if(output_data(c+(r-1)*LDT_rc%lnc(n)).ne.LDT_rc%udef) then
                     valid_data = .true.
                  endif
               enddo
            enddo
            
            if(valid_data) then 
               ftn = LDT_getNextUnitNumber()
               currTime = simGRACEJPLobs%tvals(k)*24+simGRACEJPLobs%reftime
               call LDT_julhr_date(currTime, yyyy,mm,dd,hh)

               if((simGRACEJPLobs%config.eq."default").or.&
                    (simGRACEJPLobs%config.eq."follow-on")) then 
                  write(unit=ftime,fmt='(i4.4,i2.2)') yyyy,mm
                  filename = trim(LDT_rc%odir)//'/GRACE_obs_'//&
                       trim(ftime)//'.bin'
               elseif(simGRACEJPLobs%config.eq."GRACE-2") then 
                  write(unit=ftime,fmt='(i4.4,i2.2,i2.2)') yyyy,mm,dd
                  filename = trim(LDT_rc%odir)//'/GRACE_obs_'//&
                       trim(ftime)//'.bin'
               endif
               write(LDT_logunit,*) "writing processed GRACE obs "//&
                    trim(filename)
               open(ftn,file=trim(filename), form='unformatted')
               write(ftn) output_data
               write(ftn) output_err_data
               call LDT_releaseUnitNumber(ftn)
            endif
         endif
      enddo
   endif

end subroutine readsimGRACEJPLObs

subroutine createSimGRACEJPLfilename(odir,&
     yr, mo, da, cnfg, gracefile)

  character(len=*), intent(in) :: odir
  integer         , intent(in) :: yr
  integer         , intent(in) :: mo
  integer         , intent(in) :: da
  character(len=*)             :: cnfg
  character(len=*)             :: gracefile

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  character (len=8) :: fextn

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da

  fextn = '12.bin'
  gracefile = trim(odir)//'/simGRACE_'//trim(fyr)//trim(fmo)//trim(fda)//&
       trim(fextn)
       
end subroutine createSimGRACEJPLfilename
