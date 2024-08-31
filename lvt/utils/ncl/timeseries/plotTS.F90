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
! !ROUTINE: plotTS
!  \label{plotTS}
!
program plotTS
! !USES: 
  use ESMF
!
! !DESCRIPTION: 
!  This program creates a time series plotting file for use with NCL
!  The usage is :
!   plotTS ts.config
! 
!  The ts.config file is structured as follows: 
!   
!   Output format for images:        eps # x11 or eps
!   LVT config file:                 lvt.config
!   LVT time series locations file:  TS_LOCATIONS.TXT
!   NCL output filename:         ts.ncl
!   Plot summary stats in the image: 0
!   data series 1 (model) label:     STEP0
!   data series 2 (obs) label:       STEP1
!   Number of xtics:                 5
!   Number of variables to plot:     3
!   Variable names:     Qle Qh Qg 
!   Variable units:     W/m2 W/m2 W/m2 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  20 Nov 11    Sujay Kumar  Initial Specification
! 
!EOP
  implicit none

  type(ESMF_Config)      :: ts_config
  character*50           :: model_name
  character*50           :: d1label, d2label
  integer                :: d1style, d2style
  character*100          :: tsconfigfile,ofile
  character*100          :: configfile,gfile,tsfile, statsodir
  character*100          :: rmsefile,biasfile,rcorrfile,mfile
  character*100          :: acorrfile, armsefile
  character*3            :: expcode
  type(ESMF_Config)      :: config_lvt
  integer                :: rc
  character*50           :: lsm
  integer                :: nstns
  character*100, allocatable :: varname(:)
  character*10,  allocatable :: vunits(:)
  character*100, allocatable :: stnname(:)
  character*16           :: labelx
  real, allocatable          :: labely(:)
  real,allocatable           :: labeldy(:)
  integer                :: i,k
  integer                :: wstats
  integer                :: nvars, j
  integer                :: v_index1, v_index2
  character*3            :: fid1, fid2
  real , allocatable         :: rmse(:,:)
  real,  allocatable         :: bias(:,:)
  real, allocatable          :: rcorr(:,:)
  real, allocatable          :: acorr(:,:)
  real, allocatable          :: armse(:,:)
  logical                :: file_exists,rmse_f_exists,bias_f_exists
  logical                :: rcorr_f_exists, armse_f_exists,acorr_f_exists
  character*10           :: ctmp1
  real   , allocatable       :: mv11(:),mv12(:),ov11(:),ov12(:)
  real,     allocatable      :: mv1(:),ov1(:)
  integer                :: yr,mo,da,hr,mn,ss
  logical, allocatable       :: noplot(:)
  real*8                 :: maxtime, mintime,time,dx 
  real                   :: gmt
  integer                :: doy
  integer                :: nxtics
  character*20           :: ctick
  character*20, allocatable  :: xtics(:)
  character*30           :: gformat
  integer                :: kk,ios

  ss = 0 

  k=iargc()
  if(k.lt.1) then
     print*,'Usage: '
     print*,'plotTS ts.config'
     stop
  else
     call getarg(1,tsconfigfile)

     ts_config = ESMF_ConfigCreate(rc=rc)
     call ESMF_ConfigLoadFile(ts_config,trim(tsconfigfile),rc=rc)

     call ESMF_ConfigGetAttribute(ts_config,gformat, &
          label="Output format for images:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,configfile, &
          label="LVT config file:",rc=rc)

     call ESMF_ConfigGetAttribute(ts_config,tsfile, &
          label="LVT time series locations file:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,gfile, &
          label="NCL output filename:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,wstats, &
          label="Plot summary stats in the image:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,nvars, &
          label="Number of variables to plot:",rc=rc)

     call ESMF_ConfigGetAttribute(ts_config,d1label, &
          label="data series 1 (model) label:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,d2label, &
          label="data series 2 (obs) label:",rc=rc)

!     call ESMF_ConfigGetAttribute(ts_config,d1style, &
!          label="data series 1 (model) style:",rc=rc)
!     call ESMF_ConfigGetAttribute(ts_config,d2style, &
!          label="data series 2 (obs) style:",rc=rc)

     call ESMF_ConfigGetAttribute(ts_config,nxtics, &
          label="Number of xtics:",rc=rc)
     
     allocate(xtics(nxtics))
     allocate(varname(nvars))
     allocate(vunits(nvars))

     call ESMF_ConfigFindLabel(ts_config,"Variable names:",rc=rc)
     do i=1,nvars
        call ESMF_ConfigGetAttribute(ts_config,varname(i),rc=rc)
     enddo
     call ESMF_ConfigFindLabel(ts_config,"Variable units:",rc=rc)
     do i=1,nvars
        call ESMF_ConfigGetAttribute(ts_config,vunits(i),rc=rc)
     enddo
     
     config_lvt = ESMF_ConfigCreate(rc=rc)
     call ESMF_ConfigLoadFile(config_lvt,trim(configfile),rc=rc)
     
!     call ESMF_ConfigGetAttribute(Config_lvt,lsm,label="Land surface model:",rc=rc)  
!     call ESMF_ConfigGetAttribute(config_lvt,expcode,label="Experiment code:",&
!          rc=rc)
     call ESMF_ConfigGetAttribute(config_lvt,statsodir,&
          label="Stats output directory:",&
          rc=rc)
     model_name = trim(lsm)
!     if(lsm.eq.1) then 
!        model_name = 'NOAH271'
!     elseif(lsm.eq.2) then 
!        model_name = 'CLM'
!     elseif(lsm.eq.4) then 
!        model_name = 'MOS'
!     elseif(lsm.eq.5) then 
!        model_name = 'HYSSIB'
!     elseif(lsm.eq.7) then 
!        model_name = 'CLSM'
!     elseif(lsm.eq.14) then 
!        model_name = 'TESSEL'
!     elseif(lsm.eq.21) then 
!        model_name = 'NOAH31'
!     elseif(lsm.eq.22) then 
!        model_name = 'NOAH32'
!     endif
     open(100,file=trim(tsfile),form='formatted')
     read(100,*)
     read(100,*) nstns
     read(100,*)
     read(100,*)
     read(100,*)
     allocate(stnname(nstns))
     do i=1,nstns
        read(100,*) stnname(i)
        read(100,*)
     enddo
     close(100)

     if(wstats.eq.1) then 
        rmsefile = trim(statsodir)//'/RMSE_SUMMARY_STATS.dat'
        
        inquire(file=rmsefile,exist=rmse_f_exists)
        if(rmse_f_exists) then 
           allocate(rmse(nstns,nvars))
           open(100,file=trim(rmsefile),form='formatted')
           do j=1,nvars
              read(100,*) 
              read(100,*)
              read(100,*)
              read(100,*) !all
              
              do i=1,nstns
                 read(100,fmt='(a10,E14.3)') ctmp1,rmse(i,j)
              enddo
           enddo
           close(100)
        endif
        
        biasfile = trim(statsodir)//'/BIAS_SUMMARY_STATS.dat'
        
        inquire(file=biasfile,exist=bias_f_exists)
        if(bias_f_exists) then 
           allocate(bias(nstns,nvars))
           open(100,file=trim(biasfile),form='formatted')
           do j=1,nvars
              read(100,*) 
              read(100,*)
              read(100,*)
              read(100,*) !all
              
              do i=1,nstns
                 read(100,fmt='(a10,E14.3)') ctmp1,bias(i,j)
              enddo
           enddo
           close(100)
        endif
        
        rcorrfile = trim(statsodir)//'/RCORR_SUMMARY_STATS.dat'
        
        inquire(file=rcorrfile,exist=rcorr_f_exists)
        if(rcorr_f_exists) then 
           allocate(rcorr(nstns,nvars))
           open(100,file=trim(rcorrfile),form='formatted')
           do j=1,nvars
              read(100,*) 
              read(100,*)
              read(100,*)
              read(100,*) !all
              
              do i=1,nstns
                 read(100,fmt='(a10,E14.3)') ctmp1,rcorr(i,j)
              enddo
           enddo
           close(100)
        endif
        acorrfile = trim(statsodir)//'/ACORR_SUMMARY_STATS.dat'
        
        inquire(file=acorrfile,exist=acorr_f_exists)
        if(acorr_f_exists) then 
           allocate(acorr(nstns,nvars))
           open(100,file=trim(acorrfile),form='formatted')
           do j=1,nvars
              read(100,*) 
              read(100,*)
              read(100,*)
              read(100,*) !all
              
              do i=1,nstns
                 read(100,fmt='(a10,E14.3)') ctmp1,acorr(i,j)
              enddo
           enddo
           close(100)
        endif

        armsefile = trim(statsodir)//'/ARMSE_SUMMARY_STATS.dat'
        
        inquire(file=armsefile,exist=armse_f_exists)
        if(armse_f_exists) then 
           allocate(armse(nstns,nvars))
           open(100,file=trim(armsefile),form='formatted')
           do j=1,nvars
              read(100,*) 
              read(100,*)
              read(100,*)
              read(100,*) !all
              
              do i=1,nstns
                 read(100,fmt='(a10,E14.3)') ctmp1,armse(i,j)
              enddo
           enddo
           close(100)
        endif
     endif

     allocate(noplot(nvars))
     noplot = .false.
     open(100,file=trim(gfile),form='formatted') 
     do i=1,nstns
        mfile = trim(statsodir)//'/MEAN_'//&
             trim(stnname(i))//'.dat'
        print*, 'Processing ... ',trim(mfile)
        v_index1 = 17
        do j=1,nvars
           if(.not.noplot(j)) then 

              write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"'
              write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"'
              write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"'
              write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"'
              write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"'

              write(100,*) 
              write(100,*) 'begin'
              ofile = trim(stnname(i))//'_'//&
                      trim(varname(j))//'_ts'
              if(gformat.eq."x11") then 
                 write(100,*) 'wks = gsn_open_wks("x11","'//&
                      trim(ofile)//'") '
              elseif(gformat.eq."eps") then 
                 write(100,*) 'wks = gsn_open_wks("eps","'//&
                      trim(ofile)//'") '
              endif

              write(100,*) 'fname = "'//trim(statsodir)//&
                   '/MEAN_'//trim(stnname(i))//'.dat"'

              write(100,*) 'data1 = asciiread(fname,-1,"string")'
  
              write(100,*)  'yr = stringtoint(str_get_cols(data1,0,3))'
              write(100,*)  'mo = stringtoint(str_get_cols(data1,5,6))'
              write(100,*)  'da = stringtoint(str_get_cols(data1,8,9))'
              write(100,*)  'hr = stringtoint(str_get_cols(data1,11,12))'
              write(100,*)  'mn = stringtoint(str_get_cols(data1,14,15))'
              write(100,*)  'ss = mn'
              write(100,*) 'units = "seconds since 1950-01-01 00:00:00"'
              write(100,*) 'time = cd_inv_calendar(yr,mo,da,hr,mn,ss,units,0)'
              write(100,*) 

              v_index2 = v_index1 + 14 -1
              
              write(unit=fid1,fmt='(i3.3)') v_index1
              write(unit=fid2,fmt='(i3.3)') v_index2

              write(100,*)  'var1 = stringtofloat(str_get_cols(data1,'//&
                   trim(fid1)//','//trim(fid2)//'))'

              v_index1 = v_index1+ 6*14
              v_index2 = v_index1 + 14 -1
              write(unit=fid1,fmt='(i3.3)') v_index1
              write(unit=fid2,fmt='(i3.3)') v_index2
              write(100,*)  'var2 = stringtofloat(str_get_cols(data1,'//&
                   trim(fid1)//','//trim(fid2)//'))'

              write(unit=ctick,fmt='(i3.3)') nxtics

              write(100,*) 'ntime = dimsizes(time)'
              write(100,*) 'nticks = '//trim(ctick)
              write(100,*) 'time1 = new((/nticks+1/),double)'
              write(100,*) 'dx = tointeger(ceil(1.0*ntime)/(1.0*nticks))'
              write(100,*) 'k = 0 '
              write(100,*) 'kk = 0 '
              write(100,*) 'do while(k.lt.ntime)'
              write(100,*) ' time1(kk) = time(k)'
              write(100,*) ' kk = kk + 1'
              write(100,*) '  k = k + dx'
              write(100,*) 'end do'
              write(100,*) 'stime = ut_string(time1,"%Y/%N")'
              
              write(100,*) 'res = True'
              write(100,*)  'res@gsnFrame        = False'
              write(100,*)  'res@gsnDraw        = False'
              write(100,*) 'res@vpWidthF  = 0.8'
              write(100,*) 'res@vpHeightF = 0.2'
              write(100,*) 'res@xyLineThicknesses = (/2.0,2.0/)'               
              write(100,*) 'res@xyLineColors      = (/"blue","red"/)'
!              write(100,*) 'res@xyMarkLineMode = "MarkLines"'
!              write(100,*) 'res@xyMarkers      = (/6/)'
!              write(100,*) 'res@xyMarkerColors = (/"blue"/)'
              write(100,*) 'res@tmXBMode = "Explicit"'
              write(100,*) 'res@tmXBValues = time1'
              write(100,*) 'res@tmXBLabels = stime'
              write(100,*) 'res@tmXBLabelAngleF = 70'

              write(100,*) 'res@pmLegendDisplayMode    = "Always"'
              write(100,*) 'res@pmLegendSide           = "Top"'
              write(100,*) 'res@pmLegendParallelPosF   = .9'
              write(100,*) 'res@pmLegendOrthogonalPosF = -0.4'
              
              write(100,*) 'res@pmLegendWidthF         = 0.05'
              write(100,*) 'res@pmLegendHeightF        = 0.05'
              write(100,*) 'res@lgPerimOn              = True'
              write(100,*) 'res@lgLabelFontHeightF     = .02'
              
              write(100,*) 'res@xyExplicitLegendLabels = (/"'&
                   //trim(d1label)//'","'//trim(d2label)//'"/)'
              
              write(100,*) 'res@tiYAxisString = "'//&
                   trim(varname(j))//'('//trim(vunits(j))//')"'
!              ;  res@trYMinF = 0.0       			
!              ;  res@trYMaxF = 10.0	
              write(100,*) 'var = new((/2,dimsizes(var1)/),float)'
              write(100,*) 'var(0,:) = var1'
              write(100,*) 'var(1,:) = var2'

              write(100,*) 'var@_FillValue = -9999.0'
              write(100,*) 'plot = gsn_csm_xy(wks,time,var,res)'
              write(100,*) 'draw(plot)'

              write(100,*)
              write(100,*)  'txres               = True'
              write(100,*)  'txres@txFontHeightF = 0.010'
              write(100,*)  'txres@txPerimOn = True'
              write(100,*)  'txres@txBackgroundFillColor = "Cyan"'
              write(100,*)

              write(100,*) 'minx = min(time)'
              write(100,*) 'maxx = max(time)'
              write(100,*)
              write(100,*) 'miny = min(var1)'
              write(100,*) 'maxy = max(var1)'
              write(100,*)

              if(wstats.eq.1) then
                 if(rmse_f_exists.and.rmse(i,j).ne.-0.100E+05) then 
                    write(100,*) 'v1_loc = minx + 2*(maxx-minx)/10'
                    write(100,'(a28)',advance='no') 'gsn_text(wks,plot,"RMSE: "+'
                    write(100,'(F14.3)',advance='no') rmse(i,j)
                    write(100,*) ',v1_loc,maxy,txres)'
                    write(100,*) ''
                 endif
                 if(bias_f_exists.and.bias(i,j).ne.-0.100E+05) then
                    write(100,*) 'v2_loc = maxy - 2*(maxy-miny)/10'
                    write(100,'(a28)',advance='no') 'gsn_text(wks,plot,"BIAS: "+'
                    write(100,'(F14.3)',advance='no') bias(i,j)
                    write(100,*) ',v1_loc,v2_loc,txres)'
                    write(100,*) 
                 endif
                 if(rcorr_f_exists.and.rcorr(i,j).ne.-0.100E+05) then 
                    write(100,*) 'v3_loc = maxy - 4*(maxy-miny)/10'
                    write(100,'(a25)',advance='no') 'gsn_text(wks,plot,"R: "+'
                    write(100,'(F14.3)',advance='no') rcorr(i,j)
                    write(100,*) ',v1_loc,v3_loc,txres)'
                    write(100,*) 
                 endif
                 if(acorr_f_exists.and.acorr(i,j).ne.-0.100E+05) then 
                    write(100,*) 'v3_loc = maxy - 6*(maxy-miny)/10'
                    write(100,'(a32)',advance='no') 'gsn_text(wks,plot,"Anomaly R: "+'
                    write(100,'(F14.3)',advance='no') acorr(i,j)
                    write(100,*) ',v1_loc,v3_loc,txres)'
                    write(100,*) 
                 endif
                 if(armse_f_exists.and.armse(i,j).ne.-0.100E+05) then 
                    write(100,*) 'v3_loc = maxy - 8*(maxy-miny)/10'
                    write(100,'(a35)',advance='no') 'gsn_text(wks,plot,"Anomaly RMSE: "+'
                    write(100,'(F14.3)',advance='no') armse(i,j)
                    write(100,*) ',v1_loc,v3_loc,txres)'
                    write(100,*) 
                 endif
              endif
#if 0 
              if(wstats.eq.1) then
                 if(rmse_f_exists) then 
                    write(100,*) 'v1_loc = minx + 2*(maxx-minx)/10'
                    write(100,'(a28)',advance='no') 'gsn_text(wks,plot,"RMSE: "+'
                    write(100,'(F14.3)',advance='no') rmse(i,j)
                    write(100,*) ',v1_loc,maxy,txres)'
                    write(100,*) ''
                 endif
                 if(bias_f_exists) then
                    write(100,*) 'v1_loc = v1_loc + (maxx-minx)/nticks)'
                    write(100,'(a28)',advance='no') 'gsn_text(wks,plot,"BIAS: "+'
                    write(100,'(F14.3)',advance='no') bias(i,j)
                    write(100,*) ',v1_loc,maxy,txres)'
                    write(100,*) 
                 endif
                 if(rcorr_f_exists) then 
                    write(100,*) 'v1_loc = v1_loc + (maxx-minx)/nticks)'
                    write(100,'(a25)',advance='no') 'gsn_text(wks,plot,"R: "+'
                    write(100,'(F14.3)',advance='no') rcorr(i,j)
                    write(100,*) ',v1_loc,maxy,txres)'
                    write(100,*) 
                 endif
                 if(acorr_f_exists) then 
                    write(100,*) 'v1_loc = v1_loc + (maxx-minx)/nticks)'
                    write(100,'(a32)',advance='no') 'gsn_text(wks,plot,"Anomaly R: "+'
                    write(100,'(F14.3)',advance='no') acorr(i,j)
                    write(100,*) ',v1_loc,maxy,txres)'
                    write(100,*) 
                 endif
                 if(armse_f_exists) then 
                    write(100,*) 'v1_loc = v1_loc + (maxx-minx)/nticks)'
                    write(100,'(a35)',advance='no') 'gsn_text(wks,plot,"Anomaly RMSE: "+'
                    write(100,'(F14.3)',advance='no') armse(i,j)
                    write(100,*) ',v1_loc,maxy,txres)'
                    write(100,*) 
                 endif
              endif
#endif
              write(100,*) 'frame(wks)'
              write(100,*) 'end'
              write(100,*) 
              write(100,*) 'delete(plot)'
              write(100,*) 'delete(wks)'
           endif
           v_index1 = v_index1+ 6*14
           v_index2 = v_index1 + 14 -1
        enddo
     enddo
  endif
end program plotTS

  subroutine time2date(time,doy,gmt,yr,mo,da,hr,mn)

    implicit none
! !ARGUMENTS: 
    integer :: yr,mo,da,hr,mn,ss,doy
    real*8  :: time
    real    :: gmt
! !DESCRIPTION:
! 
!  determines the value of the year, month, day of month, hour of 
!  the day, minute and second based on the specified time. This
!  method is the inverse of date2time
! 
!  The arguments are: 
!  \begin{description}
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of the month
!  \item[hr]
!   hour of day
!  \item[mn]
!   minute
!  \item[ss]
!   second
!  \item[time]
!   lvt time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    real*8 :: tmp
    integer :: yrdays,days(13)
    data days /31,28,31,30,31,30,31,31,30,31,30,31,30/
    
    yr  = dint(time)
    tmp =     (time) 
    
    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then             !correct for y2k
       yrdays=366                  
    else
       yrdays=365
    endif
    if (yrdays.eq.366) then
       days(2)=29
    else
       days(2)=28
    endif
    
    doy  = dint((tmp-yr)*float(yrdays))+1 
    tmp =      ((tmp-yr)*float(yrdays))+1 
    hr  = nint((tmp-doy)*24.d0) 
    tmp =     ((tmp-doy)*24.d0) 
    
    mn  = dint((tmp-hr)*60.d0) 
    tmp =     ((tmp-hr)*60.d0) 
    
    ss  = dint((tmp-mn)*60.d0) 
    mo=1
    do while (doy.gt.0)
       doy=doy-days(mo)
       mo=mo+1
    enddo
    mo=mo-1
    da=doy+days(mo)
    
    gmt=(((float(ss)/60.0)+float(mn))/60.0)+float(hr)
    
    if(gmt.eq.24) then
       hr = 0 
       gmt=0
       da=da+1
       if (da.gt.days(mo)) then
          da=1
          mo=mo+1
          if (mo.gt.12) then
             mo=1
             yr=yr+1
          endif
       endif
    endif
    return
  end subroutine time2date

  subroutine date2time(time,doy,gmt,yr,mo,da,hr,mn,ss)

    implicit none
! !ARGUMENTS:
    integer :: yr,mo,da,hr,mn,ss, doy
    real*8  :: time
    real    :: gmt
! !DESCRIPTION:
! 
!  determines the time, time in GMT, and the day of the year
!  based on the value of year, month, day of month, hour of 
!  the day, minute and second. This method is the inverse of 
!  time2date
! 
!  The arguments are: 
!  \begin{description}
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of the month
!  \item[hr]
!   hour of day
!  \item[mn]
!   minute
!  \item[ss]
!   second
!  \item[time]
!   lvt time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    integer :: yrdays,days(13),k
    data days /31,28,31,30,31,30,31,31,30,31,30,31,30/

    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &     !correct for leap year
         .or.(mod(yr,400).eq.0))then             !correct for y2k
       yrdays=366                  
    else
       yrdays=365
    endif
    
    doy=0
    do k=1,(mo-1)
       doy=doy+days(k)
    enddo
    doy=doy+da
    
    if(yrdays.eq.366.and.mo.gt.2)doy=doy+1
    
    time=(float(yr)+((((((float(ss)/60.d0)+float(mn))/60.d0)+ & 
         float(hr))/24.d0)+float(doy-1))/float(yrdays))
    
    gmt=( ( (float(ss)/60.0) +float(mn)) /60.0)+float(hr)
    return
  end subroutine date2time
