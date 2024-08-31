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
! !ROUTINE: plotScatter
!  \label{plotScatter}
!
program plotScatter
! !USES: 
  use ESMF
!
! !DESCRIPTION: 
!  This program creates a script for scatter plot for NCL
!  The usage is :
!   plotScatter scatter.config
! 
!  The scatter.config file is structured as follows: 
!   
!   Output format for images:        1   #1 - x11, 2-eps
!   LVT config file:                 lvt.config
!   LVT time series locations file:  TS_LOCATIONS.TXT
!   gnuplot output filename:         ts.ncl
!   Plot summary stats in the image: 0
!   data series 1 (model) label:     STEP0
!   data series 2 (obs) label:       STEP1
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
  character*100          :: tsconfigfile
  character*100          :: ofile
  character*50           :: model_name
  character*50           :: d1label, d2label
  character*100          :: configfile,gfile,tsfile, statsodir
  character*100          :: rmsefile,biasfile,rcorrfile,mfile
  character*3            :: expcode
  type(ESMF_Config)      :: config_lvt
  integer                :: rc
  character*50           :: lsm
  integer                :: nstns
  integer                :: gformat
  character*100, allocatable :: varname(:)
  character*10,  allocatable :: vunits(:)
  character*100, allocatable :: stnname(:)
  integer                :: i,k
  integer                :: wstats
  integer                :: nvars, j
  integer                :: v_index1, v_index2
  character*3            :: fid1, fid2
  real , allocatable         :: rmse(:,:)
  real,  allocatable         :: bias(:,:)
  real, allocatable          :: rcorr(:,:)
  logical                :: file_exists
  character*10           :: ctmp1
  real   , allocatable       :: mv11(:),mv12(:),ov11(:),ov12(:)
  real,     allocatable      :: mv1(:),ov1(:)
  integer                :: yr,mo,da,hr,mn
  logical, allocatable       :: noplot(:)
  integer                :: ios

  k=iargc()
  if(k.lt.1) then 
     print*,'Usage: '
     print*,'plotScatter scatter.config'
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
          label="gnuplot output filename:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,wstats, &
          label="Plot summary stats in the image:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,nvars, &
          label="Number of variables to plot:",rc=rc)

     call ESMF_ConfigGetAttribute(ts_config,d1label, &
          label="data series 1 (model) label:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,d2label, &
          label="data series 2 (obs) label:",rc=rc)

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
     
     call ESMF_ConfigGetAttribute(Config_lvt,lsm,label="Land surface model:",rc=rc)  
     call ESMF_ConfigGetAttribute(config_lvt,expcode,label="Experiment code:",&
          rc=rc)
     call ESMF_ConfigGetAttribute(config_lvt,statsodir,&
          label="Stats output directory:",&
          rc=rc)
     model_name = trim(lsm)
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
        
        inquire(file=rmsefile,exist=file_exists)
        if(file_exists) then 
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
        
        inquire(file=biasfile,exist=file_exists)
        if(file_exists) then 
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
        
        inquire(file=rcorrfile,exist=file_exists)
        if(file_exists) then 
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
     endif

     allocate(noplot(nvars))
     noplot = .false. 
     open(100,file=trim(gfile),form='formatted')
     do i=1,nstns
        mfile = trim(statsodir)//'/MEAN_'//&
             trim(stnname(i))//'.dat'
        print*, 'reading ',trim(mfile)
        ios = 0 
        
199     format(I4, 1x,I2.2, 1x,I2.2,1x, I2.2,1x, I2.2, 1x)
203     format(1x,E14.6,1x,E14.6,1x,E14.6)

        v_index1 = 17
        do j=1,nvars
           if(.not.noplot(j)) then 
              write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"'
              write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"'
              write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"'
              write(100,*) 
              write(100,*) 'begin'
              ofile = trim(stnname(i))//'_'//&
                      trim(varname(j))//'_scatter'
              if(gformat.eq.1) then 
                 write(100,*) 'wks = gsn_open_wks("x11","'//&
                      trim(ofile)//'") '
              elseif(gformat.eq.2) then 
                 write(100,*) 'wks = gsn_open_wks("eps","'//&
                      trim(ofile)//'") '
              endif

              write(100,*) 'res                   = True'
              write(100,*)  'res@gsnFrame        = False'
              write(100,*)  'res@gsnDraw        = False'
              write(100,*)  'res@xyMarkLineMode    = "Markers"'
              write(100,*)  'res@xyMarkers         =  16'
              write(100,*)  'res@xyMarkerColor     = "NavyBlue"'
              write(100,*)  'res@xyMarkerSizeF     = 0.005'
              write(100,*)  'res@tiXAxisString = "'//trim(d1label)//'"'
              write(100,*)  'res@tiYAxisString = "'//trim(d2label)//'"'
              write(100,*)  'res@tiMainString      = "'&
                   //trim(stnname(i))//'-'//&
                   trim(varname(j))//'('//trim(vunits(j))//')"'

              write(100,*)
              write(100,*)  'lnres = True'
              write(100,*)  'lnres@gsLineColor = "Red"'
              write(100,*)
              write(100,*)  'txres               = True'
              write(100,*)  'txres@txFontHeightF = 0.015'
              write(100,*)

              write(100,*) 'fname = "'//trim(statsodir)//&
                   '/MEAN_'//trim(stnname(i))//'.dat"'

              write(100,*) 'data1 = asciiread(fname,-1,"string")'
  
              write(100,*)  'yr = stringtoint(str_get_cols(data1,0,3))'
              write(100,*)  'mo = stringtoint(str_get_cols(data1,5,6))'
              write(100,*)  'da = stringtoint(str_get_cols(data1,8,9))'
              write(100,*)  'hr = stringtoint(str_get_cols(data1,11,12))'
              write(100,*)  'mn = stringtoint(str_get_cols(data1,14,15))'
              write(100,*)  'ss = mn'
  
!              write(100,*) 'set title "'//trim(stnname(i))//'-'//&
!                   trim(varname(j))//'('//trim(vunits(j))//')"'

              v_index2 = v_index1 + 14 -1
              write(unit=fid1,fmt='(i3.3)') v_index1
              write(unit=fid2,fmt='(i3.3)') v_index2

              write(100,*)  'model = stringtofloat(str_get_cols(data1,'//&
                   trim(fid1)//','//trim(fid2)//'))'

              v_index1 = v_index1+ 6*14
              v_index2 = v_index1 + 14 -1
              write(unit=fid1,fmt='(i3.3)') v_index1
              write(unit=fid2,fmt='(i3.3)') v_index2

              write(100,*)  'obs = stringtofloat(str_get_cols(data1,'//&
                   trim(fid1)//','//trim(fid2)//'))'
              write(100,*) 'model@_FillValue=-9999.0'
              write(100,*) 'obs@_FillValue=-9999.0'
              write(100,*)
              write(100,*) 'plot = gsn_csm_xy(wks,model,obs,res)'
              
              write(100,*) 'minx = min(model)'
              write(100,*) 'maxx = max(model)'
              write(100,*)
              write(100,*) 'miny = min(obs)'
              write(100,*) 'maxy = max(obs)'
              write(100,*)
              write(100,*) 'dum = new(1,graphic)'
              write(100,*) 'dum = gsn_add_polyline(wks,plot,(/minx,maxx/), (/minx,maxx/), lnres)'
              if(wstats.eq.1) then
                 write(100,*) 'v1_loc = minx + 2*(maxx-minx)/10'
                 write(100,'(a28)',advance='no') 'gsn_text(wks,plot,"RMSE: "+'
                 write(100,'(F14.3)',advance='no') rmse(i,j)
                 write(100,*) ',v1_loc,maxy,txres)'
                 write(100,*) ''
                 write(100,*) 'v2_loc = maxy - 0.5*(maxy-miny)/10'
                 write(100,'(a28)',advance='no') 'gsn_text(wks,plot,"BIAS: "+'
                 write(100,'(F14.3)',advance='no') bias(i,j)
                 write(100,*) ',v1_loc,v2_loc,txres)'
                 write(100,*) 
                 write(100,*) 'v3_loc = maxy - 1*(maxy-miny)/10'
                 write(100,'(a25)',advance='no') 'gsn_text(wks,plot,"R: "+'
                 write(100,'(F14.3)',advance='no') rcorr(i,j)
                 write(100,*) ',v1_loc,v3_loc,txres)'
                 write(100,*) 
              endif
              write(100,*) 'draw(plot)'
              write(100,*) 'frame(wks)'
              write(100,*) 'end'
              
              write(100,*) 'delete(plot)'
              write(100,*) 'delete(wks)'

           endif
           v_index1 = v_index1+ 6*14
           v_index2 = v_index1 + 14 -1
        enddo
     enddo
  endif
end program plotScatter
