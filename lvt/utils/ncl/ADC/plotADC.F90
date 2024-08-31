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
! !ROUTINE: plotADC
!  \label{plotADC}
!
program plotADC
! !USES: 
  use ESMF
!
! !DESCRIPTION: 
!  This program creates an average diurnal curve plotting file for NCL. 
!  The usage is :
!   plotADC ts.config
! 
!  The ts.config file is structured as follows: 
!   
!   Output format for images:        1
!   LVT config file:                 lvt.config
!   LVT time series locations file:  TS_LOCATIONS.TXT
!   gnuplot output filename:         ts.plt
!   data series 1 (model) label:     STEP0
!   data series 1 (model) style:     1
!   data series 2 (obs) label:       STEP1
!   data series 2 (obs) style:       1
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
  character*50           :: d1label, d2label
  integer                :: d1style, d2style
  character*100          :: tsconfigfile,ofile
  character*50           :: model_name
  character*100          :: configfile,gfile,tsfile, statsodir
  character*3            :: expcode
  type(ESMF_Config)      :: config_lvt
  integer                :: rc
  character*50           :: lsm
  integer                :: gformat
  integer                :: nstns
  character*100, allocatable :: varname(:)
  character*10,  allocatable :: vunits(:)
  character*100, allocatable :: stnname(:)
  integer                :: i,k
  integer                :: nvars, j
  integer                :: v_index1, v_index2
  character*3            :: fid1, fid2


  k=iargc()
  if(k.lt.1) then
     print*,'Usage: '
     print*,'plotTS adc.config'
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
     call ESMF_ConfigGetAttribute(ts_config,nvars, &
          label="Number of variables to plot:",rc=rc)

     call ESMF_ConfigGetAttribute(ts_config,d1label, &
          label="data series 1 (model) label:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,d2label, &
          label="data series 2 (obs) label:",rc=rc)

     call ESMF_ConfigGetAttribute(ts_config,d1style, &
          label="data series 1 (model) style:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,d2style, &
          label="data series 2 (obs) style:",rc=rc)

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

     open(100,file=trim(gfile),form='formatted')
     do i=1,nstns
        do j=1,nvars
           v_index1 = 1
           write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"'
           write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"'
           write(100,*) 'load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"'
           write(100,*) 
           write(100,*) 'begin'
           
           ofile = trim(stnname(i))//'_'//trim(varname(j))//'_adc'
           if(gformat.eq.1) then 
              write(100,*) 'wks = gsn_open_wks("x11","'//&
                   trim(ofile)//'") '
           elseif(gformat.eq.2) then 
              write(100,*) 'wks = gsn_open_wks("eps","'//&
                   trim(ofile)//'") '
           endif
           write(100,*) 'fname ="'//trim(statsodir)//'/MEAN_ADC_'//&
                trim(stnname(i))//'_'//trim(varname(j))//'.dat"'
           write(100,*) 'data1 = asciiread(fname,-1,"string")'


           v_index2 = v_index1 + 3
           write(unit=fid1,fmt='(i3.3)') v_index1-1
           write(unit=fid2,fmt='(i3.3)') v_index2
           write(100,*) 'time = stringtofloat(str_get_cols(data1,'//&
                trim(fid1)//','//trim(fid2)//'))'

           v_index1 = v_index2 + 1
           v_index2 = v_index1 + 22
           write(unit=fid1,fmt='(i3.3)') v_index1
           write(unit=fid2,fmt='(i3.3)') v_index2
           write(100,*) 'model = stringtofloat(str_get_cols(data1,'//&
                trim(fid1)//','//trim(fid2)//'))'

           v_index1 = v_index2 + 6
           v_index2 = v_index1 + 22
           write(unit=fid1,fmt='(i3.3)') v_index1
           write(unit=fid2,fmt='(i3.3)') v_index2
           write(100,*) 'obs = stringtofloat(str_get_cols(data1,'//&
                trim(fid1)//','//trim(fid2)//'))'

           write(100,*) 'res = True'
           write(100,*)  'res@gsnFrame        = False'
           write(100,*)  'res@gsnDraw        = False'
!           write(100,*) 'res@vpWidthF  = 0.8'
!           write(100,*) 'res@vpHeightF = 0.2'
           write(100,*) 'res@xyLineThicknesses = (/2.0,2.0/)'               
           write(100,*) 'res@xyLineColors      = (/"blue","red"/)'
           write(100,*) 'res@xyMarkLineMode = "MarkLines"'
           write(100,*) 'res@xyMarkers      = (/6/)'
           write(100,*) 'res@xyMarkerColors = (/"blue"/)'
           
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

           write(100,*) 'res@tiXAxisString = "Hour"'
           write(100,*) 'var = new((/2,dimsizes(time)/),float)'
           write(100,*) 'var(0,:) = model'
           write(100,*) 'var(1,:) = obs'
           
           write(100,*) 'var@_FillValue = -9999.0'
           write(100,*) 'plot = gsn_csm_xy(wks,time,var,res)'
           
           write(100,*) 'draw(plot)'
           write(100,*) 'frame(wks)'
           write(100,*) 'end'
           write(100,*) 
           write(100,*) 'delete(plot)'
           write(100,*) 'delete(wks)'
           
!           v_index1 = v_index2+3
        enddo
     enddo
  endif
end program plotADC
