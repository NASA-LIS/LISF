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
!  This program creates an average diurnal curve plotting file for gnuplot. 
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
  character*100          :: tsconfigfile  
  character*50           :: model_name
  character*100          :: configfile,gfile,tsfile, statsodir
  character*3            :: expcode
  type(ESMF_Config)      :: config_lvt
  integer                :: lsm, rc
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
          label="gnuplot output filename:",rc=rc)
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
     call ESMF_ConfigGetAttribute(config_lvt,expcode,label="Experiment code:",&
          rc=rc)
     call ESMF_ConfigGetAttribute(config_lvt,statsodir,&
          label="Stats output directory:",&
          rc=rc)
     if(lsm.eq.1) then 
        model_name = 'NOAH271'
     elseif(lsm.eq.2) then 
        model_name = 'CLM'
     elseif(lsm.eq.4) then 
        model_name = 'MOS'
     elseif(lsm.eq.5) then 
        model_name = 'HYSSIB'
     elseif(lsm.eq.7) then 
        model_name = 'CLSM'
     elseif(lsm.eq.14) then 
        model_name = 'TESSEL'
     elseif(lsm.eq.21) then 
        model_name = 'NOAH31'
     elseif(lsm.eq.22) then 
        model_name = 'NOAH32'
     endif
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
        v_index1 = 1
        do j=1,nvars
           if(gformat.eq.1) then 
              write(100,*) 'set terminal gif large font "Times-Roman" 24'
              write(100,*) 'set output "'//trim(stnname(i))//'_'//&
                   trim(varname(j))//'_adc.gif"'
           elseif(gformat.eq.2) then 
              write(100,*) 'set terminal postscript enhanced eps color "Times-Roman" 24'
              write(100,*) 'set output "'//trim(stnname(i))//'_'//&
                   trim(varname(j))//'_adc.eps"'
           endif

           write(100,*) 'set datafile missing "-0.999900E+04"'
           write(100,*) 'set key top left'
           write(100,*) 'set ylabel "'//trim(varname(j))//'"'
           write(100,*) 'set xlabel "Hour"'
           write(100,*) 'set title "'//trim(stnname(i))//'-'//&
                trim(varname(j))//'('//trim(vunits(j))//')"'

           v_index2 = v_index1 + 1

           write(unit=fid1,fmt='(i3.3)') v_index1
           write(unit=fid2,fmt='(i3.3)') v_index2

           write(100,'(a)',advance='no') 'plot "'//trim(statsodir)//'/MEAN_ADC_'//&
                trim(stnname(i))//'_'//trim(varname(j))//'.dat" using '// & 
                trim(fid1)//':'//trim(fid2)//' title "'//trim(d1label)//'" '
           if(d1style.eq.1) then 
              write(100,'(a)') trim(adjustl(" with linespoints lt 3 lw 3,\"))
           elseif(d1style.eq.2) then 
              write(100,'(a)') '" with points,\'
           endif
           

           v_index2 = v_index1 + 3
           write(unit=fid2,fmt='(i2.2)') v_index2
           write(100,'(a)',advance='no') ' "'//trim(statsodir)//'/MEAN_ADC_'//&
                trim(stnname(i))//'_'//trim(varname(j))//'.dat" using '// & 
                trim(fid1)//':'//trim(fid2)//' title "'//trim(d2label)
           if(d2style.eq.1) then 
              write(100,'(a)') '" with linespoints lt -1 lw 3'
           elseif(d2style.eq.2) then 
              write(100,'(a)') '" with points 7'
           endif

           write(100,*) 'clear'
           write(100,*) ''
!           v_index1 = v_index2+3
        enddo
     enddo
  endif
end program plotADC
