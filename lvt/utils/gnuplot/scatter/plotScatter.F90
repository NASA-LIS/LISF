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
!  This program creates a time series plotting file for gnuplot. 
!  The usage is :
!   plotScatter scatter.config
! 
!  The scatter.config file is structured as follows: 
!   
!   Output format for images:        1
!   LVT config file:                 lvt.config
!   LVT time series locations file:  TS_LOCATIONS.TXT
!   gnuplot output filename:         ts.plt
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
  character*50           :: model_name
  character*50           :: d1label, d2label
  character*100          :: configfile,gfile,tsfile, statsodir
  character*100          :: rmsefile,biasfile,rcorrfile,mfile
  character*3            :: expcode
  real  , allocatable        :: mmaxv(:),mminv(:),mdy(:)
  real  , allocatable        :: omaxv(:),ominv(:),ody(:)
  type(ESMF_Config)      :: config_lvt
  integer                :: lsm, rc
  integer                :: nstns
  integer                :: gformat
  character*100, allocatable :: varname(:)
  character*10,  allocatable :: vunits(:)
  character*100, allocatable :: stnname(:)
  real    , allocatable      :: labelx(:)
  real   , allocatable       :: labely(:)
  real    , allocatable      :: labeldy(:)
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

     if(wstats.eq.1) then
        rmsefile = trim(statsodir)//'/RMSE_'//trim(model_name)//'_'//&
             trim(expcode)//'_SUMMARY_STATS.dat'
        
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
        
        biasfile = trim(statsodir)//'/BIAS_'//trim(model_name)//'_'//&
             trim(expcode)//'_SUMMARY_STATS.dat'
        
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
        
        rcorrfile = trim(statsodir)//'/RCORR_'//trim(model_name)//'_'//&
             trim(expcode)//'_SUMMARY_STATS.dat'
        
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

     allocate(mmaxv(nvars))
     allocate(mminv(nvars))
     allocate(mdy(nvars))
     allocate(mv1(nvars))
     allocate(mv11(nvars))
     allocate(mv12(nvars))

     allocate(omaxv(nvars))
     allocate(ominv(nvars))
     allocate(ody(nvars))
     allocate(ov1(nvars))
     allocate(ov11(nvars))
     allocate(ov12(nvars))

     allocate(labelx(nvars))
     allocate(labely(nvars))
     allocate(labeldy(nvars))

     allocate(noplot(nvars))

     open(100,file=trim(gfile),form='formatted')
     do i=1,nstns
        mfile = trim(statsodir)//'/MEAN_'//trim(model_name)//'_'//&
             trim(stnname(i))//trim(expcode)//'.dat'
        print*, 'reading ',trim(mfile)
        open(200,file=trim(mfile),form='formatted')
        mmaxv = -1E20
        mminv = 1E20 
        omaxv = -1E20
        ominv = 1E20
        
        ios = 0 
        
        do while(ios.eq.0) 
           read(200,199,advance='no',iostat=ios) yr,mo,da,hr,mn
           read(200,*,iostat=ios) (mv1(j),mv11(j),mv12(j),&
                ov1(j),ov11(j),ov12(j), j=1,nvars)
            if(ios.ne.0) exit

           do j=1,nvars
              if(mv1(j).ne.-9999.0.and.mv1(j).gt.mmaxv(j)) mmaxv(j) = mv1(j)
              if(mv1(j).ne.-9999.0.and.mv1(j).lt.mminv(j)) mminv(j) = mv1(j)
              if(ov1(j).ne.-9999.0.and.ov1(j).gt.omaxv(j)) omaxv(j) = ov1(j)
              if(ov1(j).ne.-9999.0.and.ov1(j).lt.ominv(j)) ominv(j) = ov1(j)
           enddo
        enddo
        close(200)
        noplot = .true. 
        do j=1,nvars
           if(mmaxv(j).ne.-1E20.and.omaxv(j).ne.-1E20) then !nothing to plot
              mdy(j) = ((mmaxv(j)-mminv(j))/10)*0.5
              ody(j) = ((omaxv(j)-ominv(j))/10)*0.5
              
              labelx(j) = (mminv(j)+((mmaxv(j)-mminv(j))/10)*0.1 )
              labely(j) = (omaxv(j)-((omaxv(j)-ominv(j))/10)*0.1 )
              labeldy(j) = (ody(j))
              noplot(j) = .false. 
              
           endif
        enddo

199     format(I4, 1x,I2.2, 1x,I2.2,1x, I2.2,1x, I2.2, 1x)
203     format(1x,E14.6,1x,E14.6,1x,E14.6)

        v_index1 = 6
        do j=1,nvars
           if(.not.noplot(j)) then 
              if(gformat.eq.1) then 
                 write(100,*) 'set terminal gif large font "Times-Roman" 24'
                 write(100,*) 'set output "'//trim(stnname(i))//'_'//&
                      trim(varname(j))//'_scatter.gif"'
              elseif(gformat.eq.2) then 
                 write(100,*) 'set terminal postscript enhanced eps color "Times-Roman" 24'
                 write(100,*) 'set output "'//trim(stnname(i))//'_'//&
                      trim(varname(j))//'_scatter.eps"'
              endif
              write(100,*) 'set datafile missing "-0.999900E+04"'
              write(100,*) 'set ylabel "'//trim(d2label)//'"'
              !        write(100,*) 'set xlabel "'//trim(model_name)//'"'
              write(100,*) 'set xlabel "'//trim(d1label)//'"'

              write(100,*) 'set title "'//trim(stnname(i))//'-'//&
                   trim(varname(j))//'('//trim(vunits(j))//')"'

              v_index2 = v_index1 + 3 
              if(wstats.eq.1) then
                 write(100,'(a18)',advance='no') 'set label 1 "RMSE:'
                 write(100,'(F14.3)',advance='no') rmse(i,j)
                 write(100,'(a5)',advance='no') '" at '
                 write(100,'(F14.3)',advance='no') labelx(j)
                 write(100,'(a1)',advance='no') ','
                 write(100,'(F14.3)',advance='no') labely(j)
                 write(100,*) ' front'

                 write(100,'(a18)',advance='no') 'set label 2 "BIAS:'
                 write(100,'(F14.3)',advance='no') bias(i,j)
                 write(100,'(a5)',advance='no') '" at '
                 write(100,'(F14.3)',advance='no') labelx(j)
                 write(100,'(a1)',advance='no') ','
                 write(100,'(F14.3)',advance='no') labely(j)-labeldy(j)
                 write(100,*) ' front'

                 write(100,'(a18)',advance='no') 'set label 3 "R:'
                 write(100,'(F14.3)',advance='no') rcorr(i,j)
                 write(100,'(a5)',advance='no') '" at '
                 write(100,'(F14.3)',advance='no') labelx(j)
                 write(100,'(a1)',advance='no') ','
                 write(100,'(F14.3)',advance='no') labely(j)-2*labeldy(j)
                 write(100,*) ' front'
              endif
              write(unit=fid1,fmt='(i3.3)') v_index1
              write(unit=fid2,fmt='(i3.3)') v_index2
              
              write(100,'(a100)') 'plot "'//trim(statsodir)//'/MEAN_'//trim(model_name)//'_'//&
                trim(stnname(i))//trim(expcode)//'.dat" using '//&
                trim(fid1)//':'//trim(fid2)//' notitle with points 3,\'
              write(100,*) 'x notitle with lines'
              write(100,*) 'clear'
              write(100,*) ''
              v_index1 = v_index2+3
           endif
        enddo
     enddo
  endif
end program plotScatter
