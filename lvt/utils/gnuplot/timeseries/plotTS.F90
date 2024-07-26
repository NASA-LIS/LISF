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
!  This program creates a time series plotting file for gnuplot. 
!  The usage is :
!   plotTS ts.config
! 
!  The ts.config file is structured as follows: 
!   
!   Output format for images:        1
!   LVT config file:                 lvt.config
!   LVT time series locations file:  TS_LOCATIONS.TXT
!   gnuplot output filename:         ts.plt
!   Plot summary stats in the image: 0
!   data series 1 (model) label:     STEP0
!   data series 1 (model) style:     1
!   data series 2 (obs) label:       STEP1
!   data series 2 (obs) style:       1
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
  character*100          :: tsconfigfile
  character*100          :: configfile,gfile,tsfile, statsodir
  character*100          :: rmsefile,biasfile,rcorrfile,mfile
  character*3            :: expcode
  real  , allocatable        :: mmaxv(:),mminv(:),mdy(:),maxv(:),dy(:)
  real  , allocatable        :: omaxv(:),ominv(:),ody(:),minv(:)
  type(ESMF_Config)      :: config_lvt
  integer                :: lsm, rc
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
  logical                :: file_exists
  character*10           :: ctmp1
  real   , allocatable       :: mv11(:),mv12(:),ov11(:),ov12(:)
  real,     allocatable      :: mv1(:),ov1(:)
  integer                :: yr,mo,da,hr,mn,ss
  logical, allocatable       :: noplot(:)
  real*8                 :: maxtime, mintime,time,dx 
  real                   :: gmt
  integer                :: doy
  integer                :: nxtics
  character*20, allocatable  :: xtics(:)
  integer                :: gformat
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
          label="gnuplot output filename:",rc=rc)
     call ESMF_ConfigGetAttribute(ts_config,wstats, &
          label="Plot summary stats in the image:",rc=rc)
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

     open(100,file=trim(gfile),form='formatted')

     allocate(maxv(nvars))
     allocate(minv(nvars))
     allocate(dy(nvars))

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

!     allocate(labelx(nvars))
     allocate(labely(nvars))
     allocate(labeldy(nvars))

     allocate(noplot(nvars))

     open(100,file=trim(gfile),form='formatted')
!     do i=1,1
     do i=1,nstns
        mfile = trim(statsodir)//'/MEAN_'//trim(model_name)//'_'//&
             trim(stnname(i))//trim(expcode)//'.dat'
        print*, 'reading ',trim(mfile)
        open(200,file=trim(mfile),form='formatted')
        mmaxv = -1E20
        mminv = 1E20 
        omaxv = -1E20
        ominv = 1E20

        maxtime = -1E20
        mintime = 1E20

        ios = 0 
        
        do while(ios.eq.0) 
           read(200,199,advance='no',iostat=ios) yr,mo,da,hr,mn
           call date2time(time,doy,gmt,yr,mo,da,hr,mn,ss)

           if(time.gt.maxtime) maxtime = time
           if(time.lt.mintime) mintime = time

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
             
              time = mintime + ((maxtime-mintime)/10)*0.5

              call time2date(time,doy,gmt,yr,mo,da,hr,mn)

              write(labelx,'(i4.4,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2)')  &
                   yr, ' ',mo,' ', da,' ', hr,' ', mn
              maxv(j) = max(mmaxv(j),omaxv(j))
              minv(j) = min(mminv(j),ominv(j))
              dy(j) = ((maxv(j)-minv(j))/10)*0.5

              labely(j) = maxv(j)-((maxv(j)-minv(j))/10)*0.1 
              labeldy(j) = (dy(j))
              noplot(j) = .false. 

           endif
        enddo
        dx = (maxtime-mintime)/nxtics

        do j=1,nxtics
           time = mintime + (j-1)*dx
           call time2date(time,doy,gmt,yr,mo,da,hr,mn)
           
           write(xtics(j),'(i4.4,a1,I2.2,a1,I2.2,a1,I2.2,a3)')  &
                yr, ' ',mo,' ', da,' ', hr,' 00'
           
        enddo


199     format(I4, 1x,I2.2, 1x,I2.2,1x, I2.2,1x, I2.2)
203     format(1x,E14.6,1x,E14.6,1x,E14.6)


        v_index1 = 6
        do j=1,nvars
           if(.not.noplot(j)) then 
              if(gformat.eq.1) then
                 write(100,*) 'set terminal gif large font "Times-Roman" 24'
                 write(100,*) 'set output "'//trim(stnname(i))//'_'//&
                   trim(varname(j))//'_ts.gif"'
              elseif(gformat.eq.2) then 
                 write(100,*) 'set terminal postscript enhanced eps color "Times-Roman" 24'
                 write(100,*) 'set output "'//trim(stnname(i))//'_'//&
                      trim(varname(j))//'_ts.eps"'
                 write(100,*) 'set size 2,1'
              endif
              write(100,*) 'set datafile missing "-0.999900E+04"'
              write(100,*) 'set title "'//trim(stnname(i))//'-'//&
                   trim(varname(j))//'('//trim(vunits(j))//')"'
              write(100,*) 'set xdata time'
              write(100,*) 'set timefmt "%Y %m %d %H %M"'           
              write(100,*) 'set format x "%Y/%m"'
              
              v_index2 = v_index1 + 3 
              if(wstats.eq.1) then 
                 write(100,'(a18)',advance='no') 'set label 1 "RMSE:'
                 write(100,'(F14.3)',advance='no') rmse(i,j)
                 write(100,'(a6)',advance='no') '" at "'
                 write(100,'(a16)',advance='no') trim(labelx)
                 write(100,'(a2)',advance='no') '",'
                 write(100,'(F14.3)',advance='no') labely(j)
                 write(100,*) ' front'
                 
                 write(100,'(a18)',advance='no') 'set label 2 "BIAS:'
                 write(100,'(F14.3)',advance='no') bias(i,j)
                 write(100,'(a6)',advance='no') '" at "'
                 write(100,'(a16)',advance='no') trim(labelx)
                 write(100,'(a2)',advance='no') '",'
                 write(100,'(F14.3)',advance='no') labely(j)-labeldy(j)
                 write(100,*) ' front'
                 
                 write(100,'(a18)',advance='no') 'set label 3 "R:'
                 write(100,'(F14.3)',advance='no') rcorr(i,j)
                 write(100,'(a6)',advance='no') '" at "'
                 write(100,'(a16)',advance='no') trim(labelx)
                 write(100,'(a2)',advance='no') '",'
                 write(100,'(F14.3)',advance='no') labely(j)-2*labeldy(j)
                 write(100,*) ' front'
              endif
              write(unit=fid1,fmt='(i3.3)') v_index1
              write(unit=fid2,fmt='(i3.3)') v_index2
              
              write(100,'(a)',advance='no') 'set xtics ('
              do kk=1,nxtics
                 write(100,'(a)',advance='no') '"'//trim(xtics(kk))//'"'
                 if(kk.ne.nxtics) then 
                    write(100,'(a)',advance='no') ','
                 else
                    write(100,'(a)') ')'
                 endif
              enddo

              write(100,'(a)',advance='no') trim(adjustl('plot "'//trim(statsodir)//'/MEAN_'//trim(model_name)//'_'//&
                   trim(stnname(i))//trim(expcode)//'.dat" using '//&
                   '1:'//trim(fid1)//' title "'//trim(d1label)))//'" '
              if(d1style.eq.1) then 
                 write(100,'(a)') trim(adjustl(" with lines lt 3 lw 3,\"))
              elseif(d1style.eq.2) then 
                 write(100,'(a)') '" with points,\'
              endif

              write(100,'(a)',advance='no') '"'//trim(statsodir)//'/MEAN_'//trim(model_name)//'_'//&
                   trim(stnname(i))//trim(expcode)//'.dat" using '//&
                   '1:'//trim(fid2)//' title "'//trim(d2label)
              if(d2style.eq.1) then 
                 write(100,'(a)') '" with lines lt -1 lw 3'
              elseif(d2style.eq.2) then 
                 write(100,'(a)') '" with points 7'
              endif

              write(100,*) 'clear'
              write(100,*) ''
              v_index1 = v_index2+3
           endif
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
