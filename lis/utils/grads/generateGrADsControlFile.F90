!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! !ROUTINE: generateGrADsControlFile
! \label{generateGrADsControlFile}
! 
! !DESCRIPTION: 
!  This program generates a sample GrADs control file by parsing the 
!  Statistics file written out by LIS.
!  
! !REMARKS:
!  The program generates only a template. The actual descriptions of 
!  the variables need to be filled in by the user. It only works for 
!  lat/lon domains and binary data.
!
!  For generating NETCDF control files, first a file containing 
!  variable names must be generated. 
!   e.g. ncdump -h 200210290300.d01.nc | grep float  > list.dat
!   (Then parse the file to reduce it to just the list of variables).  
!  
! !REVISION HISTORY:
! 07 May 2007 : Sujay Kumar , Initial Specification
! 11 Feg 2012 : Sujay Kumar , Support for netcdf files are added. 
! 
! !ROUINE: 
program generateGrADsControlFile

  use ESMF
#if ( defined AIX )
  use xlfutility, only : iargc
#endif

  implicit none

  integer              :: nsm_levs, nst_levs 
  character*80         :: gradsfile,lisdir,modelname,configfile, listfile
  character*80         :: statsfile 
  character*40         :: lsm
  character*20         :: wform
  character*100, allocatable :: varname(:)
  character*100        :: dummy
  character*3          :: expcode 
  character*100        :: ctmp
  character*1          :: nestid(3)
  character*100        :: dset
  character*1          :: dset1(100)
  integer              :: rc
  real                 :: curr_time
  integer              :: nnest
  type(ESMF_Config)    :: config_lis
  integer              :: i,kk,c,n,ios,count
  integer              :: nc,nr,nvars
  integer              :: syr, smo, sda, shr, smn, sss
  integer              :: eyr, emo, eda, ehr, emn, ess
  real, allocatable    :: run_dd(:,:)
  integer, allocatable :: nts(:)
  integer, allocatable :: outInterval(:)
  character*3          :: cmo
  integer              :: nsm, nst
  logical              :: sm_found, st_found
  real                 :: mean,std,min,max
  integer              :: varind1,varind2
  character*10         :: temp
  character*40         :: stime
  character*19         :: var(100),var_temp, var2(100),var3(100)
  real*8  :: time
  integer :: doy
  real    :: gmt

  i=iargc()
  if(i.lt.2) then 
     print*,'Usage: '
     print*,'generateGrADsControlFile <lisconfig file> <gradsfile>'
     print*, 'for NETCDF files: '
     print*,'generateGrADsControlFile <lisconfig file> <gradsfile> <varlist> '
  else
     call getarg(1,configfile)
     call getarg(2,gradsfile)
     
     if(i.eq.3) then 
        call getarg(3,listfile)
     endif

     sm_found = .false.
     st_found = .false. 

     config_lis = ESMF_ConfigCreate(rc=rc)
     call ESMF_ConfigLoadFile(config_lis,trim(configfile),rc=rc)
     
     call ESMF_ConfigGetAttribute(config_lis,nnest,label="Number of nests:",&
          default=1)
     call ESMF_ConfigGetAttribute(config_lis,lsm,label="Land surface model:")
     call ESMF_ConfigGetAttribute(config_lis,wform,label="Output data format:")
     
     call ESMF_ConfigGetAttribute(config_lis,lisdir,label="Output directory:",&
          default='OUTPUT')
     
     if(nnest.gt.1) then 
        print*, 'Warning: This program will only create a control file '
        print*, 'for a single nest....'
     endif
     allocate(run_dd(nnest,6))
     allocate(nts(nnest))
     
     call ESMF_ConfigFindLabel(config_lis,"Run domain lower left lat:",rc=rc)
     do i=1,nnest
        call ESMF_ConfigGetAttribute(config_lis,run_dd(i,1),rc=rc)
     enddo
     call ESMF_ConfigFindLabel(config_lis,"Run domain lower left lon:",rc=rc)
     do i=1,nnest
        call ESMF_ConfigGetAttribute(config_lis,run_dd(i,2),rc=rc)
     enddo
     
     call ESMF_ConfigFindLabel(config_lis,"Run domain upper right lat:",rc=rc)
     do i=1,nnest
        call ESMF_ConfigGetAttribute(config_lis,run_dd(i,3),rc=rc)
     enddo
     
     call ESMF_ConfigFindLabel(config_lis,"Run domain upper right lon:",rc=rc)
     do i=1,nnest
        call ESMF_ConfigGetAttribute(config_lis,run_dd(i,4),rc=rc)
     enddo
     
     call ESMF_ConfigFindLabel(config_lis,"Run domain resolution (dx):",rc=rc)
     do i=1,nnest
        call ESMF_ConfigGetAttribute(config_lis,run_dd(i,5),rc=rc)
     enddo
     
     call ESMF_ConfigFindLabel(config_lis,"Run domain resolution (dy):",rc=rc)
     do i=1,nnest
        call ESMF_ConfigGetAttribute(config_lis,run_dd(i,6),rc=rc)
     enddo
     
     call ESMF_ConfigGetAttribute(config_lis,syr,label="Starting year:",&
          default=2000)
     call ESMF_ConfigGetAttribute(config_lis,smo,label="Starting month:",&
          default=1)
     call ESMF_ConfigGetAttribute(config_lis,sda,label="Starting day:",&
       default=1)
     call ESMF_ConfigGetAttribute(config_lis,shr,label="Starting hour:",&
          default=12)
     call ESMF_ConfigGetAttribute(config_lis,smn,label="Starting minute:",&
          default=0)
     call ESMF_ConfigGetAttribute(config_lis,sss,label="Starting second:",&
          default=0)
     call ESMF_ConfigGetAttribute(config_lis,eyr,label="Ending year:",&
          default=2000)
     call ESMF_ConfigGetAttribute(config_lis,emo,label="Ending month:",&
          default=1)
     call ESMF_ConfigGetAttribute(config_lis,eda,label="Ending day:",&
          default=1)
     call ESMF_ConfigGetAttribute(config_lis,ehr,label="Ending hour:",&
          default=12)
     call ESMF_ConfigGetAttribute(config_lis,emn,label="Ending minute:",&
          default=0)
     call ESMF_ConfigGetAttribute(config_lis,ess,label="Ending second:",&
          default=0)
     
     if(smo.eq.1) then 
        cmo = 'jan'
     elseif(smo.eq.2) then
        cmo = 'feb'
     elseif(smo.eq.3) then
        cmo = 'mar'
     elseif(smo.eq.4) then
        cmo = 'apr'
     elseif(smo.eq.5) then
        cmo = 'may'
     elseif(smo.eq.6) then
        cmo = 'jun'
     elseif(smo.eq.7) then
        cmo = 'jul'
     elseif(smo.eq.8) then
        cmo = 'aug'
     elseif(smo.eq.9) then
        cmo = 'sep'
     elseif(smo.eq.10) then
        cmo = 'oct'
     elseif(smo.eq.11) then
        cmo = 'nov'
     elseif(smo.eq.12) then
        cmo = 'dec'
     endif
     
     allocate(outInterval(nnest))

     if(lsm.eq."Template") then 
        modelname = 'TEMPLATE'
     elseif(lsm.eq."Noah.2.7.1") then 
        modelname = 'NOAH271'
        nsm_levs = 4
        nst_levs = 4
     elseif(lsm.eq."Noah.3.3") then 
        call ESMF_ConfigFindLabel(config_lis,&
             label="Noah.3.3 model timestep:",rc=rc)
        do i=1,nnest
           call ESMF_ConfigGetAttribute(config_lis, stime,rc=rc)
           call parseTimestring(stime, nts(i))
        enddo
        
        modelname = 'NOAH33'
        nsm_levs = 4
        nst_levs = 4
     elseif(lsm.eq."CLM2") then 
        modelname = 'CLM'
        nsm_levs = 10
        nst_levs = 10
     elseif(lsm.eq."Mosaic") then 
        modelname = 'MOS'
        nsm_levs = 3
        nst_levs = 3
     elseif(lsm.eq."HySSIB") then 
        modelname = 'HYSSIB'
     elseif(lsm.eq."CLSM") then 
        modelname = 'CLSM'
        nsm_levs = 3
        nst_levs = 5
     endif
     
     call ESMF_ConfigFindLabel(config_lis,"Surface model output interval:",&
          rc=rc)

     do i=1,nnest
        call ESMF_ConfigGetAttribute(config_lis,stime,rc=rc)
        call parseTimestring(stime, outInterval(i))
     enddo
     
     curr_time = float(shr)*3600+60*float(smn)+float(sss)

     do while(mod(curr_time,real(outInterval(1))).ne.0)  
        call tick(time,doy,gmt,syr,smo,sda,shr,smn,sss,nts(1))
        curr_time = float(shr)*3600+60*float(smn)+float(sss)
     enddo
     
     if(wform.eq."binary") then 
        do n=1,1
           write(unit=ctmp,fmt='(a1,i2.2)') 'd ',n
           read(unit=ctmp,fmt='(3a1)') nestid
           
           dset = ' DSET ^'//trim(lisdir)//'/SURFACEMODEL/'//&
                '/%y4/%y4%m2%d2/%y4%m2%d2%h2%n2.'//nestid(1)//&
                nestid(2)//nestid(3)//'.gs4r'
           
           statsfile = trim(lisdir)//'/SURFACEMODEL.d01.stats'
           
           count = 1
           open(111,file=trim(statsfile),form='formatted')
           read(111,*)
           read(111,*) 
           read(111,*)
           read(111,*)
           
           ios = 0 
           do while(ios.eq.0) 
              read(111,fmt='(a19,4F14.3)',iostat=ios) var(count),mean,std,min,max
              if(var(count).eq.'') ios = -9
              varind1 = index(var(count),'(')
              varind2 = index(var(count),')')
              var2(count) = var(count)(1:varind1-1)
              var3(count) = var(count)(varind1+1:varind2-1)
              count = count+1
           enddo
           
           
           write(unit=ctmp,fmt='(a100)') dset
           read(unit=ctmp, fmt='(100a1)') (dset1(i),i=1,100)
           
           open(100,file=trim(gradsfile),form='formatted')
           write(100,fmt='(100A1)') (dset1(i),i=1,100)
           
           write(100,*) 'OPTIONS template'
           write(100,*) 'OPTIONS sequential'
           write(100,*) 'OPTIONS big_endian'
           write(100,*) 'TITLE LIS output'
           write(100,*) 'UNDEF -9999.0'
           nc = (nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,6))) + 1
           nr = (nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,5))) + 1
           
           write(100,222) ' XDEF ',nc,' LINEAR ',run_dd(n,2),run_dd(n,5)
           write(100,222) ' YDEF ',nr,' LINEAR ',run_dd(n,1),run_dd(n,6)
           write(100,*) 'ZDEF 1 LINEAR 1 1'
           if(outInterval(n).lt.3600) then 
              outInterval(n) = outInterval(n)/60
              write(100,223) 'TDEF 1 LINEAR ',shr,':',smn,'Z',&
                   sda,cmo,syr,' ',outInterval(n),'mn'
           else 
              outInterval(n) = outInterval(n)/3600
              write(100,223) 'TDEF 1 LINEAR ',shr,':',smn,'Z',&
                   sda,cmo,syr,' ',outInterval(n),'hr'
           endif
        
           write(100,fmt='(a7,I2)') ' VARS ',count-2
           do i=1,count-2
              write(100,fmt='(a19,a11,i2,a25,a19)') var2(i),'1  99  ** ',i,&
                   ' insert description here ',var3(i)
           enddo

           write(100,*) 'ENDVARS'
           close(100)
           print*, 'Generated the control file ',trim(gradsfile),' successfully'
        end do
     elseif(wform.eq."netcdf") then !netcdf file
        
        open(20,file=trim(listfile))
        ios = 0
        nvars = 0 
        do while(ios.eq.0) 
           read(20,'(a)',iostat=ios) dummy
           nvars = nvars + 1
        enddo
        close(20)
        nvars = nvars-1

        allocate(varname(nvars))
        open(20,file=trim(listfile))
        do kk=1,nvars
           read(20,'(a)',iostat=ios) varname(kk)
        enddo
        close(20)
                
        do n=1,1
           write(unit=ctmp,fmt='(a1,i2.2)') 'd ',n
           read(unit=ctmp,fmt='(3a1)') nestid
           
           dset = ' DSET ^'//trim(lisdir)//'/SURFACEMODEL/'//&
                '/LIS_HIST_%y4%m2%d2%h2%n2.'//nestid(1)//&
                nestid(2)//nestid(3)//'.nc'
           
           statsfile = trim(lisdir)//'/SURFACEMODEL.d01.stats'
           
           count = 1
           open(111,file=trim(statsfile),form='formatted')
           read(111,*)
           read(111,*) 
           read(111,*)
           read(111,*)
           
           ios = 0 
           nsm = 0 
           nst = 0 
           var2 = "a"
           do while(ios.eq.0) 
              read(111,fmt='(a19,4F14.3)',iostat=ios) var_temp,mean,std,min,max
              if(var_temp.eq.'') ios = -9
              varind1 = index(var_temp,'(')
              varind2 = index(var_temp,')')
              var2(count) = var_temp(1:varind1-1)
              
              var3(count) = var_temp(varind1+1:varind2-1)
              
              if(trim(adjustl(var2(count))).eq."SoilMoist") then 
                 nsm = nsm + 1
              endif

              if(trim(adjustl(var2(count))).eq."SoilTemp") then 
                 nst = nst + 1
              endif
              
              if(nsm.eq.nsm_levs) then 
                 count = count + 1
                 nsm = -1
              endif
              if(nst.eq.nst_levs) then 
                 count = count + 1
                 nst = -1
              endif

              if(trim(adjustl(var2(count))).ne."SoilMoist".and.&
                   trim(adjustl(var2(count))).ne."SoilTemp".and.&
                   trim(adjustl(var2(count))).ne."a") then 
                 count = count+1
              endif
           enddo

           sm_found = .false. 
           st_found = .false. 

           write(unit=ctmp,fmt='(a100)') dset
           read(unit=ctmp, fmt='(100a1)') (dset1(i),i=1,100)
           
           open(100,file=trim(gradsfile),form='formatted')
           write(100,fmt='(100A1)') (dset1(i),i=1,100)
           write(100,*) 'DTYPE netcdf'
           write(100,*) 'OPTIONS template'
           write(100,*) 'TITLE LIS output'
           nc = (nint((run_dd(n,4)-run_dd(n,2))/run_dd(n,6))) + 1
           nr = (nint((run_dd(n,3)-run_dd(n,1))/run_dd(n,5))) + 1
           
           write(100,224) ' XDEF east_west   ',nc,' LINEAR ',run_dd(n,2),run_dd(n,5)
           write(100,224) ' YDEF north_south ',nr,' LINEAR ',run_dd(n,1),run_dd(n,6)
           write(100,*) 'ZDEF SoilMoist_profiles 4 LEVELS 0.1 0.3 0.6 1.0'
224  format (a18,I7,A9,F9.3,F7.3)
           if(outInterval(n).lt.3600) then 
              outInterval(n) = outInterval(n)/60
              write(100,225) ' TDEF  time 1 LINEAR ',shr,':',smn,'Z',&
                   sda,cmo,syr,' ',outInterval(n),'mn'
           else 
              outInterval(n) = outInterval(n)/3600
              write(100,225) ' TDEF time 1 LINEAR ',shr,':',smn,'Z',&
                   sda,cmo,syr,' ',outInterval(n),'hr'
           endif
225 format(a21,I2.2,a1,I2.2,a1,I2.2,a3,I4,a1,i2.2,a2)           
           write(100,fmt='(a7,I2)') ' VARS ',nvars
           do i=1,nvars
              if(trim(adjustl(var2(i))).eq."SoilMoist") then 
                 if(.not.sm_found) then 
                    write(temp,fmt='(i2)') nsm_levs
                    write(100,fmt='(a)') trim(adjustl(varname(i)))//'=>'//&
                         trim(adjustl(var2(i)))//'    '//trim(adjustl(temp))//&
                         '      z,y,x'   
                    sm_found = .true. 
                 endif
              elseif(trim(adjustl(var2(i))).eq."SoilTemp") then 
                 if(.not.st_found) then 
                    write(temp,fmt='(i2)') nst_levs
                    write(100,fmt='(a)') trim(adjustl(varname(i)))//'=>'//&
                         trim(adjustl(var2(i)))//'    '//trim(adjustl(temp))//&
                         '      z,y,x'             
                    st_found = .true. 
                 endif
              else
                 write(100,fmt='(a)') trim(adjustl(varname(i)))//'=>'//&
                      trim(adjustl(var2(i)))//'    0      y,x'
              endif
           enddo

           write(100,*) 'ENDVARS'
           close(100)
           print*, 'Generated the control file ',trim(gradsfile),' successfully'
        end do
     endif
  end if
222 format (a6,I7,A9,F9.3,F7.3)
223 format(a15,I2.2,a1,I2.2,a1,I2.2,a3,I4,a1,i2.2,a2)


  
end program generateGrADsControlFile

!BOP
! 
! !ROUTINE: tick
! \label{tick}
! 
! !INTERFACE:
  subroutine tick(time,doy,gmt,yr,mo,da,hr,mn,ss,ts)
    implicit none
! !ARGUMENTS:
    real*8  :: time
    integer :: yr,mo,da,hr,mn,ss,ts,doy
    real    :: gmt
! !DESCRIPTION:
! 
!  Method to advance or retract the time by the specified time
!  increment
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
!  \item[ts]
!   time increment in seconds
!  \item[time]
!   lis time
!  \item[gmt]
!   time in GMT
!  \item[doy]
!   day of the year
!  \end{description}
!EOP
    integer days(13)
    integer prvmo   !previous month
    
    data days/31,28,31,30,31,30,31,31,30,31,30,31,31/

143 format(a1,' yr',i6,' mo',i5,' dy',i5,' hr',i5, & 
         ' mn',i6,' ss',i8,' ts',i8)
    ss=ss+ts
    do while(ss.gt.59)
       ss=ss-60
       mn=mn+1
    enddo
    do while(ss.lt.0)
       ss=ss+60
       mn=mn-1
    enddo
    do while(mn.gt.59)
       mn=mn-60
       hr=hr+1
    enddo
    
    do while(mn.lt.0)
       mn=mn+60
       hr=hr-1
    enddo
    do while(hr.gt.23)
       hr=hr-24
       da=da+1
    enddo
    
    do while(hr.lt.0)
       hr=hr+24
       da=da-1
    enddo
    
    if((mod(yr,4).eq.0.and.mod(yr,100).ne.0) &    !correct for leap year
         .or.(mod(yr,400).eq.0))then               !correct for y2k
       days(2)=29                  
    else
       days(2)=28
    endif
    
    do while(da.gt.days(mo))
       da=da-days(mo)
       mo=mo+1
    enddo
    
    do while(da.lt.1)
       
       prvmo=mo-1
       if(mo.eq.1) prvmo=12
       
       da=da+days(prvmo)
       
       if(prvmo.eq.12) then
          mo=prvmo
          yr=yr-1
       else
          mo=prvmo
       endif
    enddo
    do while(mo.gt.12)
       mo=mo-12
       yr=yr+1
    enddo

    do while(mo.lt.1)
       mo=mo+12
       yr=yr-1
    enddo
    call date2time(time,doy,gmt,yr,mo,da,hr,mn,ss)
    return

  end subroutine tick
!BOP
!
! !ROUTINE: date2time
! \label{date2time} 
! 
! !INTERFACE:
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
!   lis time
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

  subroutine parseTimeString(inputstr,time_value)
    
    implicit none
    
    character(len=*)          :: inputstr
    integer                   :: time_value
    
    character*2        :: suffix
    real               :: time_input
    integer            :: iloc

    iloc  = len(trim(inputstr))-1
    suffix = inputstr(iloc:len(inputstr))
    
    read(inputstr(1:iloc-1),*) time_input

    if(suffix.eq."ss") then 
       time_value = nint(time_input)
    elseif(suffix.eq."mn") then 
       time_value = nint(time_input*60.0)
    elseif(suffix.eq."hr") then
       time_value = nint(time_input*3600.0)
    elseif(suffix.eq."da") then 
       time_value = nint(time_input*86400.0)
    elseif(suffix.eq."mo") then 
       time_value = nint(time_input*2592000.0)
    endif
    
  end subroutine ParseTimeString
