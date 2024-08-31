!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
!
! !INTERFACE:
!------------------------------------------------------------------------------

!NOTE:  Currently only 0.05 deg daily ECMWFforc data are supported.
subroutine readECMWFforcdata(source)

   ! Imports
   use ESMF
   use ECMWFforc_dataMod
   use LVT_coreMod
   use LVT_histDataMod
   use LVT_logMod
   use LVT_timeMgrMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif

   ! Defaults
   implicit none
   
   ! Arguments
   integer,intent(in) :: source

   ! Local variables
   integer             :: ftn 
   character*100       :: avgfile1,avgfile2, instfile
   logical             :: file_exists1,file_exists2
   real                :: prcp_in(ecmwfforcdata(source)%nc,ecmwfforcdata(source)%nr)
   real                :: prcp_in1(ecmwfforcdata(source)%nc*ecmwfforcdata(source)%nr)
   logical*1           :: lb(ecmwfforcdata(source)%nc*ecmwfforcdata(source)%nr)
   logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
   real                :: prcp(LVT_rc%lnc*LVT_rc%lnr)
   real                :: prcp_final(2,LVT_rc%lnc, LVT_rc%lnr)
   real                :: tmp_final(1,LVT_rc%lnc, LVT_rc%lnr)
   integer             :: t,c,r
   integer             :: iret,fret
   real                :: currTime
   logical             :: alarmCheck
   integer                      :: yr1, mo1, da1, hr1, mn1, ss1
   integer                      :: yr2, mo2, da2, hr2, mn2, ss2
   type(ESMF_Time)              :: time1
   type(ESMF_Time)              :: time2
   type(ESMF_TimeInterval)      :: lis_ts
   integer :: status
   integer :: precipid
   integer :: start(3), count(3)
   integer :: jda2

   ! Initialize variables
   prcp(:) = LVT_rc%udef
   prcp_final = LVT_rc%udef

   currTime = float(LVT_rc%dhr(source))*3600+ &
        60*LVT_rc%dmn(source) + LVT_rc%dss(source)
   ! EMK...Only read at 00Z
   alarmCheck = (mod(currtime,10800.0).eq.0)

   yr1 = LVT_rc%dyr(source)
   mo1 = LVT_rc%dmo(source)
   da1 = LVT_rc%dda(source)
   hr1 = LVT_rc%dhr(source)
   mn1 = 0
   ss1 = 0
   call ESMF_TimeSet(time1,yy=yr1, mm=mo1, dd=da1, &
       h=hr1,m=mn1,s=ss1,calendar=LVT_calendar,rc=status)
   call LVT_verify(status)

   ! If it is 00Z, use previous day's time level
   if (mod(currtime,10800.0).eq.0) then
      call ESMF_TimeIntervalSet(lis_ts, s = 10800, &
           rc=status)
      call LVT_verify(status)  
   else
      call ESMF_TimeIntervalSet(lis_ts, s = 0, &
           rc=status)
      call LVT_verify(status)  
   end if

  !  time2 = time1 - lis_ts
   time2 = time1 + lis_ts !Yeosang Yoon

   call ESMF_TimeGet(time2,yy=yr2, mm=mo2, dd=da2, &
        h=hr2,m=mn2,s=ss2,calendar=LVT_calendar, &
        dayOfYear=jda2, rc=status)
   call LVT_verify(status)

   if (alarmCheck) then
      
      call create_ECMWFforc_filename(ecmwfforcdata(source)%odir,&
           yr2,mo2,da2,hr2,avgfile1,avgfile2, instfile)

      inquire(file=trim(avgfile1),exist=file_exists1)
      inquire(file=trim(avgfile2),exist=file_exists2)
      
      if(file_exists1.and.file_exists2) then 
         write(LVT_logunit,*) '[INFO] Reading ECMWFforc data1 ',trim(avgfile1)
         write(LVT_logunit,*) '[INFO] Reading ECMWFforc data2 ',trim(avgfile2)
         
         call retrieve_inst_ecmwfvars(source, &
              instfile, &
              tmp_final, &
              fret)

         call retrieve_accum_ecmwfvars(source, &
              avgfile1, &
              avgfile2, &
              prcp_final, &
              LVT_rc%dhr(source), &
              fret)
         
            ! EMK...Use budget-bilinear interpolation if ECMWFforc is at 
            ! coarser resolution than the analysis grid; otherwise, use
            ! upscale averaging.
            
         
      else
         write(LVT_logunit,*)'[ERR] Read error with ECMWFforc file ', &
              trim(avgfile1)
         write(LVT_logunit,*)'[ERR] Read error with ECMWFforc file ', &
              trim(avgfile2)
         prcp_final = LVT_rc%udef
         tmp_final = LVT_rc%udef
      endif
   else 
      write(LVT_logunit,*)'[ERR] Missing ECMWFforc file ', trim(avgfile1)
      write(LVT_logunit,*)'[ERR] Missing ECMWFforc file ', trim(avgfile2)
      prcp_final = LVT_rc%udef
      tmp_final = LVT_rc%udef
   end if ! file_exists
   
!   open(100,file='test.bin',form='unformatted')
!   write(100) prcp_final(1,:,:)
!   close(100)
!   stop

   call LVT_logSingleDataStreamVar(LVT_MOC_TAIRFORC,&
        source,tmp_final,&
        vlevel=1,units='K')

   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final(1,:,:),&
        vlevel=1,units='kg/m2s')
   
   ! Now convert from kg/m2s to kg/m2
   do r=1,LVT_rc%lnr
      do c=1,LVT_rc%lnc
         if(prcp_final(1,c,r).ge.0) then 
            prcp_final(1,c,r) = prcp_final(1,c,r)*10800.0 !kg/m2
         else
            prcp_final(1,c,r) = LVT_rc%udef
         endif
      enddo ! c
   enddo ! r
   call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,prcp_final(1,:,:),&
        vlevel=1,units='kg/m2') 
      
end subroutine readECMWFforcdata

!------------------------------------------------------------------------------

subroutine create_ECMWFforc_filename(dir,yr,mo,da,hr,&
     avgfilename1,avgfilename2,instfilename)

  use LVT_timeMgrMod, only : LVT_tick
  ! Defaults
  implicit none
  
  character(len=*)   :: dir
  character(len=*)   :: avgfilename1
  character(len=*)   :: avgfilename2
  character(len=*)   :: instfilename
  integer            :: yr, mo, da, hr
  
  integer        :: remainder
  character      :: cyr*4,cmo*2,cda*2,chr*2,fhr*2, chr3*2
  character      :: ciyr*4, cimo*2, cida*2
  real*8         :: itime
  real           :: igmt
  integer        :: iyr,imo,ida,ihr,imn,iss,ts,idoy
  character(200) :: filename, file1, file2
  

  !instantaneous files 
  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mo
  write(cda, '(i2.2)') da
  write(chr, '(i2.2)') hr
  remainder = modulo(hr,2)
  if (remainder==0) then
     !=== Use analysis field, rather than forecast
     fhr = chr
  else if ((hr==03).or.(hr==15)) then
     write(fhr, '(i2.2)') hr-3
  else if ((hr==09).or.(hr==21)) then
     write(fhr, '(i2.2)') hr-9
  end if
  filename = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
  instfilename = trim(dir)//'/'//cyr//cmo//'/'//trim(filename)
  
  !time averaged filenames
  !=== Establish file name
  write(cyr, '(i4.4)') yr
  write(cmo, '(i2.2)') mo
  write(cda, '(i2.2)') da
  write(chr, '(i2.2)') hr
  if ((hr==09).or.(hr==21).or.(hr==06).or.(hr==18).or.(hr==12)) then
     write(chr3,'(i2.2)') hr-3
  endif
  if (hr==0) then
     !=== roll back one day for forecast initialization time
     iyr=yr;  imo=mo;  ida=da
     ihr=hr;  imn=0;   iss=0
     ts = -24*60*60
     call LVT_tick(itime,idoy,igmt,iyr,imo,ida,ihr,imn,iss,ts)
     write(ciyr, '(i4.4)') iyr
     write(cimo, '(i2.2)') imo
     write(cida, '(i2.2)') ida
     file2 = 'ecmwf.'//ciyr//cimo//cida//'12.'//cmo//cda//'00.1_4'
     file1 = 'ecmwf.'//ciyr//cimo//cida//'12.'//cimo//cida//'21.1_4'
  else if ((hr==03).or.(hr==15)) then
     write(fhr, '(i2.2)') hr-3
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'none'
  else if ((hr==06).or.(hr==18)) then
     write(fhr, '(i2.2)') hr-6
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr3//'.1_4'
  else if ((hr==09).or.(hr==21)) then
     write(fhr, '(i2.2)') hr-9
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr3//'.1_4'
  else if (hr==12) then
     write(fhr, '(i2.2)') hr-12
     file2 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr//'.1_4'
     file1 = 'ecmwf.'//cyr//cmo//cda//fhr//'.'//cmo//cda//chr3//'.1_4'
  end if
  if (hr==0) then
     avgfilename2 = trim(dir)//'/'//ciyr//cimo//'/'//trim(file2)
     avgfilename1 = trim(dir)//'/'//ciyr//cimo//'/'//trim(file1)
  else
     avgfilename2 = trim(dir)//'/'//cyr//cmo//'/'//trim(file2)
     avgfilename1 = trim(dir)//'/'//cyr//cmo//'/'//trim(file1)
  end if
end subroutine create_ECMWFforc_filename


!BOP
! !ROUTINE: retrieve_accum_ecmwfvars
! \label{retrieve_accum_ecmwfvars}
!
! !INTERFACE:
subroutine retrieve_accum_ecmwfvars(k, avgfile1, avgfile2, glbdata1,&
                                    hr, fret)
  use LVT_coreMod
  use LVT_logMod
  use ECMWFforc_dataMod

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer,   intent(in)        :: k
  character(len=*), intent(in) :: avgfile1
  character(len=*), intent(in) :: avgfile2
  real                         :: glbdata1(2,LVT_rc%lnc,LVT_rc%lnr)
  integer                      :: hr
  integer                      :: fret
!
! !DESCRIPTION: 
! This routine opens the corresponding ECMWF data file to extract
! the specified variable, which represents an accumulated value.
! Should be used for shortwave, longwave, lsp, and cp.
! 
!EOP
  integer,parameter  :: N_ACCUM_VARS=2
  integer            :: necmwf
  integer            :: iret,gbret
  integer            :: c,r,iv, ftn,igrib
  integer            :: kk,nvars
  character(len=20)  :: shortName 
  real               :: missingValue 
  logical            :: file_exists
  character(len=4)   :: varname(N_ACCUM_VARS)
  real,  allocatable     :: f(:)
  logical*1, allocatable :: lb(:)
  integer            :: var_index
  logical            :: var_found
  logical            :: var_status(N_ACCUM_VARS)
  logical            :: pcp_flag
  real               :: temp_val
  real               :: glbdata2(N_ACCUM_VARS,LVT_rc%lnc, LVT_rc%lnr)
  integer            :: rel_index(N_ACCUM_VARS)
  real,dimension(LVT_rc%lnc, LVT_rc%lnr) :: varfield
  logical            :: result1, result2  ! NaN check

  !=== set GRIB shortName 
  varname = (/ "lsp ","cp  " /)
  rel_index = (/ 1, 2 /) ! index of variable w.r.t. the
                               ! full list of 9 forcing variables

  necmwf = ecmwfforcdata(k)%nc*ecmwfforcdata(k)%nr

  varfield = 0
  var_status = .false.

#if (defined USE_GRIBAPI)
  !=== Set up to open file and retrieve specified field 
  fret = -1
  gbret = 0
  inquire(file=avgfile2,exist=file_exists) 
  if ( file_exists ) then 
     call grib_open_file(ftn,avgfile2,'r',iret)
     call LVT_verify(iret, 'failed to open file '//trim(avgfile2))

     call grib_count_in_file(ftn,nvars,iret)
     call LVT_verify(iret, 'error in grib_count_in_file in retrieve_avg_ecmwfvars')
     
     allocate(lb(necmwf))
     allocate(f(necmwf))

     do kk=1,nvars
        call grib_new_from_file(ftn,igrib,iret)
        call grib_get(igrib,'shortName',shortName)
        
        var_found = .false. 
        do iv=1,N_ACCUM_VARS
           if ( trim(shortName) == trim(varname(iv))) then
              var_found = .true. 
              var_index = rel_index(iv)
              var_status(iv) = .true. 
              exit
           endif
        enddo
        
        f = -9999.0
        call grib_get(igrib,'values',f,iret)
        call LVT_verify(iret, 'failed to get values in retrieve_accum_ecmwfvars')
        if ( iret .ne. 0 ) then  ! read not successful
           write(LVT_logunit,*) 'problem reading accum2 for ',shortName,iret
           call grib_release(igrib)
           deallocate(lb)
           deallocate(f) 
           fret = -1
           return 
        endif

        call grib_get(igrib,'missingValue',missingValue,iret)
        call LVT_verify(iret,'failed to get missingValue in retrieve_accum_ecmwfvars')

        call grib_release(igrib)
        call LVT_verify(iret,'error in grib_release in read_ecmwf')
        
        if(var_found) then
           lb = .false. 
           do c = 1, necmwf
              if ( f(c) .ne. missingValue ) then
                 lb(c) = .true.
              endif
           end do
           
           pcp_flag = .false.
           if(var_index.eq.1.or.var_index.eq.2) pcp_flag = .true.

           call interp_ecmwf(k,necmwf,f,lb,varfield)
           
           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 glbdata2(iv,c,r) = varfield(c,r)
              enddo
           enddo
        endif

     enddo
     call grib_close_file(ftn)
  
     deallocate(lb)
     deallocate(f)     
     
     do kk=1,N_ACCUM_VARS
        if ( .not. var_status(kk) ) then 
           write(LVT_logunit,*) &
                'ERR: Could not retrieve all entries in file: ',trim(avgfile2)
           fret = -1
           return
        endif
     enddo
     fret = 0 
  else
     write(LVT_logunit,*) 'ERR: Could not find file (A2) : ',avgfile2   
     fret = -1
  endif
  gbret = gbret + fret

  var_status = .false.

  if ( .not. ((hr == 03) .or. (hr == 15)) ) then

     inquire(file=avgfile1,exist=file_exists) 
     if ( file_exists ) then 
        call grib_open_file(ftn,avgfile1,'r',iret)
        call LVT_verify(iret, 'failed to open file '//trim(avgfile1))
        
        call grib_count_in_file(ftn,nvars,iret)
        call LVT_verify(iret, 'error in grib_count_in_file in retrieve_accum_ecmwfvars')
     
        allocate(lb(necmwf))
        allocate(f(necmwf))

        do kk=1,nvars
           call grib_new_from_file(ftn,igrib,iret)
           call grib_get(igrib,'shortName',shortName)
           
           var_found = .false. 
           do iv=1,N_ACCUM_VARS
              if ( trim(shortName) == trim(varname(iv))) then
                 var_found = .true. 
                 var_index = rel_index(iv)
                 var_status(iv) = .true. 
                 exit
              endif
           enddo
           
           f = -9999.0
           call grib_get(igrib,'values',f,iret)
           call LVT_verify(iret, 'failed to get values in retrieve_accum_ecmwfvars')
           if ( iret .ne. 0 ) then   ! read not successful
            write(LVT_logunit,*) 'problem reading accum1 for ',shortName,iret
            call grib_release(igrib)
            deallocate(lb)
            deallocate(f)     
            fret = -1
            return 
           endif

           call grib_get(igrib,'missingValue',missingValue,iret)
           call LVT_verify(iret,'failed to get missingValue in retrieve_accum_ecmwfvars')
           call grib_release(igrib)
           call LVT_verify(iret,'error in grib_release in read_ecmwf')

           if(var_found) then
              lb = .false. 
              do c = 1, necmwf
                 if ( f(c) .ne. missingValue ) then
                    lb(c) = .true.
                 endif
              end do

              pcp_flag = .false.
              if(var_index.eq.1.or.var_index.eq.2) pcp_flag = .true.

              call interp_ecmwf(k,necmwf,f,lb,varfield)              
              
              do r=1,LVT_rc%lnr
                 do c=1,LVT_rc%lnc
                    glbdata1(iv,c,r) = varfield(c,r)
                 enddo
              enddo
           endif

        enddo
        call grib_close_file(ftn)
  
        deallocate(lb)
        deallocate(f)     
        
        do kk=1,N_ACCUM_VARS
           if ( .not. var_status(kk) ) then 
              write(LVT_logunit,*)                                  &
                   'ERR: Could not retrieve all entries in file: ', &
                   trim(avgfile1)
              fret = -1
              return
           endif
        enddo
        fret = 0 
     else
        write(LVT_logunit,*) 'ERR: Could not find file (A1): ',avgfile1   
        fret = -1
     endif
     
  else
     do kk = 1, N_ACCUM_VARS
        iv = rel_index(kk)
        glbdata1(iv,:,:) = 0.0
     end do
  endif  ! not 03 or 15z
  gbret = gbret + fret
  
  if ( gbret .eq. 0 ) then  ! all fields valid
     do kk = 1, N_ACCUM_VARS
        iv = rel_index(kk)
        if ( iv == 1 .or. iv == 2 ) then ! lsp or cp
!=== Added because sometimes the f2 data is less than the f1 data
!=== creating negative precip. This ocurs only for small values of precip
!=== and the differences are on the order of 10E-4 so for intents and purposes 
!=== is not a substantial problem. This is an issue at the raw processing data 
!=== level issue that will need to be looked at more closely. Reprocessing 
!=== of entire dataset may be necessary if it is found that the problem can
!== be fixed. Will revisit. JG 3/10/2004.
!*** Added check and reset for NaN. H.B. 2/2/2015.
!*** Added check and reset for huge value. B.L. 2/5/2015.
           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 if(glbdata1(iv,c,r).ne.LVT_rc%udef.and.&
                      glbdata2(iv,c,r).ne.LVT_rc%udef) then 
!           result1 = ISNAN(glbdata1(iv,c))
!           result2 = ISNAN(glbdata2(iv,c))
!           if ( result1 .eqv. .true. ) glbdata1(iv,c) = 0.0
!           if ( result2 .eqv. .true. ) glbdata2(iv,c) = 0.0
                    temp_val = glbdata2(iv,c,r)-glbdata1(iv,c,r)
                    if (temp_val>1000000.) temp_val = 0.0         
                    if ( (glbdata2(iv,c,r)-glbdata1(iv,c,r)) < 0.0 ) then
                       glbdata1(iv,c,r) = 0.0
                    else
                       glbdata1(iv,c,r) = (temp_val * 1000.0)/(3.0*60*60)
                    endif
                 endif
              enddo
           enddo
        else
           glbdata1(iv,:,:) = (glbdata2(iv,:,:) - glbdata1(iv,:,:))/(3.0*60*60)
        endif
     enddo
     do r=1,LVT_rc%lnr
        do c=1,LVT_rc%lnc
           if(glbdata1(1,c,r).ne.LVT_rc%udef) then 
              glbdata1(1,c,r) = glbdata1(1,c,r) + glbdata1(2,c,r)
           endif
        enddo
     enddo
  else
   fret = gbret
  endif  ! gbret == 0
#endif
  
end subroutine retrieve_accum_ecmwfvars


!BOP
! !ROUTINE: retrieve_inst_ecmwfvars
! \label{retrieve_inst_ecmwfvars}
!
! !INTERFACE:
subroutine retrieve_inst_ecmwfvars(k, instfile, glbdata,fret)
  use LVT_coreMod
  use LVT_logMod
  use ECMWFforc_dataMod

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer,   intent(in)        :: k
  character(len=*), intent(in) :: instfile
  real                         :: glbdata(1,LVT_rc%lnc,LVT_rc%lnr)
  integer, intent(inout)       :: fret
!
! !DESCRIPTION: 
! This routine opens the corresponding ECMWF data file to extract
! the specified variable, which represents an instantaneous value.
! Should be used for near-surface temperature, specific humidity,
! winds, and surface pressure.
! 
!EOP
  integer,parameter  :: N_INST_VARS=1
  integer            :: necmwf
  integer            :: iret
  integer            :: c,r,iv,ftn,igrib
  integer            :: kk,nvars
  character(len=20)  :: shortName 
  real               :: missingValue 
  logical            :: file_exists
  character(len=2)   :: varname(N_INST_VARS)
  logical            :: var_status(N_INST_VARS)
  real,  allocatable     :: f(:)
  logical*1, allocatable :: lb(:)
  integer            :: var_index
  logical            :: var_found
  logical            :: pcp_flag
  integer            :: rel_index(N_INST_VARS)
  real, dimension(LVT_rc%lnc, LVT_rc%lnr)  :: varfield

  ! Then instantaneous files contain only 6 GRIB messages
  ! 1 = t
  ! 2 = q
  ! 3 = u
  ! 4 = v
  ! 5 = sp
  ! 6 = al <-- not used by LVT

  !=== set GRIB shortName 
!  varname = (/ "t ","q ","u ","v ","sp" /)
!  rel_index = (/ 1, 2, 5, 6, 7 /) ! index of variable w.r.t. the
                                  ! full list of 9 forcing variables

  varname = (/ "t "/)
  rel_index = (/ 1/) 

  necmwf = ecmwfforcdata(k)%nc*ecmwfforcdata(k)%nr

  varfield = 0
  var_status = .false.

#if (defined USE_GRIBAPI)
  !=== Set up to open file and retrieve specified field 
  fret = -1
  inquire(file=instfile,exist=file_exists) 
  if ( file_exists ) then 
     call grib_open_file(ftn,instfile,'r',iret)
     call LVT_verify(iret, 'failed to open file '//trim(instfile))

     call grib_count_in_file(ftn,nvars,iret)
     call LVT_verify(iret, 'error in grib_count_in_file in retrieve_inst_ecmwfvars')
     
     allocate(lb(necmwf))
     allocate(f(necmwf))

     do kk=1,nvars
        call grib_new_from_file(ftn,igrib,iret)
        call grib_get(igrib,'shortName',shortName)
        
        var_found = .false. 
        do iv=1,N_INST_VARS            
           if ( shortName == varname(iv) ) then
              var_found = .true. 
              var_index = rel_index(iv)
              var_status(iv) = .true. 
              exit
           endif
        enddo
        
        f = -9999.0
        call grib_get(igrib,'values',f,iret)
        call LVT_verify(iret, 'failed to get values in retrieve_inst_ecmwfvars')

        if ( iret .ne. 0 ) then   ! read not successful
         write(LVT_logunit,*) 'problem reading inst for ',shortName,iret
         call grib_release(igrib)
         deallocate(lb)
         deallocate(f)     
         fret = -1
         return 
        endif
 
        call grib_get(igrib,'missingValue',missingValue,iret)
        call LVT_verify(iret,'failed to get missingValue in retrieve_inst_ecmwfvars')
        call grib_release(igrib)
        call LVT_verify(iret,'error in grib_release in read_ecmwf')

        if(var_found) then
           lb = .false. 
           do c = 1, necmwf
              if ( f(c) .ne. missingValue ) then
                 lb(c) = .true.
              endif
           end do
           
           pcp_flag = .false.
           
           call interp_ecmwf(k,necmwf,f,lb,varfield)              

           do r=1,LVT_rc%lnr
              do c=1,LVT_rc%lnc
                 glbdata(kk,c,r) = varfield(c,r)
              enddo
           enddo

        endif   
     enddo
     call grib_close_file(ftn)
  
     deallocate(lb)
     deallocate(f)     
     
     do kk=1,N_INST_VARS
        if ( .not. var_status(kk) ) then 
           write(LVT_logunit,*) &
                'ERR: Could not retrieve all entries in file: ',trim(instfile)
           fret = -1
           return
        endif
     enddo
     fret = 0 
  else
     write(LVT_logunit,*) 'ERR: Could not find file (I): ',instfile   
     fret = -1
  endif
#endif
  
end subroutine retrieve_inst_ecmwfvars

!BOP
! !ROUTINE: interp_ecmwf
! \label{interp_ecmwf}
!
! !INTERFACE:
subroutine interp_ecmwf(source,necmwf,f,lb,varfield)
! !USES:

  use LVT_coreMod
  use ECMWFforc_dataMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: source
  integer, intent(in) :: necmwf
  real, intent(out)   :: f(necmwf)
  logical*1           :: lb(necmwf)
  real, intent(out)   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given ECMWF field 
!   to the LVT grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[necmwf]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmapgg
! \item[lis\_gds]
!  array description of the LVT grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LVT grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LVT grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
! 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
! \end{description}
!EOP
  integer :: iret
  integer :: count1,i,j,mo

  real, dimension(LVT_rc%lnc*LVT_rc%lnr) :: lis1d

  logical*1 :: lo(LVT_rc%lnc*LVT_rc%lnr)

!=== End variable declarations
!-----------------------------------------------------------------------
! Setting interpolation options (ip=0,bilinear)
! (km=1, one parameter, ibi=1,use undefined bitmap
! (needed for soil moisture and temperature only)
! Use budget bilinear (ip=3) for precip forcing fields
!-----------------------------------------------------------------------
  mo = LVT_rc%lnc*LVT_rc%lnr
!-----------------------------------------------------------------------
! Initialize output bitmap. Important for soil moisture and temp.
!-----------------------------------------------------------------------
  lo = .true.
!-----------------------------------------------------------------------  
! Interpolate to LVT grid
  !-----------------------------------------------------------------------  
  if (LVT_isAtAFinerResolution(ecmwfforcdata(source)%datares)) then
     call conserv_interp(LVT_rc%gridDesc,lb,f, &
          lo,lis1d, ecmwfforcdata(source)%nc*ecmwfforcdata(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr, ecmwfforcdata(source)%rlat, &
          ecmwfforcdata(source)%rlon, &
          ecmwfforcdata(source)%w112, ecmwfforcdata(source)%w122, &
          ecmwfforcdata(source)%w212, ecmwfforcdata(source)%w222, &
          ecmwfforcdata(source)%n112, ecmwfforcdata(source)%n122, &
          ecmwfforcdata(source)%n212, ecmwfforcdata(source)%n222, &
          LVT_rc%udef, iret)     
  else
     call upscaleByAveraging(&
          ecmwfforcdata(source)%nc*ecmwfforcdata(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr, LVT_rc%udef, &
          ecmwfforcdata(source)%n11, lb, &
          f, lo, lis1d)
  end if
!-----------------------------------------------------------------------    
! Create 2D array for main program. Also define a "soil" mask
! due to different geography between ECMWF & LDAS. For LDAS land 
! points not included in ECMWF geography dataset only.
!-----------------------------------------------------------------------    
  count1 = 0
  do j = 1, LVT_rc%lnr
     do i = 1, LVT_rc%lnc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + LVT_rc%lnc
  enddo

end subroutine interp_ecmwf
