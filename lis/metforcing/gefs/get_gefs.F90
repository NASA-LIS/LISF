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
! !ROUTINE: get_gefs
! \label{get_gefs}
!
!
! !REVISION HISTORY:
! 7 Mar 2013: Sujay Kumar, initial specification
! 1 Jul 2019: K. Arsenault, expand support for GEFS forecasts
! 28 Jan 2021: Sujay Kumar; Updated for GEFS operational data
! 
! !INTERFACE:
subroutine get_gefs(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod,          only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,       only : LIS_tick
  use LIS_metforcingMod,    only : LIS_forc
  use LIS_logMod
  use gefs_forcingMod,      only : gefs_struc
  use LIS_constantsMod,     only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates GEFS forecast forcing. 
!
!  At the beginning of a simulation, the code reads the most recent
!  past data (nearest 6-hour interval), and the nearest future data.
!  These two datasets are used to temporally interpolate the data to
!  the current model timestep. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    index of the forcing
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the GEFS data times
!  \item[get_reforecast_filename](\ref{get_reforecast_filename}) \newline
!    Puts together GEFS Reforecast data filenames
!  \item[get_operational_filename](\ref{get_operational_filename}) \newline
!    Puts together timestamped GEFS Operational data filenames
!  \item[read\_gefs](\ref{read_gefs}) \newline
!    Reads and Interpolates GEFS forecast data to the LIS grid
!  \end{description}
!
!EOP

  integer       :: c,f,m,v,ferror
  integer       :: order
  real*8        :: time1,time2
  integer       :: yr1,mo1,da1,hr1,mn1,ss1,doy1
  integer       :: yr2,mo2,da2,hr2,mn2,ss2,doy2
  integer       :: openfile
  integer       :: hour_cycle
  integer       :: init_yr,init_mo,init_da
  integer       :: hour_cycle1
  character(len=LIS_CONST_PATH_LEN) :: filename

!  character(10), dimension(LIS_rc%met_nf(findex)), parameter :: gefs_vars = (/ &
  character(10), dimension(8) :: gefs_vars = (/ &
       'tmp_2m    ',    &     ! (1) Inst
       'spfh_2m   ',    &     ! (2) Inst
       'dswrf_sfc ',    &     ! (3) Ave
       'dlwrf_sfc ',    &     ! (4) Ave
       'ugrd_10m  ',    &     ! (5) Inst
       'vgrd_10m  ',    &     ! (6) Inst
       'pres_sfc  ',    &     ! (7) Inst
       'apcp_sfc  '     /)    ! (8) Accum
      ! Order of orcing list found in LIS

  real            :: gmt1,gmt2,ts1,ts2
  integer         :: movetime     ! 1=move time 2 data into time 1
  integer         :: valid_hour
  integer         :: fcsthr_intv
  integer         :: fcst_hour
  type(ESMF_Time) :: currTime,initTime
  type(ESMF_TimeInterval) :: dt,dt1
  integer         :: rc
! __________________________________________________________________________________


  if(LIS_rc%ts.gt.10800) then
     write(LIS_logunit,*) '[WARN] The model timestep is > forcing data timestep ...'
     write(LIS_logunit,*) '[WARN] LIS does not support this mode currently.'
     call LIS_endrun()
  endif

  gefs_struc(n)%findtime1=0
  gefs_struc(n)%findtime2=0
  movetime=0
  openfile=0 

  if( LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1 ) then  ! beginning of run
     gefs_struc(n)%findtime1=1
     gefs_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif

  ! Read in required GEFS files:
  if( gefs_struc(n)%gefs_fcsttype .eq. "Reforecast2" ) then

    ! The 0.25deg Reforecast wind files naming convention has changed to the following:
    if(gefs_struc(n)%gefs_res == 0.25) then
      gefs_vars(5) = 'ugrd_hgt'
      gefs_vars(6) = 'vgrd_hgt'
    endif

    ! First timestep of run:
    if( LIS_rc%tscount(n).eq.1 .or. &
        LIS_rc%rstflag(n).eq.1 ) then  
      ! Bookend-time record 1
      yr1=LIS_rc%yr
      mo1=LIS_rc%mo
      da1=LIS_rc%da
      hr1=LIS_rc%hr
      mn1=0
      ss1=0
      ts1=0
      call LIS_tick(gefs_struc(n)%fcsttime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

      ! Bookend-time record 2
      yr2=LIS_rc%yr    
      mo2=LIS_rc%mo
      da2=LIS_rc%da
      hr2=6
!      hr2=3   ! If applying the 3-hourly in first 72 hours
      mn2=0
      ss2=0
      ts2=0
      call LIS_tick(gefs_struc(n)%fcsttime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

      openfile=1
    endif

    ! Determine valid times when forecasts are available to be read in:
    if( gefs_struc(n)%fcst_hour <= 69 ) then
!      fcsthr_intv = 3
!      valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
      fcsthr_intv = 6
      valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
    elseif( gefs_struc(n)%fcst_hour > 69 ) then
      fcsthr_intv = 6
      valid_hour = fcsthr_intv * (LIS_rc%hr/fcsthr_intv)
    endif

    if( (valid_hour==LIS_rc%hr .and. LIS_rc%mn==0) .or. &
         openfile == 1 ) then

      ! Forecast hour condition within each file:
      if( gefs_struc(n)%fcst_hour <= 69 ) then
        gefs_struc(n)%gribrec = gefs_struc(n)%gribrec + 2
        gefs_struc(n)%fcst_hour = gefs_struc(n)%fcst_hour + 6
!        gefs_struc(n)%fcst_hour = gefs_struc(n)%fcst_hour + 3  ! 3-hour is available up to fcst hour-66
      elseif( gefs_struc(n)%fcst_hour > 69 ) then
        gefs_struc(n)%gribrec = gefs_struc(n)%gribrec + 1
        gefs_struc(n)%fcst_hour = gefs_struc(n)%fcst_hour + 6
      endif

      ! Update bookend-time record 2:
      if( LIS_rc%tscount(n).ne.1 ) then
        gefs_struc(n)%metdata1 = gefs_struc(n)%metdata2
        gefs_struc(n)%fcsttime1 = gefs_struc(n)%fcsttime2
        yr2=LIS_rc%yr
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=valid_hour
        mn2=fcsthr_intv*60    ! Backward looking
        ss2=0
        ts2=0
        call LIS_tick(gefs_struc(n)%fcsttime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
      endif

      do v=1,LIS_rc%met_nf(findex)
        do m=1,gefs_struc(n)%max_ens_members
      
          ! Obtain filenames for different GEFS run modes: 
          if( gefs_struc(n)%gefs_runmode == "forecast" ) then

            call get_reforecast_filename(filename,gefs_struc(n)%gefs_dir, &
                 gefs_struc(n)%gefs_proj, gefs_struc(n)%gefs_res, gefs_vars(v),&
                 gefs_struc(n)%init_yr, gefs_struc(n)%init_mo, &
                 gefs_struc(n)%init_da, m)

          elseif( gefs_struc(n)%gefs_runmode == "analysis" ) then

            call get_reforecast_filename(filename,gefs_struc(n)%gefs_dir, &
                 gefs_struc(n)%gefs_proj, gefs_struc(n)%gefs_res, gefs_vars(v),&
                 yr1,mo1,da1,m)
          endif
          write(LIS_logunit,*)'[INFO] Getting GEFS file ... ',trim(filename)

          ! Read in file contents:
          if( LIS_rc%tscount(n) == 1 ) then  ! Read in first two book-ends
            order = 1
            call read_gefs_reforecast( n, m, findex, order, filename, &
                                       gefs_vars(v), ferror)
            order = 2 
            call read_gefs_reforecast( n, m, findex, order, filename, &
                                       gefs_vars(v), ferror)
          else
            order = 2
            call read_gefs_reforecast( n, m, findex, order, filename, &
                                       gefs_vars(v), ferror)
          endif
        enddo
      enddo
      openfile = 0

    endif    ! End conditional for when to open forecast files


!! GEFS OPERATIONAL FILES ARE SEPARATED BY FORECAST TIME ... BELOW APPLIES
!! THIS INCLUDES THE OPERATIONAL - ARCHIVED FILES (THOUGH THESE HAVE A 
!! A DIFFERENT NAMING CONVENTION THAN THEIR RT COUNTERPART ... )

  elseif( gefs_struc(n)%gefs_fcsttype .eq. "Operational" ) then
     
     if( LIS_rc%time.gt.gefs_struc(n)%fcsttime2 ) then 
      yr1=LIS_rc%yr
      mo1=LIS_rc%mo
      da1=LIS_rc%da
      hr1=LIS_rc%hr
      mn1=0
      ss1=0
      ts1=0
      call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
      
      yr2=LIS_rc%yr    !next hour
      mo2=LIS_rc%mo
      da2=LIS_rc%da
      hr2=LIS_rc%hr
      mn2=0
      ss2=0
      ts2=3*3600
      call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
      
      movetime = 1
      gefs_struc(n)%findtime2 = 1
      
   elseif( LIS_rc%time.eq.gefs_struc(n)%fcsttime2 ) then
      yr1=LIS_rc%yr
      mo1=LIS_rc%mo
      da1=LIS_rc%da
      hr1=LIS_rc%hr
      mn1=0
      ss1=0
      ts1=0
      call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
      
      yr2=LIS_rc%yr    !next hour
      mo2=LIS_rc%mo
      da2=LIS_rc%da
      hr2=LIS_rc%hr
      mn2=0
      ss2=0
      ts2=3*3600
      call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
      
      movetime = 1
      gefs_struc(n)%findtime2 = 1
   endif
   
   ! Open up operational -- file 1
   if( gefs_struc(n)%findtime1.eq.1 ) then
      
      ferror=0
      
      hour_cycle = 6*(int(real(LIS_rc%hr)/6.0))
      
      call ESMF_TimeSet(gefs_struc(n)%initTime, &
           yy = gefs_struc(n)%init_yr, &
           mm = gefs_struc(n)%init_mo, &
           dd = gefs_struc(n)%init_da, &
           h  = hour_cycle, &
           m  = 0,&
           s  = 0,&
           rc=rc) 
      call LIS_verify(rc,'ESMF_TimeSet failed in init_gefs')
      
      
      call ESMF_TimeSet(currTime, yy=LIS_rc%yr,&
           mm = LIS_rc%mo, &
           dd = LIS_rc%da, &
           h  = LIS_rc%hr,&
           m  = LIS_rc%mn,&
           s  = LIS_rc%ss,&
           rc=rc)
      call LIS_verify(rc,'ESMF_TimeSet failed in get_gefs')
      
      dt = (currTime - gefs_struc(n)%initTime)
      dt = ESMF_TimeIntervalAbsValue(dt)
      call ESMF_TimeIntervalGet(dt,h=gefs_struc(n)%fcst_hour,rc=rc)
      call LIS_verify(rc,'ESMF_TimeIntervalGet failed in get_gefs')
      
      do m=1,gefs_struc(n)%max_ens_members
         
         if(gefs_struc(n)%fcst_hour.eq.0) then 
! zero forecast do not contain radiation, precip fields. So special
! logic is needed
            if(hour_cycle.ne.0) then
               hour_cycle1 = hour_cycle - 6
               init_yr = gefs_struc(n)%init_yr
               init_mo = gefs_struc(n)%init_mo
               init_da = gefs_struc(n)%init_da
               fcst_hour = 9 
            else
               hour_cycle1 = 18
               call ESMF_TimeIntervalSet(dt1,h=6,rc=rc)
               call LIS_verify(rc,'ESMF_TimeIntervalSet failed in get_gefs')
               
               initTime = gefs_struc(n)%initTime - dt1
               call ESMF_TimeGet(initTime, &
                    yy = init_yr, &
                    mm = init_mo, &
                    dd = init_da, &
                    h  = hour_cycle1, &
                    rc=rc) 
               call LIS_verify(rc,'ESMF_TimeSet failed in init_gefs')
               fcst_hour = 9
               
            endif
          
            call get_operational_filename(filename,gefs_struc(n)%gefs_dir,&
                 gefs_struc(n)%gefs_res, init_yr, init_mo, &
                 init_da, hour_cycle1, &
                 fcst_hour,m)           
            
         else
            call get_operational_filename(filename,gefs_struc(n)%gefs_dir,&
                 gefs_struc(n)%gefs_res, gefs_struc(n)%init_yr, gefs_struc(n)%init_mo, &
                 gefs_struc(n)%init_da, hour_cycle, &
                 gefs_struc(n)%fcst_hour,m)
            
         endif
         
         write(LIS_logunit,*)'[INFO] getting file1.. ',trim(filename)
         order = 1
         call read_gefs_operational(n,m,findex,order,filename,ferror)
         if(ferror.ge.1) then 
            gefs_struc(n)%fcsttime1=time1
         endif
      enddo
   endif   ! end of LIS_rc%findtime=1        

   ! Open up GEFS operational -- file 2
   if(gefs_struc(n)%findtime2.eq.1) then
      hour_cycle = 6*(int(real(LIS_rc%hr)/6.0)) 
      
      call ESMF_TimeSet(gefs_struc(n)%initTime, &
           yy = gefs_struc(n)%init_yr, &
           mm = gefs_struc(n)%init_mo, &
           dd = gefs_struc(n)%init_da, &
           h  = hour_cycle, &
           m  = 0,&
           s  = 0,&
           rc=rc) 
      call LIS_verify(rc,'ESMF_TimeSet failed in init_gefs')
      
      
      call ESMF_TimeSet(currTime, yy=LIS_rc%yr,&
           mm = LIS_rc%mo, &
           dd = LIS_rc%da, &
           h  = LIS_rc%hr,&
           m  = LIS_rc%mn,&
           s  = LIS_rc%ss,&
           rc=rc)
      call LIS_verify(rc,'ESMF_TimeSet failed in get_gefs')
      
      call ESMF_TimeIntervalSet(dt1,h=3,rc=rc)
      currTime = currTime + dt1

      dt = (currTime - gefs_struc(n)%initTime)
      dt = ESMF_TimeIntervalAbsValue(dt)
      call ESMF_TimeIntervalGet(dt,h=gefs_struc(n)%fcst_hour,rc=rc)
      call LIS_verify(rc,'ESMF_TimeIntervalGet failed in get_gefs')
      
      do m=1,gefs_struc(n)%max_ens_members
         
         if(gefs_struc(n)%fcst_hour.eq.0) then 
! zero forecast do not contain radiation, precip fields. So special
! logic is needed
            if(hour_cycle.ne.0) then
               hour_cycle1 = hour_cycle - 6
               init_yr = gefs_struc(n)%init_yr
               init_mo = gefs_struc(n)%init_mo
               init_da = gefs_struc(n)%init_da
               fcst_hour = 9 
            else
               hour_cycle1 = 18
               call ESMF_TimeIntervalSet(dt1,h=6,rc=rc)
               call LIS_verify(rc,'ESMF_TimeIntervalSet failed in get_gefs')
               
               initTime = gefs_struc(n)%initTime - dt1
               call ESMF_TimeGet(initTime, &
                    yy = init_yr, &
                    mm = init_mo, &
                    dd = init_da, &
                    h  = hour_cycle1, &
                    rc=rc) 
               call LIS_verify(rc,'ESMF_TimeSet failed in init_gefs')
               fcst_hour = 9
               
            endif
          
            call get_operational_filename(filename,gefs_struc(n)%gefs_dir,&
                 gefs_struc(n)%gefs_res, init_yr, init_mo, &
                 init_da, hour_cycle1, &
                 fcst_hour,m)           
            
         else
            call get_operational_filename(filename,gefs_struc(n)%gefs_dir,&
                 gefs_struc(n)%gefs_res, gefs_struc(n)%init_yr, gefs_struc(n)%init_mo, &
                 gefs_struc(n)%init_da, hour_cycle, &
                 gefs_struc(n)%fcst_hour,m)
         endif
         
         write(LIS_logunit,*)'[INFO] getting file2.. ',trim(filename)
         order = 2
         call read_gefs_operational(n,m,findex,order,filename,ferror)

         if(ferror.ge.1) then
            gefs_struc(n)%fcsttime2=time2
         endif
     enddo
   endif
  
  endif  ! End GEFS file get product block

end subroutine get_gefs

!BOP
! !ROUTINE: get_reforecast_filename
! \label{get_reforecast_filename}
!
! !INTERFACE:
subroutine get_reforecast_filename(filename,gefsdir,gefsproj,gefsres,&
                                   var,yr,mo,da,ens_id)

  implicit none
! !ARGUMENTS: 
  character(len=*), intent(out) :: filename
  character(len=*), intent(in)  :: gefsdir
  character(len=*), intent(in)  :: gefsproj
  character(len=*), intent(in)  :: var
  real, intent(in)              :: gefsres
  integer, intent(in)           :: yr,mo,da
  integer                       :: ens_id
!  Include additional flag when adding the 8-16 lead day forecast files 

!
! !DESCRIPTION:
!  This subroutine puts together GEFS file name for 
!   Reforecast versions 
! 
!  The arguments are:
!  \begin{description}
!  \item[filename]
!   name of the GEFS Reforecast forecast file
!  \item[gefsdir]
!    Name of the GEFS Reforecast forecast data directory
!  \item[gefsproj]
!    Name of the GEFS Reforecast projection type
!  \item[var]
!    Name of the GEFS Reforecast variable, used in filename
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[ens\_id]
!   id of the ensemble member
!  \end{description}
!
!EOP
  character*3 :: month
  character*6 :: ftime0
  character*6 :: ftime1
  character*8 :: ftime2
  character*3 :: fens
  integer     :: ens2

  write(unit=ftime0, fmt='(i4.4)') yr
  write(unit=ftime1, fmt='(i4.4,i2.2)') yr,mo
  write(unit=ftime2, fmt='(i4.4,i2.2,i2.2)') yr,mo,da

  if( ens_id == 1 ) then
    ens2=ens_id-1
    write(unit=fens, fmt='(a1,i2.2)') "c",ens2
  elseif( ens_id > 1 ) then
    ens2=ens_id-1
    write(unit=fens, fmt='(a1,i2.2)') "p",ens2
  endif

  if( gefsres .eq. 0.25 ) then
    filename = trim(gefsdir)//'/'//trim(ftime0)//'/'//ftime2//&
       "00/"//trim(var)//'_'//ftime2//'00_'//fens//'.grib2'
  else
    filename = trim(gefsdir)//'/'//trim(gefsproj)//'/'//trim(ftime1)//&
         "/"//trim(var)//'_'//ftime2//'00_'//fens//'.grib2'
  endif

! ** Include 2nd GEFS Reforecast2 file here to extend to 16-day forecast
!  Need to include a conditional block to switch between the two
!   Could use the local GEFS forecast hour to be the switch at fcsthour=180?

!  filename = trim(gefsdir)//'/'//trim(gefsproj)//'/'//trim(ftime1)//&
!       "/"//trim(var)//'_'//ftime2//'00_'//fens//'_t190.grib2'
!

end subroutine get_reforecast_filename



!BOP
! !ROUTINE: get_operational_filename
! \label{get_operational_filename}
!
! !INTERFACE:
subroutine get_operational_filename(filename,gefsdir,gefsres,&
     yr,mo,da,hour_cycle,fcsthour,ens_id)
  
  implicit none
! !ARGUMENTS: 
  character(len=*), intent(out) :: filename
  character(len=*), intent(in)  :: gefsdir
  real,    intent(in)           :: gefsres
  integer, intent(in)           :: yr,mo,da,hour_cycle
  integer, intent(in)           :: fcsthour
  integer, intent(in)           :: ens_id
!
! !DESCRIPTION:
!  This subroutine puts together GEFS file name for 
!   operational products
! 
!  The arguments are:
!  \begin{description}
!  \item[name]
!   name of the timestamped GEFS forecast file
!  \item[gefsdir]
!    Name of the GEFS forecast data directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hour\_cycle]
!   cycle hour used as the baseline of the forecast (00,06,12,18)
!  \item[fcsthour]
!   forecast cycle time hour 
!  \item[ens\_id]
!   id of the ensemble member
!  \end{description}
!
!EOP
  character*2 :: fcstcycle
  character*2 :: fhr
  character*3 :: fhour
  character*3 :: month
  character*8 :: ftime
  character*3 :: fens
  integer     :: ens2

  write(unit=ftime, fmt='(i4.4,i2.2,i2.2)') yr,mo,da
  write(unit=fhr,fmt='(i2.2)') hour_cycle

  write(unit=fhour, fmt='(i3.3)') fcsthour

  if( ens_id == 1 ) then
    ens2=ens_id-1
    write(unit=fens, fmt='(a1,i2.2)') "c",ens2
  elseif( ens_id > 1 ) then
    ens2=ens_id-1
    write(unit=fens, fmt='(a1,i2.2)') "p",ens2
  endif

  if( gefsres .eq. 0.25) then
    filename = trim(gefsdir)//'/'//trim(ftime)//'/'//trim(fhr)//'/'//&
         'ge'//fens//'.t'//trim(fhr)//'z.pgrb2s.0p25.f'//fhour
  else
    filename = trim(gefsdir)//'/'//trim(ftime)//'/'//trim(fhr)//'/'//&
         'ge'//fens//'.t'//trim(fhr)//'z.pgrb2a.0p50.f'//fhour
  endif
end subroutine get_operational_filename

