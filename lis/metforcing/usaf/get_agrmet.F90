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
! !ROUTINE: get_agrmet
! \label{get_agrmet}
!
! !REVISION HISTORY:
! 29 Jul 2005; Sujay Kumar, Initial Code
! 13 Dec 2007; Marv Freimund, Do NOT abort for missing 10 files
! 14 Apr 2011; John Eylander, Added use_timestamp to create_agrmet_prod_filename
! 23 Jun 2011; Chris Franks, Use root filename from config
! 30 Apr 2012; Chris Franks, Do not update time2 data if end-of_run is true
!
! !INTERFACE:    
subroutine get_agrmet(n, findex)
! !USES:
  use LIS_coreMod, only         : LIS_rc, LIS_endofrun
  use LIS_timeMgrMod, only      : LIS_get_nstep,           &
                                  LIS_computeTimeBookEnds, &
                                  LIS_time2date,LIS_tick
  use LIS_logMod, only          : LIS_logunit
  use LIS_pluginIndices
  use AGRMET_forcingMod, only   : agrmet_struc

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!  
! !DESCRIPTION:
!  This is the entry point for calling various AGRMET analyses. 
!  At the beginning of a simulation, the code invokes call to 
!  read the most recent past data (nearest hourly interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The strategy for missing data is to go backwards up to 10 days to get
!  forcing at the same time of day. This entry point code also 
!  calls different routines based on the running mode. The AGRMET
!  mode will generate the surface analysis, merged precipitaition,
!  whereas the analysis mode will simply use the previously generated
!  data. 
!.
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_computeTimeBookEnds](\ref{LIS_computeTimeBookEnds}) \newline
!    computes the previous and past times for the specified
!    input data frequency
!  \item[readagrmetforcing](\ref{readagrmetforcing}) \newline
!    routines that perform the surface calculations, radiation
!    processing (all met forcing except precip)
!  \item[readagrmetpcpforcing](\ref{readagrmetpcpforcing}) \newline
!    routines that perform the precip analysis
!  \item[readagrmetforcinganalysis](\ref{readagrmetforcinganalysis}) \newline
!    routine to read analysis of surface fields. 
!  \item[readagrmetpcpforcinganalysis](\ref{readagrmetpcpforcinganalysis}) \newline
!    routines that reads the merged precip analysis
!  \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts time to a date format
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    computes the AGRMET read/processing times. 
!  \end{description}
!EOP

  logical :: timetoRead, timetoReadPcp
  integer :: movetime, movepcptime
  integer :: order, nforce,nstep
  real*8  :: time1,time2, timenow
  real*8  :: pcptime1,pcptime2
  integer :: ptime1,ptime2
  real    :: curr_time
  integer :: f, t
  integer :: yr1,mo1,da1,hr1,mn1,ss1,doy1,ts1
  integer :: yr2,mo2,da2,hr2,mn2,ss2,doy2,ts2
  real    :: gmt1,gmt2
  character*400 :: name 

  agrmet_struc(n)%findtime1 = 0 
  agrmet_struc(n)%findtime2 = 0 
  movetime = 0 
  movepcptime = 0 
  nforce = LIS_rc%met_nf(findex)
  ptime1 = 0 
  ptime2 = 0 
  timeToRead = .false. 
  timeToReadPcp = .false. 
  ss1 = 0

  nstep = LIS_get_nstep(LIS_rc,n)

  if(LIS_rc%runmode.ne.LIS_agrmetrunId) then
! In this mode, the code simply reads the three hour files instead of 
! doing the met forcing generation routines. 
!
     if ( LIS_get_nstep(LIS_rc,n) == 0 .or. &
          LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1) then

        agrmet_struc(n)%findtime1 = 1
        agrmet_struc(n)%findtime2 = 1
        movetime = 0 
        LIS_rc%rstflag(n) = 0 

     endif
! determine the AGRMET data times

     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr
     mn1 = LIS_rc%mn
     ss1 = 0 
     ts1 = 0 
     call LIS_tick(timenow, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
     
     yr1 = LIS_rc%yr  !previous forecast hour
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = 3*(int(real(LIS_rc%hr)/3.0))
     mn1 = 0
     ss1 = 0
     ts1 = 0
     call LIS_tick( time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, real(ts1))
     
     yr2 = LIS_rc%yr  !next forecast hour
     mo2 = LIS_rc%mo
     da2 = LIS_rc%da
     hr2 = 3*(int(real(LIS_rc%hr)/3.0))
     mn2 = 0
     ss2 = 0
     ts2 = 3*60*60
     call LIS_tick( time2, doy2, gmt2, yr2, mo2, da2, hr2, mn2, ss2, real(ts2))

     if(nstep==0 .or. nstep == 1) then 
        agrmet_struc(n)%agrmettime1 = time1
        agrmet_struc(n)%agrmettime2 = time2
     endif

     if(timenow.gt.agrmet_struc(n)%agrmettime2) then 
        movetime = 1
        agrmet_struc(n)%findtime2 = 1
     endif

     if(agrmet_struc(n)%findtime1 == 1) then 
        call create_agrmet_prod_filename(name, agrmet_struc(n)%agrmetdir, &
             agrmet_struc(n)%retroFileRoot, agrmet_struc(n)%use_timestamp, &
             yr1, mo1, da1, hr1, mn1)
        write(LIS_logunit,*) 'Getting new time1 data ', & 
             trim(name)        
        
        call readagrmetforcinganalysis(n, findex, 1, name,mo1)
     
        agrmet_struc(n)%agrmettime1 = time1
        
     endif
     
     if (movetime == 1) then 
        agrmet_struc(n)%agrmettime1 = agrmet_struc(n)%agrmettime2
        agrmet_struc(n)%findtime2 = 1
        
        do f = 1, LIS_rc%met_nf(findex)
           do t = 1, LIS_rc%ngrid(n)
              agrmet_struc(n)%metdata1(f,t) = agrmet_struc(n)%metdata2(f,t)
           end do
        end do
     end if
     
     if(agrmet_struc(n)%findtime2 == 1 .and. (.NOT. LIS_endofrun())) then 

        call create_agrmet_prod_filename(name, agrmet_struc(n)%agrmetdir, &
             agrmet_struc(n)%retroFileRoot, agrmet_struc(n)%use_timestamp, &
             yr2, mo2, da2, hr2, mn2)

        write(LIS_logunit,*) 'Getting new time2 data', & 
             trim(name)
        
        call readagrmetforcinganalysis(n, findex, 2, name,mo2)
        
        agrmet_struc(n)%agrmettime2 = time2 
     endif
     
  elseif(LIS_rc%runmode.eq.LIS_agrmetrunId) then     !YDT comment: real-time mode 
     if ( LIS_get_nstep(LIS_rc,n) == 1 .or. LIS_rc%rstflag(n) == 1) then
        curr_time = float(LIS_rc%hr)*60+float(LIS_rc%mn)+float(LIS_rc%ss/60)  
        if(mod(curr_time,60.0).eq.0) then 
           timeToRead = .true.
        endif
        curr_time = float(LIS_rc%hr)*60+float(LIS_rc%mn)+float(LIS_rc%ss/60)  
!YDT 11/27/07        if(mod(curr_time,60.0).eq.0) then 
!YDT call pcp every time step 
        if(mod(curr_time, LIS_rc%ts/60.0).eq.0) then 
           timeToReadPcp = .true.
        endif
        agrmet_struc(n)%findtime1=1
        agrmet_struc(n)%findtime2=1
        agrmet_struc(n)%metdata1 = 0
        agrmet_struc(n)%metdata2 = 0
        LIS_rc%rstflag(n) = 0
        movetime = 0 
        movepcptime = 0 
        ptime1 = 1
        ptime2 = 1
     else    
        curr_time = float(LIS_rc%hr)*60+float(LIS_rc%mn)+float(LIS_rc%ss/60)  
        if(mod(curr_time,60.0).eq.0) then 
           timeToRead = .true.
        endif
        if(timeToRead) then 
           agrmet_struc(n)%findtime2 = 1
        endif
        curr_time = float(LIS_rc%hr)*60+float(LIS_rc%mn)+float(LIS_rc%ss/60)  
!YDT 11/27/07        if(mod(curr_time,60.0).eq.0) then 
!YDT call pcp every time step 
        if(mod(curr_time, LIS_rc%ts/60.0).eq.0) then 
           timeToReadPcp = .true.
        endif
        if(timeToReadPcp) then 
           ptime2 = 1
        endif
     endif

     call LIS_computeTimeBookEnds(LIS_rc,agrmet_struc(n)%radProcessInterval,&
          time1,time2)
     call LIS_computeTimeBookEnds(LIS_rc,agrmet_struc(n)%pcpProcessInterval,&
          pcptime1,pcptime2)

! We are going to always read the previous time step for the whole hour
! no interpolation between metdata1 and metdata2 is done. 

     if(agrmet_struc(n)%findtime1.eq.1.or.agrmet_struc(n)%findtime2.eq.1) then 
        order = 1
        call readagrmetforcing(n,findex, order)
        agrmet_struc(n)%agrmettime1 = time1
     endif

     if(ptime1.eq.1.or.ptime2.eq.1) then 
        order = 1
        call readagrmetpcpforcing(n,findex, order)
        agrmet_struc(n)%agrmetpcptime1 = pcptime1
     endif
  endif
end subroutine get_agrmet

subroutine create_agrmet_prod_filename(name, rootdir, rootname, use_timestamp,&
                                       yr, mo, da, hr, mn) 

  implicit none

  character(len=*) :: name
  character(len=*) :: rootdir
  character(len=*) :: rootname
  integer, intent(in) :: use_timestamp
  integer          :: yr
  integer          :: mo
  integer          :: da
  integer          :: hr
  integer          :: mn

  character*8      :: fdate
  character*4      :: fhr

  
  write(unit=fdate,fmt='(I4.4,I2.2,I2.2)') yr, mo, da
  write(unit=fhr,fmt='(I2.2,I2.2)') hr, mn

  if (use_timestamp == 1 ) then

    name = trim(rootdir)//'/'//trim(fdate)//trim(rootname)//&
       trim(fdate)//'_DT.'//trim(fhr)//'_DF.GR1'

  else

     name = trim(rootdir)//trim(rootname)//&
        trim(fdate)//'_DT.'//trim(fhr)//'_DF.GR1'
  
  endif

end subroutine create_agrmet_prod_filename
